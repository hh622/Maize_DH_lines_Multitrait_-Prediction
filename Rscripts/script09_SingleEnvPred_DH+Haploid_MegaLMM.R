#================================================
# script09_SingleEnvPred_DH+Haploid_MegaLMM
# HH 2022-10-28
#================================================

# Clean workspace
rm(list=ls())

library(data.table)
library(rrBLUP)
library(tidyr)
library(lme4)
library(MegaLMM)
source("./Estimate_gcor_prediction.R")
library(MCMCglmm)
library(foreach)
library(doParallel)
library(BGLR)
library(tibble) #rownames_to_column()
# library(rstan)
# install.packages("rstan", dependencies = T)


runstart=1; runend=20

# set input/output dirs
DIR_gendat="../21_GenotypePreparation/"
DIR_phedat1="../41_GP_DH_Prepare_Data/"
DIR_phedat2="../51_GP_Hap_Prepare_Data/"
DIR_phedat3="../06_Phenodat_BLUPs_SingleEnv/"
DIR_output=paste0("../62_SingleEnvPred_DH_Haploid_MegaLMMgcor_Redo3bNEW_runs",
                  runstart,"-",runend,"/")
if(!dir.exists(DIR_output)) dir.create(DIR_output)

#--------------------------------------
# loading data
#--------------------------------------
list.files(DIR_phedat1)
ls1=load(file=paste0(DIR_phedat1,"Phedat_lst_SingleEnvPred.RData"))
ls1 #[1] "Phedat_lst_Env"
Phedat_lst_Env_DH = Phedat_lst_Env
sapply(Phedat_lst_Env, colnames)

list.files(DIR_phedat2)
ls2=load(file=paste0(DIR_phedat2,"Phedat_lst_SingleEnvPred.RData"))
ls2 #[1] "Phedat_lst_Env"
Phedat_lst_Env_Haploid = Phedat_lst_Env

ls3=load(file=paste0(DIR_gendat,"GenotypePrepRes.RData"))
ls3 #[1] "GenodatQC"         "GenodatAmatCoding" "Genodat012Coding"  "D" "K" "G" "Gsim"

list.files(DIR_phedat3)
ls4=load(file=paste0(DIR_phedat3,"Phenodat_BLUEs_BLUPs_SingleEnv.RData"))
ls4

# read in SingleEnvPred_DH_StalkQualityTraits_NAs
ls5=load(file=paste0(DIR_phedat1,"SingleEnvPred_DH_StalkQualityTraits_NAs.RData"))
ls5

# GRM
K[1:10,1:10]

#------------------------------------
# set up marker data matrix
# makrer matrix: {-1,0,1} coding
X=Genodat012Coding 
X[1:3,1:4] # ind in row, markers in col
X[is.na(X)] = matrix(colMeans(X,na.rm=T),nr = nrow(X),nc = ncol(X),byrow=T)[is.na(X)]
dim(X)
sum(is.na(X))

# compare Phenodat and K
sapply(list(rownames(Phedat_lst_Env[[1]]),rownames(X)), identical, rownames(K)) #TRUE TRUE

head(Phedat_lst_Env[[1]])
lengths(Phedat_lst_Env)
# BJ2013  BJ2014  QZ2013 SJZ2013 SJZ2014 
# 36      35      12      32      35

# merge Phedat_lst_Env_DH and Phedat_lst_Env_Haploid 
Phedat_lst_Env=list()
for(TraitName in names(Phedat_lst_Env_DH)){
  colnames(Phedat_lst_Env_DH[[TraitName]])=paste("DH",colnames(Phedat_lst_Env_DH[[TraitName]]),sep = "::")
  colnames(Phedat_lst_Env_Haploid[[TraitName]])=paste("Haploid",colnames(Phedat_lst_Env_Haploid[[TraitName]]),sep = "::")
  
  if(identical(rownames(Phedat_lst_Env_DH[[TraitName]]),
               rownames(Phedat_lst_Env_Haploid[[TraitName]]))){
    Phedat_lst_Env[[TraitName]]=data.frame(Phedat_lst_Env_DH[[TraitName]],
                                           Phedat_lst_Env_Haploid[[TraitName]],
                                           check.names = F)
  }
}

#------------------------------------------------------
# run single-Env Prediction 
# -----------------------------------------------------

# number of traits per Env
for(Env in names(Phedat_lst_Env)){
  print(
  paste("Env=",Env,
        "numTraits=",
        length(colnames(Phedat_lst_Env[[Env]])[grepl("DH",colnames(Phedat_lst_Env[[Env]]))]))
  )
}
# [1] "Env= BJ2013 numTraits= 36"
# [1] "Env= BJ2014 numTraits= 35"
# [1] "Env= QZ2013 numTraits= 12"
# [1] "Env= SJZ2013 numTraits= 26"
# [1] "Env= SJZ2014 numTraits= 35"

#set up parallel computation
Ncpu=50
cl <- makeCluster(Ncpu)
registerDoParallel(cl)
ptm <- proc.time() # Start the clock!

# set up foreach pars
Phedat_lst_Env0=Phedat_lst_Env
names(Phedat_lst_Env)
lengths(Phedat_lst_Env)

# run foreach
# Env = names(Phedat_lst_Env)[2];Trait = colnames(Phedat_lst_Env0[[Env]])[1];r=1
head(Phedat_lst_Env[[Env]])

# a list of agronomic traits (see MS Table 1)
AgTraits=c("InternodeLength","InternodeDiameter",           #take out "bRPR"
           "InternodeCounts",
           "FreshWeight","mRPR","DryWeight","PlantHeight",
           "EarHeight","LeafLength","LeafWidth","LeafAngle")
# a list of stalk quality traits  (see MS Table 1)
StalkqualTraits =c("CP","ADF","NDF","FAT","ASH","WSC","Lignin","IVDMD","Cellulose")

# run foreach loop
results=  foreach(Env = names(Phedat_lst_Env)[names(Phedat_lst_Env)!="QZ2013"],.combine=rbind)   %do% {
   foreach(r = runstart:runend, .combine=rbind)  %do% {
    # load library
    library(rrBLUP)
    library(BGLR)
    library(MegaLMM)
    library(tibble)
    
    # for each loop,start with the same Phedat_lst_Env
    Phedat_lst_Env=Phedat_lst_Env0
    sapply(Phedat_lst_Env,head)
    
    #--------------------------------------------------
    # set up a focal trait and secondary traits
    # focal trait is already set up from the foreach loop above
    # I used HTP to denote secondary traits 
    # HTP=NULL; HTP=Phedat_lst_Env[[Env]][,colnames(Phedat_lst_Env[[Env]])!=Trait]
    
    # update 9/28/22: only use agronomic traits to predict stalk quality traits
    HTP=NULL; HTP=Phedat_lst_Env[[Env]][,colnames(Phedat_lst_Env[[Env]])[!(((sapply(strsplit(colnames(Phedat_lst_Env[[Env]]),"::"), "[",2) %in% StalkqualTraits) &
                                                                              grepl("DH",colnames(Phedat_lst_Env[[Env]]))))]]
    colnames(HTP)
    
    # ----------------------------------------------------
    # set up training and test set (kfold partition)
    # 9/29/22
    # partition = sample(c('Training','Validation'),
    #                    size=nrow(Phedat_lst_Env[[Env]]),replace=T)# 50:50 partition
    # # note: do not include partition as one column of Phedat_lst_Env
    # nas = partition == 'Validation'
    
    # ----------------------------------------------------
    # set up training and test set (kfold partition)
    # 9/30/22 use the same partition as the DH-based prediction
    nas = lst_SingleEnv_NAs[[Env]][[r]]
    
    # ---------------------------------------------------------
    # set up / masking focal traits (i.e. stalk quality traits)
    SQT=NULL; SQT=Phedat_lst_Env[[Env]][,(sapply(strsplit(colnames(Phedat_lst_Env[[Env]]),"::"), "[",2) %in% StalkqualTraits) &
                                grepl("DH:",colnames(Phedat_lst_Env[[Env]]))]
    colnames(SQT)
    head(SQT)
    
    # masking SQT
    SQT_NA=SQT; SQT_NA[nas,] = NA
    colnames(SQT_NA)
    head(SQT_NA)
    identical(rownames(SQT), rownames(SQT_NA))
    identical(colnames(SQT), colnames(SQT_NA))
    

    #-----------------------------------
    # svd decompose of Knn
    Knn = K[nas,nas]
    sKnn = svd(Knn)
    round((sKnn$u %*% diag(sKnn$d) %*% t(sKnn$v)),6) [1:3,1:4] #  X = U D V'
    round(Knn,6)[1:3,1:4]
    
    #--------------------------------------
    # set up runID
    # runID = sprintf('../62_SingleEnvPred_DH_Haploid_MegaLMMgcor/MegaLMM_%s_%s_%d',Env,Trait,r)
    runID = sprintf(paste0(DIR_output,"MegaLMM_%s_%d"),Env,r)

    #################################################################
    # Set the parameters of the MegaLMM model with MegaLMM_control()
    #################################################################
    # There are several parameters that control the specific model that MegaLMM will construct. 
    # Most importantly, we need to specify:
    # - K , the number of latent factors 
    # - h2_divisions , the number of discrete values between 0 and 1 to evaluate each variance component 
    #   proportion for each random effect.
    
    # predict_MegaLMM = function(data,HTP_wide,K_year,runID,nas) {
    # run MegaLMM
    # lengths(Phedat_lst_Env)
    if(ncol(Phedat_lst_Env[[Env]])<20){
      run_parameters = MegaLMM_control(
        drop0_tol = 1e-10,
        scale_Y = T,   # should the columns of Y be re-scaled to have mean=0 and sd=1?
        h2_divisions = 20, # Each variance component is allowed to explain between 0% and 100% of the total variation. How many segments should the range [0,100) be divided into for each random effect?
        h2_step_size = NULL, # if NULL, all possible values of random effects are tried each iteration. If in (0,1), a new candidate set of random effect proportional variances is drawn uniformily with a range of this size
        burn = 2000,  # number of burn in samples before saving posterior samples
        K = 10 # number of factors
      )
    }else{
      run_parameters = MegaLMM_control(
        drop0_tol = 1e-10,
        scale_Y = T,   # should the columns of Y be re-scaled to have mean=0 and sd=1?
        h2_divisions = 20, # Each variance component is allowed to explain between 0% and 100% of the total variation. How many segments should the range [0,100) be divided into for each random effect?
        h2_step_size = NULL, # if NULL, all possible values of random effects are tried each iteration. If in (0,1), a new candidate set of random effect proportional variances is drawn uniformily with a range of this size
        burn = 2000,  # number of burn in samples before saving posterior samples
        K = 20 # number of factors
      )
    }
    
    ##########################################################################
    # Set the prior hyperparameters of the MegaLMM model with MegaLMM_priors()
    ##########################################################################
    
    # MegaLMM_prior
    priors = MegaLMM_priors(
      tot_Y_var = list(V = 0.5,   nu = 10),      # Prior variance of trait residuals after accounting for fixed effects and factors
      tot_F_var = list(V = 18/20, nu = 100000),     # Prior variance of factor traits. This is included to improve MCMC mixing, but can be turned off by setting nu very large
      Lambda_prior = list(
        sampler = sample_Lambda_prec_horseshoe,
        prop_0 = 0.1,
        delta = list(shape = 3, scale = 1),
        delta_iterations_factor = 100
      ),
      h2_priors_resids_fun = function(h2s,n) 1,  # Function that returns the prior density for any value of the h2s vector (ie the vector of random effect proportional variances across all random effects. 1 means constant prior. Alternative: pmax(pmin(ddirichlet(c(h2s,1-sum(h2s)),rep(2,length(h2s)+1)),10),1e-10),
      h2_priors_factors_fun = function(h2s,n) 1 # See above. Another choice is one that gives 50% weight to h2==0: ifelse(h2s == 0,n,n/(n-1))
    )
    
    ##########################################################################
    # Construct the model using 6 steps
    ##########################################################################
    
    # 0. Set up the phenotypes/response variables
    # MegaLMM-GBLUP (model 6 of MageLMM Wheat yield prediction)
    #  two partitions of the trait data 
    #  the first containing grain yield with the masked training set as described below
    #  the second containing all 620 hyperspectral bands with complete data
    if(identical(rownames(SQT_NA),rownames(HTP)))    Y = cbind(SQT_NA,HTP) 
    Y[1:10,1:6]
    HTP[1:10,1:6]
    sapply(list(Y,HTP), dim)
    Phedat=rownames_to_column(.data = Y,
                              var = "GID")
    
    # 1. First, you need to specify the model terms using setup_model_MegaLMM.
    # if(dir.exists(runID)) unlink(runID, recursive = T)
    MegaLMM_state = setup_model_MegaLMM(Y,            # n x p data matrix
                                        ~ 1 + (1|GID),
                                        data=Phedat,         # the data.frame with information for constructing the model matrices
                                        relmat = list(GID = K), # covariance matrices for the random effects. If not provided, assume uncorrelated
                                        run_parameters=run_parameters,
                                        run_ID = runID
    )
    
    
    # 2. After specifying the model, there is an optional step to divide the input matrix Y into chunks with similar patterns of missing data. 
    maps = make_Missing_data_map(MegaLMM_state,2,verbose=T)
    MegaLMM_state = set_Missing_data_map(MegaLMM_state,maps$Missing_data_map)
    # 3. Next, you specify the prior hyperparameters using the set_priors_MegaLMM funciton.
    MegaLMM_state = set_priors_MegaLMM(MegaLMM_state,priors)
    # 4. Next, you initialize all parameters of the model using initialize_variables_MegaLMM .
    MegaLMM_state = initialize_variables_MegaLMM(MegaLMM_state)
    # 5. Next, you run initialize_MegaLMM to pre-calculate many large matrices that are needed during sampling. 
    # This step can take very long to run and require large amounts of memory if the sample size is large and/or 
    # there are many random effect terms. A progress par is provided, 
    # but it is often important to check R's memory usage during this step to identify problems.
    MegaLMM_state = initialize_MegaLMM(MegaLMM_state)
    # MegaLMM_state$Posterior$posteriorSample_params = c('Lambda','U_F','U_R','Eta','F_h2')#,'Eta_mean')
    MegaLMM_state$Posterior$posteriorSample_params = c('Lambda')
    # MegaLMM_state$Posterior$posteriorFunctions = list(pred = 'U_R[,1] + U_F %*% Lambda[,1]')#only for one trait
    MegaLMM_state$Posterior$posteriorFunctions = list(pred = 'U_R[,1:ncol(SQT)] + U_F %*% Lambda[,1:ncol(SQT)]')#9/28/22
    
    # 6. Finally, you can initialize the Posterior samples data structures using clear_Posterior . 
    # This isn't critical if this is the first time running the model (with this run_ID ), 
    # but is good to include in case your directory structure isn't clean. 
    # Note, though, that if you were hoping to save posterior samples from an earlier run, this will erase them all!
    MegaLMM_state = clear_Posterior(MegaLMM_state)
    
    ################################
    # Run MCMC 
    ################################
    # Now (finally!) you're ready to run the MCMC chain.
    
    n_iter = 100;  # how many samples to collect at once?
    for(i  in 1:70) {
      print(sprintf('Run %d',i))
      MegaLMM_state = sample_MegaLMM(MegaLMM_state,n_iter)  # run MCMC chain n_iter iterations. grainSize is a paramter for parallelization (smaller = more parallelization)
      
      MegaLMM_state = save_posterior_chunk(MegaLMM_state)  # save any accumulated posterior samples in the database to release memory
      print(MegaLMM_state) # print status of current chain
      plot(MegaLMM_state) # make some diagnostic plots. These are saved in a pdf booklet: diagnostic_plots.pdf
      
      # set of commands to run during burn-in period to help chain converge
      if(MegaLMM_state$current_state$nrun < MegaLMM_state$run_parameters$burn || i <= 20) {
        MegaLMM_state = reorder_factors(MegaLMM_state,drop_cor_threshold = 0.6) # Factor order doesn't "mix" well in the MCMC. We can help it by manually re-ordering from biggest to smallest
        print(MegaLMM_state$run_parameters$burn)
      }
    }
    
    #####################################
    # Work with the Posterior samples
    #####################################
    source("./Estimate_gcor_prediction.R")
    U_MegaLMM = get_posterior_mean(load_posterior_param(MegaLMM_state,'pred'))
    rownames(U_MegaLMM)=gsub("::GID","",rownames( U_MegaLMM))
    head(U_MegaLMM)
    identical(rownames(U_MegaLMM),rownames(SQT))
    identical(colnames(U_MegaLMM),colnames(SQT))
    
    df=data.frame(Env=Env,
                  RunID=r,
                  GID=rownames(U_MegaLMM),
                  U_MegaLMM,
                  check.names = F,
                  stringsAsFactors = F)
    
    write.csv(df,
              sprintf(paste0(DIR_output,"MegaLMM_%s_%d.csv"),Env,r),
              row.names = F)
    
    
    unlink(sprintf(paste0(DIR_output,"MegaLMM_%s_%d"),Env,r),recursive=TRUE)
    
    
    return(df)
    }
}
    
#shut down clusters
stopCluster(cl)
proc.time() - ptm # Stop the clock

#output results
write.csv(results,
          paste0(DIR_output,"SingleEnvPred_DH_Haploid_MegaLMMgcor_Redo3b_runs",runstart,"_",runend,".csv"),
          row.names = F)
