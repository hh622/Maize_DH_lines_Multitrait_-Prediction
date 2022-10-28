#==========================================
# script07b_MultiEnvPred_DH_MegaLMM
# 2022-10-28
#=========================================

# Clean workspace
rm(list=ls())

library(rrBLUP)
library(MegaLMM)
library(data.table)
library(foreach)
library(doParallel)
library(lineup)

# set input/output dirs
DIR_gendat="../21_GenotypePreparation/"
DIR_phedat="../41_GP_DH_Prepare_Data/"
DIR_output="../45_MultiEnvPred_DH_MegaLMM_Redo/"
if(!dir.exists(DIR_output)) dir.create(DIR_output)

#--------------------------------------
# loading data
#--------------------------------------
list.files(DIR_phedat)
ls1=load(file=paste0(DIR_phedat,"Phedat_lst_MultiEnvPred.RData"))
ls1 #[1] "Phedat_lst_Trait"     "Phedat_lst_Trait_NAs"
ls2=load(file=paste(DIR_gendat,"GenotypePrepRes.RData"))
ls2 #[1] "GenodatQC"         "GenodatAmatCoding" "Genodat012Coding"  "D" "K" "G" "Gsim"

# GRM
K[1:10,1:10]

# compare Phenodat and K
sapply(list(rownames(Phedat_lst_Trait[[1]]),
            rownames(Phedat_lst_Trait_NAs[[1]][[1]])), 
       identical,
       rownames(K))
# TRUE TRUE
identical(names(Phedat_lst_Trait),
          names(Phedat_lst_Trait_NAs))# TRUE

# only use traits evaluated in >=3 Envs
Phedat_lst_Trait=Phedat_lst_Trait[names(Phedat_lst_Trait)[lengths(Phedat_lst_Trait)>=3]]
Phedat_lst_Trait_NAs=Phedat_lst_Trait_NAs[names(Phedat_lst_Trait)]
length(Phedat_lst_Trait)==length(Phedat_lst_Trait_NAs)
lengths(Phedat_lst_Trait);lengths(Phedat_lst_Trait_NAs)

#################################################################################
# set up megaLMM running parameters
#################################################################################
# set_MegaLMM_nthreads(RcppParallel::defaultNumThreads()-1)
# id = as.numeric(commandArgs(t=T)[1])
# if(is.na(id)) id = 101
# 
# # dataset = 0+id %/% 100
# dataset=3 #for grain yield
# seed = id %% 100

kfold=5
#################################################################
# read in Phenotypic and genotypic data and preparation 
#################################################################

#set up parallel computation
Ncpu=50
cl <- makeCluster(Ncpu)
registerDoParallel(cl)
ptm <- proc.time() # Start the clock!
runstat=1; runend=20

# a list of agronomic traits
AgTraits=c("bRPR","InternodeLength","InternodeDiameter",
           "InternodeCounts",
           "FreshWeight","mRPR","DryWeight","PlantHeight",
           "EarHeight","LeafLength","LeafWidth","LeafAngle")

# Prepare data
r=1; TraitTissue = names(Phedat_lst_Trait)[1]; #for testing scripts
# results=foreach(TraitTissue = names(Phedat_lst_Trait),.combine=rbind)  %:% 
#   foreach(r = runstart:runend, .combine=rbind)  %dopar% {
res=foreach(TraitTissue = names(Phedat_lst_Trait)[!(sapply(strsplit(names(Phedat_lst_Trait),"::"), 
                                                               "[",1) %in% AgTraits)],.combine=rbind)  %:% 
  foreach(r = runstat:runend, .combine=rbind)  %dopar% {
    #################################################################
    # set up training and test set 
    #################################################################
    
    # set up runID
    # runID = sprintf('../45_MultiEnvPred_DH_MegaLMM_Redo/MegaLMM_%s_%d',TraitTissue,r)
    runID = sprintf(paste0(DIR_output,'MegaLMM_%s_%d'),TraitTissue,r)
    
    # raw data
    Y  =as.matrix(Phedat_lst_Trait[[TraitTissue]])            #phedat without masking(rawdat)
    data = data.frame(Pedigree_taxa = rownames(Y),stringsAsFactors = F)
    
    # set up training data
    YNA=as.matrix(Phedat_lst_Trait_NAs[[TraitTissue]][[r]])   #phedat with masking
    
    # set up testing data
    YNA_TST=Y
    YNA_TST[!is.na(YNA)]=NA
    head(YNA_TST)
    
    # checking
    sapply(list(rownames(Y),rownames(YNA),rownames(YNA_TST)), identical, rownames(K))
    
    #################################################################
    # GP - megaLMM
    #################################################################
    library(MegaLMM)
    #----------------------------------------------------------------
    # Set the parameters of the MegaLMM model with MegaLMM_control()
    
    if(ncol(Y)==3){
      run_parameters = MegaLMM_control(
        # num_NA_groups = 0,
        drop0_tol = 1e-10,
        scale_Y = T,   # should the columns of Y be re-scaled to have mean=0 and sd=1?
        simulation = FALSE, # Are you running against simulated data (ex from a call to new_halfSib_simulation above)? If so, you can provide the setup list and it will make some QC plots
        h2_divisions = 20, # Each variance component is allowed to explain between 0% and 100% of the total variation. How many segments should the range [0,100) be divided into for each random effect?
        h2_step_size = NULL, # if NULL, all possible values of random effects are tried each iteration. If in (0,1), a new candidate set of random effect proportional variances is drawn uniformily with a range of this size
        burn = 2000,  # number of burn in samples before saving posterior samples
        thin = 2,
        K = 2 # number of factors
      )
    }else{
      run_parameters = MegaLMM_control(
        # num_NA_groups = 0,
        drop0_tol = 1e-10,
        scale_Y = T,   # should the columns of Y be re-scaled to have mean=0 and sd=1?
        simulation = FALSE, # Are you running against simulated data (ex from a call to new_halfSib_simulation above)? If so, you can provide the setup list and it will make some QC plots
        h2_divisions = 20, # Each variance component is allowed to explain between 0% and 100% of the total variation. How many segments should the range [0,100) be divided into for each random effect?
        h2_step_size = NULL, # if NULL, all possible values of random effects are tried each iteration. If in (0,1), a new candidate set of random effect proportional variances is drawn uniformily with a range of this size
        burn = 2000,  # number of burn in samples before saving posterior samples
        thin = 2,
        K = 3 # number of factors
      )
    }
    
    
    #-------------------------------------------------------------------------
    # Set the prior hyperparameters of the MegaLMM model with MegaLMM_priors()
    
    priors = MegaLMM_priors(
      tot_Y_var = list(V = 0.5,   nu = 10),      # Prior variance of trait residuals after accounting for fixed effects and factors
      tot_F_var = list(V = 18/20, nu = 20),     # Prior variance of factor traits. This is included to improve MCMC mixing, but can be turned off by setting nu very large
      Lambda_prior = list(
        sampler = sample_Lambda_prec_horseshoe,
        prop_0 = 0.1,
        delta = list(shape = 3, scale = 1),
        delta_iterations_factor = 100),
      h2_priors_resids_fun = function(h2s,n)  1,  # Function that returns the prior density for any value of the h2s vector (ie the vector of random effect proportional variances across all random effects. 1 means constant prior. Alternative: pmax(pmin(ddirichlet(c(h2s,1-sum(h2s)),rep(2,length(h2s)+1)),10),1e-10),
      h2_priors_factors_fun = function(h2s,n) 1 # See above. Another choice is one that gives 50% weight to h2==0: ifelse(h2s == 0,n,n/(n-1))
    )
    
    #----------------------------------------------------------------
    # Construct the model using 6 steps
    
    # 1. First, you need to specify the model terms using setup_model_MegaLMM
    MegaLMM_state = setup_model_MegaLMM(YNA,            # n x p data matrix
                                        ~ 1 + (1|Pedigree_taxa),
                                        data=data,         # the data.frame with information for constructing the model matrices
                                        relmat = list(Pedigree_taxa = K), # covariance matrices for the random effects. If not provided, assume uncorrelated
                                        run_parameters=run_parameters,
                                        run_ID = runID
    )
    
    # 2. After specifying the model, there is an optional step to divide the input matrix Y into chunks with similar patterns of missing data. 
    maps = make_Missing_data_map(MegaLMM_state,max_NA_groups = ncol(Y)+1,verbose=F)
    MegaLMM_state = set_Missing_data_map(MegaLMM_state,maps$Missing_data_map_list[[3]])
    
    
    # 3. Next, you specify the prior hyperparameters using the set_priors_MegaLMM funciton.
    MegaLMM_state = set_priors_MegaLMM(MegaLMM_state,priors)
    
    # 4. Next, you initialize all parameters of the model using initialize_variables_MegaLMM .
    MegaLMM_state = initialize_variables_MegaLMM(MegaLMM_state)
    MegaLMM_state_base = MegaLMM_state
    saveRDS(MegaLMM_state_base,file = sprintf('%s/MegaLMM_state_base.rds',runID))
    
    # 5. Next, you run initialize_MegaLMM to pre-calculate many large matrices that are needed during sampling.
    MegaLMM_state = initialize_MegaLMM(MegaLMM_state)
    MegaLMM_state$Posterior$posteriorSample_params = c('Lambda','F_h2','resid_h2','tot_Eta_prec','F','U_F','U_R')
    MegaLMM_state$Posterior$posteriorMean_params = c('Eta','Eta_mean')
    
    # 6. Finally, you can initialize the Posterior samples data structures using clear_Posterior. 
    MegaLMM_state = clear_Posterior(MegaLMM_state)
    
    #----------------------------------------------------------------
    # Run MCMC 
    n_iter = 100;  # how many samples to collect at once?
    # system(sprintf('rm %s/U_pred_samples.csv',runID))
    for(i  in 1:70) {
      print(sprintf('Run %d',i))
      MegaLMM_state = sample_MegaLMM(MegaLMM_state,n_iter)  # run MCMC chain n_iter iterations. grainSize is a paramter for parallelization (smaller = more parallelization)
      
      MegaLMM_state = save_posterior_chunk(MegaLMM_state)  # save any accumulated posterior samples in the database to release memory
      print(MegaLMM_state) # print status of current chain
      plot(MegaLMM_state) # make some diagnostic plots. These are saved in a pdf booklet: diagnostic_plots.pdf
      
      # set of commands to run during burn-in period to help chain converge
      if(MegaLMM_state$current_state$nrun < MegaLMM_state$run_parameters$burn || i < 50) {
        MegaLMM_state = reorder_factors(MegaLMM_state,drop_cor_threshold = 0.6) # Factor order doesn't "mix" well in the MCMC. We can help it by manually re-ordering from biggest to smallest
        MegaLMM_state = clear_Posterior(MegaLMM_state)
        print(MegaLMM_state$run_parameters$burn)
      }
    }
    
    #----------------------------------------------------------------
    # get Posterior means
    U_F = get_posterior_mean(MegaLMM_state,U_F,bychunk = T)
    Lambda = get_posterior_mean(MegaLMM_state,Lambda,bychunk = T)
    U_R = get_posterior_mean(MegaLMM_state,U_R,bychunk = T)
    # U_hat=U_F %*% Lambda + U_R  #predicted additive genetic value (based on mean values across iterations)
    U_hat = get_posterior_mean(MegaLMM_state,U_F %*% Lambda + U_R,bychunk = T) #based on each iteration then average
    Eta_mean = load_posterior_param(MegaLMM_state,'Eta_mean')#predicted total genetic value
    
    # ---------------------------------------------------------------
    # cross-validation
    library(lineup)

    df=data.frame(
      RunID=r,
      Kfold=kfold,
      TraitName=TraitTissue,
      Model = rep("megaLMM",ncol(Y)),
      numENV=ncol(YNA),
      TestENV=colnames(YNA),
      corA=corbetw2mat(YNA_TST,U_hat),
      corG=corbetw2mat(YNA_TST,Eta_mean)
      )
    write.csv(df,
              # sprintf('../45_MultiEnvPred_DH_MegaLMM/MegaLMM_%s_%d.csv',TraitTissue,r),
              sprintf(paste0(DIR_output,'MegaLMM_%s_%d.csv'),TraitTissue,r),
              row.names = F)
    unlink(sprintf(paste0(DIR_output,'MegaLMM_%s_%d'),TraitTissue,r),recursive=TRUE)
    return(df)

  }

#shut down clusters
stopCluster(cl)
proc.time() - ptm # Stop the clock

# dim(results)
# head(results)
#output results
write.csv(res,
          paste0(DIR_output,"MultiEnvPred_DH_MegaLMM_runs",runstat,"_",runend,".csv"),
          row.names = F)

