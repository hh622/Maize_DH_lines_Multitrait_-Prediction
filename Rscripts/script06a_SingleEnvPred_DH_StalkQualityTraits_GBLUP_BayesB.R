#=================================================================
# script06a_SingleEnvPred_DH_StalkQualityTraits_GBLUP_BayesB
# HH 2022-10-28
#=================================================================

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

runstart=1; runend=20

# set input/output dirs
DIR_gendat="../21_GenotypePreparation/"
DIR_phedat="../41_GP_DH_Prepare_Data/"
DIR_phedat2="../06_Phenodat_BLUPs_SingleEnv/"
DIR_output=paste0("../42_SingleEnvPred_DH_MegaLMMgcor_Redo3_GBLUPBayesB_runs",runstart,"-",runend,"/")
if(!dir.exists(DIR_output)) dir.create(DIR_output)

#--------------------------------------
# loading data
#--------------------------------------
list.files(DIR_phedat)
ls1=load(file=paste0(DIR_phedat,"Phedat_lst_SingleEnvPred.RData"))
ls1 #[1] "Phedat_lst_Env"
sapply(Phedat_lst_Env, colnames)

ls2=load(file=paste0(DIR_gendat,"GenotypePrepRes.RData"))
ls2 #[1] "GenodatQC"         "GenodatAmatCoding" "Genodat012Coding"  "D" "K" "G" "Gsim"
list.files(DIR_phedat2)
ls3=load(file=paste0(DIR_phedat2,"Phenodat_BLUEs_BLUPs_SingleEnv.RData"))
ls3

ls4=load(file=paste0(DIR_phedat,"SingleEnvPred_DH_StalkQualityTraits_NAs.RData"))
ls4

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
names(Phedat_lst_Env[["BJ2013"]])

#------------------------------------------------------
# run single-Env Prediction 
# -----------------------------------------------------
paste("total number of runs = ",length(Phedat_lst_Env) * sum(lengths(Phedat_lst_Env)) * 50) #36000
paste("number of runs of ",names(Phedat_lst_Env),length(Phedat_lst_Env) * lengths(Phedat_lst_Env) * 50)
# [1] "number of runs of  BJ2013 9000" 
# [2] "number of runs of  BJ2014 8750" 
# [3] "number of runs of  QZ2013 3000" 
# [4] "number of runs of  SJZ2013 6500"
# [5] "number of runs of  SJZ2014 8750"
names(Phedat_lst_Env[[1]])


# a list of agronomic traits (see MS Table 1)
AgTraits=c("InternodeLength","InternodeDiameter",           #take out "bRPR"
           "InternodeCounts",
           "FreshWeight","mRPR","DryWeight","PlantHeight",
           "EarHeight","LeafLength","LeafWidth","LeafAngle")
# a list of stalk quality traits  (see MS Table 1)
StalkqualTraits =c("CP","ADF","NDF","FAT","ASH","WSC","Lignin","IVDMD","Cellulose")

#set up parallel computation
Ncpu=50
cl <- makeCluster(Ncpu)
registerDoParallel(cl)
ptm <- proc.time() # Start the clock!

# for testing
#Env = names(Phedat_lst_Env)[1];r=1;runend=2

# set up foreach pars/assignments
Phedat_lst_Env0=Phedat_lst_Env  ##required
ResSingleEnv_UhatMegaLMM=list() #required

results=foreach(Env = names(Phedat_lst_Env),.combine = list)  %do%{
  foreach(r = runstart:runend,.combine = list)  %do% {
    #--------------------------------------------------
    # for each loop,start with the same Phedat_lst_Env
    Phedat_lst_Env=Phedat_lst_Env0
    head(Phedat_lst_Env[[Env]])
    
    # ----------------------------------------------------
    # set up training and test set (kfold partition)
    nas = lst_SingleEnv_NAs[[Env]][[r]]
    
    # ---------------------------------------------------------
    # set up / masking focal traits (i.e. stalk quality traits)
    SQT=NULL; SQT=Phedat_lst_Env[[Env]][,colnames(Phedat_lst_Env[[Env]])[sapply(strsplit(colnames(Phedat_lst_Env[[Env]]),"::"), "[",1) %in% StalkqualTraits]]
    colnames(SQT)
    head(SQT)
    
    # masking SQT
    SQT_NA=SQT; SQT_NA[nas,] = NA
    colnames(SQT_NA)
    head(SQT_NA)
    identical(rownames(SQT), rownames(SQT_NA))
    identical(colnames(SQT), colnames(SQT_NA))
    # ResSingleEnv_UhatMegaLMM[[Env]][[r]]=U_MegaLMM
    # ResSingleEnv_UhatMegaLMM[[Env]][[r]]=SQT_NA

    #-----------------------------------
    # svd decompose of Knn
    Knn = K[nas,nas]
    sKnn = svd(Knn)
    round((sKnn$u %*% diag(sKnn$d) %*% t(sKnn$v)),6) [1:3,1:4] #  X = U D V'
    round(Knn,6)[1:3,1:4]
    
    # run through each stalk quality trait
    # Trait = colnames(SQT_NA)[1]
    foreach(Trait = colnames(SQT_NA),.combine=rbind) %dopar% {
      
      # load library
      library(rrBLUP)
      library(BGLR)
      library(MegaLMM)
      library(tibble)
      
      #################################################################
      # GP - single-site GBLUP 
      #################################################################
      res_GBLUP = mixed.solve(SQT_NA[,Trait],K = K) 
      h2_rrBLUP = res_GBLUP$Vu/(res_GBLUP$Vu+res_GBLUP$Ve)
      #the h2 estimated is only using a subset (50%) of training data
      H2=hsq$hsq[hsq$TraitName==paste(unlist(strsplit(Trait, "::", fixed=TRUE))[1],"DH",Env,
                                      unlist(strsplit(Trait, "::", fixed=TRUE))[2],
                                      sep = "::")]
      # extract breeding values
      U_GBLUP=res_GBLUP$u
      
      #################################################################
      # GP - single-site BayesB 
      #################################################################
      ##BGLR - BayesB - default priors
      res_BayesB = BGLR(y = SQT_NA[,Trait],
                        ETA = list(list(X = X,model = 'BayesB')),
                        burnIn = 5000,
                        nIter = 20000,verbose=F)
      U_BayesB=X %*% res_BayesB$ETA[[1]]$b
      yhat_BayesB=res_BayesB$yHat
      # plot(yhat_BayesB,U_BayesB); abline(0,1)
      
      #####################################
      # compute gcor
      #####################################
      source("./Estimate_gcor_prediction.R")
      if(sum(is.na(SQT[,Trait])>0)>0) SQT[,Trait][is.na(SQT[,Trait])]=median(SQT[,Trait],na.rm = T)
      g_cor_GBLUP = estimate_gcor(data.frame(ID=rownames(SQT_NA)[nas],
                                             obs = SQT[,Trait][nas],
                                             pred = U_GBLUP[nas]),
                                  Knn,sKnn,method = 'MCMCglmm',normalize = T)[['g_cor']]
      
      g_cor_BayesB = estimate_gcor(data.frame(ID=rownames(SQT_NA)[nas],
                                              obs = SQT[,Trait][nas],
                                              pred = U_BayesB[nas] ),
                                   Knn,sKnn,method = 'MCMCglmm',normalize = T)[['g_cor']]
      
      #####################################
      # collect results
      #####################################
      
      df=data.frame(
        Env=Env,
        Trait=Trait,
        RunID=r,
        TraitNameFull=paste(unlist(strsplit(Trait, "::", fixed=TRUE))[1],"DH",Env,
                            unlist(strsplit(Trait, "::", fixed=TRUE))[2],sep = "::"),
        h2_rrBLUP=h2_rrBLUP,
        H2=H2,
        Model = c("GBLUP","BayesB","GBLUPgcor","BayesBgcor"),
        cor=c(cor(SQT[,Trait][nas],
                  U_GBLUP[nas],use="complete"),
              cor(SQT[,Trait][nas],
                  U_BayesB[nas],use="complete"),
              g_cor_GBLUP=g_cor_GBLUP,
              g_cor_BayesB=g_cor_BayesB
        )
      )
      
      write.csv(df,
                sprintf(paste0(DIR_output,"GBLUP_BayesB_%s_%s_%d.csv"),Env,Trait,r),
                row.names = F)
      
      
      return(df)

      
    }#foreach-Trait
    
  } #foreach-RunID/resampling
} #foreach-Env

#output results
write.csv(results,
          paste0(DIR_output,"SingleEnvPred_DH_redo3_GBLUP_BayesB_runs",runstart,"_",runend,".csv"),
          row.names = F)
