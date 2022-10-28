#==========================================
# script07a_MultiEnvPred_DH_BGLR
# 2022-10-28
#=========================================

# Clean workspace
rm(list=ls())

library(foreach)
library(BGLR)
library(tibble)
library(lineup)
library(doParallel)

# set input/output dirs
DIR_gendat="../21_GenotypePreparation/"
DIR_phedat="../41_GP_DH_Prepare_Data/"
DIR_output="../44_MultiEnvPred_DH_BGLR/"
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
K

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

#------------------------------------------------------
# run multi-Env Prediction using BGLR
# -----------------------------------------------------

#------------------------------------------------------
# create a function cat warnings and errors
# source: https://stackoverflow.com/questions/4948361/how-do-i-save-warnings-and-errors-as-output-from-a-function
myTryCatch <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error=function(e) {
      err <<- e
      NULL
    }), warning=function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    })
  list(value=value, warning=warn, error=err)
}

#set up parallel computation
Ncpu=30
cl <- makeCluster(Ncpu)
registerDoParallel(cl)
ptm <- proc.time() # Start the clock!

#------------------------------------------------------
# run multi-Env Prediction
# set up BGLR parameters
kfold=5
nIter=20000; burnIn=5000
runstart=1; runend=20

r=1; TraitTissue = names(Phedat_lst_Trait)[1]; #for testing scripts

# clean up some parallel computing going on in the background
# unregister_dopar <- function() {
#   env <- foreach:::.foreachGlobals
#   rm(list=ls(name=env), pos=env)
# }
# unregister_dopar()

results=foreach(TraitTissue = names(Phedat_lst_Trait),.combine=rbind)  %:% 
  foreach(r = runstart:runend, .combine=rbind)  %dopar% {
    #--------------------------------------------------
    # Prepare data
    Y  =as.matrix(Phedat_lst_Trait[[TraitTissue]])            #phedat without masking(rawdat)
    YNA=as.matrix(Phedat_lst_Trait_NAs[[TraitTissue]][[r]])   #phedat with masking
    sapply(list(Y,YNA), dim)
    identical(rownames(Y),rownames(YNA))
    identical(colnames(Y),colnames(YNA))
    head(Y); head(YNA)
    str(YNA)
    
    #to create a NA matrix for a pseudo U matrix of a singular fitted model
    YNA_ALL=matrix(NA,nrow = nrow(Y), ncol = ncol(Y)) 
    
    # create YNA_TST to get cor with one command line between Obs(i.e.TST) and Pred via corbetw2mat()
    YNA_TST=Y
    YNA_TST[!is.na(YNA)]=NA
    head(YNA_TST)
    
    #------------------------------------------------
    # fitting G-by-E models with GRM
    library(BGLR)
    #1) DIAG-DIAG model
    ETA <- NULL
    ETA <-list(list(K=K,model="RKHS",Cov=list(type="DIAG")))
    message=NULL;message=myTryCatch(Gkern_D_D <- Multitrait(y=YNA, 
                                                            ETA=ETA,resCov=list(type="DIAG"),
                                                            nIter=nIter,burnIn=burnIn,
                                                            verbose = F))
    if(is.null(message$error)) {
      UHat_D_D=Gkern_D_D$ETA[[1]]$u
    }else{
      UHat_D_D=YNA_ALL
    }
    
    #2) DIAG-UN model
    ETA <- NULL
    ETA <-list(list(K=K,model="RKHS",Cov=list(type="DIAG")))
    message=NULL;message=myTryCatch(Gkern_D_UN <- Multitrait(y=YNA, 
                                                             ETA=ETA, 
                                                             nIter=nIter,burnIn=burnIn,
                                                             verbose = F))
    if(is.null(message$error)) {
      UHat_D_UN=Gkern_D_UN$ETA[[1]]$u
    }else{
      UHat_D_UN=YNA_ALL
    }
    
    #3) UN-DIAG model
    ETA <- NULL
    ETA <-list(list(K=K,model="RKHS"))#the random effect is UNstructured by default
    message=NULL;message=myTryCatch(Gkern_UN_D <-Multitrait(y=YNA, 
                                                            ETA=ETA,resCov=list(type="DIAG"), 
                                                            nIter=nIter,burnIn=burnIn,
                                                            verbose = F))
    if(is.null(message$error)) {
      UHat_UN_D=Gkern_UN_D$ETA[[1]]$u
    }else{
      UHat_UN_D=YNA_ALL
    }
    
    #4) UN-UN model (# the default setting of BGLR Multitrait: UN-UN)
    # the covariance matrix of the random effect and that of model residuals are UNstructured by default.
    ETA <- NULL
    ETA <-list(list(K=K,model="RKHS"))
    message=NULL;message=myTryCatch(Gkern_UN_UN <- Multitrait(y=YNA, 
                                                              ETA=ETA,
                                                              nIter=nIter,burnIn=burnIn,
                                                              verbose = F))
    if(is.null(message$error)) {
      UHat_UN_UN=Gkern_UN_UN$ETA[[1]]$u
    }else{
      UHat_UN_UN=YNA_ALL
    }
    
    #5) FA-DIAG model
    M <- matrix(nrow = ncol(Y), ncol = 1, TRUE)
    ETA <- NULL
    ETA <-list(list(K=K,model="RKHS",Cov=list(type="FA",M=M)))
    message=NULL;message=myTryCatch(Gkern_FA_D <-Multitrait(y=YNA, 
                                                            ETA=ETA,resCov=list(type="DIAG"), 
                                                            nIter=nIter,burnIn=burnIn,
                                                            verbose = F))
    if(is.null(message$error)) {
      UHat_FA_D=Gkern_FA_D$ETA[[1]]$u
    }else{
      UHat_FA_D=YNA_ALL
    }
    
    #6) FA-UN model
    M <- matrix(nrow = ncol(Y), ncol = 1, TRUE)
    ETA <- NULL
    ETA <-list(list(K=K,model="RKHS",Cov=list(type="FA",M=M)))
    message=NULL;message=myTryCatch(Gkern_FA_UN <- Multitrait(y=YNA, 
                                                              ETA=ETA, 
                                                              nIter=nIter,burnIn=burnIn,
                                                              verbose = F))
    if(is.null(message$error)) {
      UHat_FA_UN=Gkern_FA_UN$ETA[[1]]$u
    }else{
      UHat_FA_UN=YNA_ALL
    }
    
    #-------------------------------------------------------------------------------------------
    #calculate prediction ability only for masked data points for each ENV with one command line
    library(lineup)
    data.frame(
      RunID=r,
      Kfold=kfold,
      TraitName=TraitTissue,
      Model = c(rep("D-D",ncol(Y)),
                rep("D-UN",ncol(Y)),
                rep("UN-D",ncol(Y)),
                rep("UN-UN",ncol(Y)),
                rep("FA-D",ncol(Y)),
                rep("FA-UN",ncol(Y))),
      numENV=rep(ncol(Y),6),
      TestENV=rep(colnames(Y),6),
      PredAbility=c( corbetw2mat(YNA_TST,UHat_D_D),
                     corbetw2mat(YNA_TST,UHat_D_UN),
                     corbetw2mat(YNA_TST,UHat_UN_D),
                     corbetw2mat(YNA_TST,UHat_UN_UN),
                     corbetw2mat(YNA_TST,UHat_FA_D),
                     corbetw2mat(YNA_TST,UHat_FA_UN)))
    
    # library(lineup)
    # data.frame(
    #   RunID=r,
    #   Kfold=kfold,
    #   TraitName=TraitTissue,
    #   Model = c(rep("D-D",ncol(Y))),
    #   numENV= ncol(Y),
    #   TestENV=colnames(Y),
    #   PredAbility=corbetw2mat(YNA_TST,UHat_D_D))
  }



#shut down clusters
stopCluster(cl)
proc.time() - ptm # Stop the clock
#      user   system  elapsed 
# 1178.069 1128.816 2426.173

# dim(results)
# head(results,50)


#output results
write.csv(results,
          paste0(DIR_output,"MultiEnvPred_DH_BGLR.csv"),
          row.names = F)
# tmp=subset(results,TraitName==names(Phedat_lst_Trait)[1] & RunID==1)
# library(reshape2)
# dcast(data = tmp,
#       RunID + Kfold+TraitName+Model ~TestENV,
#       value.var ="PredAbility")

