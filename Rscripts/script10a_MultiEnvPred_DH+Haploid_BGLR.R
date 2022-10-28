#==========================================
# script10a_MultiEnvPred_DH+Haploid_BGLR
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
DIR_phedat1="../41_GP_DH_Prepare_Data/"
DIR_phedat2="../51_GP_Hap_Prepare_Data/"
DIR_output="../64_MultiEnvPred_DH_Haploid_BGLR_Redo/"
if(!dir.exists(DIR_output)) dir.create(DIR_output)

#--------------------------------------
# loading data
#--------------------------------------
list.files(DIR_phedat1)
ls1=load(file=paste0(DIR_phedat1,"Phedat_lst_MultiEnvPred.RData"))
ls1 #[1] "Phedat_lst_Trait"     "Phedat_lst_Trait_NAs"
Phedat_lst_Trait_DH = Phedat_lst_Trait
Phedat_lst_Trait_NAs_DH = Phedat_lst_Trait_NAs
sapply(Phedat_lst_Trait_DH, colnames)

list.files(DIR_phedat2)
ls2=load(file=paste0(DIR_phedat2,"Phedat_lst_MultiEnvPred20220820.RData"))
ls2 #[1] "Phedat_lst_Trait" "Phedat_lst_Trait_NAs" "Phedat_lst_Trait_TST"
Phedat_lst_Trait_Haploid = Phedat_lst_Trait

ls3=load(file=paste(DIR_gendat,"GenotypePrepRes.RData"))
ls3 #[1] "GenodatQC"         "GenodatAmatCoding" "Genodat012Coding"  "D" "K" "G" "Gsim"

# GRM
K[1:10,1:10]

# merge Phedat_lst_Trait_DH and Phedat_lst_Trait_Haploid 
Phedat_lst_Trait=list()
TraitTissue = names(Phedat_lst_Trait_DH)[1]
identical(names(Phedat_lst_Trait_DH),names(Phedat_lst_Trait_Haploid))
for(TraitTissue in names(Phedat_lst_Trait_DH)){
  colnames(Phedat_lst_Trait_DH[[TraitTissue]])=paste("DH",colnames(Phedat_lst_Trait_DH[[TraitTissue]]),sep = "::")
  colnames(Phedat_lst_Trait_Haploid[[TraitTissue]])=paste("Haploid",colnames(Phedat_lst_Trait_Haploid[[TraitTissue]]),sep = "::")
  # head(Phedat_lst_Trait_DH[[TraitTissue]])
  # head(Phedat_lst_Trait_Haploid[[TraitTissue]])
  if(identical(rownames(Phedat_lst_Trait_DH[[TraitTissue]]),
               rownames(Phedat_lst_Trait_Haploid[[TraitTissue]]))){
    Phedat_lst_Trait[[TraitTissue]]=data.frame(Phedat_lst_Trait_DH[[TraitTissue]],
                                           Phedat_lst_Trait_Haploid[[TraitTissue]],
                                           check.names = F)
  }
}

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
Ncpu=50
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

# a list of agronomic traits
AgTraits=c("bRPR","InternodeLength","InternodeDiameter",
           "InternodeCounts",
           "FreshWeight","mRPR","DryWeight","PlantHeight",
           "EarHeight","LeafLength","LeafWidth","LeafAngle")
# names(Phedat_lst_Trait)[!(sapply(strsplit(names(Phedat_lst_Trait),"::"), "[",1) %in% AgTraits)]

results=foreach(TraitTissue = names(Phedat_lst_Trait)[!(sapply(strsplit(names(Phedat_lst_Trait),"::"), "[",1) %in% AgTraits)],
                .combine=rbind)  %:% 
  foreach(r = runstart:runend, .combine=rbind)  %dopar% {
    library(BGLR)
    library(lineup)
    #--------------------------------------------------
    # Prepare data
    Y  =as.matrix(Phedat_lst_Trait[[TraitTissue]])            #phedat without masking(rawdat)
    dim(Y)
    head(Y)
    
    # prepare YNA
    YNA=YNA0=NULL
    YNA0=as.matrix(Phedat_lst_Trait_NAs_DH[[TraitTissue]][[r]])   #phedat with masking
    colnames(YNA0)=paste("DH",colnames(YNA0),sep = "::")
    if(identical(rownames(YNA0),
                 rownames(Phedat_lst_Trait_Haploid[[TraitTissue]]) )){
      YNA=data.frame(YNA0,
                     Phedat_lst_Trait_Haploid[[TraitTissue]],
                     check.names = F)
      }
    YNA=as.matrix(YNA)
    
    # checking identidity of Y and YNA for rownames and colnames
    sapply(list(Y,YNA), dim)
    identical(rownames(Y),rownames(YNA))
    identical(colnames(Y),colnames(YNA))
    head(Y,50); head(YNA,50)
    corbetw2mat(Y,YNA)
    
    #to create a NA matrix for a pseudo U matrix of a singular fitted model
    YNA_ALL=matrix(NA,nrow = nrow(Y), ncol = ncol(Y)) 
    
    # create YNA_TST to get cor with one command line between Obs(i.e.TST) and Pred via corbetw2mat()
    YNA_TST=Y[,grepl("DH",colnames(Y))]
    YNA_TST[!is.na(YNA[,grepl("DH",colnames(YNA))])]=NA
    head(YNA_TST)
    colSums(!is.na(YNA_TST))
    
    sapply(list(rownames(Y),rownames(YNA),rownames(YNA_TST)), identical,rownames(K))
    sapply(list(Y,YNA,YNA_TST,K), dim)
    # [,1] [,2] [,3] [,4]
    # [1,]  187  187  187  187
    # [2,]   10   10    5  187
    #------------------------------------------------
    # fitting G-by-E models with GRM
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
    str(Gkern_FA_D)
    
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
    head(Y)
    head(YNA_TST)
    
    df=data.frame(
      RunID=r,
      Kfold=kfold,
      TraitName=TraitTissue,
      Model = c(rep("D-D",ncol(YNA_TST)),
                rep("D-UN",ncol(YNA_TST)),
                rep("UN-D",ncol(YNA_TST)),
                rep("UN-UN",ncol(YNA_TST)),
                rep("FA-D",ncol(YNA_TST)),
                rep("FA-UN",ncol(YNA_TST))),
      numENV=rep(ncol(YNA_TST),6),
      TestENV=rep(colnames(YNA_TST),6),
      PredAbility=c( corbetw2mat(YNA_TST,UHat_D_D[,grepl("DH",colnames(Y))]),
                     corbetw2mat(YNA_TST,UHat_D_UN[,grepl("DH",colnames(Y))]),
                     corbetw2mat(YNA_TST,UHat_UN_D[,grepl("DH",colnames(Y))]),
                     corbetw2mat(YNA_TST,UHat_UN_UN[,grepl("DH",colnames(Y))]),
                     corbetw2mat(YNA_TST,UHat_FA_D[,grepl("DH",colnames(Y))]),
                     corbetw2mat(YNA_TST,UHat_FA_UN[,grepl("DH",colnames(Y))]) ))
    
 
    
    write.csv(df,
              sprintf(paste0(DIR_output,'MultiEnvPred_DH_Haploid_BGLR_Redo_%s_%d.csv'),TraitTissue,r),
              row.names = F)
    return(df)

  }



#shut down clusters
stopCluster(cl)
proc.time() - ptm # Stop the clock
#      user   system  elapsed 
# 1178.069 1128.816 2426.173

dim(results)
head(results,50)

#output results
write.csv(results,
          paste0(DIR_output,"MultiEnvPred_DH_Haploid_BGLR_Redo.csv"),
          row.names = F)

