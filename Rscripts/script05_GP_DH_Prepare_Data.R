#==========================================
# script05_GP_DH_Prepare_Data
# 2022-10-28
#=========================================

# 2022-10-28
# SingleEnv Prediction 
# - for each run, 50:50 partition DH populations for 20 random samplings
# - masking 50% of stalk quality traits for model training (not masking agronomic traits)
# - this partition will be used for both GBLUP/BGLR and MegaLMM for genomic prediction

# MultiEnv Prediction
# - used 5-fold cross validation
# - predict each trait separately considering all measurements in all Envs as independent traits
# - therefore masking 20% phenotypic data in each Env for model training

# Clean workspace
rm(list=ls())

library(foreach)
library(BGLR)
library(tibble)

# set input/output dirs
DIR_gendat="../21_GenotypePreparation/"
DIR_phedat="../06_Phenodat_BLUPs_SingleEnv/"
DIR_output="../41_GP_DH_Prepare_Data/"
if(!dir.exists(DIR_output)) dir.create(DIR_output)

#--------------------------------------
# loading data
#--------------------------------------
list.files(DIR_phedat)
ls1=load(file=paste0(DIR_phedat,"Phenodat_BLUEs_BLUPs_SingleEnv.RData"))
ls1 #[1] "blue" "blup" "hsq"
ls2=load(file=paste0(DIR_gendat,"GenotypePrepRes.RData"))
ls2 #[1] "GenodatQC"         "GenodatAmatCoding" "Genodat012Coding"  "D" "K" "G" "Gsim"
dim(GenodatQC)

#---------------------------------------------
# Prepare Phenodat and G matrix
#---------------------------------------------
# an overview of factor levels
TraitNames=unique(sapply(strsplit(colnames(blup)[-1],split='::', fixed=TRUE), `[[`, 1))
Ploidy    =unique(sapply(strsplit(colnames(blup)[-1],split='::', fixed=TRUE), `[[`, 2))
Envs      =unique(sapply(strsplit(colnames(blup)[-1],split='::', fixed=TRUE), `[[`, 3))
Tissues   =unique(sapply(strsplit(colnames(blup)[-1],split='::', fixed=TRUE), `[[`, 4))
# "4"  "5"  "6"  "7"  "8"  "11" "12" "3"  "99"
Tissues=c("4","12","99")
unique(Envs) #[1] "BJ2013"  "BJ2014"  "QZ2013"  "SJZ2013" "SJZ2014"
unique(TraitNames)

# take Phedat as blup
Phedat=NULL
Phedat=blup
Phedat$GID=paste0("GID",Phedat$GID)
dim(Phedat) #193 751
blup[1:3,1:4]
Phedat[1:3,1:4]

# truncate Phedat - to make Phedat and K has the same set of lines
Phedat=Phedat[Phedat$GID %in% rownames(K),]

# sort Phedat
rownames(Phedat)=Phedat$GID
identical(rownames(Phedat),rownames(K))#before sorting
Phedat=Phedat[rownames(K),] #sort Phedat according to GID
identical(rownames(Phedat),rownames(K))#after sorting

# exclude traits with too many missings
Phedat=Phedat[,colSums(is.na(Phedat))<30]
dim(Phedat)
colnames(Phedat)

# use median to substitute NAs
Phedat[,-1] <- sweep(as.matrix(Phedat[,-1]), MARGIN = 2, #notes: mat must be matrix, df does not work
                     STATS = apply(Phedat[,-1], 2, median, na.rm=TRUE),
                     FUN =  function(x,s) ifelse(is.na(x), s, x)
)

# select for Ploidy and Tissues/Internodes
PloidySel="DH"
TissueSel=c("4","12","99")
Phedat=Phedat[,grepl(PloidySel,colnames(Phedat)) & 
                (grepl(paste0("::",TissueSel[1]),colnames(Phedat))|
                   grepl(paste0("::",TissueSel[2]),colnames(Phedat))|
                   grepl(paste0("::",TissueSel[3]),colnames(Phedat))
                )]
dim(Phedat)
colnames(Phedat)
Phedat[1:3,1:4]

# ------------------------------------------------------------------
# summary of the DH 
(TraitNames=unique(sapply(strsplit(colnames(Phedat),split='::', fixed=TRUE), `[[`, 1)))
(Ploidy    =unique(sapply(strsplit(colnames(Phedat),split='::', fixed=TRUE), `[[`, 2)))
(Envs      =unique(sapply(strsplit(colnames(Phedat),split='::', fixed=TRUE), `[[`, 3)))
(Tissues   =unique(sapply(strsplit(colnames(Phedat),split='::', fixed=TRUE), `[[`, 4)))
(TraitTissues =unique(paste(sapply(strsplit(colnames(Phedat),split='::', fixed=TRUE), `[[`, 1),
                            sapply(strsplit(colnames(Phedat),split='::', fixed=TRUE), `[[`, 4),
                            sep = "::")))
TraitTissueslst = data.frame(Trait = sapply(strsplit(TraitTissues,split='::', fixed=TRUE), `[[`, 1),
                             Tissue = sapply(strsplit(TraitTissues,split='::', fixed=TRUE), `[[`, 2),
                             stringsAsFactors = F
)
str(TraitTissueslst)
TraitTissueslst$Tissue=as.integer(TraitTissueslst$Tissue)
TraitTissueslst = TraitTissueslst[order(TraitTissueslst$Tissue,TraitTissueslst$Trait),]
identical(TraitTissueslst$Trait[TraitTissueslst$Tissue==4],TraitTissueslst$Trait[TraitTissueslst$Tissue==12])#TRUE
write.csv(TraitTissueslst, paste0(DIR_output,"Trait_names_list.csv"),row.names = F)

# ------------------------------------------------------------------
# Organize Phedat by Env into list - for SingleEnv Prediction 
# ------------------------------------------------------------------

# create an empty list
Phedat_lst_Env=list()
Env = Envs[1]
for(Env in Envs){
  print(Env)
  Phedat_lst_Env[[Env]]=Phedat[,grepl(Env,colnames(Phedat))]
  colnames(Phedat_lst_Env[[Env]])=gsub(paste0("DH::",Env,"::"),"",colnames(Phedat_lst_Env[[Env]]))
  Phedat_lst_Env[[Env]]=Phedat_lst_Env[[Env]][,colSums(Phedat_lst_Env[[Env]])!=0]
}
length(Phedat_lst_Env)
lengths(Phedat_lst_Env)
# BJ2013  BJ2014  QZ2013 SJZ2013 SJZ2014 
# 36      35      12      32      35 
# BJ2013  BJ2014  QZ2013 SJZ2013 SJZ2014  (8/17/22: delete traits with all zero values)
# 36      35      12      26      35 
range(sapply(Phedat_lst_Env, nrow)) #187
range(sapply(Phedat_lst_Env, ncol)) #5
sapply(Phedat_lst_Env, colnames)

# ------------------------------------------------------------------
# Organize Phedat by Trait into list - for MultiEnv Prediction 
# ------------------------------------------------------------------

# create an empty list
Phedat_lst_Trait=list()
# TraitTissue = TraitTissues[1]
TraitTissue = "InternodeCounts::99" #"InternodeCounts::DH::BJ2013::99"
# TraitTissue = "InternodeDiameter::4"
for(TraitTissue in TraitTissues){
  print(TraitTissue)
  Trait = unlist(strsplit(TraitTissue,split='::'))[1]
  Tissue = paste0("::",unlist(strsplit(TraitTissue,split='::'))[2])
  
  if(sum(grepl(Trait,colnames(Phedat)))==1){ #for trait evaluated in only on Env
    
    Phedat_lst_Trait[[TraitTissue]]=data.frame(Phedat[,grepl(Trait,colnames(Phedat))&
                 grepl(Tissue,colnames(Phedat))])
    colnames(Phedat_lst_Trait[[TraitTissue]])=unlist(strsplit(colnames(Phedat)[grepl(Trait,colnames(Phedat))],split='::'))[3]
    if(colSums(Phedat_lst_Trait[[TraitTissue]])==0) Phedat_lst_Trait[[TraitTissue]]=NULL
  }else{
    Phedat_lst_Trait[[TraitTissue]]=Phedat[,grepl(Trait,colnames(Phedat))&
                                       grepl(Tissue,colnames(Phedat))]
    colnames(Phedat_lst_Trait[[TraitTissue]])=sapply(strsplit(colnames(Phedat_lst_Trait[[TraitTissue]]),
                                                        split='::', fixed=TRUE), `[[`, 3)
    # if any columns of a TraitTissue are of all zero, delete that/those columns (updated 2022/08/17)
    Phedat_lst_Trait[[TraitTissue]] = Phedat_lst_Trait[[TraitTissue]][,colSums(Phedat_lst_Trait[[TraitTissue]])!=0]
    
  }
  
  
}

length(Phedat_lst_Trait) #38
lengths(Phedat_lst_Trait)
head(Phedat_lst_Trait[[1]])
range(sapply(Phedat_lst_Trait, nrow)) #187 187
table(sapply(Phedat_lst_Trait, ncol)) 
# 1  2  3  4  5 #number of Env (old results on 8/11/2022)
# 1  2  5 20 10 #number of traits
# 1  2  3  4  5 #number of Env (updated results on 8/17/2022)
# 1  2  5 26  4 #number of traits


#------------------------------------------------------
# Making NAs matrix for single-Env Prediction
# -----------------------------------------------------

# set up foreach pars
Phedat_lst_Env0=Phedat_lst_Env

# a list of agronomic traits (see MS Table 1)
AgTraits=c("InternodeLength","InternodeDiameter",           #take out "bRPR"
           "InternodeCounts",
           "FreshWeight","mRPR","DryWeight","PlantHeight",
           "EarHeight","LeafLength","LeafWidth","LeafAngle")
# a list of stalk quality traits  (see MS Table 1)
StalkqualTraits =c("CP","ADF","NDF","FAT","ASH","WSC","Lignin","IVDMD","Cellulose")

lst_SingleEnv_NAs=list() #required
runstart=1; runend=20    #required
for(Env in names(Phedat_lst_Env))  {
  lst_SingleEnv_NAs[[Env]]=vector("list",runend)
  for(r in runstart:runend){
    
    #--------------------------------------------------
    # for each loop,start with the same Phedat_lst_Env
    Phedat_lst_Env=Phedat_lst_Env0
    head(Phedat_lst_Env[[Env]])
    
    # ----------------------------------------------------
    # set up training and test set (kfold partition)
    partition = sample(c('Training','Validation'),
                       size=nrow(Phedat_lst_Env[[Env]]),replace=T)# 50:50 partition
    # note: do not include partition as one column of Phedat_lst_Env
    nas = partition == 'Validation'
    lst_SingleEnv_NAs[[Env]][[r]]=nas
  }
}

str(lst_SingleEnv_NAs)
lst_SingleEnv_NAs[[1]]
length(lst_SingleEnv_NAs)
lengths(lst_SingleEnv_NAs)

# --------------------------------------------------------------
# Making NAs matrix for multiEnv Prediction
# --------------------------------------------------------------
# - create two-level NAs matrices
# - for each of the TraitTissue, create m NA matrices
# - ref=https://github.com/MarcooLopez/Genomic-Selection/blob/master/multi_environment.md

# set up common parameters for all TraitTissue
m <- 50           # Number of replicates
percTST <- 0.2    # Percentage of the data assigned to Testing set

Phedat_lst_Trait_NAs=list()
TraitTissue = TraitTissues[1]; k=1
for(TraitTissue in TraitTissues){
  # set up TraitTissue specific parameters
  Y <- Phedat_lst_Trait[[TraitTissue]]
  n <- nrow(Y) # Number of lines (i.e., nrow(Y))
  nEnv <- ncol(Y)         # Number of Envs  (i.e., ncol(Y))
  
  
  # Creation of seed for repeated randomizations
  set.seed(123)
  seeds <- round(seq(1E3,1E6,length=m))
  
  nTST <- round(percTST*n)
  # YNA <- vector("list",m)
  Phedat_lst_Trait_NAs[[TraitTissue]]=vector("list",m)
  nNA <- nEnv*nTST
  
  for(k in 1:m){
    set.seed(seeds[k])
    YNA0 <- Y
    if(nNA<n){ indexNA <- sample(1:n,nNA,replace=FALSE) }
    # if(nNA>=n){ #my current data does not have this case, skip
    #   nRep <- floor(nNA/n)
    #   remain <- sample(1:n,nNA%%n,replace=FALSE)
    #   a0 <- sample(1:n,n,replace=FALSE)
    #   indexNA <- rep(a0,nRep)
    #   if(length(remain)>0){
    #     a1 <- floor(length(indexNA)/nTST)*nTST
    #     a2 <- nNA - a1 - length(remain)
    #     bb <- sample(a0[!a0%in%remain],a2,replace=FALSE)
    #     noInIndexNA <- c(rep(a0,nRep-1),a0[!a0%in%bb])
    #     indexNA <- c(noInIndexNA,bb,remain)
    #     }
    #   }
  indexEnv <- rep(1:nEnv,each=nTST)
  for(j in 1:nEnv) YNA0[indexNA[indexEnv==j],j] <- NA
  # YNA[[k]] <- YNA0
  Phedat_lst_Trait_NAs[[TraitTissue]][[k]] <- YNA0
  } #end of for - k in 1:m
} #end of for-TraitTissue

length(Phedat_lst_Trait_NAs) #38
lengths(Phedat_lst_Trait_NAs)#all equals 50
head(Phedat_lst_Trait_NAs[[1]][[1]])
head(Phedat_lst_Trait_NAs[[38]][[1]])
identical(Phedat_lst_Trait_NAs[[38]][[1]],
          Phedat_lst_Trait_NAs[[38]][[2]])

# --------------------------------------------------------------
# save generated lists of phenotypic data and NAs matrices
# --------------------------------------------------------------

# save lst_SingleEnv_NAs
save(Phedat_lst_Env,
     lst_SingleEnv_NAs,
     file=paste0(DIR_output,"Phedat_lst_SingleEnvPred.RData"))

save(Phedat_lst_Trait,
     Phedat_lst_Trait_NAs,
     file=paste0(DIR_output,"Phedat_lst_MultiEnvPred.RData"))

