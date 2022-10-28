#==========================================
# script08_GP_DH+Haploid_Prepare_Data
# 2022-10-28
#=========================================

# Clean workspace
rm(list=ls())

library(foreach)
library(BGLR)
library(tibble)

# set input/output dirs
DIR_gendat="../21_GenotypePreparation/"
DIR_phedat="../06_Phenodat_BLUPs_SingleEnv/"
DIR_output="../51_GP_Hap_Prepare_Data/"
if(!dir.exists(DIR_output)) dir.create(DIR_output)

#--------------------------------------
# loading data
#--------------------------------------
list.files(DIR_phedat)
ls1=load(file=paste0(DIR_phedat,"Phenodat_BLUEs_BLUPs_SingleEnv.RData"))
ls1 #[1] "blue" "blup" "hsq"
ls2=load(file=paste(DIR_gendat,"GenotypePrepRes.RData"))
ls2 #[1] "GenodatQC"         "GenodatAmatCoding" "Genodat012Coding"  "D" "K" "G" "Gsim"

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
colnames(blup)
table(sapply(strsplit(colnames(blup)[grepl("DH",colnames(blup))],split = "::"),"[[",3))
# BJ2013  BJ2014  QZ2013 SJZ2013 SJZ2014 
# 111      35     105     107      35 
table(sapply(strsplit(colnames(blup)[grepl("Haploid",colnames(blup))],split = "::"),"[[",3))
# BJ2013  BJ2014  QZ2013 SJZ2013 SJZ2014 
# 123      35      54     110      35 

# truncate Phedat - to make Phedat and K has the same set of lines
Phedat=Phedat[Phedat$GID %in% rownames(K),]

# sort Phedat
rownames(Phedat)=Phedat$GID
identical(rownames(Phedat),rownames(K))#before sorting
Phedat=Phedat[rownames(K),] #sort Phedat according to GID
identical(rownames(Phedat),rownames(K))#after sorting

# exclude traits with too many missings
Phedat=Phedat[,colSums(!is.na(Phedat))>50]
dim(Phedat)
colnames(Phedat)

# select for Tissues/Internodes
TissueSel=c("4","12","99")
Phedat=Phedat[,(grepl(paste0("::",TissueSel[1]),colnames(Phedat))|
                   grepl(paste0("::",TissueSel[2]),colnames(Phedat))|
                   grepl(paste0("::",TissueSel[3]),colnames(Phedat))
                )]
dim(Phedat)
colnames(Phedat)
Phedat[1:3,1:4]

# ------------------------------------------------------------------
# summary of the Hap 
(TraitNames=unique(sapply(strsplit(colnames(Phedat),split='::', fixed=TRUE), `[[`, 1)))
(Ploidy    =unique(sapply(strsplit(colnames(Phedat),split='::', fixed=TRUE), `[[`, 2)))
(Envs      =unique(sapply(strsplit(colnames(Phedat),split='::', fixed=TRUE), `[[`, 3)))
(Tissues   =unique(sapply(strsplit(colnames(Phedat),split='::', fixed=TRUE), `[[`, 4)))
(TraitTissues =unique(paste(sapply(strsplit(colnames(Phedat),split='::', fixed=TRUE), `[[`, 1),
                            sapply(strsplit(colnames(Phedat),split='::', fixed=TRUE), `[[`, 4),
                            sep = "::")))

# ------------------------------------------------------------------
# Organize Phedat by Env into list - for SingleEnv Prediction 
# ------------------------------------------------------------------

# create an empty list
Phedat_lst_Env=list()
Env = Envs[3]
for(Env in Envs){
  print(Env)
  Phedat_lst_Env[[Env]]=Phedat[,grepl(Env,colnames(Phedat))]
  dim(Phedat_lst_Env[[Env]])
  head(Phedat_lst_Env[[Env]])
  # modify colnames
  colnames(Phedat_lst_Env[[Env]])=gsub(paste0(Env,"::"),"",colnames(Phedat_lst_Env[[Env]]))
  # delete traits of which all data points are zero
  Phedat_lst_Env[[Env]]=Phedat_lst_Env[[Env]][,colSums(abs(Phedat_lst_Env[[Env]]),na.rm = T)!=0]
  
  head(Phedat_lst_Env[[Env]])
  
  grepl()
  
  colSums(abs(Phedat_lst_Env[[Env]]),na.rm = T)
}
length(Phedat_lst_Env)
lengths(Phedat_lst_Env)
# BJ2013  BJ2014  QZ2013 SJZ2013 SJZ2014 (Haploid, 8/19/22)
# 33      35       6      31      35 

# BJ2013  BJ2014  QZ2013 SJZ2013 SJZ2014 (8/17/22: DH, delete traits with all zero values)
# 36      35      12      26      35 

# for DH + Haploid
lengths(Phedat_lst_Env)
# BJ2013  BJ2014  QZ2013 SJZ2013 SJZ2014 (DH+Haploid, 8/19/22)
# 69      70      36      58      70 
sapply(Phedat_lst_Env, function(x) sum(grepl("Haploid",colnames(x))))
# BJ2013  BJ2014  QZ2013 SJZ2013 SJZ2014 
# 33      35       6      31      35 
sapply(Phedat_lst_Env, function(x) sum(grepl("DH",colnames(x))))
# BJ2013  BJ2014  QZ2013 SJZ2013 SJZ2014 
# 36      35      30      27      35 

range(sapply(Phedat_lst_Env, nrow)) #187
range(sapply(Phedat_lst_Env, ncol)) #36 70
sapply(Phedat_lst_Env, colnames)

# ------------------------------------------------------------------
# Organize Phedat by Trait into list - for MultiEnv Prediction 
# ------------------------------------------------------------------

# create an empty list
Phedat_lst_Trait=list()
TraitTissue = TraitTissues[1]
# TraitTissue = "InternodeCounts::99" #"InternodeCounts::DH::BJ2013::99"
# TraitTissue = "InternodeDiameter::4"
for(TraitTissue in TraitTissues){
  print(TraitTissue)
  Trait = unlist(strsplit(TraitTissue,split='::'))[1]
  Tissue = paste0("::",unlist(strsplit(TraitTissue,split='::'))[2])
  
  if(sum(grepl(Trait,colnames(Phedat)))==1){ #for trait evaluated in only one Env
    Phedat_lst_Trait[[TraitTissue]]=data.frame(Phedat[,grepl(Trait,colnames(Phedat))&
                 grepl(Tissue,colnames(Phedat))])
    colnames(Phedat_lst_Trait[[TraitTissue]])=
      paste(sapply(strsplit(colnames(Phedat_lst_Trait[[TraitTissue]]),split='::', fixed=TRUE), `[[`, 2),
            sapply(strsplit(colnames(Phedat_lst_Trait[[TraitTissue]]),split='::', fixed=TRUE), `[[`, 3),
            sep = "::")
    if(colSums(Phedat_lst_Trait[[TraitTissue]],na.rm = T)==0) Phedat_lst_Trait[[TraitTissue]]=NULL
  }else{
    Phedat_lst_Trait[[TraitTissue]]=Phedat[,grepl(Trait,colnames(Phedat))&
                                       grepl(Tissue,colnames(Phedat))]
    colnames(Phedat_lst_Trait[[TraitTissue]])=
      paste(sapply(strsplit(colnames(Phedat_lst_Trait[[TraitTissue]]),split='::', fixed=TRUE), `[[`, 2),
            sapply(strsplit(colnames(Phedat_lst_Trait[[TraitTissue]]),split='::', fixed=TRUE), `[[`, 3),
            sep = "::")
    # if any columns of a TraitTissue are of all zero, delete that/those columns (updated 2022/08/17)
    Phedat_lst_Trait[[TraitTissue]] = Phedat_lst_Trait[[TraitTissue]][,colSums(Phedat_lst_Trait[[TraitTissue]],na.rm = T)!=0]
  }
}

length(Phedat_lst_Trait) #38
lengths(Phedat_lst_Trait)
head(Phedat_lst_Trait[[1]])
range(sapply(Phedat_lst_Trait, nrow)) #187 187
table(sapply(Phedat_lst_Trait, ncol)) 
# 1  2  3  4  5 #number of Env (Haploid)
# 1  3  7 23  4 #number of traits
# for DH
# 1  2  3  4  5 #number of Env (DH)
# 1  2  5 26  4 #number of traits


sapply(Phedat_lst_Trait, colnames)

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
Phedat_lst_Trait_TST=list()
TraitTissue = TraitTissues[1]; k=1

for(TraitTissue in TraitTissues){
  # set up TraitTissue specific parameters
  Y <- Phedat_lst_Trait[[TraitTissue]]
  n <- nrow(Y[,grepl("DH",colnames(Y))])         # Number of lines (i.e., nrow(Y))
  nEnv <- ncol(Y[,grepl("DH",colnames(Y))])      # Number of Envs  (i.e., ncol(Y))
  Phedat_lst_Trait_NAs[[TraitTissue]]=vector("list",m)
  Phedat_lst_Trait_TST[[TraitTissue]]=vector("list",m)
  head(Y)
  
  # for DH, use median to substitute NAs
  Y[,grepl("DH",colnames(Y))] <- sweep(as.matrix(Y[,grepl("DH",colnames(Y))]), MARGIN = 2, #notes: mat must be matrix, df does not work
                       STATS = apply(Y[,grepl("DH",colnames(Y))], 2, median, na.rm=TRUE),
                       FUN =  function(x,s) ifelse(is.na(x), s, x)
  )
  
  # Creation of seed for repeated randomizations
  set.seed(123)
  seeds <- round(seq(1E3,1E6,length=m))
  nTST <- round(percTST*n)
  # YNA <- vector("list",m)
  Phedat_lst_Trait_NAs[[TraitTissue]]=vector("list",m)
  nNA <- nEnv*nTST

  # #------------------------------------------------------------------------------
  # #for DH only - mask percTST/20% of phenotypic values (ignoring NA in raw data)
  # for(k in 1:m){
  #   set.seed(seeds[k])
  #   YNA <- Y
  #   if(nNA<n){ indexNA <- sample(1:n,nNA,replace=FALSE) }
  #   indexEnv <- rep(1:nEnv,each=nTST)
  #   for(j in 1:nEnv) YNA[,grepl("DH",colnames(YNA))][indexNA[indexEnv==j],j] <- NA
  #   # YNA[[k]] <- YNA
  #   Phedat_lst_Trait_NAs[[TraitTissue]][[k]] <- YNA
  # } #end of for - k in 1:m
  
  #------------------------------------------------------------
  #mask percTST/20% of phenotypic values where they are not NAs
  for(k in 1:m){ #for each replication/run
    mask = matrix(F,nrow = nrow(Y),ncol = ncol(Y))#for each run, mask matrix should be the same
    for(i in 1:nEnv) { #for each Env within each run
      obs = which(!is.na(Y[,i]))
      mask[sample(obs,percTST*length(obs)),i] = T 
    }
    
    # create training data
    YNA = Y
    YNA[mask] = NA
    YNA = as.matrix(YNA)
    
    # create test data
    Y_TST = Y
    Y_TST[is.na(Y_TST)] = 999
    head(Y_TST)
    Y_TST[!mask]=NA
    
    # head(mask)
    # head(!mask)
    # colSums(mask)
    # colSums(!is.na(Y))/5
    # colSums(!is.na(Y))
    # colSums(!is.na(YNA))
    # colSums(!is.na(Y_TST))
    
    
    # add masked Phenodat to list
    Phedat_lst_Trait_NAs[[TraitTissue]][[k]] <- YNA
    Phedat_lst_Trait_TST[[TraitTissue]][[k]] <- Y_TST
  } #end of for - k in 1:m
} #end of for - TraitTissue

length(Phedat_lst_Trait_NAs) #38
lengths(Phedat_lst_Trait_NAs)#all equals 50
head(Phedat_lst_Trait_NAs[[1]][[1]])
head(Phedat_lst_Trait_NAs[[38]][[1]])
sapply(Phedat_lst_Trait_NAs[[38]], function(x) colSums(is.na(x)))
sapply(Phedat_lst_Trait_TST[[38]], function(x) colSums(!is.na(x)))

# --------------------------------------------------------------
# save generated lists of phenotypic data and NAs matrices
# --------------------------------------------------------------
# save(Phedat_lst_Env,file=paste0(DIR_output,"Phedat_lst_SingleEnvPred.RData"))
save(Phedat_lst_Trait,
     Phedat_lst_Trait_NAs,
     Phedat_lst_Trait_TST,
     file=paste0(DIR_output,"Phedat_lst_MultiEnvPred20220820.RData"))

