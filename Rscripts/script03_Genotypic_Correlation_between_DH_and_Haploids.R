# =======================================================
# script03_Genotypic_Correlation_between_DHs_and_Haploids
# =======================================================

# Clean workspace
rm(list=ls())

# load libraries
library(rrBLUP)
library(BGLR)
library(lineup)
library(tibble)
library(foreach)
library(sommer)
library(dplyr)

# set up dirs
# setwd("/home/haixiao/workdir/P4/scripts_R")
DIR_gendat="../21_GenotypePreparation/"
DIR_phedat="../06_Phenodat_BLUPs_SingleEnv/"
DIR_output="../32_GenotypicCorr/"

#-------------------------------------------------------------
# read in singleEnv phenotypic and genotypic data
#-------------------------------------------------------------
#------------------------------
# loading data
list.files(DIR_phedat)
ls1=load(file=paste0(DIR_phedat,"Phenodat_BLUEs_BLUPs_SingleEnv.RData"))
ls1 #[1] "blue" "blup" "hsq"
ls2=load(file=paste0(DIR_gendat,"GenotypePrepRes.RData"))
ls2 
# [1] "GenodatQC"         "GenodatAmatCoding" "Genodat012Coding"  "D" "K" "G" "Gsim"
K[1:3,1:4]
head(hsq)

# extract only DH phenotypic data
Phedat=blue
Phedat$GID=paste0("GID",Phedat$GID)
Phedat[1:3,1:4]

# truncate Phedat
Phedat=Phedat[Phedat$GID %in% rownames(K),]

# sort Phedat
rownames(Phedat)=Phedat$GID
identical(rownames(Phedat),rownames(K))#before sorting
Phedat=Phedat[rownames(K),]
identical(rownames(Phedat),rownames(K))#after sorting
identical(rownames(Phedat),colnames(K))
dim(Phedat)
Phedat[1:3,1:4]

# kinship matrix (G matrix)
K[1:10,1:10]
library(matrixcalc)
is.positive.semi.definite(K)

# distance matrix (Euclidean distance based on marker matrix)
D[1:3,1:4]

# makrer matrix
X=GenodatAmatCoding 
X[1:3,1:4] # ind in row, markers in col
X[is.na(X)] = matrix(colMeans(X,na.rm=T),nr = nrow(X),nc = ncol(X),byrow=T)[is.na(X)]
dim(X)

sapply(list(rownames(Phedat),rownames(K),colnames(K)), identical,rownames(X))#TRUE TRUE TRUE


#----------------------------------------------------------------------------------------
# Variance-covariance analysis between Haploid and DH lines
#-----------------------------------------------------------------------------------------

# ----------------------
# data subsetting
# select adjusted data
blue[1:3,1:4]
# dat0=blue
dat0=Phedat
GID=dat0$GID
dat0[1:3,1:4]
dim(dat0)
(EnvALL=unique(sapply(strsplit(as.character(colnames(dat0)[-1]),split='::', fixed=TRUE), `[[`, 3)))
(TissueALL=unique(sapply(strsplit(as.character(colnames(dat0)[-1]),split='::', fixed=TRUE), `[[`, 4)))
TissueALL=c("::4","::12","::99")

GRMName = "K_Amat"; Env = "BJ2014"; Tissue = "::99"; 
# trait = colnames(dat_DH)[1]

EnvSel=c("BJ2014","SJZ2014") #only the two environments, both DH and Hap data has good quality
res=foreach(Env = EnvSel, .combine = "rbind") %do% {
  foreach(Tissue = TissueALL, .combine = "rbind") %do% {
    # ----------------------------------------------------
    # subset data for a specific Env and tissue
    dat=dat0[,grepl(Env, colnames(dat0)) & grepl(Tissue, colnames(dat0))]
    # rownames(dat)=GID
    # if no phenotypic values, skip
    dim(dat)
    head(dat)
    if(dim(dat)[2]==0) return(NULL)
    
    #-----------------------------------------------------
    # subset Haploid and DH lines
    dat_DH=dat[,grepl("DH",colnames(dat))]
    dat_Hap=dat[,grepl("Hap",colnames(dat))]
    
    # imputation - using mean to replace NA in each column
    ##use median to substitute NAs
    dat_DH <- sweep(as.matrix(dat_DH), MARGIN = 2, #notes: mat must be matrix, df does not work
                    STATS = apply(dat_DH, 2, median, na.rm=TRUE),
                    FUN =  function(x,s) ifelse(is.na(x), s, x))
    dat_Hap <- sweep(as.matrix(dat_Hap), MARGIN = 2, #notes: mat must be matrix, df does not work
                     STATS = apply(dat_Hap, 2, median, na.rm=TRUE),
                     FUN =  function(x,s) ifelse(is.na(x), s, x))
    sum(is.na(dat_DH)); sum(is.na(dat_Hap))
    
    # ------------------------------------------------------
    # simplify trait names
    # colnames(dat_DH)=gsub("::DH","",colnames(dat_DH))
    # colnames(dat_Hap)=gsub("::Haploid","",colnames(dat_Hap))
    colnames(dat_DH)=sapply(strsplit(as.character(colnames(dat_DH)),split='::',
                                     fixed=TRUE), `[[`, 1)
    colnames(dat_Hap)=sapply(strsplit(as.character(colnames(dat_Hap)),split='::',
                                      fixed=TRUE), `[[`, 1)
    TraitNames=colnames(dat_DH)
    setdiff(rownames(dat_DH),rownames(K))
    setdiff(rownames(dat_Hap),rownames(K))
    
    # add ploidy to the Pheno dat
    dat_DH=data.frame(Ploidy="DH",dat_DH,stringsAsFactors = F)
    dat_Hap=data.frame(Ploidy="Hap",dat_Hap,stringsAsFactors = F)
    identical(rownames(dat_DH),rownames(K))
    identical(rownames(dat_Hap),rownames(K))
    
    # add a column of GID
    dat_DH=rownames_to_column(.data = dat_DH, var="GID")
    dat_Hap=rownames_to_column(.data = dat_Hap, var="GID")
    identical(dat_DH$GID,rownames(K))
    identical(dat_Hap$GID,rownames(K))
    
    # checking
    head(dat_DH)
    head(dat_Hap)
    dim(dat_DH)
    dim(dat_Hap)
    
    # merge dat_DH and dat_Hap
    dat_DH_Hap_wide=list()
    dat_DH_Hap_long=NULL
    if(identical(colnames(dat_DH),colnames(dat_Hap))){
      dat_DH_Hap_long=rbind(dat_DH, dat_Hap)
      for(trait in TraitNames){
        dat_DH_Hap_wide[trait]=list(left_join(dat_DH[,c("GID",trait)], dat_Hap[,c("GID",trait)],by="GID"))
        colnames(dat_DH_Hap_wide[[trait]])=c("GID","DH","Hap")
        dat_DH_Hap_wide[[trait]]=column_to_rownames(.data = dat_DH_Hap_wide[[trait]],var = "GID")
      }
    }
    
    # (trait = TraitNames[2])
    foreach(trait = TraitNames, .combine = "rbind") %do% {

      # -------------------------------------------------------
      # fitting univariate model with rrBLIP
      head(dat_DH[,trait])
      fm_rrBLUP_DH=mixed.solve(dat_DH[,trait],K=K)
      fm_rrBLUP_DH$Vu #269.3634
      fm_rrBLUP_Hap=mixed.solve(dat_Hap[,trait],K=K)
      fm_rrBLUP_Hap$Vu #93.78934
      Cor_rrBLUP_DH=cor(dat_DH[,trait],fm_rrBLUP_DH$u,use = "c")
      Cor_rrBLUP_Hap=cor(dat_Hap[,trait],fm_rrBLUP_Hap$u,use="c")
      r_P=cor(dat_DH[,trait],dat_Hap[,trait],use = "c")
      
      cor(dat_DH[,trait],fm_rrBLUP_Hap$u,use = "c")
      cor(dat_Hap[,trait],fm_rrBLUP_DH$u,use = "c")
      cor(fm_rrBLUP_Hap$u,fm_rrBLUP_DH$u)
      
      #-----------------------------------------------------
      # multitrait model fitting with BGLR (UN-DIAG model)
      nIter=20000; burnIn=5000
      ETA <- NULL
      # OLD
      # ETA <-list(G_A=list(K=K,model="RKHS"),
      #            G_nonA=list(K=K2,model="RKHS"))#the random effect is UNstructured by default
      ETA <-list(G_A=list(K=K,model="RKHS"))#the random effect is UNstructured by default
      # after meeting with DER on 8/2/22
      BGLR_UN_D <- Multitrait(y=as.matrix(dat_DH_Hap_wide[[trait]]), 
                              ETA=ETA,
                              resCov=list(type="DIAG"),
                              nIter=nIter,burnIn=burnIn)
      
      UHat_UN_D=BGLR_UN_D$ETA[[1]]$u
      head(UHat_UN_D)
      Cor_BGLR_DH=cor(dat_DH_Hap_wide[[trait]]$DH,UHat_UN_D[,1])
      Cor_BGLR_Hap=cor(dat_DH_Hap_wide[[trait]]$Hap,UHat_UN_D[,2])
      head(dat_DH_Hap_wide[[trait]])
      cor(dat_DH_Hap_wide[[trait]])
      
      # compare rrBLUP vs BGLR
      # corbetw2mat(UHat_UN_D,dat_DH_Hap_wide[[trait]])
      # cor(fm_rrBLUP_DH$u,dat_DH[,trait])
      # cor(fm_rrBLUP_Hap$u,dat_Hap[,trait])

      # extract estimated varcov
      BGLR_UN_D$ETA$G_A$Cov$Omega
      BGLR_UN_D$resCov$R
      
      gcov_D_H=BGLR_UN_D$ETA$G_A$Cov$Omega["DH","Hap"]
      gvar_D_D=BGLR_UN_D$ETA$G_A$Cov$Omega["DH","DH"]
      gvar_H_H=BGLR_UN_D$ETA$G_A$Cov$Omega["Hap","Hap"]
      r_g=gcov_D_H/(sqrt(gvar_D_D)*sqrt(gvar_H_H))
      rvar_D_D=BGLR_UN_D$resCov$R[1,1]
      rvar_H_H=BGLR_UN_D$resCov$R[2,2]
      h2_DH_BGLR=gvar_D_D/(gvar_D_D+R_D_D)
      h2_Hap_BGLR=gvar_H_H/(gvar_H_H+R_H_H)
      
      # phenotypic variation with Holland plant breeding book Eq 6.1
      rcov_D_H=0
      r_P2=(gcov_D_H+rcov_D_H)/(sqrt(gvar_D_D + rvar_D_D)*sqrt(gvar_H_H + rvar_H_H))
      r_P=cor(dat_DH[,trait],dat_Hap[,trait],use = "c")
      r_P; r_P2; r_g

      # ------------------------------------------------------
      # collect results
      data.frame(Env=Env,
                 Trait=trait,
                 Tissue=gsub("::","",Tissue),
                 BGLR_gcov_D_H=gcov_D_H,
                 BGLR_gvar_D=gvar_D_D,
                 BGLR_gvar_H=gvar_H_H,
                 BGLR_Rvar_D=rvar_D_D, #residual var, DH
                 BGLR_Rvar_H=rvar_H_H, #residual var, Hap
                 rrBLUP_varU_D=fm_rrBLUP_DH$Vu,
                 rrBLUP_varE_D=fm_rrBLUP_DH$Ve,
                 rrBLUP_varU_H=fm_rrBLUP_Hap$Vu,
                 rrBLUP_varE_H=fm_rrBLUP_Hap$Ve,
                 Cor_BGLR_DH=Cor_BGLR_DH,
                 Cor_rrBLUP_DH=Cor_rrBLUP_DH,
                 Cor_BGLR_Hap=Cor_BGLR_Hap,
                 Cor_rrBLUP_Hap=Cor_rrBLUP_Hap,
                 r_P=cor(dat_DH[,trait],dat_Hap[,trait],use = "c"),
                 r_g=r_g,
                 h2_DH_BGLR=gvar_D_D/(gvar_D_D+R_D_D),
                 h2_DH_rrBLUP=fm_rrBLUP_DH$Vu/(fm_rrBLUP_DH$Vu+fm_rrBLUP_DH$Ve),
                 h2_Hap_BGLR=gvar_H_H/(gvar_H_H+R_H_H),
                 h2_Hap_rrBLUP=fm_rrBLUP_Hap$Vu/(fm_rrBLUP_Hap$Vu+fm_rrBLUP_Hap$Ve),
                 stringsAsFactors = F 
                 ) 
      
    } #foreach-trait
  } #foreach-Tissue
}#foreach-Env
# checking results
dim(res)
write.csv(res,paste0(DIR_output,"VarCovComp_BLUP_BGLR.csv"),row.names = F)

#--------------------------------------------------------------------
# analysis of Phenotypic_and_genotypic_correlaton_between_DH_and_Hap
#--------------------------------------------------------------------
head(res)
range(res$r_g) #

res = add_column(.data = res, 
                 .after = "Tissue_Trait", 
                 Tissue=sapply(strsplit(res$Tissue_Trait,"::"), "[",1))
res = add_column(.data = res, 
                 .after = "Ploidy", 
                 EnvPloidy=paste(res$Env,res$Ploidy,sep = "::"))

res_rP=res[,colnames(res) %in% c("Env","TraitCategory","Tissue_Trait","Tissue",
                                 "r_P")]
res_rG=res[,colnames(res) %in% c("Env","TraitCategory","Tissue_Trait","Tissue",
                                 "r_g")]
colnames(res_rP)[colnames(res_rP)=="r_P"]="r"
colnames(res_rG)[colnames(res_rG)=="r_g"]="r"

res_rP= add_column(res_rP,.after = "Tissue",CorCategory="r_P")
res_rG= add_column(res_rG,.after = "Tissue",CorCategory="r_G")
lapply(list(res_rP,res_rG), head)
res_r = rbind(res_rP, res_rG)

res_r = add_column(.data = res_r, 
                   .after = "CorCategory", 
                   Env_CorCategory=paste(res_r$Env,res_r$CorCategory,sep = "::"))


