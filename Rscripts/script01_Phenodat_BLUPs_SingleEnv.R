#==================================================================
# script01_Phenodat_BLUPs_SingleEnv
# HH 2022-10-28
#==================================================================

# 2022-10-28
# compute BLUPs in each Env based on the plot mean.

# Clean workspace
rm(list=ls())

library(tibble)
library(lme4)
library(sommer)
library(dplyr)
library(lsmeans)
packageVersion("lme4") #1.1.28, the latest version
packageVersion("sommer") #4.1.5

# set dir
DIR_input="../02_Phenodat_PlotMean/"
DIR_output="../06_Phenodat_BLUPs_SingleEnv/"
if(!dir.exists(DIR_output)) dir.create(DIR_output)

# loading data
ls=load(paste0(DIR_input,"Phedat_Plotmeans_lst.RData"))
ls # "phedat_internode_lst"  "phedat_wholeplant_lst"
phedat_internode_lst[[1]]$InternodeID
names(phedat_internode_lst)

#---------------------------------------------------------------------
# SECTION 01: calculate BLUPs and hsq for each Env - Internode traits
#---------------------------------------------------------------------

# -------------------------------------
# for testing my script
phedat_lst=phedat_internode_lst
# phedat_lst=phedat_wholeplant_lst
names(phedat_lst)
i=1; ploidy="DH"; env="BJ_2014";internodeID=4 #testing data
intersect(names(phedat_internode_lst),names(phedat_wholeplant_lst))

# Empty variables going to be used 
hsq=NULL; blup=NULL; blue=NULL; use.sommer=F; use.lme4=T
blup_DH=NULL; blue_DH=NULL; blup_Hap=NULL; blue_Hap=NULL

# Run for loops for calculating BLUPS, BLURs and hsq within each Env
for(Dataset in c("internode","wholeplant")){
  # ------------------------------------
  # Selecting Data
  if(Dataset=="internode"){
    phedat_lst=phedat_internode_lst
  }
  if(Dataset=="wholeplant"){
    phedat_lst=phedat_wholeplant_lst
  }
  
  for(i in 1:length(phedat_lst)){ #3nd layer: loop through all Traits
    # for(i in 1:2){ #for testing my script
    #----------------------------------------------------------------------------------
    # take phenotypic data that contains all five Envs + both Ploidy levels for a Trait
    dfALL=phedat_lst[[i]]
    dfALL$Block=as.character(dfALL$Block)
    for(ploidy in sort(unique(dfALL$Ploidy))) { #3nd layer: loop through ploidy levels
      
      for(env in sort(unique(dfALL$Env))) { #4rd layer: loop through each Env within a ploidy level
        #take phedat that contains all tissues/internodes from a specific Env under a specific ploidy level
        dfEnv=dfALL[dfALL$Env==env & dfALL$Ploidy==ploidy,]
        # head(dfEnv)
        for(internodeID in unique(dfEnv$InternodeID)){ #loop through all internodes
          
          print(paste("########Ploidy=",ploidy,"Env=",env,"InternodeID=",internodeID, "PhenoID=",i))
          #-----------------------------------------------------------
          # Prepare data frame for model fitting
          # taking phedat from the min Unit (internode:Env:Ploidy)
          if(sum(dfEnv$InternodeID==internodeID)==0){
            next #skip if empty phenotype values there
          } 
          dfUnit=dfEnv[dfEnv$InternodeID==internodeID,]
          # dim(dfUnit) #359   8
          # str(dfUnit)
          
          # skip trait with high missing rate or only has one Rep
          if( length(unique(dfUnit$Lname)) < 30 | 
              length(levels(factor(dfUnit$Block))) < 2 ){
            next #to skip the current iteration of a loop without terminating it
          }
          
          
          # Changing all columns of metainfo to factor
          index <- 1:5
          dfUnit[ , index] <- lapply(dfUnit[ , index], as.factor)
          str(dfUnit)
          apply(dfUnit[,1:5], 2, function(x) length(levels(factor(x))))
          head(dfUnit)
          
          # ---------------------------------------------------------
          # LMM fitting
          if(use.sommer){
            fm_sommer = mmer(InternodeLength~ Block, #fixed
                             random = ~ Lname,
                             rcov = ~ units,
                             data = dfUnit)
          }
          # fitting with lme4 - lmer
          if(use.lme4){
            fm_lme4 <- lmer(dfUnit[,ncol(dfUnit)] ~ (1|Lname) + Block,data = dfUnit) 
          }
          
          # -------------------------------------------------------
          # Extract variance components and calculate heritability
          if(use.sommer){
            varcomp <- summary(fm_sommer)$varcomp
            rownames(varcomp)=sapply(strsplit(as.character(rownames(varcomp)),split='.', fixed=TRUE), `[[`, 1)
            g.var <- varcomp["Lname","VarComp"]
            e.var <- varcomp["units","VarComp"]
            NumRep=length(levels(dfUnit$Block))
            p.var <- g.var + e.var/(NumRep)
            hsq.i <- g.var/p.var
          }
          
          if(use.lme4){
            print(VarCorr(fm_lme4),comp=c("Variance"))
            varcomp=data.frame(VarCorr(fm_lme4),comp=c("Variance"))
            g.var <- varcomp[varcomp$grp=="Lname","vcov"]
            e.var <- varcomp[varcomp$grp=="Residual","vcov"]
            NumRep=length(levels(dfUnit$Block))
            p.var <- g.var + e.var/(NumRep)
            hsq.i <- g.var/p.var
            hsq.i <- data.frame(Ploidy=ploidy,
                                Env=env,
                                # TraitName=paste(names(phedat_lst)[i],internodeID,sep = "::"),
                                TraitName=paste(names(phedat_lst)[i],ploidy,gsub("_","",env),internodeID,sep = "::"),
                                numLines=length(unique(dfUnit$Lname)),
                                hsq=hsq.i,
                                stringsAsFactors = F
                                )
          }
          
          # ------------------------------------
          # Extract Line BLUPs from LMM fitting
          blup.i=NULL
          if(use.sommer){
            blup.i <- data.frame(randef(fm_sommer)$'Lname', stringsAsFactors = F)
            blup.i = data.frame(Ploidy=ploidy,
                                Env=env,
                                Lname=gsub("Lname","",rownames(blup.i)),
                                GID=gsub("[^0-9]","",rownames(blup.i)),
                                y=blup.i$y)
            # colnames(blup.i)[colnames(blup.i)=="y"]=colnames(phedat)[i]
          }
          if(use.lme4){
            blup.i=data.frame(ranef(fm_lme4)$Lname)
            colnames(blup.i)=paste(names(phedat_lst)[i],ploidy,gsub("_","",env),internodeID,sep = "::")
            blup.i = data.frame(#Ploidy=ploidy,
              # Env=env,
              #Lname=rownames(blup.i),
              GID=as.integer(gsub("[^0-9]","",rownames(blup.i))),
              blup.i,
              check.names = F,
              stringsAsFactors = F)
            rownames(blup.i)=NULL
          }
          # head(blup.i)
          # str(blup.i)
          # --------------------------------
          # Fitting LM with built-in stats
          fm_lm <- lm(dfUnit[,ncol(dfUnit)] ~ Block + Lname, data=dfUnit)
          summary(fm_lm)
          
          # -----------------------------------
          # Extract Line BLUEs from LM fitting
          blue.i=NULL
          # extract Lname and lsmean
          blue.i <- data.frame(Lname=as.vector(summary(lsmeans(fm_lm,'Lname'))$Lname),
                               y=summary(lsmeans(fm_lm,'Lname'))$lsmean,
                               stringsAsFactors = F)
          # add Ploidy and GID
          blue.i = data.frame(#Ploidy=ploidy,
            # Env=env,
            #Lname=blue.i$Lname,
            GID=as.integer(gsub("[^0-9]","",blue.i$Lname)),
            y=blue.i$y,
            stringsAsFactors = F)
          colnames(blue.i)[colnames(blue.i)=="y"]=paste(names(phedat_lst)[i],
                                                        ploidy,
                                                        gsub("_","",env),
                                                        internodeID,sep = "::")
          # head(blue.i)
          # str(blue.i)
          
          #--------------------------------------
          cmp.blue_blup=F
          if(cmp.blue_blup){
            lapply(list(blue.i,blup.i), head)
            sapply(list(blue.i,blup.i), dim)
            identical(blue.i$Lname,blup.i$Lname)
            plot(blue.i[,ncol(blue.i)],blup.i[,ncol(blup.i)]); abline(0,1)
          }
          
          #-------------------------
          #accumulate hsq
          head(hsq.i)
          if(is.null(hsq)){
            hsq <- hsq.i
          }else{
            hsq <- rbind(hsq,hsq.i)
          }
          
          #------------------------
          #accumulate BLUPs
          if(is.null(blup)){
            blup <- blup.i
            blue <- blue.i
          }else{
            blup <- full_join(blup, blup.i, by = c("GID"))
            blue <- full_join(blue, blue.i, by = c("GID"))
          }
        }#end of 5th layer: loop for internodeID
      } #end of 4th layer: loop for env
    } #end of 3rd layer: loop for ploidy
  } #end of 2nd layer: loop for phedatypes
} #end of 1st layer: loop for datasets

# --------------------------------------------------
# checking results
plot(hsq$numLines,hsq$hsq) #heritability not correlate with hsq$numLines
sapply(list(blue=blue,blup=blue,hsq=hsq), dim)
#      blue blup hsq
# [1,]  193  193 750
# [2,]  751  751   5 
# => there are 750 traits in total

lapply(list(blue,blup,hsq), function(x) x[1:10,1:5])
sapply(list(blue=blue,blup=blue,hsq=hsq), str)

sapply(list(blue=colnames(blue),
            blup=colnames(blup),
            hsq=c("GID",hsq$TraitName)), 
       identical,colnames(blue))
# blue blup  hsq 
# TRUE TRUE TRUE
# => trait names are identical among blue, blup and hsq

identical(blue$GID,blup$GID)
# GID are identical between blue and blup

library(lineup)
hist(corbetw2mat(blue,blup))
summary(corbetw2mat(blue,blup))
sum(is.na(corbetw2mat(blue,blup))) #58
# note: for 58 traits, which means trait filtering is needed.
traitnamesNA=names(corbetw2mat(blue,blup))[is.na(corbetw2mat(blue,blup))]
blue[,traitnamesNA][,1:3]
data.frame(blue[,traitnamesNA[2]],blup[,traitnamesNA[2]],check.names = F)

# -----------------------------------------------------
# save results (blue, blup, hsq)
#------------------------------------------------------
save(blue,blup,hsq, file=paste0(DIR_output,"Phenodat_BLUEs_BLUPs_SingleEnv.RData"))

