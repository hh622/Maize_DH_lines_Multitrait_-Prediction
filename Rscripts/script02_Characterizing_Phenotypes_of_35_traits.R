# ===================================================
# script02_Characterizing_Phenotypes_of_35_traits
# HH 2022-10-28
# ===================================================

# 2022-10-28
# characterize phenotypes of DH and haploid populations
#  - I computed basic statistics of DH and Haploid populations across sites(supplemental)
#  - I also plot phenotypic distribution across sites (Figure 1)

# Clean workspace
rm(list=ls())

library(tibble)
library(foreach)

# set dir
DIR_input="../02_Phenodat_PlotMean/"
DIR_output="../11_PhenotypicAnalysis/"
if(!dir.exists(DIR_output)) dir.create(DIR_output)

# loading data
ls=load(paste0(DIR_input,"Phedat_Plotmeans_lst.RData"))
ls # "phedat_internode_lst"  "phedat_wholeplant_lst"
phedat_internode_lst[[1]]$InternodeID
names(phedat_internode_lst)
names(phedat_wholeplant_lst)
head(phedat_internode_lst[[1]])
head(phedat_wholeplant_lst[[1]])

head(phedat_internode_lst[[2]])
colnames(phedat_internode_lst[[1]])

#-----------------------------------------------
# merge internode data with whole plant data
#-----------------------------------------------

i=1

# for internode data
for(i in 1:length(phedat_internode_lst)){
  phedat_internode_lst[[i]]=add_column(.data = phedat_internode_lst[[i]],
                                       .after = "Ploidy",
                                       Trait=colnames(phedat_internode_lst[[i]])[length(colnames(phedat_internode_lst[[i]]))]
                                       )
  head(phedat_internode_lst[[i]])
  colnames(phedat_internode_lst[[i]])[length(colnames(phedat_internode_lst[[i]]))]="Value"
}

# for wholeplant data
for(i in 1:length(phedat_wholeplant_lst)){
  phedat_wholeplant_lst[[i]]=add_column(.data = phedat_wholeplant_lst[[i]],
                                       .after = "Ploidy",
                                       Trait=colnames(phedat_wholeplant_lst[[i]])[length(colnames(phedat_wholeplant_lst[[i]]))]
  )
  head(phedat_wholeplant_lst[[i]])
  colnames(phedat_wholeplant_lst[[i]])[length(colnames(phedat_wholeplant_lst[[i]]))]="Value"
}


head(phedat_internode_lst[[1]])
head(phedat_wholeplant_lst[[1]])

identical(colnames(phedat_internode_lst[[1]]),colnames(phedat_wholeplant_lst[[1]]))

# merge phedat lists across traits
df_phedat_internode=foreach(i = 1:length(phedat_internode_lst),.combine = "rbind") %do% {
  phedat_internode_lst[[i]]
  # str(phedat_internode_lst[[i]])
}
df_phedat_wholeplant=foreach(i = 1:length(phedat_wholeplant_lst), .combine = "rbind") %do% {
  phedat_wholeplant_lst[[i]]
}

# merge phedat_wholeplant_lst across traits
sapply(list(df_phedat_internode,df_phedat_wholeplant), dim)
lapply(list(df_phedat_internode,df_phedat_wholeplant), head)
identical(colnames(df_phedat_internode),colnames(df_phedat_wholeplant))
df_phedat = rbind(df_phedat_internode,df_phedat_wholeplant)
head(df_phedat)

#-----------------------------------------------
# grouped histogram DH vs Haploids
#-----------------------------------------------

# Fig.cap=Phedat_ggplot_histogram_Internode_OutlierRemoval_AllEnv_DH_vs_Hap
head(df_phedat)
length(unique(df_phedat$Trait))
df_phedat=add_column(.data = df_phedat, 
           .after = "Trait",
           Internode_Trait=paste(df_phedat$InternodeID,df_phedat$Trait,sep = "::")
           )
length(unique(df_phedat$Internode_Trait))

# add trait categories
df_phedat=add_column(.data = df_phedat, .after = "Internode_Trait",TraitCategory="Stalk quality")
df_phedat$TraitCategory[df_phedat$Trait %in% c("bRPR","InternodeLength","InternodeDiameter",
                                               "InternodeCounts",
                                   "FreshWeight","mRPR","DryWeight","PlantHeight",
                                   "EarHeight","LeafLength","LeafWidth","LeafAngle")]="agronomic"
unique(df_phedat$Trait[df_phedat$TraitCategory=="Stalk quality"])


# sort df_phedat
df_phedat=df_phedat[order(df_phedat$TraitCategory,df_phedat$Trait),]

# select for internode 4,12,99 
df_phedat=df_phedat[df_phedat$InternodeID %in% c(4,12,99),]
head(df_phedat)

# further modification
df_phedat=df_phedat[df_phedat$Trait!="bRPR",]
df_phedat=df_phedat[df_phedat$Trait!="InternodeCounts",]
unique(df_phedat$Trait)
head(df_phedat)

# modify traitname before plotting
df_phedat$Internode_Trait=gsub("4::","FI::",df_phedat$Internode_Trait)
df_phedat$Internode_Trait=gsub("12::","EI::",df_phedat$Internode_Trait)
df_phedat$Internode_Trait=gsub("99::","WP::",df_phedat$Internode_Trait)
df_phedat$Internode_Trait=gsub("mRPR","RPR",df_phedat$Internode_Trait)

# sort phenotypic data: first by TraitCategory, then by Internode_Trait
df_phedat = df_phedat[order(df_phedat$TraitCategory, df_phedat$Internode_Trait),]

# ----------------------------------
# basic statistics
# ----------------------------------
head(df_phedat)
apply(df_phedat[,-ncol(df_phedat)], 2, unique)
trait = unique(df_phedat$Internode_Trait)[1]
# a function of sigLevel
sigLevel = function(pval){
  if(pval>0.05) sigLevel="NS"
  if(pval<=0.05) sigLevel="*"
  if(pval<0.01) sigLevel="**"
  if(pval<0.001) sigLevel="***"
  return(sigLevel)
}
res=foreach(trait = unique(df_phedat$Internode_Trait),.combine = "rbind") %do% {
  dat = df_phedat[df_phedat$Internode_Trait==trait,]
  head(dat)
  tail(dat)
  unique(dat$InternodeID)
  
  # distribution plot
  # ggplot(dat[dat$InternodeID == 4,], aes(x=Value,fill=Ploidy)) + 
  #   geom_histogram(bins = 30) + 
  #   facet_wrap(~factor(Internode_Trait,
  #                      levels=unique(dat[dat$InternodeID == 4,]$Internode_Trait)),
  #              scales = 'free')
  # subset for DH and Haploid
  dat_DH=subset(dat,Ploidy=="DH")
  dat_Hap=subset(dat,Ploidy=="Haploid")
  # calculate mean and sd
  mean_DH=mean(dat_DH$Value)
  mean_Hap=mean(dat_Hap$Value)
  sd_DH=sd(dat_DH$Value)
  sd_Hap=sd(dat_Hap$Value)
  
  # prep DH and Hap data 
  head(dat_DH)
  head(dat_Hap)
  dat_DH= add_column(.data = dat_DH,
                     GID=as.integer(gsub("[^0-9]","",dat_DH$Lname)),
                     .before = "Lname")
  dat_DH= add_column(.data = dat_DH,
                     newID=paste(dat_DH$Env,dat_DH$Block,dat_DH$GID,sep="_"),
                     .before = "Value")
  
  dat_Hap= add_column(.data = dat_Hap,
                     GID=as.integer(gsub("[^0-9]","",dat_Hap$Lname)),
                     .before = "Lname")
  dat_Hap= add_column(.data = dat_Hap,
                     newID=paste(dat_Hap$Env,dat_Hap$Block,dat_Hap$GID,sep="_"),
                     .before = "Value")
  colnames(dat_DH)[colnames(dat_DH)=="Value"]="Value_DH"
  colnames(dat_Hap)[colnames(dat_Hap)=="Value"]="Value_Hap"
  
  length(unique(dat_DH$newID))==dim(dat_DH)[1]
  length(unique(dat_Hap$newID))==dim(dat_Hap)[1]
  library(dplyr)
  dat2=right_join(dat_DH[,c("newID","Value_DH")],
                  dat_Hap[,c("newID","Value_Hap")],
                  by="newID")
  head(dat2)
  
  # t test
  pval=t.test(dat2$Value_DH,dat2$Value_Hap, paired=TRUE)$p.value
  
  data.frame(Internode_Trait=trait,
             TraitCategory=unique(dat$TraitCategory),
             mean_DH=round(mean_DH,2),
             mean_Hap=round(mean_Hap,2),
             mean_diff=round(mean_DH-mean_Hap,2),
             sd_DH=round(sd_DH,2),
             sd_Hap=round(sd_Hap,2),
             CV_DH=round(sd_DH/mean_DH,2),
             CV_Hap=round(sd_Hap/mean_Hap,2),
             pval=signif(pval, digits=3),
             sigLevel=sigLevel(pval)
             )
}

# output basic statistics for table S1
head(res)
output_tableS1=F
if(output_tableS1){
  write.csv(res,paste(DIR_output,"Basic_stat_DH_vs_Haplid20220923.csv"),row.names = F)
}

#-------------------------------------
# compare sd_DH vs sd_Hap
#-------------------------------------
ggplot(data=res, aes(res$sd_DH, res$sd_Hap,color=TraitCategory)) +
  geom_point()+
  # scale_color_manual(values=c("#00BFC4", "#F8766D"))
  scale_color_manual(values=c("black", "red"))+
  geom_text(aes(label=Internode_Trait))+
  geom_abline(slope=1, intercept = 0)

unique(res$TraitCategory) #[1] "agronomic"     "Stalk quality"
ggplot(data=res[res$TraitCategory=="agronomic",], aes(sd_DH, sd_Hap)) +
  geom_point()+
  # scale_color_manual(values=c("#00BFC4", "#F8766D"))
  scale_color_manual(values=c("black", "red"))+
  geom_text(aes(label=Internode_Trait))+
  geom_abline(slope=1, intercept = 0)

ggplot(data=res[res$TraitCategory=="Stalk quality",], aes(sd_DH, sd_Hap)) +
  geom_point()+
  # scale_color_manual(values=c("#00BFC4", "#F8766D"))
  scale_color_manual(values=c("black", "red"))+
  geom_text(aes(label=Internode_Trait))+
  geom_abline(slope=1, intercept = 0)

#-------------------------------------
# compare CV_DH vs CV_Hap
#-------------------------------------
head(res)
res_agr=res[res$TraitCategory=="agronomic",]
res_stalkqual=res[res$TraitCategory=="Stalk quality",]
res_cv= data.frame(
  TraitCategory=unique(res$TraitCategory),
  pval=c(
    t.test(res_agr$CV_DH,res_agr$CV_Hap, paired=TRUE)$p.value,
    t.test(res_stalkqual$CV_DH,res_stalkqual$CV_Hap, paired=TRUE)$p.value
  ))
res_cv$pval=signif(res_cv$pval,digits = 2)
head(res_cv)

library(ggrepel)
# Phenotypic_dispersion_all_InternodeTraits_CVDH_vs_CVHap
ggplot(data=res, aes(CV_DH, CV_Hap,color=TraitCategory)) +
  geom_point()+
  scale_color_manual(values=c("black", "blue"))+
  facet_wrap(~TraitCategory, ncol=1)+
  # geom_text(aes(label=Internode_Trait),show.legend = FALSE,position=position_jitter(width=0.1,height=0.1))+
  geom_text_repel(aes(label=Internode_Trait),
                  show.legend = FALSE,
                  force=1, point.padding=unit(0.1,'lines'),
                  direction = 'both',
                  nudge_x = 0.1,
                  nudge_y = 0.0, #adjust the starting y position of the text label, MUST 0 here!!!
                  segment.size=0.1)+
  geom_text(data = res_cv, aes(x = 0.05, y = 1, 
            label = paste0("p=",pval)),
            size=5,
            hjust=0, vjust=1)+
  geom_abline(slope=1, intercept = 0)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none",
        strip.text.x = element_text(size = 12))
  # theme( legend.position=c(0.95, 0.05),
  #        legend.justification = c(1, 0))

# ----------------------------------------
# res_agr
unique(res$TraitCategory) #[1] "agronomic"     "Stalk quality"
ggplot(data=res_agr, 
       aes(CV_DH, CV_Hap)) +
  geom_point()+
  # scale_color_manual(values=c("#00BFC4", "#F8766D"))
  scale_color_manual(values=c("black", "red"))+
  geom_text_repel(aes(label=Internode_Trait))+
  geom_abline(slope=1, intercept = 0)
t.test(res_agr$CV_DH,res_agr$CV_Hap, paired=TRUE)$p.value #0.01253374
sum(res_agr$CV_DH<res_agr$CV_Hap)/length(res_agr$CV_DH) #13/17=0.7647

# -----------------------------------------
# res_stalkqual
dim(res_stalkqual)
ggplot(data=res_stalkqual, aes(CV_DH, CV_Hap)) +
  geom_point()+
  # scale_color_manual(values=c("#00BFC4", "#F8766D"))
  scale_color_manual(values=c("black", "red"))+
  geom_text(aes(label=Internode_Trait))+
  geom_abline(slope=1, intercept = 0)
t.test(res_stalkqual$CV_DH,res_stalkqual$CV_Hap, paired=TRUE)$p.value #0.01253374
sum(res_stalkqual$CV_DH>res_stalkqual$CV_Hap) #11
length(res_stalkqual$CV_DH) #18

sapply(list(res_agr, res_stalkqual), nrow)#17 18

#-------------------------------
# ggplot
# ------------------------------
# Phenotypic_distribution_all_InternodeTrait
head(df_phedat)
df_phedat=df_phedat[order(df_phedat$TraitCategory,df_phedat$Internode_Trait),]

ggplot(df_phedat[df_phedat$InternodeID %in% c(4,12,99),], aes(x=Value)) + 
  # geom_boxplot(aes(fill=Ploidy))+
  geom_histogram(bins = 30,aes(fill=Ploidy)) +
  xlab("")+
  ylab("Counts")+
  facet_wrap(~ factor(Internode_Trait, levels = unique(df_phedat$Internode_Trait)), scales = 'free')+
  theme_set(theme_bw() + 
              theme(legend.position=c(1,0),
                    legend.justification = c("right", "bottom"),
                    legend.title = element_blank(),
                    axis.title=element_text(size=10,face="bold"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank()
                    ))+
  geom_text(data = res, aes(x = Inf, y = Inf, label = sigLevel), size=5,hjust=1.2, vjust=1.5)

