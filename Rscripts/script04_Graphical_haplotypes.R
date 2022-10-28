#================================
# script04_Graphical_haplotypes
# HH 2022-10-28
#================================

# 2022-10-28
# I410=C72 I412=Z58

# Clean workspace
rm(list=ls())

library(dplyr)
library(purrr)
library(tibble)
library(reshape2)
library(ggplot2)
library(plot.matrix)
options(max.print=1000000)

# set up dirs
# setwd("/home/haixiao/workdir/P4/scripts_R")
DIR_input="../00_Data/Genodat/"
DIR_gendat="../21_GenotypePreparation/"
DIR_output="../22_GraphicalHaplotypes/"
if(!dir.exists(DIR_output)) dir.create(DIR_output)

#--------------------------------------------
# read in genodat
#--------------------------------------------
# read in Genetype_Raw
Genodat=read.csv(paste0(DIR_input,"Genetype_Raw.csv"),row.names = 1,
                 header = T, stringsAsFactors = F,check.names = F)
dim(Genodat) #3072  192
Genodat[1:3,1:4]
colnames(Genodat)
sort(rownames(Genodat)[grepl("SYN",rownames(Genodat))])

# read in map file
map=read.csv(paste0(DIR_input,"Genodat_Physical_Map.csv"),row.names = 1,
             header = T, stringsAsFactors = F,check.names = F)

# modify marker names in both marker matrix and map file
# step 1: extract only letters and numbers from marker name (excluding other symbols)
# step 2: convert all letters to uppercase
map$markername=toupper(gsub("[^A-Za-z0-9]","",map$markername))
rownames(Genodat)=toupper(gsub("[^A-Za-z0-9]","",rownames(Genodat)))
# checking
sapply(list(Genodat,map), dim)
#      [,1]  [,2]
# [1,] 3072 57840
# [2,]  192     3
lapply(list(Genodat,map), function(x) x[1:3,1:3])

# loading the filtered GT data
ls=load(file=paste0(DIR_gendat,"GenotypePrepRes.RData"))
ls #[1] "GenodatQC"         "GenodatAmatCoding" "Genodat012Coding"  "D" "K" "G" "Gsim"
dim(GenodatQC) #1316  187
GenodatQC[1:3,1:4]
colnames(GenodatQC)

# select for markers included in GenodatQC
if(sum(rownames(Genodat) %in% rownames(GenodatQC))==length(rownames(GenodatQC))){
  Genodat=Genodat[rownames(Genodat) %in% rownames(GenodatQC),
                  colnames(Genodat) %in% c(gsub("GID","",colnames(GenodatQC)),"410","412")]
}
dim(Genodat) #1316  189

#---------------------------------------------
# graphical_haplotypes plotting - real data
#---------------------------------------------

dim(Genodat)
Genodat[1:3,1:4]
colnames(Genodat)
Genodat[,c("1","2","3","4","410", "412")]

# step 1: delete monomorphic and missing/-- markers between parental lines
Genodat=Genodat[Genodat$`410`!=Genodat$`412` &
                  !(Genodat[,"410"]=="--" | Genodat[,"412"]=="--"),]
unique(Genodat[,"410"]) #[1] "AA" "GG" "CC" "TT" "AC"
unique(Genodat[,"412"]) #[1] "GG" "AA" "CC" "TT"
dim(Genodat) #1314  189

# Step 2: sort Genedat by chr and position
head(map)
map=map[map$markername %in% rownames(Genodat),]
identical(sort(map$markername),sort(rownames(Genodat)))
dim(map)
str(map)
# sort map
map=map[order(map$chromosome,map$position),]
map[map$chromosome==1,]
# sort Genodat
Genodat=Genodat[map$markername,]
# checking after sorting
identical(map$markername,rownames(Genodat))
head(map)

#---------------------------------------
# Step 3: coding translation
# missing/-- = NA
# heterozygous = 2
# identical to 410(Male/Chang7-2)   = 0
# identical to 412(Female/Zheng58)  = 1

Genodat2=Genodat
unique(c(as.matrix(Genodat2))) #[1] "GG" "--" "CC" "AA" "AG" "TT" "AC" "AT" "CG"
# Genodat2[Genodat2=="AG"|
#             Genodat2=="AC"|
#             Genodat2=="AT"|
#             Genodat2=="CG"]=0.5
# Genodat2[Genodat2=="--"]=NA

# set het as missing to make graphical genotypes cleaner
Genodat2[Genodat2=="AG"|
           Genodat2=="AC"|
           Genodat2=="AT"|
           Genodat2=="CG"|
           Genodat2=="--"]=NA

Genodat2[!is.na(Genodat2) &
           Genodat2==Genodat2$`410`]=0
Genodat2[!is.na(Genodat2) &
           Genodat2==Genodat2$`412`]=1

# checking after coding translation
table(c(as.matrix(Genodat2)),useNA = "ifany")
#        0      1     CC   <NA> 
#   120963 116834      6  15829 

# set CC as NA, because likely genotyping error
# Genodat2[Genodat2=="CC" & !is.na(Genodat2)]=NA
Genodat2[Genodat2!="0" &
           Genodat2!="1" &
           Genodat2!="2" & 
           !is.na(Genodat2)]=NA

# convert chr to integer
Genodat2=apply(Genodat2, 2, as.integer)
rownames(Genodat2)=rownames(Genodat)
dim(Genodat2)

# final checking of translated genotypes
unique(c(as.matrix(Genodat2)))
# [1] 1 NA  0

# add chr and positions to Genodat2
Genodat2[1:3,1:10]
if(identical(map$markername,rownames(Genodat2))){
  Genodat2=data.frame(chr=map$chromosome,
                      position=map$position,
                      markername=map$markername,
                      Genodat2,
                      stringsAsFactors = F,
                      check.names = F)
}
str(Genodat2)
dim(Genodat2)

#----------------------------------------------------------
# for each chr, cluster individuals based on genotypic data
Genodat2chrwise=list()
plotDat=list()
numRecomb=matrix(0, nrow=ncol(Genodat2)-3-2, ncol=10)#-3 for metainfo; -2 for parental lines
i=1
for(i in 1:10){
  # order individuals chr by chr
  tmp=subset(Genodat2,chr==i)
  start=c(tmp$position[1],(tmp$position[-length(tmp$position)]+tmp$position[-1])/2)
  end=c((tmp$position[-length(tmp$position)]+tmp$position[-1])/2,tmp$position[length(tmp$position)])
  tmp[1:3,1:10]
  colnames(tmp)
  # sort DH lines based on similarity to parental lines
  # all genotypes of 410(Male/Chang7-2)   = 0
  # all genotypes of 412(Female/Zheng58)  = 1
  # coding is 0/1, het set as NA
  # so sums of genotypes closer to 0, more similar to 410(Male/Chang7-2)
  # so sums of genotypes closer to 1, more similar to 412(Female/Zheng58)
  tmp_sorted=tmp[,-1*1:3][,order(colSums(tmp[,-1*1:3],na.rm = T))]
  
  # exclude parental lines for haplotype plotting
  tmp_sorted = tmp_sorted[,!(colnames(tmp_sorted) %in% c("410","412"))]
  tmp_sorted[1:3,1:4]
  
  # to use of geom_rect to plot continous y axis, use new GID
  colnames(tmp_sorted)= 1:length(colnames(tmp_sorted))
  
  Genodat2chrwise[[i]]=data.frame(tmp[,1:3],
                                  start=start,
                                  end=end,
                                  tmp_sorted,
                                  stringsAsFactors = F,
                                  check.names = F)
  
  #convert wide to long format for ggplot
  head(Genodat2chrwise[[i]])
  dim(Genodat2chrwise[[i]])
  
  plotDat[[i]] <- melt(Genodat2chrwise[[i]],
                  id.vars = c("chr","position","markername","start","end"),
                  variable.name = "GIDnew")
  plotDat[[i]]$GIDnew=as.integer(plotDat[[i]]$GIDnew)
  str(plotDat[[i]])
  # get the number of recombinations per individuals on each chr
  Genodat2chrwise[[i]][1:3,1:10]
  resdiff=t(apply(Genodat2chrwise[[i]][,-1*1:5], 2, diff))
  numRecomb[,i]=rowSums(resdiff!=0,na.rm = T)
}

#---------------------------------------------------------------
#summary of number of recombinations per chr
rownames(numRecomb)=1:nrow(numRecomb)
colnames(numRecomb)=1:10
head(plotDat[[i]])

# plot the number of recombinations per individuals on each chr
head(numRecomb)
numRecomb_long = melt(numRecomb)
colnames(numRecomb_long)=c("GID","Chr","numRecomb")
head(numRecomb_long)
str(numRecomb_long)
# Number_of_recombinations_on_each_chr
# Number_of_recombinations_on_each_chr_outliershape_NA
# 5.6x3 for merged plot
ggplot(numRecomb_long, aes(x=as.factor(Chr), y=numRecomb)) + 
  xlab("Chromosome")+
  ylab("Number of recombinations")+
  geom_boxplot(outlier.size = 0.5)+ #outlier.shape = NA
  theme_set(theme_bw() + 
              theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.ticks.x=element_blank(),
                axis.ticks.y=element_blank(),
                # axis.ticks.length=unit(-0.15, "cm"),
                axis.text.x=element_text(size = 14),
                axis.text.y=element_text(size = 14),
                axis.title = element_text(size = 14),
                # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
              ))+
  stat_summary(
    fun.data = function(x)
      data.frame(y= Inf, 
                 label = paste(round(median(x), 2)) ),
    geom = "text",
    # aes(group = Tissue),
    hjust = 0.5, vjust = 4,
    position = position_dodge(0.75)
  ) +
  stat_summary(
    fun.data = function(x)
      data.frame(y= Inf, 
                 label = paste(round(min(x), 2)) ),
    geom = "text",
    # aes(group = Tissue),
    hjust = 0.5, vjust = 2,
    position = position_dodge(0.75)
  )+
  stat_summary(
    fun.data = function(x)
      data.frame(y= Inf, 
                 label = paste(round(max(x), 2)) ),
    geom = "text",
    # aes(group = Tissue),
    hjust = 0.5, vjust = 6,
    position = position_dodge(0.75)
  )+
  annotate(geom="text", x = 0.2, y = Inf, label="min",
           hjust = 0, vjust = 2,
           # fontface = "bold",
           color="black")+
  annotate(geom="text", x = 0.2, y = Inf, label="med",
           hjust = 0, vjust = 4,
           # fontface = "bold",
           color="black")+
  annotate(geom="text", x = 0.2, y = Inf, label="max",
           hjust = 0, vjust = 6,
           # fontface = "bold",
           color="black")




apply(numRecomb, 2, summary)
head(numRecomb)

#              Chr1      Chr2      Chr3      Chr4      Chr5      Chr6      Chr7      Chr8      Chr9     Chr10
# Min.     0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
# 1st Qu.  4.000000  2.000000  1.000000  1.000000  4.000000  1.000000  2.000000  0.000000  3.000000  1.000000
# Median   6.000000  4.000000  2.000000  2.000000  6.000000  2.000000  3.000000  2.000000  4.000000  2.000000
# Mean     9.951872  5.513369  4.010695  4.786096  9.657754  2.903743  4.042781  4.903743  5.336898  3.085561
# 3rd Qu. 12.000000  6.000000  4.000000  3.000000 15.000000  3.000000  5.000000  4.000000  7.000000  3.000000
# Max.    68.000000 54.000000 49.000000 77.000000 49.000000 37.000000 38.000000 62.000000 38.000000 33.000000

#-----------------------------------------------------------------------
# graphical genotype - plot a test on ten markers
markers10=head(plotDat[[10]]$markername,10)
tst=subset(plotDat[[10]],markername %in% markers10)
ggplot(tst, 
       aes(position,GIDnew)) +
  # xlab("position")+
  # ylab("position")+
  geom_rect(aes(NULL,NULL,xmin=start,xmax=end,fill = c("red", "blue")[ value + 1 ]),
            ymin=0,ymax=16000,colour="white")+ #colour="white" indicates boundaries of boxes
  scale_fill_identity()+
  scale_x_continuous(breaks = unique(tst$position))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#----------------------------------------------------------------------
# graphical genotpe - plot a entire chromosome
# graphical_genotpe_chr10
dim(plotDat[[10]])
head(plotDat[[10]])

ggplot(plotDat[[10]], aes(position,GIDnew)) +
  xlab("Position")+
  ylab("Lines")+
  # geom_tile(aes(width = box.w)) + 
  geom_rect(aes(NULL,NULL,xmin=start,xmax=end,
                ymin = GIDnew, ymax = GIDnew+1,
                colour="white",
                fill = c("red", "blue")[ value + 1 ])
            )+ #colour="white"
  scale_fill_identity()+
  # scale_x_continuous(breaks = unique(plotDat[[10]]$position))+
  theme_set(theme_bw() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks=element_blank(),
    # axis.ticks.length=unit(-0.25, "cm"),
    axis.text.x=element_blank(),
    axis.text.y=element_blank()
    # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
    )) +
  geom_segment(data=plotDat[[10]][!duplicated(plotDat[[10]]$markername),],
               size=0.5,
               aes(y=0,yend=-5,x=position,xend=position))

#-------------------------------------------------------
# plot graphical genotypes of all ten chrs in one plot
library(foreach)
head(plotDat[[10]])
plotDatAll=foreach(i =1:10, .combine = "rbind") %do% {
  plotDat[[i]]
}
plotDatAll$chr=paste0("Chr",plotDatAll$chr)
head(plotDatAll);tail(plotDatAll)
unique(plotDatAll$value)
# graphical_genotypes_of_all_ten_chrs_QCdat_xaxis_position
ggplot(plotDatAll, aes(position,GIDnew)) +
  xlab("Position")+
  ylab("Lines")+
  facet_wrap(~factor(chr,levels=unique(chr)),scales = "free")+
  geom_rect(aes(NULL,NULL,xmin=start,xmax=end,
                ymin = GIDnew, ymax = GIDnew+1,
                fill = c("red", "blue")[ value + 1 ])
  )+ #colour="white"
  scale_fill_identity()+
  theme_set(theme_bw() + 
              theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.ticks.x=element_blank(),
                axis.ticks.y=element_blank(),
                # axis.ticks.length=unit(-0.15, "cm"),
                axis.text.x=element_blank(),
                axis.text.y=element_blank()
                # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
              ))+
  geom_segment(data=plotDatAll[!duplicated(plotDatAll$markername),],
    aes(y=0,yend=-5,x=position,xend=position),
    size=0.3)

# graphical_genotypes_of_all_ten_chrs_QCdat_sorted_xaxis_markername 
ggplot(plotDatAll, 
       aes(markername,GIDnew, fill = c("red", "blue")[ value + 1 ])) +
  facet_wrap(~factor(chr,levels=unique(chr)),scales = "free")+
  xlab("Markers (sorted)")+
  ylab("Lines")+
  geom_tile() +
  scale_fill_identity() +
  theme_set(theme_bw() + 
              theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.ticks.x=element_blank(),
                axis.ticks.y=element_blank(),
                # axis.ticks.length=unit(-0.15, "cm"),
                axis.text.x=element_blank(),
                axis.text.y=element_blank()
                # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
              ))

