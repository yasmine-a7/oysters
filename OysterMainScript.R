######################################
##### Oyster Project #####
# Last edit 15/12/22 #
# Authurs: Yasmine, Jonna, and Alice #
######################################


# Load libraries

library(devtools)
library(dplyr)
library(mapplots)
library(shapefiles)
library(tibble)
library(ggplot2)
library(LEA)
library(quantreg)  
library(vcfR)
library(hierfstat)
library(pegas)
library(adegenet)
library(diveRsity)
library(poppr)
library(tidyverse)
library(assignPOP)
library(radiator)


# Load files

vcf <- read.vcfR("oyster_data_6169.vcf")

data <- vcfR2genind(vcf) # transform to genind object


# Investigate 

Samples_names<-rownames(data@tab)
head(Samples_names)

pop<-t(data.frame(strsplit(Samples_names, "_"))) # split names
head(pop)
data@pop<-as.factor(pop)

Pop_ID<-paste(pop[,1]) # vector with pops
head(Pop_ID)
data@pop<-as.factor(Pop_ID)

#Pop_ID<-paste(pop[,1],pop[,2],sep = "_") # vector with pops
#head(Pop_ID)
#data@pop<-as.factor(Pop_ID)


# PCA All loci

x.oyster <- tab(data, freq=TRUE, NA.method="mean")
pca.oyster <- dudi.pca(x.oyster, center=TRUE, scale=FALSE)



s.class(pca.oyster$li, fac=pop(data), col=transp(funky(15),.6), axesel=FALSE, cstar=0, cpoint=3, clabel = 0.5)
add.scatter.eig(pca.oyster$eig[1:50],3,1,2, ratio=.3)


# DAPC All loci

grp<-find.clusters(data, n.pca=50, max.n=18,scale=FALSE) # find clusters

dapc<-dapc(data, pop=grp$grp, n.pca=50)


scatter(dapc, xax=1, yax=2, grp=dapc$grp, 
        col=transp(c("forestgreen","dodgerblue4","deeppink","orange2")),
        pch=19, bg="white",cstar = 0, cellipse = 1,clabel = 1,scree.da=FALSE,
        scree.pca=FALSE)


######################################

# Outlier Selection

stat<-basic.stats(genind2hierfstat(data))$perloc
ggplot(stat,aes(Fst))+ geom_histogram(fill="steelblue4") + 
  geom_vline(xintercept = quantile(stat$Fst,0.95),
             color="forestgreen",size=1.2,linetype ="dashed") +
  geom_vline(xintercept = quantile(stat$Fst,0.99),
             color="red",size=1.2,linetype ="dashed") + theme_bw()

outliers<-stat %>%
  rownames_to_column('locus') %>%  
  filter(Fst>quantile(Fst,0.99))
outliers[,c(1,2,4,8)] # extract outliers


test<-rqss(Fst~qss(Ht,lambda = 0.02),tau=0.99,data =stat)
Fst_predict<-predict.rqss(test,stat, interval = "none", level = 0.99)
data_predict<-cbind(stat,Fst_predict)
outliers2<-data_predict %>% 
  rownames_to_column('locus') %>%  
  filter(Fst>Fst_predict)

ggplot(stat,aes(Ht,Fst))+ geom_point(color=c("grey55"),alpha=0.3) + 
  geom_point(data=outliers2,aes(Ht,Fst),color="red",size=2,alpha=0.5) +
  geom_quantile( method = "rqss", lambda = 0.02,quantiles=c(0.99),color="red")+
  theme_bw() 


# Seperating outlier and neutral loci into new datasets

## heirfstat changed loci names with - to . so...
data[loc=outliers$locus]
locNames(data)
locframe <- data.frame(loc1=locNames(data),loc2=rownames(stat)) # created loci frame with . and - comparison
outlierframe <- subset(locframe, locframe$loc2 %in% outliers$locus) # subset previous frame by outliers
## Data with just the outliers
dataoutlier<- data[loc=outlierframe$loc1] 
## Data with just neutral
neutralframe <- subset(locframe, !locframe$loc2 %in% outliers$locus)
dataneutral<- data[loc=neutralframe$loc1] 

##########

# Outlier PCA 
x.oysteroutlier <- tab(dataoutlier, freq=TRUE, NA.method="mean")
#select 5 # of axes
pca.oysteroutlier <- dudi.pca(x.oysteroutlier, center=TRUE, scale=FALSE)

s.class(pca.oysteroutlier$li, fac=pop(dataoutlier), col=transp(funky(25),.6), axesel=FALSE, cstar=0, cpoint=2, clabel = 0.5)

add.scatter.eig(pca.oysteroutlier$eig[1:50],3,1,2, ratio=.3)

## variance explained by each axes
eig.perc.outlier <- 100*pca.oysteroutlier$eig/sum(pca.oysteroutlier$eig)
head(eig.perc.outlier)

# DAPC outliers

par(mfrow=c(1,1))

set.seed(100)
grp<-find.clusters(dataoutlier, n.pca=50, max.n=18,scale=FALSE) # find clusters - 3

dapc<-dapc(dataoutlier, pop=grp$grp, n.pca=50)
dapc$posterior

set.seed(100)
scatter(dapc, xax=1, yax=2, grp=dapc$grp, 
        col=transp(c("forestgreen","dodgerblue4","deeppink","orange2")),
        pch=19, bg="white",cstar = 0, cellipse = 1,clabel = 1,scree.da=FALSE,
        scree.pca=FALSE)


###############

# Neutral PCA


x.oysterneutral <- tab(dataneutral, freq=TRUE, NA.method="mean")
pca.oysterneutral <- dudi.pca(x.oysterneutral, center=TRUE, scale=FALSE)

s.class(pca.oysterneutral$li, fac=pop(dataneutral), col=transp(funky(25),.6), axesel=FALSE, cstar=0, cpoint=3, clabel = 0.5)
add.scatter.eig(pca.oysterneutral$eig[1:50],3,1,2, ratio=.3)

## variance explained by each axes
eig.perc.neutral <- 100*pca.oysterneutral$eig/sum(pca.oysterneutral$eig)
head(eig.perc.neutral)

# DAPC neutral

set.seed(100)
grp<-find.clusters(dataneutral, n.pca=50, max.n=18,scale=FALSE) # find clusters 4

dapc<-dapc(dataneutral, pop=grp$grp, n.pca=50)
dapc$posterior

scatter(dapc, xax=1, yax=2, grp=dapc$grp, 
        col=transp(c("red", "darkorchid", "cadetblue2", "orange")),
        pch=19, bg="white",cstar = 0, cellipse = 1,clabel = 1,scree.da=FALSE,
        scree.pca=FALSE)

######################

# LEA Outlier

dataoutlier.heirf<- genind2hierfstat(dataoutlier)
write.struct(dataoutlier.heirf, ilab=indNames(dataoutlier), pop=pop(dataoutlier), fname="dataoutlier.str")

struct2geno(input.file="dataoutlier.str", ploidy=2, FORMAT=2, extra.column=2, extra.row=0)

genofile.outlier<- "dataoutlier.str.geno"


project <- snmf(genofile.outlier,K = 2:18,entropy = TRUE,
                repetitions = 2,
                project = "new")
plot(project, col = transp("steelblue4"), pch = 19)

## calculate the ancestry of each individuals for K=3
set.seed(100)
cic.snmf <- snmf(genofile.outlier, K=3, project="new")
## extract the probability matrice for K=3
qmatrix <- Q(cic.snmf, K=3)

## plot this probability matrice 

par(mfrow=c(1,1))
par(mar=c(9, 3, 0.7, 0.7))

barplot(t(qmatrix),
        col=c("forestgreen","dodgerblue4","deeppink"), 
        border=NA, space=0, xlab=NA,
        ylab = "Ancestry")
abline(v=9, lty=2, lwd=2)
abline(v=17, lty=2, lwd=2)
abline(v=29, lty=2, lwd=2)
abline(v=39, lty=2, lwd=2)
abline(v=48, lty=2, lwd=2)
abline(v=56, lty=2, lwd=2)
abline(v=66, lty=2, lwd=2)
abline(v=76, lty=2, lwd=2)
abline(v=85, lty=2, lwd=2)
abline(v=96, lty=2, lwd=2)
abline(v=108, lty=2, lwd=2)
abline(v=120, lty=2, lwd=2)
abline(v=130, lty=2, lwd=2)
abline(v=142, lty=2, lwd=2)
abline(v=148, lty=2, lwd=2)
abline(v=159, lty=2, lwd=2)
abline(v=168, lty=2, lwd=2)
abline(v=179, lty=2, lwd=2)
abline(v=191, lty=2, lwd=2)
abline(v=201, lty=2, lwd=2)
abline(v=211, lty=2, lwd=2)
abline(v=221, lty=2, lwd=2)
text(x=5,y=-0.2,cex=1.2, srt=90, "SeaSalter UK", xpd=NA)
text(x=13,y=-0.2,cex=1.2, srt=90, "Oban UK", xpd=NA)
text(x=23,y=-0.2,cex=1.2, srt=90, "Canada", xpd=NA)
text(x=35,y=-0.2,cex=1.2, srt=90, "Brest FR", xpd=NA)
text(x=45,y=-0.2,cex=1.2, srt=90, "Plymouth UK", xpd=NA)
text(x=53,y=-0.2,cex=1.2, srt=90, "Spain", xpd=NA)
text(x=62,y=-0.2,cex=1.2, srt=90, "Texel NL", xpd=NA)
text(x=71,y=-0.2,cex=1.2, srt=90, "Guernsey UK", xpd=NA)
text(x=80,y=-0.2,cex=1.2, srt=90, "Faro", xpd=NA)
text(x=91,y=-0.2,cex=1.2, srt=90, "Gothenburg", xpd=NA)
text(x=103,y=-0.2,cex=1.2, srt=90, "Ifremer FR", xpd=NA)
text(x=114,y=-0.2,cex=1.2, srt=90, "Japan", xpd=NA)
text(x=125,y=-0.2,cex=1.2, srt=90, "Italy", xpd=NA)
text(x=137,y=-0.2,cex=1.2, srt=90, "Sylt DE", xpd=NA)
text(x=146,y=-0.2,cex=1.2, srt=90, "Bangor UK", xpd=NA)
text(x=155,y=-0.2,cex=1.2, srt=90, "Denmark", xpd=NA)
text(x=166,y=-0.2,cex=1.2, srt=90, "Maldon UK", xpd=NA)
text(x=175,y=-0.2,cex=1.2, srt=90, "Norway", xpd=NA)
text(x=186,y=-0.2,cex=1.2, srt=90, "Oosterschelde NL", xpd=NA)
text(x=210,y=-0.2,cex=1.2, "France hatcheries", xpd=NA)


# LEA Neutral

dataneutral.heirf<- genind2hierfstat(dataneutral)
write.struct(dataneutral.heirf, ilab=indNames(dataneutral), pop=pop(dataneutral), fname="dataneutral.str")

struct2geno(input.file="dataneutral.str", ploidy=2, FORMAT=2, extra.column=2, extra.row=0)

genofile.neutral<- "dataneutral.str.geno"


project.neutral <- snmf(genofile.neutral,K = 2:18,entropy = TRUE,
                repetitions = 2,
                project = "new")
plot(project.neutral, col = transp("steelblue4"), pch = 19)

## calculate the ancestry of each individuals for K=4
set.seed(100)
cic.snmf.neutral <- snmf(genofile.neutral, K=4, project="new")
## extract the probability matrice for K=4
qmatrix.neutral <- Q(cic.snmf.neutral, K=4)

## plot this probability matrice 

par(mfrow=c(1,1))
par(mar=c(9, 3, 0.7, 0.7))

barplot(t(qmatrix.neutral),
        col=c("orange", "red", "cadetblue2", "darkorchid"), 
        border=NA, space=0,
        ylab = "Ancestry")
abline(v=9, lty=2, lwd=2)
abline(v=17, lty=2, lwd=2)
abline(v=29, lty=2, lwd=2)
abline(v=39, lty=2, lwd=2)
abline(v=48, lty=2, lwd=2)
abline(v=56, lty=2, lwd=2)
abline(v=66, lty=2, lwd=2)
abline(v=76, lty=2, lwd=2)
abline(v=85, lty=2, lwd=2)
abline(v=96, lty=2, lwd=2)
abline(v=108, lty=2, lwd=2)
abline(v=120, lty=2, lwd=2)
abline(v=130, lty=2, lwd=2)
abline(v=142, lty=2, lwd=2)
abline(v=148, lty=2, lwd=2)
abline(v=159, lty=2, lwd=2)
abline(v=168, lty=2, lwd=2)
abline(v=179, lty=2, lwd=2)
abline(v=191, lty=2, lwd=2)
abline(v=201, lty=2, lwd=2)
abline(v=211, lty=2, lwd=2)
abline(v=221, lty=2, lwd=2)
text(x=5,y=-0.2,cex=1.2, srt=90, "SeaSalter UK", xpd=NA)
text(x=13,y=-0.2,cex=1.2, srt=90, "Oban UK", xpd=NA)
text(x=23,y=-0.2,cex=1.2, srt=90, "Canada", xpd=NA)
text(x=35,y=-0.2,cex=1.2, srt=90, "Brest FR", xpd=NA)
text(x=45,y=-0.2,cex=1.2, srt=90, "Plymouth UK", xpd=NA)
text(x=53,y=-0.2,cex=1.2, srt=90, "Spain", xpd=NA)
text(x=62,y=-0.2,cex=1.2, srt=90, "Texel NL", xpd=NA)
text(x=71,y=-0.2,cex=1.2, srt=90, "Guernsey UK", xpd=NA)
text(x=80,y=-0.2,cex=1.2, srt=90, "Faro", xpd=NA)
text(x=91,y=-0.2,cex=1.2, srt=90, "Gothenburg", xpd=NA)
text(x=103,y=-0.2,cex=1.2, srt=90, "Ifremer FR", xpd=NA)
text(x=114,y=-0.2,cex=1.2, srt=90, "Japan", xpd=NA)
text(x=125,y=-0.2,cex=1.2, srt=90, "Italy", xpd=NA)
text(x=137,y=-0.2,cex=1.2, srt=90, "Sylt DE", xpd=NA)
text(x=146,y=-0.2,cex=1.2, srt=90, "Bangor UK", xpd=NA)
text(x=155,y=-0.2,cex=1.2, srt=90, "Denmark", xpd=NA)
text(x=166,y=-0.2,cex=1.2, srt=90, "Maldon UK", xpd=NA)
text(x=175,y=-0.2,cex=1.2, srt=90, "Norway", xpd=NA)
text(x=186,y=-0.2,cex=1.2, srt=90, "Oosterschelde NL", xpd=NA)
text(x=210,y=-0.2,cex=1.2, "France hatcheries", xpd=NA)

######################

# MAP!!! 

data2<-data.frame(cbind(Pop_ID,grp$grp))
colnames(data2)<-c("Population","group")
Count<-data.frame(data2 %>% 
                    group_by(Population,group) %>%
                    summarise(n=n()))
head(Count)


latlong<-read.table("oyster_meta_data.csv",sep = ",",header=T,dec=",") # load the latitute-longitude file csv

colnames(latlong)[2]<-"Population"  #change column name to population

Cluster_per_sites<-inner_join(latlong,Count,by="Population")
head(Cluster_per_sites) # merge the two dataset by population

xyz <- make.xyz(Cluster_per_sites$Longitude,
                Cluster_per_sites$Latitude,
                Cluster_per_sites$n,
                Cluster_per_sites$group)




#shapemap

download.file("http://thematicmapping.org/downloads/TM_WORLD_BORDERS_SIMPL-0.3.zip" , destfile="world_shape_file.zip")

shape <- read.shapefile("world_shape_file/TM_WORLD_BORDERS_SIMPL-0.3")

# Map NEUTRAL
##Format to plot the pie charts 
##Map limit
##europe save 700X700
basemap(xlim=c(-5,5),ylim=c(35,65),bg="white")
map<-draw.shape(shape, col="grey85") 
draw.pie(xyz$x, xyz$y, xyz$z, radius=1.5,
         col=c("red", "darkorchid", "cadetblue2", "orange"))
#Japan
basemap(xlim=c(125,130),ylim=c(30,50),bg="white") 
map<-draw.shape(shape, col="grey85")
draw.pie(xyz$x, xyz$y, xyz$z, radius=2,
         col=c("red", "darkorchid", "cadetblue2", "orange"))
#canada
basemap(xlim=c(-140,-125),ylim=c(40,70),bg="white") 
map<-draw.shape(shape, col="grey85")
draw.pie(xyz$x, xyz$y, xyz$z, radius=3,
         col=c("red", "darkorchid", "cadetblue2", "orange"))


## Map OUTLIER

data2<-data.frame(cbind(Pop_ID,grp$grp))
colnames(data2)<-c("Population","group")
Count<-data.frame(data2 %>% 
                    group_by(Population,group) %>%
                    summarise(n=n()))
head(Count)


latlong<-read.table("oyster_meta_data.csv",sep = ",",header=T,dec=",") # load the latitute-longitude file csv

colnames(latlong)[2]<-"Population"  #change column name to population

Cluster_per_sites<-inner_join(latlong,Count,by="Population")
head(Cluster_per_sites) # merge the two dataset by population

xyz <- make.xyz(Cluster_per_sites$Longitude,
                Cluster_per_sites$Latitude,
                Cluster_per_sites$n,
                Cluster_per_sites$group)




#shapemap

download.file("http://thematicmapping.org/downloads/TM_WORLD_BORDERS_SIMPL-0.3.zip" , destfile="world_shape_file.zip")

shape <- read.shapefile("world_shape_file/TM_WORLD_BORDERS_SIMPL-0.3")


##Format to plot the pie charts 
##Map limit
par(mfrow=c(1,1))
#europe save 700X700
basemap(xlim=c(-5,5),ylim=c(35,65),bg="white")
map<-draw.shape(shape, col="grey85") 
draw.pie(xyz$x, xyz$y, xyz$z, radius=1.5,
         col=c("forestgreen","dodgerblue4","deeppink"))
#Japan
basemap(xlim=c(125,130),ylim=c(30,50),bg="white") 
map<-draw.shape(shape, col="grey85")
draw.pie(xyz$x, xyz$y, xyz$z, radius=2,
         col=c("forestgreen","dodgerblue4","deeppink"))
#canada
basemap(xlim=c(-140,-125),ylim=c(40,70),bg="white") 
map<-draw.shape(shape, col="grey85")
draw.pie(xyz$x, xyz$y, xyz$z, radius=3,
         col=c("forestgreen","dodgerblue4","deeppink")) 

##Draw map 
map<-draw.shape(shape, col="grey85")
##Draw the camenbers 
draw.pie(xyz$x, xyz$y, xyz$z, radius=5,
         col=c("forestgreen","dodgerblue4","deeppink"))




############################


## Three plots Neutral

par(mfrow=c(3,2))
# is set the size for the first graphic
par(fig=c(0,4,3,8)/8)
#and plot it
scatter(dapc, xax=1, yax=2, grp=dapc$grp, 
        col=transp(c("forestgreen", "dodgerblue4", "deeppink", "orange2")),
        pch=19, bg="white",cstar = 0, cellipse = 1,clabel = 1,scree.da=FALSE,
        scree.pca=FALSE)
#plot(1,1)
par(new=T)
par(fig=c(4,8,3,8)/8)
#plot(1,1)
##Map limit
basemap(xlim=c(-5,5),ylim=c(35,65),bg="white")
##Draw map 
map<-draw.shape(shape, col="grey85")
##Draw map 
draw.pie(xyz$x, xyz$y, xyz$z, radius=1.5,
         col=c("red", "darkorchid", "cadetblue2", "orange"))
#then i set the size 
par(new=T)
par(fig=c(0,8,0,3)/8)
#plot(1,1)
str<-barplot(t(qmatrix.neutral), col=c("orange2", "forestgreen", "deeppink", "dodgerblue4"), 
             border=NA, space=0, 
             xlab="Individuals", 
             ylab = "Ancestry")


## Tree plots outlier

par(mfrow=c(3,2))
# is set the size for the first graphic
par(fig=c(0,4,3,8)/8)
#and plot it
scatter(dapc, xax=1, yax=2, grp=dapc$grp, 
        col=transp(c("forestgreen","dodgerblue4","deeppink")),
        pch=19, bg="white",cstar = 0, cellipse = 1,clabel = 1,scree.da=FALSE,
        scree.pca=FALSE)
#plot(1,1)
par(new=T)
par(fig=c(4,8,3,8)/8)
#plot(1,1)
##Map limit
basemap(xlim=c(-5,5),ylim=c(35,65),bg="white")
##Draw map 
map<-draw.shape(shape, col="grey85")
##Draw map 
draw.pie(xyz$x, xyz$y, xyz$z, radius = 1,
         col=c("forestgreen","dodgerblue4","deeppink"))
#then i set the size 
par(new=T)
par(fig=c(0,8,0,3)/8)
#plot(1,1)
str<-barplot(t(qmatrix), names.arg=1:232,
             col=c("forestgreen","dodgerblue4","deeppink"), 
             border=NA, space=0, 
             xlab="Individuals", 
             ylab = "Ancestry")


###################################################################
###################################################################

#FST Pairwise OUTLIER  

wc(dataneutral)

basic.stats(dataneutral)

# First we convert the genind file to a format that works with hierfstat
outlier_dat <- genind2hierfstat(dataoutlier)

matFst <- pairwise.WCfst(outlier_dat)

# Here we replace all the NA values with 0 for easier plotting
matFst <- matFst %>% 
  replace(., is.na(.), 0)

#Create label
poplabel <- c("Bangor", "BC", "Brest", "cantambria", "Faro", "FrHat1", "FrHat2", "FrHat3", "FrHat4", "Gothenburg", "Guernsey", "Ifremar", "Italy", "Japan", "Ldk", "Maldon", "Norway", "Oban", "Oosterschelde", "Plymouth", "SeaSalter", "Sylt", "Texel")

#Pop pairwise graph
table.paint(matFst, col.labels=poplabel, clegend = 0.8, csize = 0.8)

#FST Pairwise NEUTRAL

wc(dataneutral)

basic.stats(dataneutral)

# First we convert the genind file to a format that works with hierfstat
neutral_dat <- genind2hierfstat(dataneutral)

matFst <- pairwise.WCfst(neutral_dat)

# Here we replace all the NA values with 0 for easier plotting
matFst <- matFst %>% 
  replace(., is.na(.), 0)

#Create label
poplabel <- c("Bangor", "BC", "Brest", "cantambria", "Faro", "FrHat1", "FrHat2", "FrHat3", "FrHat4", "Gothenburg", "Guernsey", "Ifremar", "Italy", "Japan", "Ldk", "Maldon", "Norway", "Oban", "Oosterschelde", "Plymouth", "SeaSalter", "Sylt", "Texel")

#Pop pairwise graph
table.paint(matFst, col.labels=poplabel, clegend = 0.8, csize = 0.8)

#pop(dataneutral) 

