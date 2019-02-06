setwd("~/projects/Guittar/Reanalysis_Anders_2016/Clustering_all_ex_mock_20161208")
require(vegan)
require(nlme)
source('~/kode/R/dropRareTaxa.R')
source('~/kode/R/diversity.r')
source('~/kode/R/correlationTests.r')
source('~/kode/R/taxaplot.R')
source('~/kode/R/taxacorr.R')

md = read.csv("../metadata.csv",header=T,row.names=1)
mdt = md[md$project=="Turfs",]

mdt.c = mdt[mdt$turf.treat.group=="Control",]
mdt.t = mdt[mdt$turf.treat.group=="Transplant",]

#OTUs.ds.rl.c = OTUS.ds.rl[mdt$turf.treat.group=="Control",]
ra_all = read.table("CREST_Results/Relative_Abundance_noEuk.csv",sep=",",
                    header=T,row.names=3)
taxa_turf = ra_all[,c(T,T,md$project=="Turfs")]


# Make precipitation and temp level take intermediate levels for transplants
mdt$prec.level.comp = mdt$precip.level - as.numeric(mdt$turf.treat=="Wetter")/2
mdt$temp.level.comp = mdt$temp.level - as.numeric(mdt$turf.treat=="Warmer")/2

minAbundance = 1/3493 #Representing one read out of all classified in sample with fewest reads

taxa.ab = taxa_turf[rowMeans(taxa_turf[-c(1:2)])>minAbundance,] 
#mdt.c=mdt[mdt$turf.treat=="Control",]

# --- Taxa barplots -----------

mdt.warmer=mdt[mdt$turf.treat=="Warmer",]
mdt.wetter=mdt[mdt$turf.treat=="Wetter",]

mdt$group=rep("",53)
mdt[mdt$turf.treat=="Control",]$group = paste("T",mdt.c$temp.level," P",mdt.c$precip.level," Ctrl",sep="")
mdt[mdt$turf.treat=="Warmer",]$group = paste("T",mdt.warmer$origSite.tempLevel,"->T",mdt.warmer$temp.level," P",mdt.warmer$precip.level," Trans",sep="")
mdt[mdt$turf.treat=="Wetter",]$group = paste("T",mdt.wetter$temp.level," P",mdt.wetter$origSite.precipLevel,"->P",mdt.wetter$precip.level," Trans",sep="")
mdt$name = paste(mdt$group,mdt$rep)

for (level in c("kingdom","class","order","family", "genus")){
  levelTaxa = as.data.frame(t(taxa_turf[taxa_turf$Level==level,-c(1:2)]))
  print(paste(mean(rowSums(levelTaxa))*100,"% classified at ",level))
  row.names(levelTaxa)=mdt$name
  grouping_info<-data.frame(row.names=mdt$name,mdt$group)
  pdf(paste("../img/",level,"_barchart.pdf",sep=""),height=6,width=24)
  taxaplot(30,grouping_info,levelTaxa)
  dev.off()
}
# [1] "98.8167075471698 % classified at  domain"
# [1] "97.7753188679245 % classified at  kingdom"
# [1] "96.2931094339623 % classified at  class"
# [1] "86.9246113207547 % classified at  order"
# [1] "66.1246377358491 % classified at  family"
  
# --- Check underlying distro in taxa RA ---

orders = taxa_turf[taxa_turf$Level=="order",-c(1:2)]
orderMeans = rowMeans(orders)
hist(log(orderMeans),breaks=20,col="grey")
# So much nicer

# --- Multiple LMs for ALL soils (insetad of consensus) ------

lmBoth = function(ra,file_sign,obs){
  
  lm2 = lm(ra ~prec.level.comp+temp.level.comp, data=mdt)
  pp = summary(lm2)$coefficients[2,4]
  pt = summary(lm2)$coefficients[3,4]
  tp = summary(lm2)$coefficients[2,3]
  tt = summary(lm2)$coefficients[3,3]
  
  if (pp<0.05/obs && pt<0.05/obs) sig="TP"
  else if (pp<0.05/obs) sig="P"
  else if (pt<0.05/obs) sig="T"
  else sig="None"
  
  
  if (pp<0.05/obs){
    print(summary(lm2))
    plot(ra~prec.level.comp,data=mdt,main=o)
    abline(lm2,col="grey")
    line <- readline()
  }
  
  if (pt<0.05/obs){
    print(summary(lm2))
    plot(ra~temp.level.comp,data=mdt,main=o)
    abline(lm2,col="grey")
    line <- readline()
  }
  
  write(c(o,as.character(taxa.ab[o,]$Level),
            as.character(taxa.ab[o,]$Taxonpath),mean(ra),
            sig,pp,pp*obs,pt,pt*obs,tp,tt),
          file=file_sign,ncolumns=12,sep=",",append=T)
}

# Text output (all taxa regardless of significance)
sfile="../lm_taxa_corr_w_precip_or_temp_level_trans.csv"

write(c("Taxon","Level","Classification","Avg.abundance","Significance",
        "p.prec","p.prec.bonf","p.temp","p.temp.bonf","t.prec","t.temp"),
      file=sfile,ncolumns=12,sep=",")

for (level in c("domain","kingdom","class","order","family","genus")){
  taxa = taxa.ab[taxa.ab$Level==level,]
  obs = dim(taxa)[1]
  
  for (o in row.names(taxa)){
    ra = as.numeric(taxa[o,-c(1,2)])
    lmBoth(ra,file_sign=sfile,obs=obs)
  }
}

# ------------ 2) ANOVA model including transplants -----------------

# ---- Temp ----

mdt.temp = mdt[mdt$turf.treat!="Wetter",]
mdt.temp$temp = paste(mdt.temp$origSite.tempLevel,mdt.temp$turf.treat.group)

taxa_temp = taxa_turf[rowMeans(taxa_turf[-c(1:2)])>minAbundance,c(T,T,mdt$turf.treat!="Wetter")]

for (level in c("domain","kingdom","class","order","family","genus")){
  levelTaxa = as.data.frame(t(taxa_temp[taxa_temp$Level==level,-c(1:2)]))
  printANOVA1Factor(mdt.temp$temp,levelTaxa,0.05/dim(levelTaxa)[2])
}

 # Without Bonferroni for conf.
summary(taxa_temp$Level)
# base   class  domain  family   genus kingdom   order species 
# 2      52       2     114      91      27      81       4

for (level in c("domain","kingdom","class","order","family","genus")){
  levelTaxa = as.data.frame(t(taxa_temp[taxa_temp$Level==level,-c(1:2)]))
  printANOVA1Factor(mdt.temp$temp,levelTaxa,0.05)
}

# ---- Prec. ----

mdt.prec = mdt[mdt$turf.treat!="Warmer",]
mdt.prec$prec = paste(mdt.prec$origSite.precipLevel,mdt.prec$turf.treat.group)
taxa_prec = taxa_turf[rowSums(taxa_turf[-c(1:2)])>minAbundance,c(T,T,mdt$turf.treat!="Warmer")]
summary(taxa_prec$Level)
# base   class  domain  family   genus kingdom   order species 
# 2      70       2     214     246      38     144      31 

for (level in c("domain","kingdom","class","order","family","genus")){
  levelTaxa = as.data.frame(t(taxa_prec[taxa_prec$Level==level,-c(1:2)]))
  printANOVA1Factor(mdt.prec$prec,levelTaxa,0.05/dim(levelTaxa)[2])
}


# Without Bonferroni for confirmation
for (level in c("domain","kingdom","class","order","family","genus")){
  levelTaxa = as.data.frame(t(taxa_prec[taxa_prec$Level==level,-c(1:2)]))
  printANOVA1Factor(mdt.prec$prec,levelTaxa,0.05)
}


