setwd("~/projects/Guittar/Reanalysis_Anders_2016/Clustering_all_ex_mock_20161208")
require(vegan)
source('~/kode/R/dropRareTaxa.R')
source('~/kode/R/diversity.r')
source('~/kode/R/correlationTests.r')
source('~/kode/R/taxaplot.R')
source('~/kode/R/taxacorr.R')

# ----- New pooled parameters ----------

mdp = read.csv("../metadata_pooledTurfs.csv",header=T,row.names=1)
printAllvsAll(mdp[,c(1:3,7:13,17)],a=0.05)
# Nothing

tea1 = lm(TeabagKAvg~temperature+precipitation,data=mdp)
summary(tea1)
#No
tea1 = lm(TeabagKAvg~temperature*precipitation,data=mdp)
#Neither
require(nlme)
tea3t = lme(TeabagKAvg ~ temperature,random=~1|precip.level,data=mdp)
summary(tea3t) #Nope
tea3p = lme(TeabagKAvg ~ precipitation,random=~1|temp.level,data=mdp)
summary(tea3p)


# ---- Pool control datasets (run turfs.r first) -------

OTUs.c = OTUs[mdt$turf.treat.group=="Control",]
OTUs.c = OTUs.c[,colSums(OTUs.c)>0] #13464 / 18112 retained
mdt.c = mdt[mdt$turf.treat.group=="Control",]

OTUs.pt = data.frame(Alrust=colSums(OTUs.c[mdt.c$site=="Alrust",]),
                     Gudmedalen=colSums(OTUs.c[mdt.c$site=="Gudmedalen",]),
                     Hogsete=colSums(OTUs.c[mdt.c$site=="Hogsete",]),
                     Lavisdalen=colSums(OTUs.c[mdt.c$site=="Lavisdalen",]),
                     Rambera=colSums(OTUs.c[mdt.c$site=="Rambera",]),
                     Skjellingahaugen=colSums(OTUs.c[mdt.c$site=="Skjellingahaugen",]),
                     Ulvhaugen=colSums(OTUs.c[mdt.c$site=="Ulvhaugen",]),
                     Veskre=colSums(OTUs.c[mdt.c$site=="Veskre",]))
OTUs.p = t(OTUs.pt)
write.table(OTUs.pt,"Pooled_OTUs.tsv",quote=F,sep="\t")

OTUs.p.ra = decostand(OTUs.p,method="total")

# ---- Pool all abundant OTUs for all turf datasets ------

poolTurfs = function(otus, mdt){
  s = dim(otus)[2]
  pooled = matrix(nrow=18,ncol=s)
  for (i in c(1:18)){
    t=c(rep(1,4),rep(2,4),rep(1,3),rep(2,7))
    ot=c(rep(1,4),rep(2,4),rep(1,3),rep(2,3),rep(1,4))
    p=c(1:4,1:4,2:4,2:4,1:4)
    op=c(1:4,1:4,1:3,1:3,1:4)
    pooled[i,]=colSums(otus[mdt$temp.level==t[i] & mdt$origSite.tempLevel==ot[i] &
                            mdt$precip.level==p[i] & mdt$origSite.precipLevel==op[i],])
    
  }
  pooled = as.data.frame(pooled)
  names(pooled) = names(otus)
  row.names(pooled)=c("T1P1c","T1P2c","T1P3c","T1P4c",
                      "T2P1c","T2P2c","T2P3c","T2P4c",
                      "T1P1P2t","T1P2P3t","T1P3P4t",
                      "T2P1P2t","T2P2P3t","T2P3P4t",
                      "P1T1T2t","P2T1T2t","P3T1T2t","P4T1T2t")
  return(decostand(pooled,method="total"))
}

micPooled = poolTurfs(OTUs.ab,mdt)
micBC = as.matrix(vegdist(micPooled))
write.csv(micBC,"../micPooledBC.csv")

plantPooled = poolTurfs(plantCompT2[,c(10:180)],mdt)
plantBC = as.matrix(vegdist(plantPooled))
write.csv(plantBC,"../plantPooled.csv")

distTravelled = function(dm){
  return(c(dm[9,1],dm[10,2],dm[11,3],
           dm[12,5],dm[13,6],dm[14,7],
           dm[15,1],dm[16,2],dm[17,3],dm[18,4]))
}

distLeft = function(dm){
  return(c(dm[9,2],dm[10,3],dm[11,4],
           dm[12,6],dm[13,7],dm[14,8],
           dm[15,5],dm[16,6],dm[17,7],dm[18,8]))
}

d_vs_o = function(dm){
  return(c(dm[2,1],dm[3,2],dm[4,3],
           dm[6,5],dm[7,6],dm[8,7],
           dm[5,1],dm[6,2],dm[7,3],dm[8,4]))
}

micDT = distTravelled(micBC)
plantDT = distTravelled(plantBC)
micDL = distLeft(micBC)
plantDL = distLeft(plantBC)
micDO = d_vs_o(micBC)
plantDO = d_vs_o(plantBC)

boxplot(micDO,micDT,micDL,notch=T,col="grey",
        cex=0.8,names=c("Destination \nv origin", "Transplant \nv origin",
"Transplant v \ndestination"))
boxplot(plantDO,plantDT,plantDL,notch=T,
        cex=0.8,names=c("Destination \nv origin", "Transplant \nv origin",
                        "Transplant v \ndestination"))

#micPDTC = (micDO-micDL)/micDO

summary(micPDTC)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.08374  0.24240  0.35490  0.31250  0.46060  0.58970 
sd(micPDTC)
# 0.2260317

plantPDTC = (plantDO-plantDL)/plantDO
summary(plantPDTC)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.008735 0.139800 0.198300 0.245500 0.380900 0.476600
sd(plantPDTC)
#0.152642

boxplot(micPDTC,plantPDTC,names=(c("Prokaryotes","Plants")),notch=T,main="PDTC")

# Difference between microbes and plants site-wise
summary(micPDTC-plantPDTC)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.25210 -0.13810 -0.05463  0.06700  0.34570  0.45380 
sd(micPDTC-plantPDTC)
# 0.2695204

# Scatter plot PDTCs
lmPDTC = lm(plantPDTC~micPDTC)
summary(lmPDTC)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)  0.24014    0.09046   2.655    0.029 *
#   micPDTC      0.01711    0.23868   0.072    0.945  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1618 on 8 degrees of freedom
# Multiple R-squared:  0.0006422,	Adjusted R-squared:  -0.1243 
# F-statistic: 0.005141 on 1 and 8 DF,  p-value: 0.9446
plot(micPDTC~plantPDTC,xlab="PDTC(Plants)",ylab="PDTC(Prokaryotes)")
abline(0,1,lt=2)

plot(micDL~micDO,xlab="BC(Destination vs. origin)",ylab="BC(Transplant vs. destination)")
#,xlim=range(0.25,0.9),ylim=range(0.25,0.7))
points(plantDT~plantDO,pch=2)
# --- Diversity -------

#writeDivStats("../diversity_turfs_pooled.csv",OTUs.p)
divtemp = read.csv("../diversity_turfs_pooled.csv",header=T,row.names=1)
divp = divtemp[,c("Rarefied.richness","H","Rarefied.J")]

cor.test(divtemp$Reads,divtemp$Richness) #p=1E-3, cor=0.93
cor.test(divtemp$Reads,divp$Rarefied.richness) #0.02, cor=0.8
cor.test(divtemp$Reads,divp$H) #p=0.16 cor=0.55
cor.test(divtemp$Reads,divp$Rarefied.J) #p=7, cor=.2

printANOVA(mdp[,c("precip.level","temp.level")],divp,a=0.05)
#Nothing
printVS(divp,mdp[,c(1:4,7:17)],a=0.05)
# RR~pH: R2=0.84, tau=0.71 p=0.14
# H'~pH: R2=0.72 tau=.64 p=0.03

# --------------- 
colrain=c("lightgrey","#AABBEE","#3388FF","blue4")
raincol=colrain#(4)


nmds.s = metaMDS(OTUs.p.ra)
ordiplot(nmds.s,display="sites",type="none")
points(nmds.s,display="sites",col=raincol[mdp$precip.level],pch=(18-mdp$temp.level),cex=1.5)
text(nmds.s,display="sites",pos=1, labels=substr(row.names(mdp),start=0,stop=3))

efs1 = envfit(nmds.s,mdp[,c(1,3:4,7:8,14:17)])
efs1
# NMDS1    NMDS2     r2 Pr(>r)   
# temperature         0.91565  0.40197 0.2628  0.422   
# precipitation      -0.21427 -0.97677 0.8118  0.016 * 
#   precip.level       -0.19006 -0.98177 0.8445  0.008 **
#   temp.level          0.92281  0.38526 0.2421  0.488   
# N_Avail..pre.2012. -0.49692  0.86780 0.1682  0.642   
# pH                 -0.87612 -0.48210 0.7681  0.038 *
# TeabagK.2014       -0.97462 -0.22388 0.0170  0.955   
# TeabagK.2015        0.24638  0.96917 0.1787  0.636   
# TeabagK.2016        0.04142  0.99914 0.4029  0.252   
# TeabagKAvg          0.11909  0.99288 0.2033  0.578   

plot(efs1,p.max=0.05)

efs2 = envfit(nmds.s, mdp[,c(8:13),],na.rm=T)
efs2
#Nothing

legend("topleft",ncol=4,col=raincol[c(1:4)],legend=c(1:4),pch=(18-mdt.c$temp.level),cex=1.1,
       title="Precipitation level",box.lwd=0.5)
legend("topright",ncol=4,legend=c(1,2),pch=c(2,1),cex=1.1,
       title="Temperature level",box.lwd=0.5)

adonis(OTUs.p~mdp$precipitation+mdp$pH)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# mdp$precipitation  1   0.26893 0.26893  2.2218 0.23365  0.036 *
#   mdp$pH             1   0.27686 0.27686  2.2873 0.24054  0.035 *
#   Residuals          5   0.60521 0.12104         0.52581         
# Total              7   1.15100                 1.00000         


# Temp not


# ----------- Taxa -------

ra_all = read.table("CREST_Pooled/Relative_Abundance.tsv",sep="\t",
                    header=T,row.names=3)

ra_all = ra_all[rowSums(ra_all[,c(3:10)])>1E-2,]


printIndicators(ra_all,"kingdom",mdp[,c(7:13,17)],"Indicator_p.csv",a=0.05)
printIndicators(ra_all,"class",mdp[,c(7:13,17)],"Indicator_c.csv",a=0.05)
printIndicators(ra_all,"order",mdp[,c(7:13,17)],"Indicator_o.csv",a=0.05)
printIndicators(ra_all,"family",mdp[,c(7:13,17)],"Indicator_f.csv",a=0.05)
printIndicators(ra_all,"genus",mdp[,c(7:13,17)],"Indicator_g.csv",a=0.05)
  
printIndicatorOTUs(OTUs.p.ra,mdp[,c(7:13,17)],"Indicator_OTUs.csv")
