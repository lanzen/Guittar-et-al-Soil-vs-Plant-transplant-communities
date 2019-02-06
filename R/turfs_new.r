setwd("~/projects/Guittar/Reanalysis_Anders_2016/Clustering_all_ex_mock_20161208")
require(vegan)
source('~/kode/R/dropRareTaxa.R')
source('~/kode/R/diversity.r')
source('~/kode/R/correlationTests.r')

# ---- read data ----

minAbundance = 1/3493 #Representing one read out of all classified in sample with fewest reads


md = read.csv("../metadata.csv",header=T,row.names=1)
mdt = md[md$project=="Turfs",]

otus.t = read.table("CREST_Results/All_OTU_table.csv",sep="\t",
                    header=T,row.names=1)
otus.t = otus.t[,md$project=="Turfs"]
rs = rowSums(otus.t[,-dim(otus.t)[2]])

otus.t = otus.t[rs>0,]
  
sum(rowSums(otus.t[,-dim(otus.t)[2]])) 
# 6,507,919 reads / 35,257 OTUs for all
# 1,381,374 reads / 18,999 OTUs for turfs

otus.t.filtered = otus.t[otus.t$classification!="No hits",]
OTUs = t(otus.t.filtered[,-dim(otus.t)[2]])
names(OTUs) = row.names(otus.t.filtered)
sum(OTUs)  # 1376835 reads, 18112 OTUs
sum(rowSums(OTUs)) / sum(rowSums(otus.t[,-dim(otus.t)[2]])) # 99.7 % kept

OTUs.ra = decostand(OTUs,method="total")

OTUs.ab = OTUs[,colMeans(OTUs.ra)>minAbundance] # Leaving only 664 of 18k+ OTUs!!
OTUs.ab.ra = decostand(OTUs.ab,method="total")
OTUs.ab.ra.c = OTUs.ab.ra[mdt$turf.treat.group=="Control",]

# ---------- Plant data --------------

pcTemp = read.csv("../veg_cover_for_microbial_analysis.csv", header=T, row.names=1)
plantCompT1 = pcTemp[pcTemp$Year==2009,]
plantCompT2 = pcTemp[pcTemp$Year==2013,]
row.names(plantCompT1)=plantCompT1$turfID
row.names(plantCompT2)=plantCompT2$turfID
plantCompT1 = plantCompT1[as.character(mdt$turfID),]
plantCompT2 = plantCompT2[as.character(mdt$turfID),]
row.names(plantCompT1) == mdt$turfID
row.names(plantCompT2) == mdt$turfID
plantCompT1.ra = decostand(plantCompT1[,c(10:180)],method="total")
plantCompT2.ra = decostand(plantCompT2[,c(10:180)],method="total")

plantCompBY = pcTemp[,-c(1:9)]

# ---------- ALPHA DIVERISTY ---------

#writeDivStats("../diversity_turfs.csv",OTUs)
divtemp = read.csv("../diversity_turfs.csv",header=T,row.names=1)
div = divtemp[,c("Reads","Rarefied.richness","H","Rarefied.J")]

## ---- Sanity checks ----
cor.test(div$Reads,div$Richness) #p=1.5E-14, cor=0.8
cor.test(div$Reads,div$Rarefied.richness) #0.02, cor=0.31
cor.test(div$Reads,div$H) #p=0.01 cor=0.33
cor.test(div$Reads,div$J) #p=1E-6, cor=-.6
cor.test(div$Reads,div$Rarefied.J) #p=0.01, cor=.34

## ---- Find trend with alt. and precipitation ----

mdt$RR = RR
mdt$H = div$H
mdtp=mdt[,c(10,9,24,25)]

# B: Diversity above and including transplantation vs. not as factor ----

m0 = glm(RR~1,family="quasipoisson",data=mdt) #AIC Inf, res dev 1390
mp1 = glm(RR~precipitation,family="quasipoisson",data=mdt) 
# (Intercept)   7.147e+00  4.143e-02 172.514  < 2e-16 ***
#   precipitation 8.361e-05  2.058e-05   4.063 0.000168 ***


mtp1 = glm(RR~precipitation+temperature,family="quasipoisson",data=mdtp)
summary(mtp1)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    7.407e+00  9.794e-02  75.627  < 2e-16 ***
#   precipitation  8.331e-05  1.933e-05   4.309 7.67e-05 ***
#   temperature   -3.285e-02  1.143e-02  -2.874  0.00594 **

mtp2 = glm(RR~prec.level.comp+temp.level.comp,
           family="quasipoisson",data=mdtp)
summary(mtp2)

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      7.27371    0.06482 112.214  < 2e-16 ***
#   prec.level.comp  0.06413    0.01474   4.351 6.68e-05 ***
#   temp.level.comp -0.08996    0.03446  -2.611   0.0119 *  

mtp1lm = lm(RR~precipitation+temperature,data=mdt)
summary(mtp1lm)

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   1638.80944  146.05081  11.221 2.90e-15 ***
#   precipitation    0.12338    0.02866   4.305 7.77e-05 ***
#   temperature    -48.18294   16.90716  -2.850  0.00634 ** 


anova(m0,mp1,mtp1,mtp2,test="Chisq")
# 1        52    1389.64                          
# 2        51    1058.74  1   330.90 1.521e-05 ***
#   3        50     913.67  1   145.08   0.00418 ** <--
# 4        50     934.81  0   -21.15              

peplot(mod=mtp1lm, v="precipitation")
peplot(mod=mtp1, v="temperature")



# NMDS and envfit ----

nmds = metaMDS(OTUs.ab.ra.c)

colrain=c("lightgrey","#AABBEE","#3388FF","blue4")
raincol=colrain#(4)

ordiplot(nmds,display="sites",cex=1,lwd=0)
points(nmds,display="sites",col=raincol[mdt.c$precip.level],pch=(18-mdt.c$temp.level) ,cex=1.5)
points(nmds,display="sites",col=mdt.c$transcol,pch=(3-mdt.c$temp.level) ,cex=1.6)
legend("topleft",ncol=4,col=raincol[c(1:4)],legend=c(1:4),pch=(18-mdt.c$temp.level),cex=1.1,
       title="Precipitation level",box.lwd=0.5)
legend("topright",ncol=4,legend=c(1,2),pch=c(2,1),cex=1.1,
       title="Temperature level",box.lwd=0.5)

ef=envfit(nmds,mdt.c[,c("precipitation","temperature")])
# NMDS1    NMDS2     r2 Pr(>r)    
# precipitation -0.42654 -0.90447 0.6609  0.001 ***
#   temperature    0.53784 -0.84305 0.2303  0.063 .  

# or for DEseq2:
# precipitation -0.42653 -0.90447 0.6609  0.001 ***
#   temperature    0.53785 -0.84304 0.2303  0.057 . 

# or for removing rare <minAbundance
# NMDS1    NMDS2    r2 Pr(>r)    
# NMDS1    NMDS2     r2 Pr(>r)   
# precipitation -0.44906 -0.89350 0.5810  0.002 **
#   temperature    0.43380 -0.90101 0.3077  0.036 * 
plot(ef,p.max=.05)

# adonis ----

adonis(OTUs.ab.ra.c~precipitation+temperature,data=mdt.c)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# precipitation  1   0.46240 0.46240  4.4233 0.15208  0.001 ***
#   temperature    1   0.38277 0.38277  3.6616 0.12589  0.004 ** 
#   Residuals     21   2.19526 0.10454         0.72203           

## ---- Include also transplants ----


# ----------NMDS and envfit ----

nmds = metaMDS(OTUs.ab.ra)

ordiplot(nmds,display="sites",cex=1,lwd=0,ylim=range(-.75,.5))
points(nmds,display="sites",col=raincol[mdt$precip.level],pch=(18-mdt$temp.level) ,cex=1.5)
points(nmds,display="sites",col=as.vector(mdt$transcol),pch=(3-mdt$temp.level) ,cex=1.6)
legend(-.85,-.57,ncol=4,col=raincol[c(1:4)],legend=c(1:4),pch=(18-mdt$temp.level),cex=1.1,
       title="Precipitation level",box.lwd=0)
legend(-.85,-.72,ncol=1,legend=c("Control","Transplanted from dryer"),col=c("black","green"),pch=1,
       cex=.95,box.lwd=0)
legend(-.13,-.57, ncol=1,legend=c("1 control","2 control","Transplant 1->2"),
       pch=c(2,1,1),cex=1.1,col=c("black","black","red"),title="Temperature level",box.lwd=0)

#text(nmds,pos=1,cex=0.8)

mdt$warmer = mdt$turf.treat=="Warmer"
ef=envfit(nmds,mdt[,c("precipitation","temperature","turf.treat","warmer")])#,"origSite.tempLevel","origSite.precipLevel")])
# NMDS1    NMDS2     r2 Pr(>r)    
# precipitation -0.58119 -0.81377 0.4988  0.001 ***
#   temperature    0.91919 -0.39382 0.1962  0.006 ** 

plot(ef,p.max=.05,col="darkorange",arrow.mul=0.75)
# original site precipiation cannot be included like this because they are co-varying with the transplant site
# (all transplants done only one step!)
# 
# ord<-ordiellipse(nmds, groups=as.factor(mdt$temp.level), display = "sites", kind ="sd", 
#                  conf = 0.95, label = T, col="red",draw="lines",
#                  show.groups=1)

# ----- Model selection or adonis-------


 adonis(OTUs.ab.ra~precipitation+temperature,data=mdt)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# precipitation  1    0.9197 0.91970  9.2640 0.14018  0.001 ***
#   temperature    1    0.6772 0.67723  6.8216 0.10322  0.001 ***
#   Residuals     50    4.9638 0.09928         0.75659


adonis(OTUs.ab.ra~precip.level+temp.level,data=mdt)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# precip.level  1    0.8487 0.84866  8.3606 0.12935  0.001 ***
#   temp.level    1    0.6368 0.63676  6.2731 0.09706  0.001 ***
#   Residuals    50    5.0753 0.10151         0.77359           

adonis(OTUs.ab.ra~origSite.precipLevel+origSite.tempLevel,data=mdt)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# origSite.precipLevel  1    0.7844 0.78439  7.4941 0.11956  0.001 ***
#   origSite.tempLevel    1    0.5429 0.54292  5.1870 0.08275  0.001 ***
#   Residuals            50    5.2334 0.10467         0.79769             

# -->Connditions at destination slightly better strong predictors of community structure

# This is also significant, i.e. temperature level of original transplant
# Or warmer v not transplant. Only transplant v not is not significant i.e. precipitation
# transplants have no significant effect but are better adapted.

adonis(OTUs.ab.ra~precipitation+temperature+origSite.tempLevel,data=mdt)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# precipitation       1    0.9197 0.91970  9.5508 0.14018  0.001 ***
#   temperature         1    0.6772 0.67723  7.0328 0.10322  0.001 ***
#   origSite.tempLevel  1    0.2453 0.24532  2.5476 0.03739  0.015 *  
#   Residuals          49    4.7185 0.09630         0.71920           

adonis(OTUs.ab.ra~precipitation+temperature+warmer,data=mdt)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# precipitation  1    0.9197 0.91970  9.5235 0.14018  0.001 ***
#   temperature    1    0.6772 0.67723  7.0127 0.10322  0.001 ***
#   warmer         1    0.2318 0.23180  2.4003 0.03533  0.025 *  
#   Residuals     49    4.7320 0.09657         0.72126           
# Total         52    6.5607                 1.00000           

mdt$wetter = mdt$turf.treat=="Wetter"

adonis(OTUs.ab.ra~precipitation+temperature+wetter,data=mdt)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# precipitation  1    0.9197 0.91970  9.2731 0.14018  0.001 ***
#   temperature    1    0.6772 0.67723  6.8283 0.10322  0.001 ***
#   wetter         1    0.1041 0.10405  1.0491 0.01586  0.343    


plot(varpart(OTUs.ab.ra,mdt$temperature,mdt$precipitation,mdt$origSite.tempLevel))


adonis(OTUs.ab.ra~site,data=mdt)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# site       7    3.1025 0.44321  5.7672 0.47288  0.001 ***
#   Residuals 45    3.4583 0.07685         0.52712           

adonis(OTUs.ab.ra~turf.origSite,data=mdt)
# Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
# turf.origSite  7    2.0947 0.299240  3.0151 0.31928  0.001 ***
#   Residuals     45    4.4661 0.099246         0.68072                 

adonis(OTUs.ab.ra~site+turf.origSite,data=mdt)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# site           7    3.1025 0.44321  6.1616 0.47288  0.001 ***
#   turf.origSite  7    0.7249 0.10356  1.4397 0.11049  0.030 *  
#   Residuals     38    2.7334 0.07193         0.41662           

# --> Site as such has stronger influence than origin site but both influence




# ------ Boxplots for type of changes -----

distPartitioning = function(meta, community, col="grey"){
  
  require(vioplot)
  
  mdt.o = meta[order(meta$temp.level,meta$precip.level,meta$turf.treat),]
  otus.o = community[order(meta$temp.level,meta$precip.level,meta$turf.treat),]
  
  mdt.o$rep=c(1:53)
  mdt.o$name=paste(mdt.o$rep, "T",mdt.o$origSite.tempLevel, 
                   mdt.o$temp.level,"P",mdt.o$origSite.precipLevel, mdt.o$precip.level)
  row.names(otus.o)=mdt.o$name
  
  vd = as.matrix(vegdist(otus.o))
  
  controls = otus.o[mdt.o$turf.treat.group=="Control",]
  trans = otus.o[mdt.o$turf.treat.group=="Transplant",]
  
  cd = as.matrix(vegdist(controls))
  #write.csv(cd,"control_dists_ordered.csv")
  
  chetero = matrix(nrow=8,ncol=3)
  for (s in c(0:7)) {
    i = s*3+1
    chetero[s+1,] = c(cd[i,i+1],cd[i,i+2],cd[i+1,i+2])
  }
  chetero = as.data.frame(chetero)
  row.names(chetero)=c("T1P1","T1P2","T1P3","T1P4","T2P1","T2P2","T2P3","T2P4")
  #print(mean(unlist(chetero)),na.rm=T)
  control_hetero_sitewise = rowMeans(chetero)
  
  td = as.matrix(vegdist(trans))
  write.csv(td,"transplant_dists_ordered.csv")
  
  thetero = matrix(nrow=10,ncol=3)
  for (s in c(0:1)) {
    i = s*3+1
    thetero[s+1,] = c(td[i,i+1],td[i,i+2],td[i+1,i+2])
  }
  thetero[3,] = c(td[7,8],NaN,NaN)
  for (s in c(3:9)) {
    i = s*3
    thetero[s+1,] = c(td[i,i+1],td[i,i+2],td[i+1,i+2])
  }
  
  thetero = as.data.frame(thetero)
  row.names(thetero)=c("T1P1-2","T1P2-3","T1P3-4","T1-2P1","T1-2P2",
                       "T2P1-2","T1-2P3","T2P2-3","T1-2P4","T2P3-4")
  
  o_v_d = matrix(nrow=10,ncol=9)
  for (s in c(0:2)) {
    i = s*3+1
    o_v_d[s+1,] = c(cd[i,c((i+3):(i+5))],cd[i+1,c((i+3):(i+5))],cd[i+2,c((i+3):(i+5))])
    o_v_d[s+4,] = c(cd[i,c((i+12):(i+14))],cd[i+1,c((i+12):(i+14))],cd[i+2,c((i+12):(i+14))])
  }
  o_v_d[7,] = c(cd[10,c(22:24)],cd[11,c(22:24)],cd[12,22:24])
  
  for (s in c(4:6)) {
    i = s*3+1
    o_v_d[s+4,] = c(cd[i,c((i+3):(i+5))],cd[i+1,c((i+3):(i+5))],cd[i+2,c((i+3):(i+5))])
  }
  
  o_v_d = as.data.frame(o_v_d)
  row.names(o_v_d)=c("T1P1-2","T1P2-3","T1P3-4","T1-2P1","T1-2P2","T1-2P3","T1-2P4",
                      "T2P1-2","T2P2-3","T2P3-4")
  o_v_d_trans_wise = rowMeans(o_v_d)
  
  t_v_o = matrix(nrow=10,ncol=9)
  
  temp=matrix(nrow=3,ncol=3)
  temp = vd[c(1:3),c(7:9)]
  t_v_o[1,]=c(temp[1,],temp[2,],temp[3,])
  temp = vd[c(4:6),c(13:15)]
  t_v_o[2,]=c(temp[1,],temp[2,],temp[3,])
  temp = vd[c(19:20),c(10:12)]
  t_v_o[3,]=c(temp[1,],temp[2,],rep(NaN,3))
  temp = vd[c(24:26),c(1:3)]
  t_v_o[4,]=c(temp[1,],temp[2,],temp[3,])
  temp = vd[c(30:32),c(4:6)]
  t_v_o[5,]=c(temp[1,],temp[2,],temp[3,])
  temp = vd[c(39:41),c(10:12)]
  t_v_o[6,]=c(temp[1,],temp[2,],temp[3,])
  temp = vd[c(48:50),c(16:18)]
  t_v_o[7,]=c(temp[1,],temp[2,],temp[3,])
  temp = vd[c(33:35),c(21:23)]
  t_v_o[8,]=c(temp[1,],temp[2,],temp[3,])
  temp = vd[c(42:44),c(27:29)]
  t_v_o[9,]=c(temp[1,],temp[2,],temp[3,])
  temp = vd[c(51:53),c(38:40)]
  t_v_o[10,]=c(temp[1,],temp[2,],temp[3,])
  t_v_o = as.data.frame(t_v_o)
  row.names(t_v_o)=c("T1P1-2","T1P2-3","T1P3-4","T1-2P1","T1-2P2","T1-2P3","T1-2P4",
                     "T2P1-2","T2P2-3","T2P3-4")
  
  
  t_v_d = matrix(nrow=10,ncol=9)
  
  temp = vd[c(4:6),c(7:9)]
  t_v_d[1,]=c(temp[1,],temp[2,],temp[3,])
  temp = vd[c(10:12),c(13:15)]
  t_v_d[2,]=c(temp[1,],temp[2,],temp[3,])
  temp = vd[c(19:20),c(16:18)]
  t_v_d[3,]=c(temp[1,],temp[2,],rep(NaN,3))
  temp = vd[c(24:26),c(21:23)]
  t_v_d[4,]=c(temp[1,],temp[2,],temp[3,])
  temp = vd[c(30:32),c(27:29)]
  t_v_d[5,]=c(temp[1,],temp[2,],temp[3,])
  temp = vd[c(39:41),c(36:38)]
  t_v_d[6,]=c(temp[1,],temp[2,],temp[3,])
  temp = vd[c(48:50),c(45:47)]
  t_v_d[7,]=c(temp[1,],temp[2,],temp[3,])
  temp = vd[c(33:35),c(30:32)]
  t_v_d[8,]=c(temp[1,],temp[2,],temp[3,])
  temp = vd[c(42:44),c(39:41)]
  t_v_d[9,]=c(temp[1,],temp[2,],temp[3,])
  temp = vd[c(51:53),c(48:50)]
  t_v_d[10,]=c(temp[1,],temp[2,],temp[3,])
  t_v_d = as.data.frame(t_v_d)
  row.names(t_v_d)=c("T1P1-2","T1P2-3","T1P3-4","T1-2P1","T1-2P2","T1-2P3","T1-2P4",
                     "T2P1-2","T2P2-3","T2P3-4")
  
  t_v_d_transwise = rowMeans(t_v_d, na.rm=T)
 
  inside_control=c(chetero$V1,chetero$V2,chetero$V3)
  inside_trans = c(thetero$V1,thetero$V2,thetero$V3)
  origin_vs_dest = c(o_v_d$V1,o_v_d$V2,o_v_d$V3,o_v_d$V4,o_v_d$V5,o_v_d$V6,o_v_d$V7,o_v_d$V8,o_v_d$V9)
  origin_vs_trans = c(t_v_o$V1,t_v_o$V2,t_v_o$V3,t_v_o$V4,t_v_o$V5,t_v_o$V6,t_v_o$V7,t_v_o$V8,t_v_o$V9)
  dest_vs_trans =  c(t_v_d$V1,t_v_d$V2,t_v_d$V3,t_v_d$V4,t_v_d$V5,t_v_d$V6,t_v_d$V7,t_v_d$V8,t_v_o$V9)
  
  print("Transplant vs. dest mean:")
  m_t_v_d = mean(dest_vs_trans,na.rm=T)
  print(m_t_v_d)
  print ("origin vs. dest mean:")
  m_o_v_d = mean(origin_vs_dest,na.rm=T)
  print(m_o_v_d)
  print ("Proportional distance travelled towards convergence")
  print((m_o_v_d-m_t_v_d)/m_o_v_d)
  
  boxplot(origin_vs_dest,origin_vs_trans,dest_vs_trans,inside_control,inside_trans,notch=T,
          cex=0.8,names=c("Destination \nv origin",
                  "Transplant \nv origin","Transplant v \ndestination",
                  "Heterogeneity \ncontrols","Heterogeneity \ntransplants"),col=col)
  line <- readline()      
  
  vioplot(origin_vs_dest,as.vector(na.omit(origin_vs_trans)),
          as.vector(na.omit(dest_vs_trans)),
          inside_control,as.vector(na.omit(inside_trans)),
          names=c("Destination \nv origin","Transplant \nv origin","Transplant v \ndestination",
                  "Heterogeneity \ncontrols","Heterogeneity \ntransplants"),
          col=col)
  
  
  
  pdtc = 1 - ( (t_v_d_transwise - control_hetero_sitewise[c(2:8,6:8)]) / o_v_d_trans_wise)
  return (pdtc)
}

distPartitioning_hetreoAdjusted = function(meta, community, col="grey"){
  require(vioplot)
  
  mdt.o = meta[order(meta$temp.level,meta$precip.level,meta$turf.treat),]
  otus.o = community[order(meta$temp.level,meta$precip.level,meta$turf.treat),]
  
  mdt.o$rep=c(1:53)
  mdt.o$name=paste(mdt.o$rep, "T",mdt.o$origSite.tempLevel, 
                   mdt.o$temp.level,"P",mdt.o$origSite.precipLevel, mdt.o$precip.level)
  row.names(otus.o)=mdt.o$name
  
  vd = as.matrix(vegdist(otus.o))
  
  controls = otus.o[mdt.o$turf.treat.group=="Control",]
  trans = otus.o[mdt.o$turf.treat.group=="Transplant",]
  
  cd = as.matrix(vegdist(controls))
  #write.csv(cd,"control_dists_ordered.csv")
  
  chetero = matrix(nrow=8,ncol=3)
  for (s in c(0:7)) {
    i = s*3+1
    chetero[s+1,] = c(cd[i,i+1],cd[i,i+2],cd[i+1,i+2])
  }
  chetero = as.data.frame(chetero)
  row.names(chetero)=c("T1P1","T1P2","T1P3","T1P4","T2P1","T2P2","T2P3","T2P4")
  #print(mean(unlist(chetero)),na.rm=T)
  control_hetero_sitewise = rowMeans(chetero)
  
  td = as.matrix(vegdist(trans))
  write.csv(td,"transplant_dists_ordered.csv")
  
  thetero = matrix(nrow=10,ncol=3)
  for (s in c(0:1)) {
    i = s*3+1
    thetero[s+1,] = c(td[i,i+1],td[i,i+2],td[i+1,i+2])
  }
  thetero[3,] = c(td[7,8],NaN,NaN)
  for (s in c(3:9)) {
    i = s*3
    thetero[s+1,] = c(td[i,i+1],td[i,i+2],td[i+1,i+2])
  }
  
  thetero = as.data.frame(thetero)
  row.names(thetero)=c("T1P1-2","T1P2-3","T1P3-4","T1-2P1","T1-2P2",
                       "T2P1-2","T1-2P3","T2P2-3","T1-2P4","T2P3-4")
  
  o_v_d = matrix(nrow=10,ncol=9)
  for (s in c(0:2)) {
    i = s*3+1
    o_v_d[s+1,] = c(cd[i,c((i+3):(i+5))],cd[i+1,c((i+3):(i+5))],cd[i+2,c((i+3):(i+5))])
    o_v_d[s+4,] = c(cd[i,c((i+12):(i+14))],cd[i+1,c((i+12):(i+14))],cd[i+2,c((i+12):(i+14))])
  }
  o_v_d[7,] = c(cd[10,c(22:24)],cd[11,c(22:24)],cd[12,22:24])
  
  for (s in c(4:6)) {
    i = s*3+1
    o_v_d[s+4,] = c(cd[i,c((i+3):(i+5))],cd[i+1,c((i+3):(i+5))],cd[i+2,c((i+3):(i+5))])
  }
  
  o_v_d = as.data.frame(o_v_d)
  row.names(o_v_d)=c("T1P1-2","T1P2-3","T1P3-4","T1-2P1","T1-2P2","T1-2P3","T1-2P4",
                     "T2P1-2","T2P2-3","T2P3-4")
  o_v_d_trans_wise = rowMeans(o_v_d)
  
  t_v_o = matrix(nrow=10,ncol=9)
  
  temp=matrix(nrow=3,ncol=3)
  temp = vd[c(1:3),c(7:9)]
  t_v_o[1,]=c(temp[1,],temp[2,],temp[3,])
  temp = vd[c(4:6),c(13:15)]
  t_v_o[2,]=c(temp[1,],temp[2,],temp[3,])
  temp = vd[c(19:20),c(10:12)]
  t_v_o[3,]=c(temp[1,],temp[2,],rep(NaN,3))
  temp = vd[c(24:26),c(1:3)]
  t_v_o[4,]=c(temp[1,],temp[2,],temp[3,])
  temp = vd[c(30:32),c(4:6)]
  t_v_o[5,]=c(temp[1,],temp[2,],temp[3,])
  temp = vd[c(39:41),c(10:12)]
  t_v_o[6,]=c(temp[1,],temp[2,],temp[3,])
  temp = vd[c(48:50),c(16:18)]
  t_v_o[7,]=c(temp[1,],temp[2,],temp[3,])
  temp = vd[c(33:35),c(21:23)]
  t_v_o[8,]=c(temp[1,],temp[2,],temp[3,])
  temp = vd[c(42:44),c(27:29)]
  t_v_o[9,]=c(temp[1,],temp[2,],temp[3,])
  temp = vd[c(51:53),c(38:40)]
  t_v_o[10,]=c(temp[1,],temp[2,],temp[3,])
  t_v_o = as.data.frame(t_v_o)
  row.names(t_v_o)=c("T1P1-2","T1P2-3","T1P3-4","T1-2P1","T1-2P2","T1-2P3","T1-2P4",
                     "T2P1-2","T2P2-3","T2P3-4")
  
  
  t_v_d = matrix(nrow=10,ncol=9)
  
  temp = vd[c(4:6),c(7:9)]
  t_v_d[1,]=c(temp[1,],temp[2,],temp[3,])
  temp = vd[c(10:12),c(13:15)]
  t_v_d[2,]=c(temp[1,],temp[2,],temp[3,])
  temp = vd[c(19:20),c(16:18)]
  t_v_d[3,]=c(temp[1,],temp[2,],rep(NaN,3))
  temp = vd[c(24:26),c(21:23)]
  t_v_d[4,]=c(temp[1,],temp[2,],temp[3,])
  temp = vd[c(30:32),c(27:29)]
  t_v_d[5,]=c(temp[1,],temp[2,],temp[3,])
  temp = vd[c(39:41),c(36:38)]
  t_v_d[6,]=c(temp[1,],temp[2,],temp[3,])
  temp = vd[c(48:50),c(45:47)]
  t_v_d[7,]=c(temp[1,],temp[2,],temp[3,])
  temp = vd[c(33:35),c(30:32)]
  t_v_d[8,]=c(temp[1,],temp[2,],temp[3,])
  temp = vd[c(42:44),c(39:41)]
  t_v_d[9,]=c(temp[1,],temp[2,],temp[3,])
  temp = vd[c(51:53),c(48:50)]
  t_v_d[10,]=c(temp[1,],temp[2,],temp[3,])
  t_v_d = as.data.frame(t_v_d)
  row.names(t_v_d)=c("T1P1-2","T1P2-3","T1P3-4","T1-2P1","T1-2P2","T1-2P3","T1-2P4",
                     "T2P1-2","T2P2-3","T2P3-4")
  
  t_v_d_transwise = rowMeans(t_v_d, na.rm=T) 
  
  inside_control=c(chetero$V1,chetero$V2,chetero$V3)
  avg_ic = median(inside_control)
  print ("Control heterogeneity median:")
  print(avg_ic)
  inside_control = inside_control - avg_ic
  
  inside_trans = c(thetero$V1,thetero$V2,thetero$V3) - avg_ic
  origin_vs_dest = c(o_v_d$V1,o_v_d$V2,o_v_d$V3,o_v_d$V4,o_v_d$V5,o_v_d$V6,o_v_d$V7,o_v_d$V8,o_v_d$V9) - avg_ic
  origin_vs_trans = c(t_v_o$V1,t_v_o$V2,t_v_o$V3,t_v_o$V4,t_v_o$V5,t_v_o$V6,t_v_o$V7,t_v_o$V8,t_v_o$V9) - avg_ic
  dest_vs_trans =  c(t_v_d$V1,t_v_d$V2,t_v_d$V3,t_v_d$V4,t_v_d$V5,t_v_d$V6,t_v_d$V7,t_v_d$V8,t_v_o$V9) - avg_ic
  
  print("Transplant vs. dest median:")
  m_t_v_d = median(dest_vs_trans,na.rm=T)
  print(m_t_v_d)
  print ("origin vs. dest median:")
  m_o_v_d = median(origin_vs_dest,na.rm=T)
  print(m_o_v_d)
  
  boxplot(origin_vs_dest,origin_vs_trans,dest_vs_trans,inside_control,inside_trans,notch=T,
          cex=0.8,names=c("Destination \nv origin",
                          "Transplant \nv origin","Transplant v \ndestination",
                          "Heterogeneity \ncontrols","Heterogeneity \ntransplants"),col=col)
  abline(0,0, col="grey", lty=2)
  line <- readline()      
  
  vioplot(origin_vs_dest,as.vector(na.omit(origin_vs_trans)),
          as.vector(na.omit(dest_vs_trans)),
          inside_control,as.vector(na.omit(inside_trans)),
          names=c("Destination \nv origin","Transplant \nv origin","Transplant v \ndestination",
                  "Heterogeneity \ncontrols","Heterogeneity \ntransplants"),
          col=col)
  abline(0,0, col="grey", lty=2)
  
  
  
  
  pdtc = 1 - ( (t_v_d_transwise - control_hetero_sitewise[c(2:8,6:8)]) / o_v_d_trans_wise)
  return (pdtc)
}


#pdf("../img/16SDists.pdf",height = 6, width=8)
micPDTC = 100*distPartitioning_hetreoAdjusted(mdt, OTUs.ab.ra,col="grey")
print(paste(mean(micPDTC),sd(micPDTC),sep="±"))

plantPDTC = 100*distPartitioning_hetreoAdjusted(mdt, plantCompT2[,-c(1:9)],col="white")
print(paste(mean(plantPDTC),sd(plantPDTC),sep="±"))


plot(micPDTC~plantPDTC,
     main="PDTC (%)", col=raincol[c(2:4,1:4,2:4)], pch=c(rep(17,3),rep(16,7)),
     xlab="Prokaryotes", ylab="Plants", cex=1.5,
     xlim=range(60,105), ylim=range(60,105))
points(micPDTC~plantPDTC,pch=c(rep(2,3),rep(1,7)),
       col=c(rep("green",3),rep("red",4),rep("green",3)),cex=1.6)
abline(0,1,col="grey")
legend("bottomright",ncol=4,col=raincol[c(1:4)],legend=c(1:4),pch=(18-mdt$temp.level),
       title="Precipitation level",box.lwd=0)
legend("bottom", ncol=2,legend=c("T1","T2","T+","P+"),
       ,pch=c(17,16,1,1),col=c("black","black","red","green"),box.lwd=0)
#abline(v=0,col="grey")
#abline(h=0,col="grey")

pdf("../img/PlantDists_2009.pdf",height = 6, width=8)
distPartitioning_hetreoAdjusted(mdt, plantCompT1.ra,col="darkgrey")
dev.off()

# # Warmer or wetter only
# select=c("T1-2P1","T1-2P2","T1-2P3","T1-2P4") # warmer
# select=c("T1P1-2","T1P2-3","T1P3-4","T2P1-2","T2P2-3","T2P3-4") #wetter
# th = thetero[select,]
# inside_trans = c(th$V1,th$V2,th$V3)
# od = o_v_d[select,]
# origin_vs_dest = c(od$V1,od$V2,od$V3,od$V4,od$V5,od$V6,od$V7,od$V8,od$V9)
# ot = t_v_o[select,]
# origin_vs_trans = c(ot$V1,ot$V2,ot$V3,ot$V4,ot$V5,ot$V6,ot$V7,ot$V8,ot$V9)
# dt = t_v_d[select,]
# dest_vs_trans =  c(dt$V1,dt$V2,dt$V3,dt$V4,dt$V5,dt$V6,dt$V7,dt$V8,t_v_o$V9)


# ------- Correlation to soil parameters --------

para = read.csv("../selected.csv",row.names=1,header=T)
summary(para)
OTUs.s = OTUs.ra[row.names(para),]
mdt.s = mdt[row.names(para),]

nmds.s = metaMDS(OTUs.s)
ordiplot(nmds.s,display="sites",type="none")
points(nmds.s,display="sites",col=mdt.s$colour,pch=(18-mdt.s$temp.level),cex=1.5)

legend("bottomleft",ncol=2,col=c(1:8),legend=c("Veskre","Hogsete","Alrust","Lavisdalen",
                                               "Ulvhaugen","Skjellingahaugen","Rambera","Gudmedalen"),pch=c(16,16,16,17,17,17,16,17))
efs1 = envfit(nmds.s,para[,c(3:8)])
efs1
plot(efs1,p.max=0.05)
efs2 = envfit(nmds.s,para[,c(9:15),],na.rm=T)
efs2

adonis(OTUs.s~para$pH)


# ------- Functional results - Tax4Fun -----
mdt$prec.level.comp = mdt$precip.level - as.numeric(mdt$turf.treat=="Wetter")/2
mdt$temp.level.comp = mdt$temp.level - as.numeric(mdt$turf.treat=="Warmer")/2

t4f = read.table("../../Tax4Fun/Tax4FunProfile_Export.csv",sep="\t",header=T,row.names=1)
t4f_turf = t4f[md$project=="Turfs",]

lmBothT4F = function(ra,file_sign,obs){
  
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
  
  write(c(o,mean(ra),
          sig,pp,pp*obs,pt,pt*obs,tp,tt),
        file=file_sign,ncolumns=10,sep=",",append=T)
}

sfile="../lm_taxa4fun.csv"

write(c("Category","Avg.abundance","Significance",
        "p.prec","p.prec.bonf","p.temp","p.temp.bonf","t.prec","t.temp"),
      file=sfile,ncolumns=10,sep=",")

obs = dim(t4f_turf)[2]
  
for (o in names(t4f_turf)){
    ra = as.numeric(t4f_turf[,o])
    lmBoth(ra,file_sign=sfile,obs=obs)
  }
}


# ---- Correlation with plant compo -------

mantel(vegdist(plantCompT2.ra),vegdist(OTUs.ab.ra))
# R=.31, p<0.001

mantel(vegdist(plantCompT2.ra),vegdist(plantCompT1.ra))
# R=0.59, p<0.001

plantCompBYC = pcTemp[pcTemp$TTtreat=="TT1",]
plantCompBYT = pcTemp[pcTemp$TTtreat!="TT1",]

anosim(plantCompBYC[,-c(1:9)], plantCompBYC$Year)
# ANOSIM statistic R: -0.02587 
# Significance: 0.827 

anosim(plantCompBYT[,-c(1:9)], plantCompBYT$Year)
# NOSIM statistic R: 0.08088 
# Significance: 0.004 

plantNMDS = metaMDS(plantCompBY)

ordiplot(plantNMDS,display="sites",cex=1,lwd=0,ylim=range(-.75,.5))
points(plantNMDS,display="sites",col=raincol[pcTemp$Precipitation_level],
       pch=(18-pcTemp$Temperature_level) ,cex=1.5)
points(plantNMDS,display="sites",col=(2013-pcTemp$Year)/4+8,
        pch=(3-pcTemp$Temperature_level) ,cex=1.6)

legend("bottomleft",ncol=4,col=raincol[c(1:4)],legend=c(1:4),pch=(18-mdt$temp.level),
       title="Precipitation level",box.lwd=0)

legend("bottomright", ncol=1,legend=c("2009 T1","2009 T2", "2013 T1", "2013 T2"),
        col=c(9,9,8,8),pch=c(2,1,2,1),box.lwd=0)

