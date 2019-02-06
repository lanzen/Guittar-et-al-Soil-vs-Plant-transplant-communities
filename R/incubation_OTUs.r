setwd("~/projects/Guittar/Reanalysis_Anders_2016/Clustering_all_ex_mock_20161208")
require(vegan)
source('~/kode/R/dropRareTaxa.R')
source('~/kode/R/diversity.r')
source('~/kode/R/correlationTests.r')

md = read.csv("../metadata.csv",header=T,row.names=1)
mdi = md[md$project!="Turfs",]
mdi$tsymbol = 15
mdi[mdi$exp.temp==20,]$tsymbol=18
# Only temp level 1
mdi.dna = mdi[mdi$dna.rna=="DNA",]
mdi.rna = mdi[mdi$dna.rna=="RNA",]

otus.t = read.table("CREST_Results/All_OTU_table.csv",sep="\t",
                    header=T,row.names=1)
otus.it = otus.t[,md$project!="Turfs"]

rs = rowSums(otus.it) # remove classification field

otus.it.filtered = otus.it[otus.t$classification!="No hits",]
OTUs = t(otus.it.filtered)
names(OTUs) = row.names(otus.it.filtered)
OTUs.ra = decostand(OTUs,method="total")
OTUs.ab = OTUs[,colMeans(OTUs.ra)>1/3493] # Leaving only 619 of 30k+ OTUs!!
OTUs.ab.ra =  decostand(OTUs.ab,method="total")

divtemp = read.csv("../diversity.csv",header=T,row.names=1)
div = divtemp[md$project!="Turfs",c("Reads","Rarefied.richness","H","Rarefied.J")]

# --- 

OTUs.ab.ra.rna = OTUs.ab.ra[mdi$dna.rna=="RNA",]
OTUs.ab.ra.dna = OTUs.ab.ra[mdi$dna.rna=="DNA",] 

# --- NMDS ----

raincol=c("lightgrey","#AABBEE","#3388FF","blue4")

nmds = metaMDS(OTUs.ab.ra.dna)
ordiplot(nmds,display="sites",cex=1,lwd=0)
points(nmds,display="sites",col=raincol[mdi.dna$exp.moisture.level],pch=mdi.dna$tsymbol,cex=1.6)
points(nmds,display="sites",col=raincol[mdi.dna$origSite.precipLevel],pch=mdi.dna$tsymbol,cex=1)
text(nmds,display="sites",col=raincol[mdi.dna$exp.moisture.level],labels=mdi.dna$day,pos=1,cex=0.7)
legend("bottomleft",legend=c("10 C", "20 C"),pch=c(15,18),cex=1.1,title="Temperature",box.lwd=0.5)
ef = envfit(nmds,mdi.dna[,c("origSite.precipLevel","day","exp.temp","exp.moisture.level")])
# NMDS1    NMDS2     r2 Pr(>r)    
# origSite.precipLevel  0.49307  0.86999 0.5862  0.001 ***
plot(ef,p.max=.05)

nmds = metaMDS(OTUs.ab.ra.rna)
ordiplot(nmds,display="sites",cex=1,lwd=0)
points(nmds,display="sites",col=raincol[mdi.rna$exp.moisture.level],pch=mdi.rna$tsymbol,cex=1.6)
points(nmds,display="sites",col=raincol[mdi.rna$origSite.precipLevel],pch=mdi.rna$tsymbol,cex=1)
text(nmds,display="sites",col=raincol[mdi.rna$exp.moisture.level],labels=mdi.rna$day,pos=1,cex=0.7)
legend("bottomleft",legend=c("10 C", "20 C"),pch=c(15,18),cex=1.1,title="Temperature",box.lwd=0.5)
ef = envfit(nmds,mdi.rna[,c("origSite.precipLevel","day","exp.temp","exp.moisture.level")])

# NMDS1    NMDS2     r2 Pr(>r)    
# origSite.precipLevel  0.98024 -0.19780 0.9798  0.001 ***
#   day                  -0.02878  0.99959 0.1299  0.012 *  
#   exp.temp             -0.02462  0.99970 0.0145  0.631    
# exp.moisture.level    0.06330  0.99799 0.3159  0.001 ***
#   
plot(ef,p.max=.05)

# ---- adonis ----

adonis(OTUs.ab.ra.dna~origSite.precipLevel+day*exp.moisture.level+day*exp.temp,data=mdi.dna)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# origSite.precipLevel     1    3.2902  3.2902  41.850 0.18385  0.001 ***
#   day                      1    0.9406  0.9406  11.964 0.05256  0.001 ***
#   exp.moisture.level       1    0.4259  0.4259   5.418 0.02380  0.001 ***
#   exp.temp                 1    0.2206  0.2206   2.806 0.01233  0.015 *  
#   day:exp.moisture.level   1    0.2670  0.2670   3.396 0.01492  0.010 ** 
#   day:exp.temp             1    0.1726  0.1726   2.195 0.00964  0.037 *  
#   Residuals              160   12.5789  0.0786         0.70289           
# Total                  166   17.8958                 1.00000           

adonis(OTUs.ab.ra.rna~origSite.precipLevel+day*exp.moisture.level+day*exp.temp,data=mdi.rna)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# origSite.precipLevel    1    2.6810 2.68103  78.158 0.48480  0.001 ***
#   day                     1    0.1967 0.19673   5.735 0.03557  0.001 ***
#   exp.moisture.level      1    0.3408 0.34079   9.935 0.06162  0.001 ***
#   exp.temp                1    0.0474 0.04736   1.381 0.00856  0.204    
# day:exp.moisture.level  1    0.0895 0.08954   2.610 0.01619  0.038 *  
#   day:exp.temp            1    0.0480 0.04802   1.400 0.00868  0.210    
# Residuals              62    2.1268 0.03430         0.38457           
# Total                  68    5.5302                 1.00000           

# temperature has nothing on RNA (strange)
adonis(OTUs.ab.ra.rna~origSite.precipLevel+day*exp.moisture.level,data=mdi.rna)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# origSite.precipLevel    1    2.6810 2.68103  77.213 0.48480  0.001 ***
#   day                     1    0.1967 0.19673   5.666 0.03557  0.003 ** 
#   exp.moisture.level      1    0.3408 0.34079   9.815 0.06162  0.001 ***
#   day:exp.moisture.level  1    0.0894 0.08944   2.576 0.01617  0.036 *  
#   Residuals              64    2.2222 0.03472         0.40184           
# Total                  68    5.5302                 1.00000       


