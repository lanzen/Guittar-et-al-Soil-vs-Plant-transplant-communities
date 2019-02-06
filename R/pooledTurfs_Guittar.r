setwd("~/projects/Guittar/Reanalysis_Anders_2016/Clustering_all_ex_mock_20161208")
require(vegan)
# ----- New pooled parameters ----------

mdp = read.csv("../metadata_pooledTurfs.csv",header=T,row.names=1)

# ---- Read BC data ------

micBC = read.csv("../micPooledBC.csv", header=T, row.names = 1)
plantBC = read.csv("../plantPooled.csv", header=T, row.names = 1)

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

micPDTC = (micDO-micDL)/micDO
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

