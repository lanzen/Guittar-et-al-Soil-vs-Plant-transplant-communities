printANOVA(mdp[,c("precip.level","temp.level")],divp,a=0.05)
#Nothing

printVS(divp,mdp[,c("temperature","precipitation")],a=0.05)
#Nothing

m0 = glm(divp$Rarefied.richness~1,family="quasipoisson",data=mdp) #AIC Inf, res dev 670
mt1 = glm(divp$Rarefied.richness~temperature,family="quasipoisson",data=mdp) 
summary(mt1)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  8.89369    0.27490  32.352  5.8e-08 ***
#   temperature -0.07848    0.03628  -2.163   0.0737 .  

anova(m0,mt1,test="Chisq") # -> m1 best i.e. temperature explains richness!


lmt1 = lm(divp$Rarefied.richness~temperature,data=mdp) 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)   6436.9     1129.2    5.70  0.00126 **
#   temperature   -316.0      146.2   -2.16  0.07404 . 
# A linear model is also fine btut still not sign....


# ----

mp1 = glm(divp$Rarefied.richness~precipitation,family="quasipoisson",data=mdp) 
summary(mp1)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   8.117e+00  1.311e-01  61.918 1.19e-09 ***
#   precipitation 1.072e-04  6.731e-05   1.593    0.162    

lmp1 = lm(divp$Rarefied.richness~precipitation,data=mdp) 
summary(lmp1)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   3291.9142   517.9495   6.356 0.000712 ***
#   precipitation    0.4377     0.2755   1.589 0.163169    

# not significant

# pH also ----

lmpH = lm(divp$Rarefied.richness~pH,data=mdp)
summary(lmpH)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -6348.5     1674.8  -3.791 0.009068 ** 
#   pH            1845.5      297.2   6.210 0.000804 ***

lmboth = lm(divp$Rarefied.richness~precipitation+temperature,data=mdp)
summary(lmboth)
# Estimate Std. Error t value Pr(>|t|)   
# (Intercept)   5722.0973   867.7251   6.594   0.0012 **
#   precipitation    0.4517     0.1789   2.524   0.0529 . 
# temperature   -322.8467   106.2651  -3.038   0.0288 * 

lmpHtemp = lm(divp$Rarefied.richness~pH+temperature,data=mdp)
summary(lmpHtemp)
# only pH

# --- Plots ----

# Temp on x axis
plot(divp$Rarefied.richness~temperature,data=mdp,col=raincol[mdp$precip.level],
     pch=(18-mdp$temp.level),ylab="Rarefied OTU richness",cex=1.2)
abline(lmt1,lty=2)


# Prec on x axis

plot(divp$Rarefied.richness~precipitation,data=mdp,col=raincol[mdp$precip.level],
     pch=(18-mdp$temp.level),ylab="Rarefied OTU richness",cex=1.2)
abline(lmp1,lty=2)

# pH on x axis

plot(divp$Rarefied.richness~pH,data=mdp,col=raincol[mdp$precip.level],
     pch=(18-mdp$temp.level),ylab="Rarefied OTU richness",cex=1.2)
abline(lmpH,lty=2)


# ---- H' instead of RR

lmp1 = lm(divp$H~precipitation,data=mdp)
summary(lmp1)
# NS


# ----

lmt1 = lm(divp$H~temperature,data=mdp)
summary(lmt1)
#NS (but almost)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  7.93960    0.36574  21.709 6.24e-07 ***
#   temperature -0.11326    0.04737  -2.391   0.0539 .  


# ----

m0 = lm(divp$H~1,data=mdp) #AIC Inf, res dev 1390
mtgl = lm(divp$H~temperature+precipitation,data=mdp)
mtg = lm(divp$H~temp.level+precip.level,data=mdp)
mtgo = lm(divp$H~origSite.tempLevel+origSite.precipLevel,data=mdp)
anova(m0,mtgl,mtg,mtgo,test="Chisq")
# Res.Df    RSS Df Sum of Sq  Pr(>Chi)    
# 1     52 3.5744                           
# 2     50 2.5303  2   1.04418 3.307e-05 ***
#   3     50 2.7043  0  -0.17401              
# 4     50 2.7663  0  -0.06208                 

# Interesting that original temp. is LESS significant as opposed to for RR

summary(mtgl)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    6.968e+00  2.039e-01  34.168  < 2e-16 ***
#   temperature   -6.078e-02  2.361e-02  -2.575 0.013043 *  
#   precipitation  1.479e-04  4.002e-05   3.696 0.000544 ***
# --- Plots ----

# Temp on x axis
plot(divp$H~temperature,data=mdp,
     pch=(3-mdp$temp.level),ylab="Shannon divpersity (H')",cex=1.2)
points(divp$H~temperature,data=mdp,col=raincol[mdp$precip.level],pch=(18-mdp$temp.level))
abline(lmt1,lty=2)
points(divp.t$H~temperature,data=mdp.t,col=transcol,pch=(3-temp.level),cex=1.2)
points(divp.t$H~temperature,data=mdp.t,col=raincol[mdp$precip.level],
       pch=(18-temp.level))

# Prec on x axis

plot(divp$H~precipitation,data=mdp,
     pch=(3-mdp$temp.level),ylab="Shannon divpersity (H')",cex=1.2)
points(divp$H~precipitation,data=mdp,col=raincol[mdp$precip.level],pch=(18-mdp$temp.level))
points(divp.t$H~precipitation,data=mdp.t,col=transcol,pch=(3-temp.level),cex=1.2)
points(divp.t$H~precipitation,data=mdp.t,col=raincol[mdp$precip.level],
       pch=(18-temp.level))
abline(lmp1,lty=2)
abline(lmp1.t,lty=2,col="grey")

