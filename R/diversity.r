
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

# ---- OLD below this -----

mtp2 = glm(RR~precipitation+temperature+turf.treat.group,
           family="quasipoisson",data=mdt)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                 7.406e+00  9.925e-02  74.618  < 2e-16 ***
#   precipitation               8.363e-05  1.968e-05   4.249 9.58e-05 ***
#   temperature                -3.251e-02  1.184e-02  -2.746  0.00842 ** 
#   turf.treat.groupTransplant -4.105e-03  3.146e-02  -0.130  0.89672    

mtp3 = glm(RR~precipitation+temperature+turf.treat,
           family="quasipoisson",data=mdt)
# Intercept)       7.465e+00  1.021e-01  73.107  < 2e-16 ***
#   precipitation     8.692e-05  1.921e-05   4.525 3.98e-05 ***
#   temperature      -4.122e-02  1.251e-02  -3.296  0.00185 ** 
#   turf.treatgreen  -3.281e-02  3.453e-02  -0.950  0.34677    
# turf.treatWarmer  4.876e-02  4.178e-02   1.167  0.24896    

mtg = glm(RR~temp.level,family="quasipoisson",data=mdt)
mtgo = mtp1 = glm(RR~origSite.tempLevel,family="quasipoisson",data=mdt)
anova(m0,mtgo,mtg,test="Chisq")
# 1        52     1389.6                       
# 2        51     1294.8  1   94.807  0.04927 *
#   3        51     1298.3  0   -3.517    
# --> original site temp. is more important (same for quasipoisson)

mpg = glm(RR~precip.level,family="quasipoisson",data=mdt)
mpgo = mtp1 = glm(RR~origSite.precipLevel,family="quasipoisson",data=mdt)
anova(m0,mpgo,mpg,test="Chisq")
# Model 1: RR ~ 1
# Model 2: RR ~ origSite.precipLevel
# Model 3: RR ~ precip.level
# Resid. Df Resid. Dev Df Deviance  Pr(>Chi)    
# 1        52     1389.6                          
# 2        51     1068.8  1   320.82 6.891e-05 ***
#   3        51     1078.6  0    -9.78              
# --> original site precipitation is more important

mpgtl = glm(RR~precipitation+temperature,family="quasipoisson",data=mdt)
mpgto = mtp1 = glm(RR~origSite.precipLevel+origSite.tempLevel,family="quasipoisson",data=mdt)
mpgt = glm(RR~precip.level+temp.level,family="quasipoisson",data=mdt)
anova(m0,mpgt,mpgto,mpgtl,test="Chisq")
# 1        52    1389.64                          
# 2        50     988.91  2   400.73 2.885e-05 ***
#   3        50     983.47  0     5.45              
# 4        50     913.67  0    69.80               

# --> Destination site both is more imoirtant

mpgtl.l = lm(RR~precipitation+temperature,data=mdt)
summary(mpgto.l)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   1413.94      99.65  14.189  < 2e-16 ***
#   precip.level    88.51      22.01   4.021 0.000196 ***
#   temp.level    -103.25      47.58  -2.170 0.034788 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 167.9 on 50 degrees of freedom
# Multiple R-squared:  0.2957,	Adjusted R-squared:  0.2675 
# F-statistic: 10.49 on 2 and 50 DF,  p-value: 0.0001566

# ---- H' instead of RR


lmp1 = lm(div.c$H~precipitation,data=mdt.c)
summary(lmp1)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   6.4916508  0.1107367  58.622   <2e-16 ***
#   precipitation 0.0001558  0.0000589   2.645   0.0148 *  

# Similar to RR but slightly worse fit

lmp1.t = lm(div.t$H~precipitation,data=mdt.t)
summary(lmp1.t)
# (Intercept)   6.469e+00  1.279e-01  50.590   <2e-16 ***
#   precipitation 1.493e-04  6.296e-05   2.371   0.0251 *  

# Similar to RR but slightly worse fit

# ----

lmt1 = lm(div.c$H~temperature,data=mdt.c) 
summary(lmt1)
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  7.43359    0.27313  27.217   <2e-16 ***
#   temperature -0.08912    0.03537  -2.519   0.0195 *  

# A linear model is also fine...
# Now make separate one for transplants

lmt1.t = lm(div.t$H~temperature,data=mdt.t) 
summary(lmt1.t)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  7.07861    0.34074  20.774   <2e-16 ***
#   temperature -0.03996    0.04122  -0.969    0.341    
# 
# Also not significant like for RR

# ----

m0 = lm(div$H~1,data=mdt) #AIC Inf, res dev 1390
mtgl = lm(div$H~temperature+precipitation,data=mdt)
mtg = lm(div$H~temp.level+precip.level,data=mdt)
mtgo = lm(div$H~origSite.tempLevel+origSite.precipLevel,data=mdt)
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
plot(div.c$H~temperature,data=mdt.c,
     pch=(3-mdt.c$temp.level),ylab="Shannon diversity (H')",cex=1.2)
points(div.c$H~temperature,data=mdt.c,col=raincol[mdt.c$precip.level],pch=(18-mdt.c$temp.level))
abline(lmt1,lty=2)
points(div.t$H~temperature,data=mdt.t,col=transcol,pch=(3-temp.level),cex=1.2)
points(div.t$H~temperature,data=mdt.t,col=raincol[mdt.t$precip.level],
       pch=(18-temp.level))

# Prec on x axis

plot(div.c$H~precipitation,data=mdt.c,
     pch=(3-mdt.c$temp.level),ylab="Shannon diversity (H')",cex=1.2)
points(div.c$H~precipitation,data=mdt.c,col=raincol[mdt.c$precip.level],pch=(18-mdt.c$temp.level))
points(div.t$H~precipitation,data=mdt.t,col=transcol,pch=(3-temp.level),cex=1.2)
points(div.t$H~precipitation,data=mdt.t,col=raincol[mdt.t$precip.level],
       pch=(18-temp.level))
abline(lmp1,lty=2)
abline(lmp1.t,lty=2,col="grey")

# ---- Old stuff ----

printANOVA(mdt.c[,c("precip.level","temp.level")],div.c,a=0.05)
#Templevel 2 (out of 2) has lower RR, H, J, but weak p (near 0.05)

printVS(div.c,mdt.c[,c("temperature","precipitation")],a=0.05)
#precip.level, temp (neg.) and precipiation all correlated with RR, H' and J' (not prec as such for J)

m0 = glm(div.c$Rarefied.richness~1,family="poisson",data=mdt.c) #AIC Inf, res dev 670
m1 = glm(div.c$Rarefied.richness~temperature,family="poisson",data=mdt.c) #AIC Inf, res dev 670
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  7.653038   0.029827  256.58   <2e-16 ***
#   temperature -0.046236   0.003906  -11.84   <2e-16 ***
# 
# Null deviance: 669.85  on 23  degrees of freedom
# Residual deviance: 529.51  on 22  degrees of freedom

library(AER)
dispersiontest(m1)# > p = 2E-5 -> use quasipoisson

m0 = glm(div.c$Rarefied.richness~1,family="quasipoisson",data=mdt.c) #AIC Inf, res dev 670
mt1 = glm(div.c$Rarefied.richness~temperature,family="quasipoisson",data=mdt.c) 
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  7.65304    0.14442  52.992   <2e-16 ***
#   temperature -0.04624    0.01891  -2.445    0.023 *  
#   Null deviance: 669.85  on 23  degrees of freedom
# Residual deviance: 529.51  on 22  degrees of freedom (AIC inf)

anova(m0,mt1,test="Chisq") # -> m1 best i.e. temperature explains richness!
# Model 1: div.c$Rarefied.richness ~ 1
# Model 2: div.c$Rarefied.richness ~ temperature
# Resid. Df Resid. Dev Df Deviance Pr(>Chi)  
# 1        23     669.85                       
# 2        22     529.51  1   140.34  0.01442 *

lmt1 = lm(div.c$Rarefied.richness~temperature,data=mdt.c) 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  2007.42     216.92   9.254 4.85e-09 ***
#   temperature   -68.65      28.09  -2.444    0.023 *  
# Residual standard error: 186.8 on 22 degrees of freedom
# Multiple R-squared:  0.2135,	Adjusted R-squared:  0.1777 
# F-statistic: 5.972 on 1 and 22 DF,  p-value: 0.02302

# A linear model is also fine...
# Now make separate one for transplants

lmt1.t = lm(div.t$Rarefied.richness~temperature,data=mdt.t) 
summary(lmt1.t)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  1743.66     242.10   7.202 9.56e-08 ***
#   temperature   -32.77      29.29  -1.119    0.273

# ---> It is not significant

# ----

mp1 = glm(div.c$Rarefied.richness~precipitation,family="quasipoisson",data=mdt.c) 
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   7.143e+00  5.776e-02 123.666  < 2e-16 ***
#   precipitation 9.279e-05  2.980e-05   3.114  0.00506 ** 
# Null deviance: 669.85  on 23  degrees of freedom
# Residual deviance: 467.68  on 22  degrees of freedom

lmp1 = lm(div.c$Rarefied.richness~precipitation,data=mdt.c) 
summary(lmp1)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   1.249e+03  8.333e+01  14.990 4.98e-13 ***
#   precipitation 1.393e-01  4.432e-02   3.142  0.00474 ** 
# Residual standard error: 175 on 22 degrees of freedom
# Multiple R-squared:  0.3097,	Adjusted R-squared:  0.2784 
# F-statistic: 9.872 on 1 and 22 DF,  p-value: 0.004735

# A linear model is also fine and very similar

lmp1.t = lm(div.t$Rarefied.richness~precipitation,data=mdt.t)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   1.256e+03  8.959e+01  14.018 6.55e-14 ***
#   precipitation 1.161e-01  4.411e-02   2.632   0.0139 *  

lmboth = lm(div.c$Rarefied.richness~precipitation+temperature,data=mdt.c)
summary(lmboth)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   1782.20168  180.25335   9.887 2.36e-09 ***
#   precipitation    0.14232    0.03717   3.829 0.000978 ***
#   temperature    -70.82417   22.07456  -3.208 0.004220 ** 

# Residual standard error: 146.8 on 21 degrees of freedom
# Multiple R-squared:  0.5368,	Adjusted R-squared:  0.4927 
# F-statistic: 12.17 on 2 and 21 DF,  p-value: 0.0003094

lmboth.t = lm(div.t$Rarefied.richness~precipitation+temperature,data=mdt.t)
summary(lmboth.t)

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   1465.89850  248.07573   5.909 3.11e-06 ***
#   precipitation    0.11114    0.04458   2.493   0.0194 *  
#   temperature    -24.55013   27.01549  -0.909   0.3718    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 174.3 on 26 degrees of freedom
# Multiple R-squared:  0.2286,	Adjusted R-squared:  0.1693 
# F-statistic: 3.853 on 2 and 26 DF,  p-value: 0.03422




# --- Plots ----

# Temp on x axis
plot(div.c$Rarefied.richness~temperature,data=mdt.c,
     pch=(3-mdt.c$temp.level),ylab="Rarefied OTU richness",cex=1.2)
points(div.c$Rarefied.richness~temperature,data=mdt.c,col=raincol[mdt.c$precip.level],pch=(18-mdt.c$temp.level))
abline(lmt1,lty=2)
points(div.t$Rarefied.richness~temperature,data=mdt.t,col=transcol,pch=(3-temp.level),cex=1.2)
points(div.t$Rarefied.richness~temperature,data=mdt.t,col=raincol[mdt.t$precip.level],
       pch=(18-temp.level))


# Prec on x axis


plot(div.c$Rarefied.richness~precipitation,data=mdt.c,
     pch=(3-mdt.c$temp.level),ylab="Rarefied OTU richness",cex=1.2)
points(div.c$Rarefied.richness~precipitation,data=mdt.c,col=raincol[mdt.c$precip.level],pch=(18-mdt.c$temp.level))
points(div.t$Rarefied.richness~precipitation,data=mdt.t,col=transcol,pch=(3-temp.level),cex=1.2)
points(div.t$Rarefied.richness~precipitation,data=mdt.t,col=raincol[mdt.t$precip.level],
       pch=(18-temp.level))
abline(lmp1,lty=2)
abline(lmp1.t,lty=2,col="grey")


# ---- MULTIDIMENSIONAL MODEL WITH PREC + TEMP --------

# Result: both are significant but questionable whether this is ok with replicates!

mtp1 = glm(div.c$Rarefied.richness~precipitation+temperature,family="quasipoisson",data=mdt.c)
# (Intercept)    7.505e+00  1.223e-01  61.365  < 2e-16 ***
#   precipitation  9.570e-05  2.515e-05   3.805  0.00104 ** 
#   temperature   -4.853e-02  1.518e-02  -3.197  0.00433 ** 
# Null deviance: 669.85  on 23  degrees of freedom
# Residual deviance: 315.57  on 21  degrees of freedom

anova(m0,mp1,mtp1,test="Chisq") # --> Both temp. and precipitation!
# Analysis of Deviance Table
# 
# Model 1: div.c$Rarefied.richness ~ 1
# Model 2: div.c$Rarefied.richness ~ precipitation
# Model 3: div.c$Rarefied.richness ~ precipitation + temperature
# Resid. Df Resid. Dev Df Deviance  Pr(>Chi)    
# 1        23     669.85                          
# 2        22     467.68  1   202.17 0.0002252 ***
#   3        21     315.57  1   152.11 0.0013753 ** 


mtp2 = glm(div.c$Rarefied.richness~precipitation*temperature,family="quasipoisson",data=mdt.c) 
# (Intercept)                7.272e+00  2.791e-01  26.053   <2e-16 ***
#   precipitation              2.416e-04  1.586e-04   1.523    0.143    
# temperature               -1.789e-02  3.623e-02  -0.494    0.627    
# precipitation:temperature -1.914e-05  2.055e-05  -0.931    0.363    
# Null deviance: 669.85  on 23  degrees of freedom
# Residual deviance: 302.60  on 20  degrees of freedom

# A: Diversity (alt. and precipitation only, first alone then together, quasipoisson) ----

printANOVA(mdt[,c("site","precip.level","temp.level","turf.treat")],div,a=0.05)

printVS(div,mdt[,c("temperature","precipitation","origSite.precipLevel","origSite.tempLevel")],a=0.05)


