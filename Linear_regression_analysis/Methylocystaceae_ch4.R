## Methylocystaceae ##

library(fBasics)
library(car)
library(olsrr)
library(grid)
library(gridExtra)
library(interactions)
library(survey)
library(jtools)
library(wiqid)

# This file contains the correlation and linear modeling code for the Jordan Lake T0 
# Methylocystaceae relative abundance data

#Microbiological Data
All_16S <- read.csv(file = "16S_all.csv",header = TRUE)
T0_16S <- All_16S[All_16S$Time==0,]

Mcy <- log10(T0_16S$Methylocystaceae)

geochem <- read.csv(file = "T0_geochem.csv" ,sep = ",", header = TRUE)
row.names(geochem) <- geochem[,1] #set row names
geochem <- geochem[,-1] #remove first column
geochem <- geochem[,-5] #remove empty column at end


o2 <- geochem$Modeled_Initial_Oxygen
ch4 <- geochem$Modeled_Initial_Methane
temp <- c(20.01,20.024,20.099,27.287,27.308,27.253,28.339,28.544,28.575,27.794,
          27.93,20.268,20.153,20.127,18.287,18.276)
season <- c("O","O","O","J","J","J","J","J","J","J","J","O","O","O","O","O")

#################
## Correlation ##
#################

# Oxygen 
cor.test(o2,Mcy, alternative = "two.sided", method = "pearson")
## Negative correlation - Pearson: r = -0.4298598, t = -1.7814, and p-value = 0.09655

# Methane
cor.test(ch4,Mcy, alternative = "two.sided", method = "pearson")
## No correlation - Pearson: r = 0.199876, t = 0.76327, and p-value = 0.458

#Temperature 
cor.test(temp, Mcy, alternative = "two.sided", method = "pearson")
## Negative correlation - pearson: r = -0.9140476, t = -8.4319,  p-value = 7.386e-07

#There appears to be a negative correlation between Methylocystaceae and oxygen as well as 
#temperature. There does not appear to be any correlation between Methylocystaceae and methane.

##################
## Linear Model ##
##################

# Testing for interactions
mod.cyst.int <- lm(Mcy ~  o2* ch4 * temp)
Anova(mod.cyst.int) #Seemingly no significant interactions.
mod.cyst.int <- lm(Mcy ~  o2* temp)
Anova(mod.cyst.int) #No significant interaction
mod.cyst.int <- lm(Mcy ~  ch4* temp)
Anova(mod.cyst.int) #No significant interaction. 
mod.cyst.int <- lm(Mcy ~  o2 * ch4)
Anova(mod.cyst.int) #No significant interaction. 

#None
mod.cyst.none <- lm(Mcy ~ 1)

#All 
mod.cyst.all <- lm(Mcy ~ o2 + ch4 + temp)

#Two
mod.cyst.ct <- lm(Mcy ~ ch4 + temp)
mod.cyst.co <- lm(Mcy ~ o2 + ch4)
mod.cyst.to <- lm(Mcy ~ o2 + temp)

#One
mod.cyst.c <- lm(Mcy ~ ch4)
mod.cyst.t <- lm(Mcy ~ temp)
mod.cyst.o <- lm(Mcy ~ o2)


################
## AICc & BIC ##
################

#AICc

AICc(mod.cyst.none) #27.96169
AICc(mod.cyst.all) #7.468906 

AICc(mod.cyst.ct) #5.344672 
AICc(mod.cyst.co) #28.89513
AICc(mod.cyst.to) #5.159253 

AICc(mod.cyst.c) #30.38629
AICc(mod.cyst.t) #2.162751 -- best model
AICc(mod.cyst.o) #27.76987

#AICc indicates that the best model is temp. 


#BIC 

BIC(mod.cyst.none) #28.58379
BIC(mod.cyst.all) #5.331849

BIC(mod.cyst.ct) #4.798663 
BIC(mod.cyst.co) #28.34912
BIC(mod.cyst.to) #4.613244 

BIC(mod.cyst.c) #30.70405
BIC(mod.cyst.t) #2.480518 -- best model
BIC(mod.cyst.o) #28.08763

#BIC indicates that the best model is temp, with temp+ch4 and temp+o2 close to the best family. 

#The temperature model comes out as the best. ct and o2 interactions will also be examined. 

###########
## ANOVA ##
###########

summary(mod.cyst.ct) #p-value: 6.686e-06, adj. r^2 = 0.8155 
Anova(mod.cyst.ct) #Temperature is ***
#Response: Mcy
#            Sum Sq  Df  F value    Pr(>F)  
#ch4         0.0182  1  0.3745    0.5511    
#temp        3.1634  1 65.0479 2.046e-06 ***
#Residuals   0.6322 13          

summary(mod.cyst.to) #p-value: 6.201e-06, adj. r^2 = 0.8176
Anova(mod.cyst.to) #Temperature is extremely significant (***)
#Response: Mcy
#            Sum Sq  Df  F value    Pr(>F)    
#o2          0.02550  1  0.5304   0.4793    
#temp        2.59812  1 54.0464 5.57e-06 ***
#Residuals   0.62493 13       
 
summary(mod.cyst.t) #p-value: 7.386e-07, adj. r^2 = 0.8237
Anova(mod.cyst.t) #Temperature is highly significant (***)
#Response: Mcy
#            Sum Sq  Df  F value    Pr(>F)    
#temp         3.3032  1  71.098 7.386e-07 ***
#Residuals    0.6504 14           


#temperature is the best.  

#####################
## Quality control ##
#####################

plot(mod.cyst.t) ## qq is fine
olsrr::ols_test_normality(mod.cyst.t)

#-----------------------------------------------
#     Test              Statistic         pvalue  
#-----------------------------------------------
#Shapiro-Wilk              0.9133         0.1314 
#Kolmogorov-Smirnov         0.2           0.4835 
#Cramer-von Mises          3.4243         0.0000 
#Anderson-Darling          0.5708         0.1177 
#-----------------------------------------------

#The temperature model is normal. 


###################
## Plot of Model ##
###################

summary(mod.cyst.t)

pred.frame.temp<-data.frame(temp=seq(from=18,to=30,length=100)) #range of values to calculate CI for 

cyst.pred.ci<-predict(mod.cyst.t,int="confidence",newdata=pred.frame.temp) #actual CI values
cyst.pred.p<-predict(mod.cyst.t,int="prediction",newdata=pred.frame.temp) #predicted data

plot(temp,Mcy,pch=19,las=1, xlab = "Temperature (degrees Celsius)", 
     ylab = "Relative Abundance") + #,col = "orange2"
  text(x = 20.5, y = -3.3, labels = "P value < 0.001") +
  text(x = 20.5, y = -3.5, labels = "Adj. R^2 = 0.6543") #Basic plot of data
abline(mod.cyst.t) #adding the slope
matlines(pred.frame.temp,cyst.pred.ci,col=c("black","red","red"),lwd=3,lty=c(1,2,2)) #adding the CI lines
matlines(pred.frame.temp,cyst.pred.p,col=c("black","red","red"),lwd=3,lty=c(1,3,3))#adding the PI lines
