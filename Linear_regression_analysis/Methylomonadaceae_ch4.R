## Methylomonadaceae ##

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
# Methylococcaceae relative abundance data

All_16S <- read.csv(file = "16S_all.csv",header = TRUE)
T0_16S <- All_16S[All_16S$Time==0,]

Mm <- log10(T0_16S$Methylomonadaceae)

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
cor.test(o2,Mm, alternative = "two.sided", method = "pearson")
# Negative correlation - Pearson: r = -0.5883873, t = -2.7227, p-value = 0.0165

# Methane 

cor.test(ch4,Mm, alternative = "two.sided", method = "pearson")
## No correlation - Pearson: r = -0.3447275, t = -1.3741, and p-value = 0.191

#Oxygen significantly correlates with relative abundance, methane does not. 

##Temp
cor.test(temp,Mm, alternative = "two.sided", method = "pearson")
## No  correlation - Pearson: r = -0.1188546, t = -0.44789, p-value = 0.6611


##################
## Linear Model ##
##################

# Testing for interactions 
mod.mona.int <- lm(Mm ~ o2 + ch4 + temp + o2:ch4 + o2:temp + ch4:temp)
Anova(mod.mona.int)
# Further breakdown needed

mod.mona.to.int <- lm(Mm ~ o2 + temp * ch4)
Anova(mod.mona.to.int) #No
mod.mona.to.int <- lm(Mm ~ o2 * temp + ch4)
Anova(mod.mona.to.int) #yes
mod.mona.co.int <- lm(Mm ~ o2 * ch4 + temp)
Anova(mod.mona.co.int) #yes
# o2:ch4 and o2:temp interactions included

#None 
mod.mona.none <- lm(Mm ~ 1)

#All 
mod.mona.all.int.to <- lm(Mm ~ ch4 + o2 + temp + o2:temp)
mod.mona.all.int.co <- lm(Mm ~ ch4 + o2 + temp + o2:ch4)
mod.mona.all <- lm(Mm ~ o2 + ch4 + temp)

#Two
mod.mona.ct <- lm(Mm ~ ch4 + temp)
mod.mona.co <- lm(Mm ~ o2 + ch4)
mod.mona.to <- lm(Mm ~ o2 + temp)
mod.mona.co.int <- lm(Mm ~ o2 + ch4 + o2:ch4)
mod.mona.to.int <- lm(Mm ~ o2 + temp + o2:temp)

#One 
mod.mona.c <- lm(Mm ~ ch4)
mod.mona.t <- lm(Mm ~ temp)
mod.mona.o <- lm(Mm ~ o2)


################
## AICc & BIC ##
################

# AICc

AICc(mod.mona.none)#31.76549
AICc(mod.mona.all.int.co) #22.4862 - best
AICc(mod.mona.all.int.to) #27.37848
AICc(mod.mona.all) #34.502

AICc(mod.mona.ct) #35.4769
AICc(mod.mona.co) #30.89446
AICc(mod.mona.to) #30.20622
AICc(mod.mona.to.int) #28.19423
AICc(mod.mona.co.int) #26.24285

AICc(mod.mona.c) #32.81821
AICc(mod.mona.t) #34.61478
AICc(mod.mona.o) #28.04316

## Points to all + o2:ch4 

# BIC

BIC(mod.mona.none)#32.38759
BIC(mod.mona.all.int) #19.33666 -- family of best models 
BIC(mod.mona.all.int.co) #17.7884  -- best model
BIC(mod.mona.all.int.to) #22.68068
BIC(mod.mona.all) #32.36495

BIC(mod.mona.ct) #34.93089
BIC(mod.mona.co) #30.34845
BIC(mod.mona.to) #29.66021
BIC(mod.mona.to.int) #26.05718
BIC(mod.mona.co.int) #24.10579  

BIC(mod.mona.c) #33.13597
BIC(mod.mona.t) #34.93254
BIC(mod.mona.o) #28.36093

## According to BIC, the all + o2:ch4 interaction model is the best, with all + both interactions 
## model being the next best model for the Methylomonadaceae relative abundance. 

# Both metrics converge on the all.int and all + o2:ch4 int models as the best family of models. 

###########
## ANOVA ##
###########

summary(mod.mona.all.int) #p-value: 0.002028, adj. r^2: 0.721
Anova(mod.mona.all.int) #Interactions are insignificant
#          Sum Sq Df F value   Pr(>F)   
#ch4       0.64054  1  6.8681 0.025566 * 
#o2        1.17929  1 12.6449 0.005216 **
#temp      0.77015  1  8.2579 0.016565 * 
#o2:temp   0.07417  1  0.7952 0.393471   
#ch4:o2    0.43426  1  4.6563 0.056313 .
#Residuals 0.93262 10       


summary(mod.mona.all.int.co) #p-value: 0.0007887, adj. r^2: 0.7262
Anova(mod.mona.all.int.co) #Temp is very significant (**), o2 is significant (*) and o2:ch4 is (**)
#Response: log(mona + 5e-05)
#          Sum Sq Df F value   Pr(>F)   
#ch4       0.01265  1  0.1382 0.7170957    
#o2        1.17929  1 12.8847 0.0042471 ** 
#temp      0.77015  1  8.4145 0.0144252 *  
#ch4:o2    1.97072  1 21.5317 0.0007164 ***
#Residuals 1.00679 11                   

#I will go with the all.co.int model, as pulls out the significant interaction term and has the higher
#adj. r^2.

#####################
## Quality control ##
#####################

plot(mod.mona.all.int.co) ## Acceptable
olsrr::ols_test_normality(mod.mona.all.int.co)
#-----------------------------------------------
#   Test                  Statistic       pvalue  
#-----------------------------------------------
#Shapiro-Wilk              0.954          0.5557 
#Kolmogorov-Smirnov        0.1213         0.9500 
#Cramer-von Mises          2.9668         0.0000 
#Anderson-Darling          0.2784         0.6006 
#-----------------------------------------------

#No abnormalities and the agreement of three normality tests indicates the residuals are normal.

###################
## Plot of Model ##
###################

#Model temperature

residmod <- lm(Mm ~ ch4*o2)$resid
mod_plot_temp <- lm(residmod ~ temp)
plot(temp, residmod)
pred.frame.temp <-data.frame(temp=seq(from=18,to=29,length=100)) #range of values to calculate CI for 
mona.pred.ci<-predict(mod_plot_temp,int="confidence",newdata=pred.frame.temp) #actual CI values
mona.pred.p<-predict(mod_plot_temp,int="prediction",newdata=pred.frame.temp) #predicted data
plot(temp,residmod,pch=19,las=1, ylim = c(-0.7,0.7), xlab = "Temperature (degrees C)",
     ylab = "Residuals unexplained by O2 or CH4")  #main = "Linear Model - Methylococcaceae",
#Plot of data that temp explains, not o2 &/or ch4
abline(mod_plot_temp) #adding the slope
matlines(pred.frame.temp,mona.pred.ci,col=c("black","red","red"),lwd=3,lty=c(1,2,2)) #adding the CI lines
matlines(pred.frame.temp,mona.pred.p,col=c("black","red","red"),lwd=3,lty=c(1,3,3))#adding the PI lines
text(x = 25, y = -0.5, labels = "P = 0.006")

#Model Interactions

interactions::sim_slopes(mod.mona.all.int.co, pred = o2, modx = ch4, johnson_neyman = FALSE)

johnson_neyman(mod.mona.all.int.co, pred = o2, modx = ch4, alpha = .05)
johnson_neyman(mod.mona.all.int.co, pred = ch4, modx = o2, alpha = .05)


o2.slopes <- sim_slopes(mod.mona.all.int.co, pred = o2, modx = ch4, 
                        modx.values = c(107,200,250,300,350,378), johnson_neyman = FALSE)
ch4.slopes <- sim_slopes(mod.mona.all.int.co, pred = ch4, modx = o2, 
                        modx.values = c(70,120,170,220,270,320), johnson_neyman = FALSE)
grid.arrange(plot(o2.slopes),plot(ch4.slopes), nrow=1,ncol=2)
