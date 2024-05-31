## Methylococcaceae ##

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


#Microbiological Data
All_16S <- read.csv(file = "16S_all.csv",header = TRUE)
T0_16S <- All_16S[All_16S$Time==0,]

Mco <- log10(T0_16S$Methylococcaceae)

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

## Oxygen
cor.test(o2,Mco, alternative = "two.sided", method = "pearson")
## Negative correlation - Pearson: r = -0.827667, t = -5.5181, p-value = 7.574e-05

## Methane
cor.test(ch4,Mco, alternative = "two.sided", method = "pearson")
## No correlation - Pearson: r = -0.2773661, t= -1.0802 and p-value = 0.2983

#Temp
cor.test(temp,Mco, alternative = "two.sided", method = "pearson")
## Negative correlation - Pearson: r = -0.5567459, t = -2.5078,p = 0.02509


#There are correlations for oxygen  and temperature, not for methane.


##################
## Linear Model ##
##################

# Testing for interactions
mod.cocc.int <- lm(Mco ~ o2 + ch4 + temp + o2:temp)
Anova(mod.cocc.int) # yes
mod.cocc.int <- lm(Mco ~ o2 + ch4 + temp + o2:ch4 )
Anova(mod.cocc.int) # yes
mod.cocc.int <- lm(Mco ~ o2 + ch4 + temp + ch4:temp)
Anova(mod.cocc.int) # no

#Might be a significant interaction between oxygen and methane/temperature. 

#None
mod.cocc.none <- lm(Mco ~ 1)

#All 
mod.cocc.all <- lm(Mco ~ ch4 + o2 + temp)
mod.cocc.all.int.co <- lm(Mco ~ temp + ch4 * o2)
mod.cocc.all.int.to <- lm(Mco ~ ch4 + temp * o2)

#Two
mod.cocc.co <- lm(Mco ~ ch4 + o2)
mod.cocc.co.int <- lm(Mco ~ ch4 * o2)
mod.cocc.ct <- lm(Mco ~ ch4 + temp)
mod.cocc.to <- lm(Mco ~ o2 + temp)
mod.cocc.to.int <- lm(Mco ~ o2 * temp)

#One
mod.cocc.c <- lm(Mco ~ ch4)
mod.cocc.t <- lm(Mco ~ temp)
mod.cocc.o <- lm(Mco ~ o2)

################
## AICc & BIC ##
################

# AICc 

AICc(mod.cocc.none) #25.45115
AICc(mod.cocc.all) # 16.3453
AICc(mod.cocc.all.int.co) # 2.754088
AICc(mod.cocc.all.int.to) # 5.655261

AICc(mod.cocc.co) # 13.62973
AICc(mod.cocc.co.int) # -2.577116 -- Best model
AICc(mod.cocc.ct) # 20.41954
AICc(mod.cocc.to) # 12.83996
AICc(mod.cocc.to.int) #12.25472

AICc(mod.cocc.c) # 27.24724
AICc(mod.cocc.t) # 22.59184
AICc(mod.cocc.o) # 10.04349

## The model with int & o2*ch4 is the best model. 

# BIC 

BIC(mod.cocc.none) #26.07325
BIC(mod.cocc.all) #14.20824
BIC(mod.cocc.all.int.co) #-1.943713
BIC(mod.cocc.all.int.to) #0.9574604

BIC(mod.cocc.co) #13.08372 
BIC(mod.cocc.co.int) #-4.714172 -- Best model 
BIC(mod.cocc.ct) #19.87353
BIC(mod.cocc.to) #12.29395
BIC(mod.cocc.to.int) #10.11766

BIC(mod.cocc.c) #27.565
BIC(mod.cocc.t) #22.91559
BIC(mod.cocc.o) #10.37798

## The model with o2+ch4 and o2:ch4 is the best model.

# Based upon AIC & BIC, the o2*ch4 w/ o2:ch4 is best. 

###########
## ANOVA ##
###########

summary(mod.cocc.co.int) #p-value: 1.207e-06, adj. r^2: 0.8915 
Anova(mod.cocc.co.int)
#Response: Mco
#           Sum Sq  Df F value    Pr(>F)    
#ch4       0.00333  1  0.1362 0.7185058    
#o2        2.05838  1 84.1992 8.994e-07 ***
#ch4:o2    0.76773  1 31.4044 0.0001154 ***
#Residuals 0.29336 12             


#####################
## Quality control ##
#####################

plot(mod.cocc.co.int) ## Acceptable
olsrr::ols_test_normality(mod.cocc.co.int)
#-----------------------------------------------
#     Test               Statistic       pvalue  
#-----------------------------------------------
#Shapiro-Wilk              0.9209         0.1747 
#Kolmogorov-Smirnov        0.1746         0.6512 
#Cramer-von Mises          4.0515         0.0000 
#Anderson-Darling          0.4879         0.1916 
#-----------------------------------------------

#No abnormalities and the agreement of three normality tests indicates residuals are normal. 

###################
## Plot of Model ##
###################

interactions::sim_slopes(mod.cocc.co.int, pred = o2, modx = ch4, johnson_neyman = FALSE,cond.int = TRUE)
interactions::sim_slopes(mod.cocc.co.int, pred = ch4, modx = o2, johnson_neyman = FALSE)

o2.slopes <- sim_slopes(mod.cocc.co.int, pred = o2, modx = ch4, 
                        modx.values = c(107,200,250,300,350,378), johnson_neyman = FALSE)
ch4.slopes <- sim_slopes(mod.cocc.co.int, pred = ch4, modx = o2, 
                         modx.values = c(70,120,170,220,270,320), johnson_neyman = FALSE)
grid.arrange(plot(o2.slopes),plot(ch4.slopes), nrow=1,ncol=2)


johnson_neyman(mod.cocc.co.int, pred = o2, modx = ch4, alpha = .05)
johnson_neyman(mod.cocc.co.int, pred = ch4, modx = o2, alpha = .05)


probe_interaction(mod.cocc.co.int, pred = o2, modx = ch4, cond.int = TRUE,
                  interval = TRUE,  jnplot = TRUE)
probe_interaction(mod.cocc.co.int, pred = ch4, modx = o2, cond.int = TRUE,
                  interval = TRUE,  jnplot = TRUE)
