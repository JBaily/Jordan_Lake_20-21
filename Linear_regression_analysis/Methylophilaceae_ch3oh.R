## Methylophilaceae ##

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
# Methylophilaceae relative abundance data

#Microbiological Data
All_16S <- read.csv(file = "16S_all.csv",header = TRUE)
T0_16S <- All_16S[All_16S$Time==0,]

Mp <- log10(T0_16S$Methylophilaceae)

geochem <- read.csv(file = "T0_geochem.csv" ,sep = ",", header = TRUE)
row.names(geochem) <- geochem[,1] #set row names
geochem <- geochem[,-1] #remove first column
geochem <- geochem[,-5] #remove empty column at end


o2 <- geochem$Modeled_Initial_Oxygen
ch4 <- geochem$Modeled_Initial_Methane
temp <- c(20.01,20.024,20.099,27.287,27.308,27.253,28.339,28.544,28.575,27.794,
          27.93,20.268,20.153,20.127,18.287,18.276)
season <- c("O","O","O","J","J","J","J","J","J","J","J","O","O","O","O","O")

#palette_check(col2hex(c("black","darkgray","orange","hotpink1",
#                        "purple4","blue","green")), plot = TRUE)

#################
## Correlation ##
#################

# Oxygen
cor.test(o2,Mp, alternative = "two.sided", method = "pearson")
## No correlation - p-value = 0.1587, t = -1.489, r = -0.3697468

# Methane
cor.test(ch4,Mp, alternative = "two.sided", method = "pearson")
## No correlation - p-value = 0.3245, t = 1.0211, r = 0.2632744

# Temp 
cor.test(temp,Mp, alternative = "two.sided", method = "pearson")
## No correlation - p-value = 0.3667, t = -0.9329, r = -0.241921

# No correlations

##################
## Linear Model ##
##################

# Checking for interactions
mod.phil.int <- lm(Mp ~ o2 * ch4 * temp)
Anova(mod.phil.int) #Seems to be an interaction between oxygen and temperature
mod.phil.int <- lm(Mp ~ o2 * ch4)
Anova(mod.phil.int) #No
mod.phil.int <- lm(Mp ~ ch4 * temp)
Anova(mod.phil.int) #Yes
mod.phil.int <- lm(Mp ~ o2 * temp)
Anova(mod.phil.int) #Yes
#Potential interactions between ch4:temp and o2:temp

#None
mod.phil.none <- lm(Mp ~ 1)

#All
mod.phil.all <- lm(Mp ~ o2 + ch4 + temp)
mod.phil.all.int.to <- lm(Mp ~ ch4 + o2 * temp)
mod.phil.all.int.ct <- lm(Mp ~ ch4 * temp + o2)

#Two
mod.phil.ct <- lm(Mp ~ ch4 + temp)
mod.phil.co <- lm(Mp ~ o2 + ch4)
mod.phil.to <- lm(Mp ~ o2 + temp)
mod.phil.to.int <- lm(Mp ~ o2 * temp)
mod.phil.ct.int <- lm(Mp ~ temp * ch4)

#One 
mod.phil.c <- lm(Mp ~ ch4)
mod.phil.t <- lm(Mp ~ temp)
mod.phil.o <- lm(Mp ~ o2)


###############
## AIC & BIC ##
###############

# AICc

AICc(mod.phil.none) #-5.138527
AICc(mod.phil.all) #-0.5640082
AICc(mod.phil.all.int.to) #-18.74452 
AICc(mod.phil.all.int.ct) #0.01315979

AICc(mod.phil.ct) #-0.09793016
AICc(mod.phil.co) #-3.908485
AICc(mod.phil.to) #-0.8211958
AICc(mod.phil.to.int) #-23.96945 -- best model 
AICc(mod.phil.ct.int) #-2.698796

AICc(mod.phil.c) #-3.210927
AICc(mod.phil.t) #-3.026536
AICc(mod.phil.o) #-4.413728

## The o2:temp interaction model is the best.

# BIC

BIC(mod.phil.none) #-4.516426
BIC(mod.phil.all) #-2.701065
BIC(mod.phil.all.int.to) #-23.44232
BIC(mod.phil.all.int.ct) #-4.684641

BIC(mod.phil.ct) #-0.6439389
BIC(mod.phil.co) #-4.454494
BIC(mod.phil.to) #-1.367205
BIC(mod.phil.to.int) #-26.1065 -- best model
BIC(mod.phil.ct.int) #-4.835853

BIC(mod.phil.c) #-2.893161
BIC(mod.phil.t) #-2.70877
BIC(mod.phil.o) #-4.095962

## According to BIC, the o2:temp model is the best.

## The results indicate the o2:temp model is the best. The all + int model will be examined as 
## well for comparison as it is the second best, though not in the family of best. 

###########
## ANOVA ##
###########

summary(mod.phil.all.int.to) #p-value: 0.0001869, Adjusted R-squared:  0.7911 
Anova(mod.phil.all.int.to)
#Response: Mp # -- extremely significant interaction o2 & interaction. very sig oxygen 
#            Sum Sq  Df  F value     Pr(>F)    
#ch4       0.000520  1  0.0748  0.789561    
#o2        0.117230  1 16.8514  0.001744 ** 
#temp      0.021881  1  3.1453  0.103796    
#o2:temp   0.256168  1 36.8232 8.096e-05 ***
#Residuals 0.076524 11                  

summary(mod.phil.to.int) #p-value: 3.678e-05, Adjusted R-squared:  0.8072 
Anova(mod.phil.to.int) # -- extremely significant interaction (***), o2 barely significant (*), temp insig
#Response: Mp
#           Sum Sq Df F value   Pr(>F)    
#o2        0.04023  1  6.2667  0.02774 *  
#temp      0.00118  1  0.1837  0.67577    
#o2:temp   0.35299  1 54.9804 8.11e-06 ***
#Residuals 0.07704 12          

#Proceeding with the the mod.phil.to.int model. 

#####################
## Quality control ##
#####################

plot(mod.phil.to.int) ## Acceptable
olsrr::ols_test_normality(mod.phil.to.int)
#-----------------------------------------------
#   Test                  Statistic       pvalue  
#-----------------------------------------------
#Shapiro-Wilk              0.9273         0.2212 
#Kolmogorov-Smirnov        0.1126         0.9732 
#Cramer-von Mises          4.6446         0.0000 
#Anderson-Darling          0.3361         0.4609 
#-----------------------------------------------

#No abnormalities and the agreement of three normality tests indicates the residuals are normal. 

###################
## Plot of Model ##
###################

interactions::sim_slopes(mod.phil.to.int, pred = o2, modx = temp, johnson_neyman = FALSE)

johnson_neyman(mod.phil.to.int, pred = o2, modx = temp, alpha = .05)
johnson_neyman(mod.phil.to.int, pred = temp, modx = o2, alpha = .05) 

sim_slopes(mod.phil.to.int, pred = o2, modx = temp, 
           modx.values = c(19.6555, 22, 24, 27.87143), johnson_neyman = FALSE)
sim_slopes(mod.phil.to.int, pred = o2, modx = temp, 
           modx.values = c(min(temp),max(temp)), johnson_neyman = FALSE)

sim1 <- plot(sim_slopes(mod.phil.to.int, pred = o2, modx = temp, 
                modx.values = c(19.6555, 22, 24, 27.87143), johnson_neyman = FALSE))
sim2 <- plot(sim_slopes(mod.phil.to.int, pred = temp, modx = o2, 
                modx.values = c(70,120,170,220,270,320), johnson_neyman = FALSE))
grid.arrange(sim1,sim2, nrow=1,ncol=2)


sim_slopes(mod.phil.to.int, pred = temp, modx = o2, 
           modx.values = c(67,120,170,220,270,320), johnson_neyman = FALSE)
