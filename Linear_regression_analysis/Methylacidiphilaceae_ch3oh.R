## Methylacidiphilaceae ##

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
# Methylacidiphilaceae relative abundance data


#Microbiological Data
All_16S <- read.csv(file = "16S_all.csv",header = TRUE)
T0_16S <- All_16S[All_16S$Time==0,]

Ma <- log10(T0_16S$Methylacidiphilaceae)

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
cor.test(o2,Ma, alternative = "two.sided", method = "pearson")
## No correlation - Pearson: r = 0.3914537, t = 1.5917 and p-value = 0.1338

# Methane 
cor.test(ch4,Ma, alternative = "two.sided", method = "pearson")
## No correlation - Pearson: p-value = 0.47, t = -0.74267, r = -0.194689

#Temp 
cor.test(temp, Ma, alternative = "two.sided", method = "pearson")
## Positive correlation - pearson: r = 0.9017745, t = 7.8068, p-value = 1.818e-06

#Positive correlation with Temperature, none with others. 

##################
## Linear Model ##
##################

#Testing for interactions
mod.acid.int <- lm(Ma ~ o2 * ch4 * temp)
Anova(mod.acid.int)
mod.acid.int <- lm(Ma ~ o2 * ch4) #No
Anova(mod.acid.int)
mod.acid.int <- lm(Ma ~ ch4 * temp) #No
Anova(mod.acid.int)
mod.acid.int <- lm(Ma ~ o2 * temp) # yes
Anova(mod.acid.int)

#temp is the only thing that comes out as significant. 
#o2:temp is significant 

#All

mod.acid.none <- lm(Ma ~ 1)
mod.acid.all <- lm(Ma ~ o2 + ch4 + temp)
mod.acid.all.int.to <- lm(Ma ~ o2 + ch4 + temp + o2:temp)

#Two 
mod.acid.co <- lm(Ma ~ o2 + ch4)
mod.acid.ct <- lm(Ma ~ ch4 + temp)
mod.acid.to <- lm(Ma ~ o2 + temp)
mod.acid.to.int <- lm(Ma ~ o2 * temp)

#One 
mod.acid.c <- lm(Ma ~ ch4)
mod.acid.t <- lm(Ma ~ temp)
mod.acid.o <- lm(Ma ~ o2)


################
## AICc & BIC ##
################

#AICc

AICc(mod.acid.none) #2.334835
AICc(mod.acid.all) #-17.45949  
AICc(mod.acid.all.int.to) #-15.74845 

#Two 
AICc(mod.acid.co) #4.231961
AICc(mod.acid.ct) #-18.21524 
AICc(mod.acid.to) #-19.0372 -- Family of best models 
AICc(mod.acid.to.int) #-20.83436 -- Family of best models 

#One 
AICc(mod.acid.c) #4.793504
AICc(mod.acid.t) #-21.43147 -- Best
AICc(mod.acid.o) #2.750426

#Temperature is the unifying factor amongst the best models. Temperature only is best. 
# to.int and to are in family of best models. 


#BIC

BIC(mod.acid.none) #2.956935
BIC(mod.acid.all) #-19.59655 
BIC(mod.acid.all.int.to) #-20.44625

#Two 
BIC(mod.acid.co) #3.685953
BIC(mod.acid.ct) #-18.76125 
BIC(mod.acid.to) #-19.58321 
BIC(mod.acid.to.int) #-22.97141 -- Best Model

#One 
BIC(mod.acid.c) #5.11127
BIC(mod.acid.t) #-21.11371 -- Family of Best Models
BIC(mod.acid.o) #3.068192

#BIC also pulls temperature out as the most important. to interaction model is in family of best. 

#Both metrics converge on same two models-- temperature only and to.int. AICc also suggest to only.
#Neither agrees on best model. All three will be further analyzed.

###########
## ANOVA ##
###########

summary(mod.acid.to.int) #p-value: 7.359e-06, adj r^2: 0.853 
Anova(mod.acid.to.int)
#Response: Ma
#           Sum Sq Df F value    Pr(>F)    
#o2        0.01112  1  1.4237   0.25586    
#temp      0.53703  1 68.7610 2.599e-06 ***
#o2:temp   0.04402  1  5.6362   0.03514 *  
#Residuals 0.09372 12    

summary(mod.acid.to) #p-value: 1.109e-05, adj r^2: 0.8006 
Anova(mod.acid.to)
#           Sum Sq Df F value    Pr(>F)    
#o2        0.01112  1  1.0494    0.3243    
#temp      0.53703  1 50.6850 7.824e-06 ***
#Residuals 0.13774 13

summary(mod.acid.t) #p-value: 1.818e-06, adj. R^2: 0.7999 
Anova(mod.acid.t)
#Response: Ma
#          Sum Sq Df F value    Pr(>F)    
#temp      0.64802  1  60.945 1.818e-06 ***
#Residuals 0.14886 14           

#the o2:temp will be selected, as is has > 0.05 more adj. r^2 than the 
#temperature alone and no interactions models and no unaccounted variables.


#####################
## Quality control ##
#####################

plot(mod.acid.to.int) ## Nothing particularly stands out as unusual
olsrr::ols_test_normality(mod.acid.to.int)


#-----------------------------------------------
#   Test                  Statistic       pvalue  
#-----------------------------------------------
#Shapiro-Wilk              0.9385         0.3316 
#Kolmogorov-Smirnov        0.2133         0.4035 
#Cramer-von Mises          4.4792         0.0000 
#Anderson-Darling          0.4429         0.2504 
#-----------------------------------------------

#No abnormalities and the agreement of three normality tests indicates the residuals are normal. 
#I will proceed with the log abundance model and not proceed to the linear abundance model. 

###################
## Plot of Model ##
###################

library(survey)
library(interactions)
library(jtools)

johnson_neyman(mod.acid.to.int, pred = o2, modx = temp, alpha = .05)#For what values of temp the o2 predictor is significant.
#Oxygen is a significant predictor for temperatures outside of [21.37, 49.83].
data.frame(temp,temp <= 21.37, season)
johnson_neyman(mod.acid.to.int, pred = temp, modx = o2, alpha = .05) #[-2814.33, 79.90]
##########
ma_jn_o2_graph <- johnson_neyman(mod.acid.to.int, pred = o2, modx = temp, alpha = .05)
ma_jn_o2_graph$plot$labels$x <- expression("Temperature "~( degree*C))
ma_jn_o2_graph$plot$labels$y <- expression("Slope of O"[2] ~ "(" *mu*"M)")
ma_jn_o2_graph$plot$labels$title <- NULL
ma_jn_o2_graph$plot$theme$legend.position <- "none"
ma_jn_o2_graph_fin <- mp_jn_o2_graph$plot + annotation_custom(grid::textGrob(label = expression(bold("A")), x = unit(0.08, "npc"), 
                                                                             y = unit(0.94, "npc"), gp = gpar(color = "black", 
                                                                                                              fontsize = 17)))

                                                                                                              
ma_jn_temp_graph <- johnson_neyman(mod.acid.to.int, pred = temp, modx = o2, alpha = .05) #[183.36, 224.63]
ma_jn_temp_graph$plot$labels$x <- expression("O"[2] ~ "(" *mu*"M)")
ma_jn_temp_graph$plot$labels$y <- expression("Slope of Temperature "~( degree*C))
ma_jn_temp_graph$plot$labels$title <- NULL
ma_jn_temp_graph_fin <- mp_jn_temp_graph$plot + annotation_custom(grid::textGrob(label = expression(bold("B")), x = unit(0.08, "npc"), 
                                                                                 y = unit(0.94, "npc"), gp = gpar(color = "black", 
                                                                                                                  fontsize = 17)))

grid.arrange(mp_jn_o2_graph_fin,mp_jn_temp_graph_fin,ncol=2,nrow=1,widths=c(1.0,1.3))
