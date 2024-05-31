All_16S <- read.csv(file = "16S_all.csv",header = TRUE)
T0_16S <- All_16S[All_16S$Time==0,]

Mco <- T0_16S$Methylococcaceae
Mm <- T0_16S$Methylomonadaceae
Mcy <- T0_16S$Methylocystaceae
Ma <- T0_16S$Methylacidiphilaceae
Mp <- T0_16S$Methylophilaceae

geochem <- read.csv(file = "T0_geochem.csv" ,sep = ",", header = TRUE)
row.names(geochem) <- geochem[,1] #set row names
geochem <- geochem[,-1] #remove first column
geochem <- geochem[,-5] #remove empty column at end

o2 <- geochem$Modeled_Initial_Oxygen
ch4 <- geochem$Modeled_Initial_Methane
temp <- c(20.01,20.024,20.099,27.287,27.308,27.253,28.339,
          28.544,28.575,27.794,27.93,20.268,20.153,20.127,18.287,18.276)
rate_constant<-geochem$Rate_Constant

#Test normality 
library(fBasics)

###########################
#For Microbiological Data #
###########################

boxplot(Mco)
qqnormPlot(Mco) #bad
qqnormPlot(log(Mco)) #much better
normalTest(Mco, method = "sw") #Not normal 
normalTest(log(Mco), method = "sw")#Normal


boxplot(Mm)
qqnormPlot(Mm) #bad
qqnormPlot(log(Mm)) #much better
normalTest(Mm, method = "sw") #Not normal 
normalTest(log(Mm), method = "sw")#Normal


boxplot(Mcy)
qqnormPlot(Mcy) #bad
qqnormPlot(log(Mcy)) #much better
normalTest(Mcy, method = "sw") #Not normal 
normalTest(log(Mcy), method = "sw")#Normal


boxplot(Ma)
qqnormPlot(Ma) #Fine
qqnormPlot(log(Ma)) #About the same
normalTest(Ma, method = "sw") #Normal 
normalTest(log(Ma), method = "sw") #Still Normal


boxplot(Mp)
qqnormPlot(Mp) #bad
qqnormPlot(log(Mp)) #much better
normalTest(Mp, method = "sw") #Not normal 
normalTest(log(Mp), method = "sw")#Normal

#######################
#For Geochemical Data #
#######################

# Rate constant - probably fine
boxplot(rate_constant)
boxplot(log(rate_constant))
qqnormPlot(rate_constant) #fine
qqnormPlot(log(rate_constant)) #not much different
normalTest(rate_constant, method = "sw") #Very slightly not normal P = 0.06983  
normalTest(log(rate_constant), method = "sw")#Just under significance

# Oxygen - fine (regular) 
boxplot(o2)
qqnormPlot(o2) #fine
qqnormPlot(log(o2)) #a little worse
normalTest(o2, method = "sw")  #Just under significance 
normalTest(log(o2), method = "sw") #Normal 

#Methane -- fine (regular)
boxplot(ch4)
qqnormPlot(ch4) #fine
qqnormPlot(log(ch4)) #not much different
normalTest(ch4, method = "sw") #Normal  
normalTest(log(ch4), method = "sw") #Just under significance

#Temperature -- Temperature is usually normally distributed here, we just sampled 
#during two specific high and low periods. 

hist(temp) #Looks like was taken from a normal distribution, but with middle taken out.
boxplot(temp)
qqnormPlot(temp) #bimodal
qqnormPlot(log(temp)) #same 
normalTest(temp, method = "sw") #Not normal -- bimodal
normalTest(log(temp), method = "sw")#Same 

