library(car)

geochem <- read.csv(file = "T0_geochem.csv" ,sep = ",", header = TRUE)
row.names(geochem) <- geochem[,1] #set row names
geochem <- geochem[,-1] #remove first column
geochem <- geochem[,-5] #remove empty column at end

o2 <- geochem$Modeled_Initial_Oxygen
ch4 <- geochem$Modeled_Initial_Methane
temp <- c(20.01,20.024,20.099,27.287,27.308,27.253,28.339,28.544,28.575,27.794,
          27.93,20.268,20.153,20.127,18.287,18.276)

set.seed(12)

#testing 
consume.90 <- c(5.3,5.2,36,43.2,7.2,10.8,11.28,13.92,13.2)



#all 
test.model <- lm(rnorm(16)~o2+ch4+temp)
summary(test.model)
vif(test.model)
#       o2       ch4      temp 
# 2.098538  1.612902  2.086154 

#No temp
test.model <- lm(rnorm(16)~o2+ch4)
vif(test.model)
#       o2        ch4     
# 1.098125   1.098125

#Explaining one var as product of others & direct VIF calculations

#Temp
vif(lm(temp ~ o2 + ch4))
#o2      ch4 
#1.098125 1.098125 
summary(lm(temp ~ o2 + ch4))
1/(1-0.5206)
#2.085941 -- fine

#Oxygen
vif(lm(o2 ~ temp + ch4))
summary(lm(o2 ~ temp + ch4))
# temp      ch4 
#1.091645 1.091645 
1/(1-0.5235)
#2.098636 -- fine

#CH4
vif(lm(ch4 ~ temp + o2))
# temp       o2 
#1.420333 1.420333 
summary(lm(ch4 ~ temp + o2))
1/(1-0.38)
#1.612903 -- pretty good 

# 10 = 1/(1-X) -- 10 - 10X = 1 -- -10x = -9 -- x = 0.9

#Minimal colinearity according to VIF

#VIF direct calculation 
#1/(1-R^2)

1/(1-(0.6^2))

cor.test(o2,temp, method="pearson") #mild correlation
#t = 2.4258, df = 14, p-value = 0.02938
cor.test(o2,temp, method="spearman") #no correlation
# p = .08892 and rho = 0.4411765
