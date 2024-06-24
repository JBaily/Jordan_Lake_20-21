library(fBasics)
library(car)
library(olsrr)
library(ggplot2)

###########################
# Regular Proportion Data #
###########################

all_data_prop <- read.csv(file = "Data/16S_all.csv")
all_data_prop$Box <- as.character(all_data_prop$Box)

oct_2020_prop <- subset(all_data_prop, Experiment == "OCT_2020")
june_1_21_prop <- subset(all_data_prop, Experiment == "JL1")
june_2_21_prop <- subset(all_data_prop, Experiment == "JL2")
july_21_prop <- subset(all_data_prop, Experiment == "JL3")
oct_21_prop <- subset(all_data_prop, Experiment == "OCT_21")
nov_21_prop <- subset(all_data_prop, Experiment == "NOV_21")

time_zero_prop <- subset(all_data_prop, Time == 0)

############################
# Rarefied Proportion Data #
############################

all_data_ra <- read.csv(file = "Data/16S_all_ra.csv")
all_data_ra$Box <- as.character(all_data_ra$Box)

oct_2020_ra <- subset(all_data_ra, Experiment == "OCT_2020")
june_1_21_ra <- subset(all_data_ra, Experiment == "JL1")
june_2_21_ra <- subset(all_data_ra, Experiment == "JL2")
july_21_ra <- subset(all_data_ra, Experiment == "JL3")
oct_21_ra <- subset(all_data_ra, Experiment == "OCT_21")
nov_21_ra <- subset(all_data_ra, Experiment == "NOV_21")

time_zero_ra <- subset(all_data_ra, Time == 0)

####################################################################

#T0 correlations all look great. No R^2 less than 0.97.

#############
# Time Zero #
#############

#Testing if it matters which log is used (as long as consistent, of course)
hist(log10(time_zero_prop$Methylophilaceae)) 
hist(log10(time_zero_prop$Methylophilaceae)) 

cor.test(log10(time_zero_prop$Methylophilaceae),log10(time_zero_ra$Methylophilaceae),
         alternative = "two.sided", method = "pearson") # p = 4.701e-12, r = 0.9848588, r^2 = 0.9699469

#While looks different, ultimately does not matter in terms of correlation significance. 

########
# Phil #
########

plot(time_zero_prop$Methylophilaceae,time_zero_ra$Methylophilaceae, main = "Methylophilaceae (Amplicon)",
     xlab = "Regular Proportion", ylab = "Rarefied Proportion") #R^2 = 0.97

hist(log10(time_zero_prop$Methylophilaceae)) #Less skewed

cor.test(log10(time_zero_prop$Methylophilaceae),log10(time_zero_ra$Methylophilaceae), 
         alternative = "two.sided", method = "pearson") # p = 4.701e-12, r = 0.9848588, r^2 = 0.9699469


########
# Acid #
########

plot(time_zero_prop$Methylacidiphilaceae,time_zero_ra$Methylacidiphilaceae, main = "Methylacidiphilaceae (Amplicon)",
     xlab = "Regular Proportion", ylab = "Rarefied Proportion") #R^2 = 0.995

hist(log10(time_zero_prop$Methylacidiphilaceae)) #Less skewed

cor.test(log10(time_zero_prop$Methylacidiphilaceae),log10(time_zero_ra$Methylacidiphilaceae), 
         alternative = "two.sided", method = "pearson") # p < 2.2e-16 , r = 0.9975827, r^2 = 0.9951712


########
# Cocc #
########

plot(time_zero_prop$Methylococcaceae,time_zero_ra$Methylococcaceae, main = "Methylococcaceae (Amplicon)",
     xlab = "Regular Proportion", ylab = "Rarefied Proportion") #R^2 = 0.98

hist(log10(time_zero_prop$Methylococcaceae)) #Less skewed

cor.test(log10(time_zero_prop$Methylococcaceae),log10(time_zero_ra$Methylococcaceae), 
         alternative = "two.sided", method = "pearson") # p = 5.569e-13 , r = 0.9888528, r^2 = 0.9778299

########
# Mona #
########

plot(time_zero_prop$Methylomonadaceae,time_zero_ra$Methylomonadaceae, main = "Methylomonadaceae (Amplicon)",
     xlab = "Regular Proportion", ylab = "Rarefied Proportion") #R^2 = 0.98

hist(log10(time_zero_prop$Methylomonadaceae)) #Less skewed

cor.test(log10(time_zero_prop$Methylomonadaceae),log10(time_zero_ra$Methylomonadaceae), 
         alternative = "two.sided", method = "pearson") # p = 4.103e-13 , r = 0.9893308, r^2 = 0.9787754


########
# Cyst #
########

plot(time_zero_prop$Methylocystaceae,time_zero_ra$Methylocystaceae, main = "Methylocystaceae (Amplicon)",
     xlab = "Regular Proportion", ylab = "Rarefied Proportion") #R^2 = 0.98

hist(log10(time_zero_prop$Methylocystaceae)) 

cor.test(log10(time_zero_prop$Methylocystaceae),log10(time_zero_ra$Methylocystaceae), 
         alternative = "two.sided", method = "pearson") # p = 5.966e-14 , r = 0.9919075, r^2 = 0.9838805


#######################

###############
# Time Series #
###############

################
# October 2020 #
################


########
# Phil #
########

plot(oct_2020_prop$Methylophilaceae,oct_2020_ra$Methylophilaceae, main = "Methylophilaceae (Amplicon)",
     xlab = "Regular Proportion", ylab = "Rarefied Proportion") #R^2 = 0.986

hist(log10(oct_2020_prop$Methylophilaceae)) #Less skewed

cor.test(log10(oct_2020_prop$Methylophilaceae),log10(oct_2020_ra$Methylophilaceae), 
         alternative = "two.sided", method = "pearson") # p = 2.665e-15, r = 0.9929311, r^2 = 0.9859122


########
# Acid #
########

plot(oct_2020_prop$Methylacidiphilaceae,oct_2020_ra$Methylacidiphilaceae, main = "Methylacidiphilaceae (Amplicon)",
     xlab = "Regular Proportion", ylab = "Rarefied Proportion") #R^2 = 0.977

hist(log10(oct_2020_prop$Methylacidiphilaceae)) #Less skewed

cor.test(log10(oct_2020_prop$Methylacidiphilaceae),log10(oct_2020_ra$Methylacidiphilaceae), 
         alternative = "two.sided", method = "pearson") # p = 9.626e-14, r = 0.9885772, r^2 = 0.9772849


########
# Cocc #
########

plot(oct_2020_prop$Methylococcaceae,oct_2020_ra$Methylococcaceae, main = "Methylococcaceae (Amplicon)",
     xlab = "Regular Proportion", ylab = "Rarefied Proportion") #R^2 = 0.97

hist(log10(oct_2020_prop$Methylococcaceae)) #Less skewed

cor.test(log10(oct_2020_prop$Methylococcaceae),log10(oct_2020_ra$Methylococcaceae), 
         alternative = "two.sided", method = "pearson") # p = 5.937e-13, r = 0.9854236, r^2 = 0.9710597

########
# Mona #
########

plot(oct_2020_prop$Methylomonadaceae,oct_2020_ra$Methylomonadaceae, main = "Methylomonadaceae (Amplicon)",
     xlab = "Regular Proportion", ylab = "Rarefied Proportion") #R^2 = 0.98

hist(log10(oct_2020_prop$Methylomonadaceae)) 

cor.test(log10(oct_2020_prop$Methylomonadaceae),log10(oct_2020_ra$Methylomonadaceae), 
         alternative = "two.sided", method = "pearson") # p = 1.338e-15 , r = 0.9935532, r^2 = 0.987148

#proportion alone
ggplot() + geom_point(data = oct_2020_prop, aes(x = Time, y = Methylomonadaceae, shape = Box), 
                      color = "blue")
#rarefied
ggplot() + geom_point(data = oct_2020_ra, aes(x = Time, y = Methylomonadaceae, shape = Box), 
                      color = "orange") +theme_classic() 

#together (issue is at 1 hour for Box 1), and overall pattern is smoother in rarefied
ggplot() + geom_point(data = oct_2020_prop, aes(x = Time, y = Methylomonadaceae, shape = Box), 
                      color = "blue") + 
  geom_line(data = oct_2020_prop, aes(x = Time, y = Methylomonadaceae, linetype = Box), color = "blue") +
  geom_line(data = oct_2020_ra, aes(x = Time, y = Methylomonadaceae, linetype = Box), color = "orange") +
  geom_point(data = oct_2020_ra, aes(x = Time, y = Methylomonadaceae, shape = Box), 
             color = "orange") +theme_classic() 

#Testing w/0 B1T0 & B1T1

oct_20_mod <- oct_2020_prop[-1,]
oct_20_mod <- oct_20_mod[-1,]
oct_20_mod_ra <- oct_2020_ra[-1,]
oct_20_mod_ra <- oct_20_mod_ra[-1,]

plot(oct_20_mod$Methylomonadaceae,oct_20_mod_ra$Methylomonadaceae, main = "Methylomonadaceae (Amplicon)",
     xlab = "Regular Proportion", ylab = "Rarefied Proportion") #R^2 = 0.986

#addition to account for zeroes in dataset
cor.test(log10(oct_20_mod$Methylomonadaceae+0.00005),log10(oct_20_mod_ra$Methylomonadaceae+0.00005), 
         alternative = "two.sided", method = "pearson") # p = 1.855e-13 , r = 0.9930267, r^2 = 0.986102

########
# Cyst #
########

plot(oct_2020_prop$Methylocystaceae,oct_2020_ra$Methylocystaceae, main = "Methylocystaceae (Amplicon)",
     xlab = "Regular Proportion", ylab = "Rarefied Proportion") #R^2 = 0.985

hist(log10(oct_2020_prop$Methylocystaceae)) 

cor.test(log10(oct_2020_prop$Methylocystaceae),log10(oct_2020_ra$Methylocystaceae), 
         alternative = "two.sided", method = "pearson") # p = 2.679e-15 , r = 0.9929263, r^2 = 0.9859026

#############

#######
# JL1 #
#######

#Rest ---  R^2 above 0.954
#Methylococcaceae is 0.83 -- same broad patterns, but specifics vary. 

########
# Phil #
########

plot(june_1_21_prop$Methylophilaceae,june_1_21_ra$Methylophilaceae, main = "Methylophilaceae (Amplicon)",
     xlab = "Regular Proportion", ylab = "Rarefied Proportion") #R^2 = 0.968

hist(log10(june_1_21_prop$Methylophilaceae)) #Less skewed

cor.test(log10(june_1_21_prop$Methylophilaceae),log10(june_1_21_ra$Methylophilaceae), 
         alternative = "two.sided", method = "pearson") # p < 2.2e-16, r = 0.9843501, r^2 = 0.9689451


########
# Acid #
########

plot(june_1_21_prop$Methylacidiphilaceae,june_1_21_ra$Methylacidiphilaceae, main = "Methylacidiphilaceae (Amplicon)",
     xlab = "Regular Proportion", ylab = "Rarefied Proportion") #R^2 = 0.995

hist(log10(june_1_21_prop$Methylacidiphilaceae)) #Less skewed

cor.test(log10(june_1_21_prop$Methylacidiphilaceae),log10(june_1_21_ra$Methylacidiphilaceae), 
         alternative = "two.sided", method = "pearson") # p < 2.2e-16 , r = 0.9976678, r^2 = 0.995341


########
# Cocc #
########

plot(june_1_21_prop$Methylococcaceae,june_1_21_ra$Methylococcaceae, main = "Methylococcaceae (Amplicon)",
     xlab = "Regular Proportion", ylab = "Rarefied Proportion") #R^2 = 0.835

hist(log10(june_1_21_prop$Methylococcaceae)) #Less skewed

cor.test(log10(june_1_21_prop$Methylococcaceae),log10(june_1_21_ra$Methylococcaceae), 
         alternative = "two.sided", method = "pearson") # p = 2.74e-11, r = 0.9140375, r^2 = 0.8354646

#proportion alone
ggplot() + geom_point(data = june_1_21_prop, aes(x = Time, y = Methylococcaceae, shape = Box), 
                      color = "blue")
#rarefied
ggplot() + geom_point(data = june_1_21_ra, aes(x = Time, y = Methylococcaceae, shape = Box), 
                      color = "orange") +theme_classic() 

#Smoother with regular
ggplot() + geom_point(data = june_1_21_prop, aes(x = Time, y = Methylococcaceae, shape = Box), 
                      color = "blue") + 
  geom_line(data = june_1_21_prop, aes(x = Time, y = Methylococcaceae, linetype = Box), color = "blue") +
  geom_line(data = june_1_21_ra, aes(x = Time, y = Methylococcaceae, linetype = Box), color = "orange") +
  geom_point(data = june_1_21_ra, aes(x = Time, y = Methylococcaceae, shape = Box), 
             color = "orange") +theme_classic() 


########
# Mona #
########

plot(june_1_21_prop$Methylomonadaceae,june_1_21_ra$Methylomonadaceae, main = "Methylomonadaceae (Amplicon)",
     xlab = "Regular Proportion", ylab = "Rarefied Proportion") #R^2 = 0.95

#addition to account for zeroes in dataset
hist(log10(june_1_21_prop$Methylomonadaceae+0.00005)) #Less skewed

cor.test(log10(june_1_21_prop$Methylomonadaceae+0.00005),log10(june_1_21_ra$Methylomonadaceae+0.00005), 
         alternative = "two.sided", method = "pearson") # p < 2.2e-16, r = 0.9754358, r^2 = 0.951475

########
# Cyst #
########

plot(june_1_21_prop$Methylocystaceae,june_1_21_ra$Methylocystaceae, main = "Methylocystaceae (Amplicon)",
     xlab = "Regular Proportion", ylab = "Rarefied Proportion") #R^2 = 0.94

hist(log10(june_1_21_prop$Methylocystaceae)) 

cor.test(log10(june_1_21_prop$Methylocystaceae),log10(june_1_21_ra$Methylocystaceae), 
         alternative = "two.sided", method = "pearson") # p < 2.2e-16 , r = 0.9691773, r^2 = 0.9393046

#############

#######
# JL2 #
#######

#Methylococcaceae R^2 is 0.50. Not great, but also was at low conc. to begin with. 
# Regular seems more consistent than the rarefied.

#Methylomonadaceae R^2 is 0.74. Better than Methylococcaceae. Overall pattern broadly matches between 
#the regular and rarefied, just the specifics vary. Also generally low concentrations. 

#Methylocystaceae at 0.85

#Rest above 0.95

########
# Phil #
########

plot(june_2_21_prop$Methylophilaceae,june_2_21_ra$Methylophilaceae, main = "Methylophilaceae (Amplicon)",
     xlab = "Regular Proportion", ylab = "Rarefied Proportion") #R^2 = 0.91

hist(log10(june_2_21_prop$Methylophilaceae)) #Less skewed

cor.test(log10(june_2_21_prop$Methylophilaceae),log10(june_2_21_ra$Methylophilaceae), 
         alternative = "two.sided", method = "pearson") # p = 2.614e-12, r = 0.9524374, r^2 = 0.907137


########
# Acid #
########

plot(june_2_21_prop$Methylacidiphilaceae,june_2_21_ra$Methylacidiphilaceae, main = "Methylacidiphilaceae (Amplicon)",
     xlab = "Regular Proportion", ylab = "Rarefied Proportion") #R^2 = 0.995

hist(log10(june_2_21_prop$Methylacidiphilaceae)) #Less skewed

cor.test(log10(june_2_21_prop$Methylacidiphilaceae),log10(june_2_21_ra$Methylacidiphilaceae), 
         alternative = "two.sided", method = "pearson") # p < 2.2e-16 , r = 0.9977265, r^2 = 0.9954582


########
# Cocc #
########

plot(june_2_21_prop$Methylococcaceae,june_2_21_ra$Methylococcaceae, main = "Methylococcaceae (Amplicon)",
     xlab = "Regular Proportion", ylab = "Rarefied Proportion") #R^2 = 0.50

hist(log10(june_2_21_prop$Methylococcaceae)) #Less skewed

cor.test(log10(june_2_21_prop$Methylococcaceae),log10(june_2_21_ra$Methylococcaceae), 
         alternative = "two.sided", method = "pearson") # p = 0.0001739 , r = 0.7046921, r^2 = 0.496591

#proportion alone
ggplot() + geom_point(data = june_2_21_prop, aes(x = Time, y = Methylococcaceae, shape = Box), 
                      color = "blue")
#rarefied
ggplot() + geom_point(data = june_2_21_ra, aes(x = Time, y = Methylococcaceae, shape = Box), 
                      color = "orange") +theme_classic() 

#Smoother with regular
ggplot() + geom_point(data = june_2_21_prop, aes(x = Time, y = Methylococcaceae, shape = Box), 
                      color = "blue") + 
  geom_line(data = june_2_21_prop, aes(x = Time, y = Methylococcaceae, linetype = Box), color = "blue") +
  geom_line(data = june_2_21_ra, aes(x = Time, y = Methylococcaceae, linetype = Box), color = "orange") +
  geom_point(data = june_2_21_ra, aes(x = Time, y = Methylococcaceae, shape = Box), 
             color = "orange") +theme_classic() 


########
# Mona #
########

plot(june_2_21_prop$Methylomonadaceae,june_2_21_ra$Methylomonadaceae, main = "Methylomonadaceae (Amplicon)",
     xlab = "Regular Proportion", ylab = "Rarefied Proportion") #R^2 = 0.74

hist(log10(june_2_21_prop$Methylomonadaceae+0.00005)) #Less skewed

cor.test(log10(june_2_21_prop$Methylomonadaceae+0.00005),log10(june_2_21_ra$Methylomonadaceae+0.00005), 
         alternative = "two.sided", method = "pearson") # p = 1.479e-07, r = 0.8596352, r^2 = 0.7389727

#proportion alone
ggplot() + geom_point(data = june_2_21_prop, aes(x = Time, y = Methylomonadaceae, shape = Box), 
                      color = "blue")
#rarefied
ggplot() + geom_point(data = june_2_21_ra, aes(x = Time, y = Methylomonadaceae, shape = Box), 
                      color = "orange") +theme_classic() 

#overall pattern agrees
ggplot() + geom_point(data = june_2_21_prop, aes(x = Time, y = Methylomonadaceae, shape = Box), 
                      color = "blue") + 
  geom_line(data = june_2_21_prop, aes(x = Time, y = Methylomonadaceae, linetype = Box), color = "blue") +
  geom_line(data = june_2_21_ra, aes(x = Time, y = Methylomonadaceae, linetype = Box), color = "orange") +
  geom_point(data = june_2_21_ra, aes(x = Time, y = Methylomonadaceae, shape = Box), 
             color = "orange") +theme_classic() 


########
# Cyst #
########

plot(june_2_21_prop$Methylocystaceae,june_2_21_ra$Methylocystaceae, main = "Methylocystaceae (Amplicon)",
     xlab = "Regular Proportion", ylab = "Rarefied Proportion") #R^2 = 0.85

#addition to account for zeroes in dataset
hist(log10(june_2_21_prop$Methylocystaceae+0.005)) 

cor.test(log10(june_2_21_prop$Methylocystaceae+0.005),log10(june_2_21_ra$Methylocystaceae+0.005), 
         alternative = "two.sided", method = "pearson") # p = 4.211e-10, r = 0.9218179, r^2 = 0.8497482

ggplot() + geom_point(data = june_2_21_prop, aes(x = Time, y = Methylocystaceae, shape = Box), 
                      color = "blue") + 
  geom_line(data = june_2_21_prop, aes(x = Time, y = Methylocystaceae, linetype = Box), color = "blue") +
  geom_line(data = june_2_21_ra, aes(x = Time, y = Methylocystaceae, linetype = Box), color = "orange") +
  geom_point(data = june_2_21_ra, aes(x = Time, y = Methylocystaceae, shape = Box), 
             color = "orange") +theme_classic() 

#################

#######
# JL3 #
#######

#Methylomonadaceae is 0.81
#Methylocystaceae is 0.86

########
# Phil #
########

plot(july_21_prop$Methylophilaceae,july_21_ra$Methylophilaceae, main = "Methylophilaceae (Amplicon)",
     xlab = "Regular Proportion", ylab = "Rarefied Proportion") #R^2 = 0.93

hist(log10(july_21_prop$Methylophilaceae)) #Less skewed

cor.test(log10(july_21_prop$Methylophilaceae),log10(july_21_ra$Methylophilaceae), 
         alternative = "two.sided", method = "pearson") # p = 6.416e-09, r = 0.9648241, r^2 = 0.9308855


########
# Acid #
########

plot(july_21_prop$Methylacidiphilaceae,july_21_ra$Methylacidiphilaceae, main = "Methylacidiphilaceae (Amplicon)",
     xlab = "Regular Proportion", ylab = "Rarefied Proportion") #R^2 = 0.995

hist(log10(july_21_prop$Methylacidiphilaceae)) #Less skewed

cor.test(log10(july_21_prop$Methylacidiphilaceae),log10(july_21_ra$Methylacidiphilaceae), 
         alternative = "two.sided", method = "pearson") # p = 2.444e-16, r = 0.9974905, r^2 = 0.9949873


########
# Cocc #
########

plot(july_21_prop$Methylococcaceae,july_21_ra$Methylococcaceae, main = "Methylococcaceae (Amplicon)",
     xlab = "Regular Proportion", ylab = "Rarefied Proportion") #R^2 = 0.986

hist(log10(july_21_prop$Methylococcaceae)) #Less skewed

cor.test(log10(july_21_prop$Methylococcaceae),log10(july_21_ra$Methylococcaceae), 
         alternative = "two.sided", method = "pearson") # p = 1.57e-13, r = 0.9932039, r^2 = 0.986454

########
# Mona #
########

plot(july_21_prop$Methylomonadaceae,july_21_ra$Methylomonadaceae, main = "Methylomonadaceae (Amplicon)",
     xlab = "Regular Proportion", ylab = "Rarefied Proportion") #R^2 = 0.81

#addition to account for zeroes in dataset
hist(log10(july_21_prop$Methylomonadaceae+0.00005)) #Less skewed

cor.test(log10(july_21_prop$Methylomonadaceae+0.00005),log10(july_21_ra$Methylomonadaceae+0.00005), 
         alternative = "two.sided", method = "pearson") # p = 4.118e-06, r = 0.9026523, r^2 = 0.8147812

#together (issue is at 1 hour for Box 1), and overall pattern is smoother in rarefied
ggplot() + geom_point(data = july_21_prop, aes(x = Time, y = Methylomonadaceae, shape = Box), 
                      color = "blue") + 
  geom_line(data = july_21_prop, aes(x = Time, y = Methylomonadaceae, linetype = Box), color = "blue") +
  geom_line(data = july_21_ra, aes(x = Time, y = Methylomonadaceae, linetype = Box), color = "orange") +
  geom_point(data = july_21_ra, aes(x = Time, y = Methylomonadaceae, shape = Box), 
             color = "orange") +theme_classic() 


########
# Cyst #
########

plot(july_21_prop$Methylocystaceae,july_21_ra$Methylocystaceae, main = "Methylocystaceae (Amplicon)",
     xlab = "Regular Proportion", ylab = "Rarefied Proportion") #R^2 = 0.86

hist(log10(july_21_prop$Methylocystaceae)) 

cor.test(log10(july_21_prop$Methylocystaceae),log10(july_21_ra$Methylocystaceae), 
         alternative = "two.sided", method = "pearson") # p = 7.661e-07, r = 0.9254939, r^2 = 0.856539

ggplot() + geom_point(data = july_21_prop, aes(x = Time, y = Methylocystaceae, shape = Box), 
                      color = "blue") + 
  geom_line(data = july_21_prop, aes(x = Time, y = Methylocystaceae, linetype = Box), color = "blue") +
  geom_line(data = july_21_ra, aes(x = Time, y = Methylocystaceae, linetype = Box), color = "orange") +
  geom_point(data = july_21_ra, aes(x = Time, y = Methylocystaceae, shape = Box), 
             color = "orange") +theme_classic() 


####################

#########
# Allyl #
#########

#All are great

########
# Phil #
########

plot(oct_21_prop$Methylophilaceae,oct_21_ra$Methylophilaceae, main = "Methylophilaceae (Amplicon)",
     xlab = "Regular Proportion", ylab = "Rarefied Proportion") #R^2 = 0.986

hist(log10(oct_21_prop$Methylophilaceae)) #Less skewed

cor.test(log10(oct_21_prop$Methylophilaceae),log10(oct_21_ra$Methylophilaceae), 
         alternative = "two.sided", method = "pearson") # p < 2.2e-16, r = 0.9931612, r^2 = 0.9863692


########
# Acid #
########

plot(oct_21_prop$Methylacidiphilaceae,oct_21_ra$Methylacidiphilaceae, main = "Methylacidiphilaceae (Amplicon)",
     xlab = "Regular Proportion", ylab = "Rarefied Proportion") #R^2 = 0.99

hist(log10(oct_21_prop$Methylacidiphilaceae)) #Less skewed

cor.test(log10(oct_21_prop$Methylacidiphilaceae),log10(oct_21_ra$Methylacidiphilaceae), 
         alternative = "two.sided", method = "pearson") # p < 2.2e-16 , r = 0.9968439, r^2 = 0.9936978


########
# Cocc #
########

plot(oct_21_prop$Methylococcaceae,oct_21_ra$Methylococcaceae, main = "Methylococcaceae (Amplicon)",
     xlab = "Regular Proportion", ylab = "Rarefied Proportion") #R^2 = 0.94

#addition to account for zeroes in dataset
hist(log10(oct_21_prop$Methylococcaceae+0.00005)) #Less skewed

cor.test(log10(oct_21_prop$Methylococcaceae+0.00005),log10(oct_21_ra$Methylococcaceae+0.00005), 
         alternative = "two.sided", method = "pearson") # p = 4.605e-16, r = 0.968903, r^2 = 0.938773

########
# Mona #
########

plot(oct_21_prop$Methylomonadaceae,oct_21_ra$Methylomonadaceae, main = "Methylomonadaceae (Amplicon)",
     xlab = "Regular Proportion", ylab = "Rarefied Proportion") #R^2 = 0.91

#addition to account for zeroes in dataset
hist(log10(oct_21_prop$Methylomonadaceae+0.00005)) #Less skewed

cor.test(log10(oct_21_prop$Methylomonadaceae+0.00005),log10(oct_21_ra$Methylomonadaceae+0.00005), 
         alternative = "two.sided", method = "pearson") # p = 3.169e-14, r = 0.9555001, r^2 = 0.9129804


########
# Cyst #
########

plot(oct_21_prop$Methylocystaceae,oct_21_ra$Methylocystaceae, main = "Methylocystaceae (Amplicon)",
     xlab = "Regular Proportion", ylab = "Rarefied Proportion") #R^2 = 0.97

hist(log10(oct_21_prop$Methylocystaceae)) 

cor.test(log10(oct_21_prop$Methylocystaceae),log10(oct_21_ra$Methylocystaceae), 
         alternative = "two.sided", method = "pearson") # p < 2.2e-16 , r = 0.9873709, r^2 = 0.9749013

################

##########
# AlMeOH #
##########

#All great aside from Mona (r^2 = 0.72). Broad pattern (decreasing) matches for both, specifics differ. 

########
# Phil #
########

plot(nov_21_prop$Methylophilaceae,nov_21_ra$Methylophilaceae, main = "Methylophilaceae (Amplicon)",
     xlab = "Regular Proportion", ylab = "Rarefied Proportion") #R^2 = 0.98

hist(log10(nov_21_prop$Methylophilaceae)) #Less skewed

cor.test(log10(nov_21_prop$Methylophilaceae),log10(nov_21_ra$Methylophilaceae), 
         alternative = "two.sided", method = "pearson") # p = 5.136e-13, r = 0.9889813, r^2 = 0.978084


########
# Acid #
########

plot(nov_21_prop$Methylacidiphilaceae,nov_21_ra$Methylacidiphilaceae, main = "Methylacidiphilaceae (Amplicon)",
     xlab = "Regular Proportion", ylab = "Rarefied Proportion") #R^2 = 0.989

hist(log10(nov_21_prop$Methylacidiphilaceae)) #Less skewed

cor.test(log10(nov_21_prop$Methylacidiphilaceae),log10(nov_21_ra$Methylacidiphilaceae), 
         alternative = "two.sided", method = "pearson") # p = 4.105e-15, r = 0.9944842, r^2 = 0.9889988


########
# Cocc #
########

plot(nov_21_prop$Methylococcaceae,nov_21_ra$Methylococcaceae, main = "Methylococcaceae (Amplicon)",
     xlab = "Regular Proportion", ylab = "Rarefied Proportion") #R^2 = 0.90

#addition to account for zeroes in dataset
hist(log10(nov_21_prop$Methylococcaceae+0.0005)) #Less skewed

cor.test(log10(nov_21_prop$Methylococcaceae+0.0005),log10(nov_21_ra$Methylococcaceae+0.0005), 
         alternative = "two.sided", method = "pearson") # p = 2.079e-08, r = 0.949082, r^2 = 0.9007566

########
# Mona #
########

plot(nov_21_prop$Methylomonadaceae,nov_21_ra$Methylomonadaceae, main = "Methylomonadaceae (Amplicon)",
     xlab = "Regular Proportion", ylab = "Rarefied Proportion") #R^2 = 0.72

#addition to account for zeroes in dataset
hist(log10(nov_21_prop$Methylomonadaceae+0.00005)) #Less skewed

cor.test(log10(nov_21_prop$Methylomonadaceae+0.00005),log10(nov_21_ra$Methylomonadaceae+0.00005), 
         alternative = "two.sided", method = "pearson") # p = 3.431e-05, r = 0.8473083, r^2 = 0.7179314

#proportion alone
ggplot() + geom_point(data = nov_21_prop, aes(x = Time, y = Methylomonadaceae, shape = Box), 
                      color = "blue")
#rarefied
ggplot() + geom_point(data = nov_21_ra, aes(x = Time, y = Methylomonadaceae, shape = Box), 
                      color = "orange") +theme_classic() 

#together (issue is at 1 hour for Box 1), and overall pattern is smoother in rarefied
ggplot() + geom_point(data = nov_21_prop, aes(x = Time, y = Methylomonadaceae, shape = Box), 
                      color = "blue") + 
  geom_line(data = nov_21_prop, aes(x = Time, y = Methylomonadaceae, linetype = Box), color = "blue") +
  geom_line(data = nov_21_ra, aes(x = Time, y = Methylomonadaceae, linetype = Box), color = "orange") +
  geom_point(data = nov_21_ra, aes(x = Time, y = Methylomonadaceae, shape = Box), 
             color = "orange") +theme_classic() 


########
# Cyst #
########

plot(nov_21_prop$Methylocystaceae,nov_21_ra$Methylocystaceae, main = "Methylocystaceae (Amplicon)",
     xlab = "Regular Proportion", ylab = "Rarefied Proportion") #R^2 = 0.988

hist(log10(nov_21_prop$Methylocystaceae)) 

cor.test(log10(nov_21_prop$Methylocystaceae),log10(nov_21_ra$Methylocystaceae), 
         alternative = "two.sided", method = "pearson") # p = 6.987e-15, r = 0.9940478, r^2 = 0.988131

