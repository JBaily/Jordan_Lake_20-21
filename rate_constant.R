library(ggplot2)
library(ggtext)
library(scales)
library(grid)
library(gridExtra)

## Rate constant correlations

#Microbiological Data
All_16S <- read.csv(file = "16S_all.csv",header = TRUE)
T0_16S <- All_16S[All_16S$Time==0,]

Mco <- log10(T0_16S$Methylococcaceae)
Mm <- log10(T0_16S$Methylomonadaceae)
Mcy <- log10(T0_16S$Methylocystaceae)
Ma <- log10(T0_16S$Methylacidiphilaceae)
Mp <- log10(T0_16S$Methylophilaceae)

############
# Averages #
############

oct_2020 <- subset(All_16S, Experiment == "OCT_2020")
june_1_21 <- subset(All_16S, Experiment == "JL1")
june_2_21 <- subset(All_16S, Experiment == "JL2")
july_2021 <- subset(All_16S, Experiment == "JL3")
oct_2021 <- subset(All_16S, Experiment == "OCT_21")
nov_2021 <- subset(All_16S, Experiment == "NOV_21")

Mco_oct_2020_avg <- unname(sapply(split(oct_2020$Methylococcaceae, oct_2020$Box), mean))
Mm_oct_2020_avg <- unname(sapply(split(oct_2020$Methylomonadaceae, oct_2020$Box), mean))
Mcy_oct_2020_avg <- unname(sapply(split(oct_2020$Methylocystaceae, oct_2020$Box), mean))
Ma_oct_2020_avg <- unname(sapply(split(oct_2020$Methylacidiphilaceae, oct_2020$Box), mean))
Mp_oct_2020_avg <- unname(sapply(split(oct_2020$Methylophilaceae, oct_2020$Box), mean))

Mco_june_1_21_avg <- unname(sapply(split(june_1_21$Methylococcaceae, june_1_21$Box), mean))
Mm_june_1_21_avg <- unname(sapply(split(june_1_21$Methylomonadaceae, june_1_21$Box), mean))
Mcy_june_1_21_avg <- unname(sapply(split(june_1_21$Methylocystaceae, june_1_21$Box), mean))
Ma_june_1_21_avg <- unname(sapply(split(june_1_21$Methylacidiphilaceae, june_1_21$Box), mean))
Mp_june_1_21_avg <- unname(sapply(split(june_1_21$Methylophilaceae, june_1_21$Box), mean))

Mco_june_2_21_avg <- unname(sapply(split(june_2_21$Methylococcaceae, june_2_21$Box), mean))
Mm_june_2_21_avg <- unname(sapply(split(june_2_21$Methylomonadaceae, june_2_21$Box), mean))
Mcy_june_2_21_avg <- unname(sapply(split(june_2_21$Methylocystaceae, june_2_21$Box), mean))
Ma_june_2_21_avg <- unname(sapply(split(june_2_21$Methylacidiphilaceae, june_2_21$Box), mean))
Mp_june_2_21_avg <- unname(sapply(split(june_2_21$Methylophilaceae, june_2_21$Box), mean))

Mco_july_2021_avg <- unname(sapply(split(july_2021$Methylococcaceae, july_2021$Box), mean))
Mm_july_2021_avg <- unname(sapply(split(july_2021$Methylomonadaceae, july_2021$Box), mean))
Mcy_july_2021_avg <- unname(sapply(split(july_2021$Methylocystaceae, july_2021$Box), mean))
Ma_july_2021_avg <- unname(sapply(split(july_2021$Methylacidiphilaceae, july_2021$Box), mean))
Mp_july_2021_avg <- unname(sapply(split(july_2021$Methylophilaceae, july_2021$Box), mean))

Mco_oct_2021_avg <- unname(sapply(split(oct_2021$Methylococcaceae, oct_2021$Box), mean))
Mm_oct_2021_avg <- unname(sapply(split(oct_2021$Methylomonadaceae, oct_2021$Box), mean))
Mcy_oct_2021_avg <- unname(sapply(split(oct_2021$Methylocystaceae, oct_2021$Box), mean))
Ma_oct_2021_avg <- unname(sapply(split(oct_2021$Methylacidiphilaceae, oct_2021$Box), mean))
Mp_oct_2021_avg <- unname(sapply(split(oct_2021$Methylophilaceae, oct_2021$Box), mean))

Mco_nov_2021_avg <- unname(sapply(split(nov_2021$Methylococcaceae, nov_2021$Box), mean))
Mm_nov_2021_avg <- unname(sapply(split(nov_2021$Methylomonadaceae, nov_2021$Box), mean))
Mcy_nov_2021_avg <- unname(sapply(split(nov_2021$Methylocystaceae, nov_2021$Box), mean))
Ma_nov_2021_avg <- unname(sapply(split(nov_2021$Methylacidiphilaceae, nov_2021$Box), mean))
Mp_nov_2021_avg <- unname(sapply(split(nov_2021$Methylophilaceae, nov_2021$Box), mean))


Mco_avg <- log10(c(Mco_oct_2020_avg,Mco_june_1_21_avg,Mco_june_2_21_avg,
                      Mco_july_2021_avg,Mco_oct_2021_avg,Mco_nov_2021_avg))
Mm_avg <- log10(c(Mm_oct_2020_avg,Mm_june_1_21_avg,Mm_june_2_21_avg,
                  Mm_july_2021_avg,Mm_oct_2021_avg,Mm_nov_2021_avg))
Mcy_avg <- log10(c(Mcy_oct_2020_avg,Mcy_june_1_21_avg,Mcy_june_2_21_avg,
                   Mcy_july_2021_avg,Mcy_oct_2021_avg,Mcy_nov_2021_avg))
Ma_avg <- log10(c(Ma_oct_2020_avg,Ma_june_1_21_avg,Ma_june_2_21_avg,
                  Ma_july_2021_avg,Ma_oct_2021_avg,Ma_nov_2021_avg))
Mp_avg <- log10(c(Mp_oct_2020_avg,Mp_june_1_21_avg,Mp_june_2_21_avg,
                  Mp_july_2021_avg,Mp_oct_2021_avg,Mp_nov_2021_avg))

###############

geochem <- read.csv(file = "T0_geochem.csv" ,sep = ",", header = TRUE)
row.names(geochem) <- geochem[,1] #set row names
geochem <- geochem[,-1] #remove first column
geochem <- geochem[,-5] #remove empty column at end

rate_constant <- geochem$Rate_Constant

#########################################################
## abundance relationship between methanos and methylos##
#########################################################

methano <- log(exp(Mco) + exp(Mm) + exp(Mcy))
methylo <- log(exp(Mp) + exp(Ma))

t.test(methano, methylo, alternative = "two.sided", paired = FALSE)
#t = -5.3545, df = 19.829, p-value = 3.147e-05
# 95% confidence interval: 
# -0.8599977 -0.3775966
# mean of x      mean of y
# -1.768112      -1.149315 

t.test(methano, methylo, alternative = "two.sided", paired = TRUE)
#t = -4.6209, df = 15, p-value = 0.0003329
#95 percent confidence interval:
# -0.9042276 -0.3333667
#mean of the differences
#-0.6187972 

###########
# Pearson #
###########

#Regular 
cor.test(rate_constant,Ma, method = "pearson") 
# r = -0.6131463, t = -2.9041, p-value = 0.01155 # -- *
# 0.38
cor.test(rate_constant,Mp, method = "pearson") 
# r = 0.3691961, t = 1.4864, p-value = 0.1593 # -- NS
# 0.14
cor.test(rate_constant,Mcy, method = "pearson")
# r = 0.6689373, t = 3.3672, p-value = 0.004603 -- **
# 0.45
cor.test(rate_constant,Mm, method = "pearson")
# r = 0.5513477, t = 2.4727, p-value = 0.02684 -- *
# 0.31

cor.test(rate_constant,Mco, method = "pearson") # -- *****
# r = 0.8561366, t = 6.1991, p-value = 2.318e-05
# 0.73


#Average
cor.test(rate_constant,Ma_avg, method = "pearson") 
# r = -0.4600364, t = -1.9386, p-value = 0.07298 # -- NS
cor.test(rate_constant,Mp_avg, method = "pearson") 
# r = 0.2554762, t = 0.98871, p-value = 0.3396 # -- NS
cor.test(rate_constant,Mcy_avg, method = "pearson")
# r = 0.7214054, t = 3.8978, p-value = 0.001609 -- **
cor.test(rate_constant,Mm_avg, method = "pearson")
# r = 0.5672939, t = 2.5775, p-value = 0.02192 -- *
cor.test(rate_constant,Mco_avg, method = "pearson") 
# r = 0.68261, t = 3.495, p-value = 0.00357 -- ** 


###########
# Results #
###########

#Methylococaceae T0 has the largest significant correlation. 

###########
## Plots ##
###########

rate_Mco <- ggplot() + 
  geom_point(aes(x = 10^(Mco), y = rate_constant), color = "blue", size =3.5) +
  labs(x=bquote("Relative Abundance of"~italic("Methylococcaceae")), 
       y = expression(paste("Rate Constant (k, h", r^-1, ")"))) + 
  annotation_custom(grid::textGrob(label = "P < 0.001",
                                   x = unit(0.135, "npc"), y = unit(0.84, "npc"), 
                                   gp = gpar(color = "black", fontsize=11))) +
  annotation_custom(grid::textGrob(label = expression(R^2*" = 0.73"),
                                   x = unit(0.134, "npc"), y = unit(0.90, "npc"), 
                                   gp = gpar(color = "black", fontsize=11))) +
  scale_x_continuous(trans = "log10", breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  geom_smooth(aes(x = 10^(Mco), y = rate_constant), method = "lm", se=FALSE,
              color = "black", linewidth=0.5) +
  guides(linetype = "none", color = "none") + theme_classic() + theme(axis.text = element_text(size=rel(1.3)))

rate_Mco
