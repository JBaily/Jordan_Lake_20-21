#########################################
# Comparisons of amplicon and qPCR data #
#########################################

library(ggplot2)
library(grid)
library(ggtext)
library(gridExtra)
library(scales)

#############
# qPCR data #
#############

T0_qPCR <- read.csv(file = "qPCR_JL.csv", header = TRUE)
Oct_2020_qPCR <- read.csv(file = "Mco_time_series_relative.csv", header = TRUE)

##################
## Amplicon Data ##
###################

#Microbiological Data
All_16S <- read.csv(file = "16S_all.csv",header = TRUE)
T0_16S <- All_16S[All_16S$Time==0,]
T0_16S_select <- T0_16S[-16,]
T0_16S_select <- T0_16S_select[-15,]
T0_16S_select <- T0_16S_select[-13,]
T0_16S_select <- T0_16S_select[-1,]

Mco <- log10(T0_16S_select$Methylococcaceae)
Mcy <- log10(T0_16S_select$Methylocystaceae)

oct_2020 <- subset(All_16S, Experiment == "OCT_2020")
oct_2020 <- oct_2020[-1,]#remove B1T0
oct_2020 <- oct_2020[-1,]#remove B1T1
Mco_oct_2020 <- oct_2020$Methylococcaceae

####################
# Geochemical Data #
####################

geochem <- read.csv(file = "T0_geochem.csv" ,sep = ",", header = TRUE)
row.names(geochem) <- geochem[,1] #set row names
geochem <- geochem[,-1] #remove first column
geochem <- geochem[,-5] #remove empty column at end
geochem_qPCR <- geochem[-16,]
geochem_qPCR <- geochem_qPCR[-15,]
geochem_qPCR <- geochem_qPCR[-13,]
geochem_qPCR <- geochem_qPCR[-1,]


o2 <- geochem_qPCR$Modeled_Initial_Oxygen
ch4 <- geochem_qPCR$Modeled_Initial_Methane
temp <- c(20.024,20.099,27.287,27.308,27.253,28.339,28.544,28.575,27.794,
          27.93,20.268,20.127)
season <- c("O","O","J","J","J","J","J","J","J","J","O","O")

############
## T0 Mco ##
############

cor.test(T0_qPCR$Mco_Prop_Log10,Mco,method="pearson",alternative = "two.sided")
# p = 4.36e-06, df = 10, t = 8.9479, r = 0.9428514, r^2 = 0.89

plot(T0_qPCR$Mco_Prop_Log10,Mco)
ggplot()+geom_point(aes(x=10^T0_qPCR$Mco_Prop_Log10,y=10^Mco),color="blue",size = 3) +
  scale_x_continuous(trans = "log10", breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  scale_y_continuous(trans = "log10", breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  annotation_custom(grid::textGrob(label = "P < 0.001",
                                   x = unit(0.135, "npc"), y = unit(0.85, "npc"), 
                                   gp = gpar(color = "black", fontsize=11))) +
  annotation_custom(grid::textGrob(label = expression(R^2*" = 0.89"),
                                   x = unit(0.134, "npc"), y = unit(0.91, "npc"), 
                                   gp = gpar(color = "black", fontsize=11))) +
  geom_smooth(aes(x = 10^T0_qPCR$Mco_Prop_Log10, y = 10^(Mco)), method = "lm", se=FALSE,
              color = "black", linewidth=0.4) +
  labs(x = "qPCR Relative Abundance", y = "Amplicon Relative Abundance") + theme_classic()

###############
# Time Series #
###############

cor.test(Oct_2020_qPCR$Avg_Mco_rel_abund,Mco_oct_2020)
#t = 3.018, df = 13, p-value = 0.00989, r = 0.64, r^2 = 0.42

#for displayed data
cor.test(Oct_2020_qPCR$Avg_Mco_rel_abund[Oct_2020_qPCR$Time < 40],
         Mco_oct_2020[Oct_2020_qPCR$Time < 40])
#t = 3.2748, df = 11, p-value = 0.007402, r = 0.0.7026011, r^2 = 0.49



#######

qPCR_plot <- ggplot()+geom_point(aes(x=Oct_2020_qPCR$Time[Oct_2020_qPCR$Time < 40],
                        y=Oct_2020_qPCR$Avg_Mco_rel_abund[Oct_2020_qPCR$Time < 40]),
                        color = "blue") +
  geom_line(aes(x=Oct_2020_qPCR$Time[Oct_2020_qPCR$Time < 40],
                y=Oct_2020_qPCR$Avg_Mco_rel_abund[Oct_2020_qPCR$Time < 40],
                linetype = as.character(Oct_2020_qPCR$Box[Oct_2020_qPCR$Time < 40])),
            color = "blue") + theme_classic() +
  labs(x="Time (Hours)", y = "   qPCR 
  Relative Abundance") +
  guides(linetype = guide_legend(title = "iBag Replicate No.")) 


amplicon_plot <- ggplot() + geom_point(aes(x=Oct_2020_qPCR$Time[Oct_2020_qPCR$Time < 40],
                                  y=Mco_oct_2020[Oct_2020_qPCR$Time < 40]),
                                  color = "blue") +
  geom_line(aes(x=Oct_2020_qPCR$Time[Oct_2020_qPCR$Time < 40],
                y=Mco_oct_2020[Oct_2020_qPCR$Time < 40],
                linetype=as.character(Oct_2020_qPCR$Box[Oct_2020_qPCR$Time < 40])),
            color = "blue") + theme_classic() +
  labs(x="Time (Hours)", y = "   Amplicon 
  Relative Abundance") +
  guides(linetype = guide_legend(title = "iBag Replicate No."))

grid.arrange(qPCR_plot,amplicon_plot)

##########
Mco_abs_abund <- c(1031.143945,3771.62069,745.433083,1205.474825,726.8109836,766.9212874,
                   954.1661588,236.2442683,1753.536149,2509.215458,598.6882094,1783.791633)

plot(log10(Mco_abs_abund),Mco)
cor.test(log10(Mco_abs_abund),Mco) #R^2 = 0.39

plot(log10(Mco_abs_abund),T0_qPCR$Mco_Prop_Log10)
cor.test(log10(Mco_abs_abund),T0_qPCR$Mco_Prop_Log10)


###################
# qPCR vs geochem #
###################

Mco_qPCR <- T0_qPCR$Mco_Prop_Log10
Mcy_qPCR <- T0_qPCR$Mcy_Prop_Log10

###Mco###

# Oxygen 
cor.test(o2,Mco_qPCR, alternative = "two.sided", method = "pearson")
# Negative correlation - Pearson: r = -0.912044, t = -7.0329, p-value = 3.571e-05

# Methane 

cor.test(ch4,Mco_qPCR, alternative = "two.sided", method = "pearson")
## No correlation - Pearson: r = -0.2780195, t = -0.91526, and p-value = 0.3816

# Temp
cor.test(temp,Mco_qPCR, alternative = "two.sided", method = "pearson")
## Negative  correlation - Pearson: r = -0.7933298, t = -4.1208, p-value = 0.002074

###Mcy###

# Oxygen 
cor.test(o2,Mcy_qPCR, alternative = "two.sided", method = "pearson")
# Negative correlation - Pearson: r = -0.867105, t = -5.5047, p-value = 0.0002601

# Methane 

cor.test(ch4,Mcy_qPCR, alternative = "two.sided", method = "pearson")
## No correlation - Pearson: r = 0.02380658, t = 0.075304, and p-value = 0.9415

# Temp
cor.test(temp,Mcy_qPCR, alternative = "two.sided", method = "pearson")
## Negative  correlation - Pearson: r = -0.6818777, t = -2.9479, p-value = 0.01459

