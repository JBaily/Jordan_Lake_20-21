## Model Figures ##
library(survey)
library(interactions)
library(jtools)
library(car)
library(ggplot2)
library(grid)
library(gridExtra)

#Loading in data

#Microbiological Data
All_16S <- read.csv(file = "16S_all.csv",header = TRUE)
T0_16S <- All_16S[All_16S$Time==0,]

Mco <- log10(T0_16S$Methylococcaceae)
Mm <- log10(T0_16S$Methylomonadaceae)
Mcy <- log10(T0_16S$Methylocystaceae)
Ma <- log10(T0_16S$Methylacidiphilaceae)
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

range(temp) #18.276 28.575
range(o2) #64.40 279.03
range(ch4) #195.16 434.36

scale_3_dec <- function(x) sprintf("%.3f", x)

####################
# Methylophilaceae #
####################

# Final model 
mod.phil.to.int <- lm(Mp ~ o2 * temp)
summary(mod.phil.to.int) 
#intercept    6.338e-01  3.831e-01
#o2          -1.485e-02  1.923e-03
#temp        -1.222e-01  1.714e-02  
#o2:temp      6.005e-04  8.098e-05   
Anova(mod.phil.to.int) 


# equation: y = intercept + ax + by + c(xy)
# equation: y = 0.6338 + -0.01485(o2) + -0.1222(temp) + 0.0006005(o2*temp)
# equation: dy/(do2) = -0.01485 + 0.0006005(temp)
# equation: dy/(dtemp) = -0.1222 + 0.0006005(o2)
# equation: intercept (temp) = 0.6338 +  -0.01485(o2)
# equation: intercept (o2)   = 0.6338 + -0.1222(temp)

###############
# Temperature #
###############

interactions::sim_slopes(mod.phil.to.int, pred = temp, modx = o2, 
                         modx.values = c(65,120,180,275), johnson_neyman = FALSE)

#            Est.   S.E.   t val.      p
# for 65:   -0.08   0.01    -6.76   0.00
# for 135:  -0.05   0.01    -5.82   0.00
# for 180:  -0.01   0.01    -2.38   0.03
# for 275:   0.04   0.01     5.12   0.00


Mp_temp_line <- function(temperatures, oxygen){
  
  temperature_intercepts <- vector()
  temperature_slopes <- vector()
  
  temperature_intercepts <- 0.6338 + -0.01485*(oxygen)
  temperature_slopes <- -0.1222 + 0.0006005*(oxygen) 
  
  abundances <- 10^(temperature_intercepts+temperatures*temperature_slopes)
  
  return(abundances)
  
}

Mp_temp_plot_hl <- ggplot()+
  geom_line(aes(x = c(seq(18.3,28.3,0.1)), 
                y = Mp_temp_line(temperatures = c(seq(18.3,28.3,0.1)), oxygen = c(65)),
                linetype = "65"),color="orange2") + 
  geom_line(aes(x = c(seq(18.3,28.3,0.1)), 
                y = Mp_temp_line(temperatures = c(seq(18.3,28.3,0.1)), oxygen = c(120)),
                linetype = "120"),color="orange2") + 
  geom_line(aes(x = c(seq(18.3,28.3,0.1)), 
                  y = Mp_temp_line(temperatures = c(seq(18.3,28.3,0.1)), oxygen = c(180)),
                linetype = "180"),color="orange2") + 
  geom_line(aes(x = c(seq(18.3,28.3,0.1)), 
                y = Mp_temp_line(temperatures = c(seq(18.3,28.3,0.1)), oxygen = c(275)),
                linetype = "275"),color="orange2") + 
  annotation_custom(grid::textGrob(label = expression(bold("D")),
                                   x = unit(0.96, "npc"), y = unit(0.96, "npc"), 
                                   gp = gpar(color = "black", fontsize=17))) +
  scale_linetype_manual(breaks = c("65","120","180","275"),values = c(4,6,5,1))+
  guides(linetype = guide_legend(title = expression("O"[2] ~ "(" *mu*"M)"))) +
  scale_x_continuous(limits=c(18,29), breaks = c(18,20.75,23.50,26.25,29)) + 
  labs(x = expression("Temperature "~( degree*C)), 
       y = "Relative Abundance of"~italic("Methylophilaceae")) +
  theme_classic() + theme(axis.title.y = element_blank()) +
  scale_y_continuous(limits = c(0,0.015)) 

##########
# Oxygen # 
##########

interactions::sim_slopes(mod.phil.to.int, pred = o2, modx = temp, 
                         modx.values = c(18.3,20.8,23.3,25.8,28.3), johnson_neyman = FALSE)

#            Est.   S.E.   t val.      p
# for 18.3: -0.00   0.00    -7.44   0.00
# for 20.8: -0.00   0.00    -6.27   0.00
# for 23.3: -0.00   0.00    -2.76   0.02
# for 25.8:  0.00   0.00     1.78   0.10 -- NS
# for 28.3:  0.00   0.00     4.29   0.00

Mp_o2_line <- function(oxygen, temperatures){
  
  oxygen_intercepts <- vector()
  oxygen_slopes <- vector()
  
  oxygen_intercepts <- 0.6338 + -0.1222*(temperatures)
  oxygen_slopes <- -0.01485 + 0.0006005*(temperatures)
  
  abundances <- 10^(oxygen_intercepts+oxygen*oxygen_slopes)
  
  return(abundances)
  
}


Mp_o2_plot_hl <- ggplot()+
  geom_line(aes(x = c(seq(65,280,2.15)), 
                y = Mp_o2_line(oxygen = c(seq(65,280,2.15)), temperatures = c(18.3)),
                linetype = "18.3"),color="orange2") + 
  geom_line(aes(x = c(seq(65,280,2.15)), 
                y = Mp_o2_line(oxygen = c(seq(65,280,2.15)), temperatures = c(20.8)),
                linetype = "20.8"),color="orange2") + 
  geom_line(aes(x = c(seq(65,280,2.15)), 
                y = Mp_o2_line(oxygen = c(seq(65,280,2.15)), temperatures = c(23.3)),
                linetype = "23.3"),color="orange2") + 
  geom_line(aes(x = c(seq(65,280,2.15)), 
                y = Mp_o2_line(oxygen = c(seq(65,280,2.15)), temperatures = c(28.3)),
                linetype = "28.3"),color="orange2") + 
  annotation_custom(grid::textGrob(label = expression(bold("C")),
                                   x = unit(0.96, "npc"), y = unit(0.96, "npc"), 
                                   gp = gpar(color = "black", fontsize=17))) +
  scale_linetype_manual(breaks = c("18.3","20.8","23.3","28.3"),values = c(4,6,5,1))+
  guides(linetype = guide_legend(title = expression("Temperature "~(degree*C)))) +
  scale_x_continuous(limits=c(60,285)) + theme_classic() + 
  labs(x = expression("O"[2] ~ "(" *mu*"M)"), 
       y = "Relative Abundance of"~italic("Methylophilaceae")) +
  scale_y_continuous(limits = c(0,0.015)) + theme(legend.box.spacing = unit(0, "pt"),
                                                  legend.title.align = 0.5,
                                                  legend.box.just = "center")

##############
grid.arrange(Mp_o2_plot_hl,Mp_temp_plot_hl, nrow=1,ncol=2,widths=c(1.2,1))

########################
# Methylacidiphilaceae #
########################

# Final model 
mod.acid.to.int <- lm(Ma ~ o2 * temp)
summary(mod.acid.to.int) 
#intercept   -1.768e+00  4.225e-01
#o2          -5.376e-03  2.121e-03
#temp         1.010e-02  1.891e-02  
#o2:temp      2.121e-04  8.932e-05   
Anova(mod.acid.to.int) 


# equation: y = intercept + ax + by + c(xy)
# equation: y = -1.768 + -0.005376(o2) + 0.0101(temp) + 0.0002121(o2*temp)
# equation: dy/(do2) = -0.005376 + 0.0002121(temp)
# equation: dy/(dtemp) = 0.0101 + 0.0002121(o2)
# equation: intercept (temp) = -1.768 + -0.005376(o2)
# equation: intercept (o2)   = -1.768 + 0.0101(temp)

###############
# Temperature #
###############

interactions::sim_slopes(mod.acid.to.int, pred = temp, modx = o2, 
                         modx.values = c(95,155,215,275), johnson_neyman = FALSE)

#            Est.   S.E.   t val.      p
# for 95:    0.03   0.01     2.68   0.02
# for 155:   0.04   0.01     5.76   0.00
# for 215:   0.06   0.01     8.61   0.00
# for 275:   0.07   0.01     7.40   0.00


Ma_temp_line <- function(temperatures, oxygen){
  
  temperature_intercepts <- vector()
  temperature_slopes <- vector()
  
  temperature_intercepts <- -1.768 + -0.005376*(oxygen)
  temperature_slopes <- 0.0101 + 0.0002121*(oxygen) 
  
  abundances <- 10^(temperature_intercepts+temperatures*temperature_slopes)
  
  return(abundances)
  
}

Ma_temp_plot_hl <- ggplot()+
  geom_line(aes(x = c(seq(18.3,28.3,0.1)), 
                y = Ma_temp_line(temperatures = c(seq(18.3,28.3,0.1)), oxygen = c(95)),
                linetype = "95"),color="hotpink1") + 
  geom_line(aes(x = c(seq(18.3,28.3,0.1)), 
                y = Ma_temp_line(temperatures = c(seq(18.3,28.3,0.1)), oxygen = c(155)),
                linetype = "155"),color="hotpink1") + 
  geom_line(aes(x = c(seq(18.3,28.3,0.1)), 
                y = Ma_temp_line(temperatures = c(seq(18.3,28.3,0.1)), oxygen = c(215)),
                linetype = "215"),color="hotpink1") + 
  geom_line(aes(x = c(seq(18.3,28.3,0.1)), 
                y = Ma_temp_line(temperatures = c(seq(18.3,28.3,0.1)), oxygen = c(275)),
                linetype = "275"),color="hotpink1") + 
  annotation_custom(grid::textGrob(label = expression(bold("C")),
                                   x = unit(0.96, "npc"), y = unit(0.96, "npc"), 
                                   gp = gpar(color = "black", fontsize=17))) +
  scale_linetype_manual(breaks = c("95","155","215","275"),values = c(4,6,5,1))+
  guides(linetype = guide_legend(title = expression("O"[2] ~ "(" *mu*"M)"))) +
  scale_x_continuous(limits=c(18,29), breaks = c(18,20.75,23.50,26.25,29)) + 
  labs(x = expression("Temperature "~( degree*C)), 
       y = "Relative Abundance of"~italic("Methylacidiphilaceae")) +
  theme_classic() + 
  scale_y_continuous(limits = c(0,0.05)) 

##########
# Oxygen # 
##########

interactions::sim_slopes(mod.acid.to.int, pred = o2, modx = temp, 
                         modx.values = c(18.3,19.17,20.04,20.90), johnson_neyman = FALSE)

#            Est.   S.E.   t val.      p
# for 18.30: -0.00   0.00    -2.62   0.02
# for 19.17: -0.00   0.00    -2.56   0.02
# for 20.04: -0.00   0.00    -2.47   0.03
# for 20.90: -0.00   0.00    -2.31   0.04

Ma_o2_line <- function(oxygen, temperatures){
  
  oxygen_intercepts <- vector()
  oxygen_slopes <- vector()
  
  oxygen_intercepts <- -1.768 + 0.0101*(temperatures)
  oxygen_slopes <- -0.005376 + 0.0002121*(temperatures)
  
  abundances <- 10^(oxygen_intercepts+oxygen*oxygen_slopes)
  
  return(abundances)
  
}


Ma_o2_plot_hl <- ggplot()+
  geom_line(aes(x = c(seq(65,280,2.15)), 
                y = Ma_o2_line(oxygen = c(seq(65,280,2.15)), temperatures = c(18.3)),
                linetype = "18.3"),color="hotpink1") + 
  geom_line(aes(x = c(seq(65,280,2.15)), 
                y = Ma_o2_line(oxygen = c(seq(65,280,2.15)), temperatures = c(19.17)),
                linetype = "19.17"),color="hotpink1") + 
  geom_line(aes(x = c(seq(65,280,2.15)), 
                y = Ma_o2_line(oxygen = c(seq(65,280,2.15)), temperatures = c(20.04)),
                linetype = "20.04"),color="hotpink1") + 
  geom_line(aes(x = c(seq(65,280,2.15)), 
                y = Ma_o2_line(oxygen = c(seq(65,280,2.15)), temperatures = c(20.90)),
                linetype = "20.90"),color="hotpink1") + 
#  annotation_custom(grid::textGrob(label = expression(bold("C")),
#                                   x = unit(0.96, "npc"), y = unit(0.96, "npc"), 
#                                   gp = gpar(color = "black", fontsize=17))) +
  scale_linetype_manual(breaks = c("18.3","19.17","20.04","20.90"),values = c(4,6,5,1))+
  guides(linetype = guide_legend(title = expression("Temperature "~(degree*C)))) +
  scale_x_continuous(limits=c(60,285)) + theme_classic() + 
  labs(x = expression("O"[2] ~ "(" *mu*"M)"), 
       y = "Relative Abundance of"~italic("Methylacidiphilaceae")) +
  scale_y_continuous(limits = c(0,0.05)) + theme(legend.box.spacing = unit(0, "pt"),
                                                  legend.title.align = 0.5,
                                                  legend.box.just = "center")

##############
grid.arrange(Ma_o2_plot_hl,Ma_temp_plot_hl, nrow=1,ncol=2,widths=c(1.2,1))

####################
# Methylococcaceae #
####################

# Final model 
mod.cocc.co.int <- lm(Mco ~ ch4 * o2)
summary(mod.cocc.co.int) 
#(Intercept) -6.009e+00  7.807e-01
#ch4          1.214e-02  2.337e-03
#o2           1.936e-02  4.355e-03
#ch4:o2      -7.089e-05  1.265e-05
Anova(mod.cocc.co.int) 


# equation: y = intercept + ax + by + c(xy)
# equation: y = -6.009 + 0.01936(o2) + 0.01214(ch4) + -0.00007089(o2*ch4)
# equation: dy/(do2) = 0.01936 + -0.00007089(ch4)
# equation: dy/(dch4) = 0.01214 + -0.00007089(o2)
# equation: intercept (ch4) = -6.009 + 0.01936(o2)
# equation: intercept (o2)  = -6.009 + 0.01214(ch4)

##########
# Oxygen # 
##########

interactions::sim_slopes(mod.cocc.co.int, pred = o2, modx = ch4, 
                         modx.values = c(200,260,320,380,430), johnson_neyman = FALSE)

#            Est.   S.E.   t val.      p
# for 200:  0.01   0.00     2.77   0.02
# for 260:  0.00   0.00     0.79   0.44 - N.S. 
# for 320: -0.00   0.00    -5.58   0.00
# for 380: -0.01   0.00   -10.55   0.00
# for 430: -0.01   0.00    -9.00   0.00

Mco_o2_line <- function(oxygen, methane){
  
  oxygen_intercepts <- vector()
  oxygen_slopes <- vector()
  
  oxygen_intercepts <-  -6.009 + 0.01214*(methane)
  oxygen_slopes <-  0.01936 + -0.00007089*(methane) 
  
  abundances <- 10^(oxygen_intercepts+oxygen*oxygen_slopes)
  
  return(abundances)
  
}


Mco_o2_plot_hl <- ggplot()+
  geom_line(aes(x = c(seq(65,280,2.15)), 
                y = Mco_o2_line(oxygen = c(seq(65,280,2.15)), methane = c(200)),
                linetype = "200"),color="blue") + 
  geom_line(aes(x = c(seq(65,280,2.15)), 
                y = Mco_o2_line(oxygen = c(seq(65,280,2.15)), methane = c(320)),
                linetype = "320"),color="blue") + 
  geom_line(aes(x = c(seq(65,280,2.15)), 
                y = Mco_o2_line(oxygen = c(seq(65,280,2.15)), methane = c(380)),
                linetype = "380"),color="blue") + 
  geom_line(aes(x = c(seq(65,280,2.15)), 
                y = Mco_o2_line(oxygen = c(seq(65,280,2.15)), methane = c(430)),
                linetype = "430"),color="blue") + 
  annotation_custom(grid::textGrob(label = expression(bold("A")),
                                   x = unit(0.96, "npc"), y = unit(0.96, "npc"), 
                                   gp = gpar(color = "black", fontsize=17))) +
  scale_linetype_manual(breaks = c("200","320","380","430"),values = c(4,6,5,1))+
  guides(linetype = guide_legend(title = expression("CH"[4] ~ "(nM)"))) + theme_classic() +
  scale_x_continuous(limits=c(60,285)) + scale_y_continuous(labels = scale_3_dec) +
  labs(x = expression("O"[2] ~ "(" *mu*"M)"), 
       y = "Relative Abundance of"~italic("Methylococcaceae")) + 
  theme(legend.box.margin =  margin(0,17,0,15))

###########
# Methane # 
###########

interactions::sim_slopes(mod.cocc.co.int, pred = ch4, modx = o2, 
                         modx.values = c(65,135,205,275), johnson_neyman = FALSE)

#            Est.   S.E.   t val.      p
# for 65:   0.01   0.00     4.77   0.00
# for 135:  0.00   0.00     2.86   0.01
# for 205: -0.00   0.00    -2.87   0.01
# for 275: -0.01   0.00    -5.02   0.00

Mco_ch4_line <- function(methane, oxygen){
  
  methane_intercepts <- vector()
  methane_slopes <- vector()
  
  methane_intercepts <-  -6.009 + 0.01936*(oxygen)
  methane_slopes <-  0.01214 + -0.00007089*(oxygen) 
  
  abundances <- 10^(methane_intercepts+methane*methane_slopes)
  
  return(abundances)
  
}


Mco_ch4_plot_hl <- ggplot()+
  geom_line(aes(x = c(seq(195,435,2.4)), 
                y = Mco_ch4_line(methane = c(seq(195,435,2.4)), oxygen = c(65)),
                linetype = "65"),color="blue") + 
  geom_line(aes(x = c(seq(195,435,2.4)), 
                y = Mco_ch4_line(methane = c(seq(195,435,2.4)), oxygen = c(135)),
                linetype = "135"),color="blue") + 
  geom_line(aes(x = c(seq(195,435,2.4)), 
                y = Mco_ch4_line(methane = c(seq(195,435,2.4)), oxygen = c(205)),
                linetype = "205"),color="blue") + 
  geom_line(aes(x = c(seq(195,435,2.4)), 
                y = Mco_ch4_line(methane = c(seq(195,435,2.4)), oxygen = c(275)),
                linetype = "275"),color="blue") + 
  annotation_custom(grid::textGrob(label = expression(bold("B")),
                                   x = unit(0.96, "npc"), y = unit(0.96, "npc"), 
                                   gp = gpar(color = "black", fontsize=17))) +
  scale_linetype_manual(breaks = c("65","135","205","275"),values = c(4,6,5,1))+
  guides(linetype = guide_legend(title = expression("O"[2] ~ "(" *mu*"M)"))) +
  scale_x_continuous(limits=c(195,445)) +  theme_classic() + scale_y_continuous(labels = scale_3_dec) +
  labs(x = expression("CH"[4] ~ "(nM)"), 
       y = "Relative Abundance of"~italic("Methylococcaceae")) + 
  theme(axis.title.y = element_blank())


#####################
grid.arrange(Mco_o2_plot_hl,Mco_ch4_plot_hl, nrow=1,ncol=2, widths=c(1.07,1))


#####################
# Methylomonadaceae #
#####################

# Final model 
mod.mona.all.int.co <- lm(Mm ~ ch4 + o2 + temp + o2:ch4)
summary(mod.mona.all.int.co) 
#(Intercept) -1.144e+01  2.068e+00
#ch4          2.265e-02  5.315e-03
#o2           3.492e-02  8.740e-03
#temp         8.202e-02  2.827e-02
#ch4:o2      -1.225e-04  2.640e-05


# equation: y = intercept + ax + bz + c(xz) + e(v)
# equation: y = -11.44 + 0.03492(o2) + 0.02265(ch4) + -0.0001225(o2*ch4) + 0.08202(temp)
# equation: dy/(do2) = 0.03492 + -0.0001225(ch4)
# equation: dy/(dch4) = 0.02265 + -0.0001225(o2) 
# equation: intercept (o2)  = -11.44 + 0.02265(ch4) + 0.08202(temp)
# equation: intercept (ch4) = -11.44 + 0.03492(o2) + 0.08202(temp)

##########
# Oxygen # 
##########

interactions::sim_slopes(mod.mona.all.int.co, pred = o2, modx = ch4, 
                         modx.values = c(200,260,320,380,430), johnson_neyman = FALSE)

#            Est.   S.E.   t val.      p
# for 200:  0.01   0.00     2.87   0.02
# for 260:  0.00   0.00     1.36   0.20 - N.S. 
# for 320: -0.00   0.00    -2.99   0.01
# for 380: -0.01   0.00    -5.82   0.00
# for 430: -0.02   0.00    -5.77   0.00

Mm_o2_line <- function(oxygen, methane, temperature){
  
  oxygen_intercepts <- vector()
  oxygen_slopes <- vector()
  
  oxygen_intercepts <-  -11.144 + 0.02265*(methane) + 0.08202*(temperature)
  oxygen_slopes <-  0.03492 + -0.0001225*(methane) 
  
  abundances <- 10^(oxygen_intercepts+oxygen*oxygen_slopes)
  
  return(abundances)
  
}



Mm_o2_plot_fall_limited <- ggplot()+
  geom_line(aes(x = c(seq(65,280,2.15)), 
                y = Mm_o2_line(oxygen = c(seq(65,280,2.15)), methane = c(200), 
                               temperature = c(mean(temp[season=="O"]))),
                linetype = "200"),color="purple4") + 
  geom_line(aes(x = c(seq(65,280,2.15)), 
                y = Mm_o2_line(oxygen = c(seq(65,280,2.15)), methane = c(320),
                               temperature = mean(temp[season=="O"])),
                linetype = "320"),color="purple4") + 
  geom_line(aes(x = c(seq(65,280,2.15)), 
                y = Mm_o2_line(oxygen = c(seq(65,280,2.15)), methane = c(380),
                               temperature = mean(temp[season=="O"])),
                linetype = "380"),color="purple4") + 
  geom_line(aes(x = c(seq(65,280,2.15)), 
                y = Mm_o2_line(oxygen = c(seq(65,280,2.15)), methane = c(430),
                               temperature = mean(temp[season=="O"])),
                linetype = "430"),color="purple4") + 
  annotation_custom(grid::textGrob(label = expression(bold("A")),
                                   x = unit(0.96, "npc"), y = unit(0.96, "npc"), 
                                   gp = gpar(color = "black", fontsize=17))) +
  scale_linetype_manual(breaks = c("200","320","380","430"),values = c(4,6,5,1))+
  guides(linetype = guide_legend(title = expression("CH"[4] ~ "(nM)"))) +
  scale_x_continuous(limits=c(60,280)) + 
  labs(x = expression("O"[2] ~ "(" *mu*"M)"), 
       y = "Relative Abundance of"~italic("Methylomonadaceae")) + theme_classic() 

###########
# Methane # 
###########

interactions::sim_slopes(mod.mona.all.int.co, pred = ch4, modx = o2, 
                         modx.values = c(65,117.5,170,222.5,275), johnson_neyman = FALSE)

#             Est.   S.E.   t val.      p
# for 65:     0.01   0.00     3.93   0.00
# for 117.5:  0.01   0.00     3.18   0.01
# for 170:    0.00   0.00     1.00   0.34
# for 222.5: -0.00   0.00    -2.37   0.04
# for 275:   -0.01   0.00    -3.89   0.00

Mm_ch4_line <- function(methane, oxygen, temperature){
  
  methane_intercepts <- vector()
  methane_slopes <- vector()
  
  methane_intercepts <-  -11.44 + 0.03492*(oxygen) + 0.08202*(temperature)
  methane_slopes <-  0.02265 + -0.0001225*(oxygen)
  
  abundances <- 10^(methane_intercepts+methane*methane_slopes)
  
  return(abundances)
  
}


Mm_ch4_plot_fall_limited <- ggplot()+
  geom_line(aes(x = c(seq(190,440,2.5)), 
                y = Mm_ch4_line(methane = c(seq(190,440,2.5)), oxygen = c(65), 
                               temperature = c(mean(temp[season=="O"]))),
                linetype = "65"),color="purple4") + 
  geom_line(aes(x = c(seq(190,440,2.5)), 
                y = Mm_ch4_line(methane = c(seq(190,440,2.5)), oxygen = c(117.5),
                               temperature = mean(temp[season=="O"])),
                linetype = "117.5"),color="purple4") + 
  geom_line(aes(x = c(seq(190,440,2.5)), 
                y = Mm_ch4_line(methane = c(seq(190,440,2.5)), oxygen = c(222.5),
                               temperature = mean(temp[season=="O"])),
                linetype = "222.5"),color="purple4") + 
  geom_line(aes(x = c(seq(190,440,2.5)), 
                y = Mm_ch4_line(methane = c(seq(190,440,2.5)), oxygen = c(275),
                               temperature = mean(temp[season=="O"])),
                linetype = "275"),color="purple4") + 
  annotation_custom(grid::textGrob(label = expression(bold("B")),
                                   x = unit(0.96, "npc"), y = unit(0.96, "npc"), 
                                   gp = gpar(color = "black", fontsize=17))) +
  scale_linetype_manual(breaks = c("65","117.5","222.5","275"),values = c(4,6,5,1))+
  guides(linetype = guide_legend(title = expression("O"[2] ~ "(" *mu*"M)"))) +
  scale_x_continuous(limits=c(185,450)) + theme_classic() +
  labs(x = expression("CH"[4] ~ "(nM)"), 
       y = "Relative Abundance of"~italic("Methylomonadaceae")) + 
  theme(axis.title.y = element_text(color = "white"))

###############
grid.arrange(Mm_o2_plot_fall_limited,Mm_ch4_plot_fall_limited, 
             nrow=1,ncol=2, widths=c(1,1))

####################

#main body of paper 
grid.arrange(Mco_o2_plot_hl,Mco_ch4_plot_hl,
             Mp_o2_plot_hl,Mp_temp_plot_hl,
             nrow=2,ncol=2,widths=c(1.2,1))

#supplemental body of paper 
grid.arrange(Mm_o2_plot_fall_limited,Mm_ch4_plot_fall_limited, Ma_temp_plot_hl,
             nrow=1,ncol=3, widths=c(1,1,0.96))




