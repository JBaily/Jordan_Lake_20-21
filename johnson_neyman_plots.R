### Johnson-Neyman Plots ###

library(ggplot2)
library(grid)
library(ggtext)
library(gridExtra)
library(scales)
library(survey)
library(interactions)
library(jtools)

##########
## Data ##
##########

#Microbiological Data
All_16S <- read.csv(file = "16S_all.csv",header = TRUE)
T0_16S <- All_16S[All_16S$Time==0,]

Mco <- log10(T0_16S$Methylococcaceae)
Mm <- log10(T0_16S$Methylomonadaceae)
Mp <- log10(T0_16S$Methylophilaceae)
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

#######
# Mco #
#######
mod.cocc.co.int <- lm(Mco ~ ch4 * o2)

#Oxygen
mco_jn_o2_graph <- johnson_neyman(mod.cocc.co.int, pred = o2, modx = ch4, alpha = .05)
mco_jn_o2_graph$plot$labels$x <- expression("CH"[4] ~ "(nM)")
mco_jn_o2_graph$plot$labels$y <-  expression("Slope of O"[2] ~ "(" *mu*"M)")
mco_jn_o2_graph$plot$labels$title <- NULL
mco_jn_o2_graph$plot$theme$legend.position <- "none"
mco_jn_o2_graph_fin <- mco_jn_o2_graph$plot + annotation_custom(grid::textGrob(label = expression(bold("A")), x = unit(0.92, "npc"), 
                                                                               y = unit(0.92, "npc"), gp = gpar(color = "black", 
                                                                                                                fontsize = 17)))

#Methane
mco_jn_ch4_graph <- johnson_neyman(mod.cocc.co.int, pred = ch4, modx = o2, alpha = .05)
mco_jn_ch4_graph$plot$labels$x <- expression("O"[2] ~ "(" *mu*"M)")
mco_jn_ch4_graph$plot$labels$y <- expression("Slope of CH"[4] ~ "(nM)")
mco_jn_ch4_graph$plot$labels$title <- NULL
mco_jn_ch4_graph$plot$theme$legend.text$colour <- "white"
mco_jn_ch4_graph_fin <- mco_jn_ch4_graph$plot + annotation_custom(grid::textGrob(label = expression(bold("B")), x = unit(0.92, "npc"), 
                                                                                 y = unit(0.92, "npc"), 
                                                                                 gp = gpar(color = "black", 
                                                                                           fontsize = 17)))+
  guides(fill = guide_legend(override.aes = list(fill = "white")),
         linetype = guide_legend(override.aes = list(linetype = NA)))

######
# Mm #
######

mod.mona.all.int.co <- lm(Mm ~ ch4 + o2 + temp + o2:ch4)

#Oxygen
mm_jn_o2_graph <- johnson_neyman(mod.mona.all.int.co, pred = o2, modx = ch4, alpha = .05)
mm_jn_o2_graph$plot$labels$x <- expression("CH"[4] ~ "(nM)")
mm_jn_o2_graph$plot$labels$y <-  expression("Slope of O"[2] ~ "(" *mu*"M)")
mm_jn_o2_graph$plot$labels$title <- NULL
mm_jn_o2_graph$plot$theme$legend.position <- "none"
mm_jn_o2_graph_fin <- mm_jn_o2_graph$plot + annotation_custom(grid::textGrob(label = expression(bold("C")), x = unit(0.92, "npc"), 
                                                                             y = unit(0.92, "npc"), gp = gpar(color = "black", 
                                                                                                              fontsize = 17)))

#Methane
mm_jn_ch4_graph <- johnson_neyman(mod.mona.all.int.co, pred = ch4, modx = o2, alpha = .05)
mm_jn_ch4_graph$plot$labels$x <- expression("O"[2] ~ "(" *mu*"M)")
mm_jn_ch4_graph$plot$labels$y <- expression("Slope of CH"[4] ~ "(nM)")
mm_jn_ch4_graph$plot$labels$title <- NULL
mm_jn_ch4_graph_fin <- mm_jn_ch4_graph$plot + annotation_custom(grid::textGrob(label = expression(bold("D")), x = unit(0.92, "npc"), 
                                                                               y = unit(0.92, "npc"), gp = gpar(color = "black", 
                                                                                                                fontsize = 17)))

######
# Mp #
######

mod.phil.to.int <- lm(Mp ~ o2 * temp)

#Oxygen
mp_jn_o2_graph <- johnson_neyman(mod.phil.to.int, pred = o2, modx = temp, alpha = .05)
mp_jn_o2_graph$plot$labels$x <- expression("Temperature "~( degree*C))
mp_jn_o2_graph$plot$labels$y <- expression("Slope of O"[2] ~ "(" *mu*"M)")
mp_jn_o2_graph$plot$labels$title <- NULL
mp_jn_o2_graph$plot$theme$legend.position <- "none"
mp_jn_o2_graph_fin <- mp_jn_o2_graph$plot + 
  annotation_custom(grid::textGrob(label = expression(bold("E")), 
                                   x = unit(0.08, "npc"), y = unit(0.92, "npc"), 
                                   gp = gpar(color = "black", fontsize = 17))) 

#Temperature
mp_jn_temp_graph <- johnson_neyman(mod.phil.to.int, pred = temp, modx = o2, alpha = .05) 
mp_jn_temp_graph$plot$labels$x <- expression("O"[2] ~ "(" *mu*"M)")
mp_jn_temp_graph$plot$labels$y <- expression("Slope of Temperature "~( degree*C))
mp_jn_temp_graph$plot$labels$title <- NULL
mp_jn_temp_graph$plot$theme$legend.text$colour <- "white"
mp_jn_temp_graph_fin <- mp_jn_temp_graph$plot + 
  annotation_custom(grid::textGrob(label = expression(bold("F")), 
                                   x = unit(0.08, "npc"), y = unit(0.92, "npc"), 
                                   gp = gpar(color = "black", fontsize = 17)))+
  guides(fill = guide_legend(override.aes = list(fill = "white")),
         linetype = guide_legend(override.aes = list(linetype = NA)))

######
# Ma #
######

mod.acid.to.int <- lm(Ma ~ o2 * temp)

ma_jn_o2_graph <- johnson_neyman(mod.acid.to.int, pred = o2, modx = temp, alpha = .05)
ma_jn_o2_graph$plot$labels$x <- expression("Temperature "~( degree*C))
ma_jn_o2_graph$plot$labels$y <- expression("Slope of O"[2] ~ "(" *mu*"M)")
ma_jn_o2_graph$plot$labels$title <- NULL
ma_jn_o2_graph$plot$theme$legend.position <- "none"
ma_jn_o2_graph_fin <- ma_jn_o2_graph$plot + annotation_custom(grid::textGrob(label = expression(bold("G")), x = unit(0.08, "npc"), 
                                                                             y = unit(0.92, "npc"), gp = gpar(color = "black", 
                                                                                                              fontsize = 17)))

ma_jn_temp_graph <- johnson_neyman(mod.acid.to.int, pred = temp, modx = o2, alpha = .05) #[183.36, 224.63]
ma_jn_temp_graph$plot$labels$x <- expression("O"[2] ~ "(" *mu*"M)")
ma_jn_temp_graph$plot$labels$y <- expression("Slope of Temperature "~( degree*C))
ma_jn_temp_graph$plot$labels$title <- NULL

ma_jn_temp_graph_fin <- ma_jn_temp_graph$plot + 
  annotation_custom(grid::textGrob(label = expression(bold("H")), 
                                   x = unit(0.08, "npc"), y = unit(0.92, "npc"), 
                                   gp = gpar(color = "black", fontsize = 17)))
#######

grid.arrange(mco_jn_o2_graph_fin,mco_jn_ch4_graph_fin,
             mm_jn_o2_graph_fin,mm_jn_ch4_graph_fin,
             mp_jn_o2_graph_fin,mp_jn_temp_graph_fin,
             ma_jn_o2_graph_fin,ma_jn_temp_graph_fin,
             ncol=2,nrow=4,widths=c(1.0,1.4)) #legend positions to be adjusted in "preview"
