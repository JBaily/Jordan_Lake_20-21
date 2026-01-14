## Correlations mega-figure ##

library(ggplot2)
library(grid)
library(ggtext)
library(gridExtra)
library(scales)
library(cowplot)

##########
## Data ##
##########

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

#########################
## Oxygen Correlations ##
#########################

# Methylococcaceae #

# Methylococcaceae #

o2_Mco <- ggplot() + geom_point(size=2.3, aes(x=o2,y=10^(Mco)), col = "blue") + 
  scale_y_continuous(trans = "log10", breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  annotation_custom(grid::textGrob(label = expression(bold("A")),
                                   x = unit(0.96, "npc"), y = unit(0.915, "npc"), 
                                   gp = gpar(color = "black", fontsize=17))) +
  annotation_custom(grid::textGrob(label = "P < 0.001",
                                   x = unit(0.178, "npc"), y = unit(0.11, "npc"), 
                                   gp = gpar(color = "black", fontsize=12))) +
  annotation_custom(grid::textGrob(label = expression(R^2*" = 0.69"),
                                   x = unit(0.177, "npc"), y = unit(0.25, "npc"), 
                                   gp = gpar(color = "black", fontsize=12))) +
  geom_smooth(aes(x = o2, y = 10^(Mco)), method = "lm", se=FALSE,
              color = "black", linewidth=0.5) +xlim(c(60,300))+
  theme(axis.title.y = element_blank(), axis.line.y = element_line(color="black"),
        panel.background = element_rect(fill = "white", color = "white"), 
        axis.line.x.bottom = element_line(color="black"), axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12))


# Methylomonadaceae #

o2_Mm <- ggplot() + geom_point(size=2.3, aes(x=o2,y=10^(Mm)), col = "purple4") + 
  annotation_custom(grid::textGrob(label = expression(bold("D")),
                                   x = unit(0.96, "npc"), y = unit(0.915, "npc"), 
                                   gp = gpar(color = "black", fontsize=17))) +
  annotation_custom(grid::textGrob(label = "P = 0.017",
                                   x = unit(0.178, "npc"), y = unit(0.11, "npc"), 
                                   gp = gpar(color = "black", fontsize=12))) +
  annotation_custom(grid::textGrob(label = expression(R^2*" = 0.35"),
                                   x = unit(0.177, "npc"), y = unit(0.25, "npc"), 
                                   gp = gpar(color = "black", fontsize=12))) +
  scale_y_continuous(trans = "log10", breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  geom_smooth(aes(x = o2, y = 10^(Mm)), method = "lm", se=FALSE,
              color = "black", linewidth=0.5) + xlim(c(60,300))+
  theme(axis.title.y = element_blank(), axis.line.y = element_line(color="black",),
        panel.background = element_rect(fill = "white", color = "white"), 
        axis.line.x.bottom = element_line(color="black"), axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12))

# Methylocystaceae #

o2_Mcy <- ggplot() + geom_point(size=2.3, aes(x=o2,y=10^Mcy), col = "green3", alpha = 0.7, shape = 1, stroke = 1) + 
  annotation_custom(grid::textGrob(label = expression(bold("G")),
                                   x = unit(0.96, "npc"), y = unit(0.915, "npc"), 
                                   gp = gpar(color = "black", fontsize=17))) +
  scale_y_continuous(trans = "log10", breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + xlim(c(60,300))+ 
  theme(axis.title.y = element_blank(), axis.line.y = element_line(color="black"),
        panel.background = element_rect(fill = "white", color = "white"), 
        axis.line.x.bottom = element_line(color="black"), axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=12))

# Methylacidiphilaceae #

o2_Ma <- ggplot() + geom_point(size=2.3, aes(x=o2,y=10^Ma), col = "hotpink1",alpha = 0.7, shape = 1, stroke = 1) + 
  annotation_custom(grid::textGrob(label = expression(bold("M")),
                                   x = unit(0.96, "npc"), y = unit(0.915, "npc"), 
                                   gp = gpar(color = "black", fontsize=17))) +
  scale_y_continuous(trans = "log10", breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  annotate("text",x = 285, y = 0.01, label= "X",color="white") + xlim(c(60,300))+
  theme(axis.title.y = element_blank(), 
        axis.line.y = element_line(color="black"),
        panel.background = element_rect(fill = "white", color = "white"), 
        axis.line.x.bottom = element_line(color="black"), 
        axis.title.x = element_blank(),
        axis.text = element_text(size=12))

# Methylophilaceae # 

o2_Mp <- ggplot() + geom_point(size=2.3, aes(x=o2,y=10^Mp), col = "orange",alpha = 0.7, shape = 1, stroke = 1) + 
  annotation_custom(grid::textGrob(label = expression(bold("J")),
                                   x = unit(0.96, "npc"), y = unit(0.915, "npc"), 
                                   gp = gpar(color = "black", fontsize=17))) +
  scale_y_continuous(trans = "log10", breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)),limits=c(0.0023,0.015)) + xlim(c(60,300))+
  theme(axis.title.y = element_blank(), axis.line.y = element_line(color="black"),
        panel.background = element_rect(fill = "white", color = "white"), 
        axis.line.x.bottom = element_line(color="black"), axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=12))

##########################
## Methane Correlations ##
##########################

# Methylococcaceae #

ch4_Mco <- ggplot() + geom_point(size=2.3, aes(x=ch4,y=10^Mco), col = "blue",alpha = 0.7, shape = 1, stroke = 1) + 
  theme_classic() + 
  scale_y_continuous(trans = "log10", breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  scale_x_continuous(limits = c(195,450)) +
  annotation_custom(grid::textGrob(label = expression(bold("B")),
                                   x = unit(0.96, "npc"), y = unit(0.915, "npc"), 
                                   gp = gpar(color = "black", fontsize=17))) +
  theme(axis.title.y = element_blank(), axis.line.y = element_line(color="black",),
        panel.background = element_rect(fill = "white", color = "white"), 
        axis.line.x.bottom = element_line(color="black"), axis.title.x = element_blank(), 
        axis.text = element_blank())

# Methylomonadaceae #

ch4_Mm <- ggplot() + geom_point(size=2.3, aes(x=ch4,y=10^Mm), col = "purple4",alpha = 0.7, shape = 1, stroke = 1) +
  theme_classic() +   
  annotation_custom(grid::textGrob(label = expression(bold("E")),
                                   x = unit(0.96, "npc"), y = unit(0.915, "npc"), 
                                   gp = gpar(color = "black", fontsize=17))) +
  scale_y_continuous(trans = "log10", breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  scale_x_continuous(limits = c(195,450)) +
  theme(axis.title.y = element_blank(), axis.line.y = element_line(color="black",),
        panel.background = element_rect(fill = "white", color = "white"), 
        axis.line.x.bottom = element_line(color="black"), axis.title.x = element_blank(), 
        axis.text = element_blank())

# Methylocystaceae #

ch4_Mcy <- ggplot() + geom_point(size=2.3, aes(x=ch4,y=10^Mcy), col = "green3",alpha = 0.7, shape = 1, stroke = 1) + 
  scale_y_continuous(trans = "log10", breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  annotation_custom(grid::textGrob(label = expression(bold("H")),
                                   x = unit(0.96, "npc"), y = unit(0.915, "npc"), 
                                   gp = gpar(color = "black", fontsize=17))) +
  scale_x_continuous(limits = c(195,450)) + theme_classic() +   
  theme(axis.title.y = element_blank(), axis.line.y = element_line(color="black",),
        panel.background = element_rect(fill = "white", color = "white"), 
        axis.line.x.bottom = element_line(color="black"), axis.title.x = element_blank(), 
        axis.text = element_blank())

# Methylacidiphilaceae #

ch4_Ma <- ggplot() + geom_point(size=2.3, aes(x=ch4,y=10^Ma), col = "hotpink1",alpha = 0.7, shape = 1, stroke = 1) + 
  annotation_custom(grid::textGrob(label = expression(bold("N")),
                                   x = unit(0.96, "npc"), y = unit(0.915, "npc"), 
                                   gp = gpar(color = "black", fontsize=17))) +
  scale_y_continuous(trans = "log10", breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  scale_x_continuous(limits = c(195,450)) +
  theme(axis.title.y = element_blank(), axis.line.y = element_line(color="black"),
        panel.background = element_rect(fill = "white", color = "white"), 
        axis.line.x.bottom = element_line(color="black"), axis.title.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=12))

# Methylophilaceae # 

ch4_Mp <- ggplot() + geom_point(size=2.3, aes(x=ch4,y=10^Mp), col = "orange",alpha = 0.7, shape = 1, stroke = 1) + 
  annotation_custom(grid::textGrob(label = expression(bold("K")),
                                   x = unit(0.96, "npc"), y = unit(0.915, "npc"), 
                                   gp = gpar(color = "black", fontsize=17))) +
  scale_y_continuous(trans = "log10", breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)),limits=c(0.0023,0.015)) + theme_classic() + 
  scale_x_continuous(limits = c(195,450)) +
  theme(axis.title.y = element_blank(), axis.line.y = element_line(color="black",),
        panel.background = element_rect(fill = "white", color = "white"), 
        axis.line.x.bottom = element_line(color="black"), axis.title.x = element_blank(), 
        axis.text = element_blank())

##############################
## Temperature Correlations ##
##############################

# Methylococcaceae #

temp_Mco <- ggplot() + geom_point(size=2.3, aes(x=temp,y=10^Mco), col = "blue") + 
  annotation_custom(grid::textGrob(label = expression(bold("C")),
                                   x = unit(0.96, "npc"), y = unit(0.915, "npc"), 
                                   gp = gpar(color = "black", fontsize=17))) +
  annotation_custom(grid::textGrob(label = "P = 0.025",
                                   x = unit(0.178, "npc"), y = unit(0.11, "npc"), 
                                   gp = gpar(color = "black", fontsize=12))) +
  annotation_custom(grid::textGrob(label = expression(R^2*" = 0.31"),
                                   x = unit(0.177, "npc"), y = unit(0.25, "npc"), 
                                   gp = gpar(color = "black", fontsize=12))) +
  scale_y_continuous(trans = "log10", breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  scale_x_continuous(limits = c(18.2,29.2), breaks = c(18,20.75, 23.5,26.25,29)) +
  geom_smooth(aes(x = temp, y = 10^(Mco)), method = "lm", se=FALSE,
              color = "black", linewidth=0.5) + 
  theme(axis.title.y = element_blank(), axis.line.y = element_line(color="black",),
        panel.background = element_rect(fill = "white", color = "white"), 
        axis.line.x.bottom = element_line(color="black"), axis.title.x = element_blank(), 
        axis.text = element_blank()) + 
  guides(linetype = "none", color = "none") 

# Methylomonadaceae #

temp_Mm <- ggplot() + geom_point(size=2.3, aes(x=temp,y=10^Mm), col = "purple4",alpha = 0.7, shape = 1, stroke = 1) + 
  annotation_custom(grid::textGrob(label = expression(bold("F")),
                                   x = unit(0.96, "npc"), y = unit(0.915, "npc"), 
                                   gp = gpar(color = "black", fontsize=17))) +
  scale_y_continuous(trans = "log10", breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  scale_x_continuous(limits = c(18.2,29.2), breaks = c(18,20.75, 23.5,26.25,29)) +
  theme(axis.title.y = element_blank(), axis.line.y = element_line(color="black",),
        panel.background = element_rect(fill = "white", color = "white"), 
        axis.line.x.bottom = element_line(color="black"), axis.title.x = element_blank(), 
        axis.text = element_blank()) + 
  guides(linetype = "none", color = "none") 

# Methylocystaceae #

temp_Mcy <- ggplot() + geom_point(size=2.3, aes(x=temp,y=10^Mcy), col = "green3") + 
  annotation_custom(grid::textGrob(label = expression(bold("I")),
                                   x = unit(0.96, "npc"), y = unit(0.915, "npc"), 
                                   gp = gpar(color = "black", fontsize=17))) +
  annotation_custom(grid::textGrob(label = "P < 0.001",
                                   x = unit(0.178, "npc"), y = unit(0.11, "npc"), 
                                   gp = gpar(color = "black", fontsize=12))) +
  annotation_custom(grid::textGrob(label = expression(R^2*" = 0.84"),
                                   x = unit(0.177, "npc"), y = unit(0.25, "npc"), 
                                   gp = gpar(color = "black", fontsize=12))) +
  scale_y_continuous(trans = "log10", breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  scale_x_continuous(limits = c(18.2,29.2), breaks = c(18,20.75, 23.5,26.25,29)) +
  geom_smooth(aes(x = temp, y = 10^(Mcy)), method = "lm", se=FALSE,
              color = "black", linewidth=0.5) + 
  theme(axis.title.y = element_blank(), axis.line.y = element_line(color="black",),
        panel.background = element_rect(fill = "white", color = "white"), 
        axis.line.x.bottom = element_line(color="black"), axis.title.x = element_blank(), 
        axis.text = element_blank()) + 
  guides(linetype = "none", color = "none") 

# Methylacidiphilaceae # 

temp_Ma <- ggplot() + geom_point(size=2.3, aes(x=temp,y=10^Ma), col = "hotpink1") + 
  annotation_custom(grid::textGrob(label = expression(bold("O")),
                                   x = unit(0.96, "npc"), y = unit(0.915, "npc"), 
                                   gp = gpar(color = "black", fontsize=17))) +
  annotation_custom(grid::textGrob(label = "P < 0.001",
                                   x = unit(0.178, "npc"), y = unit(0.73, "npc"), 
                                   gp = gpar(color = "black", fontsize=12))) +
  annotation_custom(grid::textGrob(label = expression(R^2*" = 0.81"),
                                   x = unit(0.177, "npc"), y = unit(0.88, "npc"), 
                                   gp = gpar(color = "black", fontsize=12))) +
  scale_y_continuous(trans = "log10", breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) + 
  geom_smooth(aes(x = temp, y = 10^(Ma)), method = "lm", se=FALSE,
              color = "black", linewidth=0.5) + 
  scale_x_continuous(limits = c(18.2,29.2), breaks = c(18,20.75, 23.5,26.25,29)) +
  theme(axis.title.y = element_blank(), axis.line.y = element_line(color="black",),
        panel.background = element_rect(fill = "white", color = "white"), 
        axis.line.x.bottom = element_line(color="black"), axis.title.x = element_blank(), 
        axis.text.y = element_blank(),axis.text.x = element_text(size=12))

# Methylophilaceae #

temp_Mp <- ggplot() + geom_point(size=2.3, aes(x=temp,y=10^Mp), col = "orange",alpha = 0.7, shape = 1, stroke = 1) + 
  annotation_custom(grid::textGrob(label = expression(bold("L")),
                                   x = unit(0.96, "npc"), y = unit(0.915, "npc"), 
                                   gp = gpar(color = "black", fontsize=17))) +
  scale_y_continuous(trans = "log10", 
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)),limits=c(0.0023,0.015)) + 
  scale_x_continuous(limits = c(18.2,29.2), breaks = c(18,20.75, 23.5,26.25,29)) +
  theme(axis.title.y = element_blank(), axis.line.y = element_line(color="black"),
        panel.background = element_rect(fill = "white", color = "white"), 
        axis.line.x.bottom = element_line(color="black"), axis.title.x = element_blank(), 
        axis.text = element_blank()) 

#############
## Legends ##
#############

oct_con <- subset(All_16S, Experiment == "OCT_allyl")
methylo_legend <- cowplot::get_legend(ggplot() + theme_classic() +
                                        geom_point(data = oct_con, aes(x = Time, y = Methylacidiphilaceae, color="Methylacidiphilaceae"),size=3.1) + 
                                        geom_point(data = oct_con, aes(x = Time, y = Methylophilaceae, color="Methylophilaceae"),size=3.1) +
                                        scale_color_manual(breaks = c("Methylophilaceae", "Methylacidiphilaceae"), values = c("orange","hotpink1")) +
                                        theme(panel.background = element_rect(fill = "white", color = "white"),
                                              text = element_text(family = "Calibri"),
                                              legend.text = element_text(size = rel(1.0),face = "italic"), legend.title = element_text(size = rel(1.1))) +
                                        guides(color = guide_legend(title = bquote(bold("Nonmethanotrophic 
     Methylotrophs")), order = 1)))

methano_legend <- cowplot::get_legend(ggplot() + theme_classic() +
                                        geom_point(data = oct_con, aes(x = Time, y = Methylococcaceae, color = "Methylococcaceae"),size=3.1) + 
                                        geom_point(data = oct_con, aes(x = Time, y = Methylomonadaceae, color = "Methylomonadaceae"),size=3.1) +
                                        geom_point(data = oct_con, aes(x = Time, y = Methylocystaceae, color = "Methylocystaceae"),size=3.1) +  
                                        scale_color_manual(breaks = c("Methylococcaceae","Methylomonadaceae", "Methylocystaceae"),
                                                           values = c("blue","purple4","green3")) +
                                        guides(color = guide_legend(title = bquote(bold("Aerobic Methanotrophs")), order = 1)) +
                                        theme(panel.background = element_rect(fill = "white", color = "white"),
                                              legend.text = element_text(size = rel(1.0),
                                                                         face = "italic"), 
                                              text = element_text(family = "Calibri"),
                                              legend.title = element_text(size = rel(1.1)))) 


opacity_legend <- cowplot::get_legend(opacity_graph <- ggplot() + geom_point(aes(x = c(2,3,4), y = c(1,2,3),shape = "21"),stroke = 1,
                                                                             color="black",fill = "black",size=2.2) + theme_classic() +
                                        geom_point(aes(x = c(2,3,4), y = c(3,4,5), shape = "1"),color="black",size=2.2,stroke = 1) + 
                                        scale_shape_manual(values = c(21,1),breaks=c("21","1"),labels = c("P < 0.05", "P = N.S.")) +
                                        guides(shape = guide_legend(title = bquote(bold("Significance")))) +
                                        theme(panel.background = element_rect(fill = "white", color = "white"),
                                              legend.text = element_text(size = rel(1.1)), legend.title = element_text(size = rel(1.2)),
                                              text = element_text(family = "Calibri")))


legends.all <- grid.arrange(methano_legend,methylo_legend,opacity_legend,ncol=1,nrow=3)

###############
# Final Plots #
###############

library(extrafont)
font_import(pattern="[A/a]rial", prompt=FALSE)

grid.arrange(o2_Mco, ch4_Mco, temp_Mco,
             o2_Mm, ch4_Mm, temp_Mm,
             o2_Mcy, ch4_Mcy, temp_Mcy,
             o2_Mp, ch4_Mp, temp_Mp,
             o2_Ma, ch4_Ma, temp_Ma,
             nrow = 5, ncol = 3)

corr_w_xaxis <- grid.arrange(o2_Mco, ch4_Mco, temp_Mco,
                             o2_Mm, ch4_Mm, temp_Mm,
                             o2_Mcy, ch4_Mcy, temp_Mcy,
                             o2_Mp, ch4_Mp, temp_Mp,
                             o2_Ma, ch4_Ma, temp_Ma,
                             nrow = 5, ncol = 3, heights = c(1,1,1,1,1.1),widths = c(1.1,1,1),
                             bottom = textGrob(gp=gpar(fontsize = 16.5),
                                               expression("     "~~~~~~"O"[2] ~ "(" *mu*"M)"
                                                          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"CH"[4]~"(nM)"
                                                          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"Temperature "~( degree*C)~"   ")))



correlations_final <- grid.arrange(corr_w_xaxis,legends.all, nrow = 1, ncol = 2, widths=c(3.7,0.79),
                                   left =textGrob(rot=90,label="Relative Abundance",gp=gpar(fontsize=19))) #0.8

#X-axis label locations adjusted in preview