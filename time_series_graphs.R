### Data ###
library(ggplot2)
library(ggtext)
library(grid)
library(gridExtra)
library(colorblindcheck)
library(gplots)
library(scales)
library(cowplot)

#Colorblindness palette check
palette_check(col2hex(c("black","darkgray","orange","hotpink1",
                        "purple4","blue","green")), plot = TRUE)

#16S
all_data <- read.csv(file = "16S_all.csv")
all_data$Box <- as.character(all_data$Box)

oct_2020 <- subset(all_data, Experiment == "OCT_2020")
oct_2020 <- oct_2020[-1,] #remove B1T0
oct_2020 <- oct_2020[-1,] #remove B1T1
june_1_21 <- subset(all_data, Experiment == "JL1")
june_2_21 <- subset(all_data, Experiment == "JL2")
july_2021 <- subset(all_data, Experiment == "JL3")
oct_2021 <- subset(all_data, Experiment == "OCT_21")
nov_2021 <- subset(all_data, Experiment == "NOV_21")

###################
## October 2020  ##
###################

oct_20_methylo <- ggplot() + geom_point(data = oct_2020, aes(x = Time, y = Methylacidiphilaceae, color="Methylacidiphilaceae")) +
  geom_line(data = oct_2020, aes(x = Time, y = Methylacidiphilaceae, linetype = Box, color="Methylacidiphilaceae")) + 
  geom_point(data = oct_2020, aes(x = Time, y = Methylophilaceae, color="Methylophilaceae")) +
  geom_line(data = oct_2020, aes(x = Time, y = Methylophilaceae, linetype=Box, color="Methylophilaceae")) +
  labs(x="Time (Hours)", col = "Family", linetype = "Incubation Device",
       title = expression(bold("October 2020"))) +
  scale_linetype_manual(breaks = c("1","2","3"), values=c(3,1,2)) + geom_vline(xintercept = 13.6) +
  scale_color_manual(breaks = c("Methylophilaceae", "Methylacidiphilaceae"), values = c("orange","hotpink1")) +
  annotation_custom(grid::textGrob(label = expression(bold("B")),
                                   x = unit(0.94, "npc"), y = unit(0.95, "npc"), 
                                   gp = gpar(color = "black", fontsize=17))) +
  theme(axis.title.y = element_blank(), axis.line.y = element_line(color="black",),
        panel.background = element_rect(fill = "white", color = "white"), 
        axis.line.x.bottom = element_line(color="black"), 
        axis.title.x = element_text(color = "white",size = rel(1.2)), 
        legend.key = element_rect(fill = "white", color = "white"),
        title = element_text(color = "white")) + 
  guides(linetype = "none", color = "none") 


oct_20_methano <- ggplot() + geom_point(data = oct_2020, aes(x = Time, y = Methylococcaceae, color = "Methylococcaceae")) +
  geom_line(data = oct_2020, aes(x = Time, y = Methylococcaceae, color = "Methylococcaceae",linetype=Box)) + 
  geom_point(data = oct_2020, aes(x = Time, y = Methylomonadaceae, color = "Methylomonadaceae"))+
  geom_line(data = oct_2020, aes(x = Time, y = Methylomonadaceae, color = "Methylomonadaceae",linetype=Box)) +
  geom_point(data = oct_2020, aes(x = Time, y = Methylocystaceae, color = "Methylocystaceae"))+
  geom_line(data = oct_2020, aes(x = Time, y = Methylocystaceae, color = "Methylocystaceae",linetype=Box)) +  
  labs(x="Time (Hours)", y="Relative Abundance", col = "Family", linetype = "Incubation Device",
       title = expression(bold("October 2020")))+
  scale_color_manual(breaks = c("Methylomonadaceae", "Methylococcaceae", "Methylocystaceae"), values = c("purple4","blue","green3")) +
  guides(linetype = "none", color = "none") + geom_vline(xintercept = 13.6) +
  annotation_custom(grid::textGrob(label = expression(bold("A")),
                                   x = unit(0.94, "npc"), y = unit(0.95, "npc"), 
                                   gp = gpar(color = "black", fontsize=17))) +
  scale_x_continuous(limits=c(0,35)) + scale_linetype_manual(breaks = c("1","2","3"), values=c(3,1,2)) +
  theme(axis.title.x = element_text(color = "white",size = rel(1.2)), 
        axis.line.y = element_line(color="black",),
        panel.background = element_rect(fill = "white", color = "white"), 
        axis.line.x.bottom = element_line(color="black"), 
        legend.key = element_rect(fill = "white", color = "white"),
        axis.title.y = element_text(size = rel(1.2)))

#################
## Allyl 2021  ##
#################

oct_21_methylo <- ggplot() + geom_point(data = oct_2021, aes(x = Time, y = Methylacidiphilaceae, color="Methylacidiphilaceae")) +
  geom_line(data = oct_2021, aes(x = Time, y = Methylacidiphilaceae, linetype = Box, color="Methylacidiphilaceae")) + 
  geom_point(data = oct_2021, aes(x = Time, y = Methylophilaceae, color="Methylophilaceae")) +
  geom_line(data = oct_2021, aes(x = Time, y = Methylophilaceae, linetype=Box, color="Methylophilaceae")) +
  labs(x="Time (Hours)", col = "Family", linetype = "Incubation Device",
       title = expression(bold("October 2021"))) +
  guides(linetype = "none", color = "none") + 
  annotation_custom(grid::textGrob(label = expression(bold("D")),
                                   x = unit(0.94, "npc"), y = unit(0.95, "npc"), 
                                   gp = gpar(color = "black", fontsize=17))) +
  scale_color_manual(breaks = c("Methylophilaceae", "Methylacidiphilaceae"), values = c("orange","hotpink1"))  +
  scale_linetype_manual(breaks = c("1","2","3"), values=c(3,1,2)) +
  geom_segment(aes(x = 39.16, y = 0.0415, xend = 39.16, yend = 0.038), color="black",                                                                                       
               arrow = arrow(angle = 18, length = unit(0.15, "cm"))) +
  geom_segment(aes(x = 47.61, y = 0.032, xend = 47.61, yend = 0.0285), color="black",                                                                                       
               arrow = arrow(angle = 18, length = unit(0.15, "cm"))) +
  scale_x_continuous(limits=c(0,62)) + geom_vline(xintercept = 54.7) +
  theme(axis.title.x = element_text(color = "white",size = rel(1.2)), 
        axis.line.y = element_line(color="black"),
        panel.background = element_rect(fill = "white", color = "white"), 
        axis.line.x.bottom = element_line(color="black"), axis.title.y = element_blank(), 
        legend.key = element_rect(fill = "white", color = "white"),
        title = element_text(color = "white")) 


oct_21_methano <- ggplot() + geom_point(data = oct_2021, aes(x = Time, y = Methylococcaceae, color = "Methylococcaceae")) +
  geom_line(data = oct_2021, aes(x = Time, y = Methylococcaceae, color = "Methylococcaceae",linetype=Box)) + 
  geom_point(data = oct_2021, aes(x = Time, y = Methylomonadaceae, color = "Methylomonadaceae"))+
  geom_line(data = oct_2021, aes(x = Time, y = Methylomonadaceae, color = "Methylomonadaceae",linetype=Box)) +
  geom_point(data = oct_2021, aes(x = Time, y = Methylocystaceae, color = "Methylocystaceae"))+
  geom_line(data = oct_2021, aes(x = Time, y = Methylocystaceae, color = "Methylocystaceae",linetype=Box)) +  
  labs(x="Time (Hours)", y="Relative Abundance", title = expression(bold("October 2021")),
       col = "Family", linetype = "Incubation Device")+ geom_vline(xintercept = 54.7) +
  scale_color_manual(breaks = c("Methylomonadaceae", "Methylococcaceae", "Methylocystaceae"), 
                     values = c("purple4","blue","green3")) +
  guides(linetype = "none", color = "none") +
  geom_segment(aes(x = 39.16, y = 0.0044, xend = 39.16, yend = 0.0036), color="black",                                                                                       
               arrow = arrow(angle = 18, length = unit(0.15, "cm"))) +
  geom_segment(aes(x = 47.61, y = 0.0035, xend = 47.61, yend = 0.0027), color="black",                                                                                       
               arrow = arrow(angle = 18, length = unit(0.15, "cm"))) +
  annotation_custom(grid::textGrob(label = expression(bold("C")),
                                   x = unit(0.94, "npc"), y = unit(0.95, "npc"), 
                                   gp = gpar(color = "black", fontsize=17))) +
  scale_linetype_manual(breaks = c("1","2","3"), values=c(3,1,2)) +
  scale_x_continuous(limits=c(0,62)) + 
  theme(axis.title.x = element_text(color = "white",size = rel(1.2)), 
        axis.line.y = element_line(color="black",),
        panel.background = element_rect(fill = "white", color = "white"), 
        axis.line.x.bottom = element_line(color="black"), 
        legend.key = element_rect(fill = "white", color = "white"),
        axis.title.y = element_text(size = rel(1.2))) 

##################
## AlMeOH 2021  ##
##################

nov_21_methylo <- ggplot() + geom_point(data = nov_2021, aes(x = Time, y = Methylacidiphilaceae, color="Methylacidiphilaceae")) +
  geom_line(data = nov_2021, aes(x = Time, y = Methylacidiphilaceae, linetype = Box, color="Methylacidiphilaceae")) + 
  geom_point(data = nov_2021, aes(x = Time, y = Methylophilaceae, color="Methylophilaceae")) + 
  guides(linetype = "none", color = "none") + geom_vline(xintercept = 11.9) + 
  geom_line(data = nov_2021, aes(x = Time, y = Methylophilaceae, linetype=Box, color="Methylophilaceae")) +
  labs(x="Time (Hours)", col = "Family", linetype = "Incubation Device",
       title = expression(bold("November 2021"))) +
  scale_color_manual(breaks = c("Methylophilaceae", "Methylacidiphilaceae"), values = c("orange","hotpink1"))  +
  annotation_custom(grid::textGrob(label = expression(bold("F")),
                                      x = unit(0.94, "npc"), y = unit(0.95, "npc"), 
                                      gp = gpar(color = "black", fontsize=17))) + 
  scale_x_continuous(limits=c(0,18.5)) +
  geom_segment(aes(x = 0.225, y = 0.018, xend = 0.275, yend = 0.0165), color="black",                                                                                       
               arrow = arrow(angle = 18, length = unit(0.15, "cm"))) + 
  scale_linetype_manual(breaks = c("1","2","3"), values=c(3,1,2)) +
  theme(axis.line.y = element_line(color="black"),
        panel.background = element_rect(fill = "white", color = "white"), 
        axis.line.x.bottom = element_line(color="black"), axis.title.y = element_blank(),
        legend.key = element_rect(fill = "white", color = "white"),
        axis.title.x = element_text(size = rel(1.2), color = "black"),
        title = element_text(color = "white")) 

nov_21_methano <- ggplot() + geom_point(data = nov_2021, aes(x = Time, y = Methylococcaceae, color = "Methylococcaceae")) +
  geom_line(data = nov_2021, aes(x = Time, y = Methylococcaceae, linetype=Box, color = "Methylococcaceae")) + 
  geom_point(data = nov_2021, aes(x = Time, y = Methylomonadaceae, color = "Methylomonadaceae"))+
  geom_line(data = nov_2021, aes(x = Time, y = Methylomonadaceae, linetype=Box, color = "Methylomonadaceae")) +
  geom_point(data = nov_2021, aes(x = Time, y = Methylocystaceae, color = "Methylocystaceae"))+
  geom_line(data = nov_2021, aes(x = Time, y = Methylocystaceae, linetype=Box, color = "Methylocystaceae")) +  
  labs(x="Time (Hours)", y="Relative Abundance", col = "Family", linetype = "Incubation Device",
       title = expression(bold("November 2021")))+
  scale_color_manual(breaks = c("Methylomonadaceae", "Methylococcaceae", "Methylocystaceae"), 
                     values = c("purple4","blue","green3"))  + 
  guides(linetype = "none", color = "none") + geom_vline(xintercept = 11.9) + 
  annotation_custom(grid::textGrob(label = expression(bold("E")),
                                   x = unit(0.94, "npc"), y = unit(0.95, "npc"), 
                                   gp = gpar(color = "black", fontsize=17))) + 
  scale_x_continuous(limits=c(0,18.5)) +
  geom_segment(aes(x = 0.225, y = 0.0132, xend = 0.275, yend = 0.012), color="black",                                                                                       
               arrow = arrow(angle = 18, length = unit(0.15, "cm")))+
  scale_linetype_manual(breaks = c("1","2","3"), values=c(3,1,2)) +
  theme(axis.line.y = element_line(color="black",),
        panel.background = element_rect(fill = "white", color = "white"), 
        axis.line.x.bottom = element_line(color="black"),
        legend.key = element_rect(fill = "white", color = "white"),
        axis.title.x = element_text(size = rel(1.2)),
        axis.title.y = element_text(size = rel(1.2))) 


##################
## June 1 2021  ##
##################

june_1_methylo <- ggplot() + geom_point(data = june_1_21, aes(x = Time, y = Methylacidiphilaceae, color="Methylacidiphilaceae")) +
  geom_line(data = june_1_21, aes(x = Time, y = Methylacidiphilaceae, linetype = Box, color="Methylacidiphilaceae")) + 
  geom_point(data = june_1_21, aes(x = Time, y = Methylophilaceae, color="Methylophilaceae")) +
  geom_line(data = june_1_21, aes(x = Time, y = Methylophilaceae, linetype=Box, color="Methylophilaceae")) +
  labs(x="Time (Hours)", y="Relative Abundance",title = expression(bold("June 2021 - 1"))) + 
  geom_vline(xintercept = 67.7) + 
  annotation_custom(grid::textGrob(label = expression(bold("B")),
                                   x = unit(0.94, "npc"), y = unit(0.95, "npc"), 
                                   gp = gpar(color = "black", fontsize=17))) +
  scale_color_manual(breaks = c("Methylophilaceae", "Methylacidiphilaceae"), values = c("orange","hotpink1")) +
  scale_linetype_manual(breaks = c("1","2","3"), values=c(3,1,2)) +
  guides(linetype = "none", color = "none") + scale_x_continuous(limits=c(0,76)) +
  theme(axis.title.x = element_text(color = "white",size = rel(1.2)), 
        axis.line.y = element_line(color="black",),
        panel.background = element_rect(fill = "white", color = "white"), 
        axis.line.x.bottom = element_line(color="black"), axis.title.y = element_blank(), 
        legend.key = element_rect(fill = "white", color = "white"),
        title = element_text(color = "white")) 


june_1_methano <- ggplot() + geom_point(data = june_1_21, aes(x = Time, y = Methylococcaceae, color = "Methylococcaceae")) +
  geom_line(data = june_1_21, aes(x = Time, y = Methylococcaceae, color = "Methylococcaceae",linetype=Box)) + 
  geom_point(data = june_1_21, aes(x = Time, y = Methylomonadaceae, color = "Methylomonadaceae"))+
  geom_line(data = june_1_21, aes(x = Time, y = Methylomonadaceae, color = "Methylomonadaceae",linetype=Box)) +
  geom_point(data = june_1_21, aes(x = Time, y = Methylocystaceae, color = "Methylocystaceae"))+
  geom_line(data = june_1_21, aes(x = Time, y = Methylocystaceae, color = "Methylocystaceae",linetype=Box)) +  
  labs(x="Time (Hours)", y="Relative Abundance", col = "Family", linetype = "Incubation Device",
       title = expression(bold("June 2021 - 1")))+
  scale_linetype_manual(breaks = c("1","2","3"), values=c(3,1,2)) + geom_vline(xintercept = 67.7) + 
  scale_color_manual(breaks = c("Methylomonadaceae", "Methylococcaceae", "Methylocystaceae"), values = c("purple4","blue","green3")) +
  guides(linetype = "none", color = "none") + scale_x_continuous(limits=c(0,75)) +
  annotation_custom(grid::textGrob(label = expression(bold("A")),
                                   x = unit(0.94, "npc"), y = unit(0.95, "npc"), 
                                   gp = gpar(color = "black", fontsize=17))) +
  theme(axis.title.x = element_text(color = "white",size = rel(1.2)), 
        axis.line.y = element_line(color="black",),
        panel.background = element_rect(fill = "white", color = "white"), 
        axis.line.x.bottom = element_line(color="black"), 
        legend.key = element_rect(fill = "white", color = "white"),
        axis.title.y = element_text(size = rel(1.2))) 


##################
## June 2 2021  ##
##################

scale_Y_axis <- function(x) sprintf("%.3f", x)

june_2_methylo <- ggplot() + geom_point(data = june_2_21, aes(x = Time, y = Methylacidiphilaceae, color="Methylacidiphilaceae")) +
  geom_line(data = june_2_21, aes(x = Time, y = Methylacidiphilaceae, linetype = Box, color="Methylacidiphilaceae")) + 
  geom_point(data = june_2_21, aes(x = Time, y = Methylophilaceae, color="Methylophilaceae")) +
  geom_line(data = june_2_21, aes(x = Time, y = Methylophilaceae, linetype=Box, color="Methylophilaceae")) +
  labs(x="Time (Hours)", y="Relative Abundance", col = "Family", linetype = "Incubation Device",
       title = expression(bold("June 2021 - 2"))) +
  scale_color_manual(breaks = c("Methylophilaceae", "Methylacidiphilaceae"), values = c("orange","hotpink1")) +
  scale_linetype_manual(breaks = c("1","2","3"), values=c(3,1,2)) +
  annotation_custom(grid::textGrob(label = expression(bold("D")),
                                   x = unit(0.94, "npc"), y = unit(0.95, "npc"), 
                                   gp = gpar(color = "black", fontsize=17))) + 
  scale_y_continuous(labels=scale_Y_axis) + geom_vline(xintercept = 51.7) + 
  guides(linetype = "none", color = "none") + scale_x_continuous(limits=c(0,60)) + #was 68
  theme(axis.title.x = element_text(color = "white",size = rel(1.2)), 
        axis.line.y = element_line(color="black",),
        panel.background = element_rect(fill = "white", color = "white"), 
        axis.line.x.bottom = element_line(color="black"), axis.title.y = element_blank(), 
        legend.key = element_rect(fill = "white", color = "white"),
        title = element_text(color = "white")) 


june_2_methano <- ggplot() + geom_point(data = june_2_21, aes(x = Time, y = Methylococcaceae, color = "Methylococcaceae")) +
  geom_line(data = june_2_21, aes(x = Time, y = Methylococcaceae, linetype=Box, color = "Methylococcaceae")) + 
  geom_point(data = june_2_21, aes(x = Time, y = Methylomonadaceae, color = "Methylomonadaceae"))+
  geom_line(data = june_2_21, aes(x = Time, y = Methylomonadaceae, linetype=Box, color = "Methylomonadaceae")) +
  geom_point(data = june_2_21, aes(x = Time, y = Methylocystaceae, color = "Methylocystaceae"))+
  geom_line(data = june_2_21, aes(x = Time, y = Methylocystaceae, linetype=Box, color = "Methylocystaceae")) +  
  labs(x="Time (Hours)", y="Relative Abundance", col = "Family", linetype = "Incubation Device",
       title = expression(bold("June 2021 - 2")))+
  scale_linetype_manual(breaks = c("1","2","3"), values=c(3,1,2)) + geom_vline(xintercept = 51.7) + 
  scale_color_manual(breaks = c("Methylomonadaceae", "Methylococcaceae", "Methylocystaceae"), values = c("purple4","blue","green3")) +
  guides(linetype = "none", color = "none") + scale_x_continuous(limits=c(0,60)) + 
  annotation_custom(grid::textGrob(label = expression(bold("C")),
                                   x = unit(0.94, "npc"), y = unit(0.95, "npc"), 
                                   gp = gpar(color = "black", fontsize=17))) + 
  scale_y_continuous(labels = label_comma(accuracy = 0.0001)) +
  theme(axis.title.x = element_text(color = "white",size = rel(1.2)), 
        axis.line.y = element_line(color="black",),
        panel.background = element_rect(fill = "white", color = "white"), 
        axis.line.x.bottom = element_line(color="black"), 
        legend.key = element_rect(fill = "white", color = "white"),
        axis.title.y = element_text(size = rel(1.2)),
        axis.text.y.left = element_text(size=7.2)) 

###############
## July 2021  ##
###############

july_21_methylo <- ggplot() + geom_point(data = july_2021, aes(x = Time, y = Methylacidiphilaceae, color="Methylacidiphilaceae")) +
  geom_line(data = july_2021, aes(x = Time, y = Methylacidiphilaceae, linetype = Box, color="Methylacidiphilaceae")) + 
  geom_point(data = july_2021, aes(x = Time, y = Methylophilaceae, color="Methylophilaceae")) +
  geom_line(data = july_2021, aes(x = Time, y = Methylophilaceae, linetype=Box, color="Methylophilaceae")) +
  labs(x="Time (Hours)", y="Relative Abundance", col = "Family", linetype = "Incubation Device",
       title = expression(bold("July 2021"))) +
  scale_color_manual(breaks = c("Methylophilaceae", "Methylacidiphilaceae"), values = c("orange","hotpink1")) +
  scale_linetype_manual(breaks = c("1","2","3"), values=c(3,1,2)) + geom_vline(xintercept = 25.6) + 
  annotation_custom(grid::textGrob(label = expression(bold("F")),
                                   x = unit(0.94, "npc"), y = unit(0.95, "npc"), 
                                   gp = gpar(color = "black", fontsize=17))) +
  guides(linetype = "none", color = "none") + scale_x_continuous(limits=c(0,50)) +
  theme(axis.line.y = element_line(color="black",),
        panel.background = element_rect(fill = "white", color = "white"), 
        axis.line.x.bottom = element_line(color="black"), axis.title.y = element_blank(),
        legend.key = element_rect(fill = "white", color = "white"),
        axis.title.x = element_text(size = rel(1.2), color = "black"),
        title = element_text(color = "white")) 

july_21_methano <- ggplot() + geom_point(data = july_2021, aes(x = Time, y = Methylococcaceae, color = "Methylococcaceae")) +
  geom_line(data = july_2021, aes(x = Time, y = Methylococcaceae, color = "Methylococcaceae",linetype=Box)) + 
  geom_point(data = july_2021, aes(x = Time, y = Methylomonadaceae, color = "Methylomonadaceae"))+
  geom_line(data = july_2021, aes(x = Time, y = Methylomonadaceae, color = "Methylomonadaceae",linetype=Box)) +
  geom_point(data = july_2021, aes(x = Time, y = Methylocystaceae, color = "Methylocystaceae"))+
  geom_line(data = july_2021, aes(x = Time, y = Methylocystaceae, color = "Methylocystaceae",linetype=Box)) +  
  labs(x="Time (Hours)", y="Relative Abundance", col = "Family", linetype = "Incubation Device",
       title = expression(bold("July 2021")))+
  scale_color_manual(breaks = c("Methylomonadaceae", "Methylococcaceae", "Methylocystaceae"), values = c("purple4","blue","green3")) +
  scale_linetype_manual(breaks = c("1","2","3"), values=c(3,1,2)) +
  annotation_custom(grid::textGrob(label = expression(bold("E")),
                                   x = unit(0.94, "npc"), y = unit(0.95, "npc"), 
                                   gp = gpar(color = "black", fontsize=17))) +
  guides(linetype = "none", color = "none") + geom_vline(xintercept = 25.6) +
  scale_x_continuous(limits=c(0,50)) + theme_classic() + 
  theme(axis.title.x = element_text(size = rel(1.2)), axis.title.y = element_text(size = rel(1.2)))


###########
# Legends #
###########


boxes_legend <- cowplot::get_legend(ggplot() + 
  geom_line(data = oct_2021, aes(x = Time, y = Methylacidiphilaceae, linetype = Box)) + 
  geom_line(data = oct_2021, aes(x = Time, y = Methylophilaceae, linetype=Box)) +
  labs(x="Time (Hours)", y="Relative Abundance", col = "Family", linetype = "Incubation Device") +
  scale_linetype_manual(breaks = c("1","2","3"), values=c(3,1,2)) +
  guides(linetype = guide_legend(title = bquote(bold("iBag Replicate No.")), order = 1))+
  theme_classic() + theme(panel.background = element_rect(fill = "white", color = "white"),
                          legend.text = element_text(size = rel(1.0)), 
                          legend.title = element_text(size = rel(1.2))))

grid.newpage()
grid.draw(boxes_legend)
grid.newpage()

methylo_legend <- cowplot::get_legend(ggplot() + 
  geom_line(data = oct_2021, aes(x = Time, y = Methylacidiphilaceae, color="Methylacidiphilaceae")) + 
  geom_line(data = oct_2021, aes(x = Time, y = Methylophilaceae, color="Methylophilaceae")) +
  scale_color_manual(breaks = c("Methylophilaceae", "Methylacidiphilaceae"), values = c("orange","hotpink1"))+ 
  guides(color = guide_legend(title = bquote(bold("Nonmethanotrophic 
     Methylotrophs")), order = 1)) + theme_classic() +
    theme(panel.background = element_rect(fill = "white", color = "white"),
          legend.text = element_text(size = rel(1.0),face = "italic"), legend.title = element_text(size = rel(1.2))))

grid.newpage()
grid.draw(methylo_legend)
grid.newpage()


methano_legend <- cowplot::get_legend(ggplot() + theme_classic() +
  geom_line(data = oct_2021, aes(x = Time, y = Methylococcaceae, color = "Methylococcaceae")) + 
  geom_line(data = oct_2021, aes(x = Time, y = Methylomonadaceae, color = "Methylomonadaceae")) +
  geom_line(data = oct_2021, aes(x = Time, y = Methylocystaceae, color = "Methylocystaceae")) +  
  scale_color_manual(breaks = c("Methylomonadaceae","Methylococcaceae",  "Methylocystaceae"), 
                     values = c("purple4","blue","green3")) +
  guides(color = guide_legend(title = bquote(bold("Aerobic Methanotrophs")), order = 1))+ theme_classic() +
    theme(panel.background = element_rect(fill = "white", color = "white"), #legend.text.align = 0.5,
          legend.text = element_text(size = rel(1.0),face = "italic"), legend.title = element_text(size = rel(1.2))))

grid.newpage()
grid.draw(methano_legend)
grid.newpage()

grid.arrange(methano_legend,methylo_legend,boxes_legend,ncol=1,nrow=3)

###########
# Arrange #
###########

legends.all <- grid.arrange(methano_legend,methylo_legend,boxes_legend,ncol=1,nrow=3)

#Fall 
fall.graphs <- grid.arrange(oct_20_methano, oct_20_methylo, 
                            oct_21_methano, oct_21_methylo,
                            nov_21_methano, nov_21_methylo,
                            nrow=3, ncol=2, widths=c(5,4.74))

fall.all <- grid.arrange(fall.graphs,legends.all, nrow = 1, ncol = 2, widths=c(3.5,1))

#Summer 
summer.graphs <- grid.arrange(june_1_methano, june_1_methylo,
                            june_2_methano, june_2_methylo,
                            july_21_methano, july_21_methylo,
                            nrow=3, ncol=2, widths=c(5,4.74))


summer.all <- grid.arrange(summer.graphs,legends.all, nrow = 1, ncol = 2, widths=c(3.5,1))
