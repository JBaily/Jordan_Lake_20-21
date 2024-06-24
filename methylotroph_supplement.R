library(ggplot2)
library(ggnewscale)
library(grid)
library(gridExtra)

### Load Data ##

all_data <- read.csv(file = "16S_all.csv")
all_data$Box <- as.character(all_data$Box)

JL1 <- subset(all_data, Experiment == "July_2021_25")
JL2 <- subset(all_data, Experiment == "July_2021_50")
JL3 <- subset(all_data, Experiment == "July_2021_BES")

############
## Graphs ##
############

#########
## JL1 ##
#########

JL1_acid <- ggplot() + geom_point(data = JL1, aes(x = Time, y = Methylacidiphilaceae, 
                                                  color="Methylacidiphilaceae")) +
  geom_line(data = JL1, aes(x = Time, y = Methylacidiphilaceae, linetype = Box, 
                            color="Methylacidiphilaceae")) + 
  labs(x="Time (Hours)", y="Relative Abundance",title = expression(bold("June 2021 - 1"))) + 
  geom_vline(xintercept = 67.7) + 
  annotation_custom(grid::textGrob(label = expression(bold("A")),
                                   x = unit(0.94, "npc"), y = unit(0.95, "npc"), 
                                   gp = gpar(color = "black", fontsize=17))) +
  scale_color_manual(breaks = c("Methylacidiphilaceae"), values = c("hotpink1")) +
  scale_linetype_manual(breaks = c("1","2","3"), values=c(3,1,2)) +
  guides(linetype = "none", color = "none") + scale_x_continuous(limits=c(0,76)) +
  theme(axis.title.x = element_text(color = "white",size = rel(1.2)), 
        axis.line.y = element_line(color="black",),
        panel.background = element_rect(fill = "white", color = "white"), 
        axis.line.x.bottom = element_line(color="black"), 
        axis.title.y = element_text(color = "black"), 
        legend.key = element_rect(fill = "white", color = "white"),
        title = element_text(color = "black")) 

JL1_phil <- ggplot() + 
  geom_point(data = JL1, aes(x = Time, y = Methylophilaceae, color="Methylophilaceae")) +
  geom_line(data = JL1, aes(x = Time, y = Methylophilaceae, linetype=Box, color="Methylophilaceae")) +
  labs(x="Time (Hours)", y="Relative Abundance",title = expression(bold("June 2021 - 1"))) + 
  geom_vline(xintercept = 67.7) + 
  annotation_custom(grid::textGrob(label = expression(bold("B")),
                                   x = unit(0.94, "npc"), y = unit(0.95, "npc"), 
                                   gp = gpar(color = "black", fontsize=17))) +
  scale_color_manual(breaks = c("Methylophilaceae"), values = c("orange")) +
  scale_linetype_manual(breaks = c("1","2","3"), values=c(3,1,2)) +
  guides(linetype = "none", color = "none") + scale_x_continuous(limits=c(0,76)) +
  theme(axis.title.x = element_text(color = "white",size = rel(1.2)), 
        axis.line.y = element_line(color="black",),
        panel.background = element_rect(fill = "white", color = "white"), 
        axis.line.x.bottom = element_line(color="black"), axis.title.y = element_blank(), 
        legend.key = element_rect(fill = "white", color = "white"),
        title = element_text(color = "white"))

#########
## JL2 ##
#########

scale_Y_axis <- function(x) sprintf("%.3f", x)

JL2_acid <- ggplot() + geom_point(data = JL2, aes(x = Time, y = Methylacidiphilaceae, color="Methylacidiphilaceae")) +
  geom_line(data = JL2, aes(x = Time, y = Methylacidiphilaceae, linetype = Box, color="Methylacidiphilaceae")) + 
  labs(x="Time (Hours)", y="Relative Abundance", col = "Family", linetype = "Incubation Device",
       title = expression(bold("June 2021 - 2"))) +
  scale_color_manual(breaks = c("Methylacidiphilaceae"), values = c("hotpink1")) +
  scale_linetype_manual(breaks = c("1","2","3"), values=c(3,1,2)) +
  annotation_custom(grid::textGrob(label = expression(bold("C")),
                                   x = unit(0.94, "npc"), y = unit(0.95, "npc"), 
                                   gp = gpar(color = "black", fontsize=17))) + 
  scale_y_continuous(labels=scale_Y_axis) + geom_vline(xintercept = 51.7) + 
  guides(linetype = "none", color = "none") + scale_x_continuous(limits=c(0,60)) + #was 68
  theme(axis.title.x = element_text(color = "white",size = rel(1.2)), 
        axis.line.y = element_line(color="black",),
        panel.background = element_rect(fill = "white", color = "white"), 
        axis.line.x.bottom = element_line(color="black"), 
        axis.title.y = element_text(color = "black"), 
        legend.key = element_rect(fill = "white", color = "white"),
        title = element_text(color = "black")) 

JL2_phil <- ggplot() + geom_point(data = JL2, aes(x = Time, y = Methylophilaceae, color="Methylophilaceae")) +
  geom_line(data = JL3, aes(x = Time, y = Methylophilaceae, linetype = Box, color="Methylophilaceae")) + 
  labs(x="Time (Hours)", y="Relative Abundance", col = "Family", linetype = "Incubation Device",
       title = expression(bold("June 2021 - 2"))) +
  scale_color_manual(breaks = c("Methylophilaceae"), values = c("orange")) +
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

#########
## JL3 ##
#########

JL3_acid <- ggplot() + geom_point(data = JL3, aes(x = Time, y = Methylacidiphilaceae, color="Methylacidiphilaceae")) +
  geom_line(data = JL3, aes(x = Time, y = Methylacidiphilaceae, linetype = Box, color="Methylacidiphilaceae")) + 
  labs(x="Time (Hours)", y="Relative Abundance", col = "Family", linetype = "Incubation Device",
       title = expression(bold("July 2021"))) +
  scale_color_manual(breaks = c("Methylacidiphilaceae"), values = c("hotpink1")) +
  scale_linetype_manual(breaks = c("1","2","3"), values=c(3,1,2)) + geom_vline(xintercept = 25.6) + 
  annotation_custom(grid::textGrob(label = expression(bold("E")),
                                   x = unit(0.94, "npc"), y = unit(0.95, "npc"), 
                                   gp = gpar(color = "black", fontsize=17))) + 
  guides(linetype = "none", color = "none") + scale_x_continuous(limits=c(0,50)) +
  theme(axis.title.x = element_text(color = "black",size = rel(1.2)), 
        axis.line.y = element_line(color="black",),
        panel.background = element_rect(fill = "white", color = "white"), 
        axis.line.x.bottom = element_line(color="black"), 
        axis.title.y = element_text(color = "black"), 
        legend.key = element_rect(fill = "white", color = "white"),
        title = element_text(color = "black")) 

JL3_phil <- ggplot() + geom_point(data = JL3, aes(x = Time, y = Methylophilaceae, color="Methylophilaceae")) +
  geom_line(data = JL3, aes(x = Time, y = Methylophilaceae, linetype = Box, color="Methylophilaceae")) + 
  labs(x="Time (Hours)", y="Relative Abundance", col = "Family", linetype = "Incubation Device",
       title = expression(bold("July 2021"))) +
  scale_color_manual(breaks = c("Methylophilaceae"), values = c("orange")) +
  scale_linetype_manual(breaks = c("1","2","3"), values=c(3,1,2)) +
  annotation_custom(grid::textGrob(label = expression(bold("F")),
                                   x = unit(0.94, "npc"), y = unit(0.95, "npc"), 
                                   gp = gpar(color = "black", fontsize=17))) + 
  scale_y_continuous(labels=scale_Y_axis) + geom_vline(xintercept = 25.6) + 
  guides(linetype = "none", color = "none") + scale_x_continuous(limits=c(0,50)) + #
  theme(axis.title.x = element_text(color = "black",size = rel(1.2)), 
        axis.line.y = element_line(color="black",),
        panel.background = element_rect(fill = "white", color = "white"), 
        axis.line.x.bottom = element_line(color="black"), axis.title.y = element_blank(), 
        legend.key = element_rect(fill = "white", color = "white"),
        title = element_text(color = "white"))
#################

###########
# Legends #
###########

library(cowplot)


boxes_legend <- cowplot::get_legend(ggplot() + 
                                      geom_line(data = JL1, aes(x = Time, y = Methylacidiphilaceae, linetype = Box)) + 
                                      geom_line(data = JL1, aes(x = Time, y = Methylophilaceae, linetype=Box)) +
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
                                        geom_line(data = JL1, aes(x = Time, y = Methylacidiphilaceae, color="Methylacidiphilaceae")) + 
                                        geom_line(data = JL1, aes(x = Time, y = Methylophilaceae, color="Methylophilaceae")) +
                                        scale_color_manual(breaks = c("Methylophilaceae","Methylacidiphilaceae"), values = c("orange","hotpink1"))+ 
                                        guides(color = guide_legend(title = bquote(bold("Nonmethanotrophic 
     Methylotrophs")), order = 1)) + theme_classic() +
                                        theme(panel.background = element_rect(fill = "white", color = "white"),
                                              legend.text = element_text(size = rel(1.0),face = "italic"), legend.title = element_text(size = rel(1.2))))

grid.newpage()
grid.draw(methylo_legend)
grid.newpage()


legends.all <- grid.arrange(methylo_legend,boxes_legend,ncol=1,nrow=2)


##############
## Combined ##
##############

methylo.graphs <- grid.arrange(JL1_acid, JL1_phil, 
                               JL2_acid, JL2_phil,
                               JL3_acid, JL3_phil,
                            nrow=3, ncol=2, widths=c(5,4.74))

methylo.all <- grid.arrange(methylo.graphs,legends.all, nrow = 1, ncol = 2, widths=c(3.5,1))
methylo.all


##############