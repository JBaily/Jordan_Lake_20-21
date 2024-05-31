## Figures of Spirulina DNA reads.
library(ggplot2)
library(ggtext)
library(grid)
library(gridExtra)

# Load Data
spir_diff <- read.csv(file = "Spirulina_diff.csv", sep = ",", header = TRUE)
spir_diff$Replicate <- as.character(spir_diff$Replicate)

# Graphs
my.colors <- c("1" = "#00BA38", "2" = "#F8766D", "3" = "#619CFF")

##########
# Oct 21 #
##########

# Open Arrow, lines
allyl <- ggplot() + geom_point(data = spir_diff[spir_diff$Experiment == "A",], 
                      aes(x=Time, y=Proportion, color = Replicate)) + 
  geom_line(data = spir_diff[spir_diff$Experiment == "A",], 
            aes(x=Time, y=Proportion, color = Replicate)) +
  geom_segment(aes(x = 39.16, y = 0.090, xend = 39, yend = 0.079), color="black", 
               arrow = arrow(angle = 13, length = unit(0.25, "cm"))) + #two capsules
  geom_segment(aes(x = 47.6, y = 0.090, xend = 47.6, yend = 0.079), color="black", 
               arrow = arrow(angle = 13, length = unit(0.25, "cm"))) + #ten capsules
  scale_color_manual(values = my.colors) +
  ggplot2::annotate(geom = "richtext", label = "<b>A</b>", x = 120, y = 0.085, size = 10, 
                    fill = NA, label.color = NA) + ylab("Relative Abundance") + xlab("Time (Hours)") + 
  ylim(0,0.09) + theme_classic()



##########
# Nov 21 #
##########

almeoh <- ggplot() + geom_point(data = spir_diff[spir_diff$Experiment == "M",], 
                      aes(x=Time, y=Proportion, color = Replicate)) + 
  geom_line(data = spir_diff[spir_diff$Experiment == "M",], 
             aes(x=Time, y=Proportion, color = Replicate)) +
  geom_segment(aes(x = 0.275, y = 0.1, xend = 0.275, yend = 0.080), color="black",                                                                                       
               arrow = arrow(angle = 13, length = unit(0.3, "cm")))+
  ggplot2::annotate(geom = "richtext", label = "<b>B</b>", x = 30, y = 0.095, size = 10, 
                    fill = NA, label.color = NA) + ylab("Relative Abundance") + xlab("Time (Hours)") + 
  scale_color_manual(values = my.colors) +
  theme(title = element_text(color = "black"),axis.title.x = element_text(color = "black"), 
        axis.title.y = element_text(color = "white"), plot.background = element_rect(color = "white"),
        panel.background = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), legend.background = element_rect(color = "white"),
        axis.line = element_line(colour = "black"), 
        legend.key = element_rect(colour = "transparent", fill = "white"))

##########

grid.arrange(allyl, almeoh, nrow = 1, ncol = 2)
