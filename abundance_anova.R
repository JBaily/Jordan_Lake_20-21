library(ggplot2)
library(ggtext)
library(scales)
library(ggsignif)
library(multcompView)
library(tidyverse)
library(dplyr)


All_16S <- read.csv(file = "16S_all.csv",header = TRUE)
T0_16S <- All_16S[All_16S$Time==0,]

Mco <- log10(T0_16S$Methylococcaceae)
Mm <- log10(T0_16S$Methylomonadaceae)
Mcy <- log10(T0_16S$Methylocystaceae)
Ma <- log10(T0_16S$Methylacidiphilaceae)
Mp <- log10(T0_16S$Methylophilaceae)

rel_abund <- c(Ma,Mp,Mcy,Mco,Mm)

families <- c(rep(c("Methylacidiphilaceae"),16,each=1),
              rep(c("Methylophilaceae"),16,each=1),
              rep(c("Methylocystaceae"),16,each=1),
              rep(c("Methylococcaceae"),16,each=1),
              rep(c("Methylomonadaceae"),16,each=1))

abundances <- data.frame(families,rel_abund)
colnames(abundances) <- c("groups","abund")
abundances$groups <- factor(abundances$groups , levels=c("Methylacidiphilaceae",
                                                         "Methylophilaceae", 
                                                         "Methylocystaceae",
                                                         "Methylococcaceae",
                                                         "Methylomonadaceae"))


res.aov <- aov(abund ~ groups, data = abundances)
summary(res.aov)
TukeyHSD(res.aov)
#                                             diff         lwr         upr     p adj
#Methylococcaceae-Methylacidiphilaceae  -1.3034356 -1.72495462 -0.88191656 0.0000000
#Methylocystaceae-Methylacidiphilaceae  -1.1750697 -1.59658870 -0.75355064 0.0000000
#Methylomonadaceae-Methylacidiphilaceae -1.6503418 -2.07186084 -1.22882278 0.0000000
#Methylophilaceae-Methylacidiphilaceae  -0.7202440 -1.14176305 -0.29872499 0.0000831
#Methylocystaceae-Methylococcaceae       0.1283659 -0.29315311  0.54988495 0.9134649
#Methylomonadaceae-Methylococcaceae     -0.3469062 -0.76842525  0.07461281 0.1561955
#Methylophilaceae-Methylococcaceae       0.5831916  0.16167254  1.00471060 0.0021099
#Methylomonadaceae-Methylocystaceae     -0.4752721 -0.89679117 -0.05375311 0.0191986
#Methylophilaceae-Methylocystaceae       0.4548257  0.03330662  0.87634468 0.0279748
#Methylophilaceae-Methylomonadaceae      0.9300978  0.50857876  1.35161682 0.0000003


tukey_results <- TukeyHSD(res.aov)
letters_tukey <- data.frame(multcompLetters2(formula = abund ~ groups, data = abundances,
                                            x = tukey_results$groups[,4])$Letters)
colnames(letters_tukey)[1] <- "Letter"
letters_tukey$Groups <- rownames(letters_tukey)
rownames(letters_tukey) <- NULL

placement <- abundances %>% 
  group_by(groups) %>%
  summarise(quantile(abund)[4])

colnames(placement)[2] <- "Placement.Value"
colnames(placement)[1] <- "Groups"
letters_tukey <- left_join(letters_tukey, placement)

letters_tukey[,3] <- 10^(letters_tukey[,3])


########
# Plot #
########

set.seed(11)
abund_plot <- ggplot(abundances, aes(x=groups, y=10^(abund)), color = black) + #put color= groups in aes for edges
  labs(y = "Relative Abundance", x = "Families") +
  geom_boxplot(outlier.shape = NA, aes(fill=groups), alpha = 0.5) +
  scale_y_continuous(trans = "log10", breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)),
                     limits = 10^c(-4.1,-1)) +
  geom_jitter(size=1, shape=21, alpha=1, aes(fill=groups),position = position_jitter(0.25)) + #alpha is 0.5 for most versions
  #geom_jitter(size=1, pch=21, alpha=1, aes(color=groups),position = position_jitter(0.25)) +
  geom_text(data = letters_tukey, aes(x = Groups, y = Placement.Value, label = Letter), 
            size = 4, color = "black", vjust = c(-3.2,-4,-6,-9.6,-6.7)) +
  scale_color_manual(values = c("hotpink1","orange","green","blue","purple4"),
                     breaks = c("Methylacidiphilaceae","Methylophilaceae",
                                "Methylocystaceae","Methylococcaceae",
                                "Methylomonadaceae")) + guides(color = "none", fill = "none") +
  scale_fill_manual(values = c("hotpink1","orange","green","blue","purple4"),
                     breaks = c("Methylacidiphilaceae","Methylophilaceae",
                                "Methylocystaceae","Methylococcaceae",
                                "Methylomonadaceae")) + 
  theme_classic() + theme(axis.text.x.bottom = element_text(face = "italic"))

abund_plot


############
## Means  ##
############

#Since data are log10-transformed, the means will be back-transformed to correct units. 

10^(mean(Mco)) #0.001389472
10^(mean(Mm)) #0.0006250916
10^(mean(Mcy)) #0.001867307
10^(mean(Ma)) #0.0279438
10^(mean(Mp)) #0.005321591


