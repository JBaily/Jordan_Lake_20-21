#Load packages 
library(phyloseq)
library(Biostrings)
library(ggplot2)
library(microbiome)
library(RColorBrewer)
library(rioja)
library(microbiomeutilities)

theme_set(theme_classic())

#################################
#      Making data tables       #
#################################

OTUs <- read.table(file = 'JL_asvs_count.tsv', sep = '\t', header = TRUE)
row.names(OTUs) <- OTUs[,1]
OTUs <- OTUs[,-1]
OTUs <- subset(OTUs, select = -c("Set to positions of blanks")) #Remove blanks

TAXA <- read.table(file = 'JL_asvs_tax.tsv', sep = '\t', header = TRUE)
tax_row <- TAXA[,1]
row.names(TAXA) <- TAXA[,1]
TAXA <- TAXA[,-1]

SAMPLE <- read.csv(file = 'JL_asvs_sample_data.csv', sep = ',', header = TRUE)
SAMPLE$Replicate <- as.character(SAMPLE$Replicate)
SAMPLE$Experiment <- as.character(SAMPLE$Experiment)
row.names(SAMPLE) <- SAMPLE[,1]
SAMPLE <- SAMPLE[,-1]

##########################################
#  Convert to phyloseq-friendly objects  #
##########################################

OTU_table <- otu_table(OTUs, taxa_are_rows = TRUE)
TAXA_table <- tax_table(TAXA)
row.names(TAXA_table) <- tax_row
colnames(TAXA_table) <- c("Kingdom", "Phylum", "Class", "Order", "Family", 
                          "Genus", "Species")

SAMPLE_table <- sample_data(SAMPLE)

############################
### Make phyloseq object ###
############################

ps <- phyloseq(OTU_table, SAMPLE_table, TAXA_table)

############################
#   Subset Samples by Exp  #
############################

#Regular 
Oct_20 <- subset_samples(ps, Experiment=="Oct_20")
JL1 <- subset_samples(ps, Experiment=="JL1")
JL2 <- subset_samples(ps, Experiment=="JL2")
JL3 <- subset_samples(ps, Experiment=="JL3")
Allyl <- subset_samples(ps, Experiment=="Allyl")
AlMeOH <- subset_samples(ps, Experiment=="AlMeOH")

#Rarefy
Oct_20.ra <- subset_samples(ps.ra, Experiment=="Oct_20")
JL1.ra <- subset_samples(ps.ra, Experiment=="JL1")
JL2.ra <- subset_samples(ps.ra, Experiment=="JL2")
JL3.ra <- subset_samples(ps.ra, Experiment=="JL3")
Allyl.ra <- subset_samples(ps.ra, Experiment=="Allyl")
AlMeOH.ra <- subset_samples(ps.ra, Experiment=="AlMeOH")

##############################
## Proportion and Rarefying ##
##############################

#Relative abundance - Regular 
Oct_20.prop <- transform_sample_counts(Oct_20, function(x){x / sum(x)})
JL1.prop <- transform_sample_counts(JL1, function(x){x / sum(x)})
JL2.prop <- transform_sample_counts(JL2, function(x){x / sum(x)})
JL3.prop <- transform_sample_counts(JL3, function(x){x / sum(x)})
Allyl.prop <- transform_sample_counts(Allyl, function(x){x / sum(x)})
AlMeOH.prop <- transform_sample_counts(AlMeOH, function(x){x / sum(x)})

#Relative abundance - Rarefying 
Oct_20.ra.prop <- transform_sample_counts(Oct_20.ra, function(x){x / sum(x)})
JL1.ra.prop <- transform_sample_counts(JL1.ra, function(x){x / sum(x)})
JL2.ra.prop <- transform_sample_counts(JL2.ra, function(x){x / sum(x)})
JL3.ra.prop <- transform_sample_counts(JL3.ra, function(x){x / sum(x)})
Allyl.ra.prop <- transform_sample_counts(Allyl.ra, function(x){x / sum(x)})
AlMeOH.ra.prop <- transform_sample_counts(AlMeOH.ra, function(x){x / sum(x)})

###########################################
#     Subsetting Methylo/Methano taxa            #
###########################################

export_group <- function(ps, filename) {
  temp <- as(otu_table(ps), "matrix")
  row.names(temp) <- NULL
  tempdf <- colSums(temp)
  write.csv(tempdf, paste0("Data/Groups_noarthro/",filename))
}

#Oct_20 
Oct_20_cocc <- subset_taxa(Oct_20.prop, Family == "Methylococcaceae")
Oct_20_mona <- subset_taxa(Oct_20.prop, Family == "Methylomonadaceae")
Oct_20_beij <- subset_taxa(Oct_20.prop, Genus == "Methylocystis")
Oct_20_acid <- subset_taxa(Oct_20.prop, Order == "Methylacidiphilales")
Oct_20_phil <- subset_taxa(Oct_20.prop, Family == "Methylophilaceae")

export_group(Oct_20_acid, "Oct_20_acid_noarthro.csv")
export_group(Oct_20_beij, "Oct_20_beij_noarthro.csv")
export_group(Oct_20_cocc, "Oct_20_cocc_noarthro.csv")
export_group(Oct_20_mona, "Oct_20_mona_noarthro.csv")
export_group(Oct_20_phil, "Oct_20_phil_noarthro.csv")

Oct_20_cocc_ra <- subset_taxa(Oct_20.ra.prop, Family == "Methylococcaceae")
Oct_20_mona_ra <- subset_taxa(Oct_20.ra.prop, Family == "Methylomonadaceae")
Oct_20_beij_ra <- subset_taxa(Oct_20.ra.prop, Genus == "Methylocystis")
Oct_20_acid_ra <- subset_taxa(Oct_20.ra.prop, Order == "Methylacidiphilales")
Oct_20_phil_ra <- subset_taxa(Oct_20.ra.prop, Family == "Methylophilaceae")

export_group(Oct_20_acid_ra, "Oct_20_acid_noarthro_ra.csv")
export_group(Oct_20_beij_ra, "Oct_20_beij_noarthro_ra.csv")
export_group(Oct_20_cocc_ra, "Oct_20_cocc_noarthro_ra.csv")
export_group(Oct_20_mona_ra, "Oct_20_mona_noarthro_ra.csv")
export_group(Oct_20_phil_ra, "Oct_20_phil_noarthro_ra.csv")


#JL1 
JL1_cocc <- subset_taxa(JL1.prop, Family == "Methylococcaceae")
JL1_mona <- subset_taxa(JL1.prop, Family == "Methylomonadaceae")
JL1_beij <- subset_taxa(JL1.prop, Genus == "Methylocystis")
JL1_acid <- subset_taxa(JL1.prop, Order == "Methylacidiphilales")
JL1_phil <- subset_taxa(JL1.prop, Family == "Methylophilaceae")

export_group(JL1_acid, "JL1_acid_noarthro.csv")
export_group(JL1_beij, "JL1_beij_noarthro.csv")
export_group(JL1_cocc, "JL1_cocc_noarthro.csv")
export_group(JL1_mona, "JL1_mona_noarthro.csv")
export_group(JL1_phil, "JL1_phil_noarthro.csv")

JL1_cocc_ra <- subset_taxa(JL1.ra.prop, Family == "Methylococcaceae")
JL1_mona_ra <- subset_taxa(JL1.ra.prop, Family == "Methylomonadaceae")
JL1_beij_ra <- subset_taxa(JL1.ra.prop, Genus == "Methylocystis")
JL1_acid_ra <- subset_taxa(JL1.ra.prop, Order == "Methylacidiphilales")
JL1_phil_ra <- subset_taxa(JL1.ra.prop, Family == "Methylophilaceae")

export_group(JL1_acid_ra, "JL1_acid_noarthro_ra.csv")
export_group(JL1_beij_ra, "JL1_beij_noarthro_ra.csv")
export_group(JL1_cocc_ra, "JL1_cocc_noarthro_ra.csv")
export_group(JL1_mona_ra, "JL1_mona_noarthro_ra.csv")
export_group(JL1_phil_ra, "JL1_phil_noarthro_ra.csv")


#JL2 
JL2_cocc <- subset_taxa(JL2.prop, Family == "Methylococcaceae")
JL2_mona <- subset_taxa(JL2.prop, Family == "Methylomonadaceae")
JL2_beij <- subset_taxa(JL2.prop, Genus == "Methylocystis")
JL2_acid <- subset_taxa(JL2.prop, Order == "Methylacidiphilales")
JL2_phil <- subset_taxa(JL2.prop, Family == "Methylophilaceae")

export_group(JL2_acid, "JL2_acid_noarthro.csv")
export_group(JL2_beij, "JL2_beij_noarthro.csv")
export_group(JL2_cocc, "JL2_cocc_noarthro.csv")
export_group(JL2_mona, "JL2_mona_noarthro.csv")
export_group(JL2_phil, "JL2_phil_noarthro.csv")

JL2_cocc_ra <- subset_taxa(JL2.ra.prop, Family == "Methylococcaceae")
JL2_mona_ra <- subset_taxa(JL2.ra.prop, Family == "Methylomonadaceae")
JL2_beij_ra <- subset_taxa(JL2.ra.prop, Genus == "Methylocystis")
JL2_acid_ra <- subset_taxa(JL2.ra.prop, Order == "Methylacidiphilales")
JL2_phil_ra <- subset_taxa(JL2.ra.prop, Family == "Methylophilaceae")

export_group(JL2_acid_ra, "JL2_acid_noarthro_ra.csv")
export_group(JL2_beij_ra, "JL2_beij_noarthro_ra.csv")
export_group(JL2_cocc_ra, "JL2_cocc_noarthro_ra.csv")
export_group(JL2_mona_ra, "JL2_mona_noarthro_ra.csv")
export_group(JL2_phil_ra, "JL2_phil_noarthro_ra.csv")


#JL3
JL3_cocc <- subset_taxa(JL3.prop, Family == "Methylococcaceae")
JL3_mona <- subset_taxa(JL3.prop, Family == "Methylomonadaceae")
JL3_beij <- subset_taxa(JL3.prop, Genus == "Methylocystis")
JL3_acid <- subset_taxa(JL3.prop, Order == "Methylacidiphilales")
JL3_phil <- subset_taxa(JL3.prop, Family == "Methylophilaceae")

export_group(JL3_acid, "JL3_acid_noarthro.csv")
export_group(JL3_beij, "JL3_beij_noarthro.csv")
export_group(JL3_cocc, "JL3_cocc_noarthro.csv")
export_group(JL3_mona, "JL3_mona_noarthro.csv")
export_group(JL3_phil, "JL3_phil_noarthro.csv")

JL3_cocc_ra <- subset_taxa(JL3.ra.prop, Family == "Methylococcaceae")
JL3_mona_ra <- subset_taxa(JL3.ra.prop, Family == "Methylomonadaceae")
JL3_beij_ra <- subset_taxa(JL3.ra.prop, Genus == "Methylocystis")
JL3_acid_ra <- subset_taxa(JL3.ra.prop, Order == "Methylacidiphilales")
JL3_phil_ra <- subset_taxa(JL3.ra.prop, Family == "Methylophilaceae")

export_group(JL3_acid_ra, "JL3_acid_noarthro_ra.csv")
export_group(JL3_beij_ra, "JL3_beij_noarthro_ra.csv")
export_group(JL3_cocc_ra, "JL3_cocc_noarthro_ra.csv")
export_group(JL3_mona_ra, "JL3_mona_noarthro_ra.csv")
export_group(JL3_phil_ra, "JL3_phil_noarthro_ra.csv")


#Allyl 
Allyl_cocc <- subset_taxa(Allyl.prop, Family == "Methylococcaceae")
Allyl_mona <- subset_taxa(Allyl.prop, Family == "Methylomonadaceae")
Allyl_beij <- subset_taxa(Allyl.prop, Genus == "Methylocystis")
Allyl_acid <- subset_taxa(Allyl.prop, Order == "Methylacidiphilales")
Allyl_phil <- subset_taxa(Allyl.prop, Family == "Methylophilaceae")

export_group(Allyl_acid, "Allyl_acid_noarthro.csv")
export_group(Allyl_beij, "Allyl_beij_noarthro.csv")
export_group(Allyl_cocc, "Allyl_cocc_noarthro.csv")
export_group(Allyl_mona, "Allyl_mona_noarthro.csv")
export_group(Allyl_phil, "Allyl_phil_noarthro.csv")

Allyl_cocc_ra <- subset_taxa(Allyl.ra.prop, Family == "Methylococcaceae")
Allyl_mona_ra <- subset_taxa(Allyl.ra.prop, Family == "Methylomonadaceae")
Allyl_beij_ra <- subset_taxa(Allyl.ra.prop, Genus == "Methylocystis")
Allyl_acid_ra <- subset_taxa(Allyl.ra.prop, Order == "Methylacidiphilales")
Allyl_phil_ra <- subset_taxa(Allyl.ra.prop, Family == "Methylophilaceae")

export_group(Allyl_acid_ra, "Allyl_acid_noarthro_ra.csv")
export_group(Allyl_beij_ra, "Allyl_beij_noarthro_ra.csv")
export_group(Allyl_cocc_ra, "Allyl_cocc_noarthro_ra.csv")
export_group(Allyl_mona_ra, "Allyl_mona_noarthro_ra.csv")
export_group(Allyl_phil_ra, "Allyl_phil_noarthro_ra.csv")


#AlMeOH
AlMeOH_cocc_ra <- subset_taxa(AlMeOH.ra.prop, Family == "Methylococcaceae")
AlMeOH_mona_ra <- subset_taxa(AlMeOH.ra.prop, Family == "Methylomonadaceae")
AlMeOH_beij_ra <- subset_taxa(AlMeOH.ra.prop, Genus == "Methylocystis")
AlMeOH_acid_ra <- subset_taxa(AlMeOH.ra.prop, Order == "Methylacidiphilales")
AlMeOH_phil_ra <- subset_taxa(AlMeOH.ra.prop, Family == "Methylophilaceae")

export_group(AlMeOH_acid_ra, "AlMeOH_acid_noarthro_ra.csv")
export_group(AlMeOH_beij_ra, "AlMeOH_beij_noarthro_ra.csv")
export_group(AlMeOH_cocc_ra, "AlMeOH_cocc_noarthro_ra.csv")
export_group(AlMeOH_mona_ra, "AlMeOH_mona_noarthro_ra.csv")
export_group(AlMeOH_phil_ra, "AlMeOH_phil_noarthro_ra.csv")

AlMeOH_cocc <- subset_taxa(AlMeOH.prop, Family == "Methylococcaceae")
AlMeOH_mona <- subset_taxa(AlMeOH.prop, Family == "Methylomonadaceae")
AlMeOH_beij <- subset_taxa(AlMeOH.prop, Genus == "Methylocystis")
AlMeOH_acid <- subset_taxa(AlMeOH.prop, Order == "Methylacidiphilales")
AlMeOH_phil <- subset_taxa(AlMeOH.prop, Family == "Methylophilaceae")

export_group(AlMeOH_acid, "AlMeOH_acid_noarthro.csv")
export_group(AlMeOH_beij, "AlMeOH_beij_noarthro.csv")
export_group(AlMeOH_cocc, "AlMeOH_cocc_noarthro.csv")
export_group(AlMeOH_mona, "AlMeOH_mona_noarthro.csv")
export_group(AlMeOH_phil, "AlMeOH_phil_noarthro.csv")

################################
#    Abundance Bar plots       #
################################
plot_bar(Allyl_methano, fill = "Family", x = "Hours.Elapsed", facet_grid = ~Box, 
         title = "Allyl - Methanotrophs")

plot_bar(Allyl_cocc, fill = "Genus", x = "Hours.Elapsed", facet_grid = ~Box, 
         title = "Allyl - Methylococcaceae")
plot_bar(Allyl_mona, fill = "Genus", x = "Hours.Elapsed", facet_grid = ~Box, 
         title = "Allyl - Methylomonadaceae")
plot_bar(Allyl_beij, fill = "Genus", x = "Hours.Elapsed", facet_grid = ~Box, 
         title = "Allyl - Methylocystis")


plot_bar(Allyl_methylo, fill = "Family", x = "Hours.Elapsed", facet_grid = ~Box, 
         title = "Allyl - Methylotrophs")

plot_bar(Allyl_acid, fill = "Family", x = "Hours.Elapsed", facet_grid = ~Box, 
         title = "Allyl - Methylacidiphilaceae")
plot_bar(Allyl_phil, fill = "Family", x = "Hours.Elapsed", facet_grid = ~Box, 
         title = "Allyl - Methylophilaceae")
plot_bar(Allyl_ros, fill = "Family", x = "Hours.Elapsed", facet_grid = ~Box, 
         title = "Allyl - Methylorosula")
plot_bar(Allyl_vers, fill = "Family", x = "Hours.Elapsed", facet_grid = ~Box, 
         title = "Allyl - Methyloversatilis")


plot_bar(AlMeOH_methano, fill = "Family", x = "Hours.Elapsed", facet_grid = ~Box, 
         title = "AlMeOH - Methanotrophs")
plot_bar(AlMeOH_cocc, fill = "Genus", x = "Hours.Elapsed", facet_grid = ~Box, 
         title = "AlMeOH - Methylococcaceae")
plot_bar(AlMeOH_mona, fill = "Genus", x = "Hours.Elapsed", facet_grid = ~Box, 
         title = "AlMeOH - Methylomonadaceae")
plot_bar(AlMeOH_beij, fill = "Genus", x = "Hours.Elapsed", facet_grid = ~Box, 
         title = "AlMeOH - Methylocystis")

plot_bar(AlMeOH_methylo, fill = "Family", x = "Hours.Elapsed", facet_grid = ~Box, 
         title = "AlMeOH - Methylotrophs")
plot_bar(AlMeOH_acid, fill = "Family", x = "Hours.Elapsed", facet_grid = ~Box, 
         title = "AlMeOH - Methylacidiphilaceae")
plot_bar(AlMeOH_phil, fill = "Family", x = "Hours.Elapsed", facet_grid = ~Box, 
         title = "AlMeOH - Methylophilaceae")
plot_bar(AlMeOH_ros, fill = "Family", x = "Hours.Elapsed", facet_grid = ~Box, 
         title = "AlMeOH - Methylorosula")
plot_bar(AlMeOH_vers, fill = "Family", x = "Hours.Elapsed", facet_grid = ~Box, 
         title = "AlMeOH - Methyloversatilis")


plot_bar(all_methano, fill = "Family", x = "Hours.Elapsed", facet_grid = ~Box, 
         title = "all - Methanotrophs")
plot_bar(all_methylo, fill = "Family", x = "Hours.Elapsed", facet_grid = ~Box, 
         title = "all - Methylotrophs")
