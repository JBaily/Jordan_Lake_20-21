library(phyloseq)

setwd("~/Lloyd Lab Dropbox/Jennifer Baily/AMOR/Total_ASVs/Analysis")

###################
## Rate Constant ##
###################

############
## Family ##
############

T0_ps <- subset_samples(ps.prop.fam, Hours_Elapsed == "0")

temp_taxtable <- data.frame(as(tax_table(T0_ps), "matrix"))
temp_counttable <- t(data.frame(as(otu_table(T0_ps), "matrix")))

all_relabund <- temp_counttable

all_relabund <- all_relabund[,colSums(all_relabund != 0) == nrow(all_relabund)]

rate_constant <- c(0.327,0.462,0.448,
                   0.087,0.075,0.051,
                   0.029,0.042,0.038,
                   0.238,0.193,
                   0.236,0.22,0.231,
                   0.246,0.242) #Reordered values to match phyloseq sample reordering. 

correlations <- t(as.data.frame(cor(rate_constant, log10(all_relabund))))
colnames(correlations) <- c("Correlation_coef")
correlations <- as.data.frame(correlations)

corr_sort <- correlations[order(-correlations$Correlation_coef),,drop=FALSE]

#Sorting and filling in Family names 

fam_corr_sort <- data.frame(matrix(nrow = nrow(corr_sort), ncol = 8))
row.names(fam_corr_sort) <- row.names(corr_sort)
colnames(fam_corr_sort) <- c("Phylum","Family","Pos_Neg","Correlation_coef","t","df","P-value",
                                 "Significant")

i <- 1

for(i in 1:nrow(fam_corr_sort)){
  
  ASV <- row.names(fam_corr_sort)[i]
  
  fam_corr_sort[i,4] <- corr_sort[i,1] #Correlation_coef
  fam_corr_sort[i,3] <- ifelse(corr_sort[i,1] == 0, "Zero", 
                                   ifelse(corr_sort[i,1] > 0, "Positive","Negative"))
  
  temp_cor_test <- cor.test(rate_constant,log10(all_relabund[,ASV]))
  
  fam_corr_sort[i,5] <- temp_cor_test$statistic
  fam_corr_sort[i,6] <- temp_cor_test$parameter
  fam_corr_sort[i,7] <- temp_cor_test$p.value
  fam_corr_sort[i,8] <- ifelse(temp_cor_test$p.value <= 0.05, "Yes","No")
  
  
  fam_corr_sort[i,1] <- temp_taxtable[ASV,2]
  fam_corr_sort[i,2] <- temp_taxtable[ASV,5]
  
  if(is.na(fam_corr_sort[i,2]) == TRUE){
    
    fam_corr_sort[i,2] <- "NA"
    
  }
  
  if(is.na(fam_corr_sort[i,1]) == TRUE){
    
    fam_corr_sort[i,1] <- "NA"
    
  }

}


###########
## Genus ##
###########

T0_ps_gen <- subset_samples(ps.prop, Hours_Elapsed == "0")

temp_taxtable_gen <- data.frame(as(tax_table(T0_ps_gen), "matrix"))
temp_counttable_gen <- t(data.frame(as(otu_table(T0_ps_gen), "matrix")))

all_relabund_gen <- temp_counttable_gen

all_relabund_gen <- all_relabund_gen[,colSums(all_relabund_gen != 0) == nrow(all_relabund_gen)]

rate_constant_gen <- c(0.327,0.462,0.448,
                   0.087,0.075,0.051,
                   0.029,0.042,0.038,
                   0.238,0.193,
                   0.236,0.22,0.231,
                   0.246,0.242) #Reordered values to match phyloseq sample reordering. 

correlations_gen <- t(as.data.frame(cor(rate_constant_gen, log10(all_relabund_gen))))
colnames(correlations_gen) <- c("Correlation_coef")
correlations_gen <- as.data.frame(correlations_gen)

corr_sort_gen <- correlations_gen[order(-correlations_gen$Correlation_coef),,drop=FALSE]

#Sorting and filling in Family names 

i <- 1

fam_corr_sort_gen <- data.frame(matrix(nrow = nrow(corr_sort_gen), ncol = 9))
row.names(fam_corr_sort_gen) <- row.names(corr_sort_gen)
colnames(fam_corr_sort_gen) <- c("Phylum","Family","Genus","Pos_Neg","Correlation_coef","t","df","P-value",
                                 "Significant")

for(i in 1:nrow(fam_corr_sort_gen)){
  
  ASV_gen <- row.names(fam_corr_sort_gen)[i]
  
  fam_corr_sort_gen[i,5] <- corr_sort_gen[i,1] #Correlation_coef
  fam_corr_sort_gen[i,4] <- ifelse(corr_sort_gen[i,1] == 0, "Zero", 
                                   ifelse(corr_sort_gen[i,1] > 0, "Positive","Negative"))
  
  temp_cor_test_gen <- cor.test(rate_constant_gen,log10(all_relabund_gen[,ASV_gen]))
  
  fam_corr_sort_gen[i,6] <- temp_cor_test_gen$statistic
  fam_corr_sort_gen[i,7] <- temp_cor_test_gen$parameter
  fam_corr_sort_gen[i,8] <- temp_cor_test_gen$p.value
  fam_corr_sort_gen[i,9] <- ifelse(temp_cor_test_gen$p.value <= 0.05, "Yes","No")
  
  
  fam_corr_sort_gen[i,1] <- temp_taxtable_gen[ASV_gen,2]
  fam_corr_sort_gen[i,2] <- temp_taxtable_gen[ASV_gen,5]
  fam_corr_sort_gen[i,3] <- temp_taxtable_gen[ASV_gen,6]
  
  if(is.na(fam_corr_sort_gen[i,3]) == TRUE){
    
    fam_corr_sort_gen[i,3] <- "NA"
    
  }
  
  if(is.na(fam_corr_sort_gen[i,2]) == TRUE){
    
    fam_corr_sort_gen[i,2] <- "NA"
    
  }
  
  if(is.na(fam_corr_sort_gen[i,1]) == TRUE){
    
    fam_corr_sort_gen[i,1] <- "NA"
    
  }
  
  
  
}



###########

write.csv(fam_corr_sort,"JL_Family_correlations.csv")
write.csv(fam_corr_sort_gen,"JL_Genera_correlations.csv")

##################
## Initial Rate ##
##################

############
## Family ##
############

T0_ps_Rate <- subset_samples(ps.prop.fam, Hours_Elapsed == "0")

temp_taxtable_Rate <- data.frame(as(tax_table(T0_ps_Rate), "matrix"))
temp_counttable_Rate <- t(data.frame(as(otu_table(T0_ps_Rate), "matrix")))

all_relabund_Rate <- temp_counttable_Rate

all_relabund_Rate <- all_relabund_Rate[,colSums(all_relabund_Rate != 0) == nrow(all_relabund_Rate)]


rate_constant <- c(0.327,0.462,0.448,
                   0.087,0.075,0.051,
                   0.029,0.042,0.038,
                   0.238,0.193,
                   0.236,0.22,0.231,
                   0.246,0.242) #Reordered values to match phyloseq sample reordering. 

inital_methane <- c(334.07, 334.60, 369.88,
                    301.60, 331.44, 318.03,
                    362.09, 378.98, 434.36,
                    254.92, 195.16,
                    327.93,333.83,347.38,
                    379.15,408.75) #Reordered values to match phyloseq sample reordering. 

rate_AMOR <- rate_constant * inital_methane


correlations_Rate <- t(as.data.frame(cor(rate_AMOR, log10(all_relabund_Rate))))
colnames(correlations_Rate) <- c("Correlation_coef")
correlations_Rate <- as.data.frame(correlations_Rate)

corr_sort_Rate <- correlations_Rate[order(-correlations_Rate$Correlation_coef),,drop=FALSE]

#Sorting and filling in Family names 

fam_corr_sort_Rate <- data.frame(matrix(nrow = nrow(corr_sort_Rate), ncol = 8))
row.names(fam_corr_sort_Rate) <- row.names(corr_sort_Rate)
colnames(fam_corr_sort_Rate) <- c("Phylum","Family","Pos_Neg","Correlation_coef","t","df","P-value",
                             "Significant")

i <- 1

for(i in 1:nrow(fam_corr_sort_Rate)){
  
  ASV <- row.names(fam_corr_sort_Rate)[i]
  
  fam_corr_sort_Rate[i,4] <- corr_sort_Rate[i,1] #Correlation_coef
  fam_corr_sort_Rate[i,3] <- ifelse(corr_sort_Rate[i,1] == 0, "Zero", 
                               ifelse(corr_sort_Rate[i,1] > 0, "Positive","Negative"))
  
  temp_cor_test_Rate <- cor.test(rate_AMOR,log10(all_relabund_Rate[,ASV]))
  

  fam_corr_sort_Rate[i,5] <- temp_cor_test_Rate$statistic
  fam_corr_sort_Rate[i,6] <- temp_cor_test_Rate$parameter
  fam_corr_sort_Rate[i,7] <- temp_cor_test_Rate$p.value
  fam_corr_sort_Rate[i,8] <- ifelse(temp_cor_test_Rate$p.value <= 0.05, "Yes","No")
  
  
  fam_corr_sort_Rate[i,1] <- temp_taxtable_Rate[ASV,2]
  fam_corr_sort_Rate[i,2] <- temp_taxtable_Rate[ASV,5]
  
  if(is.na(fam_corr_sort_Rate[i,2]) == TRUE){
    
    fam_corr_sort_Rate[i,2] <- "NA"
    
  }
  
  if(is.na(fam_corr_sort_Rate[i,1]) == TRUE){
    
    fam_corr_sort_Rate[i,1] <- "NA"
    
  }
  
}


###########
## Genus ##
###########

T0_ps_gen_Rate <- subset_samples(ps.prop, Hours_Elapsed == "0")

temp_taxtable_gen_Rate <- data.frame(as(tax_table(T0_ps_gen_Rate), "matrix"))
temp_counttable_gen_Rate <- t(data.frame(as(otu_table(T0_ps_gen_Rate), "matrix")))

all_relabund_gen_Rate <- temp_counttable_gen_Rate

all_relabund_gen_Rate <- all_relabund_gen_Rate[,colSums(all_relabund_gen_Rate != 0) == nrow(all_relabund_gen_Rate)]

rate_constant_gen <- c(0.327,0.462,0.448,
                   0.087,0.075,0.051,
                   0.029,0.042,0.038,
                   0.238,0.193,
                   0.236,0.22,0.231,
                   0.246,0.242) #Reordered values to match phyloseq sample reordering. 

initial_methane_gen <- c(334.07, 334.60, 369.88,
     301.60, 331.44, 318.03,
     362.09, 378.98, 434.36,
     254.92, 195.16,
     327.93,333.83,347.38,
     379.15,408.75) #Reordered values to match phyloseq sample reordering. 

rate_AMOR_gen <- rate_constant_gen * initial_methane_gen 

correlations_gen_Rate <- t(as.data.frame(cor(rate_AMOR_gen, log10(all_relabund_gen_Rate))))
colnames(correlations_gen_Rate) <- c("Correlation_coef")
correlations_gen_Rate <- as.data.frame(correlations_gen_Rate)

corr_sort_gen_Rate <- correlations_gen_Rate[order(-correlations_gen_Rate$Correlation_coef),,drop=FALSE]

#Sorting and filling in Family names 

i <- 1

fam_corr_sort_gen_Rate <- data.frame(matrix(nrow = nrow(corr_sort_gen_Rate), ncol = 9))
row.names(fam_corr_sort_gen_Rate) <- row.names(corr_sort_gen_Rate)
colnames(fam_corr_sort_gen_Rate) <- c("Phylum","Family","Genus","Pos_Neg","Correlation_coef","t","df","P-value",
                                 "Significant")

for(i in 1:nrow(fam_corr_sort_gen_Rate)){
  
  ASV_gen_Rate <- row.names(fam_corr_sort_gen_Rate)[i]
  
  fam_corr_sort_gen_Rate[i,5] <- corr_sort_gen_Rate[i,1] #Correlation_coef
  fam_corr_sort_gen_Rate[i,4] <- ifelse(corr_sort_gen_Rate[i,1] == 0, "Zero", 
                                   ifelse(corr_sort_gen_Rate[i,1] > 0, "Positive","Negative"))
  
  temp_cor_test_gen_Rate <- cor.test(rate_AMOR_gen,log10(all_relabund_gen_Rate[,ASV_gen_Rate]))
  
  fam_corr_sort_gen_Rate[i,6] <- temp_cor_test_gen_Rate$statistic
  fam_corr_sort_gen_Rate[i,7] <- temp_cor_test_gen_Rate$parameter
  fam_corr_sort_gen_Rate[i,8] <- temp_cor_test_gen_Rate$p.value
  fam_corr_sort_gen_Rate[i,9] <- ifelse(temp_cor_test_gen_Rate$p.value <= 0.05, "Yes","No")
  
  
  fam_corr_sort_gen_Rate[i,1] <- temp_taxtable_gen_Rate[ASV_gen_Rate,2]
  fam_corr_sort_gen_Rate[i,2] <- temp_taxtable_gen_Rate[ASV_gen_Rate,5]
  fam_corr_sort_gen_Rate[i,3] <- temp_taxtable_gen_Rate[ASV_gen_Rate,6]
  
  if(is.na(fam_corr_sort_gen_Rate[i,3]) == TRUE){
    
    fam_corr_sort_gen_Rate[i,3] <- "NA"
    
  }
  
  if(is.na(fam_corr_sort_gen_Rate[i,2]) == TRUE){
    
    fam_corr_sort_gen_Rate[i,2] <- "NA"
    
  }
  
  if(is.na(fam_corr_sort_gen_Rate[i,1]) == TRUE){
    
    fam_corr_sort_gen_Rate[i,1] <- "NA"
    
  }
  
}



###########

write.csv(fam_corr_sort_Rate,"JL_Family_correlations_Rate.csv")
write.csv(fam_corr_sort_gen_Rate,"JL_Genera_correlations_Rate.csv")


#############################################
# qPCR pmoA vs qPCR bacterial 16S vs Oxygen #
#############################################

qPCR_metadata <- read.csv("JL_extra_metadata.csv")
qPCR_metadata

cor.test(qPCR_metadata$pmoA_Mco_avg,qPCR_metadata$Rate.Constant.hr..1)
#t = 2.6343, df = 10, p-value = 0.02498, R = 0.6400487, R^2 0.4096623

cor.test(qPCR_metadata$pmoA_Mco_avg,qPCR_metadata$Modeled.Initial.Oxygen..µM.)
#t = -2.6339, df = 10, p-value = 0.02499, R = -0.6399948, R^2 0.4095933

cor.test(qPCR_metadata$pmoA_Mco_avg,qPCR_metadata$Modeled.Initial.Methane..nM.)
#t = -0.91316, df = 10, p-value = 0.3826, R = -0.2774311, R^2 0.07696802


cor.test(qPCR_metadata$X16S_bacterial_avg,qPCR_metadata$Rate.Constant.hr..1)
#t = -4.9719, df = 12, p-value = 0.0003242, R = -0.8204908, R^2 0.6732052

cor.test(qPCR_metadata$X16S_bacterial_avg,qPCR_metadata$Modeled.Initial.Oxygen..µM.)
#t = 3.8664, df = 12, p-value = 0.002243, R = 0.7447893, R^2 0.5547111

cor.test(qPCR_metadata$X16S_bacterial_avg,qPCR_metadata$Modeled.Initial.Methane..nM.)
#t = 0.2693, df = 12, p-value = 0.7923, R = 0.07750727, R^2 0.006007377


cor.test(log10(qPCR_metadata$X16S_Mco_bacterial_rel),log10(qPCR_metadata$pmoA_rel_bac))

