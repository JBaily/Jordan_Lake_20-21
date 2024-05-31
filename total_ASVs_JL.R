library(dada2); packageVersion("dada2")

#original

fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
plotQualityProfile(fnFs[1:20]) 
plotQualityProfile(fnRs[1:20]) 


#Filtering 
filtFs <- file.path("filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
#cutting after 230 for forward and 200 for reverse (due to summer datasets)
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=c(19,20), 
                     truncLen = c(230,200), maxN=0, maxEE=c(2,2), truncQ=2, 
                     rm.phix=TRUE, compress=TRUE, multithread=TRUE)
head(out)

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers)

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 248:258]
dim(seqtab2)
table(nchar(getSequences(seqtab2)))

seqtab2.nochim <- removeBimeraDenovo(seqtab2, method="consensus", 
                                     multithread=TRUE, verbose=TRUE)
dim(seqtab2.nochim)
sum(seqtab2.nochim)/sum(seqtab2)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), 
               sapply(mergers, getN), rowSums(seqtab2), rowSums(seqtab2.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "cutoff", "nonchim")
rownames(track) <- sample.names
head(track)

section.size <- 1000

sections <- split(c(1:nrow(seqtab2.nochim)),
                  sort(c(1:nrow(seqtab2.nochim))%%ceiling(nrow(seqtab2.nochim)/section.size)))
sections.split <- lapply(sections,
                            function(x){return(assignTaxonomy(seqtab2.nochim[x,], 
                                                              refFasta="taxa_ref/silva_nr99_v138.1_wSpecies_train_set.fa.gz",
                                                              verbose = TRUE))})
taxa <- do.call(rbind, sections.split)

#Removing chloro and mito from seqtab2.nochim

is.chloro <- taxa[,"Order"] %in% "Chloroplast"
is.mito <- taxa[,"Family"] %in% "Mitochondria"

print(dim(seqtab2.nochim))
seqtab2.nochloro.nomito <- seqtab2.nochim[,!(is.chloro | is.mito)]
print(dim(seqtab2.nochloro.nomito))

#Removing chloro and mito from taxa
print(dim(taxa))
taxa.nochloro.nomito <- taxa[!(is.chloro | is.mito),]
print(dim(taxa.nochloro.nomito))

#Removing arthrospira from seqtab2.nochim

is.chloro <- taxa[,"Order"] %in% "Chloroplast"
is.mito <- taxa[,"Family"] %in% "Mitochondria"
is.arthro <- taxa[,"Genus"] %in% "Arthrospira PCC-7345"

print(dim(seqtab2.nochim))
seqtab2.noarthro <- seqtab2.nochim[,!(is.chloro | is.mito | is.arthro)]
print(dim(seqtab2.noarthro))

#Removing arthrospira from taxa
print(dim(taxa))
taxa.noarthro <- taxa[!(is.chloro | is.mito | is.arthro),]
print(dim(taxa.noarthro))


#Save relevant variables for other analyses
save(taxa.noarthro, seqtab2.noarthro, 
     file = "JL_asvs.rda")

#Make data more informative

taxa.print <- taxa.noarthro
rownames(taxa.print) <- NULL
head(taxa.print)

# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab2.noarthro)
asv_headers <- vector(dim(seqtab2.noarthro)[2], mode="character")

for (i in 1:dim(seqtab2.noarthro)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "JL_asvs.fa")

# count table:
asv_tab <- t(seqtab2.noarthro)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "JL_asvs_count.tsv", sep="\t", quote=F, col.names=NA)

#taxonomy table:
row.names(taxa.print) <- sub(">", "", asv_headers)
write.table(taxa.print, "JL_asvs_tax.tsv", sep="\t", quote=F, col.names=NA)  
