#Microbiological Data
All_16S <- read.csv(file = "16S_all.csv",header = TRUE)
T0_16S <- All_16S[All_16S$Time==0,]

Mco <- log10(T0_16S$Methylococcaceae)
Mm <- log10(T0_16S$Methylomonadaceae)
Mcy <- log10(T0_16S$Methylocystaceae)
Ma <- log10(T0_16S$Methylacidiphilaceae)
Mp <- log10(T0_16S$Methylophilaceae)

plot(Mm,Mco)

plot(Mm,Mp)
cor.test(Mm,Mp)
# p = 0.019, R^2 = 0.33
plot(Mm,Ma)
plot(Mm,log10(10^Ma+10^Mp))

plot(Mco,Mp)
cor.test(Mco,Mp)
# p = 0.02, R^2 = 0.30
plot(Mco,Ma)
cor.test(Mco,Ma)
#p = 0.11
plot(Mco,log10(10^Ma+10^Mp))
