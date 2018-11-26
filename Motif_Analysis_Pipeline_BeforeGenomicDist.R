#For motif File without annotation

setwd("~/Desktop/Paramecium_Genome_Data")

dfmotif = read.csv("12nt_MEME_All_Info.csv", header=T)

#Compare to yeast database
dfyeast_tfbs = read.table("Yeast_TFBS.tab", header=F)
colnames(dfyeast_tfbs) = c("Element", "Transcription-Factor", "Regex", "Something")

vhits = c()
for(i in 1:nrow(dfyeast_tfbs)){
  hit = gregexpr(dfyeast_tfbs$Regex[i], dfmotif$Regex, ignore.case = T)
  hits = which(as.numeric(hit) != -1)
  vhits = append(vhits, hits)
}

dfhits = dfyeast_tfbs[vhits,]
dfhits = dfhits[which(is.na(dfhits$Element) == F),]

#Single variant motifs
single_variant = dfmotif[which(dfmotif$MotifMismatch == 0),]
nrow(single_variant)
median(single_variant$EValue)
hist(log(as.numeric(format(single_variant$EValue, scientific=T)), base=10), breaks=200, main="Distribution of E-values for Single Variant Motifs", xlab="log(E-value)")

#Motifs found multiple times
dup_motifs = dfmotif[which(duplicated(dfmotif$Consensus)),]
median(as.numeric(format(dup_motifs$EValue, scientific=T)))
hist(log(as.numeric(format(dup_motifs$EValue, scientific=T)), base=10), breaks=100, main="Distribution of E-values for Duplicated Motifs", xlab="log(E-value)")

#Motifs close to start codon
utr_motifs = dfmotif[which(dfmotif$DistanceToStart < 15),]
nrow(utr_motifs)
median(utr_motifs$EValue)
hist(log(as.numeric(format(utr_motifs$EValue, scientific=T)), base=10), breaks=100, main="Distribution of E-values for UTR Motifs", xlab="log(E-value)")
dfmotif[which(dfmotif$DistanceToStart == 12),]
twleve = (dfmotif[which(dfmotif$DistanceToStart == 12),])


#High E-value
high_evale = dfmotif[as.numeric(which(dfmotif$EValue)) < 1e-140,]
