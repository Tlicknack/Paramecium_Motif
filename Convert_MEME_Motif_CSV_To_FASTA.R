#DONE
#This converts a table of motifs into a fasta file of motifs


#source("http://bioconductor.org/biocLite.R")
#biocLite(seqinr)
#library(seqinr)

motifs= read.csv("MEMExmlParaOrtholog.csv", header=T, sep=",")
flist = list()
for(i in 1:nrow(motifs)){
  motif = as.character(motifs[i,2])
  flist[[i]] = motif
}
len = length(flist)

write.fasta(sequences=flist, names=1:len, file.out="Motifs_ParaOrtholog.fasta")
