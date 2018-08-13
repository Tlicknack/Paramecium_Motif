#DONE
#This program simply takes as an INPUT a modified GFF containing only genes and the fasta file of a genome assembly of interest
#It OUTPUTS a fasta file containing only intergenic sequences


#source("http://bioconductor.org/biocLite.R")
#biocLite("seqinr")
library("seqinr")

species_fasta = read.fasta("/N/u/tlicknac/Carbonate/Paramecium_FASTA/biaurelia_V1-4_assembly_v1.fa", as.string=TRUE)
species_intergenic_gff = read.table("/N/u/tlicknac/Carbonate/Paramecium_GFF/pbi-intergenic.tab", skip=2)   #gff's ended up with first row of NAs
len_intergenic_gff = nrow(species_intergenic_gff)

for(k in 1:len_intergenic_gff){                                          #iterate through that number
  scaf = species_intergenic_gff[[1]][k]                                  #get scaffold name of k'th element in gff
  if(species_intergenic_gff[[2]][k] < species_intergenic_gff[[3]][k]){              #if gene is on + strand (left to right)
    start_seq = species_intergenic_gff[[2]][k]
    end_seq = species_intergenic_gff[[3]][k]
  }
  else{                                                       #else, gene is on the - strand
    start_seq = species_intergenic_gff[[3]][k]
    end_seq = species_intergenic_gff[[2]][k]
  }
  annotation = paste(">", scaf, "|5'+:", species_intergenic_gff[[4]][k], "_5'-:", species_intergenic_gff[[5]][k], "_3'+:", species_intergenic_gff[[6]][k], "_3'-:", species_intergenic_gff[[7]][k], sep="")
  scaf_seq = species_fasta[[scaf]][1]                                                                                                                       #get sequence of entire scaf
  intergenic_seq = as.character(substr(scaf_seq, start_seq, end_seq))                                                                                                     #get only intergenic sequence of scaf
  cat(annotation, "\n", intergenic_seq, "\n", file="/N/u/tlicknac/Carbonate/Paramecium_FASTA/pbi-intergenic.fasta", append=TRUE)
  cat("Row completed: ", k)
}
