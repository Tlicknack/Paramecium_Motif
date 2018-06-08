#DONE
#This program simply takes as an INPUT a modified GFF containing only genes and the fasta file of a genome assembly of interest
#It OUTPUTS a fasta file containing only intergenic sequences


#source("http://bioconductor.org/biocLite.R")
#biocLite("seqinr")
library("seqinr")

makeIntergenicFASTA = function(species_fasta, species_gff){       #uses intergenic gff to make intergenic fasta
  
  cat("", file="/home/tlicknac/Desktop/Paramecium_Genome_Data/Paramecium_FASTA/ptet-intergenic.fasta")   #initiate file, clears it if it already exists. HARDCODED
  #cat("", file="/N/u/tlicknac/Carbonate/Paramecium_GFF/ptet-intergenic.fasta")                            #if running on server
  len_gff = nrow(species_gff)                                   #get number of rows in genic-gff
  
  for(k in 1:len_gff){                                          #interate through that number
    
    scaf = species_gff[[1]][k]                                  #get scaffold name of k'th element in gff
    
    if(species_gff[[2]][k] < species_gff[[3]][k]){              #if gene is on + strand (left to right)
      start_seq = species_gff[[2]][k]
      end_seq = species_gff[[3]][k]
    }
    else{                                                       #else, gene is on the - strand
      start_seq = species_gff[[3]][k]
      end_seq = species_gff[[2]][k]
    }
    
    annotation = paste(">", scaf, "|5'+:", species_gff[[4]][k], "_5'-:", species_gff[[5]][k], "_3'+:", species_gff[[6]][k], "_3'-:", species_gff[[7]][k], sep="")   #create string for >description
    scaf_seq = species_fasta[[scaf]][1]                                                                                                                       #get sequence of entire scaf
    intergenic_seq = substr(scaf_seq, start_seq, end_seq)                                                                                                     #get only intergenic sequence of scaf
    cat(annotation, "\n", intergenic_seq, "\n", file="/home/tlicknac/Desktop/Paramecium_Genome_Data/Paramecium_FASTA/ptet-intergenic.fasta", append=TRUE) #write out fasta file.... HARD CODED OUTFILE    #CHANGE
    #cat(annotation, "\n", intergenic_seq, "\n", file="/N/u/tlicknac/Carbonate/Paramecium_FASTA/ptet-intergenic.fasta", append=TRUE)                 #if running on server
    #intergenic_fasta = paste(annotation, intergenic_seq, sep="\n", append=TRUE)
  }
  #return(intergenic_fasta)
}


species_fasta = read.fasta("/home/tlicknac/Desktop/Paramecium_Genome_Data/Paramecium_FASTA/ptetraurelia_mac_51.fa", as.string=TRUE)         #fasta of 1 species   INPUT              #CHANGE
species_gff = read.table("/home/tlicknac/Desktop/Paramecium_Genome_Data/Paramecium_GFF/ptet-intergenic.tab", sep="\t", as.is=T, header=TRUE)#intergenic gff of 1 species  INPUT      #CHANGE
makeIntergenicFASTA(species_fasta, species_gff)   #Create intergenic fasta file that can be downloaded
