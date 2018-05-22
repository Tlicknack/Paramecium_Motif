#This program will compare multiple sequence allignments (fasta) and calculate "distances' between them
#These distances will determine if this block of sequences in similar across the ~200nts
#This will be used to determine which fastas will be run through MEME for motif discovery

#source("http://bioconductor.org/biocLite.R")
#biocLite("seqinr")
library("seqinr")


outputDist = function(input_list, seq_dist){
  
}



calcDist = function(input_fasta, len_fasta){
  
  sum_of_mismatches = 0
  
  for(k in 1:len_fasta){
    counter = 1
    for(l in 1:len_fasta){
      
      k_letter = substr(input_fasta[[k]], start=counter, stop=counter)
      l_letter = substr(input_fasta[[l]], start=counter, stop=counter)
      
      if(k_letter != l_letter){
        
        sum_of_mismatches = sum_of_mismatches + 1
        print(sum_of_mismatches)
      }
      
      if(k_letter == "-" & l_letter ==  "-"){
        sum_of_mismatches = sum_of_mismatches + 1
      }
      counter = counter + 1
    }
  }
  return(sum_of_mismatches)
}


#MAIN
#input_list = list.files("/N/u/tlicknac/Carbonate/Upstream_Alligned_Test/", "*.fasta", all.files=FALSE, recursive=FALSE)  #get all fasta files in a directory

input_list = list.files("/home/tlicknac/Desktop/Paramecium_Genome_Data/Upstream_Alligned_Test", "*.fasta", all.files=FALSE, recursive=FALSE)  #get all fasta files in a directory

for(i in 1:length(input_list)) {                                                      #iterate through fasta file names
  
  seq_dist = 0
  
  #path = "/N/u/tlicknac/Carbonate/Upstream_Alligned_Test/"                            #save path to files as a string
  path = "/home/tlicknac/Desktop/Paramecium_Genome_Data/Upstream_Alligned_Test/"
  path_i = paste(path, input_list[[i]], sep="")                                       #concatenate path and i'th filename
  input_fasta = read.alignment(path_i, format="fasta", forceToLower = TRUE)           #input the fasta files using allignment reader, which saves elements in the following way:
  #[[1]] = number of seqs in the current fasta file,    [[2]] = list of geneIDs,                          [[3]] = list of sequence strings
  #[[1]] has 1 dimension,                               [[2]][i] allows you to go acress i'th geneID,     [[3]][i] allows you to access the i'th seq
  
  len_fasta = length(input_fasta[[3]])
  #print(input_fasta[[3]])
  
  seq_dist = calcDist(input_fasta[[3]], len_fasta)
  #ouputDist(input_list[[i]], seq_dist)
  print("Next file:")
}



