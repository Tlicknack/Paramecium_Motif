#This program will compare multiple sequence allignments (fasta) and calculate "distances' between them
#These distances will determine if this block of sequences in similar across the ~200nts
#This will be used to determine which fastas will be run through MEME for motif discovery

#source("http://bioconductor.org/biocLite.R")
#biocLite("seqinr")
library("seqinr")


calcDist = function(input_fasta, len_fasta, input_names){
  
  sum_of_mismatches = 0                                                     #Initiate variable for counting mismatches between sequences
  
  for(k in 1:len_fasta){                                                    #Iterate k through length of fasta
    counter = 1                                                             #Initiate counter to go through all 201nt in each upstream region
    for(l in 1:len_fasta){                                                  #Do the same for l .... k=1-l=1, k=1-l=2, k=1-l=3 ... will go through all combinations
      
      k_letter = substr(input_fasta[[k]], start=counter, stop=counter)      #Get 1st, 2nd, 3rd, ..., etc. letter from the k'th sequence in the fasta
      l_letter = substr(input_fasta[[l]], start=counter, stop=counter)      #Do the same for the l'th sequence in the fasta
      
      if(k_letter != l_letter){                                             #If the whatever'th letter in each sequence is different, then add one to the mismatch counter
        
        sum_of_mismatches = sum_of_mismatches + 1 
        #print(sum_of_mismatches)
      }
    }
  }
  name_dist = basename(input_names)
  name_dist = paste("/N/u/tlicknac/Carbonate/Paramecium_Aligned_Distances/", name_dist, sep='')      #Change path to determine where it OUTPUTS to                CHANGE
  cat(sum_of_mismatches, file=name_dist)                                              #Print this out to a file of the same name as the input file
  print("File number:", k)
}




#MAIN

input_list = list.files("/N/u/tlicknac/Carbonate/Paramecium_Upstream_Alligned/", "*.fasta", all.files=FALSE, recursive=FALSE)  #get all fasta files in a directory          CHANGE
#input_list = list.files("/home/tlicknac/Desktop/Paramecium_Genome_Data/Upstream_Alligned_Test", "*.fasta", all.files=FALSE, recursive=FALSE)  #get all fasta files in a directory

for(i in 1:length(input_list)) {                                                      #iterate through fasta file names
  
  seq_dist = 0
  
  path = "/N/u/tlicknac/Carbonate/Paramecium_Upstream_Alligned/"                            #save path to files as a string                                                 CHANGE
  #path = "/home/tlicknac/Desktop/Paramecium_Genome_Data/Upstream_Alligned_Test/"
  path_i = paste(path, input_list[[i]], sep="")                                       #concatenate path and i'th filename
  input_fasta = read.alignment(path_i, format="fasta", forceToLower = TRUE)           #input the fasta files using allignment reader, which saves elements in the following way:
  #[[1]] = number of seqs in the current fasta file,    [[2]] = list of geneIDs,                          [[3]] = list of sequence strings
  #[[1]] has 1 dimension,                               [[2]][i] allows you to go acress i'th geneID,     [[3]][i] allows you to access the i'th seq
  
  input_names = input_list[[i]]
  input_seqs = input_fasta[[3]]
  len_fasta = length(input_seqs)
  
  calcDist(input_seqs, len_fasta, input_names)
}

