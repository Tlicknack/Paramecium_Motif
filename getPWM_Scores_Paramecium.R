#This program will use Jeff's C++ code to scan a genome of interest, and give scores to each position based upon a given PWM for a DNA motif
#These PWMs come from the Paraortholog study
#With these scores, we can look at distribution of scores across a scaffold for an organism, compare between organisms, see what genes these motifs lie near
#INPUT:
  #C++ code- getPWM_Scores.cpp
  #PWMs, or list of PWMs
  #FASTA files for analysis
#OUTPUT:
  #Dataframe with consensus motif, species, scaffold, position, score ONLY for positions where score>0

library(Rcpp)
sourceCpp("/N/u/tlicknac/Carbonate/getPWM_Scores.cpp")
library(seqinr)


consensusFromPWM = function(pwm_tab){
  motif = c()
  i=1
  max_cols = max.col(pwm_tab, ties.method="first")
  for(biggest in max_cols){
    i=i+1
    motif[i] = names(pwm_tab[1,biggest]) 
  }
  motif = paste(motif, collapse="")
  motif = substr(motif, 3, nchar(motif))
  return(motif)
}


fasta_directory = "/N/u/tlicknac/Carbonate/Paramecium_FASTA/"                                                         #location of PWMs
FASTAs = list.files(fasta_directory, recursive=F)                                                                     #List all files in directory
FASTAs = FASTAs[grep("*.fa", FASTAs)]                                                                                 #Rm all files that arent fasta files *this will capture .fasta and .fa

pwm_directory = "/N/dc2/scratch/tlicknac/All_Aurelias_PWMs/"                                                          #location of PWMs
#PWMs = list.files(pwm_directory, recursive=F)                                                                         #List all files in directory
PWMs = c("pbi_7605.fasta_PWM.tab", "pbi_736.fasta_PWM.tab", "pbi_4750.fasta_PWM.tab", "pbi_16238.fasta_PWM.tab", "pbi_676.fasta_PWM.tab")              #first attempt at this... 5 best motif

firstloop=TRUE

for(assembly in FASTAs){                                                                                              #Iterate through fasta file names
  cat("Starting with species: ", assembly, "\n")
  a = paste(fasta_directory, assembly, sep="")                                 
  fasta = read.fasta(a, as.string=T)                                                                                  #read in fasta 1 at a time... must preserve memory
  
  for(scaf in fasta){                                                                                                 #Iterate through each scaffold
    cat("Onto scaffold: ", getName(scaf), " of species: ", assembly, "\n")
    scaf_seq = as.character(scaf[1])                                                                                  #get sequence
    
    for(pwm in PWMs){
      pwm_file_path = paste(pwm_directory, pwm, sep="") 
      t = read.table(pwm_file_path, header=T, as.is=T, sep=",")                                                       #read in PWM as table
      pwm_tab = as.matrix(t[ , 2:5])                                      #####Change to c(A,  T ... )######
      vscores = getScores(pwm_tab, scaf_seq) #C++ code
      vgood = which(vscores>-100)                                #POSITION ####CHANGE CUTOFF IF NEEDED. 0 CHOSEN AS DEFAULT... only a handful of positions will meet this criterion####
      
      if(length(vgood)>0){
        
        for(position in vgood){
        
          if(firstloop==TRUE){
            finaldf = data.frame(c( consensusFromPWM(pwm_tab), pwm, assembly, getName(scaf), position, vscores[position] ))
            finaldf = t(finaldf)
            colnames(finaldf) = c("Consensus_Motif", "POFF_Row", "Species", "Scaffold", "Position", "PWM_Score")
            firstloop=FALSE
          
          } else{
            tmpdf = data.frame(c( consensusFromPWM(pwm_tab), pwm, assembly, getName(scaf), position, vscores[position] ))  #all the info we need
            tmpdf = t(tmpdf)
            colnames(tmpdf) = c("Consensus_Motif", "POFF_Row", "Species", "Scaffold", "Position", "PWM_Score")
            finaldf = rbind(finaldf, tmpdf)
            cat("Wrote to finaldf", "\n")
          }
        }
      }
    }
  }
  #rm(setdiff(ls(), c("FASTAs", "PWMs", "firstloop",  "getScores", "fasta_directory", "pwm_directory")))               #Being frugle, rm'ing unessesary objects
}
write.csv(finaldf, file="test_5motifs_15spp_Rownames.csv")
rownames(finaldf) = NULL
write.csv(finaldf, file="test_5motifs_15spp.csv")


#Create data frame with 3 columns
  #scaffold number, position, intergenic?
