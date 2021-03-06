#This will take in a table containing positions with high scores from PWM across genomes of interest
#INPUT:
  #DF with "Consensus_Motif", "POFF_Row", "Species", "Scaffold", "Position", "PWM_Score"
  #Intergenic GFFs with: "scaffold"      "start_position"        "end_position"  "5'+strand"     "5'-strand"     "3'+strand"     "3'-strand"

#OUTPUT
  #DF with the above + genic features, position to start

#MAIN

#Load intergenic gffs
lgff = list()
lgff[["pbi"]] = read.table("/N/u/tlicknac/Carbonate/Paramecium_GFF/pbi-intergenic.tab", header=T, sep="\t")
lgff[["psex"]] = read.table("/N/u/tlicknac/Carbonate/Paramecium_GFF/psex-intergenic.tab", header=T, sep="\t")
lgff[["ptet"]] = read.table("/N/u/tlicknac/Carbonate/Paramecium_GFF/ptet-intergenic.tab", header=T, sep="\t")


#Load pwm scoring data frame
input_file = "pwm-distribution-duplicate_8450-2.csv"                                  ##CHange this for each new PWM file
file_directory_path = "/N/dc2/scratch/tlicknac/Data-PWM_Distribution/"
input = paste(file_directory_path, input_file, sep="")
pwm_table = read.csv(input, header=T, as.is=T)
new_column = data.frame(matrix(ncol=5, nrow=nrow(pwm_table)))
colnames(new_column) = c("Intergenic", "FivePrimeGenePlusStrand", "FivePrimeGeneMinusStrand", "DistancetoFivePlus", "DistanceToFiveMinus")
pwm_table = cbind(pwm_table, new_column)

for(i in 1:nrow(pwm_table)){ 
  #Load PWMs and GFFs
  pwm_row = pwm_table[i,]
  
  if(pwm_row$Species == "ptetraurelia_mac_51.fa"){
    gff_tab = lgff[["ptet"]]
  } 
  if(pwm_row$Species == "biaurelia_V1-4_assembly_v1.fa"){
    gff_tab = lgff[["pbi"]]
  } 
  if(pwm_row$Species == "sexaurelia_AZ8-4_assembly_v1.fasta"){
    gff_tab = lgff[["psex"]]
  }
  
  #For line i, only proceed if the scaffold of the PWM hit matches the scaffold in the GFF... not all scaffolds are in both annotation and assembly
  if(pwm_row$Scaffold %in% gff_tab$scaffold){ 
    #If the scaffold is in the annotation, get the first row for which end of the intergenic is larger than the position of the high PWM score
    gff_current_scaf = gff_tab[which(gff_tab$scaffold ==  pwm_row$Scaffold),]                       
    firstrow = gff_current_scaf[which(gff_current_scaf$end_position > pwm_row$Position),][1,]   
    
    #Some hits won't be smaller than any end position for that scaffold because it'll be between the last annotated gene and the telomere... get the final row of the GFF for current scaffold
    if(is.na(firstrow$scaffold)){ 
      lastrow = gff_current_scaf[nrow(gff_current_scaf),]  
      difference = lastrow$end_position
      
      #If the PWM hit is within 200bps
      if(difference < 200){
      
        if(is.na(gff_current_scaf$scaffold[1]) == F){
          intergenic = "Y"
          five_plus = as.character(lastrow$X5..strand)
          five_minus = as.character(lastrow$X5..strand.1)
        }
      } else{
        intergenic = "N"
        five_plus = NA
        five_minus = NA
      }
    } else{ 
    
      if(pwm_row$Position > firstrow$start_position){
        intergenic = "Y"
        five_plus =  as.character(firstrow[4]$X5..strand)
        five_minus = as.character(firstrow[5]$X5..strand)
      } else{
        intergenic = "N"
        five_plus = NA
        five_minus = NA
      }
     } 
    if(is.na(five_plus) == FALSE){
      dist_to_plus = firstrow$end_position - pwm_row$Position
    } else{
      dist_to_plus = NA
    }
  
    if(is.na(five_minus) == FALSE){
      dist_to_minus = pwm_row$Position- firstrow$start_position 
    } else{
      dist_to_minus = NA
    }
      cat("Writing to line: ", i, "\n")
      pwm_table$Intergenic[i] = intergenic
      pwm_table$FivePrimeGenePlusStrand[i] = five_plus
      pwm_table$FivePrimeGeneMinusStrand[i] = five_minus
      pwm_table$DistancetoFivePlus[i] = dist_to_plus
      pwm_table$DistanceToFiveMinus[i] = dist_to_minus
  } 
}
pwm_table$X = NULL                                                  #Might need this if old row names are added as a column
length(which(is.na(pwm_table$Intergenic)))                          #gives the number of PWM hits nowhere near any annotation
pwm_table = pwm_table[-which(is.na(pwm_table$Intergenic)),]         #this will remove all rows that didn't receive any annotation
output_file_name = gsub(".csv", "_With-Annotation.csv", input_file) #modify file name
write.csv(pwm_table, file=output_file_name, row.names=FALSE)  
