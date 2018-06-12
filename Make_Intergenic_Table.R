#DONE
#This script will convert the modified GFF (containing only genes) into a GFF with only intergenic regions
#INPUT- paramecium-gene.tab
  #scaf genestart geneend geneID  
#OUPUT- Paramecium-intergenic.tab
  #scaf start end 5'+ 3'+ 5'- 3'-
#There are 154 scaffolds in the Ptet genome with 1 gene... according to our annotation... this should get them too
  #Bash one liner: cat Paramecium_GFF/ptet-gene.tab | cut -f1 | sort -r | uniq -c | sort -n | grep " 1 " | wc -l

#source("http://bioconductor.org/biocLite.R")
#biocLite("seqinr")
library("seqinr")

species_gff = read.table("/home/tlicknac/Desktop/Paramecium_Genome_Data/Paramecium_GFF/ptet-gene.tab", sep="\t", as.is=T)                                           #CHANGE

vscaf = unique(species_gff$V1)                                                #If I want to go scaffold by scaffold
final_table = matrix(ncol=7)                                                  #Initiate final table that will be exported
colnames(final_table) = c("scaffold", "start_position", "end_position", "5'+strand", "5'-strand", "3'+strand", "3'-strand")        #Column names for final table 
#completed_rows = 0

for(scaf in vscaf){                     #iterate through scaffolds
  #print(scaf)
  ts = species_gff[which(species_gff$V1==scaf), ]                             #set variable ts equal to all lines in the gff that match the scaffold of interest (going 1 by 1)
  row_count = nrow(ts)                                                        #get count of rows for each scaffold
  output_table = matrix(nrow=nrow(ts), ncol=7)                                #Initiate table and values for genes that go into them ... will have to hard-wire nrow !
  
  #completed_rows = completed_rows + 1
  
  #First gene in each scaffold
  five_minus_1 = ""                                                           #Set values for each geneID position holder... not doing this will cause bug in lines 45-48
  five_plus_1 = ""
  three_minus_1 = ""
  three_plus_1 = ""
  
  row1 = ts[1, 1:5]                                                           #Get the first row of each scaffold; columns 1-5
  start_row1 = row1[,2]                                                       #Get the start of each gene position on row1
  end_row1 = row1[,3]                                                         #End position
  row1_orient = row1[,4]                                                      #Get gene's orientation
  row1_id = row1[,5]                                                          #Get geneID
  if(row1_orient == "+"){                                                     #If gene is on + strand, intergenic is 5' of plus. if on - strand, then it is the 3' of minus
    five_plus_1 = row1_id
  } else{
    three_minus_1 = row1_id
  }
  
  
  output_table[1,1] = ts$V1[1]                                                # Plug in values for row1
  output_table[1,2] = as.numeric(start_row1)
  output_table[1,3] = as.numeric(end_row1)
  output_table[1,4] = five_plus_1
  output_table[1,5] = five_minus_1
  output_table[1,6] = three_plus_1
  output_table[1,7] = three_minus_1
 
#----------------------------------------------------------------------------------------------------------------------------------------
  #Middle genes
  i = 2                                                             #start i at 2 ... will iterate through middle genes two rows at a time
  while(i < nrow(ts)){                                                                  
  
    five_minus_i = ""
    five_plus_i = ""
    three_minus_i = ""
    three_plus_i = ""
      
    rowi = ts[i-1, 1:5]                                             #All of row 1 for position i-1 (starts at 1)
    rowi2 = ts[i, 1:5]                                              #All of row 2 for position i (starts at 2)
    start_rowi2 = rowi2[,2]                                         #gene start position is 2nd column
    end_rowi = rowi[,3]                                             #gene end position is 3rd column ... middle is intergenic
    rowi2_orient = rowi2[,4]                                        # + or - strand of row i
    rowi_orient = rowi[,4]                                          # + or - strand of row i-1
    rowi2_id = rowi2[,5]                                            #geneID of row i
    rowi_id = rowi[,5]                                              #geneID of row i-1
    
    if(rowi_orient == "-"){
      five_minus_i = rowi_id
    }
    if(rowi_orient == "+"){
      three_plus_i = rowi_id 
    }
    if(rowi2_orient == "-"){
      three_minus_i = rowi2_id
    }
    if(rowi2_orient == "+"){
      five_plus_i = rowi2_id
    }
    #Head to Head
      # - +
    #Head to Tail
      # - - 
    #Tail to Tail
      # + -
    #Tail to Head
      # + +
    #Output values into table
    output_table[i,1] = ts$V1[i]
    output_table[i,2] = as.numeric(start_rowi2)
    output_table[i,3] = as.numeric(end_rowi)
    output_table[i,4] = five_plus_i 
    output_table[i,5] = five_minus_i
    output_table[i,6] = three_plus_i
    output_table[i,7] = three_minus_i
    
    i = i+1
  }
  #-------------------------------------------------------------------------------------------------------------------------------
  #Last gene
  five_minus_n = ""
  five_plus_n = ""
  three_minus_n = ""
  three_plus_n = ""
                                                                              #Do the same for the last row in each scaffold
  rown = ts[row_count, 1:5]
  start_rown = rown[,2]
  end_rown = rown[,3]
  rown_orient = rown[,4]
  rown_id = rown[,5]
  if(rown_orient == "+") {
    five_plus_n = rown_id
  } else {
    three_minus_n = rown_id
  }
  
  output_table[row_count,1] = ts$V1[row_count]
  output_table[row_count,2] = as.numeric(start_rown)
  output_table[row_count,3] = as.numeric(end_rown)
  output_table[row_count,4] = five_plus_n
  output_table[row_count,5] = three_plus_n
  output_table[row_count,6] = five_minus_n
  output_table[row_count,7] = three_minus_n
  
  final_table = rbind(final_table, output_table)
}

write.table(final_table, sep="\t", row.names=FALSE, file = "ptet-intergenic.tab")
