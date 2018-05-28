#This script will convert the modified GFF (containing only genes) into a GFF with only intergenic regions
#INPUT- paramecium-gene.tab
  #scaf genestart geneend geneID  
#OUPUT- Paramecium-intergenic.tab
  #scaf start end 5'+ 3'+ 5'- 3'-
#There are 154 scaffolds in the Ptet genome with 1 gene... according to our annotation... this should get them too
  #Bash one liner: cat Paramecium_GFF/ptet-gene.tab | cut -f1 | sort -r | uniq -c | sort -n | grep " 1 " | wc -l


species_gff = read.table("/home/tlicknac/Desktop/Paramecium_Genome_Data/Paramecium_GFF/ptet-gene.tab", sep="\t", as.is=T)                                           #CHANGE

vscaf = unique(species_gff$V1)                                                         #If I want to go scaffold by scaffold
uniq_scaf = length(vscaf)
#print(uniq_scaf)
                                       
 
output_table = data.frame(matrix(nrow=50000, ncol=7))                       #Initiate table and values for genes that go into them ... will have to hard-wire nrow !

for(scaf in vscaf){
  #print(scaf)
  ts = species_gff[which(species_gff$V1==scaf), ]                                         #set variable ts equal to all lines in the gff that match the scaffold of interest (going 1 by 1)
  row_count = nrow(ts)  
  completed_rows = 
  
  #First gene in each scaffold
  five_minus_1 = ""
  five_plus_1 = ""
  three_minus_1 = ""
  three_plus_1 = ""
  
  row1 = ts[1, 1:5]
  start_row1 = row1[,2]
  end_row1 = row1[,3]
  row1_orient = row1[,4]
  row1_id = row1[,5]
  if(row1_orient == "+"){
    five_plus_1 = row1_id
  } else{
    three_minus_1 = row1_id
  }
  
  colnames(output_table) = c("scaffold", "start_position", "end_position", "5'+strand", "5'-strand", "3'+strand", "3'-strand")
  output_table[1,1] = ts$V1[1]
  output_table[1,2] = start_row1
  output_table[1,3] = end_row1
  output_table[1,4] = five_plus_1
  output_table[1,5] = five_minus_1
  output_table[1,6] = three_plus_1
  output_table[1,7] = three_minus_1
 
#----------------------------------------------------------------------------------------------------------------------------------------
  #Middle genes
  i = 2
  while(i < nrow(ts)){                                                                  
  
    five_minus_i = ""
    five_plus_i = ""
    three_minus_i = ""
    three_plus_i = ""
      
    rowi = ts[i-1, 1:5]                                             #All of row 2
    rowi2 = ts[i, 1:5]                                              #All of row 1
    start_rowi2 = rowi2[,2]                                         #gene start position is 2nd column
    end_rowi = rowi[,3]                                             #gene end position is 3rd column
    rowi2_orient = rowi2[,4]                                        # + or - strand 
    rowi_orient = rowi[,4]
    rowi2_id = rowi2[,5]
    rowi_id = rowi[,5]
    
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
    output_table[i,2] = start_rowi2
    output_table[i,3] = end_rowi
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
  output_table[row_count,2] = start_rown
  output_table[row_count,3] = end_rown
  output_table[row_count,4] = five_plus_n
  output_table[row_count,5] = three_plus_n
  output_table[row_count,6] = five_minus_n
  output_table[row_count,7] = three_minus_n
}

write.table( output_table, sep="\t", row.names=FALSE, file = "ptet-intergenic.tab")
