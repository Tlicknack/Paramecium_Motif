#This script will convert the modified GFF (containing only genes) into a GFF with only intergenic regions
#INPUT- paramecium-gene.tab
  #scaf genestart geneend geneID  
#OUPUT- Paramecium-intergenic.tab
  #scaf start end 5'+ 3'+ 5'- 3'-


species_gff = read.table("/home/tlicknac/Desktop/Paramecium_Genome_Data/Paramecium_GFF/ptet-gene.tab", sep="\t", as.is=T)                                           #CHANGE

#vscaf = unique(species_gff$V1)                                                         #If I want to go scaffold by scaffold
#uniq_scaf = length(vscaf)
#print(uniq_scaf)
#for(scaf in vscaf){
#  ts = species_gff[which(species_gff$V1==scaf), ]                                        #set variable ts equal to all lines in the gff that match the scaffold of interest (going 1 by 1)
#}

for(i in 1:nrow(species_gff)){
  
  output_table = data.frame(matrix(nrow=nrow(species_gff), ncol=7))

  if(species_gff$V1[i] == species_gff$V1[i+1]){
    five_minus = ""
    five_plus = ""
    three_minus = ""
    three_plus = ""
    row2 = species_gff[i+1, 1:5]                                  #All of row 2
    row1 = species_gff[i, 1:5]                                    #All of row 1
    start_row2 = row2[,2]                                          #gene start position is 2nd column
    end_row1 = row1[,3]                                            #gene end position is 3rd column
    row2_orient = row2[,4]                                        # + or - strand 
    row1_orient = row1[,4]
    row2_id = row2[,5]
    row1_id = row1[,5]
    
    if(row1_orient == "-"){
      five_minus = row1_id
    }
    if(row1_orient == "+"){
      three_plus = row1_id 
    }
    if(row2_orient == "-"){
      three_minus = row2_id
    }
    if(row2_orient == "+"){
      five_plus = row2_id
    }
    #Head to Head
      # - +
    #Head to Tail
      # - - 
    #Tail to Tail
      # + -
    #Tail to Head
      # + +

    colnames(output_table) = c("scaffold", "start_position", "end_position", "5'+strand", "5'-strand", "3'+strand", "3'-strand")
    output_table[i,1] = species_gff$V1[i]
    output_table[i,2] = start_row2
    output_table[i,3] = end_row1
    output_table[i,4] = five_plus 
    output_table[i,5] = five_minus
    output_table[i,6] = three_plus
    output_table[i,7] = three_minus
  } else {
    break
  }
}
