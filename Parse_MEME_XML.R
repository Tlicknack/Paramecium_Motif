#This program will parse MEME xml outputs to retrieve important information:
  #Consensus sequence/Regex
  #E_value
  #PWM
  #Background AT content
  #Motif position within promoters
#INPUT:
  #meme.xml (from a directory containing)


install.packages("XML")
library("XML")


getATcontent = function(xml_data){
  lAT = list()                                                                #initiate list
  len = length(xml_data$model$background_frequencies)                         #get number of promoters
  for(j in 1:len){                                                             
    a = xml_data$model$background_frequencies$alphabet_array[[1]][j]              #get frequency of As in each sequence
    t = xml_data$model$background_frequencies$alphabet_array[[4]][j]              #get frequency of Ts in each sequence
    at = as.integer(a)+as.integer(t)
  }
}

getConsensusSeq = function(xml_data){
  lseq_evalues = list()                                                       #initiate list
  len = length(xml_data$motifs)                                               #get number of motifs
  for(i in 1:len){                                                            #interate through that number
    consensus_seq = as.character(xml_data$motifs[[i]]$.attrs[2])              #get consensus sequence; located in motifs[1:len] then the 2nd attribute
    e_value = as.character(xml_data$motifs[[i]]$.attrs[9])                    #get e_value; located in motifs[1:len] then the 9th attribute
    seq_evalue = c(consensus_seq, e_value)                                    #concatenate the consensus and e_value
    lseq_evalues[[i]] = seq_evalue                                            #add them to a list for storage
  }
  return(lseq_evalues)                                                        #return a list, which should have 5 tuples
}

#MAIN ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
directory = "/N/dc2/scratch/tlicknac/All_Aurelias_MEME_Results/"                                                      #location of meme outputs
meme_files = list.files(directory, recursive=F)                                                                       #List all subdirectories containing meme outputs

firstOne = T                                                                                                          #Used to initialize data frame 
final_dataframe = 0                                                                                                   #^

for(subdir in meme_files){
  cat("Working on file: ", subdir, "\n")
  meme = paste(directory, subdir, "/meme.xml", sep="")                                                                #Concatenate path to file name
  data = xmlParse(meme)                                                                                               #Get xml file
  xml_data = xmlToList(data)                                                                                          #Convert xml format to an R list
  #Get Consensus Sequence
  lseq_evalues = getConsensusSeq(xml_data)                                                                            #Function to get list of tuples: each tuple is consensus sequence and e_value
  df_seq_evalues = data.frame(matrix(unlist(lseq_evalues), nrow=length(lseq_evalues), byrow=T), stringsAsFactors=F)   #convert the list to a df
  colnames(df_seq_evalues) = c("motif", "e_value")                                                                    #Give the df colnames
  if(firstOne==T){                                                                                                    #Our strategy initiates final df using the first motif in the dir/
    final_dataframe = df_seq_evalues                                                                                  #fina_dataframe will now have 5 occupied rows
    firstOne = F                                                                                                      #set this boolean to false so we can continue appending to final_dataframe
  } else {
    final_dataframe = rbind(final_dataframe, df_seq_evalues)                                                          #append until completion
  }
  
  #Get AT Content
  #lAT = getATcontent(xml_data)
}
write.csv(final_dataframe, file="MEMExmlParaOrtholog.csv")














