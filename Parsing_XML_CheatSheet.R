#This program will parse MEME outputs and organize them into a large table of relevant information, allowing for easier sorting of motifs
#There are several output formats from MEME- txt and xml



install.packages("XML")
library("XML")

data <- xmlParse("/N/dc2/scratch/tlicknac/All_Aurelias_MEME_Results/pbi_10193.fasta/meme.xml")   #24 sequences

xml_data = xmlToList(data)

length(xml_data)  #return 5: training set, model, motifs, scanned_sites_summary, attributes
#Background A/T
xml_data$model$background_frequencies$alphabet_array  #contains frequency of each letter in entire sequence
xml_data$model$background_frequencies$alphabet_array[[1]][1]  #freq(A)
xml_data$model$background_frequencies$alphabet_array[[2]][1]  #freq(C)
xml_data$model$background_frequencies$alphabet_array[[3]][1]  #freq(G)
xml_data$model$background_frequencies$alphabet_array[[4]][1]  #freq(T)
#Position Score/Freq Matrix
length(xml_data$motifs$motif) #returns 5: scores, probabilities, regexs, contributing sites, attributes
length(xml_data$motifs$motif$scores$alphabet_matrix)  #returns length of motif ... 18 for us
xml_data$motifs$motif$scores$alphabet_matrix[[1]][1]$value$text  #score(A) for position 1
xml_data$motifs$motif$scores$alphabet_matrix[[1]][2]$value$text  #score(C) for position 1
xml_data$motifs$motif$scores$alphabet_matrix[[1]][3]$value$text  #score(G) ... 
xml_data$motifs$motif$scores$alphabet_matrix[[1]][4]$value$text  #score(T) ... 
xml_data$motifs$motif$scores$alphabet_matrix[[2]][1]$value$text  #score(A) for position 2
length(xml_data$motifs$motif$probabilities$alphabet_matrix)  #same as scores
xml_data$motifs$motif$probabilities$alphabet_matrix[[1]][1]$value$text  #Freq(A) for position 1
xml_data$motifs$motif$probabilities$alphabet_matrix[[1]][2]$value$text  #Freq(C) for position 1
xml_data$motifs$motif$probabilities$alphabet_matrix[[1]][3]$value$text  #Freq(G) ...
xml_data$motifs$motif$probabilities$alphabet_matrix[[1]][4]$value$text  #Freq(T) ..
xml_data$motifs$motif$probabilities$alphabet_matrix[[2]][1]$value$text  #Freq(A) for position 2
#Regular Expression
xml_data$motifs$motif$regular_expression  #Returns regex AS CHARACTER with flanking \n
#Flanking Region
length(xml_data$motifs$motif$contributing_sites)  #24 ... number of total sequences
xml_data$motifs$motif$contributing_sites[[1]][1]  #left flank for 1st seq
xml_data$motifs$motif$contributing_sites[[1]][3]  #right flank for 1st seq
xml_data$motifs$motif$contributing_sites[[2]][1]  #left flank for 2nd seq
xml_data$motifs$motif$contributing_sites[[2]][3]  #right flank for 2nd seq
#Consensus Sequence
xml_data$motifs[[1]]$.attrs[2]  #1st consensus motif
xml_data$motifs[[2]]$.attrs[2]  #2nd consensus motif
xml_data$motifs[[3]]$.attrs[2]  #...
xml_data$motifs[[4]]$.attrs[2]
xml_data$motifs[[5]]$.attrs[2]
#E_Value
xml_data$motifs[[1]]$.attrs[9]  #e_value for first motif
xml_data$motifs[[2]]$.attrs[9]  #...
xml_data$motifs[[3]]$.attrs[9]
xml_data$motifs[[4]]$.attrs[9]
xml_data$motifs[[5]]$.attrs[9]
#Additional Data on Motif Position in Each Sequence
length(xml_data$scanned_sites_summary)  #Returns 28 ... all of different lengths ranging from 1 to 6
length(xml_data$scanned_sites_summary[[1]])  #Returns 4: scanned site, scanned site, scanned site, attributes
length(xml_data$scanned_sites_summary[[2]])  #Returns 3: scanned site, scanned site, attributes ..... number of scanned sites varies based on motif hits in that sequence
xml_data$scanned_sites_summary[[1]]  #FIRST SEQUENCE
xml_data$scanned_sites_summary[[2]]  #SECOND SEQUENCE
xml_data$scanned_sites_summary[[1]][1]  #first motif hit (in first sequence): motif_id, strand, position, pvalue
xml_data$scanned_sites_summary[[1]][2]  #second motif hit (in first sequence)
xml_data$scanned_sites_summary[[2]][1]  #first motif hit (in second sequence): motif_id, strand, position, pvalue
xml_data$scanned_sites_summary[[2]][2]  #second motif hit (in second sequence)
xml_data$scanned_sites_summary[[1]][1]$scanned_site[1]  #Returns first motif_id in first sequence
xml_data$scanned_sites_summary[[1]][1]$scanned_site[2]  #Returns strand 
xml_data$scanned_sites_summary[[1]][1]$scanned_site[3]  #Returns position of first motif in first sequence
xml_data$scanned_sites_summary[[1]][1]$scanned_site[4]  #Returns p_value for first sequence
xml_data$scanned_sites_summary[[1]][2]$scanned_site[1]  #Returns second motif_id in first sequence
xml_data$scanned_sites_summary[[1]][2]$scanned_site[2]  #Returns strand
xml_data$scanned_sites_summary[[1]][2]$scanned_site[3]  
xml_data$scanned_sites_summary[[1]][2]$scanned_site[4]  
xml_data$scanned_sites_summary[[2]][1]$scanned_site[1]  #Returns first motif_id in second sequence



















