#Goal: To extract 5' and 3' intergenic DNA of groups of OrthoParalogous genes in the Paramecium Aurelia Species Complex
#Due to the nature of their shared evolutionary history, we must group together species based on their location in 1 of 3 clades
	#Tet= Tetraurelia, Octaurelia, Septaurelia, Decaurelia, Dodecaurelia
	#Bi= Biaurelia, Novaurelia, Primaurelia, Pentaurelia, Tredecaurelia, Quadecaurelia
	#Sex= Sexaurelia, Jenningsi, Sonnonborni
#POFF files contain information about genes that are orthoparalogs of each other within each clade
	#All POFF files have 3 columns for "# of species", "Genes", "Alg.-Conn" + more columns containing GeneIDs for each species
		#Bi clade will have 6 more columns; 1 for each species
		#Tet clade will have 5 more columns
		#Sex will have 3 more columns
#TAB files are modified from GFF to get only genes, their positions on each scaffold, their strand orientation, and their GeneID
	#This GeneID is the same as the one in the POFF, so we have to match them and retrieve their location in the genome
#FASTA files contain scaffold name and sequence 
	#Using the scaffoldID, positions, and orientation, we must extract the upstream and downstream sequence of each gene and put them into files
	#These intergenic regions will be used for MEME software to search for motifs


Download_Files = function(

	
#Main function which will run and loop all other functions
main = function(args){
	
}	
	
args = commandArgs(trailingOnly = TRUE)
	main(args)