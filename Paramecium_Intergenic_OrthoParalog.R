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

#Before you execute, make sure you're in my Carbonate directory (or just change the location when importing files) and install seqinr with: install.packages("seqinr")

library("seqinr")

extractPoffIDs = function(bi_poff, tet_poff, sex_poff){
	
	bi_homologs = 
	sex_homologs = 
	tet_homologs = 
}

#Main function which will run and loop all other functions
main = function(){
	# Working Directory: /N/u/tlicknac/Carbonate/
	bi_poff = read.table("~/Paramecium_POFF/biaurelia_group.poff", header=TRUE, sep="\t", row.names="id")
	tet_poff = read.table("tetraurelia_group.poff", header=TRUE, sep="\t", row.names="id")
	sex_poff = read.table("sexaurelia_group.poff", header=TRUE, sep="\t", row.names="id")
	
	pbi_gff = read.table("~/Paramecium_GFF/pbi-gene.tab", sep="\t", row.names="id")
	pdec_gff = read.table("pdec-gene.tab", sep="\t", row.names="id")
	pdodec_gff = read.table("pdodec-gene.tab", sep="\t", row.names="id")
	pjen_gff = read.table("pjen-gene.tab", sep="\t", row.names="id")
	pnov_gff = read.table("pnov-gene.tab", sep="\t", row.names="id")
	poct_gff = read.table("poct-gene.tab", sep="\t", row.names="id")
	ppent_gff = read.table("ppent-gene.tab", sep="\t", row.names="id")
	pprim_gff = read.table("pprim-gene.tab", sep="\t", row.names="id") 
	pquad_gff = read.table("pquad-gene.tab", sep="\t", row.names="id")
	psept_gff = read.table("psept-gene.tab", sep="\t", row.names="id")
	psex_gff = read.table("psex-gene.tab", sep="\t", row.names="id")
	pson_gff = read.table("pson-gene.tab", sep="\t", row.names="id")
	ptet_gff = read.table("ptet-gene.tab", sep="\t", row.names="id")
	ptre_gff = read.table("ptre-gene.tab", sep="\t", row.names="id")
	
	pbi_fa = read.fasta("~/Paramecium_FASTA/biaurelia_V1-4_assembly_v1.fa", as.string=TRUE)
	pdec_fa = read.fasta("pdecaurelia_mac_223_v1.0.fa", as.string=TRUE)
	pdodec_fa = read.fasta("pdodecaurelia_mac_274_v1.0.fa", as.string=TRUE)
	pjen_fa = read.fasta("pjenningsi_mac_M_v1.0.fa", as.string=TRUE)
	pnov_fa = read.fasta("pnovaurelia_mac_TE_v1.0.fa", as.string=TRUE)
	poct_fa = read.fasta("poctaurelia.fasta", as.string=TRUE)
	ppent_fa = read.fasta("ppentaurelia_mac_87_v1.0.fa", as.string=TRUE)
	pprim_fa = read.fasta("primaurelia_Ir4-2_assembly_v1.fasta", as.string=TRUE)
	pquad_fa = read.fasta("pquadecaurelia_mac_NiA_v1.0.fa", as.string=TRUE)
	psept_fa =read.fasta("pseptaurelia_mac_38_v1.0.fa", as.string=TRUE)
	psex_fa = read.fasta("sexaurelia_AZ8-4_assembly_v1.fasta", as.string=TRUE)
	pson_fa = read.fasta("psonneborni_mac_ATCC_30995_v1.0.fa", as.string=TRUE)
	ptet_fa = read.fasta("ptetraurelia_mac_51.fa", as.string=TRUE)
	ptre_fa = read.fasta("ptredecaurelia_209_AP38_filtered.fa", as.string=TRUE)
	
}	
	
