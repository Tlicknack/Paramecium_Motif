#This program will group together the DNA sequence upstream of ortho-parologs and output them into a fasta file. This file can then be used for motif searching (e.g. MEME)
#The program requires:
  #FASTA files for each species 
  #modified GFFs which contains only the scaffold/chromosome number, gene start, gene end, strand, and geneID for each species
  #POFF tables of orthoparalogs which contains number of species,number of genes, Alg.-Conn, and a variable number of species names under which geneIDs will be located
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#source("http://bioconductor.org/biocLite.R")
#biocLite("seqinr")
#install.packages("stringi")
library("seqinr")


extractUpstream = function(sp, geneID, lGFF, lfasta, fname){                          #Pass species name, geneID, list of GFFs and list of FASTAs for that family

  grep_geneID = which(lGFF[[sp]]$V5 == geneID, arr.ind=TRUE)
  cat("Grep= ", grep_geneID, "\n")
  #grep_geneID = as.integer(grep(geneID, lGFF[[sp]][,5]))                             #search each gff for geneID, return column number
  row_geneID = lGFF[[sp]][grep_geneID, ]
  if(row_geneID[,4] == "+"){                                      #get Gene's starting position, depending upon its strand 
    gene_start = row_geneID[,2]   
  } else{
    gene_start = row_geneID[,3]
  }
  scaf = row_geneID[,1]                                             #find the scaffold holding the gene 
  cat("scaf= ", scaf, ".....", "gene start= ", gene_start, ".....", "strand= ", row_geneID[,4], "\n")
  scafSeq=as.character(lfasta[[sp]][scaf])
  if(row_geneID[,4] == "+"){                                       #extract the correct upstream sequence. if +, then its just 200bp from start of gene, if -, then its 200 from end of gene
    upstream = substr(scafSeq, gene_start, gene_start-200)
  } else {
    upstream = substr(scafSeq, gene_start, gene_start+200)
    upstream = c2s(rev(comp(s2c(upstream))))
  }
  return(upstream)
  #cat("Upstream sequence for ", geneID, "has been done successfully. \n")
}
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


getPOFFid = function(vsp, POFF, lGFF, lFasta) {                                     #Pass POFF and list of species names
  
  i = 1                                                               #counter for iterating through poff
  for(i in 1:nrow(POFF)) {                                            #iterate through POFF
    fname = paste(names(POFF[4]), "_", i, ".fasta", sep="")       #file name that resets each row
    cat("Parsing line ", i, "\n")           #fun printed message
    
    vupstream = c()
    for(sp in vsp) {                                                  #iterate through species names, which are identifiers
      sp_family = as.character(POFF[i,sp])                            #get proteinIDs from each column in each row
      vsplit = strsplit(sp_family, ",")[[1]]                          #split each proteinID as it comes out, some will have commas, some won't. Will do this for each species as we iterate
      
      for(geneID in vsplit){                                          #iterate through 1 or 2 proteinIDs (depending on if they needed to be split)

        if(geneID != "*"){
          cat("Now extracting sequence upstream of ", geneID, " (which is in", sp, ")\n")                     
          upstreamSeq = extractUpstream(sp, geneID, lGFF, lFasta, fname)
          cat("Sequence extracted: '", upstreamSeq, "'\n")
          vupstream[geneID] = upstreamSeq                     #changed from append to set each position, will be able to iterate through names in writeVectorAsFasta fxn
        }
      }
    }
    writeVectorAsFasta(vupstream, fname)
    #cat("line", i, "done.\n-----------------------\n")      #print that line one is complete
  }
}
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


writeVectorAsFasta <- function(vseq, fname){                          #Pass list of sequences from extractUpstream and name of file. This will write to fasta file as output
  #cat("FNAME: ", fname, "\n")
  #cat("", file=fname)                                         # Just to reset the file to empty if it already exists
  print(vseq)
  j = 0
  for(id in names(vseq)){
    j = j + 1
    #cat("Writing ", id ," to file '", fname, "'")
    cat(">", id, "\n", vseq[[j]], "\n", sep="", file=fname, append=T)       #dirty fasta format. will write to current directory
  }
}
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


#MAIN:                             #WD= Desktop/Paramecium_Genome_Data/Paramecium_Upstream_FASTA. ALl other files are ../Paramecium_whatever

  bi_poff = read.table("../Paramecium_POFF/head_bi_poff", header=TRUE, sep="\t", comment.char="$")  #read POFFs, change comment character so first line is read as header
  tet_poff = read.table("../Paramecium_POFF/head_tet_poff", header=TRUE, sep="\t", comment.char="$")
  sex_poff = read.table("../Paramecium_POFF/head_sex_poff", header=TRUE, sep="\t", comment.char="$")
  colnames(bi_poff) = gsub(".faa", "", colnames(bi_poff))     #fix column names
  colnames(tet_poff) = gsub(".faa", "", colnames(tet_poff)) 
  colnames(sex_poff) = gsub(".faa", "", colnames(sex_poff)) 
  
  vsp_bi = colnames(bi_poff)[4:ncol(bi_poff)]                 #get list of spcies names for each family... just in case
  vsp_tet = colnames(tet_poff)[4:ncol(tet_poff)] 
  vsp_sex = colnames(sex_poff)[4:ncol(sex_poff)] 
#----------------------------------------
  lbi_GFF = list()
  lsex_GFF = list()
  ltet_GFF = list()
  
  for(sp_bi in vsp_bi){
    bi_gff_name = paste("../Paramecium_GFF/", sp_bi, "-gene.tab", sep="")             #name of GFFs all look the same, so use list of names to create file names
    bi_gff_load = read.table(bi_gff_name, sep="\t")                       #load each file given its name
    lbi_GFF[[sp_bi]] = bi_gff_load                                         #add each file to list of GFFs with name of species as identifier
  }
  for(sp_tet in vsp_tet){
    tet_gff_name = paste("../Paramecium_GFF/", sp_tet, "-gene.tab", sep="") 
    tet_gff_load = read.table(tet_gff_name, sep="\t")
    ltet_GFF[[sp_tet]] = tet_gff_load
  }
  for(sp_sex in vsp_sex){
    sex_gff_name = paste("../Paramecium_GFF/", sp_tet, "-gene.tab", sep="") 
    sex_gff_load = read.table(sex_gff_name, sep="\t")
    ltet_GFF[[sp_sex]] = sex_gff_load    
  }
#---------------------------------------
  lbi_FASTA = list()
  lbi_FASTA[["pbi"]] = read.fasta("../Paramecium_FASTA/biaurelia_V1-4_assembly_v1.fa", as.string=TRUE)
  lbi_FASTA[["ppent"]] = read.fasta("../Paramecium_FASTA/ppentaurelia_mac_87_v1.0.fa", as.string=TRUE)
  lbi_FASTA[["pprim"]] = read.fasta("../Paramecium_FASTA/primaurelia_Ir4-2_assembly_v1.fasta", as.string=TRUE)
  lbi_FASTA[["pquad"]] = read.fasta("../Paramecium_FASTA/pquadecaurelia_mac_NiA_v1.0.fa", as.string=TRUE)
  lbi_FASTA[["ptre"]] = read.fasta("../Paramecium_FASTA/ptredecaurelia_209_AP38_filtered.fa", as.string=TRUE)
  lbi_FASTA[["pnov"]] = read.fasta("../Paramecium_FASTA/pnovaurelia_mac_TE_v1.0.fa", as.string=TRUE)
  getPOFFid(vsp_bi, bi_poff, lbi_GFF, lbi_FASTA)
  lsex_FASTA = list()
  lsex_FASTA[["psex"]] = read.fasta("../Paramecium_FASTA/sexaurelia_AZ8-4_assembly_v1.fasta", as.string=TRUE)
  lsex_FASTA[["pjen"]] = read.fasta("../Paramecium_FASTA/pjenningsi_mac_M_v1.0.fa", as.string=TRUE)
  lsex_FASTA[["pson"]] = read.fasta("../Paramecium_FASTA/psonneborni_mac_ATCC_30995_v1.0.fa", as.string=TRUE)
  getPOFFid(vsp_sex, sex_poff, lsex_GFF, lsex_FASTA)
  ltet_FASTA = list()
  ltet_FASTA[["pdec"]] = read.fasta("../Paramecium_FASTA/pdecaurelia_mac_223_v1.0.fa", as.string=TRUE)
  ltet_FASTA[["pdodec"]] = read.fasta("../Paramecium_FASTA/pdodecaurelia_mac_274_v1.0.fa", as.string=TRUE)
  ltet_FASTA[["poct"]] = read.fasta("../Paramecium_FASTA/poctaurelia.fasta", as.string=TRUE)
  ltet_FASTA[["psept"]] =read.fasta("../Paramecium_FASTA/pseptaurelia_mac_38_v1.0.fa", as.string=TRUE)
  ltet_FASTA[["ptet"]] = read.fasta("../Paramecium_FASTA/ptetraurelia_mac_51.fa", as.string=TRUE)
  getPOFFid(vsp_tet, tet_poff, ltet_GFF, ltet_FASTA)
  