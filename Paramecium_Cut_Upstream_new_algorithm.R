#This program will group together the DNA sequence upstream of ortho-parologs and output them into a fasta file. This file can then be used for motif searching (e.g. MEME)
#The program requires:
  #FASTA files for each species 
  #modified GFFs which contains only the scaffold/chromosome number, gene start, gene end, strand, and geneID for each species
  #POFF tables of orthoparalogs which contains number of species,number of genes, Alg.-Conn, and a variable number of species names under which geneIDs will be located
#FILES ARE OUTPUTTED IN THE FOLLOWING WAY:
  #pbi_i.fasta
  #pjen_i.fasta
  #pdec_i.fasta
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#source("http://bioconductor.org/biocLite.R")
#biocLite("seqinr")
library("seqinr")


extractUpstream = function(sp, geneID, lGFF, lfasta, fname){      #Pass species name, geneID, list of (modified) GFFs and list of FASTAs for that family

  grep_geneID = which(lGFF[[sp]]$V5 == geneID, arr.ind=TRUE)      #go into the 5th column of the GFF for the current species, return the row number which matches the current geneID
  #cat("Grep= ", grep_geneID, "\n")
  row_geneID = lGFF[[sp]][grep_geneID, ]                          #return the contents of the entire row where the above match was found
  if(row_geneID[,4] == "+"){                                      #get Gene's starting position, depending upon its strand 
    gene_start = row_geneID[,2]   
  } else{
    gene_start = row_geneID[,3]
  }
  scaf = as.character(row_geneID[,1])                             #return the name of the scaffold holding the gene 
  #cat("scaf= ", scaf, ".....", "gene start= ", gene_start, ".....", "strand= ", row_geneID[,4], "\n")
  scafSeq = as.character(lfasta[[sp]][scaf])                      #return the entire sequence of the scaffold which holds the current gene
  if(row_geneID[,4] == "+"){                                      #extract the correct upstream sequence. if +, then its just 200bp from start of gene, if -, then its 200 from end of gene
    upstream = substr(scafSeq, gene_start-200, gene_start)
  } else {
    upstream = substr(scafSeq, gene_start, gene_start+200)
    upstream = c2s(rev(comp(s2c(upstream))))                      #must reverse complement the gene on the minus strand. Seqinr has a strange way of doing this
  }
  return(upstream)
  #cat("Upstream sequence for ", geneID, "has been done successfully. \n")
}
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


getPOFFid = function(vsp, POFF, lGFF, lfasta) {                           #Pass vector of species names, POFF of family of interest, list of GFFs and FASTAs for that family
  
  i = 1                                                                   #counter for iterating through poff
  for(i in 1:nrow(POFF)) {                                                #use i to iterate through each row number
    fname = paste(names(POFF[4]), "_", i, ".fasta", sep="")               #return file name based on the species name in column 4 (all POFFs have at least 6 columns)
    #cat("Parsing line ", i, "\n")           #fun printed message
    vupstream = c()                                                       #create empty vector to store all upstream sequences for genes in this row
    for(sp in vsp) {                                                      #iterate through species names
      sp_family = as.character(POFF[i,sp])                                #get proteinIDs from each row (based on iterator) and column (based on species name which is the header)
      vsplit = strsplit(sp_family, ",")[[1]]                              #some columns have two geneIDs, must split them on the comma that separates them 
      for(geneID in vsplit){                                              #iterate through 1 or 2 geneIDs (if there's only 1, then it's fine too)
        if(geneID != "*"){                                                #only proceed if the geneID != '*' ... *'s are used when no ortho-paralog were found. It returns no upstream sequence
          #cat("Now extracting sequence upstream of ", geneID, " (which is in", sp, ")\n")     #might as well remove*'s before running proceeding function                
          upstreamSeq = extractUpstream(sp, geneID, lGFF, lfasta, fname)
          #cat("Sequence extracted: '", upstreamSeq, "'\n")
          vupstream[geneID] = upstreamSeq                                 #add the upstream sequence into a vector with the position given by the geneID. can easily call that position w/ geneID
        }
      }
    }
    if(length(vupstream) >= 4) {                                          #only create a fasta of these sequences if there are at least 4. MEME will not be effective with 3 or less
      writeVectorAsFasta(vupstream, fname)
      cat("line", i, "done.\n-----------------------\n")                  #print that line row is complete ... I'll keep this print-out
    }
    }
}
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


writeVectorAsFasta <- function(vseq, fname){                          #Pass vector of sequences from extractUpstream and name of file. This will write to fasta file as output
  #cat("FNAME: ", fname, "\n")
  cat("", file=fname)                                                 #reset the file to empty if it already exists
  #print(vseq)
    for(id in names(vseq)){                                           #iterate through geneIDs which are names for each corresponding sequence
      #cat("Writing ", id ," to file '", fname, "'")
      seq_len = nchar(vseq[[id]])                                     #return length of sequence for each geneID
      if(seq_len > 10){                                               #only write to fasta if the length is greater than 10 (just in case something slipped by)
        cat(">", id, "\n", vseq[[id]], "\n", sep="", file=fname, append=T)       
        #write to fasta file. >geneID \n "sequence of that geneID". name the file after the species family and row number, append until this row is finished (fname will change next row)
      }
    }
}
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


#MAIN:                             #WD= Desktop/Paramecium_Genome_Data/Paramecium_Upstream_FASTA. ALl other files are ../Paramecium_whatever

  bi_poff = read.table("../Paramecium_POFF/biaurelia_group.poff", header=TRUE, sep="\t", comment.char="$", as.is=T)  #read POFFs, change comment character so first line is read as header
  tet_poff = read.table("../Paramecium_POFF/tetraurelia_group.poff", header=TRUE, sep="\t", comment.char="$", as.is=T)
  sex_poff = read.table("../Paramecium_POFF/sexaurelia_group.poff", header=TRUE, sep="\t", comment.char="$", as.is=T)
  colnames(bi_poff) = gsub(".faa", "", colnames(bi_poff))                 #fix column names. should just be a species identifier
  colnames(tet_poff) = gsub(".faa", "", colnames(tet_poff)) 
  colnames(sex_poff) = gsub(".faa", "", colnames(sex_poff)) 
  
  vsp_bi = colnames(bi_poff)[4:ncol(bi_poff)]                             #get vector of spcies names for each family
  vsp_tet = colnames(tet_poff)[4:ncol(tet_poff)] 
  vsp_sex = colnames(sex_poff)[4:ncol(sex_poff)] 
#----------------------------------------
  lbi_GFF = list()                                                        #make lists for GFFs
  lsex_GFF = list()
  ltet_GFF = list()
  
#---------------------------------------
  
  groupToDo = c("bi_group", "tet_group", "sex_group")[1]                  #do one group at a time. easier for troubleshooting
  
  if( groupToDo=="bi_group"){

    for(sp_bi in vsp_bi){
      bi_gff_name = paste("../Paramecium_GFF/", sp_bi, "-gene.tab", sep="")                                     #name of GFFs all look the same, so use list of names to call file names
      bi_gff_load = read.table(bi_gff_name, sep="\t", as.is=T)                                                  #load each file given its name
      lbi_GFF[[sp_bi]] = bi_gff_load                                                                            #add each file to list of GFFs with name of species as identifier
  }
  
  lbi_FASTA = list()                                                                                            #make list for fastas
  lbi_FASTA[["pbi"]] = read.fasta("../Paramecium_FASTA/biaurelia_V1-4_assembly_v1.fa", as.string=TRUE)          
  lbi_FASTA[["ppent"]] = read.fasta("../Paramecium_FASTA/ppentaurelia_mac_87_v1.0.fa", as.string=TRUE)
  lbi_FASTA[["pprim"]] = read.fasta("../Paramecium_FASTA/primaurelia_Ir4-2_assembly_v1.fasta", as.string=TRUE)
  lbi_FASTA[["pquad"]] = read.fasta("../Paramecium_FASTA/pquadecaurelia_mac_NiA_v1.0.fa", as.string=TRUE)
  lbi_FASTA[["ptre"]] = read.fasta("../Paramecium_FASTA/ptredecaurelia_209_AP38_filtered.fa", as.string=TRUE)
  lbi_FASTA[["pnov"]] = read.fasta("../Paramecium_FASTA/pnovaurelia_mac_TE_v1.0.fa", as.string=TRUE)
  getPOFFid(vsp_bi, bi_poff, lbi_GFF, lbi_FASTA)                                                                #call first function which will call all other and write fasta output 
  
  }
 
if(groupToDo=="sex_group"){

  for(sp_sex in vsp_sex){
    sex_gff_name = paste("../Paramecium_GFF/", sp_sex, "-gene.tab", sep="") 
    sex_gff_load = read.table(sex_gff_name, sep="\t", as.is=T)
    lsex_GFF[[sp_sex]] = sex_gff_load    
  }
  
  lsex_FASTA = list()
  lsex_FASTA[["psex"]] = read.fasta("../Paramecium_FASTA/sexaurelia_AZ8-4_assembly_v1.fasta", as.string=TRUE)
  lsex_FASTA[["pjen"]] = read.fasta("../Paramecium_FASTA/pjenningsi_mac_M_v1.0.fa", as.string=TRUE)
  lsex_FASTA[["pson"]] = read.fasta("../Paramecium_FASTA/psonneborni_mac_ATCC_30995_v1.0.fa", as.string=TRUE)
  getPOFFid(vsp_sex, sex_poff, lsex_GFF, lsex_FASTA)
  
}

if(groupToDo=="tet_group"){  
  
  for(sp_tet in vsp_tet){
    tet_gff_name = paste("../Paramecium_GFF/", sp_tet, "-gene.tab", sep="") 
    tet_gff_load = read.table(tet_gff_name, sep="\t", as.is=T)
    ltet_GFF[[sp_tet]] = tet_gff_load
  }
  
  ltet_FASTA = list()
  ltet_FASTA[["pdec"]] = read.fasta("../Paramecium_FASTA/pdecaurelia_mac_223_v1.0.fa", as.string=TRUE)
  ltet_FASTA[["pdodec"]] = read.fasta("../Paramecium_FASTA/pdodecaurelia_mac_274_v1.0.fa", as.string=TRUE)
  ltet_FASTA[["poct"]] = read.fasta("../Paramecium_FASTA/poctaurelia.fasta", as.string=TRUE)
  ltet_FASTA[["psept"]] =read.fasta("../Paramecium_FASTA/pseptaurelia_mac_38_v1.0.fa", as.string=TRUE)
  ltet_FASTA[["ptet"]] = read.fasta("../Paramecium_FASTA/ptetraurelia_mac_51.fa", as.string=TRUE)
  getPOFFid(vsp_tet, tet_poff, ltet_GFF, ltet_FASTA)
  
}
  