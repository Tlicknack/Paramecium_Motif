extractUpstream = function(sp, geneID, lGFF, lFasta){                          #make gff and fasta into list with id as sp

  
}



getPOFFid = function(vsp, POFF) {                                     #Pass POFF and list of species names

  i = 1                                                               #counter for iterating through poff
  for(i in 1:nrow(POFF)) {                                            #iterate through POFF
    cat("Parsing line ", i, "\n") 
    for(sp in vsp) {                                                  #iterate through species names, which are identifiers
      sp_family = as.character(POFF[i,sp])                            #get proteinIDs from each column in each row
      vsplit = strsplit(sp_family, ",")[[1]]                          #split each proteinID as it comes out, some will have commas, some won't. Will do this for each species as we iterate
      
      for(geneID in vsplit){                                          #iterate through
        cat("Now extracting sequence upstream of ", geneID, " (which is in", sp, ")\n")
        extractUpstream(sp, geneID, lGFF, lFasta)                     #call function that extracts sequence given species name, the proteinID, a list of GFFs, and a list of FASTAs
      }
    }
      cat("line", i, "done.\n-----------------------\n")
  }
}


writeVectorAsFasta <- function(vseq, fname){                          #Pass list of sequences from extractUpstream and name of file. This will write to fasta file as output
  cat("", file=fname) # Just to reset the file to empty if it already exists
  for(id in names(vseq)){
    cat(">", id, "\n", seq, "\n", sep="", file=fname, append=T)       #dirty fasta format. will write to current directory
  }
}



main = function(){
  
  bi_poff = read.table("~/Paramecium_POFF/biaurelia_group.poff", header=TRUE, sep="\t", comment.char="$")  #read POFFs, change comment character so first line is read as header
  tet_poff = read.table("~/Paramecium_POFF/tetraurelia_group.poff", header=TRUE, sep="\t", comment.char="$")
  sex_poff = read.table("~/Paramecium_POFF/sexaurelia_group.poff", header=TRUE, sep="\t", comment.char="$")
  
  for(poff in POFF){                                        #get list of species names
    colnames(poff) = gsub(".faa", "", colnames(poff))                   #modify column names to remove .faa .... CHANGE FOR DIFFERENT COLUMN NAMES
    vsp = colnames(poff)[4:ncol(poff)]                                  #capture proteinIDs in each line using column name as identifier
  }
  #----------------------------------------
  
  lGFF = list()                                                         #Make list of GFFs. 
  for(sp in vsp){                                     
    gff_name = paste("~/", sp, "-gene.tab", sep="")             #name of GFFs all look the same, so use list of names to create file names
    gff_load = read.table(gff_name, sep="\t")                   #load each file given its name
    lGFF[sp] = gff_load                                         #add each file to list of GFFs with name of species as identifier
  }
  bi_GFF = lGFF["pbi", "ppent", "pprim", "pquad", "ptre", "pnov"]
  sex_GFF = lGFF["psex", "pjen", "pson"]
  tet_GFF = lGFF["ptet", "psept", "poct", "pdodec", "pdec"]
  #---------------------------------------
  
  lFASTA = list()                                                             #Make list of FASTA. Must manually plug in each due to different file names
  lFASTA["pbi"] = "biaurelia_V1-4_assembly_v1.fa"
  lFASTA["pdec"] = "pdecaurelia_mac_223_v1.0.fa" 
  lFASTA["dodec"] = "pdodecaurelia_mac_274_v1.0.fa" 
  lFASTA["pjen"] = "pjenningsi_mac_M_v1.0.fa"
  lFASTA["pnov"] = "pnovaurelia_mac_TE_v1.0.fa"
  lFASTA["poct"] = "poctaurelia.fasta"
  lFASTA["ppent"] = "ppentaurelia_mac_87_v1.0.fa"   
  lFASTA["pprim"] = "primaurelia_Ir4-2_assembly_v1.fasta"
  lFASTA["pquad"] = "pquadecaurelia_mac_NiA_v1.0.fa"
  lFASTA["psept"] = "pseptaurelia_mac_38_v1.0.fa"
  lFASTA["pson"] = "psonneborni_mac_ATCC_30995_v1.0.fa"
  lFASTA["ptre"] = "ptredecaurelia_209_AP38_filtered.fa"
  lFASTA["ptet"] = "ptetraurelia_mac_51.fa"
  
  bi_FASTA = list()
  for(bi_f in lFASTA["pbi", "ppent", "pprim", "pquad", "ptre", "pnov"]){
    for(sp in vsp){
      bi_FASTA[sp] = bi_f
    }
  }
  sex_FASTA = list()
  for(sex_f in lFASTA["psex", "pjen", "pson"]){
    sex_FASTA = append(sex_FASTA, sex_g)
  }
  tet_FASTA = list()
  for(tet_f in lFASTA["ptet", "psept", "poct", "pdodec", "pdec"]){
    tet_FASTA = append(tet_FASTA, tet_g)
  }
}  
#---------------------------------------------------------------------------------------------------------------------------------------  
 














  
pbi_gff = read.table("~/Paramecium_GFF/pbi-gene.tab", sep="\t")
pdec_gff = read.table("pdec-gene.tab", sep="\t")
pdodec_gff = read.table("pdodec-gene.tab", sep="\t")
pjen_gff = read.table("pjen-gene.tab", sep="\t")
pnov_gff = read.table("pnov-gene.tab", sep="\t")
poct_gff = read.table("poct-gene.tab", sep="\t")
ppent_gff = read.table("ppent-gene.tab", sep="\t")
pprim_gff = read.table("pprim-gene.tab", sep="\t") 
pquad_gff = read.table("pquad-gene.tab", sep="\t")
psept_gff = read.table("psept-gene.tab", sep="\t")
psex_gff = read.table("psex-gene.tab", sep="\t")
pson_gff = read.table("pson-gene.tab", sep="\t")
ptet_gff = read.table("ptet-gene.tab", sep="\t")
ptre_gff = read.table("ptre-gene.tab", sep="\t")

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


pbi_gff = "~/Paramecium_GFF/pbi-gene.tab"
pdec_gff = "~/pdec-gene.tab"
pdodec_gff = "~/pdodec-gene.tab "
pjen_gff = "~/pjen-gene.tab "
pnov_gff = "~/pnov-gene.tab "
poct_gff = "~/poct-gene.tab "
ppent_gff = "~/ppent-gene.tab "
pprim_gff = "~/pprim-gene.tab "
pquad_gff = "~/pquad-gene.tab "
psept_gff = "~/psept-gene.tab "
psex_gff = "~/psex-gene.tab "
pson_gff = "~/pson-gene.tab "
ptet_gff = "~/ptet-gene.tab "
ptre_gff = "~/ptre-gene.tab "
lGFF = c(pbi_gff, pdec_gff, pdodec_gff, pjen_gff, pnov_gff, poct_gff, ppent_gff, pprim_gff, pquad_gff, psept_gff, psex_gff, pson_gff, ptet_gff, ptre_gff)
names(lGFF) = c("pbi", "pdec", "pdodec", "pjen", "pnov", "poct", "ppent", "pprim", "pquad", "psept", "psex", "pson", "ptet", "ptre")







