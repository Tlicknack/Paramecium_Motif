


fasta = read.fasta("/home/tlicknac/Desktop/Paramecium_Genome_Data/Paramecium_FASTA/pdecaurelia_mac_223_v1.0.fa", as.string=T, forceDNAtolower=T)      #CHANGE
#vfasta = list.files("/home/tlicknac/Desktop/Paramecium_Genome_Data/Paramecium_GFF", recursive=F)
gff = read.table("/home/tlicknac/Desktop/Paramecium_Genome_Data/Paramecium_GFF/pdec-full.gff", header=F, sep="\t")
#matched_scafs = list()
matched_counter = 1
vseqs = c()
vscafs = c()

for(i in 1:length(fasta)) {
  scaf_i = getName(fasta[[i]])
  for(j in 1:nrow(gff)){
    scaf_j = as.character(gff[j,1])
    
    if(scaf_i == scaf_j){
      vscafs[matched_counter] = scaf_i
      vseqs[matched_counter] = as.character(getSequence(fasta[[i]], as.string=T))
                                        
      #matched_scafs[[matched_counter]]["scaf"] = scaf_i
      #matched_scafs[[matched_counter]]["seq"] =  as.character(getSequence(fasta[[i]], as.string=T))
    }
  }
  matched_counter=matched_counter+1 
}

new_vscafs = vscafs[which(is.na(vscafs)==F)]
new_vseqs = vseqs[which(is.na(vseqs)==F)]

write.fasta(new_vseqs, new_vscafs, file.out="pdec_mod.fa")

#Move fasta to the server if doing it locally
  #Do this with scp /home/user/tlicknac/.....  tlicknac@carbonate.uits.iu.edu:/N/u/Carbonate/....
