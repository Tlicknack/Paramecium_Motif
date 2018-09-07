#DONE
#This script will take a csv which has motif, evalue, wgd_family and another table attached to it that has file and mismatch per site... can change to include two separate tables
#It will combine them on the basis of the file from which the motif/evalue and mismatch count came

mismatch_table = read.csv("All_Aurelia_12ntMotif_MismatchPerSite.csv", sep=",", header=T)
  
motif_table = read.csv("All_Aurelias_12nt_Regexs_Evalues.csv", sep=",", header=T)

wgd_families = as.character(unique(motif_table$File))
final_df = data.frame(matrix(ncol = 4))
colnames(final_df) = c("Regex", "E_Value", "File", "rep_mismatch")

for(family in wgd_families){
  family_df = motif_table[which(motif_table$File == family),]
  file_df = mismatch_table[which(mismatch_table$File == family),]
  
  if(nrow(file_df) > 0){
    rep_mismatch =  rep(file_df$Mismatch.Site, times=nrow(family_df))
    new_df = as.data.frame(cbind(family_df, rep_mismatch))
    final_df = as.data.frame(rbind(final_df, new_df))
  }
}
final_df = final_df[-1,]
write.csv(final_df, file="All_Aurelias_12nt_Motif_Evalue_Mismatch.csv", sep=",", row.names=F)
