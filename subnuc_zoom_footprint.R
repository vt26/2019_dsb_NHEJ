#subnucleosomal zoom in of WT/KU/MRE11

titles = c("Pre-Induction",
           
           "15 Mins",
           "30 Mins",
           "60 Mins",
           "90 Mins",
           "120 Mins"
           
)


sacCer3_files2 = list("/data/home/vt26/DM1129/DM1129_sacCer3_m1_2019-06-19-01-05.bam",
                      "/data/home/vt26/DM1130/DM1130_sacCer3_m1_2019-06-19-01-21.bam",
                      "/data/home/vt26/DM1131/DM1131_sacCer3_m1_2019-06-19-01-36.bam",
                      "/data/home/vt26/DM1132/DM1132_sacCer3_m1_2019-06-19-01-51.bam",
                      "/data/home/vt26/DM1133/DM1133_sacCer3_m1_2019-06-19-02-07.bam",
                      "/data/home/vt26/DM1134/DM1134_sacCer3_m1_2019-06-19-02-20.bam"
)
sacCer3_files3 = list("/data/home/vt26/DM1147/DM1147_sacCer3_m1_2019-07-30-12-12.bam",
                      "/data/home/vt26/DM1148/DM1148_sacCer3_m1_2019-07-30-12-18.bam",
                      "/data/home/vt26/DM1149/DM1149_sacCer3_m1_2019-07-30-12-33.bam",
                      "/data/home/vt26/DM1150/DM1150_sacCer3_m1_2019-07-30-12-40.bam",
                      "/data/home/vt26/DM1151/DM1151_sacCer3_m1_2019-07-30-12-56.bam",
                      "/data/home/vt26/DM1152/DM1152_sacCer3_m1_2019-07-30-13-05.bam"
)



chr2_files2 = list("/data/home/vt26/DM1129/DM1129_chr2_pho5matindel_m1_2019-06-18-22-32.bam",
                   "/data/home/vt26/DM1130/DM1130_chr2_pho5matindel_m1_2019-06-18-22-45.bam",
                   "/data/home/vt26/DM1131/DM1131_chr2_pho5matindel_m1_2019-06-18-22-58.bam",
                   "/data/home/vt26/DM1132/DM1132_chr2_pho5matindel_m1_2019-06-18-23-10.bam",
                   "/data/home/vt26/DM1133/DM1133_chr2_pho5matindel_m1_2019-06-18-23-22.bam",
                   "/data/home/vt26/DM1134/DM1134_chr2_pho5matindel_m1_2019-06-18-23-34.bam"
)
chr2_files3 = list("/data/home/vt26/DM1147/DM1147_chr2_pho5matindel_m1_2019-07-29-15-26.bam",
                   "/data/home/vt26/DM1148/DM1148_chr2_pho5matindel_m1_2019-07-29-15-35.bam",
                   "/data/home/vt26/DM1149/DM1149_chr2_pho5matindel_m1_2019-07-29-15-44.bam",
                   "/data/home/vt26/DM1150/DM1150_chr2_pho5matindel_m1_2019-07-29-15-53.bam",
                   "/data/home/vt26/DM1151/DM1151_chr2_pho5matindel_m1_2019-07-29-16-03.bam",
                   "/data/home/vt26/DM1152/DM1152_chr2_pho5matindel_m1_2019-07-29-16-12.bam"
)

depth.v = vector()

for (i in 1:6){
  
  #read in the bam files to data frames
  chr = "chrII"
  df2 = get_dot_mat(as.character(chr2_files2[i]), chr, 1, get_chr_length(as.character(chr2_files2[i]), chr))
  df3 = get_dot_mat(as.character(chr2_files3[i]), chr, 1, get_chr_length(as.character(chr2_files3[i]), chr))
  
  df_2 = rbind(df2, df3)
  
  depth.v = c(depth.v, dim(df_2)[1])
}

sampling_depth = min(depth.v)



ho_start = 431525

ho_end = 431641

start = ho_start-500

end = ho_end+500

chr = "chrII"


file_name = paste(start, "_", end, "_pho5-wt_zoomed.png", sep="")
png(file_name, width = 2.5, height = 4, units = "in", res = 300)
#par(bg=NA)
par(mfcol=c(3,1))

#this function plots the gene bodies on top of the typhoon plots in gray boxes. if you only want annotated protein coding genes set the proteinCoding = T, otherwise if you set it to F (as i have done here) every ORF will be represented
par(mar=c(1,3,3.5,4))
#MakeArrowSchematic_ho("2", start, end, cex_title = 1.0, proteinCoding = F)

text(x= (431878+431888+117*2)/2,
     y = .25,
     labels = expression(bold("Sum1")),
     srt = 90,
     cex = 1.2
)

text(x= (ho_start+ho_end)/2,
     y = .25,
     labels = expression(bold("HOcs")),
     #las = 2
     srt = 90,
     cex = 1.2
)

par(mar=c(3,3,1,4))

nuc_occupancy.l = list()
nuc_pos.l = list()
nuc_fuzz.l = list()

for (i in seq(1,6,3)){
  
  #plot(density(merged_df$fsize, bw = 1))
  set.seed(9)
  
  #read in the bam files to data frames
  chr = "chrII"
  df2 = get_dot_mat(as.character(chr2_files2[i]), chr, 1, get_chr_length(as.character(chr2_files2[i]), chr))
  df3 = get_dot_mat(as.character(chr2_files3[i]), chr, 1, get_chr_length(as.character(chr2_files3[i]), chr))
  
  df_2 = rbind(df2, df3)
  
  norm.df = df_2[sample(nrow(df_2), sampling_depth, replace = FALSE),]
  
  
  df2 = norm.df
  
  window.df = df2[which(df2$mpoint<(end) & df2$mpoint>start),]
  
  
  dcolor = densColsDM(window.df$mpoint, window.df$fsize,
                      nbin=c(1024,1024), 
                      bandwidth=c(36,16), 
                      transformation = function(x) x^.5,
                      colramp = colorRampPalette(brewer.pal(9, "Oranges")),
                      z_factor = 1
  )
  
 
  
  if (i == 4){
    plot(window.df$mpoint, window.df$fsize, 
         col=dcolor, 
         cex=0.25, pch=20, 
         main=titles[i], 
         xlab='Chr II Position (bp)', 
         ylab='Fragment Size',
         #xaxt = "n",
         mgp = c(2, 1, 0),
         cex.axis = 0.8,
         xaxs = "i"
    )
  } else{
    
    plot(window.df$mpoint, window.df$fsize, 
         col=dcolor, 
         cex=0.25, pch=20, 
         main=titles[i], 
         xlab='', 
         ylab='Fragment Size',
         #xaxt = "n",
         mgp = c(2, 1, 0),
         cex.axis = 0.8,
         xaxs = "i"
    )
    
  }
  
  abline(v = ho_start, col = "blue", lty = "dotted")
  abline(v = ho_end, col = "blue", lty = "dotted")

}



dev.off()


####### YKU70
set.seed(9)

titles = c("Pre-Induction",
           "15 Mins",
           "30 Mins",
           "60 Mins",
           "90 Mins",
           "120 Mins"
           
)


chr2_files2 = list("/data/home/vt26/DM1135/DM1135_chr2_pho5matindel_m1_2019-06-18-23-47.bam",
                   "/data/home/vt26/DM1136/DM1136_chr2_pho5matindel_m1_2019-06-19-00-00.bam",
                   "/data/home/vt26/DM1137/DM1137_chr2_pho5matindel_m1_2019-06-19-00-12.bam",
                   "/data/home/vt26/DM1138/DM1138_chr2_pho5matindel_m1_2019-06-19-00-25.bam",
                   "/data/home/vt26/DM1139/DM1139_chr2_pho5matindel_m1_2019-06-19-00-37.bam",
                   "/data/home/vt26/DM1140/DM1140_chr2_pho5matindel_m1_2019-06-19-00-50.bam"
)
chr2_files3 = list("/data/home/vt26/DM1153/DM1153_chr2_pho5matindel_m1_2019-07-29-16-21.bam",
                   "/data/home/vt26/DM1154/DM1154_chr2_pho5matindel_m1_2019-07-29-16-45.bam",
                   "/data/home/vt26/DM1155/DM1155_chr2_pho5matindel_m1_2019-07-29-16-56.bam",
                   "/data/home/vt26/DM1156/DM1156_chr2_pho5matindel_m1_2019-07-29-17-29.bam",
                   "/data/home/vt26/DM1157/DM1157_chr2_pho5matindel_m1_2019-07-30-01-44.bam",
                   "/data/home/vt26/DM1158/DM1158_chr2_pho5matindel_m1_2019-07-30-01-54.bam"
)


sacCer3_files2 = list("/data/home/vt26/DM1135/DM1135_sacCer3_m1_2019-06-19-02-37.bam",
                      "/data/home/vt26/DM1136/DM1136_sacCer3_m1_2019-06-19-02-51.bam",
                      "/data/home/vt26/DM1137/DM1137_sacCer3_m1_2019-06-19-03-06.bam",
                      "/data/home/vt26/DM1138/DM1138_sacCer3_m1_2019-06-19-03-20.bam",
                      "/data/home/vt26/DM1139/DM1139_sacCer3_m1_2019-06-19-03-35.bam",
                      "/data/home/vt26/DM1140/DM1140_sacCer3_m1_2019-06-19-03-49.bam"
)
sacCer3_files3 = list("/data/home/vt26/DM1153/DM1153_sacCer3_m1_2019-07-30-13-13.bam",
                      "/data/home/vt26/DM1154/DM1154_sacCer3_m1_2019-07-30-13-23.bam",
                      "/data/home/vt26/DM1155/DM1155_sacCer3_m1_2019-07-30-13-35.bam",
                      "/data/home/vt26/DM1156/DM1156_sacCer3_m1_2019-07-30-13-48.bam",
                      "/data/home/vt26/DM1157/DM1157_sacCer3_m1_2019-07-30-13-57.bam",
                      "/data/home/vt26/DM1158/DM1158_sacCer3_m1_2019-07-30-14-07.bam"
)

depth.v = vector()

for (i in 1:6){
  
  #read in the bam files to data frames
  chr = "chrII"
  df2 = get_dot_mat(as.character(chr2_files2[i]), chr, 1, get_chr_length(as.character(chr2_files2[i]), chr))
  df3 = get_dot_mat(as.character(chr2_files3[i]), chr, 1, get_chr_length(as.character(chr2_files3[i]), chr))
  
  df_2 = rbind(df2, df3)
  
  depth.v = c(depth.v, dim(df_2)[1])
}

sampling_depth = min(depth.v)


ho_start = 431525

ho_end = 431641

start = ho_start-500

end = ho_end+500

chr = "chrII"


file_name = paste(start, "_", end, "_pho5-yku70_zoomed.png", sep="")
png(file_name, width = 2.5, height = 4, units = "in", res = 300)
#par(bg=NA)
par(mfcol=c(3,1))

#this function plots the gene bodies on top of the typhoon plots in gray boxes. if you only want annotated protein coding genes set the proteinCoding = T, otherwise if you set it to F (as i have done here) every ORF will be represented
par(mar=c(1,3,3.5,4))
MakeArrowSchematic_ho("2", start, end, cex_title = 1.0, proteinCoding = F)


text(x= (431878+431888+117*2)/2,
     y = .25,
     labels = expression(bold("Sum1")),
     srt = 90,
     cex = 1.2
)

text(x= (ho_start+ho_end)/2,
     y = .25,
     labels = expression(bold("HOcs")),
     #las = 2
     srt = 90,
     cex = 1.2
)

par(mar=c(3,3,1,4))

nuc_occupancy.l = list()
nuc_pos.l = list()
nuc_fuzz.l = list()



for (i in c(1,4)){
  

  set.seed(9)
  
  #read in the bam files to data frames
  chr = "chrII"
  df2 = get_dot_mat(as.character(chr2_files2[i]), chr, 1, get_chr_length(as.character(chr2_files2[i]), chr))
  df3 = get_dot_mat(as.character(chr2_files3[i]), chr, 1, get_chr_length(as.character(chr2_files3[i]), chr))
  
  df_2 = rbind(df2, df3)
  
  norm.df = df_2[sample(nrow(df_2), sampling_depth, replace = FALSE),]
  
  
  df2 = norm.df
  
  window.df = df2[which(df2$mpoint<(end) & df2$mpoint>start),]
  
  
  dcolor = densColsDM(window.df$mpoint, window.df$fsize,
                      nbin=c(1024,1024), 
                      bandwidth=c(36,16), 
                      transformation = function(x) x^.5,
                      colramp = colorRampPalette(brewer.pal(9, "Oranges")),
                      z_factor = 1
  )
  
  
  
  if (i == 4){
    plot(window.df$mpoint, window.df$fsize, 
         col=dcolor, 
         cex=0.25, pch=20, 
         main=titles[i], 
         xlab='Chr II Position (bp)', 
         ylab='Fragment Size',
         #xaxt = "n",
         mgp = c(2, 1, 0),
         cex.axis = 0.8,
         xaxs = "i"
    )
  } else{
    
    plot(window.df$mpoint, window.df$fsize, 
         col=dcolor, 
         cex=0.25, pch=20, 
         main=titles[i], 
         xlab='', 
         ylab='Fragment Size',
         #xaxt = "n",
         mgp = c(2, 1, 0),
         cex.axis = 0.8,
         xaxs = "i"
    )
    
  }
  
  abline(v = ho_start, col = "blue", lty = "dotted")
  abline(v = ho_end, col = "blue", lty = "dotted")
  
  
}



dev.off()






#### MRE11

set.seed(9)

titles = c("Pre-Induction",
           "15 Mins",
           "30 Mins",
           "60 Mins",
           "90 Mins",
           "120 Mins"
           
)


chr2_files2 = list("/data/home/vt26/DM1141/DM1141_chr2_pho5matindel_m1_2019-06-19-01-02.bam",
                   "/data/home/vt26/DM1142/DM1142_chr2_pho5matindel_m1_2019-06-19-01-15.bam",
                   "/data/home/vt26/DM1143/DM1143_chr2_pho5matindel_m1_2019-06-19-01-26.bam",
                   "/data/home/vt26/DM1144/DM1144_chr2_pho5matindel_m1_2019-06-19-01-39.bam",
                   "/data/home/vt26/DM1145/DM1145_chr2_pho5matindel_m1_2019-06-19-01-50.bam",
                   "/data/home/vt26/DM1146/DM1146_chr2_pho5matindel_m1_2019-06-19-02-01.bam"
)
chr2_files3 = list("/data/home/vt26/DM1159/DM1159_chr2_pho5matindel_m1_2019-07-30-02-04.bam",
                   "/data/home/vt26/DM1160/DM1160_chr2_pho5matindel_m1_2019-07-30-02-13.bam",
                   "/data/home/vt26/DM1161/DM1161_chr2_pho5matindel_m1_2019-07-30-02-22.bam",
                   "/data/home/vt26/DM1162/DM1162_chr2_pho5matindel_m1_2019-07-30-02-32.bam",
                   "/data/home/vt26/DM1163/DM1163_chr2_pho5matindel_m1_2019-07-30-02-42.bam",
                   "/data/home/vt26/DM1164/DM1164_chr2_pho5matindel_m1_2019-07-30-02-51.bam"
)

sacCer3_files2 = list("/data/home/vt26/DM1141/DM1141_sacCer3_m1_2019-06-19-04-03.bam",
                      "/data/home/vt26/DM1142/DM1142_sacCer3_m1_2019-06-19-04-18.bam",
                      "/data/home/vt26/DM1143/DM1143_sacCer3_m1_2019-06-19-04-30.bam",
                      "/data/home/vt26/DM1144/DM1144_sacCer3_m1_2019-06-19-04-45.bam",
                      "/data/home/vt26/DM1145/DM1145_sacCer3_m1_2019-06-19-04-57.bam",
                      "/data/home/vt26/DM1146/DM1146_sacCer3_m1_2019-06-19-05-08.bam"
)
sacCer3_files3 = list("/data/home/vt26/DM1159/DM1159_sacCer3_m1_2019-07-30-14-17.bam",
                      "/data/home/vt26/DM1160/DM1160_sacCer3_m1_2019-07-30-14-28.bam",
                      "/data/home/vt26/DM1161/DM1161_sacCer3_m1_2019-07-30-14-42.bam",
                      "/data/home/vt26/DM1162/DM1162_sacCer3_m1_2019-07-30-14-53.bam",
                      "/data/home/vt26/DM1163/DM1163_sacCer3_m1_2019-07-30-15-02.bam",
                      "/data/home/vt26/DM1164/DM1164_sacCer3_m1_2019-07-30-15-15.bam"
)

depth.v = vector()

for (i in 1:6){
  
  #read in the bam files to data frames
  chr = "chrII"
  df2 = get_dot_mat(as.character(chr2_files2[i]), chr, 1, get_chr_length(as.character(chr2_files2[i]), chr))
  df3 = get_dot_mat(as.character(chr2_files3[i]), chr, 1, get_chr_length(as.character(chr2_files3[i]), chr))
  
  df_2 = rbind(df2, df3)
  
  depth.v = c(depth.v, dim(df_2)[1])
}

sampling_depth = min(depth.v)



ho_start = 431525

ho_end = 431641

start = ho_start-500

end = ho_end+500

chr = "chrII"

file_name = paste(start, "_", end, "_pho5-mre11_zoomed.png", sep="")
png(file_name, width = 2.5, height = 4, units = "in", res = 300)
#par(bg=NA)
par(mfcol=c(3,1))

#this function plots the gene bodies on top of the typhoon plots in gray boxes. if you only want annotated protein coding genes set the proteinCoding = T, otherwise if you set it to F (as i have done here) every ORF will be represented
par(mar=c(1,3,3.5,4))
MakeArrowSchematic_ho("2", start, end, cex_title = 1.0, proteinCoding = F)

text(x= (431878+431888+117*2)/2,
     y = .25,
     labels = expression(bold("Sum1")),
     srt = 90,
     cex = 1.2
)

text(x= (ho_start+ho_end)/2,
     y = .25,
     labels = expression(bold("HOcs")),
     #las = 2
     srt = 90,
     cex = 1.2
)

par(mar=c(3,3,1,4))

nuc_occupancy.l = list()
nuc_pos.l = list()
nuc_fuzz.l = list()



for (i in c(1,4)){
  
  set.seed(9)
  
  #read in the bam files to data frames
  chr = "chrII"
  df2 = get_dot_mat(as.character(chr2_files2[i]), chr, 1, get_chr_length(as.character(chr2_files2[i]), chr))
  df3 = get_dot_mat(as.character(chr2_files3[i]), chr, 1, get_chr_length(as.character(chr2_files3[i]), chr))
  
  df_2 = rbind(df2, df3)
  
  norm.df = df_2[sample(nrow(df_2), sampling_depth, replace = FALSE),]
  
  
  df2 = norm.df
  
  window.df = df2[which(df2$mpoint<(end) & df2$mpoint>start),]
  
  
  dcolor = densColsDM(window.df$mpoint, window.df$fsize,
                      nbin=c(1024,1024), 
                      bandwidth=c(36,16), 
                      transformation = function(x) x^.5,
                      colramp = colorRampPalette(brewer.pal(9, "Oranges")),
                      z_factor = 1
  )
  
  
  if (i == 4){
    plot(window.df$mpoint, window.df$fsize, 
         col=dcolor, 
         cex=0.25, pch=20, 
         main=titles[i], 
         xlab='Chr II Position (bp)', 
         ylab='Fragment Size',
         #xaxt = "n",
         mgp = c(2, 1, 0),
         cex.axis = 0.8,
         xaxs = "i"
    )
  } else{
    
    plot(window.df$mpoint, window.df$fsize, 
         col=dcolor, 
         cex=0.25, pch=20, 
         main=titles[i], 
         xlab='', 
         ylab='Fragment Size',
         #xaxt = "n",
         mgp = c(2, 1, 0),
         cex.axis = 0.8,
         xaxs = "i"
    )
    
  }
  
  abline(v = ho_start, col = "blue", lty = "dotted")
  abline(v = ho_end, col = "blue", lty = "dotted")
  
  
  
  
}



dev.off()


