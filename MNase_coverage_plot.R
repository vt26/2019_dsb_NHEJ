#Plot coverage trace of mnase time course, wild type samples showen here
#replace the bam file paths with YKU70-, MRE11- or NHEJ data to generate repsective plots

#list of paths to the bam files aligned to the sacCer3 Genome

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

#list of paths to the bam files aligned to the custom ChrII with the HO cut site added

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


#calculate sampling depth for all samples to be at the same read depth

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


#set boundaries of coverage window

ho_start = 431525

ho_end = 431641

start = ho_start-25000

end = ho_end+25000

chr = "chrII"


#intialize an empty list and then calculate a vector of coverage values for each time point and append to this list
wt_mnase_coverage.l = list()

for (i in 1:6){
  
  set.seed(9)
  
  #read in the bam files to data frames
  chr = "chrII"
  df2 = get_dot_mat(as.character(chr2_files2[i]), chr, 1, get_chr_length(as.character(chr2_files2[i]), chr))
  df3 = get_dot_mat(as.character(chr2_files3[i]), chr, 1, get_chr_length(as.character(chr2_files3[i]), chr))
  
  df_2 = rbind(df2, df3)
  
  norm.df = df_2[sample(nrow(df_2), sampling_depth, replace = FALSE),]
  coverage.v = vector()
  for (n in seq(start,(end-500),10)){
    
    mod_start = n
    mod_end = n + 500
    
    occupancy_of_window = which(norm.df$mpoint > mod_start & 
                                  norm.df$mpoint< mod_end)
    
    
    reads_in_window = length(occupancy_of_window)
    
    coverage.v = c(coverage.v, reads_in_window)
    
    
    # create progress bar
    pb <- txtProgressBar(min = 0, max = (end-500-start), style = 3)
    #    Sys.sleep(0.1)
    # update progress bar
    setTxtProgressBar(pb, (n-start))
    
    
  }
  
  wt_mnase_coverage.l[[i]] = coverage.v
  
  print(i)
  
}

close(pb)

#plot the coverage

plot(x=seq(start,(end-500),10),
     y= log2(wt_mnase_coverage.l[[1]]/wt_mnase_coverage.l[[1]]),
     type = "l",
     col = blues9[4],
     lty = 1,
     main = "Relative (Pre-Induction) MNase Coverage",
     ylab = "Log2(FPKM Normalized Ratio)",
     xlab = "Position (ChrII)",
     ylim = c(-3, 1)
)

for (i in 2:6){
  
  lines(x = seq(start,(end-500),10),
        y= log2(wt_mnase_coverage.l[[i]]/wt_mnase_coverage.l[[1]]),
        type = "l",
        col = blues9[i+3],
        lty = 2
  )
  
}

#add a legend

legend(453500, y = -.65, legend = c(0,15,30,60,90,120),
       border = "gray",
       col = blues9[4:9],
       pch = "-",
       lwd = 3,
       cex = 0.4
)

#mark the HO cut site
abline(v = ho_end, col = "red")
abline(v = ho_start, col = "red")


#Plotting Smoothed MNase coverage

plot(x=seq(start,(end-500),10),
     y = smooth.spline(log2(wt_mnase_coverage.l[[6]]/wt_mnase_coverage.l[[1]]), spar = .65)$y,
     type = "l",
     col = "blue",
     lty = 1,
     main = "Smoothed MNase Coverage (120 min)",
     ylab = "Log2(FPKM Normalized Ratio)",
     xlab = "ChrII Position (bp)",
     ylim = c(-2, 1),
     mgp = c(2,1,0),
     cex.axis = 0.8
)
abline(v = ho_end+500, col = "red")
abline(v = ho_start+500, col = "red")

lines(x=seq(start,(end-500),10),
      y = smooth.spline(log2(yku70_mnase_coverage.l[[6]]/yku70_mnase_coverage.l[[1]]), spar = .65)$y,
      type = "l",
      col = "blue",
      lty = 2,
      main = "WT Relative (Pre-Induction) MNase Coverage",
      ylab = "Log2(FPKM Normalized Ratio)",
      xlab = "ChrII Position (bp)",
      ylim = c(-2, 1),
      mgp = c(2,1,0),
      cex.axis = 0.8
)

lines(x=seq(start,(end-500),10),
      y = smooth.spline(log2(mre11_mnase_coverage.l[[6]]/mre11_mnase_coverage.l[[1]]), spar = .65)$y,
      type = "l",
      col = "blue",
      lty = 3,
      main = "Smoothed MNase Coverage (120 min)",
      ylab = "Log2(FPKM Normalized Ratio)",
      xlab = "ChrII Position (bp)",
      ylim = c(-2, 1),
      mgp = c(2,1,0),
      cex.axis = 0.8
)

abline(h = 0, lwd = 0.5)


legend("bottomright", 
       legend = c("WT", 
                  expression(italic("yku70"*Delta)),
                  expression(italic("mre11"*Delta))),
       border = "gray",
       col = "blue",
       lty = c(1,2,3),
       lwd = 1,
       cex = 1
)

dev.off()