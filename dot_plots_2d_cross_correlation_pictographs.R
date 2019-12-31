#new WT Data

set.seed(9)

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

#this data is derived from the broggard analysis (separate code)
#this is how the idealized nucleosome kernel is generated for each experiment

#WT parameters
std_x = sqrt(2539.909)
nucleosome_size = 159
std_y = sqrt(598.262)

# #YKU70- parameters
# nucleosome_size = 158
# std_x = sqrt(2503.144)
# std_y = sqrt(558.9836)
# 
# #MRE11- parameters
# nucleosome_size = 160
# std_x = sqrt(2581.28)
# std_y = sqrt(625.9666)
# 
# #WT NHEJ Parameters
# nucleosome_size = 166
# std_x = sqrt(2499.701)
# std_y = sqrt(527.0342)
# 
# #DNL4- NHEJ Parameters
# nucleosome_size = 166
# std_x = sqrt(2595.372)
# std_y = sqrt(532.3775)


#find cross correlation score to kernel density across the region of interest

library(reshape)
library(mvtnorm)
library(lattice)

# all lengths and positions
lengths <- seq(1, 250, by = 1)
pos <- seq(-75, 75, by = 1)

# combinatorial combination of lengths and positions
dens.df <- as.data.frame(expand.grid(lengths, pos))
colnames(dens.df) <- c('length', 'pos')



# function to feed into apply
#sigma is variance not std
get.density <- function(row)
{
  length.mean <- 159
  pos.mean <- 0
  length.std <- 598.262/4
  pos.std <- 2539.909/4
  
  covmat <- matrix(c(length.std, 0, 0, pos.std), ncol=2)
  return(dmvnorm(c(row[1], row[2]), mean=c(length.mean, pos.mean), sigma=covmat))
}

# get the density for each pair of length and pos
dens.df$density <- apply(dens.df, 1, get.density)

# pivot narrow length, pos dataframe into a kernel matrix: 
# each row is a length, each column is a position
kernel.df <- cast(dens.df, length ~ pos, value='density')
kernel.mat <- as.matrix(kernel.df)

#this function plots the 2d kernel based on the parameters above

png(filename = "2d_kernel_wt.png", width = 5, height = 10, units = "in",  res = 300)
levelplot(t(kernel.mat),
          main = "Idealized Nucleosome (WT) 2D Kernel",
          xlab = "Pos From Nucleosome Center",
          ylab = "Fragment Size",
          colorkey = list(labels = list(col="white")
          ),
          scales=
            list(
              x= list(at = seq(0,151,25),
                      labels = seq(-75,75,25)),
              y = list(at = seq(0,250,50),
                       labels = seq(0,250,50))
            ),
          col.regions = colorRampPalette(c("white","blue"))(100)
)

dev.off()




###calculate cor dist on chrIV for broggard nucs and calibrate fuzziness to this

#do this for only the first preinduction wt sample
i = 1

#read in the bam files to data frames
chr = "chrIV"
df2 = get_dot_mat(as.character(sacCer3_files2[i]), chr, 1, get_chr_length(as.character(sacCer3_files2[i]), chr))
df3 = get_dot_mat(as.character(sacCer3_files3[i]), chr, 1, get_chr_length(as.character(sacCer3_files3[i]), chr))

df_2 = rbind(df2, df3)

set.seed(9)

norm.df = df_2[sample(nrow(df_2), (sampling_depth/813184)*1531933, replace = FALSE),]


reads.m = get_matrix_from_df(norm.df, 
                             1, get_chr_length(as.character(sacCer3_files2[i]), chr))


score = signal$correlate2d(reads.m,kernel.mat, mode='valid')
chrIV_cor.v = as.vector(score[1,])

#shift the chrIV_cor.v vector by padding 75bp

chrIV_cor.v = c(rep(0,75),chrIV_cor.v)

#normalize to read depth

scaling_factor = (dim(norm.df)[1]/1531933)*2016

chrIV_cor.v = (chrIV_cor.v)*10000/scaling_factor
#cor.v contains all the correlation values for ChrIV
#find the cor values for the nucs on broggard positions
#use wget to get the unique nuc positions from the broggard published data set
#https://www.nature.com/articles/nature11142
#https://media.nature.com/original/nature-assets/nature/journal/v486/n7404/extref/nature11142-s2.txt


broggard_unique.df <- fread("https://media.nature.com/original/nature-assets/nature/journal/v486/n7404/extref/nature11142-s2.txt")

colnames(broggard_unique.df) = c("chrom", "pos", "NCPscore", "NCPscorenoiseratio")

#make a vector with all the x positions of the nucleosomes
chr4_nucs_pos = broggard_unique.df$pos[which(broggard_unique.df$chrom == "chrIV")]

#subset the correlation vector for all of chr4 on the nuc position vector
chr4_nucs_cor.v = chrIV_cor.v[chr4_nucs_pos]

#remove 0 values
chr4_nucs_cor.v = chr4_nucs_cor.v[chr4_nucs_cor.v!=0]

#get the statistics for this distribution
summary(chr4_nucs_cor.v)

#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.0000001 0.0332500 0.0565300 0.0713100 0.0876400 1.1030000

quantile(chr4_nucs_cor.v, c(0.01, 0.025, 0.05, 0.1, .90, .95, .975))
#1%        2.5%          5%         10%         90%         95%       97.5% 
#0.003943144 0.007622462 0.011578938 0.018158418 0.125505977 0.157078384 0.206973987 

#set the threshold value equal to the 10% quantile

threshold_value = 0.018158418


ho_start = 431525

ho_end = 431641

start = ho_start-900

end = ho_end+1000

chr = "chrII"

nucleosome_size = 159

file_name = paste(start, "_", end, "_pho5-wt_duplicates_new_merged_blue-bar_relative-lin_scaled_fuzz.png", sep="")
png(file_name, width = 4, height = 10, units = "in", res = 300)
#par(bg=NA)
par(mfcol=c(7,1))

#this function plots the gene bodies on top of the typhoon plots in gray boxes. if you only want annotated protein coding genes set the proteinCoding = T, otherwise if you set it to F (as i have done here) every ORF will be represented
par(mar=c(1,3,3.5,4))
MakeArrowSchematic_ho("2", start, end, cex_title = 1.5, proteinCoding = F)

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


for (i in 1:6){
  
  set.seed(9)
  
  #read in the bam files to data frames
  chr = "chrII"
  df2 = get_dot_mat(as.character(chr2_files2[i]), chr, 1, get_chr_length(as.character(chr2_files2[i]), chr))
  df3 = get_dot_mat(as.character(chr2_files3[i]), chr, 1, get_chr_length(as.character(chr2_files3[i]), chr))
  
  df_2 = rbind(df2, df3)
  
  norm.df = df_2[sample(nrow(df_2), sampling_depth, replace = FALSE),]
  
  reads.m = get_matrix_from_df(norm.df, 
                               (start-76), (end+151))
  
  cor.v = vector()
  
  for (n in 1:((end)-(start-76))){
    
    mod_start = n
    mod_end = n + 150
    
    temp_window = reads.m[1:250,mod_start:mod_end]
    
    score = signal$correlate2d(temp_window,kernel.mat, mode='valid')
    cor.v = c(cor.v, score)
    
  }
  #subet on the window we want to plot
  cor.v = cor.v[1:((end)-start)]
  
  df2 = norm.df
  
  window.df = df2[which(df2$mpoint<(end) & df2$mpoint>start),]
  
  cor.v = cor.v*10000/dim(window.df)[1]
  
  
  dcolor = densColsDM(window.df$mpoint, window.df$fsize,
                      nbin=c(1024,1024), 
                      bandwidth=c(36,16), 
                      transformation = function(x) x^.5,
                      colramp = colorRampPalette(brewer.pal(9, "Oranges")),
                      z_factor = 1
  )
  
  
  if (i == 6){
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
  
  
  peaks.df = data.frame(findpeaks(smooth.spline(cor.v, spar = .001)$y, 
                                  minpeakdistance = 75, 
                                  minpeakheight = threshold_value, 
                                  sortstr= TRUE,
                                  nups= 1))
 
  colnames(peaks.df) = c("height", "center", "start", "end")
  
  peaks.v = sort(peaks.df$center)
  
  
  
  nuc_pos.l[[i]] = peaks.v+start
  nuc_fuzz.l[[i]] = cor.v[peaks.v]
  nuc_pos.v = vector()
  
  y = nucleosome_size
  x = peaks.v+start
  
 
  
  for (p in 1:length(x)){
    nuc_window.df = df2[which(df2$mpoint<=(x[p]+std_x) &
                                df2$mpoint>=(x[p]-std_x) &
                                df2$fsize>=(y-std_y) &
                                df2$fsize<=(y+std_y)),]
    
    
    nuc_pos.v = c(nuc_pos.v, dim(nuc_window.df)[1])
    
    rect(xleft = x[p]-std_x,
         ybottom = 0,
         ytop = 250,
         xright = x[p]+std_x,
         border = NA,
         col = rgb(0,0,1.0,alpha=.15*nuc_fuzz.l[[i]][p]/max(unlist(nuc_fuzz.l))))
    
  }
  
  nuc1L = nuc_pos.l[[1]][5]
  
  nuc3R = nuc_pos.l[[1]][8]
  
  nfr = (nuc_pos.l[[1]][4] + nuc_pos.l[[1]][5])/2
  
  
  points(x= nuc1L,
         y = 235,
         pch = 25,
         cex = 1,
         col = "deepskyblue3",
         bg="deepskyblue3")
  
  nuc_window.df = df2[which(df2$mpoint<=(x[nuc1L]+std_x) &
                              df2$mpoint>=(x[nuc1L]-std_x) &
                              df2$fsize>=(y-std_y) &
                              df2$fsize<=(y+std_y)),]
  
  #comment or uncomment this line below to get the occupancy (read) counts for each nucleosome printed on the graph
  #this is what is plotted for nucleosome kinetics in the line plots
  text(x[nuc1L], y = 35, labels = dim(nuc_window.df)[1], cex = 0.75)
  
  points(x= nuc3R,
         y = 235,
         pch = 25,
         cex = 1,
         col = "darkgoldenrod3",
         bg="darkgoldenrod3")
  
  nuc_window.df = df2[which(df2$mpoint<=(x[nuc3R]+std_x) &
                              df2$mpoint>=(x[nuc3R]-std_x) &
                              df2$fsize>=(y-std_y) &
                              df2$fsize<=(y+std_y)),]
  #comment or uncomment this line below to get the occupancy (read) counts for each nucleosome printed on the graph
  #this is what is plotted for nucleosome kinetics in the line plots
  text(x[nuc3R], y = 35, labels = dim(nuc_window.df)[1], cex = 0.75)
  
  points(x= nfr,
         y = 235,
         pch = 25,
         cex = 1,
         col = "red",
         bg="red")
  
  nuc_window.df = df2[which(df2$mpoint<=(x[nfr]+std_x) &
                              df2$mpoint>=(x[nfr]-std_x) &
                              df2$fsize>=(y-std_y) &
                              df2$fsize<=(y+std_y)),]
  #comment or uncomment this line below to get the occupancy (read) counts for each nucleosome printed on the graph
  #this is what is plotted for nucleosome kinetics in the line plots
  text(x[nfr], y = 35, labels = dim(nuc_window.df)[1], cex = 0.75)
  
  
  
  nuc_occupancy.l[[i]] = nuc_pos.v
  #abline(v = peaks.v, col = "deeppink")
  
  par(new=TRUE)
  
  plot(cor.v, 
       type = "l",
       ylim = c(0, 0.2),
       axes = FALSE,
       ann=FALSE,
       col = "darkslategray",
       xaxs = "i"
  )
  
  
  print(titles[i])
  
}



dev.off()




file_name = paste(start, "_", end, "_pho5-wt_new_nuc_pos_time.png", sep="")
png(file_name, width = 6, height = 6, units = "in", res = 300)
#par(bg=NA)
par(mfcol=c(2,1))

#this function plots the gene bodies on top of the typhoon plots in gray boxes. if you only want annotated protein coding genes set the proteinCoding = T, otherwise if you set it to F (as i have done here) every ORF will be represented
par(mar=c(1,3.5,8,1))
MakeArrowSchematic_ho("2", start, end, cex_title = 1.5, proteinCoding = F)
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

par(mar=c(3.5,3.5,1,1))

i = 1

heights = c(13, 11.5, 10, 7, 4, 1)

x <- nuc_pos.l[[i]]

plot(x = x,
     y = rep(heights[i], length(x)),
     pch = 19, 
     ylab = "Time Post Induction (Mins)", 
     xlab = "Chr II Position (bp)", 
     cex = 2*sqrt(nuc_occupancy.l[[i]]/max(unlist(nuc_occupancy.l))),
     col = rgb(0,0,1, alpha = nuc_fuzz.l[[i]]/max(unlist(nuc_fuzz.l))),
     ylim = c(0,14),
     yaxt='n',
     main = "WT Pho5 Nucleosome Positioning",
     mgp = c(2.5, 1, 0),
     cex.axis = 0.8,
     xlim = c(start,end),
     xaxs = "i"
)

abline(v = ho_start, col = "blue", lty = "dotted")
abline(v = ho_end, col = "blue", lty = "dotted")

axis(2, at = c(13, 11.5, 10, 7, 4, 1), 
     labels = c("Pre", "15", "30", "60", "90", "120"),
     las=1,
     cex = 0.8)

for (i in 2:6){
  
  
  x <- nuc_pos.l[[i]]
  
  points(x = x,
         y = rep(heights[i], length(x)),
         pch = 19, 
         ylab = "Time Post Induction", 
         xlab = "Position (ChrII)", 
         cex = 2*sqrt(nuc_occupancy.l[[i]]/max(unlist(nuc_occupancy.l))),
         col = rgb(0,0,1, alpha = nuc_fuzz.l[[i]]/max(unlist(nuc_fuzz.l))),
         ylim = c(0,14),
         yaxt='n',
         main = "WT Pho5 Nucleosome Positioning",
         mgp = c(2, 1, 0),
         cex.axis = 0.8
  )
  
}

dev.off()

