sacCer3_files = list("/data/home/vt26/DM1129/DM1129_sacCer3_m1_2019-06-19-01-05.bam",
                     "/data/home/vt26/DM1130/DM1130_sacCer3_m1_2019-06-19-01-21.bam",
                     "/data/home/vt26/DM1131/DM1131_sacCer3_m1_2019-06-19-01-36.bam",
                     "/data/home/vt26/DM1132/DM1132_sacCer3_m1_2019-06-19-01-51.bam",
                     "/data/home/vt26/DM1133/DM1133_sacCer3_m1_2019-06-19-02-07.bam",
                     "/data/home/vt26/DM1134/DM1134_sacCer3_m1_2019-06-19-02-20.bam",
                     "/data/home/vt26/DM1147/DM1147_sacCer3_m1_2019-07-30-12-12.bam",
                     "/data/home/vt26/DM1148/DM1148_sacCer3_m1_2019-07-30-12-18.bam",
                     "/data/home/vt26/DM1149/DM1149_sacCer3_m1_2019-07-30-12-33.bam",
                     "/data/home/vt26/DM1150/DM1150_sacCer3_m1_2019-07-30-12-40.bam",
                     "/data/home/vt26/DM1151/DM1151_sacCer3_m1_2019-07-30-12-56.bam",
                     "/data/home/vt26/DM1152/DM1152_sacCer3_m1_2019-07-30-13-05.bam"
)

chr2_files = list("/data/home/vt26/DM1129/DM1129_chr2_pho5matindel_m1_2019-06-18-22-32.bam",
                  "/data/home/vt26/DM1130/DM1130_chr2_pho5matindel_m1_2019-06-18-22-45.bam",
                  "/data/home/vt26/DM1131/DM1131_chr2_pho5matindel_m1_2019-06-18-22-58.bam",
                  "/data/home/vt26/DM1132/DM1132_chr2_pho5matindel_m1_2019-06-18-23-10.bam",
                  "/data/home/vt26/DM1133/DM1133_chr2_pho5matindel_m1_2019-06-18-23-22.bam",
                  "/data/home/vt26/DM1134/DM1134_chr2_pho5matindel_m1_2019-06-18-23-34.bam",
                  "/data/home/vt26/DM1147/DM1147_chr2_pho5matindel_m1_2019-07-29-15-26.bam",
                  "/data/home/vt26/DM1148/DM1148_chr2_pho5matindel_m1_2019-07-29-15-35.bam",
                  "/data/home/vt26/DM1149/DM1149_chr2_pho5matindel_m1_2019-07-29-15-44.bam",
                  "/data/home/vt26/DM1150/DM1150_chr2_pho5matindel_m1_2019-07-29-15-53.bam",
                  "/data/home/vt26/DM1151/DM1151_chr2_pho5matindel_m1_2019-07-29-16-03.bam",
                  "/data/home/vt26/DM1152/DM1152_chr2_pho5matindel_m1_2019-07-29-16-12.bam"
)

std_x = 50.39751

nucleosome_size = 159

#var(nuc_y_pos) 598.262

std_y = 24.45939

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

#levelplot(t(kernel.mat))



chroms.v = c("chrI",
             "chrII",
             "chrIII",
             "chrIV",
             "chrV",
             "chrVI",
             "chrVII",
             "chrVIII",
             "chrIX",
             "chrX",
             "chrXI",
             "chrXII",
             "chrXIII",
             "chrXIV",
             "chrXV",
             "chrXVI"
)

nuc_pos1.l = list()
nuc_occupancy1.l = list()
nuc_fuzz1.l = list()

nuc_pos6.l = list()
nuc_occupancy6.l = list()
nuc_fuzz6.l = list()

for (c in 2:2){
  if (c==2){
    files = chr2_files
  }else{
    files = sacCer3_files
  }
  
  depth.v = vector()
  for (i in seq(1,6,5)){
    
    #read in the bam files to data frames
    chr = chroms.v[c]
    df1 = get_dot_mat(as.character(files[i]), chr, 1, get_chr_length(as.character(files[i]), chr))
    df2 = get_dot_mat(as.character(files[i+6]), chr, 1, get_chr_length(as.character(files[i+6]), chr))
    
    merged.df = rbind(df1, df2)
    
    depth.v = c(depth.v, dim(merged.df)[1])
  }
  
  sampling_depth = min(depth.v)
  
  set.seed(9)
  
  chr = chroms.v[c]
  df1 = get_dot_mat(as.character(files[1]), chr, 1, get_chr_length(as.character(files[1]), chr))
  df2 = get_dot_mat(as.character(files[7]), chr, 1, get_chr_length(as.character(files[7]), chr))
  
  merged.df = rbind(df1, df2)
  
  
  norm.df = merged.df[sample(nrow(merged.df), sampling_depth, replace = FALSE),]
  
  reads.m = get_matrix_from_df(norm.df, 
                               1, get_chr_length(as.character(files[1]), chr))
  
  
  score = signal$correlate2d(reads.m,kernel.mat, mode='valid')
  chr_cor.v = as.vector(score[1,])
  
  #shift the chrIV_cor.v vector by padding 75bp
  
  chr_cor.v = c(rep(0,75),chr_cor.v)
  
  peaks.df = data.frame(findpeaks(chr_cor.v, 
                                  minpeakdistance = 75, 
                                  minpeakheight = 0.011578938, 
                                  sortstr= TRUE,
                                  nups= 1))
  colnames(peaks.df) = c("height", "center", "start", "end")
  
  peaks.v = sort(peaks.df$center)
  
  #abline(v =peaks.v)
  
  nuc_pos1.l[[c]] = peaks.v
  nuc_fuzz1.l[[c]] = chr_cor.v[peaks.v]
  nuc_pos.v = vector()
  
  y = nucleosome_size
  x = nuc_pos1.l[[c]]
  
  for (p in 1:length(x)){
    nuc_window.df = df2[which(norm.df$mpoint<=(x[p]+std_x) &
                                norm.df$mpoint>=(x[p]-std_x) &
                                norm.df$fsize>=(y-std_y) &
                                norm.df$fsize<=(y+std_y)),]
    
    
    nuc_pos.v = c(nuc_pos.v, dim(nuc_window.df)[1])
    
    
  }
  
  
  
  nuc_occupancy1.l[[c]] = nuc_pos.v
  
  
  chr = chroms.v[c]
  df1 = get_dot_mat(as.character(files[6]), chr, 1, get_chr_length(as.character(files[6]), chr))
  df2 = get_dot_mat(as.character(files[12]), chr, 1, get_chr_length(as.character(files[12]), chr))
  
  merged.df = rbind(df1, df2)
  
  
  norm.df = merged.df[sample(nrow(merged.df), sampling_depth, replace = FALSE),]
  
  reads.m = get_matrix_from_df(norm.df, 
                               1, get_chr_length(as.character(files[6]), chr))
  
  
  score = signal$correlate2d(reads.m,kernel.mat, mode='valid')
  chr_cor.v = as.vector(score[1,])
  
  #shift the chrIV_cor.v vector by padding 75bp
  
  chr_cor.v = c(rep(0,75),chr_cor.v)
  
  
  peaks.df = data.frame(findpeaks(chr_cor.v,
                                  minpeakdistance = 75,
                                  minpeakheight = 0.011578938,
                                  sortstr= TRUE,
                                  nups= 1))
  
  colnames(peaks.df) = c("height", "center", "start", "end")
  
  peaks.v = sort(peaks.df$center)
  
  #abline(v =peaks.v)
  
  nuc_pos6.l[[c]] = peaks.v
  
  nuc_fuzz6.l[[c]] = chr_cor.v[peaks.v]
  nuc_pos.v = vector()
  
  for (p in 1:length(x)){
    nuc_window.df = df2[which(norm.df$mpoint<=(x[p]+std_x) &
                                norm.df$mpoint>=(x[p]-std_x) &
                                norm.df$fsize>=(y-std_y) &
                                norm.df$fsize<=(y+std_y)),]
    
    
    nuc_pos.v = c(nuc_pos.v, dim(nuc_window.df)[1])
    
    
  }
  
  
  
  nuc_occupancy6.l[[c]] = nuc_pos.v
  
  
  total <- 16
  # create progress bar
  pb <- txtProgressBar(min = 1, max = total, style = 3)
  #    Sys.sleep(0.1)
  # update progress bar
  setTxtProgressBar(pb, c)
  
}

close(pb)



#calculate min distance from nucs in t120 to nucs in t0


c = 2

#for (c in 1:16){
nuc_pos.l = list()
nuc_dist.l = list()
if (c==2){
  files = chr2_files
}else{
  files = sacCer3_files
}

depth.v = vector()
for (i in seq(1,6,5)){
  
  #read in the bam files to data frames
  chr = chroms.v[c]
  df1 = get_dot_mat(as.character(files[i]), chr, 1, get_chr_length(as.character(files[i]), chr))
  df2 = get_dot_mat(as.character(files[i+6]), chr, 1, get_chr_length(as.character(files[i+6]), chr))
  
  merged.df = rbind(df1, df2)
  
  depth.v = c(depth.v, dim(merged.df)[1])
}

sampling_depth = min(depth.v)

chr = chroms.v[c]
df1 = get_dot_mat(as.character(files[1]), chr, 1, get_chr_length(as.character(files[1]), chr))
df2 = get_dot_mat(as.character(files[7]), chr, 1, get_chr_length(as.character(files[7]), chr))

merged.df = rbind(df1, df2)


norm.df = merged.df[sample(nrow(merged.df), sampling_depth, replace = FALSE),]

reads.m = get_matrix_from_df(norm.df, 
                             1, get_chr_length(as.character(files[1]), chr))


score = signal$correlate2d(reads.m,kernel.mat, mode='valid')
chr_cor.v = as.vector(score[1,])

#shift the chrIV_cor.v vector by padding 75bp

chr_cor.v = c(rep(0,75),chr_cor.v)



peaks.df = data.frame(findpeaks(chr_cor.v, 
                                minpeakdistance = 75, 
                                minpeakheight = 0.011578938, 
                                sortstr= TRUE,
                                nups= 1))
colnames(peaks.df) = c("height", "center", "start", "end")


peaks.v = sort(peaks.df$center)

#abline(v =peaks.v)

nuc_pos.l[[1]] = peaks.v
nuc_dist.l[[1]] = rep(0, length(peaks.v))


set.seed(9)
#  for (i in 2:6){
i = 6   
chr = chroms.v[c]
df1 = get_dot_mat(as.character(files[i]), chr, 1, get_chr_length(as.character(files[i]), chr))
df2 = get_dot_mat(as.character(files[i+6]), chr, 1, get_chr_length(as.character(files[i+6]), chr))

merged.df = rbind(df1, df2)


norm.df = merged.df[sample(nrow(merged.df), sampling_depth, replace = FALSE),]

reads.m = get_matrix_from_df(norm.df, 
                             1, get_chr_length(as.character(files[i]), chr))


score = signal$correlate2d(reads.m,kernel.mat, mode='valid')
chr_cor.v = as.vector(score[1,])

#shift the chrIV_cor.v vector by padding 75bp

chr_cor.v = c(rep(0,75),chr_cor.v)


peaks.df = data.frame(findpeaks(chr_cor.v,
                                minpeakdistance = 75,
                                minpeakheight = 0.011578938,
                                sortstr= TRUE,
                                nups= 1))

colnames(peaks.df) = c("height", "center", "start", "end")


peaks.v = sort(peaks.df$center)

#  abline(v = peaks.v)

nuc_pos.l[[i]] = peaks.v

#abline(v =peaks.v)

nuc_dist.v = vector()
for (n in 1:length(nuc_pos1.l[[2]])){
  min_nuc_dist = min(abs(nuc_pos1.l[[2]][n]-nuc_pos6.l[[2]]))
  
  nuc_dist.v = c(nuc_dist.v, min_nuc_dist)
}


nuc_dist.l[[i]] = nuc_dist.v


file_name = paste("WT_chrII_nuc_occupancy_pos_0vsAll_120_5p_threshold.png", sep="")
png(file_name, width = 8, height = 4, units = "in", res = 300)
#par(bg=NA)
par(mar=c(3,3,2,3))
options(scipen = 999)
plot(x=nuc_pos1.l[[2]],
     # y=log2(nuc_dist.l[[i]]+1),
     y=nuc_dist.v,
     pch = 1,
     cex = .5,
     lwd = 2,
     col= "blue",
     ylim = c(0,500),
     xlim = c((ho_start-25000), (ho_end+25000)),
     #ylim = c(0,max(log2(unlist(nuc_dist.l))+1)),
     xlab = "Chr II Position (bp)",
     ylab = "Min Distance to Pre-Ind N.N." ,
     main = "Nuclosome Positioning vs. Occupancy",
     mgp = c(2,1,0),
     cex.axis = 0.8
)

par(new = TRUE)

plot(x=nuc_pos1.l[[2]],
     y= log2(nuc_occupancy6.l[[2]]/nuc_occupancy1.l[[2]]),
     pch = 2,
     cex = .5,
     lwd = 2,
     col= "red",
     ylim = c(-3,2),
     xlim = c((ho_start-25000), (ho_end+25000)),
     xlab = "",
     ylab = "",
     axes = F,
     #xlab = "Chr II Position (bp)",
     #ylab = "log2(Relative Nuc. Occupancy)",
     #main = "Change in Nuclosome Occupancy",
     cex.axis = 0.8
)

axis(4,
     at = seq(-3,2,1),
     labels = T,
     cex.axis = 0.8
)
mtext("log2(Relative Nuc. Occupancy)", side = 4, line = 2)

abline(v = c((ho_start-1000), (ho_end+1000)),
       col = "black",
       lty = "dotted",
       lwd = 1.5)

abline(v = c((ho_start), (ho_end)),
       col = "red",
       lty = "solid",
       lwd = 1)

abline(v= (ho_start-8000),
       col = "green",
       lty = "dotted",
       lwd = 1.5)
abline(v= (ho_end+8000),
       col = "green",
       lty = "dotted",
       lwd = 1.5)

legend("topright",
       legend = c("Position", "Occupancy"),
       pch = c(1,2),
       lwd = 2,
       lty = 0,
       col = c("blue", "red"),
       cex = .8)



dev.off()

