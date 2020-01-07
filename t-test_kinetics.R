#t test for nuc kinetics

#WT Data

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
  
  #depth.v = c(depth.v, dim(df2)[1], dim(df3)[1])
  df_2 = rbind(df2, df3)
  
  depth.v = c(depth.v, dim(df_2)[1])
}

sampling_depth = round(min(depth.v)/2)


#these kernel parameters are derived in the idealized nuc profiling and used in the dot plotting
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

#nuc positions and the linker (nfr) are derived in the dot plot code from the pre-induction sample

nfr = 431325
nuc1L = 431444
nuc3R = 432141

nuc1L_rep1.v = vector()
nuc1L_rep2.v = vector()

nuc3R_rep1.v = vector()
nuc3R_rep2.v = vector()

nfr_rep1.v = vector()
nfr_rep2.v = vector()

for (i in 1:6){
  
  #plot(density(merged_df$fsize, bw = 1))
  set.seed(9)
  
  #read in the bam files to data frames
  chr = "chrII"
  df2 = get_dot_mat(as.character(chr2_files2[i]), chr, 1, get_chr_length(as.character(chr2_files2[i]), chr))
  df3 = get_dot_mat(as.character(chr2_files3[i]), chr, 1, get_chr_length(as.character(chr2_files3[i]), chr))
  
  norm.df2 = df2[sample(nrow(df2), sampling_depth, replace = T),]
  norm.df3 = df3[sample(nrow(df3), sampling_depth, replace = T),]
  
  y = nucleosome_size

  
#nuc1L
nuc1L_rep1.v = c(nuc1L_rep1.v,
  length(
    which(norm.df2$mpoint<=(nuc1L+std_x) &
            norm.df2$mpoint>=(nuc1L-std_x) &
            norm.df2$fsize>=(y-std_y) &
            norm.df2$fsize<=(y+std_y))
  )
)

nuc1L_rep2.v = c(nuc1L_rep2.v,
                 length(
                   which(norm.df3$mpoint<=(nuc1L+std_x) &
                           norm.df3$mpoint>=(nuc1L-std_x) &
                           norm.df3$fsize>=(y-std_y) &
                           norm.df3$fsize<=(y+std_y))
                 )
)


#nuc 3R

nuc3R_rep1.v = c(nuc3R_rep1.v,
                 length(
                   which(norm.df2$mpoint<=(nuc3R+std_x) &
                           norm.df2$mpoint>=(nuc3R-std_x) &
                           norm.df2$fsize>=(y-std_y) &
                           norm.df2$fsize<=(y+std_y))
                 )
)

nuc3R_rep2.v = c(nuc3R_rep2.v,
                 length(
                   which(norm.df3$mpoint<=(nuc3R+std_x) &
                           norm.df3$mpoint>=(nuc3R-std_x) &
                           norm.df3$fsize>=(y-std_y) &
                           norm.df3$fsize<=(y+std_y))
                 )
)

#linker
nfr_rep1.v = c(nfr_rep1.v,
                 length(
                   which(norm.df2$mpoint<=(nfr+std_x) &
                           norm.df2$mpoint>=(nfr-std_x) &
                           norm.df2$fsize>=(y-std_y) &
                           norm.df2$fsize<=(y+std_y))
                 )
)

nfr_rep2.v = c(nfr_rep2.v,
                 length(
                   which(norm.df3$mpoint<=(nfr+std_x) &
                           norm.df3$mpoint>=(nfr-std_x) &
                           norm.df3$fsize>=(y-std_y) &
                           norm.df3$fsize<=(y+std_y))
                 )
)
    
    
    
  
  print(titles[i])
  
}



#MRE11


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
  
  depth.v = c(depth.v, dim(df2)[1], dim(df3)[1])
}

sampling_depth = min(depth.v)

#these kernel parameters are derived in the idealized nuc profiling and used in the dot plotting

std_x = 50.8063

nucleosome_size = 160

#var(nuc_y_pos) 625.9666
#std(nuc_y_pos) 25.01932

std_y = 25.01932

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
  length.mean <- 158
  pos.mean <- 0
  length.std <- 625.9666/4
  pos.std <- 2581.28/4
  
  covmat <- matrix(c(length.std, 0, 0, pos.std), ncol=2)
  return(dmvnorm(c(row[1], row[2]), mean=c(length.mean, pos.mean), sigma=covmat))
}

# get the density for each pair of length and pos
dens.df$density <- apply(dens.df, 1, get.density)

# pivot narrow length, pos dataframe into a kernel matrix: 
# each row is a length, each column is a position
kernel.df <- cast(dens.df, length ~ pos, value='density')
kernel.mat <- as.matrix(kernel.df)


#nuc positions and the linker (nfr) are derived in the dot plot code from the pre-induction sample
nfr = 431300.5
nuc1L = 431410
nuc3R = 432143

mre11_nuc1L_rep1.v = vector()
mre11_nuc1L_rep2.v = vector()

mre11_nuc3R_rep1.v = vector()
mre11_nuc3R_rep2.v = vector()

mre11_nfr_rep1.v = vector()
mre11_nfr_rep2.v = vector()

for (i in 1:6){
  
  #plot(density(merged_df$fsize, bw = 1))
  set.seed(9)
  
  #read in the bam files to data frames
  chr = "chrII"
  df2 = get_dot_mat(as.character(chr2_files2[i]), chr, 1, get_chr_length(as.character(chr2_files2[i]), chr))
  df3 = get_dot_mat(as.character(chr2_files3[i]), chr, 1, get_chr_length(as.character(chr2_files3[i]), chr))
  
  norm.df2 = df2[sample(nrow(df2), sampling_depth, replace = FALSE),]
  norm.df3 = df3[sample(nrow(df3), sampling_depth, replace = FALSE),]
  
  y = nucleosome_size
  
  
  #nuc1L
  mre11_nuc1L_rep1.v = c(mre11_nuc1L_rep1.v,
                         length(
                           which(norm.df2$mpoint<=(nuc1L+std_x) &
                                   norm.df2$mpoint>=(nuc1L-std_x) &
                                   norm.df2$fsize>=(y-std_y) &
                                   norm.df2$fsize<=(y+std_y))
                         )
  )
  
  mre11_nuc1L_rep2.v = c(mre11_nuc1L_rep2.v,
                         length(
                           which(norm.df3$mpoint<=(nuc1L+std_x) &
                                   norm.df3$mpoint>=(nuc1L-std_x) &
                                   norm.df3$fsize>=(y-std_y) &
                                   norm.df3$fsize<=(y+std_y))
                         )
  )
  
  
  #nuc 3R
  
  mre11_nuc3R_rep1.v = c(mre11_nuc3R_rep1.v,
                         length(
                           which(norm.df2$mpoint<=(nuc3R+std_x) &
                                   norm.df2$mpoint>=(nuc3R-std_x) &
                                   norm.df2$fsize>=(y-std_y) &
                                   norm.df2$fsize<=(y+std_y))
                         )
  )
  
  mre11_nuc3R_rep2.v = c(mre11_nuc3R_rep2.v,
                         length(
                           which(norm.df3$mpoint<=(nuc3R+std_x) &
                                   norm.df3$mpoint>=(nuc3R-std_x) &
                                   norm.df3$fsize>=(y-std_y) &
                                   norm.df3$fsize<=(y+std_y))
                         )
  )
  
  #linker
  mre11_nfr_rep1.v = c(mre11_nfr_rep1.v,
                       length(
                         which(norm.df2$mpoint<=(nfr+std_x) &
                                 norm.df2$mpoint>=(nfr-std_x) &
                                 norm.df2$fsize>=(y-std_y) &
                                 norm.df2$fsize<=(y+std_y))
                       )
  )
  
  mre11_nfr_rep2.v = c(mre11_nfr_rep2.v,
                       length(
                         which(norm.df3$mpoint<=(nfr+std_x) &
                                 norm.df3$mpoint>=(nfr-std_x) &
                                 norm.df3$fsize>=(y-std_y) &
                                 norm.df3$fsize<=(y+std_y))
                       )
  )
  
  
  
  
  print(titles[i])
  
}

for (n in 1:6){

  print(
     
t.test(
  c(log2(nuc1L_rep1.v[n]/nuc1L_rep1.v[1]), log2(nuc1L_rep2.v[n]/nuc1L_rep2.v[1])), 
       c(log2(mre11_nuc1L_rep1.v[n]/mre11_nuc1L_rep1.v[1]), log2(mre11_nuc1L_rep2.v[n]/mre11_nuc1L_rep2.v[1]))
       , alternative = "less")$p.value
  ) 
  
}
