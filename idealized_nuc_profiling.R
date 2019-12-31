#extract parameters for idealized nucleosome kernel from all nucs on chrIV

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

depth.v = vector()

for (i in 1:6){
  
  #read in the bam files to data frames
  chr = "chrIV"
  df2 = get_dot_mat(as.character(sacCer3_files2[i]), chr, 1, get_chr_length(as.character(sacCer3_files2[i]), chr))
  df3 = get_dot_mat(as.character(sacCer3_files3[i]), chr, 1, get_chr_length(as.character(sacCer3_files3[i]), chr))
  
  df_2 = rbind(df2, df3)
  
  depth.v = c(depth.v, dim(df_2)[1])
}

sampling_depth = min(depth.v)

#get the unique nuc positions from the broggard published data set
#https://www.nature.com/articles/nature11142
#https://media.nature.com/original/nature-assets/nature/journal/v486/n7404/extref/nature11142-s2.txt

broggard_unique.df <- fread("https://media.nature.com/original/nature-assets/nature/journal/v486/n7404/extref/nature11142-s2.txt")

colnames(broggard_unique.df) = c("chrom", "pos", "NCPscore", "NCPscorenoiseratio")

chr4_nucs_pos = broggard_unique.df$pos[which(broggard_unique.df$chrom == "chrIV")]

my_nuc.m = matrix(0, nrow = 250, ncol = 401)

df2 = get_dot_mat(as.character(sacCer3_files2[1]), chr, 1, get_chr_length(as.character(sacCer3_files2[1]), chr))
df3 = get_dot_mat(as.character(sacCer3_files3[1]), chr, 1, get_chr_length(as.character(sacCer3_files3[1]), chr))

df_2 = rbind(df2, df3)

reads.df = df_2[sample(nrow(df_2), sampling_depth, replace = FALSE),]

for (i in 1:length(chr4_nucs_pos)){
  start = chr4_nucs_pos[i]-200
  end = chr4_nucs_pos[i] + 200

reads.m = get_matrix_from_df(reads.df, start, end)

my_nuc.m = my_nuc.m+reads.m

total <- length(chr4_nucs_pos)
# create progress bar
pb <- txtProgressBar(min = 0, max = total, style = 3)
#    Sys.sleep(0.1)
# update progress bar
setTxtProgressBar(pb, i)

}

close(pb)

#subset the matrix on what we want to plot


levelplot(t(my_nuc.m),
          main = "ChrIV Aggregate Nuc Plot",
          xlab = "Pos From Nucleosome Center",
          ylab = "Fragment Size",
          colorkey = list(labels = list(col="white")
            ),
          scales=
            list(
              x= list(at = seq(0,400,50),
                      labels = seq(-200,200,50)),
              y = list(at = seq(0,250,50),
                       labels = seq(0,250,50))
            )
)
              
plot(x= (-200:200), 
     y = colSums(my_nuc.m),
     type = "l",
     ylim = c(min(colSums(my_nuc.m)), max(colSums(my_nuc.m))),
     xlim = c(-200,200),
     xaxs = "i"
     )

#subset the matrix on what we want to plot

subset.m = my_nuc.m[,101:301]


levelplot(t(subset.m),
          main = "ChrIV Aggregate Nuc Plot",
          xlab = "Pos From Nucleosome Center",
          ylab = "Fragment Size",
          colorkey = list(labels = list(col="white")
          ),
          scales=
            list(
              x= list(at = seq(0,200,25),
                      labels = seq(-100,100,25)),
              y = list(at = seq(0,250,50),
                       labels = seq(0,250,50))
            )
)

plot(x= (-100:100), 
     y = colSums(subset.m),
     type = "l",
     lwd = 2,
     col = "blue",
     ylim = c(min(colSums(subset.m)), max(colSums(subset.m))),
     xlim = c(-200,200),
     xaxs = "i",
     xlab = "Distance From Nucleosome Center",
     ylab = "Density",
     main = "1D Nucleosome Positional Density (ChrIV)",
     cex = 1.5,
     cex.axis = 1.5,
     cex.lab=1.5
)

plot(x= (1:250), 
     y = rowSums(subset.m),
     type = "l",
     lwd = 2,
     col = "red",
     ylim = c(min(rowSums(subset.m)), max(rowSums(subset.m))),
     xlim = c(0,250),
     xaxs = "i",
     xlab = "Fragment Size (MNase)",
     ylab = "Density",
     main = "1D Nucleosomal Fragment Size Density (ChrIV)",
     cex = 1.5,
     cex.axis = 1.5,
     cex.lab=1.5
)

#we now need fragment size mean for the nucleosomes and var and same for the positioning

#use rep to rep the index by the number of times it appears and then use this to calculate var

x_pos = colSums(subset.m)
index = seq(-100,100,1)

nuc_x_pos = rep(index, x_pos)

var(nuc_x_pos)
#use the sqrt of this to get std dev

y_pos = rowSums(subset.m)[Mode(nuc_y_pos):250]
index = seq(Mode(nuc_y_pos),250,1)

nuc_y_pos = rep(index, y_pos)
var(nuc_y_pos)
#use the sqrt of this to get stdev

nucleosome_size = Mode(nuc_y_pos)