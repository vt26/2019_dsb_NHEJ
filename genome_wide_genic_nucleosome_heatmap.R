#make a genome-wide heatmap of nucleosome fuzziness for every single genic nucleosome

set.seed(9)

titles = c("Pre-Induction",
           "15 Mins",
           "30 Mins",
           "60 Mins",
           "90 Mins",
           "120 Mins"
           
)


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

nuc_fuzz1.l = list()

for (c in 1:16){
  if (c==2){
    files = chr2_files
  }else{
    files = sacCer3_files
  }
  
  depth.v = vector()
  for (i in seq(1,6,1)){
    
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
  
  #normalize to read depth
  
  #chr_cor.v = (chr_cor.v)*10000/dim(norm.df)[1]
  
  peaks.df = data.frame(findpeaks(chr_cor.v, 
                                  minpeakdistance = 75, 
                                  minpeakheight = 0.011578938, 
                                  sortstr= TRUE,
                                  nups= 1))
  #plot(smooth.spline(cor.v, spar = .001)$y)
  colnames(peaks.df) = c("height", "center", "start", "end")
  
  peaks.v = sort(peaks.df$center)
  
  nuc_fuzz1.l[[c]] = chr_cor.v[peaks.v]
  
  nuc_pos1.l[[c]] = peaks.v
}

#now we have all the nuc positions on all the chromosomes in a list

#look at each chromosome and each gene

all_chroms_corr.l = list()

for (c in 1:16){
  
  genes_on_chr.df = yeast_gene.df2[which(yeast_gene.df2$chr==c),]
  
  genes_on_chr.df$t0 = 0
  genes_on_chr.df$t15 = 0
  genes_on_chr.df$t30 = 0
  genes_on_chr.df$t60 = 0
  genes_on_chr.df$t90 = 0
  genes_on_chr.df$t120 = 0
  
  if (c==2){
    files = chr2_files
  }else{
    files = sacCer3_files
  }
  
  depth.v = vector()
  for (i in seq(1,6,1)){
    
    #read in the bam files to data frames
    chr = chroms.v[c]
    df1 = get_dot_mat(as.character(files[i]), chr, 1, get_chr_length(as.character(files[i]), chr))
    df2 = get_dot_mat(as.character(files[i+6]), chr, 1, get_chr_length(as.character(files[i+6]), chr))
    
    merged.df = rbind(df1, df2)
    
    depth.v = c(depth.v, dim(merged.df)[1])
  }
  
  sampling_depth = min(depth.v)
  
  
  for (a in 1:6){
    
    set.seed(9)
    
    chr = chroms.v[c]
    df1 = get_dot_mat(as.character(files[a]), chr, 1, get_chr_length(as.character(files[a]), chr))
    df2 = get_dot_mat(as.character(files[a+6]), chr, 1, get_chr_length(as.character(files[a+6]), chr))
    
    merged.df = rbind(df1, df2)
    
    
    norm.df = merged.df[sample(nrow(merged.df), sampling_depth, replace = FALSE),]
    
    reads.m = get_matrix_from_df(norm.df, 
                                 1, get_chr_length(as.character(files[a]), chr))
    
    
    score = signal$correlate2d(reads.m,kernel.mat, mode='valid')
    chr_cor.v = as.vector(score[1,])
    
    #shift the chrIV_cor.v vector by padding 75bp
    
    chr_cor.v = c(rep(0,75),chr_cor.v)
    
    for (n in 1:(dim(genes_on_chr.df)[1])){
      
      genic_nucs.idx = 
        which(nuc_pos1.l[[c]] > genes_on_chr.df$start[n] & 
                nuc_pos1.l[[c]] < (genes_on_chr.df$end[n]))
      
      if (length(genic_nucs.idx)>0){
        nuc_cor.idx = nuc_pos1.l[[c]][genic_nucs.idx]
        genes_on_chr.df[n,a+7] = sum(chr_cor.v[nuc_cor.idx])
      }
      
    }
    
  }
  all_chroms_corr.l[[c]] = genes_on_chr.df
  print(c)
}

all_chroms_corr.df = rbind(all_chroms_corr.l[[1]],
                           all_chroms_corr.l[[2]],
                           all_chroms_corr.l[[3]],
                           all_chroms_corr.l[[4]],
                           all_chroms_corr.l[[5]],
                           all_chroms_corr.l[[6]],
                           all_chroms_corr.l[[7]],
                           all_chroms_corr.l[[8]],
                           all_chroms_corr.l[[9]],
                           all_chroms_corr.l[[10]],
                           all_chroms_corr.l[[11]],
                           all_chroms_corr.l[[12]],
                           all_chroms_corr.l[[13]],
                           all_chroms_corr.l[[14]],
                           all_chroms_corr.l[[15]],
                           all_chroms_corr.l[[16]])

means.v = colSums(all_chroms_corr.df[,8:13])/dim(all_chroms_corr.df)[1]

all_genes_mean_norm2.m = sweep(all_chroms_corr.df[,8:13], 2, means.v, "/")

file_name = paste("genome_wide_nuc_corr_chr_order_retest.png", sep="")
png(file_name, width = 4, height = 12, units = "in", res = 300)
image(t(
  log2(
    (all_genes_mean_norm2.m+.01)/(all_genes_mean_norm2.m[,1]+.01)
  )
),
col = colorRampPalette(c("blue", "white", "red"))(n=51),
zlim = c(-2,2),
axes = F,
xlab = "Time (Minutes)",
#ylab = "Location"
)

axis(1,
     at = seq(0,1,.2),
     labels = c(0, 15, 30, 60, 90, 120))

chroms_length.v = vector()
for (c in 1:16){
  chroms_length.v = c(chroms_length.v, get_chr_length(as.character(sacCer3_files[1]), chroms.v[c]))
}

chrom_marks.v = vector()
for(i in 1:16){
  chrom_marks.v = c(chrom_marks.v, sum(chroms_length.v[1:i])/sum(chroms_length.v))
}

center_marks.v = vector()
for(i in 1:15){
  center_marks.v = c(center_marks.v, (chrom_marks.v[i] + chrom_marks.v[i+1])/2)
}

center_marks.v = c((0+chrom_marks.v[1])/2, center_marks.v)

axis(2,
     at= center_marks.v,
     labels = chroms.v,
     las = 2,
     tick = F)

axis(2,
     at = c(0,chrom_marks.v),
     labels = F,
     tick = T
)

dev.off()

chr2_genes = which(all_chroms_corr.df$chr==2)

file_name = paste("chrII_nuc_corr_chr_order_test.png", sep="")
png(file_name, width = 4, height = 10, units = "in", res = 300)
image(
  t(
    log2(
      (all_genes_mean_norm2.m[chr2_genes,]+.01)/(all_genes_mean_norm2.m[chr2_genes,1]+.01)
    )
  ),
col = colorRampPalette(c("blue", "white", "red"))(n=51),
zlim = c(-2,2),
axes = F,
xlab = "Time (Minutes)"
#ylab = "Location"
)

axis(1,
     at = seq(0,1,.2),
     labels = c(0, 15, 30, 60, 90, 120))

axis(2,
     at= .95*276253/813184,
     labels = expression(italic("GAL")),
     las = 2,
     tick = F)

axis(2,
     at = .95*((ho_start+ho_end)/2)/813184,
     labels = "HOcs",
     tick = F,
     las = 2
)

dev.off()

#make legend
file_name = paste("heatmap_legend.png", sep="")
png(file_name, width = 3, height = 12, units = "in", res = 300)
par(mar=c(1,1,1,4))
library(autoimage)

legend.scale(c(-2,2), col = colorRampPalette(c("blue", "white", "red"))(n=51), horizontal = F,
             axis.args = list(cex.axis = 3)
             )
dev.off()

