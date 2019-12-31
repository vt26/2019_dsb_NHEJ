#all the packages and functiones needed to create plots

library(Rsamtools)
library(graphics)
library(repr)
library(reshape)
library(reticulate)
library(shape)
library(pracma)
library(lattice)
library(stats)
library(KernSmooth)
library(mvtnorm)
library(RColorBrewer)
library(reticulate)
library(Hmisc)

use_python('/usr/lib/python2.7/')
signal<-import('scipy.signal')

library(repr)
library(KernSmooth)
library(RColorBrewer)

library(data.table)
yeast_gene.df <- fread('http://macalpine.vm.duke.edu/~vt26/yeast_gene.df')

densColsDM <- function(x, y = NULL, nbin = 128, bandwidth, transformation = function(x) x^1, colramp = colorRampPalette(blues9), z_factor = 1) 
{
  xy <- xy.coords(x, y)
  select <- is.finite(xy$x) & is.finite(xy$y)
  x <- cbind(xy$x, xy$y)[select, ]
  map <- grDevices:::.smoothScatterCalcDensity(x, nbin, bandwidth)
  mkBreaks <- function(u) u - diff(range(u))/(length(u) - 1)/2
  xbin <- cut(x[, 1], mkBreaks(map$x1), labels = FALSE)
  ybin <- cut(x[, 2], mkBreaks(map$x2), labels = FALSE)
  dens <- map$fhat[cbind(xbin, ybin)]
  dens[is.na(dens)] <- 0
  dens[] <- transformation(dens)
  #print(length(dens))
  colpal <- cut(dens, length(dens), labels = FALSE)
  colpal<-ceiling(as.integer(colpal/z_factor)+1)
  #print(length(colpal))
  cols <- rep(NA_character_, length(select))
  #print(length(cols))
  cols[select] <- colramp(length(dens))[colpal]
  #print(length(cols[select]))
  #print(length(colramp(length(dens))[colpal]))
  cols
}


get_dot_mat <-
  function(bam_file_name, chr, start_pos, end_pos){
    
    # Load the bam_files
    chr.gr = convert_paired_read_bam_to_gr(bam_file_name, chr, start_pos, end_pos)
    
    # Subset on reads less than 250 bp
    chr.gr = chr.gr[width(chr.gr) <= 250]
    
    # Get a new data frame that will contain the positions of the midpoints reads
    mid_points = data.frame(fsize = width(chr.gr), mpoint = start(chr.gr) + floor(width(chr.gr)/2))
    
    # Return the mat.m
    return(mid_points)
    
  }

get_matrix_from_df <-
  function(df, start_pos, end_pos){
    
    
    # Set the query
    
    my_reads = which(df$mpoint<end_pos & df$mpoint>start_pos)
    
    # Subset on the indices that overlap with the query
    
    my_reads.df = df[my_reads,]
    
    # Set up the matrix
    mat.m = matrix(0, nrow = 250, ncol = end_pos - start_pos + 1)
    colnames(mat.m) = start_pos:end_pos
    
    # Iterate through each fragment width
    for(i in 20:250){
      
      # Get the reads that have a particular fragment width
      idx = which(my_reads.df$fsize== i)
      
      if(any(idx)){
        x_positions = my_reads.df[idx,]$mpoint-start_pos
        # Count the overlaps with the mat.gr-
        mat.m[i,x_positions] = mat.m[i,x_positions]+1
        
      }
      
    }
    
    # Return the mat.m
    return(mat.m)
    
  }


MakeArrowSchematic2 <-
  function(feature_chr, feature_start, feature_end,
           cex_title = 1, bg_type = "white",
           fwd_gene_col = "gray", rev_gene_col = "gray",
           proteinCoding = T, geneName = T, omit_genes = NA
  ){
    
    # Set up the plot
    plot(0, 0, type = "n", bty = "n", bg = bg_type,
         xlim = c(feature_start, feature_end), xaxs = "i", xaxt = "n",
         ylim = c(0, 1), yaxs = "i", yaxt = "n",
         ann = F
    )
    
    # Subset only on protein coding if selected
    if(proteinCoding){
      
      idx = which(as.character(yeast_gene.df$name) != as.character(yeast_gene.df$sgd_name))
      
      yeast_gene.df = yeast_gene.df[idx,]
      
    }
    
    # Omit any gene if necessary
    if(any(!is.na(omit_genes))){
      
      yeast_gene.df = yeast_gene.df[-which(yeast_gene.df$sgd_name %in% omit_genes),]
      
    }
    
    # Convert to a GenomicRanges object
    gene.gr = GenomicRanges::GRanges(seqnames = yeast_gene.df$chr,
                                     ranges = IRanges::IRanges(start = yeast_gene.df$start, end = yeast_gene.df$end),
                                     strand = yeast_gene.df$strand
    )
    names(gene.gr) = yeast_gene.df$name
    
    # Create the feature gr
    feature.gr = GenomicRanges::GRanges(seqnames = feature_chr,
                                        ranges = IRanges::IRanges(start = feature_start, end = feature_end)
    )
    
    # Find the overlaps
    overlaps.hits = GenomicRanges::findOverlaps(feature.gr, gene.gr)
    
    if(any(S4Vectors::subjectHits(overlaps.hits))){
      
      # Get the subjectHits
      subject_hits.v = S4Vectors::subjectHits(overlaps.hits)
      
      # Enter in the genes
      for(i in 1:length(subject_hits.v)){
        PlotGene2(yeast_gene.df[subject_hits.v[i],], y_low = 0, y_high = 1,
                  feature_start, feature_end, cex_title, geneName, x_pos_title = 50, fwd_gene_col, rev_gene_col)
      }
      
    }
    
  }



PlotGene2 <-
  function(gene.v, y_low, y_high, x_start, x_end,
           cex_title, geneName, x_pos_title = 50,
           fwd_gene_color, rev_gene_color
  ){
    
    # Get y_mid
    y_mid = (y_high + y_low) / 2
    
    # Add in the text
    if(gene.v$strand == "+"){
      
      # Make the rectangle
      rect(gene.v$start, y_mid + 0.1, gene.v$end, y_high - 0.15, col = fwd_gene_color)
      
      segments(gene.v$start, y_mid + 0.1, gene.v$start, y_low + 0.4,
               col = "black", lty = "solid", lwd = 1
      )
      segments(gene.v$start, y_low + 0.4, gene.v$start+150, y_low + 0.4,
               col = "black", lty = "solid", lwd = 1
      )
      
      if (x_end<gene.v$end){
        Arrowhead(gene.v$start+150, y_low + 0.4, 
                  arr.length = 0.25,
                  arr.width = .25,
                  arr.adj = 1,
                  arr.type = "triangle", 
                  arr.col = "black",
                  arr.lwd = 1,
                  lcol = "black")
      }else{
        Arrowhead(gene.v$start+150, y_low + 0.4, 
                  arr.length = 0.25,
                  arr.width = 0.25,
                  arr.adj = 1,
                  arr.type = "triangle",
                  arr.col = "black",
                  arr.lwd = 1,
                  lcol = "black")
      }
      
      if(geneName){
        if(gene.v$start >= x_start){
          text(x = gene.v$start + x_pos_title, y = y_high - 0.2, adj = c(0, 1),
               labels = gene.v$sgd_name, font = 3, cex = cex_title)
        }else{
          text(x = gene.v$end - x_pos_title, y = y_high - 0.2, adj = c(1, 1),
               labels = gene.v$sgd_name, font = 3, cex = cex_title)
        }
      }
    }else{
      
      # Make the rectangle
      rect(gene.v$start, y_low + 0.15, gene.v$end, y_mid - 0.1, col = rev_gene_color)
      
      segments(gene.v$end, y_mid - 0.1, gene.v$end, y_high - 0.4,
               col = "black", lty = "solid", lwd = 1
      )
      
      segments(gene.v$end, y_high - 0.4, gene.v$end-150, y_high - 0.4,
               col = "black", lty = "solid", lwd = 1
      )
      
      if (x_start>gene.v$start){
        Arrowhead(gene.v$end-150, y_high - 0.4, 
                  arr.length = 0.25,
                  arr.width = .25,
                  arr.adj = 1,
                  arr.type = "triangle", 
                  arr.col = "black",
                  arr.lwd = 1,
                  lcol = "black",
                  angle=180)
      }else{
        Arrowhead(gene.v$end-150, y_high - 0.4, 
                  arr.length = 0.25,
                  arr.width = 0.25,
                  arr.adj = 1,
                  arr.type = "triangle",
                  arr.col = "black",
                  arr.lwd = 1,
                  lcol = "black",
                  angle=180)
      }
      
      
      if(geneName){
        if(gene.v$end <= x_end){
          text(x = gene.v$end - x_pos_title, y = y_mid - 0.15, adj = c(1, 1),
               labels = gene.v$sgd_name, srt = 0, font = 3, cex = cex_title)
        }else{
          text(x = gene.v$start + x_pos_title, y = y_mid - 0.15, adj = c(1, 1),
               labels = gene.v$sgd_name, srt = 0, font = 3, cex = cex_title)
        }
      }
    }
    
  }



DensDotPlotVT <-
  function(dot.m, z_min = 0, z_max = NA,
           low_col = "white", med_col = "", high_col = "blue", num_colors = 100,
           plot_title = "", plot_title_line = NA,
           x_label = "", y_label = "",
           x_axt = "s", y_axt = "s",
           plot_box = TRUE){
    
    # Check that dot.m is a matrix
    if(!is.matrix(dot.m)){
      stop("dot.m is not a matrix, DensDotPlot requires a matrix!")
    }
    
    if(is.na(z_max)){
      # Make the z_max equivalent to the 95th percentile if z_max is not specified
      z_max = quantile(as.vector(dot.m), probs = 0.95)
      
    }
    
    # For points that are either above or below z_max or z_min respectively, set them to
    # the z_min and z_max (otherwise, plot shows arbitrary colors
    dot.m[which(dot.m >= z_max)] = z_max
    dot.m[which(dot.m <= z_min)] = z_min
    
    # Get the current column names
    if(!is.null(colnames(dot.m))){
      
      position.v = as.numeric(colnames(dot.m))
      
    }else{
      
      position.v = 1:ncol(dot.m)
      
    }
    
    # Get the row values
    row.v = 1:nrow(dot.m)
    
    # Make the colorpanel
    if(med_col == ""){
      make_colorpanel = gplots::colorpanel(n = num_colors, low = low_col, high = high_col)
    }else{
      make_colorpanel = gplots::colorpanel(n = num_colors, low = low_col, med = med_col, high = high_col)
    }
    
    # Make the heatmap utilizing the parameters specified above
    image(position.v, row.v, t(dot.m[nrow(dot.m):1,]), col = make_colorpanel, zlim = c(z_min, z_max),
          xlab = x_label, ylab = y_label, xaxt = x_axt, yaxt = y_axt, bty = "n", axes = F
    )
    
    # Set the title
    title(main = plot_title, line = plot_title_line)
    
    # Add a box around the plot
    if(plot_box){
      box(which = "plot", lty = "solid")
    }
    
  }

#make new gene df with +117 for all genes beyond ho-site

yeast_gene.df2 = yeast_gene.df

yeast_gene.df2$start[which(yeast_gene.df2$chr==2 & yeast_gene.df2$start>ho_end)] = 
  yeast_gene.df2$start[which(yeast_gene.df2$chr==2 & yeast_gene.df2$start>ho_end)]+117 

yeast_gene.df2$end[which(yeast_gene.df2$chr==2 & yeast_gene.df2$start>ho_end)] = 
  yeast_gene.df2$end[which(yeast_gene.df2$chr==2 & yeast_gene.df2$start>ho_end)]+117 


MakeArrowSchematic_ho <-
  function(feature_chr, feature_start, feature_end,
           cex_title = 1, bg_type = "white",
           fwd_gene_col = "gray", rev_gene_col = "gray",
           proteinCoding = T, geneName = T, omit_genes = NA
  ){
    
    # Set up the plot
    plot(0, 0, type = "n", bty = "n", bg = bg_type,
         xlim = c(feature_start, feature_end), xaxs = "i", xaxt = "n",
         ylim = c(0, 1), yaxs = "i", yaxt = "n",
         ann = F
    )
    
    # Subset only on protein coding if selected
    if(proteinCoding){
      
      idx = which(as.character(yeast_gene.df2$name) != as.character(yeast_gene.df2$sgd_name))
      
      yeast_gene.df2 = yeast_gene.df2[idx,]
      
    }
    
    # Omit any gene if necessary
    if(any(!is.na(omit_genes))){
      
      yeast_gene.df2 = yeast_gene.df2[-which(yeast_gene.df2$sgd_name %in% omit_genes),]
      
    }
    
    # Convert to a GenomicRanges object
    gene.gr = GenomicRanges::GRanges(seqnames = yeast_gene.df2$chr,
                                     ranges = IRanges::IRanges(start = yeast_gene.df2$start, end = yeast_gene.df2$end),
                                     strand = yeast_gene.df2$strand
    )
    names(gene.gr) = yeast_gene.df2$name
    
    # Create the feature gr
    feature.gr = GenomicRanges::GRanges(seqnames = feature_chr,
                                        ranges = IRanges::IRanges(start = feature_start, end = feature_end)
    )
    
    # Find the overlaps
    overlaps.hits = GenomicRanges::findOverlaps(feature.gr, gene.gr)
    
    if(any(S4Vectors::subjectHits(overlaps.hits))){
      
      # Get the subjectHits
      subject_hits.v = S4Vectors::subjectHits(overlaps.hits)
      
      # Enter in the genes
      for(i in 1:length(subject_hits.v)){
        PlotGene2(yeast_gene.df2[subject_hits.v[i],], y_low = 0, y_high = 1,
                  feature_start, feature_end, cex_title, geneName, x_pos_title = 50, fwd_gene_col, rev_gene_col)
      }
      
    }
    
  }




# Function to convert a BAM file to a GR file
# By default reads in all the reads for the entire chromosome
# 	Otherwise, specify a particular start_pos and end_pos for just a specific region
convert_paired_read_bam_to_gr = function(bam_file_name, chr, start_pos = 1, end_pos = -1){
  
  # Create the BAM File object
  bf = BamFile(bam_file_name, index = paste(bam_file_name, ".bai", sep = ""))
  
  # Get the chr list
  chr_length.v = scanBamHeader(bf)$targets
  
  # Update the end_pos if necessary
  if(end_pos == -1){
    
    # Update the end_pos
    end_pos = chr_length.v[chr]
    
  }
  
  # Make a GR file for the chromosome
  chr.gr = GRanges(seqnames = chr, 
                   ranges = IRanges(start = max(start_pos - 250, 1), end = end_pos), 
                   strand = "*"
  )
  
  # Specify the scan bam paramaeters
  p = ScanBamParam(what = c("pos", "isize"), which = chr.gr, flag = scanBamFlag(isMinusStrand = FALSE))
  
  # Get the reads that meet these conditions
  reads.l = scanBam(bf, param = p)
  
  if(length(reads.l[[1]][["pos"]]) > 0){
    
    # Convert these reads to a GR object
    IP.gr = GRanges(seqnames = factor(chr, levels = names(chr_length.v)),
                    ranges = IRanges(start = reads.l[[1]][["pos"]], width = reads.l[[1]][["isize"]]),
                    strand = "*"
    )
    seqlengths(IP.gr) = chr_length.v
    
  }else{
    
    IP.gr = GRanges()
    
  }
  
  return(IP.gr)
}	

get_chr_length = function(bam_file, chr){
  
  # Get the summary of the .bai file
  bai_file.v = system(command = paste("samtools idxstats ", bam_file, sep = ""), intern = T)
  
  # Convert to list
  bai_file.l = strsplit(bai_file.v, split = "\t")
  
  # Convert to a dataframe
  bai_file.df = as.data.frame(matrix(unlist(bai_file.l), ncol = 4, byrow = T))
  colnames(bai_file.df) = c("chr", "length", "read_num", "unaligned")
  
  # Convert the read number into a numeric
  bai_file.df$length = as.numeric(as.character(bai_file.df$length))
  
  # Get the total count number over the called chr
  return(bai_file.df$length[which(bai_file.df$chr == chr)])
  
}