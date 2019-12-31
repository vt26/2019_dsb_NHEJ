library(flowCore)
#make an empty list
f= list()

#loop through all the fcs files and add them to the list

f[[1]] = read.FCS(paste("S1_NHEJ_P2_001.fcs", sep=""), transformation = FALSE)

for (i in 2:7){
  
  f[[i]] = read.FCS(paste("S1_NHEJ_00", i-1, "_P2_00", i, ".fcs", sep=""), transformation = FALSE)
  
}

titles = c("Pre", expression("3.5h "*alpha*"-F"), "1h Gal", "1h Dex", "2h Dex", "4h Dex", "6h Dex")

png(filename= "NHEJ_wt_repeated_08272019_time-course_FACScanto.png", width = 12, height = 2, units = "in", res = 300)
par(mfrow=c(1,7))
for (i in 1:7){
  plot(density(exprs(f[[i]])[,3], bw = 20), xlab = "DNA Content", ylab = "Frequency", 
       main = titles[i],
       #sub = expression(italic(dnl4)*italic(Delta))
  )
}

dev.off()

#make an empty list
f= list()

#loop through all the fcs files and add them to the list

f[[8]] = read.FCS("S1_NHEJ_007_P2_008.fcs", transformation = FALSE)
f[[9]] = read.FCS("S1_NHEJ_008_P2_009.fcs", transformation = FALSE)
f[[10]] = read.FCS("S1_NHEJ_009_P2_010.fcs", transformation = FALSE)
f[[11]] = read.FCS("S1_NHEJ_010_P2_011.fcs", transformation = FALSE)
f[[12]] = read.FCS("S1_NHEJ_011_P2_012.fcs", transformation = FALSE)
f[[13]] = read.FCS("S1_NHEJ_012_P2_013.fcs", transformation = FALSE)
f[[14]] = read.FCS("S1_NHEJ_013_P2_014.fcs", transformation = FALSE)

titles = c("Pre", expression("3.5h "*alpha*"-F"), "1h Gal", "1h Dex", "2h Dex", "4h Dex", "6h Dex")

png(filename= "NHEJ_dnl4_repeated_08272019_time-course_FACScanto.png", width = 12, height = 2, units = "in", res = 300)
par(mfrow=c(1,7))
for (i in 8:14){
  plot(density(exprs(f[[i]])[,3], bw = 20), xlab = "DNA Content", ylab = "Frequency", 
       main = titles[i-7],
       sub = expression(italic(dnl4)*italic(Delta))
        )
}

dev.off()
