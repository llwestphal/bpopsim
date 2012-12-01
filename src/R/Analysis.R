#/usr/bin/R

#This is just some code I wrote to do analysis on the output of bpopsim
#You will probably need to modify it heavily to get it to work with your particular data

#Genotypes_dat <- as.matrix(read.table("Location/Number_Unique_Genotypes.dat", 
#                 sep="\t", header=T))

library(plyr)

location <- "/Users/austin/Desktop/output_test"

#parallelMuts_dat <- (as.list(read.table(paste(location, "/SignificantParallelMutations.dat", sep=""),
#                      sep="\n")))[[1]]
                      
#clumpiness_dat <- (as.list(read.table(paste(location, "/SweepClumpiness.dat", sep=""),
#                      sep="\n")))[[1]]
                      
#time_to_sweep_dat <- ((as.list(read.table(paste(location, "/TimeToSweep.dat", sep=""),
#                      sep="\n")))[[1]]*6.64)


genotype_frequency_dat <- data.frame(as.matrix(read.table(paste(location, "/Genotype_Frequencies.dat", sep=""), fill=TRUE, sep=" ", header=T)))
genotype_frequency_dat <- genotype_frequency_dat[,-length(genotype_frequency_dat)]

colors = c(hsv(runif(1000), 1, .65))
yticks <- seq(-1,1,by=.1)
xticks <- seq(-50000,50000,by=100)

# png(paste(location, "/Genotype_Frequencies.png", sep=""), 
# width=1280, 
# height=1024, 
# units = "px")
    
par(mar=c(6,6,6,4))
par(mgp=c(4,1,0))

plot(genotype_frequency_dat[,1],
     col=colors[1],
     pch=20,
     axes=FALSE,
     ylim=c(0,1),
     yaxt="n",
     #log="y",
     xlab="Time (Generations)", 
     ylab="Genotype Frequency",
     cex.lab=2)

axis(1, at=xticks,
     labels=xticks*6.64,
     col.axis="black", 
     adj=0, tck=.01,
     lwd=2, cex.axis=2)
  
axis(2, at=yticks, 
     col.axis="black", 
     adj=0, tck=.01, lwd=2,
     las=1, cex.axis=2)
  
title(main='Bpopsim Genotype Frequencies', 
      col.main="black", 
      cex.main = 3)

for(i in 2:length(genotype_frequency_dat)) {
  lines(genotype_frequency_dat[,i], 
        col=colors[i],
        pch=20,
        lwd=4)
}

legend_names <- colnames(genotype_frequency_dat)

#legend(0,.99, legend_names, cex=1.5, pch=20, col=colors)

# dev.off()

#This is just an easy way to comment out several lines
if(FALSE) {
totaltime <- (dim(Genotypes_dat))

averageVals = c()
for (time in 3:totaltime[1]) {
  averageAtTime <- mean(Genotypes_dat[time,][3:totaltime[[2]]])
  averageVals[time] <- averageAtTime
}

lenxvals <- dim(Genotypes_dat)[1]
a <- ((lenxvals+1)/10*6.64)

yvals <- (max(averageVals)-min(averageVals))/10

xticks <- (c(-2*a,
             0,
             2*a,
             4*a,
             6*a,
             8*a,
             10*a,
             12*a)/6.64)
             
xticksnames <- (xticks*6.64)
             
yticks <- c(0, 
            2*yvals+min(averageVals), 
            4*yvals+min(averageVals), 
            6*yvals+min(averageVals), 
            8*yvals+min(averageVals),
            10*yvals+min(averageVals),
            12*yvals+min(averageVals))
            
png(paste(location, "/Number_Unique_Genotypes.png", sep=""), 
width=1280, 
height=1024, 
units = "px")

par(mar=c(8,8,8,4))
par(mgp=c(6,1,0))

plot(averageVals,
     pch=20,
     xlab="Time (Transfers)", 
     ylab="Unique Genotypes",
     axes=FALSE,
     yaxt="n",
     xaxt="n",
     cex.lab=2)
  
axis(1,
     at=xticks,
     labels=xticksnames,
     col.axis="black", 
     adj=0,
     lwd=2,
     tck=.01,
     cex.axis=2)
  
axis(2,
     col.axis="black", 
     adj=0,
     lwd=2,
     tck=.01,
     las=1,
     cex.axis=2)
  
title(main='Number of Unique Genotypes', 
     col.main="black", 
     cex.main = 3)

dev.off()
}

png(paste(location, "/SignificantParallelMutations.png", sep=""), 
width=1280, 
height=1024, 
units = "px")

par(mar=c(8,8,8,4))
par(mgp=c(6,1,0))

xticks1 <- c(seq(from=-0.2, to=1, by=.2))
yticks1 <- c(seq(from=-10, to=100, by=1))

hist(parallelMuts_dat,
     breaks=c(seq(from=0, to=1, by=0.05)),
     col=heat.colors(11),
     main='Sweeping Neighbors', 
     col.main="black", 
     cex.main = 3,
     xaxt="n",
     yaxt="n",
     axes=FALSE,
     xlab="Max Difference", 
     ylab="Raw Counts",
     #probability=T,
     cex.lab=2)
     
axis(1,
     col.axis="black", 
     at=xticks1,
     labels=xticks1,
     adj=0,
     lwd=2,
     tck=.02,
     cex.axis=2)
  
axis(2,
     col.axis="black", 
     at=yticks1,
     labels=yticks1,
     adj=0,
     lwd=2,
     tck=.02,
     las=1,
     cex.axis=2)
  
#lines(density(parallelMuts_dat, bw=.009), lwd=2, col="black")

dev.off()

png(paste(location, "/SweepClumpiness.png", sep=""), 
width=1280, 
height=1024, 
units = "px")

par(mar=c(8,8,8,4))
par(mgp=c(6,1,0))

xticks2 <- c(seq(from=-1, to=10, by=1))
yticks2 <- c(seq(from=-20, to=120, by=1))

hist(clumpiness_dat, 
     breaks=c(.5:6.5),
     col=heat.colors(11),
     main='Number of Simultaneous Sweeps', 
     col.main="black", 
     cex.main = 3,
     xaxt="n",
     yaxt="n",
     axes=FALSE,
     xlab="Simultaneous Sweeps", 
     ylab="Raw Counts",
     cex.lab=2,
     #probability=T
     )
     
axis(1,
     col.axis="black", 
     at=xticks2,
     labels=xticks2,
     adj=0,
     lwd=2,
     tck=.02,
     cex.axis=2)
  
axis(2,
     col.axis="black", 
     at=yticks2,
     adj=0,
     lwd=2,
     tck=.02,
     las=1,
     cex.axis=2)
  
dev.off()

png(paste(location, "/TimeToSweep.png", sep=""), 
width=1280, 
height=1024, 
units = "px")

par(mar=c(8,8,8,4))
par(mgp=c(6,1,0))

xticks3 <- c(seq(from=-2000, to=10000, by=2000))
yticks3 <- c(seq(from=-20, to=120, by=1))

hist(time_to_sweep_dat,
     col=heat.colors(11),
     main='Total Time to Sweep', 
     col.main="black", 
     cex.main = 3,
     xaxt="n",
     yaxt="n",
     axes=FALSE,
     xlab="Time to Sweep (Generations)", 
     ylab="Raw Counts",
     cex.lab=2,
     #probability=T
     )
     
axis(1,
     at=xticks3,
     col.axis="black",
     adj=0,
     lwd=2,
     tck=.02,
     cex.axis=2)
  
axis(2,
     at=yticks3,
     col.axis="black", 
     adj=0,
     lwd=2,
     tck=.02,
     las=1,
     cex.axis=2)
  
#lines(density(time_to_sweep_dat), lwd=2, col="black")
  
dev.off()