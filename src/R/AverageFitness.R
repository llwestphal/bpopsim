#/usr/bin/R

#The comments below this can be uncommented to output a png graphic

average_fit <- (as.matrix(read.table("Location/AverageFitness.dat", sep="\t", header=F)))[1,]

lenxvals <- length(average_fit)
a <- ((lenxvals+1)/10*6.64)

xticks <- round(c(-1000,
             0,
             2*a,
             4*a,
             6*a,
             8*a,
             10*a,
             12*a,
             14*a,
             16*a)/6.64)
             
xticksnames <- round(xticks*6.64*3/1000)

yticks <- c(1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0)

#png("Location/RelativeFitness.png", 
#width=1280, 
#height=1024, 
#units = "px")

#par(mar=c(8,8,8,3))
#par(mgp=c(5,1,0))

plot(average_fit,
  pch=20,
  xlab="Generations (In Thousands)", 
  ylab="Relative Fitness",
  axes=FALSE,
  xaxt="n",
  yaxt="n",
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
  at=yticks,
  labels=yticks,
  adj=0,
  lwd=2,
  tck=.01,
  las=1,
  cex.axis=2)
  
title(main='Average Fitness', 
  col.main="black", 
  cex.main = 3)
  
#dev.off()