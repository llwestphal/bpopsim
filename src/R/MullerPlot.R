#/usr/bin/R

## Set the working directory
setwd('.')

##Get filename from command line
muller_matrix_file_name = "muller_matrix.dat";

############################################################################
# This script is for plotting the muller_matrix.dat file
############################################################################

require(grDevices)

############################################################################
# Matrix manipulation methods
############################################################################
# Flip matrix (upside-down)
flip.matrix <- function(x) {
  mirror.matrix(rotate180.matrix(x))
}

# Mirror matrix (left-right)
mirror.matrix <- function(x) {
  xx <- as.data.frame(x);
  xx <- rev(xx);
  xx <- as.matrix(xx);
  xx;
}

# Rotate matrix 90 clockworks
rotate90.matrix <- function(x) {
  t(mirror.matrix(x))
}

# Rotate matrix 180 clockworks
rotate180.matrix <- function(x) { 
  xx <- rev(x);
  dim(xx) <- dim(x);
  xx;
}

# Rotate matrix 270 clockworks
rotate270.matrix <- function(x) {
  mirror.matrix(t(x))
}

############################################################################
# Draw methods
############################################################################
# 'colorTable' instead of 'col' to be more explicit.
image.matrix <- function(x, colorTable=NULL, xlab="x", ylab="y", ...) {
  image(x=1:ncol(x),y=1:nrow(x),z=rotate270.matrix(x), col=colorTable, xlab=xlab, ylab=ylab, ...)
}

muller_mat <- t(as.matrix(read.table(muller_matrix_file_name, header=F)))*10

png("Muller.png", 
    width=2000,
    height=500)

#opar <- par(ask=TRUE)

colors = c(hsv(runif(10000), 1, .65))

par(mar=c(6,1,1,1))
par(mgp=c(4.75,2,0))

xticksnames=seq(0,20000, by=1000)
xticks=xticksnames/6.64

image.matrix(muller_mat,
             colorTable = colors,
             ylab="",
             xlab="Time (Bacterial Generations)",
             axes=FALSE,
             xaxt="n",
             yaxt="n",
             cex.lab=2.75)

axis(1,
     at=xticks,
     labels=xticksnames,
     col.axis="black", 
     adj=0,
     lwd=5,
     tck=-0.05,
     cex.axis=2)

# title(main='Bpopsim Muller Matrix', 
#       col.main="black", 
#       cex.main = 5)

dev.off()
