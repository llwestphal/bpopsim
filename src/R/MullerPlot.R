#/usr/bin/R

options <- commandArgs(trailingOnly = TRUE)

##Get set commandline arguments
input_file_name = options[1];
output_file_name = options[2]
output_file_type = options[3]

## Set the working directory
setwd('.')

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

is.adjacent <- function(x, y, vec) {
  pos_x <- which(vec == x)
  pos_y <- which(vec == y)
  
  for(i in pos_y) {
    if(1 %in% abs(i - pos_x)) 
      return(T)
  }
  
  return(F)
}

check.colors <- function(color_vec, adj, color_pallete) {
  count <- 1
  for(i in 1:(length(color_vec) - 1)) {
    for(j in (i+1):length(color_vec)) {
      if(adj[count] & color_vec[i] == color_vec[j] & color_vec[i] != 0 & color_vec[j] != 0 
         & color_vec[i] != color_pallete[length(color_pallete)] 
         & color_vec[j] != color_pallete[length(color_pallete)]) {
        return(T)
      }
      count <- count + 1
    }
  }
  
  return(F)
}

muller_mat <- t(as.matrix(read.table(input_file_name, header=F)))*10

if( output_file_type == "pdf") {
  pdf(paste(output_file_name, ".pdf", sep = ''), 
      width = 20,
      height = 5)
}

if( output_file_type == "svg") {
  svg(paste(output_file_name, ".svg", sep = ''), 
      width = 20,
      height = 5)
}

if( output_file_type == "png") {
  png(paste(output_file_name, ".png", sep = ''), 
      width = 20,
      height = 5, 
      units = "in",
      res = 300)
}

num_pixels = length(muller_mat[, 1])
num_times = length(muller_mat[1, ])
small_genotypes = c()
grown_up = c()

## I go through the muller matrix looking for genotypes that never rise above 1% frequency
## Then, I reset those that do not rise high enough to the number 5 and thus gray
## There is probably a more efficient way to do this
for( i in 1:num_times ) {
  this_time <- muller_mat[, i]
  freq_frame <- as.data.frame(table(this_time) / num_pixels)
  colnames(freq_frame) <- c("Genotype", "Freq")
  freq_frame$Genotype = as.numeric(as.character(freq_frame$Genotype))
  
  small_genotypes = append(small_genotypes, freq_frame$Genotype)
  grown_up  = append(grown_up, freq_frame$Genotype[freq_frame$Freq > 0.01])
}

grown_up <- unique(grown_up)
small_genotypes <- unique(small_genotypes)
small_genotypes = setdiff(small_genotypes, grown_up)

print(sort(grown_up))
#print(sort(small_genotypes))

## Here I figure out which genotypes are adjacent to each other
## Again there is probably a little more efficient way to do it
## For loops in R are generally not a good idea

num <- 0
for( i in 1:(length(grown_up) - 1) ) {
  num <- num + i
}

adjacency <- rep(F, num)

for( i in 1:length(muller_mat[1, ]) ) {
  this_time <- muller_mat[, i]
  
  counter <- 1
  
  for( j in 1:(length(grown_up) - 1) ){
    for( k in (j+1):length(grown_up) ) {
      if( !adjacency[counter] ) {
        adjacency[counter] = is.adjacent(grown_up[j], grown_up[k], this_time)
      }
      counter <- counter + 1
    }
  }
}
## Now that we have colored all small genotypes gray we will start coloring the rest of the plot

#colors = hsv(seq(0, 0.5, length.out = length(grown_up)), s = 0.75, v = 0.65)
#colors <- colors()[1:100]
cols <- c('black', 'green', 'darkred', 'tan', 'lightgreen', 'orange', 'purple',
          'red', 'gold', 'darkgreen', 'brown')

cols = append(cols, 'darkgray')

col_vec <- rep(0, length(grown_up))
col_vec[length(col_vec)] <- length(cols)

this_color <- 1
counter <- 1

col_vec[1] = 1

for( i in 1:(length(grown_up) - 1 )) {
  for( j in (i+1):length(grown_up) ) {
    if( !adjacency[counter] & col_vec[i] == 0) {
      col_vec[j] = this_color
    }
    
    else if( col_vec[j] == 0 ) {
      keep_trying <- T
      permutations <- 0
      
      while(keep_trying) {
        this_color <- this_color + 1
        if( this_color >= length(cols) )
          this_color = 1
        
        col_vec[j] = this_color
        if( !check.colors(col_vec, adjacency, cols)) {
          break
        }
        else{
          col_vec[j] = 0
        }
        
        if(permutations >= length(cols)) {
          stop("Color possibilities exhausted... You will have to add more colors to prevent color adjacency.")
        }
        
        print(this_color)
        permutations <- permutations + 1
      }
    }
  }
}


muller_mat[which(muller_mat %in% small_genotypes)] = length(cols)
for( i in 1:length(grown_up)) {
  muller_mat[muller_mat == grown_up[i]] = col_vec[i]
}

print(check.colors(col_vec, adjacency, cols))

print(col_vec)

par(mar=c(6,1,1,1))
par(mgp=c(4.75,2,0))

xticksnames=seq(0,20000, by=1000)
xticks=xticksnames/6.64

image.matrix(muller_mat,
             colorTable = cols,
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
