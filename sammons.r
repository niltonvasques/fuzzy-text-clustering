library(MASS)
cosineDist <- function(x){
  as.dist(1 - x%*%t(x)/(sqrt(rowSums(x^2) %*% t(rowSums(x^2))))) 
}

drawCluster <- function(label, norm,  method, frequencys, sam){
  points <- as.data.frame(sam$points)
  filename = paste(method, "clusters", sep=".")
  if(file.exists(filename)){
    clusters = read.table(filename)
    clusters <- clusters[frequencys$id,]
    points$class <- clusters

    print("Output image to output.pdf")
    pdf(file=paste(label, method, norm, "pdf", sep="."))

    print("Plot results")
    plot(points$V1, points$V2, col=points$class)
    #text(points, labels = as.character(1:nrow(frequencys)))
    dev.off()
  }
}

applySammonMapping <- function(frequencys, distances, norm, label){
  print("Performing matrix reduction with sammons mapping")
  sam <- sammon(distances)

  drawCluster(label, norm, "soft-fdcl", frequencys, sam)
  drawCluster(label, norm, "mixed-pdcl", frequencys, sam)
  drawCluster(label, norm, "pdcl", frequencys, sam)
}

args <- commandArgs(trailingOnly = TRUE)
freq_filename <- args[1]
title <- args[2]

print("Reading data")
table = read.table(freq_filename)
table <- as.data.frame(table)
table$id <- c(1:nrow(table))

print("Checking for duplicates")
table <- table[!duplicated(table[,c(1:ncol(table)-1)]), ]

print("Computing all distances with cosine measure")
d = cosineDist(as.matrix(table))
applySammonMapping(table, d, "cosine", title)

print("Computing all distances with euclidian measure")
d = dist(as.matrix(table))
applySammonMapping(table, d, "euclidian", title)

print("Process finished")
