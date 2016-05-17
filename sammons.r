library(MASS)
cosineDist <- function(x){
  as.dist(1 - x%*%t(x)/(sqrt(rowSums(x^2) %*% t(rowSums(x^2))))) 
}

drawCluster <- function(label, norm,  method, frequencys, sam){
  points <- as.data.frame(sam$points)
  filename = paste(method, "clusters", sep=".")
  filename = paste(path, filename, sep="/")
  if(file.exists(filename)){
    clusters = read.table(filename)
    clusters <- clusters[frequencys$id,] + 1
    points$class <- clusters

    pdf_path <- paste(label, method, norm, "pdf", sep=".")
    pdf_path <- paste(path, pdf_path, sep="/")
    print(paste("Output image to", pdf_path, sep=" "))
    pdf(file=pdf_path)

    print("Plot results")
    plot(points$V1, points$V2, col=clusters)
    dev.off()

    # How to plot each cluster separated and respecting fuzzy partitions
    #plot(table.sam$points, col=grey(1-m$V6), pch=16)
    #plot(table.sam$points, col=grey(1-m$V1), pch=16)
    #plot(table.sam$points, col=grey(1-m$V2), pch=16)
    #plot(table.sam$points, col=grey(1-m$V3), pch=16)
    #plot(table.sam$points, col=grey(1-m$V4), pch=16)
    #m <- read.table("memberships.matrix")
    #text(points, labels = as.character(1:nrow(frequencys)))
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
path <- args[3]

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
