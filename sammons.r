library(MASS)
cosineDist <- function(x){
  as.dist(1 - x%*%t(x)/(sqrt(rowSums(x^2) %*% t(rowSums(x^2))))) 
}

drawCluster <- function(title, method, frequencys, sam){
  points <- as.data.frame(sam$points)
  filename = paste(method, "clusters", sep=".")
  if(file.exists(filename)){
    clusters = read.table(filename)
    clusters <- clusters[frequencys$id,]
    points$class <- clusters
    points$color <- "black"
    points$color[points$class==1] <- "green"
    points$color[points$class==2] <- "red"
    points$color[points$class==3] <- "blue"
    points$color[points$class==4] <- "cyan"
    points$color[points$class==5] <- "magenta"
    points$color[points$class==6] <- "grey"

    print("Output image to output.pdf")
    pdf(file=paste(method, title, "pdf", sep="."))

    print("Plot results")
    plot(points$V1, points$V2, col=points$color)
    #text(points, labels = as.character(1:nrow(frequencys)))
    dev.off()
  }
}

applySammonMapping <- function(frequencys, distances, title){
  print("Performing matrix reduction with sammons mapping")
  sam <- sammon(distances)

  drawCluster(title, "soft-fdcl", frequencys, sam)
  drawCluster(title, "mixed-pdcl", frequencys, sam)
  drawCluster(title, "pdcl", frequencys, sam)
}

args <- commandArgs(trailingOnly = TRUE)
freq_filename <- args[1]

print("Reading data")
table = read.table(freq_filename)
table <- as.data.frame(table)
table$id <- c(1:nrow(table))

print("Checking for duplicates")
table <- table[!duplicated(table[,c(1:ncol(table)-1)]), ]

print("Computing all distances with cosine measure")
d = cosineDist(as.matrix(table))
applySammonMapping(table, d, "cosine")

print("Computing all distances with euclidian measure")
d = dist(as.matrix(table))
applySammonMapping(table, d, "euclidian")

print("Process finished")
