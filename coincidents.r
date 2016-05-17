library(MASS)
cosineDist <- function(x){
  as.dist(1 - x%*%t(x)/(sqrt(rowSums(x^2) %*% t(rowSums(x^2))))) 
}

if(file.exists("prototypes.matrix")){
	print("Evaluate coincident prototypes")
	tip <- as.matrix(read.table("prototypes.matrix"))
	dist(tip)
}
if(file.exists("tipicalities.matrix")){
	print("Evaluate coincident tipicalities")
	tip <- as.matrix(read.table("tipicalities.matrix"))
	dist(t(tip))
}
if(file.exists("memberships.matrix")){
	print("Evaluate coincident memberships")
	tip <- as.matrix(read.table("memberships.matrix"))
	dist(t(tip))
}

