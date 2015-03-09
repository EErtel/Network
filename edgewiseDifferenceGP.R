#### Revised 09 March 2015 by Ethan Ertel ####
rm(list=ls())

network_vector <- c("~/Documents/Corrada/cases.tsv","~/Documents/Corrada/controls.tsv")
label_vector <- c("cases","controls")
perms <- 10
#inputs <- commandArgs(trailingOnly = TRUE)
#network_vector <- c(inputs[1],inputs[2])
#label_vector <- c(inputs[3],inputs[4])
#perms <- as.integer(inputs[5])

### Global Permutation (GP) of networks - samples two new networks from the pooled set 
### of all edges in the two observation-based networks
EdgewiseDifferenceGP <- function(network_vector, labels, permutations = 100){
  networks <- list()
  for(i in 1:2){networks[[i]] <- read.table(file=network_vector[i], sep="\t", header=TRUE)}
  ### Removes edges that are NA in either observation-based network
  goodRows<-as.logical((!is.na(networks[[1]]$correlation))&(!is.na(networks[[2]]$correlation))) 
  
  edgeLabels <- networks[[1]]$Taxa.and.Samples[goodRows]
  population1 <- networks[[1]]$correlation[goodRows]
  population2 <- networks[[2]]$correlation[goodRows]
  n <- length(population1)
  edgeDifferences <- c(rep(0,n))
  for(idx in 1:n){edgeDifferences[idx] = population2[idx] - population1[idx]}
  pool = c(population1,population2)
  difference_vector <- rep(0,n)  
  for(i in 1:permutations){
    shuffled_labels <- sample(1:(2*n), (2*n), replace = FALSE, prob = NULL)
    population1 <- pool[shuffled_labels[1:n]]
    population2 <- pool[shuffled_labels[(n+1):(2*n)]]
    for(idx in 1:n){difference_vector[idx] <- population1[idx]-population2[idx]}
    edgeDifferences <- rbind(edgeDifferences, difference_vector)
  }
  ### compares the absolute value of calculated difference in observed edge weights with permuted edge weights
  perm_results <- colMeans(abs(edgeDifferences[2:(permutations+1),])<abs(edgeDifferences[1,]))
  return(data.frame(Edge_Labels=edgeLabels,permutationStatistics=perm_results))
}
EdgewiseDifferenceGP(network_vector,label_vector,permutations)
