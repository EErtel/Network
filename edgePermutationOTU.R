#### Revised 09 March 2015 by Ethan Ertel ####
network_vector <- c("~/Documents/Corrada/cases.tsv","~/Documents/Corrada/controls.tsv")
label_vector <- c("cases","controls")
perms <- 10
#inputs <- commandArgs(trailingOnly = TRUE)
#network_vector <- c(inputs[1],inputs[2])
#label_vector <- c(inputs[3],inputs[4])
#perms <- as.integer(inputs[5])

### Calculates t-statistic to compare two distributions
### ??? The result differs the t-statistic reported in R's base package by a factor of sqrt(2) ???
t_statistic <- function(distribution1, distribution2){  # NA values will result in error
  return(  (mean(distribution1)-mean(distribution2)) / ( sqrt( var(distribution1) / ((length(distribution1))) ) + sqrt( var(distribution2) / ((length(distribution2))) ) )  )
}

###########################################################################################

edgesByOTU <- function(edgeList){  ### parses adjacency list, outputs n*n matrix (n*1 OTU names, n*(n-1) edges)
  ### calculates number of edges based on length of list
  a<-length(edgeList[,1])
  n<-sqrt(2*a+1/4)+1/2
  ### splits edge name entry into two separate columns
  edgeList[,1] = as.character(edgeList[,1])
  edgeList$list1 <- sapply(strsplit(edgeList[,1],"-"), function(i){i[1]})
  edgeList$list2 <- sapply(strsplit(edgeList[,1],"-"), function(i){i[2]})
  edgeList[,1] <- NULL
  ### summarize OTU names, populate matrix with edges 
  OTU_List=c(rep(0,n))
  count = 1
  edges<-matrix(0,nrow = n,ncol = n)
  for(i in 1:n-1){
    if(i==1){count=1}
    OTU_List[i]=edgeList$list1[count]
    
    for(j in (i+1):n){
      edges[i,j]<-edgeList$correlation[count]
      edges[j,i]<-edgeList$correlation[count]
      count = count + 1
    }
  }
  OTU_List[n-1]<-edgeList$list2[a-2] #### corner cases for unlabeled final OTU
  OTU_List[n]<-edgeList$list2[a-1] #### corner cases for unlabeled final OTU
  ### creates matrix - edge for each OTU
  for(i in 1:n){
    if(i==1){edgeMat <- matrix(edges[1,c(2:n)],nrow = 1,ncol=n-1)}
    else if(i<n){edgeMat <- rbind(edgeMat,edges[i,c(1:(i-1),(i+1):n)])}
    else{edgeMat <- rbind(edgeMat,edges[n,c(1:n-1)])}  
  }
  return(cbind(as.integer(OTU_List), edgeMat))
}

###########################################################################################

edgePermutationOTU <- function(network_vector, labels, permutations = 1000){
  if(length(network_vector) != 2){stop("Input vector must contain two networks")}
  if(length(labels) != 2){stop("Label vector must contain two labels")}
  len <- length(network_vector)
  networks <- list()
  ### Reads network files to generate matrix of edges for each OTU
  for(i in 1:len){networks[[i]] <- edgesByOTU(read.table(file=network_vector[i], sep="\t", header=TRUE))}
  n <- length(networks[[1]][,1])
  OTU.labels <- networks[[1]][,1]
  permutation_results <- rep(0,n)
  for(i in 1:n){
    population1 <- networks[[1]][i,2:n] #Edges for OTU i; network 1
    population2 <- networks[[2]][i,2:n] #Edges for OTU i; network 2
    notNA <- (!is.na(population1))&(!is.na(population2))
    tStatistic1 <- t_statistic(population1[notNA],population2[notNA])
    pool <- c(population1[notNA],population2[notNA])  ### pool of all edges
    permutation_t_statistics <- rep(0,permutations)
    nGood = sum(notNA)
    for(t in 1:permutations){
        shuffled_labels <- sample(1:(2*nGood), (2*nGood), replace = FALSE, prob = NULL)
        population1 <- pool[shuffled_labels[1:nGood]]
        population2 <- pool[shuffled_labels[(nGood+1):(2*nGood)]]
        permutation_t_statistics[t] <- t_statistic(population1,population2)
    }
    #permutationPvalues <- rep(0,permutations)
    permutation_results[i] <- mean(permutation_t_statistics>tStatistic1)
  }
  return(list(data.frame(OTU.labels,permutation_results),labels))
}
edgePermutationOTU(network_vector = network_vector,labels=label_vector,permutations)
