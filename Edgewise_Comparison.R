#### Created 30 January 2015 by Ethan Ertel ####
#### Revised 23 February 2015 by Ethan Ertel ####
rm(list=ls())

network_vector <- c("~/Documents/Corrada/cases.tsv","~/Documents/Corrada/controls.tsv")
#network_vector <- for(thing in commandArgs(trailingOnly = TRUE)){}

###########################################################################################
### parses adjacency list, outputs n*n matrix (n*1 OTU names, n*(n-1) edges)
edgesByOTU <- function(edgeList){
  ###########################################################################################
  ### calculates number of edges based on length of list
  a<-length(edgeList[,1])
  n<-sqrt(2*a+1/4)+1/2
  ###########################################################################################
  ### splits edge name entry into two separate columns
  edgeList[,1] = as.character(edgeList[,1])
  edgeList$list1 <- sapply(strsplit(edgeList[,1],"-"), function(i){i[1]})
  edgeList$list2 <- sapply(strsplit(edgeList[,1],"-"), function(i){i[2]})
  #caseList$list1[a] <- caseList$list2[(a-1)]
  #caseList$list1[(a-1)] <- caseList$list2[(a-2)]
  edgeList[,1] <- NULL
  #caseList <- caseList[order(caseList$list1),]
  ###########################################################################################
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
  
  ###########################################################################################
  ### creates matrix - edge for each OTU
  for(i in 1:n){
    if(i==1){edgeMat <- matrix(edges[1,c(2:n)],nrow = 1,ncol=n-1)}
    else if(i<n){edgeMat <- rbind(edgeMat,edges[i,c(1:(i-1),(i+1):n)])}
    else{edgeMat <- rbind(edgeMat,edges[n,c(1:n-1)])}  
  }
  return(cbind(as.integer(OTU_List), edgeMat))
}
###########################################################################################

#### This script takes as inputs two network file locations (weighted adjacency lists of coexpressed OTUs) and a label vector (phenotypic state for each network);
#### optional inputs are (or, may be) 
#### the function and outputs a list of edge significance (by pt)
edgeSignificance <- function(network_vector, labels, permutations = 10){
  if(length(network_vector) != 2){stop("Input vector must contain two networks")}
  if(length(labels) != 2){stop("Label vector must contain two labels")}
  len <- length(network_vector)
  networks <- list()
  ### Reads network files to generate matrix of edges for each OTU
  for(i in 1:len){ 
    networks[[i]] <- edgesByOTU(read.table(file=network_vector[i], sep="\t", header=TRUE))
  }
  #if(sum(apply(X=cbind(a[[1]][,1],a[[2]][,1]),MARGIN=1,FUN=function(i){i[[1]]!=i[[2]]}))!=0){stop("Networks must consist of same OTUs")} ### Compares OTU names to confirm networks are comparable
  n <- length(a[[1]][,1])
  pval_results <- matrix(ncol = 2,nrow=n*(n-1)/2) #matrix(nrow=n,ncol=2)
  permutationPvalues <- rep(0,permutations)
  count = 0
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      dist1 <- c(networks[[1]][i,2:n],networks[[1]][j,2:n]) #Pooled edges for OTUs i,j; network 1
      dist2 <- c(networks[[2]][i,2:n],networks[[2]][j,2:n]) #Pooled edges for OTUs i,j; network 2
      ### edge i,j = networks[[1]][i,j] will be present in duplicate; is this a problem?
      pvalue1 <- t.test(dist1,dist2, alternative = "two.sided")[[3]]
      pool <- c(dist1,dist2)  ### pool of all edges
      for(t in 1:permutations){
        samps <- sample(1:(4*(n-1)), (4*(n-1)), replace = FALSE, prob = NULL)
        dist1 <- pool[samps[1:(2*(n-1))]]
        dist2 <- pool[samps[(2*(n-1)+1):(4*(n-1))]]
        permutationPvalues[t]<-t.test(dist1,dist2, alternative = "two.sided")[[3]]
        ### t.test scales (about) linearly?
      }
      count = count + 1
      pval_results[count,1] <- paste(as.character(i),as.character(j),sep='-')
      pval_results[count,2] <- as.numeric(sum(permutationPvalues<=pvalue1))/permutations
      ### This is a raw value; Holm's step-wise method may be more appropriate
    }
  }
  return(pval_results)
}
