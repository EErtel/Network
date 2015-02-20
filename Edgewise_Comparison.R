#### Created 30 January 2015 by Ethan Ertel ####
#### Revised 20 February 2015 by Ethan Ertel ####
rm(list=ls())
source("http://bioconductor.org/biocLite.R")
biocLite()
setwd("~/Documents/Corrada/")
caseList=read.table(file="cases.tsv", sep="\t", header=TRUE)
controlList=read.table("controls.tsv", sep="\t", header=TRUE)
globalList=read.table("global.tsv", sep="\t", header=TRUE)
###########################################################################################
### calculates number of edges based on length of list
a<-length(caseList[,1])
n<-sqrt(2*a+1/4)+1/2


###########################################################################################
### splits edge name entry into two separate columns
caseList[,1] = as.character(caseList[,1])
caseList$list1 <- sapply(strsplit(caseList[,1],"-"), function(i){i[1]})
caseList$list2 <- sapply(strsplit(caseList[,1],"-"), function(i){i[2]})
#caseList$list1[a] <- caseList$list2[(a-1)]
#caseList$list1[(a-1)] <- caseList$list2[(a-2)]
caseList[,1] <- NULL
#caseList <- caseList[order(caseList$list1),]
###########################################################################################
### summarize OTU names, populate matrix with edges 
caseList[(a-5):a,]
OTUs=list()
count = 1
edges<-matrix(0,nrow = n,ncol = n)
for(i in 1:n-1){
  if(i==1){count=1}
  OTUs[i]=caseList$list1[count]
  
  for(j in (i+1):n){
    edges[i,j]<-caseList$correlation[count]
    edges[j,i]<-caseList$correlation[count]
    count = count + 1
  }
}
OTUs[n-1]<-caseList$list2[a-2] #### corner cases for unlabeled final OTU
OTUs[n]<-caseList$list2[a-1] #### corner cases for unlabeled final OTU
