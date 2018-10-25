###   5/23/17 Script to Compare GCN clusters of Mok networks (correlation based)

###   Set working directory
#setwd("~/ferdig_rotation/pB_data/my_results/DTWMIC/rewired_nets/rewiring_code/")
library(nettools)
library(igraph)

################ read in data, mix up, make networks ################
PB_data <- read.delim("Gene_Data_BC_all.csv",header=TRUE,sep=",",as.is=TRUE)

### Get all NF54 data and PB-58 data
NF54_data <- PB_data[which(PB_data[,"pBLine"]=="NF54"),]
PB58_data <- PB_data[which(PB_data[,"pBLine"]=="PB-58"),]

### Make datasets with replicates averaged together
NF54_sum <- matrix(NA,ncol=5540,nrow=12)
colnames(NF54_sum) <- colnames(NF54_data)[6:5545]
for(i in 6:length(NF54_data[1,])){
  NF54_sum[,i-5] <- aggregate(as.formula(paste(paste(colnames(NF54_data)[i],collapse=""),"~ hpi")),data=NF54_data,mean)[,2]
}

PB58_sum <- matrix(NA,ncol=5540,nrow=12)
colnames(PB58_sum) <- colnames(PB58_data)[6:5545]
for(i in 6:length(PB58_data[1,])){
  PB58_sum[,i-5] <- aggregate(as.formula(paste(paste(colnames(PB58_data)[i],collapse=""),"~ hpi")),data=PB58_data,mean)[,2]
}

#scamble data, for each column swap all of the values
NF54_permmat <- NF54_sum
PB58_permmat <- PB58_sum
for(ii in 1:ncol(NF54_sum)){
  # for each column, swap all the values
  NF54_permmat[,ii] <- sample(NF54_sum[,ii])
  PB58_permmat[,ii] <- sample(PB58_sum[,ii])
}

#make the network
NF54_network <- mat2adj(NF54_permmat,method="DTWMIC")
PB58_network <- mat2adj(PB58_permmat,method="DTWMIC")

#write.csv(PB57_AM,"PB57_DTQMIC_mat.csv")


############### rewiring code ####################

###   Read in DTWMIC matrix 
#NF54_network <- read.delim("NF54_DTWMIC_mat.csv",sep=",",header=TRUE,as.is=TRUE,row.names=1)
#PB58_network <- read.delim("PB58_DTQMIC_mat.csv",sep=",",header=TRUE,as.is=TRUE,row.names=1)


###   Turn DTWMIC into edge lists
### Function flattenSquareMatrix
### Input: Square matrix with p-value in top diagonal and cor in bottom diagonal
### Output: dataframe with node pairs in columns 1 and 2, cor in column 3 and p-value in column 4
flattenSquareMatrix <- function(m) {
  if( (class(m) != "matrix") | (nrow(m) != ncol(m))) stop("Must be a square matrix.") 
  if(!identical(rownames(m), colnames(m))) stop("Row and column names must be equal.")
  ut <- upper.tri(m)
  data.frame(i = rownames(m)[row(m)[ut]],
             j = rownames(m)[col(m)[ut]],
             cor=t(m)[ut])
}

NF54_EL <- flattenSquareMatrix(data.matrix(NF54_network))
PB58_EL <- flattenSquareMatrix(data.matrix(PB58_network))

colnames(NF54_EL) <- c("k","j","cor")
colnames(PB58_EL) <- c("k","j","cor")

###   Rewiring of 2 matrixes
rewiring_score <- function(i) {
  cor_r <- NF54_EL[i,"cor"]
  cor_s <- PB58_EL[i,"cor"]
  Fr <- (1/2)*log((1+cor_r)/(1-cor_r))
  Fs <- (1/2)*log((1+cor_s)/(1-cor_s))
  score <- abs(Fr-Fs)/sqrt((1/(12-3))+(1/(12-3)))
  p <- pnorm(score,lower.tail=FALSE)
  #if(i%%100000==0){
  #  write.csv(i,paste("out",i,".csv"))
  #}
  return(c(levels(NF54_EL$k[i])[as.numeric(NF54_EL$k[i])],levels(NF54_EL$j[i])[as.numeric(NF54_EL$j[i])],cor_r,cor_s,score,p))
}

require(parallel)
cores=4
#15343030 p*(p-1)/2 #i:number of rows in edgelist
PB58_rewiring <- matrix(unlist(mclapply(1:nrow(NF54_EL), FUN=rewiring_score, mc.cores=cores)),ncol=6,byrow=TRUE)
colnames(PB58_rewiring) <- c("Gene_i","Gene_j","cor_i","cor_j","score","p")
adj.p <- p.adjust(PB58_rewiring[,"p"],method="fdr")

#write.csv(cbind(PB54_rewiring,adj.p),"PB58_rewiring_sig.csv")
PB58_rewiring <- cbind(PB58_rewiring,adj.p)
PB58_rewired <- PB58_rewiring[which(PB58_rewiring[,"p"]<0.05),]
#write.csv(PB58_rewired,"PB58_rewiring_p0.05.csv")

#get degree of nodes
PB58_rew_graph <- graph_from_edgelist(PB58_rewired[,1:2])
PB58_rew_degree <- degree(PB58_rew_graph)
degree_mat <- as.matrix(PB58_rew_degree)
#write.table(degree_mat, file = "PB58_rewire_perm_degrees.csv", append=T, sep=",", row.names = T, col.names = F)
lgd <- format(Sys.time(), "%Y_%m_%d_%H_%M_%S")
write.table(degree_mat, file = paste("degree_output_files/PB58_rw_deg", lgd, ".csv", sep=''), append=T, sep=",", row.names = T, col.names = F)

rm(list=ls())
gc()
quit()
n

