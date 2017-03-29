# This project is created by:
# Name: Aditya Bhardwaj
# Unity ID: abhardw2

library(igraph)
library(lsa)

graph=read_graph("/Users/ADITYA/BI/Projects/Project6/data/fb_caltech_small_edgelist.txt",format= c("edgelist"))
attribute_data <- read.csv("/Users/ADITYA/BI/Projects/Project6/data/fb_caltech_small_attrlist.csv",header = TRUE)

setwd("/Users/ADITYA/BI/Projects/Project6/")

# Cosine function for finding the modularity
cosine_update <- function(attrList, h, membership, values)
{
  indices <- which(values == membership)
  similar <- 0
  for(i in indices)
  {
    similar <- similar + cosine(as.numeric(attrList[h,]), as.numeric(attrList[i,]))
  }
  similar <- similar/length(indices)
}


#Phase 1 of Sac1 Algorithm

phase1 <- function(graph, mapped_values, alpha, attributes){
  for(h in 1:15)
  {
    x <- mapped_values
    for(i in 1:vcount(graph))
    {
      
      index <- 0
      max <- 0
      
      n <- neighbors(graph, i)
      for(j in unique(mapped_values[n]))
      {
        
        membership_prime <- mapped_values
        oldMod <- modularity(graph,membership_prime)
        membership_prime[i] <- j
        newMod <- modularity(graph,membership_prime)
        
        delta_Q_newman <- newMod - oldMod
        delta_Q_attr <- cosine_update(attributes, i, j, mapped_values)
        delta_Q <- (1-alpha)*delta_Q_attr + (alpha)*delta_Q_newman
        
        if(i!=j && delta_Q > max){
          index <- j
          max <- delta_Q
        }
      }
      if(index !=0){
        mapped_values[i] <- index
      }
      
    }
    if(isTRUE(all.equal(x, mapped_values)))
    {
      break
    }
    x <- mapped_values
    
  }
  mapped_values  
}

#Phase2 of sac1 algorithm
phase2 <- function(graph, mapped_values, alpha, attributes){
  
  x <- mapped_values
  for(i in 1:15)
  {
    graph_prime <- contract.vertices(graph, mapped_values)
    
    graph_prime_1 <- simplify(graph_prime, remove.multiple = TRUE, remove.loops = TRUE)
    
    mapped_values <- phase1(graph_prime_1, mapped_values, alpha, attributes)
    
    if(isTRUE(all.equal(x, mapped_values)))
    {
      break
    }
    
    x <- mapped_values
  }
  
  mapped_values
}


#By changing the value of alpha we can get different communities
sac1 <- function(alpha, attributes = attribute_data){
  mapped_communities <- phase1(graph, alpha=alpha, mapped_values = c(1:324), attributes)
  
  mapped_communities <- phase2(graph, alpha=alpha, mapped_values = mapped_communities, attributes)
  
  return(mapped_communities)
}

fileWrite <- function(mappedCommunity, alpha){
  
  fileName<-paste("communities",alpha,sep="_")
  fileName<-paste(fileName,"txt",sep=".")
  fileConnection<-file(fileName,"w")
  
  for(i in 1:length(unique(mappedCommunity)))
  {
    community <- vector("numeric")
    for(j in 1:324)
    {
      if(mappedCommunity[j]==unique(mappedCommunity)[i]){
        
        community <- append(community, j-1, after = length(community))
        
      }
    }
    cat(as.character(community), file=fileConnection, sep = ",")
    cat("\n", file=fileConnection)
  }
  
  close(fileConnection)
  
}


args <- commandArgs(trailingOnly = TRUE)
alpha = as.numeric(args[1])

# Running SAC1 Algorithm
mappedComunity <- sac1(alpha = alpha)
fileWrite(mappedComunity, alpha = alpha)

# Discussed with Anshuman Goel and Ayush Kumar for project
