# This project is created by:
# Name: Aditya Bhardwaj
# Unity ID: abhardw2

library(igraph)
library(lsa)

# Reading the data from the files:
g <- read.graph(file = "/Users/ADITYA/BI/Projects/Project6/data/fb_caltech_small_edgelist.txt", format = c("edgelist"))
attrList <- read.csv("/Users/ADITYA/BI/Projects/Project6/data/fb_caltech_small_attrlist.csv", header = TRUE)

# Taking Command Line Input
args <- commandArgs(trailingOnly = TRUE)
alpha=as.numeric(args[1])

cosine_update <- function(k, memebership, values, attribute)
{
  indices=which(values==memebership)
  sim=0
  for(i in indices)
  {
    sim=sim+cosine(as.numeric(attribute[k,]),as.numeric(attribute[i,]))
  }
  
  sim <- sim/length(indices)
  
}

phase1 <- function(g, alpha){
  
  mapped_values <- 1:vcount(g)
  
  for(k in 1:15)
  {
    x = mapped_values
    for(i in 1:vcount(g))
    {
      
      index <- 0
      maxValue <- 0
      n <- neighbors(g, i)
      for(j in unique(mapped_values[n]))
      {
        community <- mapped_values
        oldModularity <- modularity(g, community)
        community[i]=j
        newModularity <- modularity(g, community)
        
        modularity_gain <- (alpha)*(newModularity - oldModularity) + (1-alpha)*(cosine_update(i, j, mapped_values, attrList))
        
        if(i!=j && modularity_gain > maxValue){
          index <- j
          maxValue <- modularity_gain
        }
      }
      if(index !=0){
        mapped_values[i] <- index
      }
      
    }
    if(isTRUE(all.equal(x,mapped_values)))
    {
      break
    }
    x <- mapped_values
    
  }
  return(mapped_values)
}


phase2 <- function(alphaValue, community_mapped){
  oldMember <- community_mapped
  for(h in 1:15)
  {
    g2 <- contract.vertices(g, oldMember)
    g3 <- simplify(g2, remove.multiple = TRUE, remove.loops = TRUE)
    community_mapped <- phase1(g3, alphaValue, community_mapped)
    if(isTRUE(all.equal(oldMember, community_mapped)))
    {
      break
    }
    oldMember <- community_mapped
  }
  return(community_mapped)
}

sac1 <- function(alphaVal){
  mapped_communities <- phase1(g,alpha=alphaVal)
  mapped_communities <- phase2(alphaValue = alphaVal, mapped_communities)
  
  # Writing Data into the files
  
  fileName<-paste("communities",alphaVal,sep="_")
  fileName<-paste(fileName,"txt",sep=".")
  fileDes<-file(fileName,"w")
  
  for(i in 1:length(unique(mapped_communities))-1)
  {
    finalComm <- vector("numeric")
    for(j in 1:vcount(g))
    {
      if(mapped_communities[j]==unique(mapped_communities)[i]){
        finalComm <- append(finalComm,j,after = length(finalComm))
      }
    }
    cat(as.character(finalComm), file=fileDes,sep = ",")
    cat("\n", file=fileDes)
  }
  close(fileDes)
  
}

sac1(alphaVal = alpha)