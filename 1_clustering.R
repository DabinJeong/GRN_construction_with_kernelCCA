#!/usr/bin/env python
# coding: utf-8
# # Multilevel clustering _ R
library(igraph)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
file_name = args[1] 
out_file = args[2] 
out_file2 = args[3]
## Data import and cleaning
print('START')
parsed_genes = read.csv(file_name, sep='\t',header=TRUE,stringsAsFactors = FALSE)

## Creating igraph object

edgelist = as.matrix(parsed_genes[,1:2])
g = graph.edgelist(edgelist, directed=FALSE)

## Community detection : Multilevel algorithm
mlc <- multilevel.community(g)

## Writing the results to a file

mlc_community_list <- as.data.frame(as.matrix(membership(mlc)))
mlc_community_list$gene <- rownames(mlc_community_list)
mlc_community_member <- arrange(mlc_community_list, V1) %>%
  select(gene, V1)
colnames(mlc_community_member)[2] <- 'community'

community_info = mlc_community_member

# Lookup Table
vals <- community_info[,2]
keys <- community_info[,1]
lookup <- setNames(vals, keys)

index = which(lookup[parsed_genes[,1]] != lookup[parsed_genes[,2]])
filtered_edgelist <- parsed_genes[-index,]

write.table(mlc_community_member, row.names=FALSE, col.names=TRUE, file=out_file, sep='\t', quote=FALSE)
write.table(filtered_edgelist, row.names=FALSE, col.names=TRUE, file=paste(out_file2,sep=''), sep='\t',quote=FALSE)

