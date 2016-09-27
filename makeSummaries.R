library(synapseClient)
synapseLogin()

rankNetwork <- list()

modulesObj <- synGet('syn5700963')
enrichmentsObj <- synGet('syn5701314')
networkObj <- synGet('syn5700531')

library(data.table)
library(dplyr)
load(networkObj@filePath)
rankNetwork$network <- sparseNetwork
rankNetwork$modules <- data.frame(fread(modulesObj@filePath,stringsAsFactors = F))
rankNetwork$enrichments <- data.frame(fread(enrichmentsObj@filePath,stringsAsFactors=F))
rankNetwork$enrichments <- arrange(rankNetwork$enrichments,fdr)


