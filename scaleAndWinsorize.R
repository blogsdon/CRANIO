library(synapseClient)
synapseClient::synapseLogin()

winsorize <- function(x,per=.99){
  up <- quantile(x,per,na.rm=T)
  low <- quantile(x,1-per,na.rm=T)
  x[x>=up] <- up
  x[x<=low] <- low
  return(x)
}


winsorizeAndScale <- function(x){
  library(dplyr)
  replaceNaMean <- function(x){
    if(sum(is.na(x))>0){
      y <- x
      y[is.na(x)] <- mean(x,na.rm=T)
      return(y)
    }else{
      return(x)
    }
  }
  
  x <- t(x)
  x <- apply(x,2,winsorize)
  x <- apply(x,2,replaceNaMean)
  x <- scale(x)
  return(x)
}
#######query private project and pull data
#dataType=mRNA
#foo <- synapseClient::synQuery('select name,id from file where projectId==\'syn2370594\' and dataType==\'mRNA\'')
#foo <- foo[1:7,]
#foo <- synapseClient::synTableQuery("SELECT * FROM syn8681664 where dataType = 'mRNA' and assay = 'RNAseq' and columnScaled = 'FALSE'")@values
#get data for just rnaseq studies
foo <- c('syn8555302','syn8545679')
bar <- lapply(foo,synGet)

#function to load data quickly
loadData <- function(x){
  library(data.table)
  testData <- read.table(x@filePath,header=T,row.names=1)
  #testData <- data.table::fread(x@filePath,data.table=F)
  #testData <- t(testData)
  return(testData)
}

#load all 7 data sets into a lis
listOfData <- lapply(bar,loadData)


listOfData2 <- lapply(listOfData, function(x){
  library(dplyr)
  foob <- apply(x,1,sd,na.rm=T)
  naTest <- foob %>%
    is.na %>%
    which
  if(length(naTest)>0){
    x <- x[-naTest,]
  }
  return(x)
})

#winsorize and scale
reformattedData <- lapply(listOfData2,winsorizeAndScale)

library(dplyr)

#fileNames <- paste0('Scaled_Winsorized_',foo$name)
fileNames <- sapply(bar,function(x){return(paste0('scaled and winsorized ',x@properties$name))})
parentIds <- sapply(bar,function(x){return(x@properties$parentId)})
annotations <- sapply(bar,function(x){return(synGetAnnotations(x))})
comments <- rep('Winsorizing and scaling',2)
useds <- sapply(bar,function(x){return(x@properties$id)})

library(githubr)

permLink1 =githubr::getPermlink(repository = 'blogsdon/CRANIO',
                                ref = 'branch',
                                refName = 'master',
                                repositoryPath = 'scaleAndWinsorize.R')

permLink2 =githubr::getPermlink(repository = 'Sage-Bionetworks/rSynapseUtilities',
                                ref = 'branch',
                                refName = 'dev',
                                repositoryPath = 'R/pushToSynapseWrapper.R')

executeds <- vector('list',2)
for(i in 1:2){
  executeds[[i]]<-as.list(c(permLink1,permLink2))
}

activityNames <- rep('Post Process CRANIO RNAseq data',2)

activityDescriptions <- rep('Transpose, Winsorize, set Nas to mean, then Scale each data-set',2)

#push back to synapse
synObjs<-mapply(rSynapseUtilities::pushToSynapseWrapper,
                df=reformattedData,
                fileName=fileNames,
                synapseFolderId=parentIds,
                annos=annotations,
                comment=comments,
                usedVector=useds,
                executedVector=executeds,
                activityName1=activityNames,
                activityDescription1=activityDescriptions)