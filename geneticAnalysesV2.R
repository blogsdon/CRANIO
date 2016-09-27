require(synapseClient)
synapseLogin()


foo <- read.csv('genotypes.csv',stringsAsFactors=F,row.names=NULL)



#############reformat into genotype matrix + annotations matrix

table(foo$Proband)
#identify unique samples
#identify unique snps
#populate genotype matrix
#populate genotype annotaiton data.frame
#push to synapse

uniqueSamples <- unique(foo$Proband)

buildGenoIdentifier = function(x,y,z){
  if(x=='.'){
    return(paste0('chr',y,':',z))
  }else{
    return(x)
  }
}

genoIdentifiers <- mapply(buildGenoIdentifier,foo$snp138,foo$Chr,foo$Start)

foo$genoIdentifiers <- genoIdentifiers
uniqueGenotypes <- unique(genoIdentifiers)

genotypeMatrix <- matrix(0,length(uniqueSamples),length(uniqueGenotypes))
rownames(genotypeMatrix) <- uniqueSamples
colnames(genotypeMatrix) <- uniqueGenotypes
for(i in 1:nrow(foo)){
  if(foo$Zygosity[i]=='het'){
    genotypeMatrix[as.character(foo$Proband[i]),foo$genoIdentifiers[i]] <- 1
  }else{
    genotypeMatrix[as.character(foo$Proband[i]),foo$genoIdentifiers[i]] <- 2
  }
}

pairs(svd(scale(genotypeMatrix))$u[,1:5])

###make annotation matrix

#which(duplicated(foo$genoIdentifiers))[1]
#View(dplyr::filter(foo,genoIdentifiers==foo$genoIdentifiers[345]))

foo2 <- dplyr::filter(foo,!duplicated(genoIdentifiers))
foo2 <- dplyr::select(foo2,-Proband,-Info)
View(foo2)


#download gene expression data
exprMat <- read.csv('cranioRNAseq.csv',stringsAsFactors=F)




#mappingFileObj <- synGet('syn2823605')
mappingFile <- synTableQuery('SELECT * FROM syn2823605')
fooBar <- mappingFile@values


fooBar2 <- dplyr::select(fooBar,`Px Code`, SeqSampleName)

code1 <- fooBar2[,1]
names(code1) <- fooBar2[,2]

genoToExpr <- rownames(genotypeMatrix)

code2 <- code1[genoToExpr]
rownames(genotypeMatrix) <- code2


rownames(exprMat) = exprMat$X
exprMat = exprMat[,-1]

idsToKeep <- intersect(rownames(genotypeMatrix),rownames(exprMat))
genotypeMatrix <- genotypeMatrix[idsToKeep,]
exprMat <- exprMat[idsToKeep,]

library(SKAT)

######write utility funciton to extract genotypes for a gene set and a set of filters
######
