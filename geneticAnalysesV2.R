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

#exprMatObj = synGet('syn7300582')
#exprMat = read.delim(exprMatObj@filePath,stringsAsFactors=F)


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





####first type of filter: grab all variants for a set of genes.
library(dplyr)
getGeneSetVariantDataFrame = function(geneSet,genotypeMatrix,variantAnnotation){
  variants = dplyr::filter(variantAnnotation, Gene.refGene%in%geneSet) %>%
    dplyr::select(genoIdentifiers)
  return(variants)
}


#####analyses to be run
###1) Control vs CASE
###2) gene (snp) to gene via SKAT-O at CADD filters
###3) gene (snp) x mod to gene via SKAT-O at CADD filters
###4) gene (snp) x mod to gene x mod via SKAT-O at CADD filters


###cadd 
caddRaw = foo2$CADD_raw %>%as.numeric
caddThresholds = quantile(caddRaw,na.rm=T,c(.8,.9,.95,.99))

fxn1 = function(threshold,
                genotypeMatrix,
                annotationMatrix){
  
  variants = dplyr::filter(annotationMatrix, CADD_raw>threshold) %>%
    dplyr::select(genoIdentifiers)
  return(genotypeMatrix[,variants$genoIdentifiers])
  
}

filteredGenotypes = lapply(caddThresholds,fxn1,genotypeMatrix,foo2)

###gene by gene analysis






foo3 = dplyr::filter(foo2,CADD_raw > 10)
foobar2=getGeneSetVariantDataFrame(test1,genotypeMatrix,foo3)


bar10=genotypeMatrix[,foobar2$genoIdentifiers]
bar23 = rowSums(bar10)
test1=modList[[1]]

map234 = genes$ensembl_gene_id
names(map234) = genes$external_gene_name
greenModExpr = exprMat[,map234[test1]]
eigenGenes = svd(scale(greenModExpr))
summary(lm(bar23 ~ eigenGenes$u[,1:4]))
