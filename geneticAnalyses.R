#grab genetic data
#pkg <- c("package:IRanges","package:XML","package:digest","package:bitops","package:DBI","package:stats4","package:RSQLite","package:S4Vectors","package:rjson","package:tools","package:Biobase","package:RCurl","package:parallel","package:BiocGenerics","package:AnnotationDbi","package:biomaRt")
#lapply(pkg, detach, character.only = TRUE, unload = TRUE,force=TRUE)


library(synapseClient)
synapseLogin()

foo <- synQuery('select name,id from file where parentId==\'syn5752525\'')

bar <- lapply(foo$file.id,synGet)
library(data.table)
dataMatrices <- lapply(bar,function(x){return(read.delim(x@filePath,stringsAsFactors=F,sep='\t'))})

esp6500 <- dataMatrices[[1]]
cases <- strsplit(esp6500$cases,'_')
cases2 <- strsplit(esp6500$cases,'_')
cases <- c(unlist(cases))
uniqueNames <- unique(cases)

esp6500individual <- matrix(0,length(uniqueNames),nrow(esp6500))
rownames(esp6500individual) <- uniqueNames
colnames(esp6500individual) <- paste0('chr',esp6500$Chr,':',esp6500$Pos)

for (i in 1:ncol(esp6500individual)){
  specificCases <- cases2[[i]]
  if(esp6500$Zygosity[i]=='het'){
    esp6500individual[specificCases,i]<-1
  }else{
    esp6500individual[specificCases,i]<-2
  }
}

w1 <- which(p.adjust(esp6500$p.value,method = 'BH')<=0.05)
geneList <- esp6500$Gene.refGene[w1]

cranioModsObj <- synGet('syn5700963')
cranioMods <- read.delim(cranioModsObj@filePath,stringsAsFactors=F)
ensembl=biomaRt::useMart('ENSEMBL_MART_ENSEMBL',
                         dataset = 'hsapiens_gene_ensembl',
                         host='www.ensembl.org')
library(biomaRt)
genes<-getBM(attributes = c('ensembl_gene_id','external_gene_name'),
             filters='ensembl_gene_id',
             values=cranioMods$GeneIDs,
             mart=ensembl)

#genes<-getBM(attributes = c('ensembl_gene_id','external_gene_name'),
#             filters='ensembl_gene_id',
#             values=n1,
#             mart=ensembl)

cranioMods2 <- merge(cranioMods,
                     genes,
                     by.x='GeneIDs',
                     by.y='ensembl_gene_id')

keep20 <-names(which(table(cranioMods2$modulelabels)>20))
cranioMods2 <- dplyr::filter(cranioMods2,modulelabels%in%keep20)

modList <- sapply(unique(cranioMods2$modulelabels),
                  utilityFunctions::listify,
                  cranioMods2$external_gene_name,
                  cranioMods2$modulelabels)


names(modList) <- unique(cranioMods2$modulelabels)

pvals<-sapply(modList,utilityFunctions::fisherWrapperPval,geneList,genes$external_gene_name)
ors <- sapply(modList,utilityFunctions::fisherWrapperOR,geneList,genes$external_gene_name)

exprMat <- read.csv('cranioRNAseq.csv',stringsAsFactors=F)



#################################SKAT ANALYSES


####merge expression sets and genotype sets







greenGene <- filter(cranioMods2,modulelabels=='green')$GeneIDs
greenCranio <- exprMat[,greenGene]
baz <- hclust(dist(t(scale(greenCranio))))
fob <- cutree(baz,4)
newMatrix <- data.frame(cbind(names(fob),fob),stringsAsFactors=F)
colnames(newMatrix) <- c('ensemblId','subModule')
newMatrix <- merge(newMatrix,genes,by.x='ensemblId',by.y='ensembl_gene_id',stringsAsFactors=F)

greenModList <- sapply(unique(newMatrix$subModule),utilityFunctions::listify,newMatrix$external_gene_name,newMatrix$subModule)

pvals2<-sapply(greenModList,utilityFunctions::fisherWrapperPval,geneList,genes$external_gene_name)
ors2 <- sapply(greenModList,utilityFunctions::fisherWrapperOR,geneList,genes$external_gene_name)


greenModAnnoObj <- synGet('syn5758862')
greenModAnno <- read.delim(greenModAnnoObj@filePath,sep='\t',stringsAsFactors=F)

genesInGeneList <- genes$external_gene_name%in%geneList
genes2 <- cbind(genes,genesInGeneList)
genes2 <- filter(genes2,ensembl_gene_id%in%greenGene)
write.csv(genes2,file='greenanno2.csv',quote=F)
write.csv(newMatrix,file='submodules.csv',quote=F)
