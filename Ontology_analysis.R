
############################################################################
# This program performs Deseq2 analysis and ontology analysis on count     #
# data. To use this script you will need count data from an RNA-seq        #
# experiment                                                               #
############################################################################


setwd("~/Dropbox/Documents/Rfiles/ZCC_STAR/")

# First establish some background data sets/functions

source("http://bioconductor.org/biocLite.R")
biocLite()

library( DESeq2 )     # The DESeq tools
library(org.Hs.eg.db) # gene annotation file for human genes
library(org.Rn.eg.db) # gene annotation file for rat
library(reactome.db)  # CSH based ontology file; this is the default for this program
library(GO.db)        # GO consortium ontology file; you can also use this one if desired
library(Hmisc)

load("ZCCcounts.RData")


#Load a data set, with the following variables:
# countdata : a table with the read counts, with technical replicates summed up,
# coldata : a table with metadata on the count table’s columns, i.e., on the samples,
# rowdata : a table with metadata on the count table’s rows, i.e., on the genes

dds <- DESeqDataSetFromMatrix(
  countData = countdata,
  colData = coldata,
  design = ~ condition)

dds = DESeq(dds)
res = results(dds)

# Define a count threshold, below which we throw out values
filterThreshold = 2
keep = rowMeans( counts( dds, normalized=TRUE ) ) > filterThreshold
table( keep )

# Now filter the results so that these low-count genes are eliminated
resSig = res[keep, ] 

# Now re-calculate the adjusted p-value and order based on logfold change
resSig$padj = p.adjust( resSig$pvalue, method="BH" ) 
resSig$bonf = p.adjust( resSig$pvalue, method="bonferroni" ) 
resSig = resSig[ order( resSig$log2FoldChange ), ]



##########################################################################################
# The next part of code calculates ontology categories based on p-values and fold change #
##########################################################################################


catagory_filter = 5 # we will ignore ontology catagories with less than this number of genes

# If you need to convert one set of genes identifiers into another set of gene identifiers, 
# this is a subfunction that allows conversion of keys from the gene annotation file
convertIDs = function( ids, fromKey, toKey, db, ifMultiple=c( "putNA", "useFirst" ) ) {
  stopifnot( inherits( db, "AnnotationDb" ) )
  ifMultiple = match.arg( ifMultiple )
  selRes = AnnotationDbi::select( db, keys=ids, keytype=fromKey, columns=c(fromKey,toKey) )
  if( ifMultiple == "putNA" ) {
    duplicatedIds = selRes[ duplicated( selRes[,1] ), 1 ]
    selRes = selRes[ ! selRes[,1] %in% duplicatedIds, ] }
  return( selRes[ match( ids, selRes[,1] ), 2 ] )
}


# first, put your genenames into a column for ease of reference
resSig$symbol = rownames(resSig)

# reactome.db works with ENTREZIDs; convert if needed (sample code below)
resSig$entrez = convertIDs( resSig$symbol, "SYMBOL", "ENTREZID", org.Rn.eg.db )

# Now filter results so that you only have genes that are in the reactome.db; also filter for non-NA p-values
res2 = resSig[ resSig$entrez %in% keys( reactome.db, "ENTREZID" ) & !is.na( resSig$pvalue) , ]

# Create a reactomeTable that matches ENTREZIDs to REACTOMEIDs
reactomeTable = AnnotationDbi::select( reactome.db, keys=res2$entrez, keytype="ENTREZID", columns=c("ENTREZID","REACTOMEID") )

# Create a boolean matrix (incm) with one row for each REACTOMEID and one column for ENTREZID, 1 if gene in reactome path, 0 if not
incm = do.call( rbind, with(reactomeTable, tapply( ENTREZID, factor(REACTOMEID), function(x) res2$entrez %in% x ) ))
colnames(incm) <- res2$entrez

# Remove ontology groups with less than the "catagory_filter" number of genes (see above)
incm = incm[ rowSums(incm) >= catagory_filter, ]

# Finally, define a function that will assemble a data frame of reactomeIDs that gives a p-value for 
# whether their constituents are statistically changed from 0.

testCategory = function( reactomeID ) {
  isMember = incm[ reactomeID, ]
  reactomeName=reactomePATHID2NAME[[reactomeID]]
  reactomeName=ifelse(is.null(reactomeName),reactomeID,reactomeName)
  data.frame(
    reactomeID = reactomeID,
    numGenes = sum( isMember ),
    avgLFC = mean( res2$log2FoldChange[isMember] ),
    strength = sum( res2$log2FoldChange[isMember] ) / sqrt(sum(isMember)),
    pvalue = t.test( res2$log2FoldChange ~ isMember )$p.value,
    reactomeName = reactomeName ) }

# Now, recursively apply this program to all of the reactomeID rows in incm
reactomeResult = do.call( rbind, lapply( rownames(incm), testCategory ) )

# Finally, apply an FDR adjusted p-value to each catagory
reactomeResult$padjust = p.adjust( reactomeResult$pvalue, "BH" )
reactomeResult$zscore = qnorm(1-reactomeResult$padjust/2)*sign(reactomeResult$avgLFC)
reactomeResult_ordered = reactomeResult[order(reactomeResult$padjust),]
reactomeResult_ordered = reactomeResult_ordered[reactomeResult_ordered$padjust<0.05,]
#write.csv(reactomeResult_ordered,'~/Desktop/reactomeResult_ordered.csv')




