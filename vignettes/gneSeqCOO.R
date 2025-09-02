## ----eval=FALSE---------------------------------------------------------------
# 
# library(DESeq2)
# cds = DESeqDataSet(strandedCtSummary,design=~1)
# 

## ----eval=FALSE---------------------------------------------------------------
# ## counts is a matrix, with rows of features and columns of samples.
# cds = DESeqDataSetFromMatrix(counts,
#   colData=data.frame(ID=colnames(counts)),
#   rowData=data.frame(ID=rownames(counts)),
#   design=~1)

## ----echo=FALSE,message=FALSE-------------------------------------------------
library(gneSeqCOO)
themodel=coo_GOYA_21gene_igis3$coefficients
names(themodel)[names(themodel)!="(Intercept)"]

## ----eval=FALSE---------------------------------------------------------------
# 
# pred = coo_rnaseq(cds)
# 

## ----eval=FALSE---------------------------------------------------------------
# norm = coo_normalize(cds)

## ----eval=FALSE---------------------------------------------------------------
# pred = coo_predict(norm_main)

## ----eval=FALSE---------------------------------------------------------------
# dzsig_pred = dzsig_rnaseq(cds)

## ----eval=FALSE---------------------------------------------------------------
# refined_pred = refined_coo_rnaseq(cds)

