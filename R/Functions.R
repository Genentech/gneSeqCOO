

#' Combined COO function
#' @description wrapper for running coo_normalize and coo_predict together
#' @param cds - a DESeqDataSet object, with raw counts in the "counts" slot.
#' @details This function is a wrapper for the two core COO prediction functions.
#'     This function intentionally forces a standard procedure for both 
#'     normalization and modeling, as this should give us the best chance of 
#'     the model matching up with our training data.
#' @return a data frame of COO predictions. See coo_predict for more details
#' @export

coo_rnaseq <- function(cds){
	cat("Normalizing count matrix...")
	cds_norm = coo_normalize(cds)
	cat("Done!\n")
	result = coo_predict(cds_norm)
	return(result)
}

#' Combined DZsig function
#' @description wrapper for running coo_normalize and dzsig_predict together
#' @param cds - a DESeqDataSet object, with raw counts in the "counts" slot.
#' @details This function is a wrapper for the core normalization and prediction functions.
#' @return a data frame of DZsig predictions. See dzsig_predict for more details
#' @export

dzsig_rnaseq <- function(cds){
	cat("Normalizing count matrix...")
	cds_norm = coo_normalize(cds)
	cat("Done!\n")
	result = dzsig_predict(cds_norm)
	return(result)
}

#' Refined COO prediction function
#' @description Wrapper for running coo_normalize, coo_predict, and dzsig_predict
#' @param cds - a DESeqDataSet object, with raw counts in the "counts" slot.
#' @details This function is a wrapper that runs the normalization and both COO
#'    and DZsig predictions, and then combines the results into a single data frame, 
#'    including the refined COO class.
#' @return a data frame with six columns: 
#'    * Sample names
#'    * COO LPS score prediction
#'    * COO class
#'    * DZsig LPS score prediction
#'    * DZsig class
#'    * Refined COO class
#' @export

refined_coo_rnaseq <- function(cds){
	
	cat("Normalizing count matrix...")
	cds_norm = coo_normalize(cds)
	cat("Done!\n")
	
	cat("Predicting COO...")
	coo = coo_predict(cds_norm)
 	colnames(coo) <- c("Sample", "LPS_COO", "COO")
	cat("Predicting DZsig...")
 	dzsig = dzsig_predict(cds_norm)
 	colnames(dzsig) <- c("Sample", "LPS_DZsig", "DZsig")
 	cat("Done!\n")

	result = merge(coo,dzsig,by="Sample")
	result$refined_COO <-
		ifelse(result$COO == "ABC", "ABC",
			ifelse(result$DZsig == "DZsig-POS", "DZsig-POS", as.character(result$COO))
		)
	return(result)
}





#' COO Normalization function
#' @description The standard normalization function for RNASeq data
#' @param cds - a DESeqDataSet object, with raw counts in the "counts" slot. 
#' 
#' @details this function applies a robust library size normalization procedure
#'    to the raw count matrix provided. The goal is to specifically ensure that 
#'    the expression values for the data are on a consistent scale, regardless 
#'    of the source dataset. For more details on the robust library size 
#'    normalization function, refer to the "COO_Classifier" supplemental 
#'    vignette.
#'
#'    Please note that rownames of the raw count matrix must be either refseq 
#'    IDs or ENSEMBL gene IDs. Gene symbols will not work. 
#' 
#'    Notably, for the purposes of machine learning, we do not care about the 
#'    absolute expression of genes, but only their expression relative to other
#'    samples. This, combined with the fact that people may not know how to get
#'    read lengths easily for new data means that we're skipping any
#'    normalization step involving read length.
#' @importFrom stats median
#' @importFrom methods new
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom DESeq2 counts
#' @importFrom Biobase ExpressionSet
#' @return A normalized DESeqDataSet
#' @export

coo_normalize <- function(cds){

	if(nrow(cds)<2000){
		warning(sprintf("Only %i genes present in the CDS.\ncoo_normalize uses the entire expression profile for normalization, so it is not recommended to subset the count matrix before running."))
	}

	## figure out the IGIS version, if not provided
	ident = get_id_type(cds)

	## Step 1: remove low count genes (because they affect libSize calculation)
	## I'm gonna use a simple heuristic of "minimum 1 count per sample"
	cts <- rowSums(counts(cds))
	ind <- (cts < ncol(cds))
	# ind <- (cts == 0)
	cds2 <- cds[!ind,]

	## add a flag to the fData matrix
	fData = data.frame(rowData(cds))
	## If the user hasn't initialized the fData matrix, initialize with just IDs
	if(is.null(fData) || ncol(fData)==0){  
		fData = data.frame(ID=rownames(cds))
		rownames(fData) = rownames(cds)
	}
	fData$isLowCountGene = ind 
	if(is.null(rownames(fData))){
		rownames(fData) = fData$ID  ## Not sure how these get lost, but they do
	}

  ## This is the updated "library size" calculation, optimized for use with
  ## a single sample. This calculates a DESeq2-esque size factor compared with
  ## a 'reference' sample from GOYA, and 
	libSizes = refLibrarySize(cds2,ident)

  ## now actually normalize the data
  normMat = t( t(log(counts(cds)+1)) + log(10^6/libSizes) )
	
	## I'm just going to send a warning about samples with very low library size
	n = sum(libSizes<5*10^5)
	if(n>0){
		warning(sprintf("%i samples with low library size detected. 
			COO estimates may be unreliable on these samples.",n))
	}

	## add some of our data to the pData matrices
	pData = data.frame(colData(cds))
	pData$LibrarySize  = libSizes
	pData$isLowCountSample = libSizes<5*10^5

	## Wrap up the normalized matrix in a new expressionSet and return
	ExpressionSet(assayData = normMat,
		phenoData   = new("AnnotatedDataFrame",data=pData),
		featureData = new("AnnotatedDataFrame",data=fData)
	)
}


#' refLibrarySize
#' @description Function for library ref size
#' @param cds2 - a DESeqDataSet object, with raw counts in the "counts" slot. 
#' @param ident - Identifiers used in the dataset.
#'        Should be one of c("refseq","ensembl")
#'        If not provided, the IGIS version will be inferred from the data.
#'
#' @return vector describing the robust library size for the sample
#' @export
#'

refLibrarySize <- function( cds2, ident ){
	goyaRefSample = list("refseq"=goyaRefSample3, "ensembl"=goyaRefSample4)[[ident]]
  ## match up the rownames of cds2 with our reference dataset
  keep = rownames(cds2) %in% names(goyaRefSample)
  cds2 = cds2[keep,]
  i = match(rownames(cds2), names(goyaRefSample))

  ## For each sample in "cds2", calculate a library size
  apply(counts(cds2), 2, function(x){
    sizeFactor = median((x / goyaRefSample[i])[x!=0])
    goyaLibSize * sizeFactor
  })
}







#' COO prediction function
#' @description Function for predicting COO using a normalized RNASeq dataset
#' @param norm_eset A normalized ExpressionSet, as returned by coo_normalize
#' @param model     The fitted COO model.
#' @details This function uses the model generated on the GOYA dataset to
#'      predict Cell of Origin on a normalized expression dataset. No subsetting
#'      needs to be performed, but the expression matrix must have either refseq
#'      (e.g. "geneID:7900") or ENSEMBL (e.g. ENSG....) gene IDs as the rownames
#' @return a data frame with three columns:
#'    * Sample names
#'    * LPS score prediction
#'    * COO class
#' @importFrom stats predict
#' @importFrom Biobase fData exprs
#' @export

coo_predict <- function(norm_eset, model) {
    ## figure out the IGIS version, if not provided
    ids = get_id_type(norm_eset)

    ## Check the model, and if not provided, load the default model
    if (missing(model)) {
        model = list("refseq" = get("coo_GOYA_21gene_igis3"), "ensembl" = get("coo_GOYA_21gene_igis4"))[[ids]]
        cat("Model not specified. Using the 21-gene GOYA model...\n")
    }

    coefs = model$coefficients[-1]
    coefNames = gsub("\\`", "", names(coefs))

    ## Check that all the genes are present in the dataset
    if (!all(coefNames %in% rownames(norm_eset))) {
        missing = coefNames[!(coefNames %in% rownames(norm_eset))]
        stop(
            "Some required genes are missing.\n",
            "Missing genes: ", paste(missing, collapse = " ")
        )
    }

    ## Check if any of the genes are flagged as "low expression"
    if (any(fData(norm_eset[coefNames, ])$isLowCountGene)) {
        warning(
            "Some genes in the classifer were flagged as low count. ",
            "Predictions may be lower confidence."
        )
    }

    ## pull out the data
    x = as.data.frame(t(exprs(norm_eset)[coefNames, ]))

    lps = predict(model, newdata = x)
    coo = cut(lps, c(-Inf, 1920, 2430, Inf),
        labels = c("GCB", "UNCLASSIFIED", "ABC")
    )

    out = data.frame(Sample = names(lps), LPS = lps, COO = coo)
    rownames(out) = NULL
    return(out)
}


#' DZsig prediction function
#' @description Function for predicting DZsig using a normalized RNASeq dataset
#' @param norm_eset A normalized ExpressionSet, as returned by coo_normalize
#' @param model     The fitted DZsig model.
#' @details This function uses the model generated on the GAMBL dataset to
#'      predict Dark Zone Signature (DZsig) on a normalized expression dataset. No subsetting
#'      needs to be performed, but the expression matrix must have either refseq
#'      (e.g. "geneID:7900") or ENSEMBL (e.g. ENSG....) gene IDs as the rownames. 
#'      Thresholds used for classification are -15.6 and -6.3, as per the original DLBCL90 publication 
#'      (Ennishi et al, *J Clin Oncol* 2019).
#' @return a data frame with three columns: 
#'    * Sample names
#'    * LPS score prediction
#'    * DZsig class
#' @importFrom stats predict
#' @importFrom Biobase fData exprs
#' @export

dzsig_predict <- function(norm_eset, model){
  
	## figure out the IGIS version, if not provided
	ids = get_id_type(norm_eset)

	## Check the model, and if not provided, load the default model
	if(missing(model)){
		model = list("refseq"=get("dzsig_GAMBL_29gene_igis3"), "ensembl"=get("dzsig_GAMBL_29gene_igis4"))[[ids]]
		cat("Model not specified. Using the 29-gene GAMBL model...\n")
	}

	coefs = model$coefficients[-1]
	coefNames = gsub("\\`","",names(coefs))

	## Check that all the genes are present in the dataset
	if(!all(coefNames %in% rownames(norm_eset))){
		missing = coefNames[!(coefNames %in% rownames(norm_eset))]
		stop("Some required genes are missing.\n",
			"Missing genes: ",paste(missing,collapse=" "))
	}

	## Check if any of the genes are flagged as "low expression"
	if(any(fData(norm_eset[coefNames,])$isLowCountGene)){
		warning("Some genes in the classifer were flagged as low count. ",
			"Predictions may be lower confidence.")
	}

	## pull out the data
	x = as.data.frame(t(exprs(norm_eset)[coefNames,]))

	lps = predict(model, newdata=x)
	dzsig = cut(lps, c(-Inf, -15.6, -6.3, Inf), 
		labels=c("DZsig-NEG","DZsig-IND","DZsig-POS"))

	out = data.frame(Sample=names(lps), LPS=lps, DZsig=dzsig)
	rownames(out) = NULL
	return(out)
}


#' Attempt to guess the identifiers used in the model
#' @description Function for attempting to guess identifier type
#' @param cds - a DESeqDataSet object, with raw counts in the "counts" slot.
#' @return string of IGIS version 3.0 or 4.0
#' @export

get_id_type <- function(cds){
  i3 = sum(rownames(cds) %in% names(goyaRefSample3))
  i4 = sum(rownames(cds) %in% names(goyaRefSample4))
  
  if(i3 > i4 & i3 > 2000){
    return("refseq")
  }
  if(i4 > i3 & i4 > 2000){
    return("ensembl")
  }
  stop("Could not identify identifier type")
}




