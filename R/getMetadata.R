#' saveReadTable
#' 
#' save version of read.table, which can adjust for tabs that have been replaced by 4 or more spaces
#' @param fn file name
#' @param ... arguments passed on to read.table
#' @return result of read.table(fn, ...), after possible error adjustment
#' @author Fabian Mueller
#' @noRd
saveReadTable <- function(fn, ...){
	dots <- list(...)
	sep <- dots[["sep"]]
	res <- tryCatch(
		read.table(fn, ...),
		error=function(err) {
			if (grepl("^line.*did not have.*elements", err$message)){
				num.elems <- as.numeric(gsub("^line.*did not have(.*)elements", "\\1", err$message))
				if (sep=="\t"){
					warning(paste0("Trying to repair tab/space separation (", fn, ") due to error in reading file: ", err$message))
					lls <- readLines(fn)
					# replace 4 or more whitespaces by tab
					lls <- gsub(" {4,}", "\t", lls)
					# collapse too many or too few tabs
					lls.split <- strsplit(lls, "\t")
					num.elems <- length(lls.split[[1]])
					lls.adj <- unlist(lapply(lls.split, FUN=function(x){
						if (length(x) == num.elems){
							paste(x, collapse="\t")
						} else if (length(x) > num.elems){
							# more tabs than expected --> use regular tab separation for the first elements
							# and collapse superfluent elements into one element
							paste(paste(x[1:(num.elems-1)], collapse="\t"), paste(x[num.elems:length(x)], collapse=""), sep="\t")
						} else {
							# fewer tabs than expected --> append tabs
							paste(c(x, rep("", num.elems-length(x))), collapse="\t")
						}
					}))

					tmpFn <- tempfile()
					writeLines(lls, tmpFn)
					read.table(tmpFn, ...)
				} else {
					stop(paste0("Error in reading table file (", fn, "): ", err$message))
				}
			} else {
				stop(paste0("Error in reading table file (", fn, "): ", err$message))
			}
		}
	)
	return(res)
}

#' getDataFrameFromTabFile
#' 
#' given a two column, tab-separated file, extract relevant information as a data frame
#' @param fn file name
#' @return a dataframe with one row containing relevant sample information
#' @author Fabian Mueller
#' @export 
getDataFrameFromTabFile <- function(fn){
	tt <- saveReadTable(fn, sep="\t", colClasses=rep("character",2), stringsAsFactors=FALSE, comment.char="")
	res <- data.frame(t(matrix(tt[,2])),stringsAsFactors=FALSE)
	colnames(res) <- tt[,1]
	return(res)
}

#' getSampleMdFromFile
#' 
#' given a sample metadata file, extract relevant information as a data frame
#' @param fn sample metadata file name
#' @param lenient do not enforce required fields, just return a warning if one of them is missing
#' @return a dataframe with one row containing relevant sample information
#' @author Fabian Mueller
#' @export 
getSampleMdFromFile <- function(fn, lenient=FALSE){
	requiredFields <- c("DEEP_SAMPLE_ID","ORGANISM","POOL","STATUS","DISEASE","BIOMATERIAL_PROVIDER","BIOMATERIAL_TYPE","CELL_TYPE","MARKERS","DONOR_ID","DONOR_AGE","DONOR_HEALTH_STATUS","DONOR_SEX","DONOR_ETHNICITY","COMMENT","SORTING_PARAMETERS")
	fdata <- getDataFrameFromTabFile(fn)
	presentFields <- colnames(fdata)
	missingFields <- setdiff(requiredFields, presentFields)
	additionalFields <- setdiff(presentFields, requiredFields)
	if (length(missingFields)>0){
		msg <- paste0("Missing fields in sample metadata: ",paste(missingFields,collapse=","))
		if (lenient){
			warning(msg)
		} else {
			stop(msg)
		}
	}
	cols <- c(requiredFields, additionalFields)
	cols.pres <- c(intersect(requiredFields, presentFields), additionalFields)
	res <- data.frame(rep(list(NA),length(cols)))
	colnames(res) <- cols
	res[,cols.pres] <- fdata[,cols.pres]
	return(res)
}

#' getExperimentMdFromFile
#' 
#' given a experiment metadata file, extract relevant information as a data frame
#' @param fn sample metadata file name
#' @param assay assay type. Possible values: 'wgbs','chip','rna','chracc'
#' @param lenient do not enforce required fields, just return a warning if one of them is missing
#' @return a dataframe with one row containing relevant sample information
#' @author Fabian Mueller
#' @export 
getExperimentMdFromFile <- function(fn, assay=NULL, lenient=FALSE){
	baseFields <- c("EXPERIMENT_ID","SAMPLE_ID","TECHNICAL_REPLICATE_ID","BIOLOGICAL_REPLICATE_ID","CONTACT_PERSON","CONTACT_LAB","EXPERIMENT_DATE","EXPERIMENT_TYPE","LIBRARY_GENERATION_PROTOCOL","LIBRARY_GENERATION_INITIAL_INPUT_QNTY","LIBRARY_GENERATION_PCR_F_PRIMER_SEQUENCE","LIBRARY_GENERATION_PCR_NUMBER_CYCLES","LIBRARY_GENERATION_PCR_POLYMERASE_TYPE","LIBRARY_GENERATION_PCR_PRIMER_CONC","LIBRARY_GENERATION_PCR_PRODUCT_ISOLATION_PROTOCOL","LIBRARY_GENERATION_PCR_R_PRIMER_SEQUENCE","LIBRARY_GENERATION_PCR_TEMPLATE_CONC","LIBRARY_GENERATION_PCR_THERMOCYCLING_PROGRAM","LIBRARY_GENERATION_PCR_PRODUCT_LENGTH","LIBRARY_GENERATION_FRAGMENT_SIZE_SELECTION","LIBRARY_GENERATION_cDNA_FRAGMENTATION")
	fields.wgbs <- c("EXTRACTION_PROTOCOL_SONICATION_FRAGMENT_SIZE_RANGE","BISULFITE_CONVERSION_PROTOCOL","LIBRARY_GENERATION_ADAPTOR_LIGATION_PROTOCOL","LIBRARY_GENERATION_ADAPTOR_SEQUENCE")
	fields.chip <- c("CHIP_PROTOCOL","CHIP_PROTOCOL_CHROMATIN_AMOUNT","CHIP_ANTIBODY","CHIP_ANTIBODY_CATALOG","CHIP_ANTIBODY_LOT","CHIP_ANTIBODY_PROVIDER","CHIP_PROTOCOL_ANTIBODY_AMOUNT","CHIP_PROTOCOL_BEAD_AMOUNT","CHIP_PROTOCOL_BEAD_TYPE","CHIP_PROTOCOL_PERC_CH2O_CROSSLINKING","CHIP_PROTOCOL_DURATION_CROSSLINKING","CHIP_PROTOCOL_SHEARING_METHOD","CHIP_PROTOCOL_SHEARING_CYCLES","CHIP_PROTOCOL_SHEARING_FREQUENCY","CHIP_PROTOCOL_IPSTAR_CONDITIONS")
	fields.chracc <- c("ChromatinAccessibility_METHOD","ChromatinAccessibility_PROTOCOL","ChromatinAccessibility_CHROMATIN_AMOUNT","ChromatinAccessibility_ENZYME","ChromatinAccessibility_ENZYME_PROVIDER","ChromatinAccessibility_ENZYME_LOT","ChromatinAccessibility_CHROMATIN_STATUS","ChromatinAccessibility_ENZYME_AMOUNT")
	fields.rna <- c("RNA_PROTOCOL_ENRICHMENT","RNA_cDNA_PREPARATION_INITIAL_RNA_QNTY","RNA_cDNA_PREPARATION_PCR_NUMBER_CYCLES","RNA_cDNA_PREPARATION_REVERSE_TRANSCRIPTION_PROTOCOL","RNA_EXTRACTION_PROTOCOL_FRAGMENTATION","RNA_PREPARATION_FRAGMENT_SIZE_RANGE","RNA_PREPARATION_INITIAL_QNTY","RNA_EXTRACTION_PROTOCOL","RNA_PREPARATION_3'_RNA_ADAPTER_LIGATION_PROTOCOL","RNA_PREPARATION_3'_RNA_ADAPTER_SEQUENCE","RNA_PREPARATION_5'_DEPHOSPHORYLATION","RNA_PREPARATION_5'_PHOSPHORYLATION","RNA_PREPARATION_5'_RNA_ADAPTER_LIGATION_PROTOCOL","RNA_PREPARATION_5'_RNA_ADAPTER_SEQUENCE","RNA_PREPARATION_REVERSE_TRANSCRIPTION_PRIMER_SEQUENCE","RNA_PREPARATION_REVERSE_TRANSCRIPTION_PROTOCOL")
	requiredFields <- baseFields
	if (is.character(assay)) {
		if (assay == "wgbs"){
			requiredFields <- c(requiredFields, fields.wgbs)
		} else if (assay == "chip"){
			requiredFields <- c(requiredFields, fields.chip)
		} else if (assay == "rna"){
			requiredFields <- c(requiredFields, fields.rna)
		} else if (assay == "chracc"){
			requiredFields <- c(requiredFields, fields.chracc)
		} else {
			stop("Unknown assay type")
		}
	}
	fdata <- getDataFrameFromTabFile(fn)
	presentFields <- colnames(fdata)
	missingFields <- setdiff(requiredFields, presentFields)
	additionalFields <- setdiff(presentFields, requiredFields)
	if (length(missingFields)>0){
		msg <- paste0("Missing fields in experiment metadata: ",paste(missingFields,collapse=","))
		if (lenient){
			warning(msg)
		} else {
			stop(msg)
		}
	}
	cols <- c(requiredFields, additionalFields)
	cols.pres <- c(intersect(requiredFields, presentFields), additionalFields)
	res <- data.frame(rep(list(NA),length(cols)))
	colnames(res) <- cols
	res[,cols.pres] <- fdata[,cols.pres]
	return(res)
}
