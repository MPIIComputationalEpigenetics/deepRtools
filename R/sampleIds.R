################################################################################
# Globals
################################################################################
DEEP_ORGANS <- unlist(getDeepVocab()[["DEEP"]][["ORGAN"]])
DEEP_CELL_TYPES <- unlist(getDeepVocab()[["DEEP"]][["CELL_TYPE"]])
DEEP_DISEASES <- unlist(getDeepVocab()[["DEEP"]][["DISEASE"]])
DEEP_ASSAYS <- unlist(getDeepVocab()[["DEEP"]][["ASSAY"]])
DEEP_CENTERS <- unlist(getDeepVocab()[["DEEP"]][["CENTER"]])
DEEP_CELL_LINES <- getDeepVocab()[["DEEP"]][["CELL_LINE"]]


REGEX_DEEP_SAMPLE_PREFIX <- paste0(
	"(?P<sp>([45][1-4]|0[01]))_",
	"((?P<organism>[HM])(?P<sex>[mf])(?P<donor>[0-9a-z]{2})|(?P<cellline>",paste(names(DEEP_CELL_LINES),collapse="|"),"))_",
	"(?P<organ>",paste(names(DEEP_ORGANS),collapse="|"),")","(?P<celltype>",paste(names(DEEP_CELL_TYPES),collapse="|"),")_",
	"(?P<disease>",paste(names(DEEP_DISEASES),collapse="|"),")(?P<breplicate>[0-9]?)"
)
REGEX_DEEP_SAMPLE_SUFFIX <- paste0(
	"(?P<assay>",paste(names(DEEP_ASSAYS),collapse="|"),")(_",
	"(?P<center>",paste(names(DEEP_CENTERS),collapse="|"),")",
	"(_(?P<treplicate>[1-9]))?)?"
)
REGEX_DEEP_SAMPLE_ID <- paste0("^",REGEX_DEEP_SAMPLE_PREFIX,"$")
REGEX_DEEP_SAMPLE_ID_FULL <- paste0("^",REGEX_DEEP_SAMPLE_PREFIX,"_",REGEX_DEEP_SAMPLE_SUFFIX,"$")

################################################################################
# Helpers
################################################################################
#regex parsing (from https://stat.ethz.ch/pipermail/r-devel/2011-February/060086.html)
parse.one <- function(string,result){
	m <- do.call(rbind,lapply(seq_along(string),function(i){
	st <- attr(result,"capture.start")[i,]
	substring(string[i],st,st+attr(result,"capture.length")[i,]-1)
	}))
	colnames(m) <- attr(result,"capture.names")
	m
}
################################################################################
# Methods
################################################################################
#' deepSampleIds2dataFrame
#' 
#' parses a vector of sample IDs and returns avector with the corresponding fields
#' @param ids a vector of sample ids as strings
#' @param ignore.case flag indicating whether to ignore casing in parsing the sample ids
#' @return A dataframe with columns corresponding to the fields that can be inferred from the names
#' @author Fabian Mueller
#' @export 
#' @examples 
#' deepSampleIds2dataFrame(c("43_Hm03_BlMa_TO1_WGBS_E_1", "43_Hm03_BlMa_TO1_WGBS","43_Hm05_BlMa_Ct","43_Hm05_BlMa_Ct_NOMe_S_2","01_HepaRG_LiHR_D32","41_Hf01_LiHe_Ct1_H3K4me1_F_1"))
deepSampleIds2dataFrame <- function(ids, ignore.case=FALSE){
	res <- do.call("rbind",lapply(ids,FUN=function(x){
		re <- NA
		fullId <- FALSE
		if (grepl(REGEX_DEEP_SAMPLE_ID,x,perl=TRUE)) re <- REGEX_DEEP_SAMPLE_ID
		if (grepl(REGEX_DEEP_SAMPLE_ID_FULL,x,perl=TRUE)){
			re <- REGEX_DEEP_SAMPLE_ID_FULL
			fullId <- TRUE
		}
		df <- data.frame(
			id=x,
			validId=FALSE,
			isCellLine=FALSE,
			subproject=NA,
			organism=NA,
			sex=NA,
			donor=NA,
			organ=NA,
			cellType=NA,
			disease=NA,
			breplicate=NA,
			assay=NA,
			center=NA,
			treplicate=NA,
			subproject_short=NA,
			organ_short=NA,
			cellType_short=NA,
			disease_short=NA,
			assay_short=NA,
			center_short=NA
		)
		if (!is.na(re)){
			rem <- regexpr(re,x,perl=TRUE, ignore.case=ignore.case)
			parsed <- parse.one(x,rem)
			
			df[["validId"]] <- TRUE
			df[["subproject"]] <- paste0("SP",paste(strsplit(parsed[1,"sp"],"")[[1]],collapse="."))
			df[["subproject_short"]] <- parsed[1,"sp"]
			df[["organism"]]   <- ifelse(parsed[1,"organism"]=="H","human",ifelse(parsed[1,"organism"]=="M","mouse",NA))
			df[["sex"]]   <- ifelse(parsed[1,"sex"]=="m","male",ifelse(parsed[1,"sex"]=="f","female",NA))
			if (nchar(parsed[1,"donor"])>0) df[["donor"]] <- parsed[1,"donor"]
			df[["organ"]] <- DEEP_ORGANS[parsed[1,"organ"]]
			df[["organ_short"]] <- parsed[1,"organ"]
			df[["cellType"]] <- DEEP_CELL_TYPES[parsed[1,"celltype"]]
			df[["cellType_short"]] <- parsed[1,"celltype"]
			df[["disease"]] <- DEEP_DISEASES[parsed[1,"disease"]]
			df[["disease_short"]] <- parsed[1,"disease"]
			if (nchar(parsed[1,"breplicate"])>0) df[["breplicate"]] <- parsed[1,"breplicate"]

			cln <- parsed[1,"cellline"]
			if (nchar(cln)>0) {
				df[["isCellLine"]] <- TRUE
				df[["organism"]] <- DEEP_CELL_LINES[[cln]]$organism
				df[["sex"]] <- DEEP_CELL_LINES[[cln]]$sex
			}
		}
		if (fullId){
			df[["assay"]] <- DEEP_ASSAYS[parsed[1,"assay"]]
			df[["assay_short"]] <- parsed[1,"assay"]
			df[["center"]] <- DEEP_CENTERS[parsed[1,"center"]]
			df[["center_short"]] <- parsed[1,"center"]
			trep <- parsed[1,"treplicate"]
			if (nchar(trep)>0) df[["treplicate"]] <- trep
		}
		return(df)
	}))
	return(res)
}
#' deepSampleBaseIds
#' 
#' parses a vector of sample IDs and returns a vector with the base ids (without assay, center and technical replicate fields). Returns the same string if it is not a valid ID.
#' @param ids a vector of sample ids as strings
#' @param ignore.case flag indicating whether to ignore casing in parsing the sample ids
#' @return a vector of sample ids as strings
#' @author Fabian Mueller
#' @export 
#' @examples 
#' deepSampleBaseIds(c("43_Hm03_BlMa_TO1_WGBS_E_1","43_Hm05_BlMa_Ct","43_Hm05_BlMa_Ct_NOMe_S_2","01_HepaRG_LiHR_D32","41_Hf01_LiHe_Ct1_H3K4me1_F_1"))
deepSampleBaseIds <- function(ids, ignore.case=FALSE){
	res <- gsub(paste0("(",REGEX_DEEP_SAMPLE_PREFIX,")","(_",REGEX_DEEP_SAMPLE_SUFFIX,")?"),"\\1",ids, perl=TRUE, ignore.case=ignore.case)
	return(res)
}
#' getDeepSampleIdsFromString
#' 
#' looks for occurrences of DEEP IDs in a vector of strings and returns them as a list
#' @param x string vector in which to look for DEEP Ids
#' @param ignore.case flag indicating whether to ignore casing in parsing the sample ids
#' @param fullId only retrieve full ids
#' @return a list containing the found DEEP IDs for each element in the string vector
#' @author Fabian Mueller
#' @export 
#' @examples 
#' getDeepSampleIdsFromString(c("blubbsad43_Hm03_BlMa_TO1_WGBS_E_1XXX43_Hm05_BlMa_CtNOI43_Hm05_BlMa_Ct","43_Hm05_BlMa_Ct_NOMe_S_2uiadyf01_HepaRG_LiHR_D32","blubb"))
getDeepSampleIdsFromString <- function(x, fullId=FALSE, ignore.case=FALSE){
	pp <- REGEX_DEEP_SAMPLE_PREFIX
	if (fullId) pp <- paste0(REGEX_DEEP_SAMPLE_PREFIX,"_",REGEX_DEEP_SAMPLE_SUFFIX)
	m <- gregexpr(pp, x, perl=TRUE, ignore.case=ignore.case)
	res <- regmatches(x, m)
	return(res)
}
