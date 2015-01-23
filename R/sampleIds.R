################################################################################
# Globals
################################################################################
DEEP_ORGANS <- unlist(getDeepVocab()[["DEEP"]][["ORGAN"]])
DEEP_CELL_TYPES <- unlist(getDeepVocab()[["DEEP"]][["CELL_TYPE"]])
DEEP_DISEASES <- unlist(getDeepVocab()[["DEEP"]][["DISEASE"]])
DEEP_ASSAYS <- unlist(getDeepVocab()[["DEEP"]][["ASSAY"]])
DEEP_CENTERS <- unlist(getDeepVocab()[["DEEP"]][["CENTER"]])


REGEX_DEEP_SAMPLE_PREFIX <- paste0(
	"(?P<sp>([45][1-4]|0[01]))_",
	"(?P<organism>[HM])(?P<sex>[mf])(?P<donor>[0-9]{2})_",
	"(?P<organ>",paste(names(DEEP_ORGANS),collapse="|"),")","(?P<celltype>",paste(names(DEEP_CELL_TYPES),collapse="|"),")_",
	"(?P<disease>",paste(names(DEEP_DISEASES),collapse="|"),")(?P<breplicate>[0-9]?)"
)
REGEX_DEEP_SAMPLE_SUFFIX <- paste0(
	"(?P<assay>",paste(names(DEEP_ASSAYS),collapse="|"),")_",
	"(?P<center>",paste(names(DEEP_CENTERS),collapse="|"),")",
	"(_(?P<treplicate>[1-9]))?",
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
#' @return A dataframe with columns corresponding to the fields that can be inferred from the names
#' @author Fabian Mueller
#' @export 
#' @examples 
#' deepSampleIds2dataFrame(c("43_Hm03_BlMa_TO1_WGBS_E","43_Hm05_BlMa_Ct","43_Hm05_BlMa_Ct_NOMe_S","01_HepaRG_LiHR_D32","41_Hf01_LiHe_Ct1_H3K4me1_F"))
deepSampleIds2dataFrame <- function(ids){
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
			subproject=NA,
			organism=NA,
			sex=NA,
			donor=NA,
			organ=NA,
			cellType=NA,
			disease=NA,
			replicate=NA,
			assay=NA,
			center=NA
		)
		if (!is.na(re)){
			rem <- regexpr(re,x,perl=TRUE)
			parsed <- parse.one(x,rem)
			df[["validId"]] <- TRUE
			df[["subproject"]] <- paste0("SP",paste(strsplit(parsed[1,"sp"],"")[[1]],collapse="."))
			df[["organism"]]   <- ifelse(parsed[1,"organism"]=="H","human",ifelse(parsed[1,"organism"]=="M","mouse",NA))
			df[["sex"]]   <- ifelse(parsed[1,"sex"]=="m","male",ifelse(parsed[1,"sex"]=="f","female",NA))
			df[["donor"]] <- parsed[1,"donor"]
			df[["organ"]] <- DEEP_ORGANS[parsed[1,"organ"]]
			df[["cellType"]] <- DEEP_CELL_TYPES[parsed[1,"celltype"]]
			df[["disease"]] <- DEEP_DISEASES[parsed[1,"disease"]]
			if (nchar(parsed[1,"breplicate"])>0) df[["breplicate"]] <- parsed[1,"breplicate"]
		}
		if (fullId){
			df[["assay"]] <- DEEP_ASSAYS[parsed[1,"assay"]]
			df[["center"]] <- DEEP_CENTERS[parsed[1,"center"]]
			if (nchar(parsed[1,"treplicate"])>0) df[["treplicate"]] <- parsed[1,"treplicate"]
		}
		return(df)
	}))
	return(res)
}

