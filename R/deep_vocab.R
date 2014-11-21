getDeepVocab <- function(){
	require(jsonlite)
	fn <- system.file(file.path("extdata", "deepVocab.json"), package = "deepRtools")
	vocab <- fromJSON(fn)
	return(vocab)
}

