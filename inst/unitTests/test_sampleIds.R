test_sampleId_tab <- function(){
	fn <- system.file(file.path("extdata", "sampleIds_20150123.txt"), package = "deepRtools")
	ids <- readLines(fn)
	df <- deepSampleIds2dataFrame(ids)
	# df[sort(sample(1:nrow(df),10)),] #show some examples
	# checkEquals(sum(!df$validId),0)
	checkTrue(all(df$validId))
}

test_sampleIds <- function(){
	require(RUnit)
	require(deepRtools)

	test_sampleId_tab()
}

test_sampleIds()
