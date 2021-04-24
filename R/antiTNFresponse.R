#' ATRPred: Anti-TNF Treatment Response Predictor
#'
#' This function predicts the response to anti-TNF therapy in Rheumatoid Arthritis patients at baseline
#' It requires an excel input for 17 proteins along with gender and baseline DAS
#' Any values not present will be imputed
#' Please don't change the order of proteins in the sample excel file
#'
#' @export
antiTNFresponse <- function(){
  options(warn=-1, verbose = FALSE) # Comment to not suppress the warnings
	print("##### ATRPred: Anti-TNF Treatment Response Predictor #####")
	print("Please select the excel file having pateint's NPX values:")
	patient <- readxl::read_excel(file.choose())
	gen=patient[2,2]
	newRow <- patient[c(1,3:19),2]
	rownames(newRow)=colnames(spl)
	consolidated <- rbind(spl,t(newRow))

#	require('RANN')
	hpc <- caret::preProcess(consolidated, method = c("knnImpute","center", "scale"))
	transformed <- stats::predict(hpc, newdata = consolidated)

	val <- cbind(gen,transformed[90,])
	wt <- c(0.116,2.133,-2.126,-2.068,0.421,2.488,-2.595,-0.960,-0.651,2.557,-0.830,-0.758,1.281,0.744,-0.421,2.661,-0.243,2.990,-2.574)
	S <- sum(wt*val)+3.8
	t <- 0.7136
	p <- 1/(1+exp(-S+t))
	print(paste0(sprintf('Calculated probabilty of response = %.2f',100*p),"%"))
	if(S>t) print('The patient is a RESPONDER.') else print('The patient is a NON-RESPONDER.')
}
