rMSIprocPeakMatrix <- rMSI2::LoadPeakMatrix("D:/DATOMA/MSIData/MultipleBrainsPeakMatrix.pkmat")
params <- rMSI2::ProcessingParameters()
params$peakAnnotation$isotopeLikelihoodScoreThreshold <- 0.7
params$peakAnnotation$ppmMassTolerance <- 30
ann  <- rMSI2::peakAnnotation(rMSIprocPeakMatrix = rMSIprocPeakMatrix, params = params)
rMSI2::plotAnnotatedSpectra(rMSIprocPeakMatrix, ann)
