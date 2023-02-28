rMSIprocPeakMatrix <- rMSI2::LoadPeakMatrix("D:/DATOMA/MSIData/MultipleBrainsPeakMatrix.pkmat")
params <- rMSI2::ProcessingParameters()
params$peakAnnotation$ppmMassTolerance <- 10
params$peakAnnotation$isotopeLikelihoodScoreThreshold <- 0.7
rMSI2::peakAnnotation(rMSIprocPeakMatrix = rMSIprocPeakMatrix, params = params)

for (i in 1:50)
{
  params$peakAnnotation$ppmMassTolerance <- params$peakAnnotation$ppmMassTolerance + 1
  rMSI2::peakAnnotation(rMSIprocPeakMatrix = rMSIprocPeakMatrix, params = params)
}






