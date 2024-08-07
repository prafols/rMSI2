// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// CNormalizationsAndMeans
List CNormalizationsAndMeans(Rcpp::List rMSIObj_list, int numOfThreads, double memoryPerThreadMB, Rcpp::NumericVector commonMassAxis);
RcppExport SEXP _rMSI2_CNormalizationsAndMeans(SEXP rMSIObj_listSEXP, SEXP numOfThreadsSEXP, SEXP memoryPerThreadMBSEXP, SEXP commonMassAxisSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type rMSIObj_list(rMSIObj_listSEXP);
    Rcpp::traits::input_parameter< int >::type numOfThreads(numOfThreadsSEXP);
    Rcpp::traits::input_parameter< double >::type memoryPerThreadMB(memoryPerThreadMBSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type commonMassAxis(commonMassAxisSEXP);
    rcpp_result_gen = Rcpp::wrap(CNormalizationsAndMeans(rMSIObj_list, numOfThreads, memoryPerThreadMB, commonMassAxis));
    return rcpp_result_gen;
END_RCPP
}
// CparseBrukerXML
List CparseBrukerXML(String xml_path);
RcppExport SEXP _rMSI2_CparseBrukerXML(SEXP xml_pathSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type xml_path(xml_pathSEXP);
    rcpp_result_gen = Rcpp::wrap(CparseBrukerXML(xml_path));
    return rcpp_result_gen;
END_RCPP
}
// testingimzMLBinWriteSequential
Rcpp::DataFrame testingimzMLBinWriteSequential(const char* ibdFname, Rcpp::String mz_dataTypeString, Rcpp::String int_dataTypeString, Rcpp::String str_uuid, Rcpp::NumericMatrix mzArray, Rcpp::NumericMatrix intArray);
RcppExport SEXP _rMSI2_testingimzMLBinWriteSequential(SEXP ibdFnameSEXP, SEXP mz_dataTypeStringSEXP, SEXP int_dataTypeStringSEXP, SEXP str_uuidSEXP, SEXP mzArraySEXP, SEXP intArraySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const char* >::type ibdFname(ibdFnameSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type mz_dataTypeString(mz_dataTypeStringSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type int_dataTypeString(int_dataTypeStringSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type str_uuid(str_uuidSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type mzArray(mzArraySEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type intArray(intArraySEXP);
    rcpp_result_gen = Rcpp::wrap(testingimzMLBinWriteSequential(ibdFname, mz_dataTypeString, int_dataTypeString, str_uuid, mzArray, intArray));
    return rcpp_result_gen;
END_RCPP
}
// CimzMLBinCreateNewIBD
void CimzMLBinCreateNewIBD(const char* ibdFname, Rcpp::String str_uuid);
RcppExport SEXP _rMSI2_CimzMLBinCreateNewIBD(SEXP ibdFnameSEXP, SEXP str_uuidSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const char* >::type ibdFname(ibdFnameSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type str_uuid(str_uuidSEXP);
    CimzMLBinCreateNewIBD(ibdFname, str_uuid);
    return R_NilValue;
END_RCPP
}
// CimzMLBinAppendMass
uint64_t CimzMLBinAppendMass(const char* ibdFname, Rcpp::String mz_dataTypeString, Rcpp::NumericVector mzNew);
RcppExport SEXP _rMSI2_CimzMLBinAppendMass(SEXP ibdFnameSEXP, SEXP mz_dataTypeStringSEXP, SEXP mzNewSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const char* >::type ibdFname(ibdFnameSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type mz_dataTypeString(mz_dataTypeStringSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mzNew(mzNewSEXP);
    rcpp_result_gen = Rcpp::wrap(CimzMLBinAppendMass(ibdFname, mz_dataTypeString, mzNew));
    return rcpp_result_gen;
END_RCPP
}
// CimzMLBinAppendIntensity
uint64_t CimzMLBinAppendIntensity(const char* ibdFname, Rcpp::String int_dataTypeString, Rcpp::NumericVector intNew);
RcppExport SEXP _rMSI2_CimzMLBinAppendIntensity(SEXP ibdFnameSEXP, SEXP int_dataTypeStringSEXP, SEXP intNewSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const char* >::type ibdFname(ibdFnameSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type int_dataTypeString(int_dataTypeStringSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type intNew(intNewSEXP);
    rcpp_result_gen = Rcpp::wrap(CimzMLBinAppendIntensity(ibdFname, int_dataTypeString, intNew));
    return rcpp_result_gen;
END_RCPP
}
// CimzMLBinWriteModifyMass
void CimzMLBinWriteModifyMass(const char* ibdFname, unsigned int NPixels, Rcpp::String mz_dataTypeString, Rcpp::String int_dataTypeString, bool continuous, Rcpp::NumericVector mzNew, uint64_t mzOffset);
RcppExport SEXP _rMSI2_CimzMLBinWriteModifyMass(SEXP ibdFnameSEXP, SEXP NPixelsSEXP, SEXP mz_dataTypeStringSEXP, SEXP int_dataTypeStringSEXP, SEXP continuousSEXP, SEXP mzNewSEXP, SEXP mzOffsetSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const char* >::type ibdFname(ibdFnameSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type NPixels(NPixelsSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type mz_dataTypeString(mz_dataTypeStringSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type int_dataTypeString(int_dataTypeStringSEXP);
    Rcpp::traits::input_parameter< bool >::type continuous(continuousSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mzNew(mzNewSEXP);
    Rcpp::traits::input_parameter< uint64_t >::type mzOffset(mzOffsetSEXP);
    CimzMLBinWriteModifyMass(ibdFname, NPixels, mz_dataTypeString, int_dataTypeString, continuous, mzNew, mzOffset);
    return R_NilValue;
END_RCPP
}
// CimzMLBinReadMass
Rcpp::NumericVector CimzMLBinReadMass(const char* ibdFname, unsigned int NPixels, unsigned int N, uint64_t offset, Rcpp::String dataTypeString, bool continuous);
RcppExport SEXP _rMSI2_CimzMLBinReadMass(SEXP ibdFnameSEXP, SEXP NPixelsSEXP, SEXP NSEXP, SEXP offsetSEXP, SEXP dataTypeStringSEXP, SEXP continuousSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const char* >::type ibdFname(ibdFnameSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type NPixels(NPixelsSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type N(NSEXP);
    Rcpp::traits::input_parameter< uint64_t >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type dataTypeString(dataTypeStringSEXP);
    Rcpp::traits::input_parameter< bool >::type continuous(continuousSEXP);
    rcpp_result_gen = Rcpp::wrap(CimzMLBinReadMass(ibdFname, NPixels, N, offset, dataTypeString, continuous));
    return rcpp_result_gen;
END_RCPP
}
// CimzMLBinReadIntensity
Rcpp::NumericVector CimzMLBinReadIntensity(const char* ibdFname, unsigned int NPixels, unsigned int N, uint64_t offset, Rcpp::String dataTypeString, bool continuous);
RcppExport SEXP _rMSI2_CimzMLBinReadIntensity(SEXP ibdFnameSEXP, SEXP NPixelsSEXP, SEXP NSEXP, SEXP offsetSEXP, SEXP dataTypeStringSEXP, SEXP continuousSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const char* >::type ibdFname(ibdFnameSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type NPixels(NPixelsSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type N(NSEXP);
    Rcpp::traits::input_parameter< uint64_t >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type dataTypeString(dataTypeStringSEXP);
    Rcpp::traits::input_parameter< bool >::type continuous(continuousSEXP);
    rcpp_result_gen = Rcpp::wrap(CimzMLBinReadIntensity(ibdFname, NPixels, N, offset, dataTypeString, continuous));
    return rcpp_result_gen;
END_RCPP
}
// CimzMLReadPeakList
Rcpp::List CimzMLReadPeakList(const char* ibdFname, Rcpp::List imzML_peakList_descriptor, unsigned int PixelID);
RcppExport SEXP _rMSI2_CimzMLReadPeakList(SEXP ibdFnameSEXP, SEXP imzML_peakList_descriptorSEXP, SEXP PixelIDSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const char* >::type ibdFname(ibdFnameSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type imzML_peakList_descriptor(imzML_peakList_descriptorSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type PixelID(PixelIDSEXP);
    rcpp_result_gen = Rcpp::wrap(CimzMLReadPeakList(ibdFname, imzML_peakList_descriptor, PixelID));
    return rcpp_result_gen;
END_RCPP
}
// overwriteIbdUUid
void overwriteIbdUUid(const char* ibdFname, Rcpp::String newUUID);
RcppExport SEXP _rMSI2_overwriteIbdUUid(SEXP ibdFnameSEXP, SEXP newUUIDSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const char* >::type ibdFname(ibdFnameSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type newUUID(newUUIDSEXP);
    overwriteIbdUUid(ibdFname, newUUID);
    return R_NilValue;
END_RCPP
}
// Cload_imzMLSpectra
Rcpp::NumericMatrix Cload_imzMLSpectra(Rcpp::List rMSIobj, Rcpp::IntegerVector pixelIDs, Rcpp::NumericVector commonMassAxis, unsigned int number_of_threads);
RcppExport SEXP _rMSI2_Cload_imzMLSpectra(SEXP rMSIobjSEXP, SEXP pixelIDsSEXP, SEXP commonMassAxisSEXP, SEXP number_of_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type rMSIobj(rMSIobjSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type pixelIDs(pixelIDsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type commonMassAxis(commonMassAxisSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type number_of_threads(number_of_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(Cload_imzMLSpectra(rMSIobj, pixelIDs, commonMassAxis, number_of_threads));
    return rcpp_result_gen;
END_RCPP
}
// CimzMLParse
List CimzMLParse(String xml_path);
RcppExport SEXP _rMSI2_CimzMLParse(SEXP xml_pathSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type xml_path(xml_pathSEXP);
    rcpp_result_gen = Rcpp::wrap(CimzMLParse(xml_path));
    return rcpp_result_gen;
END_RCPP
}
// CimzMLStore
bool CimzMLStore(String fname, List imgInfo, const char* mass_spectrometer_file_format);
RcppExport SEXP _rMSI2_CimzMLStore(SEXP fnameSEXP, SEXP imgInfoSEXP, SEXP mass_spectrometer_file_formatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type fname(fnameSEXP);
    Rcpp::traits::input_parameter< List >::type imgInfo(imgInfoSEXP);
    Rcpp::traits::input_parameter< const char* >::type mass_spectrometer_file_format(mass_spectrometer_file_formatSEXP);
    rcpp_result_gen = Rcpp::wrap(CimzMLStore(fname, imgInfo, mass_spectrometer_file_format));
    return rcpp_result_gen;
END_RCPP
}
// AlignSpectrumToReference
Rcpp::List AlignSpectrumToReference(NumericVector mass, NumericVector ref, NumericVector spectrumInterpolated, NumericVector massProcessedMode, NumericVector intensityProcessedMode, bool bilinear, double lagRefLow, double lagRefMid, double lagRefHigh, int iterations, double lagLimitppm, int fftOverSampling, double winSizeRelative);
RcppExport SEXP _rMSI2_AlignSpectrumToReference(SEXP massSEXP, SEXP refSEXP, SEXP spectrumInterpolatedSEXP, SEXP massProcessedModeSEXP, SEXP intensityProcessedModeSEXP, SEXP bilinearSEXP, SEXP lagRefLowSEXP, SEXP lagRefMidSEXP, SEXP lagRefHighSEXP, SEXP iterationsSEXP, SEXP lagLimitppmSEXP, SEXP fftOverSamplingSEXP, SEXP winSizeRelativeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type mass(massSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ref(refSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type spectrumInterpolated(spectrumInterpolatedSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type massProcessedMode(massProcessedModeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type intensityProcessedMode(intensityProcessedModeSEXP);
    Rcpp::traits::input_parameter< bool >::type bilinear(bilinearSEXP);
    Rcpp::traits::input_parameter< double >::type lagRefLow(lagRefLowSEXP);
    Rcpp::traits::input_parameter< double >::type lagRefMid(lagRefMidSEXP);
    Rcpp::traits::input_parameter< double >::type lagRefHigh(lagRefHighSEXP);
    Rcpp::traits::input_parameter< int >::type iterations(iterationsSEXP);
    Rcpp::traits::input_parameter< double >::type lagLimitppm(lagLimitppmSEXP);
    Rcpp::traits::input_parameter< int >::type fftOverSampling(fftOverSamplingSEXP);
    Rcpp::traits::input_parameter< double >::type winSizeRelative(winSizeRelativeSEXP);
    rcpp_result_gen = Rcpp::wrap(AlignSpectrumToReference(mass, ref, spectrumInterpolated, massProcessedMode, intensityProcessedMode, bilinear, lagRefLow, lagRefMid, lagRefHigh, iterations, lagLimitppm, fftOverSampling, winSizeRelative));
    return rcpp_result_gen;
END_RCPP
}
// MergeMassAxisAutoBinSize
List MergeMassAxisAutoBinSize(NumericVector mz1, NumericVector mz2);
RcppExport SEXP _rMSI2_MergeMassAxisAutoBinSize(SEXP mz1SEXP, SEXP mz2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type mz1(mz1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mz2(mz2SEXP);
    rcpp_result_gen = Rcpp::wrap(MergeMassAxisAutoBinSize(mz1, mz2));
    return rcpp_result_gen;
END_RCPP
}
// COverallAverageSpectrum
NumericVector COverallAverageSpectrum(Rcpp::List rMSIObj_list, int numOfThreads, double memoryPerThreadMB, Rcpp::NumericVector commonMassAxis, double minTIC, double maxTic);
RcppExport SEXP _rMSI2_COverallAverageSpectrum(SEXP rMSIObj_listSEXP, SEXP numOfThreadsSEXP, SEXP memoryPerThreadMBSEXP, SEXP commonMassAxisSEXP, SEXP minTICSEXP, SEXP maxTicSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type rMSIObj_list(rMSIObj_listSEXP);
    Rcpp::traits::input_parameter< int >::type numOfThreads(numOfThreadsSEXP);
    Rcpp::traits::input_parameter< double >::type memoryPerThreadMB(memoryPerThreadMBSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type commonMassAxis(commonMassAxisSEXP);
    Rcpp::traits::input_parameter< double >::type minTIC(minTICSEXP);
    Rcpp::traits::input_parameter< double >::type maxTic(maxTicSEXP);
    rcpp_result_gen = Rcpp::wrap(COverallAverageSpectrum(rMSIObj_list, numOfThreads, memoryPerThreadMB, commonMassAxis, minTIC, maxTic));
    return rcpp_result_gen;
END_RCPP
}
// CcommonMassAxis
List CcommonMassAxis(List rMSIObj_list, int numOfThreads, double memoryPerThreadMB);
RcppExport SEXP _rMSI2_CcommonMassAxis(SEXP rMSIObj_listSEXP, SEXP numOfThreadsSEXP, SEXP memoryPerThreadMBSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type rMSIObj_list(rMSIObj_listSEXP);
    Rcpp::traits::input_parameter< int >::type numOfThreads(numOfThreadsSEXP);
    Rcpp::traits::input_parameter< double >::type memoryPerThreadMB(memoryPerThreadMBSEXP);
    rcpp_result_gen = Rcpp::wrap(CcommonMassAxis(rMSIObj_list, numOfThreads, memoryPerThreadMB));
    return rcpp_result_gen;
END_RCPP
}
// CRunFillPeaks
void CRunFillPeaks(Rcpp::List rMSIObj_list, int numOfThreads, double memoryPerThreadMB, Rcpp::Reference preProcessingParams, Rcpp::NumericVector commonMassAxis, Rcpp::List peakMatrix);
RcppExport SEXP _rMSI2_CRunFillPeaks(SEXP rMSIObj_listSEXP, SEXP numOfThreadsSEXP, SEXP memoryPerThreadMBSEXP, SEXP preProcessingParamsSEXP, SEXP commonMassAxisSEXP, SEXP peakMatrixSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type rMSIObj_list(rMSIObj_listSEXP);
    Rcpp::traits::input_parameter< int >::type numOfThreads(numOfThreadsSEXP);
    Rcpp::traits::input_parameter< double >::type memoryPerThreadMB(memoryPerThreadMBSEXP);
    Rcpp::traits::input_parameter< Rcpp::Reference >::type preProcessingParams(preProcessingParamsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type commonMassAxis(commonMassAxisSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type peakMatrix(peakMatrixSEXP);
    CRunFillPeaks(rMSIObj_list, numOfThreads, memoryPerThreadMB, preProcessingParams, commonMassAxis, peakMatrix);
    return R_NilValue;
END_RCPP
}
// CInternalReferenceSpectrum
List CInternalReferenceSpectrum(Rcpp::List rMSIObj_list, int numOfThreads, double memoryPerThreadMB, Rcpp::NumericVector referenceSpectrum, Rcpp::NumericVector commonMassAxis);
RcppExport SEXP _rMSI2_CInternalReferenceSpectrum(SEXP rMSIObj_listSEXP, SEXP numOfThreadsSEXP, SEXP memoryPerThreadMBSEXP, SEXP referenceSpectrumSEXP, SEXP commonMassAxisSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type rMSIObj_list(rMSIObj_listSEXP);
    Rcpp::traits::input_parameter< int >::type numOfThreads(numOfThreadsSEXP);
    Rcpp::traits::input_parameter< double >::type memoryPerThreadMB(memoryPerThreadMBSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type referenceSpectrum(referenceSpectrumSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type commonMassAxis(commonMassAxisSEXP);
    rcpp_result_gen = Rcpp::wrap(CInternalReferenceSpectrum(rMSIObj_list, numOfThreads, memoryPerThreadMB, referenceSpectrum, commonMassAxis));
    return rcpp_result_gen;
END_RCPP
}
// CRunPeakPicking
Rcpp::List CRunPeakPicking(Rcpp::List rMSIObj_list, int numOfThreads, double memoryPerThreadMB, Rcpp::Reference preProcessingParams, Rcpp::StringVector uuid, Rcpp::String outputDataPath, Rcpp::StringVector imzMLoutFnames, Rcpp::NumericVector commonMassAxis);
RcppExport SEXP _rMSI2_CRunPeakPicking(SEXP rMSIObj_listSEXP, SEXP numOfThreadsSEXP, SEXP memoryPerThreadMBSEXP, SEXP preProcessingParamsSEXP, SEXP uuidSEXP, SEXP outputDataPathSEXP, SEXP imzMLoutFnamesSEXP, SEXP commonMassAxisSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type rMSIObj_list(rMSIObj_listSEXP);
    Rcpp::traits::input_parameter< int >::type numOfThreads(numOfThreadsSEXP);
    Rcpp::traits::input_parameter< double >::type memoryPerThreadMB(memoryPerThreadMBSEXP);
    Rcpp::traits::input_parameter< Rcpp::Reference >::type preProcessingParams(preProcessingParamsSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type uuid(uuidSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type outputDataPath(outputDataPathSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type imzMLoutFnames(imzMLoutFnamesSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type commonMassAxis(commonMassAxisSEXP);
    rcpp_result_gen = Rcpp::wrap(CRunPeakPicking(rMSIObj_list, numOfThreads, memoryPerThreadMB, preProcessingParams, uuid, outputDataPath, imzMLoutFnames, commonMassAxis));
    return rcpp_result_gen;
END_RCPP
}
// CRunPreProcessing
List CRunPreProcessing(Rcpp::List rMSIObj_list, int numOfThreads, double memoryPerThreadMB, Rcpp::Reference preProcessingParams, Rcpp::NumericVector reference, Rcpp::StringVector uuid, Rcpp::String outputDataPath, Rcpp::StringVector imzMLoutFnames, Rcpp::NumericVector commonMassAxis);
RcppExport SEXP _rMSI2_CRunPreProcessing(SEXP rMSIObj_listSEXP, SEXP numOfThreadsSEXP, SEXP memoryPerThreadMBSEXP, SEXP preProcessingParamsSEXP, SEXP referenceSEXP, SEXP uuidSEXP, SEXP outputDataPathSEXP, SEXP imzMLoutFnamesSEXP, SEXP commonMassAxisSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type rMSIObj_list(rMSIObj_listSEXP);
    Rcpp::traits::input_parameter< int >::type numOfThreads(numOfThreadsSEXP);
    Rcpp::traits::input_parameter< double >::type memoryPerThreadMB(memoryPerThreadMBSEXP);
    Rcpp::traits::input_parameter< Rcpp::Reference >::type preProcessingParams(preProcessingParamsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type reference(referenceSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type uuid(uuidSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type outputDataPath(outputDataPathSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type imzMLoutFnames(imzMLoutFnamesSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type commonMassAxis(commonMassAxisSEXP);
    rcpp_result_gen = Rcpp::wrap(CRunPreProcessing(rMSIObj_list, numOfThreads, memoryPerThreadMB, preProcessingParams, reference, uuid, outputDataPath, imzMLoutFnames, commonMassAxis));
    return rcpp_result_gen;
END_RCPP
}
// NoiseEstimationFFTCosWin
NumericVector NoiseEstimationFFTCosWin(NumericVector x, int filWinSize);
RcppExport SEXP _rMSI2_NoiseEstimationFFTCosWin(SEXP xSEXP, SEXP filWinSizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type filWinSize(filWinSizeSEXP);
    rcpp_result_gen = Rcpp::wrap(NoiseEstimationFFTCosWin(x, filWinSize));
    return rcpp_result_gen;
END_RCPP
}
// NoiseEstimationFFTExpWin
NumericVector NoiseEstimationFFTExpWin(NumericVector x, int filWinSize);
RcppExport SEXP _rMSI2_NoiseEstimationFFTExpWin(SEXP xSEXP, SEXP filWinSizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type filWinSize(filWinSizeSEXP);
    rcpp_result_gen = Rcpp::wrap(NoiseEstimationFFTExpWin(x, filWinSize));
    return rcpp_result_gen;
END_RCPP
}
// NoiseEstimationFFTCosWinMat
NumericMatrix NoiseEstimationFFTCosWinMat(NumericMatrix x, int filWinSize);
RcppExport SEXP _rMSI2_NoiseEstimationFFTCosWinMat(SEXP xSEXP, SEXP filWinSizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type filWinSize(filWinSizeSEXP);
    rcpp_result_gen = Rcpp::wrap(NoiseEstimationFFTCosWinMat(x, filWinSize));
    return rcpp_result_gen;
END_RCPP
}
// NoiseEstimationFFTExpWinMat
NumericMatrix NoiseEstimationFFTExpWinMat(NumericMatrix x, int filWinSize);
RcppExport SEXP _rMSI2_NoiseEstimationFFTExpWinMat(SEXP xSEXP, SEXP filWinSizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type filWinSize(filWinSizeSEXP);
    rcpp_result_gen = Rcpp::wrap(NoiseEstimationFFTExpWinMat(x, filWinSize));
    return rcpp_result_gen;
END_RCPP
}
// C_adductAnnotation
Rcpp::List C_adductAnnotation(int numMonoiso, int numAdducts, int tolerance, int numMass, NumericVector R_monoisitopeMassVector, NumericVector R_adductMassVector, List R_isotopes, NumericVector R_isotopeListOrder, NumericVector R_massAxis, NumericMatrix R_peakMatrix, int numPixels, NumericVector R_labelAxis, NumericVector R_monoisotopicIndexVector);
RcppExport SEXP _rMSI2_C_adductAnnotation(SEXP numMonoisoSEXP, SEXP numAdductsSEXP, SEXP toleranceSEXP, SEXP numMassSEXP, SEXP R_monoisitopeMassVectorSEXP, SEXP R_adductMassVectorSEXP, SEXP R_isotopesSEXP, SEXP R_isotopeListOrderSEXP, SEXP R_massAxisSEXP, SEXP R_peakMatrixSEXP, SEXP numPixelsSEXP, SEXP R_labelAxisSEXP, SEXP R_monoisotopicIndexVectorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type numMonoiso(numMonoisoSEXP);
    Rcpp::traits::input_parameter< int >::type numAdducts(numAdductsSEXP);
    Rcpp::traits::input_parameter< int >::type tolerance(toleranceSEXP);
    Rcpp::traits::input_parameter< int >::type numMass(numMassSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type R_monoisitopeMassVector(R_monoisitopeMassVectorSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type R_adductMassVector(R_adductMassVectorSEXP);
    Rcpp::traits::input_parameter< List >::type R_isotopes(R_isotopesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type R_isotopeListOrder(R_isotopeListOrderSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type R_massAxis(R_massAxisSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type R_peakMatrix(R_peakMatrixSEXP);
    Rcpp::traits::input_parameter< int >::type numPixels(numPixelsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type R_labelAxis(R_labelAxisSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type R_monoisotopicIndexVector(R_monoisotopicIndexVectorSEXP);
    rcpp_result_gen = Rcpp::wrap(C_adductAnnotation(numMonoiso, numAdducts, tolerance, numMass, R_monoisitopeMassVector, R_adductMassVector, R_isotopes, R_isotopeListOrder, R_massAxis, R_peakMatrix, numPixels, R_labelAxis, R_monoisotopicIndexVector));
    return rcpp_result_gen;
END_RCPP
}
// C_isotopeAnnotator
Rcpp::List C_isotopeAnnotator(int massPeaks, int massChannels, int numPixels, int numIso, NumericMatrix PeakMtx, NumericVector massVec, NumericVector massChanVec, int tolerance, double scoreThreshold, bool ToleranceInScans, int charge);
RcppExport SEXP _rMSI2_C_isotopeAnnotator(SEXP massPeaksSEXP, SEXP massChannelsSEXP, SEXP numPixelsSEXP, SEXP numIsoSEXP, SEXP PeakMtxSEXP, SEXP massVecSEXP, SEXP massChanVecSEXP, SEXP toleranceSEXP, SEXP scoreThresholdSEXP, SEXP ToleranceInScansSEXP, SEXP chargeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type massPeaks(massPeaksSEXP);
    Rcpp::traits::input_parameter< int >::type massChannels(massChannelsSEXP);
    Rcpp::traits::input_parameter< int >::type numPixels(numPixelsSEXP);
    Rcpp::traits::input_parameter< int >::type numIso(numIsoSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type PeakMtx(PeakMtxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type massVec(massVecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type massChanVec(massChanVecSEXP);
    Rcpp::traits::input_parameter< int >::type tolerance(toleranceSEXP);
    Rcpp::traits::input_parameter< double >::type scoreThreshold(scoreThresholdSEXP);
    Rcpp::traits::input_parameter< bool >::type ToleranceInScans(ToleranceInScansSEXP);
    Rcpp::traits::input_parameter< int >::type charge(chargeSEXP);
    rcpp_result_gen = Rcpp::wrap(C_isotopeAnnotator(massPeaks, massChannels, numPixels, numIso, PeakMtx, massVec, massChanVec, tolerance, scoreThreshold, ToleranceInScans, charge));
    return rcpp_result_gen;
END_RCPP
}
// CRunPeakBinning
List CRunPeakBinning(Rcpp::List rMSIObj_list, int numOfThreads, double memoryPerThreadMB, Rcpp::Reference preProcessingParams);
RcppExport SEXP _rMSI2_CRunPeakBinning(SEXP rMSIObj_listSEXP, SEXP numOfThreadsSEXP, SEXP memoryPerThreadMBSEXP, SEXP preProcessingParamsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type rMSIObj_list(rMSIObj_listSEXP);
    Rcpp::traits::input_parameter< int >::type numOfThreads(numOfThreadsSEXP);
    Rcpp::traits::input_parameter< double >::type memoryPerThreadMB(memoryPerThreadMBSEXP);
    Rcpp::traits::input_parameter< Rcpp::Reference >::type preProcessingParams(preProcessingParamsSEXP);
    rcpp_result_gen = Rcpp::wrap(CRunPeakBinning(rMSIObj_list, numOfThreads, memoryPerThreadMB, preProcessingParams));
    return rcpp_result_gen;
END_RCPP
}
// DetectPeaks_C
NumericMatrix DetectPeaks_C(NumericVector mass, NumericVector intensity, double SNR, int WinSize, int UpSampling);
RcppExport SEXP _rMSI2_DetectPeaks_C(SEXP massSEXP, SEXP intensitySEXP, SEXP SNRSEXP, SEXP WinSizeSEXP, SEXP UpSamplingSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type mass(massSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type intensity(intensitySEXP);
    Rcpp::traits::input_parameter< double >::type SNR(SNRSEXP);
    Rcpp::traits::input_parameter< int >::type WinSize(WinSizeSEXP);
    Rcpp::traits::input_parameter< int >::type UpSampling(UpSamplingSEXP);
    rcpp_result_gen = Rcpp::wrap(DetectPeaks_C(mass, intensity, SNR, WinSize, UpSampling));
    return rcpp_result_gen;
END_RCPP
}
// TestPeakInterpolation_C
NumericVector TestPeakInterpolation_C(NumericVector mass, NumericVector intensity, int peakIndex, int WinSize, int UpSampling, bool useHanning, int Iterations);
RcppExport SEXP _rMSI2_TestPeakInterpolation_C(SEXP massSEXP, SEXP intensitySEXP, SEXP peakIndexSEXP, SEXP WinSizeSEXP, SEXP UpSamplingSEXP, SEXP useHanningSEXP, SEXP IterationsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type mass(massSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type intensity(intensitySEXP);
    Rcpp::traits::input_parameter< int >::type peakIndex(peakIndexSEXP);
    Rcpp::traits::input_parameter< int >::type WinSize(WinSizeSEXP);
    Rcpp::traits::input_parameter< int >::type UpSampling(UpSamplingSEXP);
    Rcpp::traits::input_parameter< bool >::type useHanning(useHanningSEXP);
    Rcpp::traits::input_parameter< int >::type Iterations(IterationsSEXP);
    rcpp_result_gen = Rcpp::wrap(TestPeakInterpolation_C(mass, intensity, peakIndex, WinSize, UpSampling, useHanning, Iterations));
    return rcpp_result_gen;
END_RCPP
}
// TestHanningWindow
NumericVector TestHanningWindow(NumericVector mass, int WinSize, int UpSampling);
RcppExport SEXP _rMSI2_TestHanningWindow(SEXP massSEXP, SEXP WinSizeSEXP, SEXP UpSamplingSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type mass(massSEXP);
    Rcpp::traits::input_parameter< int >::type WinSize(WinSizeSEXP);
    Rcpp::traits::input_parameter< int >::type UpSampling(UpSamplingSEXP);
    rcpp_result_gen = Rcpp::wrap(TestHanningWindow(mass, WinSize, UpSampling));
    return rcpp_result_gen;
END_RCPP
}
// TestAreaWindow
NumericVector TestAreaWindow(NumericVector mass, int WinSize, int UpSampling);
RcppExport SEXP _rMSI2_TestAreaWindow(SEXP massSEXP, SEXP WinSizeSEXP, SEXP UpSamplingSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type mass(massSEXP);
    Rcpp::traits::input_parameter< int >::type WinSize(WinSizeSEXP);
    Rcpp::traits::input_parameter< int >::type UpSampling(UpSamplingSEXP);
    rcpp_result_gen = Rcpp::wrap(TestAreaWindow(mass, WinSize, UpSampling));
    return rcpp_result_gen;
END_RCPP
}
// ReduceDataPointsC
List ReduceDataPointsC(NumericVector mass, NumericVector intensity, double massMin, double massMax, int npoints);
RcppExport SEXP _rMSI2_ReduceDataPointsC(SEXP massSEXP, SEXP intensitySEXP, SEXP massMinSEXP, SEXP massMaxSEXP, SEXP npointsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type mass(massSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type intensity(intensitySEXP);
    Rcpp::traits::input_parameter< double >::type massMin(massMinSEXP);
    Rcpp::traits::input_parameter< double >::type massMax(massMaxSEXP);
    Rcpp::traits::input_parameter< int >::type npoints(npointsSEXP);
    rcpp_result_gen = Rcpp::wrap(ReduceDataPointsC(mass, intensity, massMin, massMax, npoints));
    return rcpp_result_gen;
END_RCPP
}
// Ccreate_rMSIXBinData
List Ccreate_rMSIXBinData(List rMSIobj, int number_of_threads);
RcppExport SEXP _rMSI2_Ccreate_rMSIXBinData(SEXP rMSIobjSEXP, SEXP number_of_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type rMSIobj(rMSIobjSEXP);
    Rcpp::traits::input_parameter< int >::type number_of_threads(number_of_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(Ccreate_rMSIXBinData(rMSIobj, number_of_threads));
    return rcpp_result_gen;
END_RCPP
}
// Cload_rMSIXBinData
List Cload_rMSIXBinData(String path, String fname);
RcppExport SEXP _rMSI2_Cload_rMSIXBinData(SEXP pathSEXP, SEXP fnameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type path(pathSEXP);
    Rcpp::traits::input_parameter< String >::type fname(fnameSEXP);
    rcpp_result_gen = Rcpp::wrap(Cload_rMSIXBinData(path, fname));
    return rcpp_result_gen;
END_RCPP
}
// Cload_rMSIXBinIonImage
NumericMatrix Cload_rMSIXBinIonImage(List rMSIobj, unsigned int ionIndex, unsigned int ionCount, NumericVector normalization_coefs, int number_of_threads);
RcppExport SEXP _rMSI2_Cload_rMSIXBinIonImage(SEXP rMSIobjSEXP, SEXP ionIndexSEXP, SEXP ionCountSEXP, SEXP normalization_coefsSEXP, SEXP number_of_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type rMSIobj(rMSIobjSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type ionIndex(ionIndexSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type ionCount(ionCountSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type normalization_coefs(normalization_coefsSEXP);
    Rcpp::traits::input_parameter< int >::type number_of_threads(number_of_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(Cload_rMSIXBinIonImage(rMSIobj, ionIndex, ionCount, normalization_coefs, number_of_threads));
    return rcpp_result_gen;
END_RCPP
}
// Smoothing_SavitzkyGolay
NumericVector Smoothing_SavitzkyGolay(NumericVector x, int sgSize);
RcppExport SEXP _rMSI2_Smoothing_SavitzkyGolay(SEXP xSEXP, SEXP sgSizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type sgSize(sgSizeSEXP);
    rcpp_result_gen = Rcpp::wrap(Smoothing_SavitzkyGolay(x, sgSize));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rMSI2_CNormalizationsAndMeans", (DL_FUNC) &_rMSI2_CNormalizationsAndMeans, 4},
    {"_rMSI2_CparseBrukerXML", (DL_FUNC) &_rMSI2_CparseBrukerXML, 1},
    {"_rMSI2_testingimzMLBinWriteSequential", (DL_FUNC) &_rMSI2_testingimzMLBinWriteSequential, 6},
    {"_rMSI2_CimzMLBinCreateNewIBD", (DL_FUNC) &_rMSI2_CimzMLBinCreateNewIBD, 2},
    {"_rMSI2_CimzMLBinAppendMass", (DL_FUNC) &_rMSI2_CimzMLBinAppendMass, 3},
    {"_rMSI2_CimzMLBinAppendIntensity", (DL_FUNC) &_rMSI2_CimzMLBinAppendIntensity, 3},
    {"_rMSI2_CimzMLBinWriteModifyMass", (DL_FUNC) &_rMSI2_CimzMLBinWriteModifyMass, 7},
    {"_rMSI2_CimzMLBinReadMass", (DL_FUNC) &_rMSI2_CimzMLBinReadMass, 6},
    {"_rMSI2_CimzMLBinReadIntensity", (DL_FUNC) &_rMSI2_CimzMLBinReadIntensity, 6},
    {"_rMSI2_CimzMLReadPeakList", (DL_FUNC) &_rMSI2_CimzMLReadPeakList, 3},
    {"_rMSI2_overwriteIbdUUid", (DL_FUNC) &_rMSI2_overwriteIbdUUid, 2},
    {"_rMSI2_Cload_imzMLSpectra", (DL_FUNC) &_rMSI2_Cload_imzMLSpectra, 4},
    {"_rMSI2_CimzMLParse", (DL_FUNC) &_rMSI2_CimzMLParse, 1},
    {"_rMSI2_CimzMLStore", (DL_FUNC) &_rMSI2_CimzMLStore, 3},
    {"_rMSI2_AlignSpectrumToReference", (DL_FUNC) &_rMSI2_AlignSpectrumToReference, 13},
    {"_rMSI2_MergeMassAxisAutoBinSize", (DL_FUNC) &_rMSI2_MergeMassAxisAutoBinSize, 2},
    {"_rMSI2_COverallAverageSpectrum", (DL_FUNC) &_rMSI2_COverallAverageSpectrum, 6},
    {"_rMSI2_CcommonMassAxis", (DL_FUNC) &_rMSI2_CcommonMassAxis, 3},
    {"_rMSI2_CRunFillPeaks", (DL_FUNC) &_rMSI2_CRunFillPeaks, 6},
    {"_rMSI2_CInternalReferenceSpectrum", (DL_FUNC) &_rMSI2_CInternalReferenceSpectrum, 5},
    {"_rMSI2_CRunPeakPicking", (DL_FUNC) &_rMSI2_CRunPeakPicking, 8},
    {"_rMSI2_CRunPreProcessing", (DL_FUNC) &_rMSI2_CRunPreProcessing, 9},
    {"_rMSI2_NoiseEstimationFFTCosWin", (DL_FUNC) &_rMSI2_NoiseEstimationFFTCosWin, 2},
    {"_rMSI2_NoiseEstimationFFTExpWin", (DL_FUNC) &_rMSI2_NoiseEstimationFFTExpWin, 2},
    {"_rMSI2_NoiseEstimationFFTCosWinMat", (DL_FUNC) &_rMSI2_NoiseEstimationFFTCosWinMat, 2},
    {"_rMSI2_NoiseEstimationFFTExpWinMat", (DL_FUNC) &_rMSI2_NoiseEstimationFFTExpWinMat, 2},
    {"_rMSI2_C_adductAnnotation", (DL_FUNC) &_rMSI2_C_adductAnnotation, 13},
    {"_rMSI2_C_isotopeAnnotator", (DL_FUNC) &_rMSI2_C_isotopeAnnotator, 11},
    {"_rMSI2_CRunPeakBinning", (DL_FUNC) &_rMSI2_CRunPeakBinning, 4},
    {"_rMSI2_DetectPeaks_C", (DL_FUNC) &_rMSI2_DetectPeaks_C, 5},
    {"_rMSI2_TestPeakInterpolation_C", (DL_FUNC) &_rMSI2_TestPeakInterpolation_C, 7},
    {"_rMSI2_TestHanningWindow", (DL_FUNC) &_rMSI2_TestHanningWindow, 3},
    {"_rMSI2_TestAreaWindow", (DL_FUNC) &_rMSI2_TestAreaWindow, 3},
    {"_rMSI2_ReduceDataPointsC", (DL_FUNC) &_rMSI2_ReduceDataPointsC, 5},
    {"_rMSI2_Ccreate_rMSIXBinData", (DL_FUNC) &_rMSI2_Ccreate_rMSIXBinData, 2},
    {"_rMSI2_Cload_rMSIXBinData", (DL_FUNC) &_rMSI2_Cload_rMSIXBinData, 2},
    {"_rMSI2_Cload_rMSIXBinIonImage", (DL_FUNC) &_rMSI2_Cload_rMSIXBinIonImage, 5},
    {"_rMSI2_Smoothing_SavitzkyGolay", (DL_FUNC) &_rMSI2_Smoothing_SavitzkyGolay, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_rMSI2(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
