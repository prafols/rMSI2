// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// CparseBrukerXML
List CparseBrukerXML(String xml_path);
RcppExport SEXP _rMSI_CparseBrukerXML(SEXP xml_pathSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type xml_path(xml_pathSEXP);
    rcpp_result_gen = Rcpp::wrap(CparseBrukerXML(xml_path));
    return rcpp_result_gen;
END_RCPP
}
// testingimzMLBinRead
Rcpp::NumericVector testingimzMLBinRead(const char* ibdFname, unsigned int NPixels, unsigned int N, unsigned long offset, Rcpp::String dataTypeString, bool read_mz, bool continuous);
RcppExport SEXP _rMSI_testingimzMLBinRead(SEXP ibdFnameSEXP, SEXP NPixelsSEXP, SEXP NSEXP, SEXP offsetSEXP, SEXP dataTypeStringSEXP, SEXP read_mzSEXP, SEXP continuousSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const char* >::type ibdFname(ibdFnameSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type NPixels(NPixelsSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type N(NSEXP);
    Rcpp::traits::input_parameter< unsigned long >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type dataTypeString(dataTypeStringSEXP);
    Rcpp::traits::input_parameter< bool >::type read_mz(read_mzSEXP);
    Rcpp::traits::input_parameter< bool >::type continuous(continuousSEXP);
    rcpp_result_gen = Rcpp::wrap(testingimzMLBinRead(ibdFname, NPixels, N, offset, dataTypeString, read_mz, continuous));
    return rcpp_result_gen;
END_RCPP
}
// CimzMLParse
List CimzMLParse(String xml_path);
RcppExport SEXP _rMSI_CimzMLParse(SEXP xml_pathSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type xml_path(xml_pathSEXP);
    rcpp_result_gen = Rcpp::wrap(CimzMLParse(xml_path));
    return rcpp_result_gen;
END_RCPP
}
// CimzMLStore
bool CimzMLStore(String fname, List imgInfo);
RcppExport SEXP _rMSI_CimzMLStore(SEXP fnameSEXP, SEXP imgInfoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type fname(fnameSEXP);
    Rcpp::traits::input_parameter< List >::type imgInfo(imgInfoSEXP);
    rcpp_result_gen = Rcpp::wrap(CimzMLStore(fname, imgInfo));
    return rcpp_result_gen;
END_RCPP
}
// CalcMassAxisBinSize
NumericVector CalcMassAxisBinSize(NumericVector mass, NumericVector intensity);
RcppExport SEXP _rMSI_CalcMassAxisBinSize(SEXP massSEXP, SEXP intensitySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type mass(massSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type intensity(intensitySEXP);
    rcpp_result_gen = Rcpp::wrap(CalcMassAxisBinSize(mass, intensity));
    return rcpp_result_gen;
END_RCPP
}
// MergeMassAxis
List MergeMassAxis(NumericVector mz1, NumericVector bins1, NumericVector mz2, NumericVector bins2);
RcppExport SEXP _rMSI_MergeMassAxis(SEXP mz1SEXP, SEXP bins1SEXP, SEXP mz2SEXP, SEXP bins2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type mz1(mz1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bins1(bins1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mz2(mz2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bins2(bins2SEXP);
    rcpp_result_gen = Rcpp::wrap(MergeMassAxis(mz1, bins1, mz2, bins2));
    return rcpp_result_gen;
END_RCPP
}
// ReduceDataPointsC
List ReduceDataPointsC(NumericVector mass, NumericVector intensity, double massMin, double massMax, int npoints);
RcppExport SEXP _rMSI_ReduceDataPointsC(SEXP massSEXP, SEXP intensitySEXP, SEXP massMinSEXP, SEXP massMaxSEXP, SEXP npointsSEXP) {
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
RcppExport SEXP _rMSI_Ccreate_rMSIXBinData(SEXP rMSIobjSEXP, SEXP number_of_threadsSEXP) {
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
RcppExport SEXP _rMSI_Cload_rMSIXBinData(SEXP pathSEXP, SEXP fnameSEXP) {
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
RcppExport SEXP _rMSI_Cload_rMSIXBinIonImage(SEXP rMSIobjSEXP, SEXP ionIndexSEXP, SEXP ionCountSEXP, SEXP normalization_coefsSEXP, SEXP number_of_threadsSEXP) {
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
// Cload_imzMLSpectra
NumericMatrix Cload_imzMLSpectra(List rMSIobj, IntegerVector pixelIDs);
RcppExport SEXP _rMSI_Cload_imzMLSpectra(SEXP rMSIobjSEXP, SEXP pixelIDsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type rMSIobj(rMSIobjSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type pixelIDs(pixelIDsSEXP);
    rcpp_result_gen = Rcpp::wrap(Cload_imzMLSpectra(rMSIobj, pixelIDs));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rMSI_CparseBrukerXML", (DL_FUNC) &_rMSI_CparseBrukerXML, 1},
    {"_rMSI_testingimzMLBinRead", (DL_FUNC) &_rMSI_testingimzMLBinRead, 7},
    {"_rMSI_CimzMLParse", (DL_FUNC) &_rMSI_CimzMLParse, 1},
    {"_rMSI_CimzMLStore", (DL_FUNC) &_rMSI_CimzMLStore, 2},
    {"_rMSI_CalcMassAxisBinSize", (DL_FUNC) &_rMSI_CalcMassAxisBinSize, 2},
    {"_rMSI_MergeMassAxis", (DL_FUNC) &_rMSI_MergeMassAxis, 4},
    {"_rMSI_ReduceDataPointsC", (DL_FUNC) &_rMSI_ReduceDataPointsC, 5},
    {"_rMSI_Ccreate_rMSIXBinData", (DL_FUNC) &_rMSI_Ccreate_rMSIXBinData, 2},
    {"_rMSI_Cload_rMSIXBinData", (DL_FUNC) &_rMSI_Cload_rMSIXBinData, 2},
    {"_rMSI_Cload_rMSIXBinIonImage", (DL_FUNC) &_rMSI_Cload_rMSIXBinIonImage, 5},
    {"_rMSI_Cload_imzMLSpectra", (DL_FUNC) &_rMSI_Cload_imzMLSpectra, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_rMSI(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
