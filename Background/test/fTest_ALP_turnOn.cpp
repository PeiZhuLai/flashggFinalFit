#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include "boost/program_options.hpp"
#include "boost/lexical_cast.hpp"
#include "TFile.h"
#include "TMath.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooHist.h"
#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooArgSet.h"
#include "RooFitResult.h"
#include "RooMinimizer.h"
#include "RooMsgService.h"
#include "RooDataHist.h"
#include "RooExtendPdf.h"
#include "TRandom3.h"
#include "TLatex.h"
#include "TMacro.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TArrow.h"
#include "TKey.h"
#include "RooCategory.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooMultiPdfCombine.h"
#include "../interface/PdfModelBuilder.h"
#include <Math/PdfFuncMathCore.h>
#include <Math/ProbFunc.h>
#include <iomanip>
#include "boost/program_options.hpp"
#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string/classification.hpp"
#include "boost/algorithm/string/predicate.hpp"
#include "../../tdrStyle/tdrstyle.C"
#include "../../tdrStyle/CMS_lumi.C"
#include "../interface/PlotStyle.h"
#include "TStyle.h"
#include "TColor.h"
#include <algorithm>
#include <cctype>
#include <cstring>
#include <cmath>

using namespace std;
using namespace RooFit;
using namespace boost;

namespace po = program_options;

static void transferMacros(TFile* in, TFile* out);
static void ensureOrderSuffix(RooAbsPdf* pdf, int order);
static int getBestFitFunction(RooMultiPdf *mpdf, RooAbsData *data, RooCategory *catIndex, bool silent);
static bool checkPdfDataObservables(RooAbsPdf* pdf, RooAbsData* data, bool verboseDiag = true);
static void syncMassRangeToData(RooRealVar* mass, RooAbsData* data, bool verboseDiag = true);
static void sanitizePdfParams(RooAbsPdf* pdf, RooAbsData* data, bool verboseDiag = false);
static void randomizePdfParamsUniform(RooAbsPdf* pdf, RooAbsData* data);

bool BLIND = true;
bool runFtestCheckWithToys=false;
bool PLOT_ONLY = false;

int FTEST_NTOYS = 500; // was 5000; lower = much faster (override with --ftoys)
int GOF_NTOYS   = 200; // was 500;  lower = faster (override with --gtoys)
int MIN_ENVELOPE_PDFS = 0; // No artificial floor: let envelope size be driven by the F-test
                           // (looser upperEnvThreshold) + one representative per family.
                           // Padding with sub-threshold "next-best" orders introduced jumps in
                           // the limits because the filler PDF was inconsistent across mA.

float mgglow_ =95.;//FIXME
float mgghigh_ =180;//FIXME
float mggblindlow_ =115;//FIXME
float mggblindhigh_ =135;//FIXME

float mgg_low =95.;//FIXME
float mgg_high =180.;//FIXME
float nBinsForMass = 1.*(mgg_high-mgg_low);
float mgg_blind_low =115;//FIXME
float mgg_blind_high =135;//FIXME

RooRealVar *intLumi_ = new RooRealVar("IntLumi","hacked int lumi", 1000.);

TRandom3 *RandomGen = new TRandom3();

// --- Helpers to decouple *plot* range from *fit* range ---
// We sometimes clamp the mass range to the dataset range for fit stability.
// But for plotting, we still want to show the user-requested window (--mhLow/--mhHigh).
static inline int binsForRange(double lo, double hi) {
  const double w = hi - lo;
  if (!(w > 0.0)) return 1;
  int nb = (int)std::lround(w);
  if (nb < 1) nb = 1;
  return nb;
}
struct MassStateGuard {
  RooRealVar* v = nullptr;
  double oldMin = 0.0;
  double oldMax = 0.0;
  int oldBins = 0;
  bool active = false;
  explicit MassStateGuard(RooRealVar* vv) : v(vv) {
    if (!v) return;
    oldMin = v->getMin();
    oldMax = v->getMax();
    oldBins = v->getBins();
    active = true;
  }
  ~MassStateGuard() {
    if (!active || !v) return;
    v->setRange(oldMin, oldMax);
    v->setBins(oldBins);
  }
};
static inline void setMassPlotState(RooRealVar* v) {
  if (!v) return;
  v->setRange((double)mgglow_, (double)mgghigh_);
  v->setBins(binsForRange((double)mgglow_, (double)mgghigh_));
}

struct EnvelopeCandidate {
  RooAbsPdf* pdf = nullptr;
  std::string family;
  std::string name;
  int order = -1;
  double gof = -1.0;
  double score = 1e300;
  bool isTruth = false;
  bool fromFallback = false;
  bool fromFamilyFallback = false;
};

static bool betterEnvelopeCandidate(const EnvelopeCandidate& lhs, const EnvelopeCandidate& rhs) {
  const bool lhsHasValidGof = std::isfinite(lhs.gof) && lhs.gof >= 0.0;
  const bool rhsHasValidGof = std::isfinite(rhs.gof) && rhs.gof >= 0.0;
  if (lhsHasValidGof != rhsHasValidGof) return lhsHasValidGof;
  if (lhsHasValidGof && rhsHasValidGof && std::fabs(lhs.gof - rhs.gof) > 1e-12) return lhs.gof > rhs.gof;
  if (std::fabs(lhs.score - rhs.score) > 1e-12) return lhs.score < rhs.score;
  if (lhs.family != rhs.family) return lhs.family < rhs.family;
  return lhs.order < rhs.order;
}

static bool lowerScoreEnvelopeCandidate(const EnvelopeCandidate& lhs, const EnvelopeCandidate& rhs) {
  if (std::fabs(lhs.score - rhs.score) > 1e-12) return lhs.score < rhs.score;
  const bool lhsHasValidGof = std::isfinite(lhs.gof) && lhs.gof >= 0.0;
  const bool rhsHasValidGof = std::isfinite(rhs.gof) && rhs.gof >= 0.0;
  if (lhsHasValidGof != rhsHasValidGof) return lhsHasValidGof;
  if (lhsHasValidGof && rhsHasValidGof && std::fabs(lhs.gof - rhs.gof) > 1e-12) return lhs.gof > rhs.gof;
  if (lhs.family != rhs.family) return lhs.family < rhs.family;
  return lhs.order < rhs.order;
}

static std::string envelopeCandidateKey(const EnvelopeCandidate& cand) {
  return cand.family + "::" + std::to_string(cand.order) + "::" + cand.name;
}

RooAbsPdf* getPdf(PdfModelBuilder &pdfsModel, string type, int order, const char* ext="", int mass_ALP=1)
{
  if (type=="Bernstein") return pdfsModel.getBernsteinStepxGau("Bern",order, mass_ALP);//PZ
  else if (type=="Chebychev") return pdfsModel.getChebychev("Che",order);
  else if (type=="Exponential") return pdfsModel.getExponentialStepxGau("Exp",order,2, mass_ALP);//PZ
  else if (type=="PowerLaw") return pdfsModel.getPowerLawStepxGau("Pow",order,2, mass_ALP);//PZ
  else if (type=="Laurent") return pdfsModel.getLaurentStepxGau("Lau",order,2, mass_ALP);//PZ

  else
  {
    cerr << "[ERROR] -- getPdf() -- type " << type << " not recognised." << endl;
    return NULL;
  }
}

// Exclusive upper bound on the order scanned per family (loop runs while order < this).
// Bias study (2026-06-13): with Bern4 in the envelope the bias study failed badly at high mass
// (Bern4 worst at 11 of mA4-30, up to -0.488 at mA27) because the 4th-order Bernstein is flexible
// enough to absorb signal. fTest never selects Bern4 anyway -> cap Bernstein at order 3 for ALL
// masses (envMaxOrder=4 -> scans orders 1,2,3). This also lets mA1-3 reach Bern3 (Bern2 under-fit
// gave mA3 = -0.213). At mA=1 the steep turn-on makes Pow3 (+0.52)/Exp3 (+0.30) badly biased ->
// keep exponential/power-law capped at order 2 there.
// Per-mA optimal Bernstein order (2026-06-14 bias study). For each mA the order
// is the one with the smallest |bias| among the GOF-passing (p>0.01) Bernstein
// orders, requiring the bias study to pass (|median|<0.2). Returns 0 where no
// Bernstein order passes -> Bernstein dropped from the envelope at that mA
// (mA 1: all functions fail GOF; mA 6: best Bern4=-0.244 still biases; mA 28: no
// GOF-passing Bernstein). The hung high-order points (12,21,24) land on a stable
// lower order that passes (B5/B4/B4), so this also avoids the Bern6 instability.
static int bestBernOrder(int m)
{
  // Per-mA Bernstein order chosen from ALL bias runs (Bern3/4/5/6 + validation v1):
  // among orders that pass GOF (p>0.01) AND have |median bias|<0.2 in EVERY available
  // measurement (robust against the Bernstein irreproducibility seen in single runs),
  // pick the one with the smallest max|bias| (largest margin). Returns 0 -> drop:
  //   mA 1, 28: no Bernstein order passes GOF;  mA 5, 6: no order robustly passes bias.
  switch (m) {
    case 2: return 2; case 3: return 3; case 4: return 5; case 9: return 6;
    case 10: return 5; case 11: return 5; case 12: return 5; case 13: return 4;
    case 14: return 5; case 15: return 5; case 16: return 6; case 17: return 5;
    case 18: return 5; case 19: return 5; case 20: return 6; case 22: return 6;
    case 23: return 6; case 24: return 5; case 25: return 6; case 26: return 6;
    case 27: return 6; case 29: return 6; case 30: return 6;
    default: return 0;  // mA 1, 5, 6, 7, 8, 21, 28: Bernstein dropped
  }
}

static int envMaxOrder(const std::string& funcType, int mass_ALP)
{
  // Bernstein: cap at the per-mA optimal order (exclusive upper bound = order+1).
  // 0 -> return 1 so the scan loop body never runs (Bernstein also excluded from
  // functionClasses below, so this is belt-and-braces).
  if (funcType=="Bernstein") {
    int o = bestBernOrder(mass_ALP);
    return (o > 0) ? (o + 1) : 1;
  }
  if (mass_ALP==1 && (funcType=="Exponential" || funcType=="PowerLaw")) return 3;  // excludes Exp3/Pow3
  return 7;
}

void runFit(RooAbsPdf *pdf, RooAbsData *data, double *NLL, int *stat_t, int MaxTries){
  if (PLOT_ONLY) {
    if (stat_t) *stat_t = 0;
    if (NLL) *NLL = 0.;
    return;
  }
  if (!pdf) {
    std::cerr << "[ERROR] runFit called with null pdf. Skip fit." << std::endl;
    if (stat_t) *stat_t = 5;
    if (NLL) *NLL = 1e12;
    return;
  }
  if (!data) {
    std::cerr << "[ERROR] runFit called with null data. Skip fit." << std::endl;
    if (stat_t) *stat_t = 5;
    if (NLL) *NLL = 1e12;
    return;
  }

  if (!checkPdfDataObservables(pdf, data, /*verboseDiag=*/true)) {
    std::cerr << "[ERROR] Observable mismatch between pdf and data. Skip fit to avoid RooFit abort." << std::endl;
    if (stat_t) *stat_t = 5;
    if (NLL) *NLL = 1e12;
    return;
  }


  // Keep a fixed fit window across all mass points. Do not shrink to the populated
  // data range, otherwise the effective lower edge becomes mass-dependent.
  // We only align the pdf/data observables to the user-requested window.
  RooRealVar* mData = dynamic_cast<RooRealVar*>(data && data->get() ? data->get()->find("CMS_hza_mass") : nullptr);
  std::unique_ptr<RooArgSet> _obsSet(pdf ? pdf->getObservables(*data) : nullptr);
  RooRealVar* mPdf = nullptr;
  if (_obsSet) mPdf = dynamic_cast<RooRealVar*>(_obsSet->find("CMS_hza_mass"));
  if (!mPdf && _obsSet && mData) mPdf = dynamic_cast<RooRealVar*>(_obsSet->find(mData->GetName()));
  if (!mPdf) mPdf = mData;

  struct RangeGuard {
    RooRealVar* v = nullptr;
    double oldMin = 0.0;
    double oldMax = 0.0;
    bool active = false;
    explicit RangeGuard(RooRealVar* vv) : v(vv) {}
    void clamp(double lo, double hi) {
      if (!v) return;
      oldMin = v->getMin();
      oldMax = v->getMax();
      v->setRange(lo, hi);
      active = true;
    }
    ~RangeGuard() {
      if (active && v) v->setRange(oldMin, oldMax);
    }
  } _rgPdf(mPdf), _rgData((mData && mData!=mPdf) ? mData : nullptr);

  if (mData && mPdf) {
    const double lo = std::max<double>((double)mgglow_, std::max(mPdf->getMin(), mData->getMin()));
    const double hi = std::min<double>((double)mgghigh_, std::min(mPdf->getMax(), mData->getMax()));
    if (hi > lo) {
      _rgPdf.clamp(lo, hi);
      if (mData!=mPdf) _rgData.clamp(lo, hi);
    }
  }

  int tries = 0;
  int status = 1;
  double bestNll = 1e12;

  const bool fastMode = (MaxTries<=1);

  sanitizePdfParams(pdf, data, /*verboseDiag=*/false);

  while (status!=0 && tries<MaxTries) {
    std::unique_ptr<RooAbsReal> nll(pdf->createNLL(
      *data,
      RooFit::Offset(true),
      RooFit::Optimize(2)
      // RooFit::PrintEvalErrors(0)   // 只出摘要(每個component計數)，不逐條狂刷
      // RooFit::PrintEvalErrors(-1) // 若你想完全不印
    ));
    RooMinimizer minim(*nll);
    minim.setPrintLevel(-1);
    minim.setStrategy(fastMode ? 0 : 1);
    minim.setOffsetting(true);
    minim.optimizeConst(1);
    minim.setEps(fastMode ? 1e-4 : 1e-6);
    minim.setMaxFunctionCalls(fastMode ? 2000 : 20000);
    minim.setMaxIterations(fastMode ? 2000 : 20000);

    status = minim.minimize("Minuit2","migrad");

    // Light fallback for real-data fits: try higher strategy once.
    if (!fastMode && status!=0) {
      minim.setStrategy(2);
      status = minim.minimize("Minuit2","migrad");
    }

    // if (status!=0) {
    //   minim.setStrategy(2);
    //   status = minim.minimize("Minuit2","minimize");
    // }
    // if (status!=0) {
    //   status = minim.minimize("Minuit2","simplex");
    //   if (status==0) {
    //     status = minim.minimize("Minuit2","migrad");
    //   }
    // }

    std::unique_ptr<RooFitResult> res(minim.save());
    double nllVal = nll->getVal();
    if (nllVal < bestNll) bestNll = nllVal;

    if (!fastMode && status!=0) {
      // If the fit failed, restart from a *sane* random point inside allowed parameter ranges.
      // (RooFitResult::randomizePars() relies on a valid covariance and may yield NaNs when the fit is unstable.)
      randomizePdfParamsUniform(pdf, data);
      sanitizePdfParams(pdf, data, /*verboseDiag=*/false);
    }
    tries++;
  }

  if (stat_t) *stat_t = status;
  if (NLL) *NLL = bestNll;
}

double getProbabilityFtest(double chi2, int ndof,RooAbsPdf *pdfNull, RooAbsPdf *pdfTest, RooRealVar *mass, RooAbsData *data, std::string name){
  if (PLOT_ONLY) return TMath::Prob(chi2,ndof);

  double prob_asym = TMath::Prob(chi2,ndof);
  if (!data || !pdfNull || !pdfTest || !mass) {
    std::cerr << "[WARN] getProbabilityFtest called with null input. Return asymptotic prob." << std::endl;
    return prob_asym;
  }

  syncMassRangeToData(mass, data, /*verboseDiag=*/true);

  if (!checkPdfDataObservables(pdfNull, data, /*verboseDiag=*/true) ||
      !checkPdfDataObservables(pdfTest, data, /*verboseDiag=*/true)) {
    std::cerr << "[WARN] getProbabilityFtest: observable mismatch. Return asymptotic prob (no-toys)." << std::endl;
    return prob_asym;
  }

  if (!runFtestCheckWithToys) return prob_asym;

  int ndata = data->sumEntries();

  // Speed/stability: for F-test (esp. with toys), fit a binned clone of the data.
  mass->setBins(nBinsForMass);
  // Some ROOT versions do NOT provide RooAbsData::binnedClone().
  // Only RooDataSet has binnedClone(), so do it conditionally.
  std::unique_ptr<RooAbsData> data_binned;
  RooAbsData* data_for_fit = data;
  if (auto* ds = dynamic_cast<RooDataSet*>(data)) {
    data_binned.reset(ds->binnedClone());
    if (data_binned) data_for_fit = data_binned.get();
  }

  // Ensure we keep a single pointer used for the fits.
  data_for_fit = data_binned ? data_binned.get() : data;

  double nllNullData = 1e12, nllTestData = 1e12;
  int statNullData = 1, statTestData = 1;
  runFit(pdfNull, data_for_fit, &nllNullData, &statNullData, /*MaxTries=*/3);
  runFit(pdfTest, data_for_fit, &nllTestData, &statTestData, /*MaxTries=*/3);

  if (statNullData != 0 || statTestData != 0) {
    std::cerr << "[WARN] getProbabilityFtest: fit to data failed (null=" << statNullData
              << ", test=" << statTestData << "). Return asymptotic prob." << std::endl;
    return prob_asym;
  }

  RooArgSet *params_null = pdfNull->getParameters((const RooArgSet*)(0));
  RooArgSet preParams_null;
  params_null->snapshot(preParams_null);
  RooArgSet *params_test = pdfTest->getParameters((const RooArgSet*)(0));
  RooArgSet preParams_test;
  params_test->snapshot(preParams_test);

  int ntoys = FTEST_NTOYS;
  TCanvas *can = new TCanvas();
  can->SetLogy();
  TH1F toyhist(Form("toys_fTest_%s.pdf",pdfNull->GetName()),";Chi2;",60,-2,10);
  TH1I toyhistStatN(Form("Status_%s.pdf",pdfNull->GetName()),";FitStatus;",8,-4,4);
  TH1I toyhistStatT(Form("Status_%s.pdf",pdfTest->GetName()),";FitStatus;",8,-4,4);

  TGraph *gChi2 = new TGraph();
  gChi2->SetLineColor(kGreen+2);
  double w = toyhist.GetBinWidth(1);

  int ipoint=0;

  for (int b=0;b<toyhist.GetNbinsX();b++){
	double x = toyhist.GetBinCenter(b+1);
	if (x>0){
	  gChi2->SetPoint(ipoint,x,(ROOT::Math::chisquared_pdf(x,ndof)));
	  ipoint++;
	}
  }
  int npass =0; int nsuccesst =0;
  mass->setBins(nBinsForMass);
  for (int itoy = 0 ; itoy < ntoys ; itoy++){

        params_null->assignValueOnly(preParams_null);
        params_test->assignValueOnly(preParams_test);
	RooDataHist *binnedtoy = pdfNull->generateBinned(RooArgSet(*mass),ndata,0,1);

    if (binnedtoy && binnedtoy->get()) {
      if (auto* xtoy = dynamic_cast<RooRealVar*>(binnedtoy->get()->find(mass->GetName()))) {
        xtoy->setMin(mass->getMin());
        xtoy->setMax(mass->getMax());
        xtoy->setRange(mass->getMin(), mass->getMax());
        xtoy->setBins(mass->getBins());
      }
    }

    int stat_n=1;
    int stat_t=1;
    int ntries = 0;
    double nllNull = 1e12, nllTest = 1e12;
    int MaxTries = 2;

    while (stat_n!=0){
      if (ntries>=MaxTries) break;
      runFit(pdfNull, binnedtoy, &nllNull, &stat_n, /*MaxTries=*/1);
      if (stat_n!=0) params_null->assignValueOnly(preParams_null);
      ntries++;
    }

    ntries = 0;
    while (stat_t!=0){
      if (ntries>=MaxTries) break;
      runFit(pdfTest, binnedtoy, &nllTest, &stat_t, /*MaxTries=*/1);
      if (stat_t!=0) params_test->assignValueOnly(preParams_test);
      ntries++;
    }

    toyhistStatN.Fill(stat_n);
    toyhistStatT.Fill(stat_t);

    if (stat_t !=0 || stat_n !=0) { delete binnedtoy; continue; }
    nsuccesst++;
    double chi2_t = 2*(nllNull-nllTest);
    if (chi2_t >= chi2) npass++;
    toyhist.Fill(chi2_t);

    delete binnedtoy;
  }

  double prob=0;
  if (nsuccesst!=0)  prob = (double)npass / nsuccesst;
  toyhist.Scale(1./(w*toyhist.Integral()));
  toyhist.Draw();
  TArrow lData(chi2,toyhist.GetMaximum(),chi2,0);
  lData.SetLineWidth(2);
  lData.Draw();
  gChi2->Draw("L");
  TLatex *lat = new TLatex();
  lat->SetNDC();
  lat->SetTextFont(42);
  lat->DrawLatex(0.1,0.91,Form("Prob (asymptotic) = %.4f (%.4f)",prob,prob_asym));
  can->SaveAs(name.c_str());

  TCanvas *stas =new TCanvas();
  toyhistStatN.SetLineColor(2);
  toyhistStatT.SetLineColor(1);
  TLegend *leg = new TLegend(0.2,0.6,0.4,0.87); leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->AddEntry(&toyhistStatN,"Null Hyp","L");
  leg->AddEntry(&toyhistStatT,"Test Hyp","L");
  toyhistStatN.Draw();
  toyhistStatT.Draw("same");
  leg->Draw();
  stas->SaveAs(Form("%s_fitstatus.pdf",name.c_str()));

  params_null->assignValueOnly(preParams_null);
  params_test->assignValueOnly(preParams_test);

  delete can; delete stas;
  delete gChi2;
  delete leg;
  delete lat;

  return prob_asym;

}

double getGoodnessOfFit(RooRealVar *mass, RooAbsPdf *mpdf, RooAbsData *data, std::string name){
  if (PLOT_ONLY) return 1.0;

  if (!data || !mpdf || !mass) {
    std::cerr << "[WARN] getGoodnessOfFit called with null input. Return prob=0." << std::endl;
    return 0.;
  }

  syncMassRangeToData(mass, data, /*verboseDiag=*/false);

  double prob;
  int ntoys = GOF_NTOYS;
  name+="_gofTest.pdf";
  const double nData = data->sumEntries();
  const double nInit = (nData > 0.) ? nData : 1.0;
  const double nMax  = (nData > 0.) ? 10E6 : 10.0;
  RooRealVar norm("norm","norm", nInit, 1e-6, nMax);

  RooExtendPdf *pdf = new RooExtendPdf("ext","ext",*mpdf,norm);

  RooPlot *plot_chi2 = mass->frame();
  data->plotOn(plot_chi2,
              RooFit::Name("data"),
              Binning(nBinsForMass),
              RooFit::DataError(RooAbsData::SumW2));
  pdf->plotOn(plot_chi2, RooFit::Name("pdf"));
  int np = pdf->getParameters(*data)->getSize();

  double chi2 = plot_chi2->chiSquare("pdf","data",np);
  std::cout << "[INFO] Calculating GOF for pdf " << pdf->GetName() << ", using " <<np << " fitted parameters" <<std::endl;

  if ((double)data->sumEntries()/nBinsForMass < 5 ){

    std::cout << "[INFO] Running toys for GOF test " << std::endl;
    RooArgSet *params = pdf->getParameters(*data);
    RooArgSet preParams;
    params->snapshot(preParams);
    int ndata = data->sumEntries();

    int npass =0;
    std::vector<double> toy_chi2;
    for (int itoy = 0 ; itoy < ntoys ; itoy++){
      params->assignValueOnly(preParams);
      int nToyEvents = RandomGen->Poisson(ndata);
      RooDataHist *binnedtoy = pdf->generateBinned(RooArgSet(*mass),nToyEvents,0,1);

      if (binnedtoy && binnedtoy->get()) {
        if (auto* xtoy = dynamic_cast<RooRealVar*>(binnedtoy->get()->find(mass->GetName()))) {
          xtoy->setMin(mass->getMin());
          xtoy->setMax(mass->getMax());
          xtoy->setRange(mass->getMin(), mass->getMax());
          xtoy->setBins(mass->getBins());
        }
      }

      double tmpNll = 1e12;
      int tmpStat = 1;
      runFit(pdf, binnedtoy, &tmpNll, &tmpStat, /*MaxTries=*/2);

      RooPlot *plot_t = mass->frame();
      binnedtoy->plotOn(plot_t, RooFit::Name("data"), RooFit::DataError(RooAbsData::SumW2));
      pdf->plotOn(plot_t, RooFit::Name("pdf"));

      double chi2_t = plot_t->chiSquare("pdf","data",np);
      if( chi2_t>=chi2) npass++;
      toy_chi2.push_back(chi2_t*(nBinsForMass-np));
      delete plot_t;
      delete binnedtoy;
    }
    std::cout << "[INFO] complete" << std::endl;
    prob = (double)npass / ntoys;

    TCanvas *can = new TCanvas();
    double medianChi2 = toy_chi2[(int)(((float)ntoys)/2)];
    double rms = TMath::Sqrt(medianChi2);

    TH1F toyhist(Form("gofTest_%s.pdf",pdf->GetName()),";Chi2;",50,medianChi2-5*rms,medianChi2+5*rms);
    for (std::vector<double>::iterator itx = toy_chi2.begin();itx!=toy_chi2.end();itx++){
      toyhist.Fill((*itx));
    }
    toyhist.Draw();

    TArrow lData(chi2*(nBinsForMass-np),toyhist.GetMaximum(),chi2*(nBinsForMass-np),0);
    lData.SetLineWidth(2);
    lData.Draw();
    can->SaveAs(name.c_str());

    params->assignValueOnly(preParams);
  } else {
    prob = TMath::Prob(chi2*(nBinsForMass-np),nBinsForMass-np);
  }
  std::cout << "[INFO] Chi2 in Observed =  " << chi2*(nBinsForMass-np) << std::endl;
  std::cout << "[INFO] p-value  =  " << prob << std::endl;
  delete pdf;
  return prob;

}

//Each Function Plot
void eachFunc_plot(RooRealVar *mass, RooAbsPdf *pdf, RooAbsData *data, string name,vector<string> flashggCats_, int status, double *prob){
  if (!pdf || !mass) return;
  if (!data) {
    std::cerr << "[WARN] eachFunc_plot(Each function pdf) called with null data. Skip plotting." << std::endl;
    if (prob) *prob = 0.;
    return;
  }

  RooPlot *plot_chi2 = mass->frame();
  data->plotOn(plot_chi2,
    Binning(nBinsForMass),
    RooFit::DataError(RooAbsData::SumW2),
    RooFit::Name("data"));

  pdf->plotOn(plot_chi2);

  int np = pdf->getParameters(*data)->getSize()+1;
  double chi2 = plot_chi2->chiSquare(np);

  // *prob = getGoodnessOfFit(mass,pdf,data,name);
  std::string probLabel = "N/A";
  if (prob) *prob = -1.;
  if (status == 0) {
    const double p = getGoodnessOfFit(mass, pdf, data, name);
    if (prob) *prob = p;
    probLabel = Form("%.2f", p);
  } else {
    std::cout << "[INFO] Skip GOF because fitStatus=" << status
              << " for pdf " << pdf->GetName() << std::endl;
  }
  
  // For output plots, force x-axis to the user-requested window (--mhLow/--mhHigh).
  // Keep the (possibly clamped) fit range for subsequent fits in the caller.
  MassStateGuard _fitState(mass);
  setMassPlotState(mass);
  const int plotBins = binsForRange((double)mgglow_, (double)mgghigh_);

  RooPlot *plot = mass->frame();
  mass->setRange("unblindReg_1", mgglow_,       mggblindlow_);
  mass->setRange("unblindReg_2", mggblindhigh_, mgghigh_);

  plot->GetXaxis()->SetTitle("m_{ll#gamma#gamma} (GeV)");
  plot->GetXaxis()->SetTitleSize(0.05);
  plot->GetXaxis()->SetLabelSize(0.04);
  plot->GetYaxis()->SetTitle("Events / 1 GeV");
  plot->GetYaxis()->SetTitleSize(0.05);
  plot->GetYaxis()->SetLabelSize(0.04);

  if (BLIND) {
    data->plotOn(plot,Binning(plotBins),CutRange("unblindReg_1"),RooFit::DataError(RooAbsData::SumW2));
    data->plotOn(plot,Binning(plotBins),CutRange("unblindReg_2"),RooFit::DataError(RooAbsData::SumW2));
  }
  else data->plotOn(plot,Binning(plotBins),RooFit::DataError(RooAbsData::SumW2));

  float leftMargion = 0.14;
  float bottomMargion = 0.13;
  TCanvas *canv = new TCanvas("","",800,800);
  canv->SetLeftMargin(leftMargion);
  canv->SetBottomMargin(bottomMargion);
  canv->SetRightMargin(0.05);
  canv->SetTopMargin(0.07);

  pdf->plotOn(plot);
  pdf->paramOn(plot,RooFit::Layout(0.14,0.96,0.89),RooFit::Format("NEA",AutoPrecision(1)));
  if (BLIND) plot->SetMinimum(0.0001);
  plot->SetTitle("");
  plot->GetYaxis()->SetTitleOffset(leftMargion*10.);
  plot->GetXaxis()->SetTitleOffset(bottomMargion*10.-0.05);

  plot->Draw();
  TLatex *lat = new TLatex();
  lat->SetNDC();
  lat->SetTextFont(42);
  lat->SetTextSize(0.04);
  lat->DrawLatex(0.14,0.94,Form("#chi^{2} = %.3f, Prob = %.2f, Fit Status = %d ",chi2*(nBinsForMass-np),*prob,status));
  canv->SaveAs(name.c_str());

  delete canv;
  delete lat;
}

// MultiPdf
void multipdf_plot(RooRealVar *mass, RooMultiPdf *pdfs, RooCategory *catIndex, RooAbsData *data, string name, vector<string> flashggCats_, int cat, int bestFitPdf=-1){
  if (!pdfs || !mass || !catIndex) return;
  if (!data) {
    std::cerr << "[WARN] multipdf_plot(multipdf) called with null data. Skip plotting." << std::endl;
    return;
  }

  // Save current (fit) mass state; we may have clamped it to the dataset range for stability.
  MassStateGuard _fitState(mass);
  const double fitMin  = mass->getMin();
  const double fitMax  = mass->getMax();
  const int    fitBins = mass->getBins();

  // For plotting, still show the user-requested window (--mhLow/--mhHigh).
  setMassPlotState(mass);
  const int plotBins = binsForRange((double)mgglow_, (double)mgghigh_);


  TCanvas *canv = new TCanvas("","",800,800);
  canv->SetLeftMargin(PlotStyleCfg::canvasLeftMargin);
  canv->SetRightMargin(PlotStyleCfg::canvasRightMargin);
  canv->SetTopMargin(PlotStyleCfg::canvasTopMargin);
  canv->SetBottomMargin(PlotStyleCfg::canvasBottomMargin);
  gStyle->SetOptStat(PlotStyleCfg::showStatBox ? 1 : 0);
  gStyle->SetPadTickX(PlotStyleCfg::axisTickX);
  gStyle->SetPadTickY(PlotStyleCfg::axisTickY);

  gStyle->SetLineScalePS(4.0);
  gStyle->SetEndErrorSize(0);
  gStyle->SetLineStyleString(2,"[16 12] 0");
  gStyle->SetLineStyleString(3,"[8 12] 0");
  gStyle->SetLineStyleString(4,"[24 12] 0");
  gStyle->SetLineStyleString(5,"[4 8] 0");

  TLegend *leg = new TLegend(PlotStyleCfg::multipdfLegendX1,
                             PlotStyleCfg::multipdfLegendY1,
                             PlotStyleCfg::multipdfLegendX2,
                             PlotStyleCfg::multipdfLegendY2);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.035);
  leg->SetTextFont(42);

  RooPlot *plot = mass->frame();

  mass->setRange("unblindReg_1", mgglow_,       mggblindlow_);
  mass->setRange("unblindReg_2", mggblindhigh_, mgghigh_);
  if (BLIND) {
    data->plotOn(plot,Binning(plotBins),CutRange("unblindReg_1"),RooFit::DataError(RooAbsData::SumW2));
    data->plotOn(plot,Binning(plotBins),CutRange("unblindReg_2"),RooFit::DataError(RooAbsData::SumW2));
    data->plotOn(plot,Binning(plotBins),Invisible());
  }
  else data->plotOn(plot,Binning(plotBins),RooFit::DataError(RooAbsData::SumW2));
  RooHist *plotdata = (RooHist*)plot->getObject(plot->numItems()-1);
  bool doRatioPlot_= PlotStyleCfg::enableRatio;
  TPad *pad1 = new TPad("pad1","pad1",0,PlotStyleCfg::pad2Height,1,1);
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,PlotStyleCfg::pad2Height);
  pad1->SetTopMargin(PlotStyleCfg::pad1TopMargin);
  pad1->SetBottomMargin(PlotStyleCfg::pad1BottomMargin);
  pad2->SetTopMargin(PlotStyleCfg::pad2TopMargin);
  pad2->SetBottomMargin(PlotStyleCfg::pad2BottomMargin);

  pad1->SetLeftMargin(PlotStyleCfg::canvasLeftMargin);
  pad1->SetRightMargin(PlotStyleCfg::canvasRightMargin);
  pad2->SetLeftMargin(PlotStyleCfg::canvasLeftMargin);
  pad2->SetRightMargin(PlotStyleCfg::canvasRightMargin);
  pad1->SetTicks(1,1);
  pad2->SetTicks(1,1);
  pad2->Draw();
  pad1->Draw();
  pad1->cd();

  int currentIndex = catIndex->getIndex();
  TObject *datLeg = plot->getObject(int(plot->numItems()-1));
  leg->AddEntry(datLeg,"Data","LEP");
  int style=1;
  RooAbsPdf *pdf = nullptr;
  RooCurve *nomBkgCurve = nullptr;
  int bestcol= -1;

  std::vector<RooCurve*> pdfCurves;

  for (int icat=0;icat<catIndex->numTypes();icat++){
    int col = PlotStyleCfg::colorForIndex(icat);
    if (col==kWhite || col==kYellow) col = kOrange+7;

    if (icat>6) { col=kBlack; style++; }
    catIndex->setIndex(icat);
    if (!PLOT_ONLY) {
      // Fit in the (possibly clamped) fit window to avoid RooFit range/data mismatches.
      mass->setRange(fitMin, fitMax);
      mass->setBins(fitBins);
      // [PZ-FIX] Use our lightweight RooMinimizer-based fit (no Hessian) to avoid Minuit2 warnings.
      runFit(pdfs->getCurrentPdf(), data, nullptr, nullptr, /*MaxTries=*/2);
    }
    // Plot in the user-requested window.
    setMassPlotState(mass);

    std::string curveName = Form("bkg_curve_%d",icat);
    pdfs->getCurrentPdf()->plotOn(
      plot,
      RooFit::Binning(plotBins),
      LineColor(col),
      LineStyle(style),
      LineWidth(4),
      Name(curveName.c_str())
    );

    RooCurve *thisCurve = dynamic_cast<RooCurve*>(plot->findObject(curveName.c_str()));
    if (thisCurve) {
      thisCurve->SetFillStyle(0);
      thisCurve->SetLineWidth(4);
      pdfCurves.push_back(thisCurve);
    }

    TObject *pdfLeg = thisCurve ? (TObject*)thisCurve : plot->getObject(int(plot->numItems()-1));
    std::string ext = "";
    if (bestFitPdf==icat) {
      ext=" (Best Fit) ";
      pdf= pdfs->getCurrentPdf();
      nomBkgCurve = thisCurve ? thisCurve : (RooCurve*)plot->getObject(plot->numItems()-1);
      bestcol = col;
    }
    leg->AddEntry(pdfLeg,Form("%s%s",pdfs->getCurrentPdf()->GetName(),ext.c_str()),"L");
  }
  if (!pdf) pdf = pdfs->getCurrentPdf();
  bool canDoRatio = (pdf && nomBkgCurve);

  plot->GetYaxis()->SetTitle("Events / 1 GeV");
  plot->GetYaxis()->SetTitleFont(PlotStyleCfg::axisTitleFont);
  plot->GetYaxis()->SetLabelFont(PlotStyleCfg::axisLabelFont);
  plot->GetYaxis()->SetTitleSize(PlotStyleCfg::pad1axisTitleSizeY);
  plot->GetYaxis()->SetLabelSize(PlotStyleCfg::pad1axisLabelSizeY);
  plot->GetYaxis()->SetTitleOffset(PlotStyleCfg::pad1axisTitleOffsetY);
  plot->GetXaxis()->SetTitle("m_{ll#gamma#gamma} (GeV)");


  plot->SetTitle(Form("Category %s",flashggCats_[cat].c_str()));
  if (BLIND) plot->SetMinimum(0.0001);
  plot->Draw();
  CMS_lumi(canv, 22, 0);

  leg->Draw("same");

  for (auto* c : pdfCurves) {
    if (!c) continue;
    c->SetFillStyle(0);
    c->SetLineWidth(3);
    c->Draw("L same");
  }

  pad1->RedrawAxis();

  int npoints = plotdata->GetN();
  double xtmp,ytmp;
  int point =0;
  TGraphAsymmErrors *hdatasub = new TGraphAsymmErrors(npoints);
  for (int ipoint=0; ipoint<npoints; ++ipoint) {
  plotdata->GetPoint(ipoint, xtmp,ytmp);
  double bkgval = 0.;
  if (canDoRatio) {
    bkgval = nomBkgCurve->interpolate(xtmp);
  } else {
    continue;
  }
  if (BLIND) {
   if ((xtmp > mggblindlow_) && (xtmp < mggblindhigh_)) continue;
  }
 double errhi = plotdata->GetErrorYhigh(ipoint);
 double errlow = plotdata->GetErrorYlow(ipoint);

 bool drawZeroBins_ =1;
 if (!drawZeroBins_) if(fabs(ytmp)<1e-5) continue;
 hdatasub->SetPoint(point,xtmp,ytmp-bkgval);
 hdatasub->SetPointError(point,0.,0.,errlow,errhi );
 point++;
  }
  pad2->cd();
  TH1 *hdummy = new TH1D("hdummyweight","",plotBins,mgglow_,mgghigh_);
  hdummy->SetMaximum(hdatasub->GetHistogram()->GetMaximum()+1);
  hdummy->SetMinimum(hdatasub->GetHistogram()->GetMinimum()-1);

  const double kScale = 1.0/PlotStyleCfg::pad2Height;
  hdummy->GetXaxis()->SetTitle("m_{ll#gamma#gamma} (GeV)");
  hdummy->GetYaxis()->SetTitle("Data - Best Fit");

  hdummy->GetXaxis()->SetTitleFont(PlotStyleCfg::axisTitleFont);
  hdummy->GetYaxis()->SetTitleFont(PlotStyleCfg::axisTitleFont);
  hdummy->GetXaxis()->SetLabelFont(PlotStyleCfg::axisLabelFont);
  hdummy->GetYaxis()->SetLabelFont(PlotStyleCfg::axisLabelFont);

  hdummy->GetYaxis()->SetTitleSize(PlotStyleCfg::pad2axisTitleSizeY);
  hdummy->GetYaxis()->SetLabelSize(PlotStyleCfg::pad2axisLabelSizeY);
  hdummy->GetYaxis()->SetTitleOffset(PlotStyleCfg::pad2axisTitleOffsetY);

  hdummy->GetXaxis()->SetTitleSize(PlotStyleCfg::pad2axisTitleSizeX);
  hdummy->GetXaxis()->SetLabelSize(PlotStyleCfg::pad2axisLabelSizeX);
  hdummy->GetXaxis()->SetTitleOffset(PlotStyleCfg::pad2axisTitleOffsetX);

  if (PlotStyleCfg::axisCenterTitle) {
    hdummy->GetXaxis()->CenterTitle(true);
    hdummy->GetYaxis()->CenterTitle(true);
  }
  hdummy->GetYaxis()->SetNdivisions(808);
  hdummy->Draw("HIST");

  if (canDoRatio) {
    TLine *line3 = new TLine(mgglow_,0.,mgghigh_,0.);
    line3->SetLineColor(bestcol);
    line3->SetLineWidth(PlotStyleCfg::zeroLineWidth);
    line3->Draw();
    hdatasub->Draw("PESAME");
  }

  pad1->Modified(); pad1->Update();
  pad2->Modified(); pad2->Update();
  canv->Modified(); canv->Update();

  canv->SaveAs(Form("%s.pdf",name.c_str()));
  canv->SaveAs(Form("%s.png",name.c_str()));
  catIndex->setIndex(currentIndex);
  delete canv;
}

void truth_plot(RooRealVar *mass,
  std::map<std::string,RooAbsPdf*> pdfs,
  RooAbsData *data,std::string name,
  std::vector<std::string> flashggCats_,
  int cat,int bestFitPdf = -1){
  if (!mass) return;
  if (!data) {
    std::cerr << "[WARN] truth_plot(map) called with null data. Skip plotting." << std::endl;
    return;
  }

  // Force plotting window to the user-requested range (--mhLow/--mhHigh).
  MassStateGuard _fitState(mass);
  setMassPlotState(mass);
  const int plotBins = binsForRange((double)mgglow_, (double)mgghigh_);


  TCanvas *canv = new TCanvas("","",800,800);
  canv->SetLeftMargin(PlotStyleCfg::canvasLeftMargin);
  canv->SetRightMargin(PlotStyleCfg::canvasRightMargin);
  canv->SetTopMargin(0.07);
  canv->SetBottomMargin(PlotStyleCfg::canvasBottomMargin);

  gStyle->SetOptStat(PlotStyleCfg::showStatBox ? 1 : 0);
  gStyle->SetPadTickX(PlotStyleCfg::axisTickX);
  gStyle->SetPadTickY(PlotStyleCfg::axisTickY);
  gStyle->SetLineScalePS(3.0);
  gStyle->SetEndErrorSize(0);
  gStyle->SetLineStyleString(2,"[16 12] 0");
  gStyle->SetLineStyleString(3,"[8 12] 0");
  gStyle->SetLineStyleString(4,"[24 12] 0");
  gStyle->SetLineStyleString(5,"[4 8] 0");

  TLegend *leg = new TLegend(PlotStyleCfg::truthLegendX1,
                PlotStyleCfg::truthLegendY1,
                PlotStyleCfg::truthLegendX2,
                PlotStyleCfg::truthLegendY2);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  leg->SetTextFont(42);

  // ---- RooPlot ----
  RooPlot *plot = mass->frame();

  mass->setRange("unblindReg_1", mgg_low,        mgg_blind_low);
  mass->setRange("unblindReg_2", mgg_blind_high, mgg_high);
  if (BLIND) {
  data->plotOn(plot, Binning(mgg_high - mgg_low), CutRange("unblindReg_1"), RooFit::DataError(RooAbsData::SumW2));
  data->plotOn(plot, Binning(mgg_high - mgg_low), CutRange("unblindReg_2"), RooFit::DataError(RooAbsData::SumW2));
  data->plotOn(plot, Binning(mgg_high - mgg_low), Invisible()); // 占位，便于取 legend 对象
  } else {
  data->plotOn(plot, Binning(mgg_high - mgg_low), RooFit::DataError(RooAbsData::SumW2));
  }

  TObject *datLeg = plot->getObject(int(plot->numItems()-1));
  if (datLeg) leg->AddEntry(datLeg, "Data", "LEP");

  int i = 0, style = 1;
  RooCurve* nomBkgCurve = nullptr;
  int bestcol = -1;
  std::vector<RooCurve*> pdfCurves;

  for (auto it = pdfs.begin(); it != pdfs.end(); ++it, ++i) {
  RooAbsPdf* p = it->second;
  if (!p) {
  std::cerr << "[WARN] truth_plot(map): skip null pdf entry " << it->first << std::endl;
  continue;
  }

  int col = PlotStyleCfg::colorForIndex(i);
  if (col == kWhite || col == kYellow) col = kOrange+7;
  if (i > 6) { col = kBlack; style++; }

  std::string curveName = Form("bkg_curve_%d", i);
  p->plotOn(plot,
  RooFit::Binning(plotBins),
  LineColor(col),
  LineStyle(style),
  LineWidth(4),
  Name(curveName.c_str()));

  RooCurve *thisCurve = dynamic_cast<RooCurve*>(plot->findObject(curveName.c_str()));
  if (thisCurve) {
  thisCurve->SetFillStyle(0);
  thisCurve->SetLineWidth(4);
  pdfCurves.push_back(thisCurve);
  }

  std::string ext = "";
  if (bestFitPdf == i) {
  ext = " (Best Fit) ";
  nomBkgCurve = thisCurve ? thisCurve : (RooCurve*)plot->getObject(plot->numItems()-1);
  bestcol = col;
  }
  TObject *pdfLeg = thisCurve ? (TObject*)thisCurve : plot->getObject(int(plot->numItems()-1));
  if (pdfLeg) leg->AddEntry(pdfLeg, Form("%s%s", it->first.c_str(), ext.c_str()), "L");
  }

  // ---- 坐标轴样式 ----
  plot->GetYaxis()->SetTitle("Events / 1 GeV");
  plot->GetYaxis()->SetTitleFont(PlotStyleCfg::axisTitleFont);
  plot->GetYaxis()->SetLabelFont(PlotStyleCfg::axisLabelFont);
  plot->GetYaxis()->SetTitleSize(0.05);
  plot->GetYaxis()->SetLabelSize(0.04);
  plot->GetYaxis()->SetTitleOffset(1.2);

  plot->GetXaxis()->SetTitle("m_{ll#gamma#gamma} (GeV)");
  plot->GetXaxis()->SetTitleFont(PlotStyleCfg::axisTitleFont);
  plot->GetXaxis()->SetLabelFont(PlotStyleCfg::axisLabelFont);
  plot->GetXaxis()->SetTitleSize(0.05);
  plot->GetXaxis()->SetLabelSize(0.04);
  plot->GetXaxis()->SetTitleOffset(1.2);

  plot->SetTitle(Form("Category %s", flashggCats_[cat].c_str()));
  if (BLIND) plot->SetMinimum(0.0001);

  // ---- 绘制到画布 ----
  plot->Draw();
  // CMS_lumi(canv, 22, 0);
  TLatex *lat = new TLatex();
  lat->SetNDC();
  lat->SetTextFont(42);
  lat->SetTextSize(0.045);
  lat->DrawLatex(PlotStyleCfg::canvasLeftMargin, 0.94, "#bf{CMS} #it{Preliminary}");
  // 減越多，字越靠近左邊
  lat->DrawLatex(1.-PlotStyleCfg::canvasRightMargin*11.-0.03, 0.94, "172.13 fb^{-1} (13.6 TeV)");

  // 确保曲线在最上层
  for (auto* c : pdfCurves) {
  if (!c) continue;
  c->SetFillStyle(0);
  c->SetLineWidth(3);
  c->Draw("L same");
  }

  leg->Draw("same");
  canv->RedrawAxis();

  // ---- 输出 ----
  canv->SaveAs(Form("%s.pdf", name.c_str()));
  canv->SaveAs(Form("%s.png", name.c_str()));

  // ---- 清理 ----
  delete leg;
  delete canv;
}

int main(int argc, char* argv[]){

  setTDRStyle();
  writeExtraText = true;
  extraText  = "Preliminary";
  lumi_8TeV  = "19.1 fb^{-1}";
  lumi_7TeV  = "4.9 fb^{-1}";
  lumi_sqrtS = "13.6 TeV";

  string fileName;
  string channelName;
  int ncats, mass_ALP;
  int singleCategory;
  string datfile;
  string outDir;
  string outfilename;
  bool is2011=false;
  bool verbose=false;
  bool saveMultiPdf=false;
	int isFlashgg_ =1;
  string flashggCatsStr_;
  vector<string> flashggCats_;
  bool isData_ =1;
  std::string datasetName_;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",                                                                                  "Show help")
    ("infilename,i", po::value<string>(&fileName),                                              "In file name")
    ("channel", po::value<string>(&channelName),                                                "channel name")
    ("ncats,c", po::value<int>(&ncats)->default_value(5),                                       "Number of categories")
    ("mass_ALP", po::value<int>(&mass_ALP)->default_value(1),                                   "Mass of ALP")
    ("singleCat", po::value<int>(&singleCategory)->default_value(-1),                           "Run A single Category")
    ("datfile,d", po::value<string>(&datfile)->default_value("dat/fTest.dat"),                  "Right results to datfile for BiasStudy")
    ("outDir,D", po::value<string>(&outDir)->default_value("plots/fTest"),                      "Out directory for plots")
    ("saveMultiPdf", po::value<string>(&outfilename),         					"Save a MultiPdf model with the appropriate pdfs")
    ("runFtestCheckWithToys", 									"When running the F-test, use toys to calculate pvals (and make plots) ")
    ("ftoys", po::value<int>(&FTEST_NTOYS)->default_value(500), "Number of toys for F-test p-values (default 500; was 5000)")
    ("gtoys", po::value<int>(&GOF_NTOYS)->default_value(200), "Number of toys for Goodness-of-Fit (default 200; was 500)")
    ("is2011",                                                                                  "Run 2011 config")
    ("is2012",                                                                                  "Run 2012 config")
    ("unblind",  									                                                              "Dont blind plots")
    ("isFlashgg",  po::value<int>(&isFlashgg_)->default_value(1),  								    	        "Use Flashgg output ")
    ("isData",  po::value<bool>(&isData_)->default_value(1),  								    	            "Use Data not MC ")
		("flashggCats,f", po::value<string>(&flashggCatsStr_)->default_value("UntaggedTag_0,UntaggedTag_1,UntaggedTag_2,UntaggedTag_3,UntaggedTag_4,VBFTag_0,VBFTag_1,VBFTag_2,TTHHadronicTag,TTHLeptonicTag,VHHadronicTag,VHTightTag,VHLooseTag,VHEtTag"),       "Flashgg category names to consider")
    ("datasetName", po::value<std::string>(&datasetName_)->default_value("Data_13p6TeV"),                                   "Override dataset name in workspace")
    ("plotOnly",                                                                         "Skip all fits; build minimal pdf set and only plot styles")
    ("verbose,v",                                                                               "Run with more output")
    ("mhLow,L", po::value<float>(&mgglow_)->default_value(95.),                                 "Low ALP mass point")
    ("mhHigh,H", po::value<float>(&mgghigh_)->default_value(180.),                              "High ALP mass point")
    ("mhLowBlind,LB", po::value<float>(&mggblindlow_)->default_value(115.),                     "Low ALP blind mass point")
    ("mhHighBlind,HB", po::value<float>(&mggblindhigh_)->default_value(135.),                   "High ALP blind mass point")
  ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc,argv,desc),vm);
  po::notify(vm);
  if (vm.count("help")) { cout << desc << endl; exit(1); }
  if (vm.count("is2011")) is2011=true;
	if (vm.count("unblind")) BLIND=false;
  saveMultiPdf = vm.count("saveMultiPdf");
  if (vm.count("plotOnly")) {
    PLOT_ONLY = true;
    runFtestCheckWithToys = false;
    if (verbose) std::cout << "[INFO] plotOnly mode enabled: all fits/toys skipped." << std::endl;
  }

  if (vm.count("verbose")) verbose=true;
  if (vm.count("runFtestCheckWithToys")) runFtestCheckWithToys=true;

  mgg_low =mgglow_;
  mgg_high =mgghigh_;
  nBinsForMass = 1.*(mgg_high-mgg_low);
  mgg_blind_low =mggblindlow_;
  mgg_blind_high =mggblindhigh_;

  if (!verbose) {
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    RooMsgService::instance().setSilentMode(false);
    gErrorIgnoreLevel=kWarning;
  }
	split(flashggCats_,flashggCatsStr_,boost::is_any_of(","));

	int startingCategory=0;
  if (singleCategory >-1){
	ncats=singleCategory+1;
	startingCategory=singleCategory;
  }
	if (isFlashgg_==1){

	ncats= flashggCats_.size();

	}

  if(verbose) std::cout << "[INFO] SaveMultiPdf? " << saveMultiPdf << std::endl;
  TFile *outputfile;
  RooWorkspace *outputws;

  if (saveMultiPdf){
	outputfile = new TFile(outfilename.c_str(),"RECREATE");
	outputws = new RooWorkspace(); outputws->SetName("multipdf");
  }

  system(Form("mkdir -p %s",outDir.c_str()));
  TFile *inFile = TFile::Open(fileName.c_str());
  if (!inFile || inFile->IsZombie()) {
    std::cerr << "[FATAL] Cannot open input file: " << fileName << std::endl;
    return 1;
  }
  RooWorkspace *inWS;
	if(isFlashgg_){
		if (isData_){
			inWS = (RooWorkspace*)inFile->Get("CMS_hza_workspace");
		} else {
			inWS = (RooWorkspace*)inFile->Get("CMS_hza_workspace");
		}
	} else {
		inWS = (RooWorkspace*)inFile->Get("CMS_hza_workspace");
	}
	if (verbose) std::cout << "[INFO]  inWS open " << inWS << std::endl;
	if (!inWS) {
    std::cerr << "[FATAL] Workspace not found in file: " << fileName << std::endl;
    return 1;
  }
	if (saveMultiPdf){
		transferMacros(inFile,outputfile);

		RooRealVar *intL;
		RooRealVar *sqrts;

		if (isFlashgg_){
			intL  = intLumi_;
			sqrts = (RooRealVar*)inWS->var("SqrtS");
			if (!sqrts){ sqrts = new RooRealVar("SqrtS","SqrtS",13.6); }
		  std::cout << "[INFO] got intL and sqrts " << intL << ", " << sqrts << std::endl;


		} else {
			intL  = intLumi_;
			sqrts = (RooRealVar*)inWS->var("Sqrts");
		}
		outputws->import(*intL);
		std::cout << "[INFO] got intL and sqrts " << intL << ", " << sqrts << std::endl;
	}

  // function switch
	// Per-mA Bernstein (2026-06-14 bias study): include Bernstein only at the masses
	// where some GOF-passing order also passes the bias study; bestBernOrder() caps
	// it at that optimal order. Dropped at mA 1, 6, 28 (bestBernOrder==0). Exp/Pow/Lau
	// are always included (they pass GOF wherever any function does and stay unbiased).
	vector<string> functionClasses;
	if (bestBernOrder(mass_ALP) > 0) functionClasses.push_back("Bernstein");
	functionClasses.push_back("Exponential");
	functionClasses.push_back("PowerLaw");
	functionClasses.push_back("Laurent");

	map<string,string> namingMap;
	namingMap.insert(pair<string,string>("Bernstein","pol"));
	namingMap.insert(pair<string,string>("Exponential","exp"));
	namingMap.insert(pair<string,string>("PowerLaw","pow"));
	namingMap.insert(pair<string,string>("Laurent","lau"));

	FILE *resFile ;
	if  (singleCategory >-1) resFile = fopen(Form("%s/fTestResults_%s.txt",outDir.c_str(),flashggCats_[singleCategory].c_str()),"w");
	else resFile = fopen(Form("%s/fTestResults.txt",outDir.c_str()),"w");
  FILE *logFile = fopen(Form("%s/EnvelopeResults.txt",outDir.c_str()),"w");
	vector<map<string,int> > choices_vec;
	vector<map<string,std::vector<int> > > choices_envelope_vec;
	vector<map<string,RooAbsPdf*> > pdfs_vec;

	PdfModelBuilder pdfsModel;
	RooRealVar *mass = (RooRealVar*)inWS->var("CMS_hza_mass");
	std:: cout << "[INFO] Got mass from ws " << mass << std::endl;
	if (!mass) {
    std::cerr << "[FATAL] CMS_hza_mass not found in workspace. Abort." << std::endl;
    return 1;
  }
  mass->setRange(mgg_low, mgg_high);
  mass->setBins((int)nBinsForMass);
	pdfsModel.setObsVar(mass);
	double upperEnvThreshold = 0.10; // Looser F-test cut for envelope construction; truth-model
	                                  // selection above stays at the stricter 0.05 (hardcoded).
	const double minEnvelopeGof = 0.01;

	fprintf(resFile,"Truth Model & d.o.f & $\\Delta NLL_{N+1}$ & $p(\\chi^{2}>\\chi^{2}_{(N\\rightarrow N+1)})$ \\\\\n");
	fprintf(resFile,"\\hline\n");

  std::string ext = "13p6TeV";
	if (isFlashgg_) ext = "13p6TeV";
	for (int cat=startingCategory; cat<ncats; cat++){
		map<string,int> choices;
		map<string,std::vector<int> > choices_envelope;
		map<string,RooAbsPdf*> pdfs;
		string catname;
		if (isFlashgg_){
			catname = Form("%s",flashggCats_[cat].c_str());
		} else {
			catname = Form("cat%d",cat);
		}
		RooAbsData *dataFull = nullptr;
    std::vector<std::string> triedNames;

    if (!datasetName_.empty()) {
      triedNames.push_back(datasetName_);
      dataFull = inWS->data(datasetName_.c_str());
    }
    if (!dataFull) {
      if (isData_) {
        std::vector<std::string> candidates = {
          "Data_13p6TeV",
          "Data_13TeV",
          Form("Data_%s",catname.c_str()),
          Form("Data_%s_%s",catname.c_str(),ext.c_str()),
          "Data"
        };
        for (auto const& n : candidates) {
          triedNames.push_back(n);
          if ((dataFull = inWS->data(n.c_str()))) { if (verbose) std::cout << "[INFO] Using dataset: " << n << std::endl; break; }
        }
      } else {
        std::vector<std::string> candidates = {
          Form("data_mass_%s",catname.c_str()),
          "data_mass"
        };
        for (auto const& n : candidates) {
          triedNames.push_back(n);
          if ((dataFull = inWS->data(n.c_str()))) { if (verbose) std::cout << "[INFO] Using dataset: " << n << std::endl; break; }
        }
      }
    }

    std::cout << "[INFO] Processing for mA = " << mass_ALP << std::endl;

    if (!dataFull) {
      std::cerr << "[ERROR] Could not retrieve dataset for category " << catname << ". Tried names:" << std::endl;
      for (auto const& n : triedNames) std::cerr << "  - " << n << std::endl;
      std::cerr << "[HINT] Provide --datasetName <name> if your dataset has a custom name." << std::endl;
      continue;
    }
    if (dataFull->numEntries() <= 0) {
      std::cerr << "[WARN] Dataset has zero entries for category " << catname << ". Skip." << std::endl;
      continue;
    }

    syncMassRangeToData(mass, dataFull, verbose);

		mass->setBins(nBinsForMass);
		RooAbsData *data = dataFull;
		// Guard against empty datasets: extended fits with n=0 lead to NaN NLL.
		if (!data || data->sumEntries() <= 0.) {
		  std::cerr << "[WARN] Category " << catname << " has zero entries in the fit range. Skip this category." << std::endl;
		  continue;
		}

				RooArgList storedPdfs("store");
				EnvelopeCandidate globalFallback;
				std::vector<EnvelopeCandidate> acceptedEnvelopeCandidates;
				std::map<std::string, EnvelopeCandidate> bestAcceptedByFamily;
				std::map<std::string, EnvelopeCandidate> bestSeenByFamily;

			fprintf(resFile,"\\multicolumn{4}{|c|}{\\textbf{Category %d}} \\\\\n",cat);
			fprintf(resFile,"\\hline\n");

		double MinimimNLLSoFar=1e10;
		int simplebestFitPdfIndex = 0;

    if (PLOT_ONLY) {
      map<string,RooAbsPdf*> pdfs;
      for (auto funcType : functionClasses) {
        int order = 1;
        RooAbsPdf *qpdf = getPdf(pdfsModel, funcType, order, Form("quick_pdf_cat%d_%s",cat,ext.c_str()), mass_ALP);
        if (qpdf) {
          ensureOrderSuffix(qpdf, order);
          pdfs.insert(std::make_pair(Form("%s%d",funcType.c_str(),order), qpdf));
          if (saveMultiPdf) storedPdfs.add(*qpdf);
        }
        choices.insert(std::make_pair(funcType, order-1));
      }

      truth_plot(mass,pdfs,data,Form("%s/truths_cat%d",outDir.c_str(),cat),flashggCats_,cat);

      if (saveMultiPdf && storedPdfs.getSize()>0) {
        string catindexname, catname;
        if (isFlashgg_){
          catindexname = Form("pdfindex_%s_%s",flashggCats_[cat].c_str(),ext.c_str());
          catname = Form("%s",flashggCats_[cat].c_str());
        } else {
          catindexname = Form("pdfindex_cat%d_%s",cat,ext.c_str());
          catname = Form("cat%d",cat);
        }
        RooCategory catIndex(catindexname.c_str(),"c");
        RooMultiPdf *pdf = new RooMultiPdf(Form("CMS_hgg_%s_%s_bkgshape",catname.c_str(),ext.c_str()),"all pdfs",catIndex,storedPdfs);
        const double nData = data->sumEntries();
        const double nInit = (nData > 0.) ? nData : 1.0;
        const double nMax  = (nData > 0.) ? 3.0*nData : 10.0;
        RooRealVar nBackground(Form("CMS_hgg_%s_%s_bkgshape_norm",catname.c_str(),ext.c_str()),"nbkg", nInit, 1e-6, nMax);
        mass->setBins(nBinsForMass);
        RooDataHist dataBinned(Form("roohist_data_mass_%s",catname.c_str()),"data",*mass,*dataFull);

        outputws->import(*pdf);
        outputws->import(nBackground);
        outputws->import(catIndex);
        outputws->import(dataBinned);
        outputws->import(*data);

        int bestFitPdfIndex = 0;
        catIndex.setIndex(bestFitPdfIndex);
        multipdf_plot(mass,pdf,&catIndex,data,Form("%s/multipdf_%s",outDir.c_str(),catname.c_str()),flashggCats_,cat,bestFitPdfIndex);
      }

      choices_vec.push_back(choices);
      choices_envelope_vec.push_back(map<string,std::vector<int> >());
      pdfs_vec.push_back(map<string,RooAbsPdf*>());
      continue;
    }

		for (vector<string>::iterator funcType=functionClasses.begin();
				funcType!=functionClasses.end(); funcType++){

			double thisNll=0.; double prevNll=0.; double chi2=0.; double prob=0.;
			int order=1; int prev_order=0; int cache_order=0;

			RooAbsPdf *prev_pdf=NULL;
			RooAbsPdf *cache_pdf=NULL;
				int counter =0;
				while (prob<0.05 && order < envMaxOrder(*funcType, mass_ALP) ){
				
        RooAbsPdf *bkgPdf = getPdf(pdfsModel,*funcType,order,"", mass_ALP);
				if (!bkgPdf){
					order++;
				}

				else {
					int fitStatus = 0;
					bkgPdf->Print();
					runFit(bkgPdf,data,&thisNll,&fitStatus,/*max iterations*/3);

          if (fitStatus!=0) std::cout << "[WARNING] Warning -- Fit status for " << bkgPdf->GetName() << " at " << fitStatus <<std::endl;
          
					chi2 = 2.*(prevNll-thisNll);
					if (chi2<0. && order>1) chi2=0.;
					if (prev_pdf!=NULL){
						prob = getProbabilityFtest(chi2,order-prev_order,prev_pdf,bkgPdf,mass,data
								,Form("%s/Ftest_from_%s%d_cat%d.pdf",outDir.c_str(),funcType->c_str(),order,cat));
						std::cout << "[INFO]  F-test Prob(chi2>chi2(data)) == " << prob << std::endl;
					} else {
						prob = 0;
					}
					double gofProb=0;
					if (!saveMultiPdf && fitStatus != 5) eachFunc_plot(mass,bkgPdf,data,Form("%s/%s%d_cat%d.pdf",outDir.c_str(),funcType->c_str(),order,cat),flashggCats_,fitStatus,&gofProb);
          cout << "[INFO]\t funcType: " << *funcType << " order: " << order << " prevNll: " << prevNll << " thisNll: " << thisNll << " chi2: " << chi2 << " prob: " << prob << endl;
					prevNll=thisNll;
					cache_order=prev_order;
					cache_pdf=prev_pdf;
					prev_order=order;
					prev_pdf=bkgPdf;
					order++;
				}
				counter++;
			}

			fprintf(resFile,"%15s & %d & %5.2f & %5.2f \\\\\n",funcType->c_str(),cache_order+1,chi2,prob);
			choices.insert(pair<string,int>(*funcType,cache_order));

			if (cache_pdf) {
        pdfs.insert(pair<string,RooAbsPdf*>(Form("%s%d",funcType->c_str(),cache_order),cache_pdf));
      } else {
        std::cerr << "[WARN] No valid cached pdf for family " << *funcType << " in " << catname << ". Skip adding to truth set." << std::endl;
      }

			int truthOrder = cache_order;

			if (saveMultiPdf){
				chi2=0.;
				thisNll=0.;
				prevNll=0.;
				prob=0.;
				order=1;
				prev_order=0.;
				cache_order=0.;
				std::cout << "[INFO] Determining Envelope Functions for Family " << *funcType << ", cat " << cat << std::endl;
				std::cout << "[INFO] Upper end Threshold for highest order function " << upperEnvThreshold <<std::endl;


				while (prob<upperEnvThreshold && order < envMaxOrder(*funcType, mass_ALP) ){
					RooAbsPdf *bkgPdf = getPdf(pdfsModel,*funcType,order,"", mass_ALP);

          if (!bkgPdf ){
						if (order >6) { std::cout << " [WARNING] could not add ] " << std::endl; break ;}
						order++;
					}
					else {
            ensureOrderSuffix(bkgPdf, order);

						int fitStatus=0;
						runFit(bkgPdf,data,&thisNll,&fitStatus,/*max iterations*/3);
						if (fitStatus!=0) std::cout << "[WARNING] Warning -- Fit status for " << bkgPdf->GetName() << " at " << fitStatus <<std::endl;
						double myNll = 2.*thisNll;
						chi2 = 2.*(prevNll-thisNll);
						if (chi2<0. && order>1) chi2=0.;
						prob = TMath::Prob(chi2,order-prev_order);

						cout << "[INFO] \t funcType: " << *funcType << " order: " << order << " prevNll: " << prevNll << " thisNll: " << thisNll << " chi2: " << chi2 << " prob: " << prob << endl;
						prevNll=thisNll;
						cache_order=prev_order;
						cache_pdf=prev_pdf;

								double gofProb =0;

	            if(fitStatus != 5)
	            {
							  eachFunc_plot(mass,bkgPdf,data,Form("%s/%s%d_cat%d.pdf",outDir.c_str(),funcType->c_str(),order,cat),flashggCats_,fitStatus,&gofProb);
	            }

								const int nvars = bkgPdf->getVariables() ? bkgPdf->getVariables()->getSize() : 0;
								const double candidateScore = myNll + nvars;
								const bool candidateHasValidGof = (gofProb >= 0.0);
									if (fitStatus != 5 && std::isfinite(candidateScore)) {
										EnvelopeCandidate fallbackCandidate;
									fallbackCandidate.pdf = bkgPdf;
									fallbackCandidate.family = *funcType;
									fallbackCandidate.name = bkgPdf->GetName();
									fallbackCandidate.order = order;
									fallbackCandidate.gof = gofProb;
									fallbackCandidate.score = candidateScore;
									fallbackCandidate.isTruth = (order == truthOrder);
									const bool takeFallback =
										(!globalFallback.pdf) ||
										betterEnvelopeCandidate(fallbackCandidate, globalFallback);
										if (takeFallback) {
											globalFallback = fallbackCandidate;
										}
										auto seenIt = bestSeenByFamily.find(*funcType);
										if (seenIt == bestSeenByFamily.end() || betterEnvelopeCandidate(fallbackCandidate, seenIt->second)) {
											bestSeenByFamily[*funcType] = fallbackCandidate;
										}
									}
		            
									if ((prob < upperEnvThreshold) ) {

									if (gofProb > minEnvelopeGof) {
										std::cout << "[INFO] Adding to Envelope " << bkgPdf->GetName() << " "<< gofProb
											<< " 2xNLL + c is " << myNll + bkgPdf->getVariables()->getSize() <<  std::endl;
										EnvelopeCandidate acceptedCandidate;
										acceptedCandidate.pdf = bkgPdf;
										acceptedCandidate.family = *funcType;
										acceptedCandidate.name = bkgPdf->GetName();
										acceptedCandidate.order = order;
										acceptedCandidate.gof = gofProb;
											acceptedCandidate.score = myNll + bkgPdf->getVariables()->getSize();
											acceptedCandidate.isTruth = (order == truthOrder);
											acceptedEnvelopeCandidates.push_back(acceptedCandidate);
											auto bestAccIt = bestAcceptedByFamily.find(*funcType);
											if (bestAccIt == bestAcceptedByFamily.end() || lowerScoreEnvelopeCandidate(acceptedCandidate, bestAccIt->second)) {
												bestAcceptedByFamily[*funcType] = acceptedCandidate;
											}
										} else {
										std::cout << "[INFO] Rejecting from Envelope " << bkgPdf->GetName()
										          << " because gof=" << gofProb
									          << " <= " << minEnvelopeGof << std::endl;
								}
							}

						prev_order=order;
						prev_pdf=bkgPdf;
						order++;
					}
				}

					fprintf(resFile,"%15s & %d & %5.2f & %5.2f \\\\\n",funcType->c_str(),cache_order+1,chi2,prob);
				}
			}

			fprintf(resFile,"\\hline\n");
			choices_vec.push_back(choices);
			choices_envelope_vec.push_back(choices_envelope);
			pdfs_vec.push_back(pdfs);

		truth_plot(mass,pdfs,data,Form("%s/truths_cat%d",outDir.c_str(),cat),flashggCats_,cat);

		if (saveMultiPdf){

		      std::vector<EnvelopeCandidate> finalEnvelopeCandidates;
		      std::map<std::string, bool> seenEnvelopeKeys;
		      auto addUniqueEnvelopeCandidate = [&](const EnvelopeCandidate& cand) {
		        if (!cand.pdf) return;
		        const std::string key = envelopeCandidateKey(cand);
		        if (seenEnvelopeKeys[key]) return;
		        seenEnvelopeKeys[key] = true;
		        finalEnvelopeCandidates.push_back(cand);
		      };

		      for (const auto& family : functionClasses) {
		        auto bestAcceptedIt = bestAcceptedByFamily.find(family);
		        if (bestAcceptedIt != bestAcceptedByFamily.end()) {
		          addUniqueEnvelopeCandidate(bestAcceptedIt->second);
		          continue;
		        }
		        auto bestSeenIt = bestSeenByFamily.find(family);
		        if (bestSeenIt != bestSeenByFamily.end()) {
		          EnvelopeCandidate familyRep = bestSeenIt->second;
		          familyRep.fromFallback = true;
		          familyRep.fromFamilyFallback = true;
		          addUniqueEnvelopeCandidate(familyRep);
		        }
		      }

		      std::stable_sort(acceptedEnvelopeCandidates.begin(), acceptedEnvelopeCandidates.end(), lowerScoreEnvelopeCandidate);
		      const int minEnvelopeSize = std::max(MIN_ENVELOPE_PDFS, (int)functionClasses.size());
		      const int targetEnvelopeSize = std::max(minEnvelopeSize, (int)finalEnvelopeCandidates.size());
		      for (const auto& cand : acceptedEnvelopeCandidates) {
		        addUniqueEnvelopeCandidate(cand);
		        if ((int)finalEnvelopeCandidates.size() >= targetEnvelopeSize) break;
		      }

		      if (finalEnvelopeCandidates.empty()) {
		        if (globalFallback.pdf) {
		          std::cerr << "[WARN] storedPdfs is empty for " << catname
		                    << ". Forcing fallback pdf " << globalFallback.name
		                    << " with score=" << globalFallback.score
		                    << " and gof=" << globalFallback.gof << std::endl;
		          globalFallback.fromFallback = true;
		          addUniqueEnvelopeCandidate(globalFallback);
		        } else {
		          std::cerr << "[WARN] storedPdfs is empty for " << catname << ". Skip MultiPdf for this category." << std::endl;
		          continue;
		        }
		      }

		      if ((int)acceptedEnvelopeCandidates.size() > targetEnvelopeSize) {
		        std::cout << "[INFO] Pruning envelope in " << catname << " from "
		                  << acceptedEnvelopeCandidates.size() << " to protected size "
		                  << finalEnvelopeCandidates.size() << " (target " << targetEnvelopeSize
		                  << ", with at least one representative per family)" << std::endl;
		      }

		      std::stable_sort(finalEnvelopeCandidates.begin(), finalEnvelopeCandidates.end(), lowerScoreEnvelopeCandidate);
		      for (const auto& cand : finalEnvelopeCandidates) {
		        if (!cand.pdf) continue;
		        storedPdfs.add(*cand.pdf);
		        choices_envelope[cand.family].push_back(cand.order);
		        if (cand.score < MinimimNLLSoFar) {
		          simplebestFitPdfIndex = storedPdfs.getSize()-1;
		          MinimimNLLSoFar = cand.score;
		        }
		        if (logFile) {
		          fprintf(logFile,
		                  "category : %d , pdf : %s , gof : %f, isTruth : %d%s%s\n ",
		                  cat,
		                  cand.name.c_str(),
		                  cand.gof,
		                  cand.isTruth ? 1 : 0,
		                  cand.fromFallback ? " [fallback-empty-envelope]" : "",
		                  cand.fromFamilyFallback ? " [family-fallback]" : "");
		        }
		      }

		      choices_envelope_vec.back() = choices_envelope;

				string catindexname;
				string catname;
			if (isFlashgg_){
        catindexname = Form("pdfindex_%s_%s",flashggCats_[cat].c_str(),ext.c_str());
				catname = Form("%s",flashggCats_[cat].c_str());
			} else {
        catindexname = Form("pdfindex_cat%d_%s",cat,ext.c_str());
				catname = Form("cat%d",cat);
			}

			RooCategory catIndex(catindexname.c_str(),"c");
			RooMultiPdf *pdf = new RooMultiPdf(Form("CMS_hgg_%s_%s_bkgshape",catname.c_str(),ext.c_str()),"all pdfs",catIndex,storedPdfs);
			const double nData = data->sumEntries();
			const double nInit = (nData > 0.) ? nData : 1.0;
			const double nMax  = (nData > 0.) ? 3.0*nData : 10.0;
			RooRealVar nBackground(Form("CMS_hgg_%s_%s_bkgshape_norm",catname.c_str(),ext.c_str()),"nbkg", nInit, 1e-6, nMax);
			int bestFitPdfIndex = getBestFitFunction(pdf,data,&catIndex,!verbose);
			catIndex.setIndex(bestFitPdfIndex);
			std::cout << "// ------------------------------------------------------------------------- //" <<std::endl;
			std::cout << "[INFO] Created MultiPdf " << pdf->GetName() << ", in Category " << cat << " with a total of " << catIndex.numTypes() << " pdfs"<< std::endl;
			storedPdfs.Print();
			std::cout << "[INFO] Best Fit Pdf = " << bestFitPdfIndex << ", " << storedPdfs.at(bestFitPdfIndex)->GetName() << std::endl;
			std::cout << "// ------------------------------------------------------------------------- //" <<std::endl;
			std::cout << "[INFO] Simple check of index "<< simplebestFitPdfIndex <<std::endl;

			mass->setBins(nBinsForMass);
      RooDataHist dataBinned(Form("roohist_data_mass_%s",catname.c_str()),"data",*mass,*dataFull);

			outputws->import(*pdf);
			outputws->import(nBackground);
			outputws->import(catIndex);
			outputws->import(dataBinned);
			outputws->import(*data);
			multipdf_plot(mass,pdf,&catIndex,data,Form("%s/multipdf_%s",outDir.c_str(),catname.c_str()),flashggCats_,cat,bestFitPdfIndex);
      
		}

		}
		if (saveMultiPdf){
			outputfile->cd();
			outputws->Write();
			outputfile->Close();
		}

		FILE *dfile = fopen(datfile.c_str(),"w");
		cout << "[RESULT] Recommended options" << endl;

		for (int cat=startingCategory; cat<ncats; cat++){
			cout << "Cat " << cat << endl;
			fprintf(dfile,"cat=%d\n",cat);
			for (map<string,int>::iterator it=choices_vec[cat-startingCategory].begin(); it!=choices_vec[cat-startingCategory].end(); it++){
				cout << "\t" << it->first << " - " << it->second << endl;
				fprintf(dfile,"truth=%s:%d:%s%d\n",it->first.c_str(),it->second,namingMap[it->first].c_str(),it->second);
			}
			for (map<string,std::vector<int> >::iterator it=choices_envelope_vec[cat-startingCategory].begin(); it!=choices_envelope_vec[cat-startingCategory].end(); it++){
				std::vector<int> ords = it->second;
				for (std::vector<int>::iterator ordit=ords.begin(); ordit!=ords.end(); ordit++){
					fprintf(dfile,"paul=%s:%d:%s%d\n",it->first.c_str(),*ordit,namingMap[it->first].c_str(),*ordit);
				}
			}
			fprintf(dfile,"\n");
		}
		inFile->Close();

    if (logFile) fclose(logFile);

		return 0;
}
static void transferMacros(TFile* in, TFile* out) {
  if (!in || !out) return;
  out->cd();
  TIter next(in->GetListOfKeys());
  TKey* key = nullptr;
  while ((key = (TKey*)next())) {
    const char* className = key->GetClassName();
    if (!className) continue;
    if (std::string(className) == "TMacro") {
      TObject* obj = key->ReadObj();
      if (obj) {
        obj->Write();
        delete obj;
      }
    }
  }
}
static void ensureOrderSuffix(RooAbsPdf* pdf, int order) {
  if (!pdf) return;
  std::string name(pdf->GetName() ? pdf->GetName() : "");
  if (!name.empty() && std::isdigit(name.back())) return;
  pdf->SetName(Form("%s%d", name.c_str(), order));
}

static void sanitizePdfParams(RooAbsPdf* pdf, RooAbsData* data, bool verboseDiag) {
  if (!pdf || !data) return;
  std::unique_ptr<RooArgSet> params(pdf->getParameters(*data));
  if (!params) return;

  std::unique_ptr<TIterator> it(params->createIterator());
  for (RooAbsArg* a = (RooAbsArg*)it->Next(); a; a = (RooAbsArg*)it->Next()) {
    RooRealVar* v = dynamic_cast<RooRealVar*>(a);
    if (!v) continue;
    if (v->isConstant()) continue;

    // Reset NaN/Inf values to something sane.
    const double val = v->getVal();
    if (!std::isfinite(val)) {
      double newv = 0.0;
      if (v->hasMin() && v->hasMax() && std::isfinite(v->getMin()) && std::isfinite(v->getMax()) && v->getMax() > v->getMin()) {
        newv = 0.5*(v->getMin() + v->getMax());
      } else if (v->hasMin() && std::isfinite(v->getMin())) {
        newv = v->getMin() + 1.0;
      } else if (v->hasMax() && std::isfinite(v->getMax())) {
        newv = v->getMax() - 1.0;
      }
      v->setVal(newv);
      if (verboseDiag) {
        std::cerr << "[WARN] sanitizePdfParams: reset non-finite " << v->GetName()
                  << " to " << newv << std::endl;
      }
    }

	    // Nudge away from hard boundaries (helps Minuit).
	    if (v->hasMin() && std::isfinite(v->getMin()) && v->getVal() <= v->getMin()) {
	      const double span = (v->hasMax() && std::isfinite(v->getMax())) ? (v->getMax() - v->getMin()) : 1.0;
	      v->setVal(v->getMin() + 1e-3*span);
    }
	    if (v->hasMax() && std::isfinite(v->getMax()) && v->getVal() >= v->getMax()) {
	      const double span = (v->hasMin() && std::isfinite(v->getMin())) ? (v->getMax() - v->getMin()) : 1.0;
	      v->setVal(v->getMax() - 1e-3*span);
	    }

	    if (v->hasMin() && v->hasMax() && std::isfinite(v->getMin()) && std::isfinite(v->getMax()) && v->getMax() > v->getMin()) {
	      const double span = v->getMax() - v->getMin();
	      std::string vname = v->GetName() ? v->GetName() : "";
	      double guardFrac = 1e-3;
	      if (vname.find("stepWidth") != std::string::npos || vname.find("_width_") != std::string::npos) {
	        guardFrac = 0.15;
	      } else if (vname.find("gsigma") != std::string::npos || vname.find("_sigma_") != std::string::npos) {
	        guardFrac = 0.12;
	      } else if (vname.find("turnon") != std::string::npos || vname.find("_step") != std::string::npos) {
	        guardFrac = 0.10;
	      }
	      const double guardedLow = v->getMin() + guardFrac*span;
	      const double guardedHigh = v->getMax() - guardFrac*span;
	      if (guardedHigh > guardedLow) {
	        if (v->getVal() < guardedLow) v->setVal(guardedLow);
	        if (v->getVal() > guardedHigh) v->setVal(guardedHigh);
	      }
	    }

    // Ensure step size is finite and non-zero.
    const double err = v->getError();
    if (!std::isfinite(err) || err <= 0.0) {
      double step = 1e-2;
      if (v->hasMin() && v->hasMax() && std::isfinite(v->getMin()) && std::isfinite(v->getMax()) && v->getMax() > v->getMin()) {
        step = 0.1*std::fabs(v->getMax() - v->getMin());
      } else {
        const double aval = std::fabs(v->getVal());
        step = std::max(1e-2, 0.1*aval);
      }
      if (!std::isfinite(step) || step <= 0.0) step = 1e-2;
      v->setError(step);
    }
  }
}

static void randomizePdfParamsUniform(RooAbsPdf* pdf, RooAbsData* data) {
  if (!pdf || !data || !RandomGen) return;

  std::unique_ptr<RooArgSet> params(pdf->getParameters(*data));
  if (!params) return;

  std::unique_ptr<TIterator> it(params->createIterator());
  for (RooAbsArg* a = (RooAbsArg*)it->Next(); a; a = (RooAbsArg*)it->Next()) {
    RooRealVar* v = dynamic_cast<RooRealVar*>(a);
    if (!v) continue;
    if (v->isConstant()) continue;

    double newv = v->getVal();

    if (v->hasMin() && v->hasMax() && std::isfinite(v->getMin()) && std::isfinite(v->getMax()) && v->getMax() > v->getMin()) {
      newv = RandomGen->Uniform(v->getMin(), v->getMax());
    } else {
      double step = v->getError();
      if (!std::isfinite(step) || step <= 0.0) {
        const double aval = std::fabs(v->getVal());
        step = std::max(1e-2, 0.1*aval);
      }
      newv = v->getVal() + RandomGen->Gaus(0.0, step);
    }

    if (!std::isfinite(newv)) {
      if (v->hasMin() && v->hasMax() && std::isfinite(v->getMin()) && std::isfinite(v->getMax()) && v->getMax() > v->getMin()) {
        newv = 0.5*(v->getMin() + v->getMax());
      } else {
        newv = 0.0;
      }
    }
    v->setVal(newv);
  }
}

static int getBestFitFunction(RooMultiPdf *mpdf, RooAbsData *data, RooCategory *catIndex, bool silent) {
  if (!mpdf || !data || !catIndex) return 0;
  if (PLOT_ONLY) return 0;
  int n = catIndex->numTypes();
  if (n<=0) return 0;
  int bestIndex = 0;
  double bestScore = 1e300;
  int current = catIndex->getIndex();
  for (int i=0; i<n; ++i) {
    catIndex->setIndex(i);
    RooAbsPdf* pdf = mpdf->getCurrentPdf();
    if (!pdf) continue;
    double nll = 0.;
    int status = 0;
    runFit(pdf, data, &nll, &status, 3);
    if (status!=0 && !silent) {
      std::cout << "[WARN] getBestFitFunction: fit status " << status << " for index " << i << std::endl;
    }
    RooArgSet* pars = pdf->getParameters(*data);
    int npar = pars ? pars->getSize() : 0;
    double score = 2.0*nll + npar;
    if (!std::isfinite(score)) score = 1e300;
    if (score < bestScore) {
      bestScore = score;
      bestIndex = i;
    }
  }
  catIndex->setIndex(current);
  return bestIndex;
}
static bool checkPdfDataObservables(RooAbsPdf* pdf, RooAbsData* data, bool verboseDiag) {
  if (!pdf || !data) return false;

  std::unique_ptr<RooArgSet> pdfObs(pdf->getObservables(*data));
  const RooArgSet* dataVars = data->get(); // data 的當前 row vars（包含 observables）

  if (!pdfObs || !dataVars) {
    if (verboseDiag) {
      std::cerr << "[ERROR] checkPdfDataObservables: null pdfObs or dataVars"
                << " pdfObs=" << pdfObs.get() << " dataVars=" << dataVars << std::endl;
    }
    return false;
  }

  std::unique_ptr<TIterator> it(pdfObs->createIterator());
  for (RooAbsArg* a = (RooAbsArg*)it->Next(); a; a = (RooAbsArg*)it->Next()) {
    if (!a) continue;
    const char* n = a->GetName();
    if (!n || !(*n)) continue;

    if (!dataVars->find(n)) {
      if (verboseDiag) {
        std::cerr << "[ERROR] Observable '" << n << "' required by pdf '" << pdf->GetName()
                  << "' is missing in data '" << data->GetName() << "'" << std::endl;
        std::cerr << "  pdf observables: "; pdfObs->Print("V");
        std::cerr << "  data vars:       "; dataVars->Print("V");
      }
      return false;
    }
  }
  return true;
}
static void syncMassRangeToData(RooRealVar* mass, RooAbsData* data, bool verboseDiag) {
  if (!mass) return;

  double dmin = 0., dmax = 0.;
  if (data) data->getRange(*mass, dmin, dmax);

  mgg_low  = mgglow_;
  mgg_high = mgghigh_;
  nBinsForMass = 1.0f * (mgg_high - mgg_low);

  mgg_blind_low  = std::min(std::max(mgg_blind_low,  mgg_low),  mgg_high);
  mgg_blind_high = std::min(std::max(mgg_blind_high, mgg_low),  mgg_high);
  if (mgg_blind_high < mgg_blind_low) std::swap(mgg_blind_low, mgg_blind_high);

  if (verboseDiag) {
    std::cout << "[INFO] syncMassRangeToData: data range [" << dmin << "," << dmax << "], "
              << "forcing fixed window [" << mgg_low << "," << mgg_high << "]\n";
  }

  mass->setRange(mgg_low, mgg_high);
  mass->setBins((int)nBinsForMass);
}
