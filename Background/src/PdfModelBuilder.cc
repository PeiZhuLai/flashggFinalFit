#include "TCanvas.h"

#include "RooPlot.h"
#include "RooBernstein.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooGenericPdf.h"
#include "RooExponential.h"
#include "../interface/RooPowerLaw.h"
#include "../interface/RooPowerLawSum.h"
#include "RooKeysPdf.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooConstVar.h"
#include "RooFitResult.h"
#include "RooRandom.h"
#include "RooGaussian.h"
#include "RooFFTConvPdf.h"
#include "RooProdPdf.h"
#include "RooNumConvPdf.h"
#include "RooGaussModel.h"
#include "TMath.h"
#include "RooAbsData.h"
#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string/classification.hpp"
#include "boost/algorithm/string/predicate.hpp"
#include "../interface/PdfModelBuilder.h"
#include "HiggsAnalysis/CombinedLimit/interface/HGGRooPdfs.h"
#include "HiggsAnalysis/CombinedLimit/interface/HZGRooPdfs.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooBernsteinFast.h"

using namespace std;
using namespace RooFit;
using namespace boost;

namespace {
constexpr int kBkgFftBins = 1024;
constexpr int kBkgCacheBins = 1024;
constexpr double kBkgBufferFraction = 0.25;

struct StepGausWindow {
  double turnon;
  double turnonLo;
  double turnonHi;
  double sigma;
  double sigmaLo;
  double sigmaHi;
  double width;
  double widthLo;
  double widthHi;
};

inline double clampDouble(double v, double lo, double hi) {
  return std::max(lo, std::min(v, hi));
}

// RooFit normalizes these shapes, so one overall coefficient scale is unidentifiable.
inline void fixUnitShapeCoeff(RooRealVar* coeff) {
  if (!coeff) return;
  coeff->setRange(0.0, 10.0);
  coeff->setVal(1.0);
  coeff->setConstant(true);
}

inline StepGausWindow makeStableStepGausWindow(int massALP, double turnonShift, double sigmaScale, double widthScale) {
  const double m = clampDouble((double)massALP, 1.0, 30.0);
  const double baseTurnon = clampDouble(104.5 + 0.32*m + turnonShift, 103.0, 117.0);
  const double baseSigma = clampDouble((2.2 + 0.14*m)*sigmaScale, 1.2, 8.5);
  const double baseWidth = clampDouble((0.9 + 0.08*m)*widthScale, 0.8, 5.0);

  StepGausWindow cfg;
  cfg.turnon = baseTurnon;
  cfg.turnonLo = std::max(100.0, baseTurnon - 3.5);
  cfg.turnonHi = std::min(125.0, baseTurnon + 3.5);
  if (cfg.turnonHi <= cfg.turnonLo) cfg.turnonHi = cfg.turnonLo + 1.0;

  cfg.sigma = baseSigma;
  cfg.sigmaLo = std::max(0.8, 0.45*baseSigma);
  cfg.sigmaHi = std::min(15.0, std::max(baseSigma + 2.0, 1.9*baseSigma));
  if (cfg.sigmaHi <= cfg.sigmaLo) cfg.sigmaHi = cfg.sigmaLo + 1.0;

  cfg.width = baseWidth;
  cfg.widthLo = std::max(0.5, 0.45*baseWidth);
  cfg.widthHi = std::min(20.0, std::max(baseWidth + 1.5, 2.6*baseWidth));
  if (cfg.widthHi <= cfg.widthLo) cfg.widthHi = cfg.widthLo + 1.0;

  return cfg;
}

inline void setFftObsBinning(RooRealVar* obs) {
  if (!obs) return;
  obs->setBins(kBkgFftBins, "fft");
  obs->setBins(kBkgCacheBins, "cache");
  obs->setBins(kBkgCacheBins);
}

inline void setCacheObsBinning(RooRealVar* obs) {
  if (!obs) return;
  obs->setBins(kBkgCacheBins, "cache");
  obs->setBins(kBkgCacheBins);
}
}

PdfModelBuilder::PdfModelBuilder():
  obs_var_set(false),
  signal_modifier_set(false),
  signal_set(false),
  bkgHasFit(false),
  sbHasFit(false),
  keysPdfAttributesSet(false),
  verbosity(0)
{
  
  recognisedPdfTypes.push_back("Bernstein");
  recognisedPdfTypes.push_back("Exponential");
  recognisedPdfTypes.push_back("PowerLaw");
  recognisedPdfTypes.push_back("Laurent");
  recognisedPdfTypes.push_back("KeysPdf");
  recognisedPdfTypes.push_back("File");
  //bing
  recognisedPdfTypes.push_back("BernsteinStepxGau");
  recognisedPdfTypes.push_back("PowerLawStepxGau");
  recognisedPdfTypes.push_back("LaurentStepxGau");
  recognisedPdfTypes.push_back("ExponentialStepxGau");

  wsCache = new RooWorkspace("PdfModelBuilderCache");
};

PdfModelBuilder::~PdfModelBuilder(){};

void PdfModelBuilder::setObsVar(RooRealVar *var){
  obs_var=var;
  obs_var_set=true;
}

void PdfModelBuilder::setSignalModifier(RooRealVar *var){
  signalModifier=var;
  signal_modifier_set=true;
}

void PdfModelBuilder::setSignalModifierVal(float val){
  signalModifier->setVal(val);
}

void PdfModelBuilder::setSignalModifierConstant(bool val){
  signalModifier->setConstant(val);
}

RooAbsPdf* PdfModelBuilder::getChebychev(string prefix, int order){
  
  RooArgList *coeffList = new RooArgList();
  for (int i=0; i<order; i++){
    string name = Form("%s_p%d",prefix.c_str(),i);
    RooRealVar *param = new RooRealVar(name.c_str(),name.c_str(),0.01,-10.,10.);
    params.insert(pair<string,RooRealVar*>(name,param));
    coeffList->add(*params[name]);
  }
  RooPolynomial *cheb = new RooPolynomial(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  return cheb;
}

RooAbsPdf* PdfModelBuilder::getBernstein(string prefix, int order){
  
  RooArgList *coeffList = new RooArgList();
  for (int i=0; i<order; i++){
    string name = Form("%s_p%d",prefix.c_str(),i);
    RooRealVar *param = new RooRealVar(name.c_str(),name.c_str(),0.1*(i+1),-5.,5.);
    RooFormulaVar *form = new RooFormulaVar(Form("%s_sq",name.c_str()),Form("%s_sq",name.c_str()),"@0*@0",RooArgList(*param));
    params.insert(pair<string,RooRealVar*>(name,param));
    prods.insert(pair<string,RooFormulaVar*>(name,form));
    coeffList->add(*prods[name]);
  }
  RooRealVar *fgaus = new RooRealVar("fgaus", "gaus fraction",0.5,0.,1.) ;
  RooRealVar *mean = new RooRealVar(Form("%s_mean",prefix.c_str()),Form("%s_mean",prefix.c_str()),0.0,-5.0,5.0) ;
  mean->setConstant(true);
  RooRealVar *sigma = new RooRealVar(Form("%s_sigma",prefix.c_str()),Form("%s_sigma",prefix.c_str()),8,0.2,30.0) ;
  RooGaussian *gaus = new RooGaussian(Form("%s_gaus",prefix.c_str()),Form("%s_gaus",prefix.c_str()),*obs_var,*mean,*sigma) ;
  RooRealVar *step_value = new RooRealVar(Form("%s_step",prefix.c_str()), Form("%s_step",prefix.c_str()), 105., 100., 110.);
  RooRealVar *step_width = new RooRealVar(Form("%s_stepWidth",prefix.c_str()), Form("%s_stepWidth",prefix.c_str()), 2.0, 0.2, 10.0);
  if (order==1) {
	 RooBernsteinFast<1> *bern = new RooBernsteinFast<1>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
    RooGenericPdf *soft_step = new RooGenericPdf(
      Form("%s_softstep",prefix.c_str()), Form("%s_softstep",prefix.c_str()),
      "1e-20 + 0.5*(1.0 + TMath::Erf((@0-@1)/(@2*sqrt(2.)))) * @3",
      RooArgList(*obs_var, *step_value, *step_width, *bern)
    );
    obs_var->setRange(-200, 500);
    setFftObsBinning(obs_var);
    RooFFTConvPdf *bern_gaus = new RooFFTConvPdf(Form("%s_berngaus",prefix.c_str()), Form("%s_berngaus",prefix.c_str()), *obs_var, *gaus, *soft_step);
    setCacheObsBinning(obs_var);
    bern_gaus->setBufferFraction(kBkgBufferFraction);
    obs_var->setRange(95, 180);
    return bern_gaus;
  } else if (order==2) {
	RooBernsteinFast<2> *bern = new RooBernsteinFast<2>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
    RooGenericPdf *soft_step = new RooGenericPdf(
      Form("%s_softstep",prefix.c_str()), Form("%s_softstep",prefix.c_str()),
      "1e-20 + 0.5*(1.0 + TMath::Erf((@0-@1)/(@2*sqrt(2.)))) * @3",
      RooArgList(*obs_var, *step_value, *step_width, *bern)
    );
    obs_var->setRange(-200, 500);
    setFftObsBinning(obs_var);
    RooFFTConvPdf *bern_gaus = new RooFFTConvPdf(Form("%s_berngaus",prefix.c_str()), Form("%s_berngaus",prefix.c_str()), *obs_var, *gaus, *soft_step);
    setCacheObsBinning(obs_var);
    bern_gaus->setBufferFraction(kBkgBufferFraction);
    obs_var->setRange(95, 180);
    return bern_gaus;
  } else if (order==3) {
	RooBernsteinFast<3> *bern = new RooBernsteinFast<3>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
     RooGenericPdf *soft_step = new RooGenericPdf(
      Form("%s_softstep",prefix.c_str()), Form("%s_softstep",prefix.c_str()),
      "1e-20 + 0.5*(1.0 + TMath::Erf((@0-@1)/(@2*sqrt(2.)))) * @3",
      RooArgList(*obs_var, *step_value, *step_width, *bern)
    );
    obs_var->setRange(-200, 500);
    setFftObsBinning(obs_var);
    RooFFTConvPdf *bern_gaus = new RooFFTConvPdf(Form("%s_berngaus",prefix.c_str()), Form("%s_berngaus",prefix.c_str()), *obs_var, *gaus, *soft_step);
    setCacheObsBinning(obs_var);
    bern_gaus->setBufferFraction(kBkgBufferFraction);
    obs_var->setRange(95, 180);
    return bern_gaus;
  } else if (order==4) {
	RooBernsteinFast<4> *bern = new RooBernsteinFast<4>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
     RooGenericPdf *soft_step = new RooGenericPdf(
      Form("%s_softstep",prefix.c_str()), Form("%s_softstep",prefix.c_str()),
      "1e-20 + 0.5*(1.0 + TMath::Erf((@0-@1)/(@2*sqrt(2.)))) * @3",
      RooArgList(*obs_var, *step_value, *step_width, *bern)
    );
    obs_var->setRange(-200, 500);
    setFftObsBinning(obs_var);
    RooFFTConvPdf *bern_gaus = new RooFFTConvPdf(Form("%s_berngaus",prefix.c_str()), Form("%s_berngaus",prefix.c_str()), *obs_var, *gaus, *soft_step);
    setCacheObsBinning(obs_var);
    bern_gaus->setBufferFraction(kBkgBufferFraction);
    obs_var->setRange(95, 180);
    return bern_gaus;
  } else if (order==5) {
	RooBernsteinFast<5> *bern = new RooBernsteinFast<5>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
     RooGenericPdf *soft_step = new RooGenericPdf(
      Form("%s_softstep",prefix.c_str()), Form("%s_softstep",prefix.c_str()),
      "1e-20 + 0.5*(1.0 + TMath::Erf((@0-@1)/(@2*sqrt(2.)))) * @3",
      RooArgList(*obs_var, *step_value, *step_width, *bern)
    );
    obs_var->setRange(-200, 500);
    setFftObsBinning(obs_var);
    RooFFTConvPdf *bern_gaus = new RooFFTConvPdf(Form("%s_berngaus",prefix.c_str()), Form("%s_berngaus",prefix.c_str()), *obs_var, *gaus, *soft_step);
    setCacheObsBinning(obs_var);
    bern_gaus->setBufferFraction(kBkgBufferFraction);
    obs_var->setRange(95, 180);
    return bern_gaus;
  } else if (order==6) {
	RooBernsteinFast<6> *bern = new RooBernsteinFast<6>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
     RooGenericPdf *soft_step = new RooGenericPdf(
      Form("%s_softstep",prefix.c_str()), Form("%s_softstep",prefix.c_str()),
      "1e-20 + 0.5*(1.0 + TMath::Erf((@0-@1)/(@2*sqrt(2.)))) * @3",
      RooArgList(*obs_var, *step_value, *step_width, *bern)
    );
    obs_var->setRange(-200, 500);
    setFftObsBinning(obs_var);
    RooFFTConvPdf *bern_gaus = new RooFFTConvPdf(Form("%s_berngaus",prefix.c_str()), Form("%s_berngaus",prefix.c_str()), *obs_var, *gaus, *soft_step);
    setCacheObsBinning(obs_var);
    bern_gaus->setBufferFraction(kBkgBufferFraction);
    obs_var->setRange(95, 180);
    return bern_gaus;
  } else {
	return NULL;
  }
}

RooAbsPdf* PdfModelBuilder::getBernsteinStepxGau(string prefix, int order, int mass_ALP){
  
  auto stableBern = makeStableStepGausWindow(mass_ALP, -1.0, 0.85, 0.75);
  double turnon_bern = stableBern.turnon, turnon_lbern = stableBern.turnonLo, turnon_hbern = stableBern.turnonHi;
  double sigma_bern  = stableBern.sigma,  sigma_lbern  = stableBern.sigmaLo,  sigma_hbern  = stableBern.sigmaHi;
  double width_bern  = stableBern.width,  width_lbern  = stableBern.widthLo,  width_hbern  = stableBern.widthHi;

  // Recent fits for these masses converged in narrow regions or pinned to bounds.
  // Use mass-specific seeds/ranges to stabilize the minimizer.
  if (mass_ALP == 14) {
    turnon_bern = 106.8; turnon_lbern = 103.0; turnon_hbern = 110.0;
    sigma_bern  = 4.2;   sigma_lbern  = 1.5;   sigma_hbern  = 7.5;
    width_bern  = 0.8;   width_lbern  = 0.3;   width_hbern  = 4.5;
  } else if (mass_ALP == 21) {
    turnon_bern = 106.9; turnon_lbern = 103.0; turnon_hbern = 111.0;
    sigma_bern  = 1.8;   sigma_lbern  = 0.8;   sigma_hbern  = 4.0;
    width_bern  = 0.8;   width_lbern  = 0.3;   width_hbern  = 4.5;
  } else if (mass_ALP == 23) {
    turnon_bern = 110.5; turnon_lbern = 108.0; turnon_hbern = 113.0;
    sigma_bern  = 4.5;   sigma_lbern  = 1.5;   sigma_hbern  = 9.0;
    width_bern  = 0.8;   width_lbern  = 0.3;   width_hbern  = 5.0;
  } else if (mass_ALP == 29) {
    turnon_bern = 110.5; turnon_lbern = 107.0; turnon_hbern = 114.0;
    sigma_bern  = 7.0;   sigma_lbern  = 3.0;   sigma_hbern  = 12.0;
    width_bern  = 2.0;   width_lbern  = 0.5;   width_hbern  = 6.0;
  }

  RooRealVar *g_mean  = new RooRealVar(Form("%s_gmean",prefix.c_str()),Form("%s_gmean",prefix.c_str()),0.);
  g_mean->setConstant(true);
  RooRealVar *g_sigma = new RooRealVar(Form("%s_gsigma",prefix.c_str()),Form("%s_gsigma",prefix.c_str()), sigma_bern, sigma_lbern, sigma_hbern);
  RooGaussian *gaus   = new RooGaussian(Form("%s_gaus",prefix.c_str()),Form("%s_gaus",prefix.c_str()),*obs_var,*g_mean,*g_sigma);

  RooRealVar *step_value = new RooRealVar(Form("%s_step",prefix.c_str()), Form("%s_step",prefix.c_str()), turnon_bern, turnon_lbern, turnon_hbern);
  RooRealVar *step_width = new RooRealVar(Form("%s_stepWidth",prefix.c_str()), Form("%s_stepWidth",prefix.c_str()), width_bern, width_lbern, width_hbern);

  RooArgList coeffList;
  std::vector<RooRealVar*> rawPars;
  std::vector<RooFormulaVar*> sqPars;
  for (int i=0; i<order; i++){
    std::string pname = Form("%s_b%02d",prefix.c_str(),i);
    RooRealVar* raw = new RooRealVar(pname.c_str(),pname.c_str(), 0.1*(i+1), -3., 3.);
    if (i == 0) fixUnitShapeCoeff(raw);
    rawPars.push_back(raw);
    std::string sname = Form("%s_sq_b%02d",prefix.c_str(),i);
    RooFormulaVar* sq = new RooFormulaVar(sname.c_str(),sname.c_str(),"@0*@0",RooArgList(*raw));
    sqPars.push_back(sq);
    coeffList.add(*sq);
  }

  RooAbsPdf* bern = nullptr;
  if      (order==1) bern = new RooBernsteinFast<1>(Form("%s_bern",prefix.c_str()),Form("%s_bern",prefix.c_str()),*obs_var,coeffList);
  else if (order==2) bern = new RooBernsteinFast<2>(Form("%s_bern",prefix.c_str()),Form("%s_bern",prefix.c_str()),*obs_var,coeffList);
  else if (order==3) bern = new RooBernsteinFast<3>(Form("%s_bern",prefix.c_str()),Form("%s_bern",prefix.c_str()),*obs_var,coeffList);
  else if (order==4) bern = new RooBernsteinFast<4>(Form("%s_bern",prefix.c_str()),Form("%s_bern",prefix.c_str()),*obs_var,coeffList);
  else if (order==5) bern = new RooBernsteinFast<5>(Form("%s_bern",prefix.c_str()),Form("%s_bern",prefix.c_str()),*obs_var,coeffList);
  else if (order==6) bern = new RooBernsteinFast<6>(Form("%s_bern",prefix.c_str()),Form("%s_bern",prefix.c_str()),*obs_var,coeffList);
  else return NULL;

  RooGenericPdf *soft_step_times_bern = new RooGenericPdf(
    Form("%s_softstepxbern",prefix.c_str()), Form("%s_softstepxbern",prefix.c_str()),
    "1e-30 + 0.5*(1.0 + TMath::Erf((@0-@1)/(@2*sqrt(2.)))) * @3",
    RooArgList(*obs_var, *step_value, *step_width, *bern)
  );

  setFftObsBinning(obs_var);
  RooFFTConvPdf *conv = new RooFFTConvPdf(
    Form("%s",prefix.c_str()), Form("%s",prefix.c_str()),
    *obs_var, *soft_step_times_bern, *gaus
  );
  setCacheObsBinning(obs_var);
  conv->setBufferFraction(kBkgBufferFraction);
  conv->setBufferFraction(0.15);

  return conv;
}

RooAbsPdf* PdfModelBuilder::getPowerLawStepxGau(string prefix, int order, int cat, int mass_ALP){
  if(order%2==0) return NULL;
  RooRealVar *mean = new RooRealVar(Form("%s_mean",prefix.c_str()),Form("%s_mean",prefix.c_str()),0.);
  mean->setConstant(true);
  double sigma_pow,sigma_lpow,sigma_hpow;
  double turnon_pow,turnon_lpow,turnon_hpow;
  double width_pow,width_lpow,width_hpow;
  double par1_pow1, par1_pow3, par3_pow3, par1_pow5, par3_pow5, par5_pow5;
  double par1_hpow1, par1_hpow3, par3_hpow3, par1_hpow5, par3_hpow5, par5_hpow5;
  double par1_lpow1, par1_lpow3, par3_lpow3, par1_lpow5, par3_lpow5, par5_lpow5;
  double coeff1_pow1, coeff1_pow3, coeff3_pow3, coeff1_pow5, coeff3_pow5, coeff5_pow5;
  double coeff1_hpow1, coeff1_hpow3, coeff3_hpow3, coeff1_hpow5, coeff3_hpow5, coeff5_hpow5;
  double coeff1_lpow1, coeff1_lpow3, coeff3_lpow3, coeff1_lpow5, coeff3_lpow5, coeff5_lpow5;

  coeff1_pow1 = 1.0;        coeff1_lpow1 = 0.;    coeff1_hpow1 = 10.;
  auto stablePow = makeStableStepGausWindow(mass_ALP, 0.6, 1.00, 1.00);
  par1_pow1 = -7.00;        par1_lpow1 = -12.;    par1_hpow1 = -5.;
  sigma_pow = stablePow.sigma;            sigma_lpow = stablePow.sigmaLo;      sigma_hpow = stablePow.sigmaHi;
  turnon_pow = stablePow.turnon;          turnon_lpow = stablePow.turnonLo;    turnon_hpow = stablePow.turnonHi;
  width_pow = stablePow.width;            width_lpow = stablePow.widthLo;      width_hpow = stablePow.widthHi;

  coeff1_pow3 = 1.0;      coeff1_lpow3 = 0.;    coeff1_hpow3 = 10.;
  coeff3_pow3 = 0.10;     coeff3_lpow3 = 1.e-4; coeff3_hpow3 = 1.;
  par1_pow3 = -6.80;      par1_lpow3 = -11.;    par1_hpow3 = -5.;
  par3_pow3 = -4.80;      par3_lpow3 = -8.;     par3_hpow3 = -2;
  sigma_pow = stablePow.sigma;          sigma_lpow = stablePow.sigmaLo;      sigma_hpow = stablePow.sigmaHi;
  turnon_pow = stablePow.turnon;        turnon_lpow = stablePow.turnonLo;    turnon_hpow = stablePow.turnonHi;

  coeff1_pow5 = 1.0;      coeff1_lpow5 = 0.;    coeff1_hpow5 = 10.;
  coeff3_pow5 = 0.10;     coeff3_lpow5 = 1.e-4; coeff3_hpow5 = 1.;
  coeff5_pow5 = 0.05;     coeff5_lpow5 = 1.e-4; coeff5_hpow5 = 1.;
  par1_pow5 = -6.80;      par1_lpow5 = -11.;    par1_hpow5 = -5.;
  par3_pow5 = -4.80;      par3_lpow5 = -8.;     par3_hpow5 = -2.;
  par5_pow5 = -5.80;      par5_lpow5 = -8.5;    par5_hpow5 = -1.5;
  sigma_pow = stablePow.sigma;          sigma_lpow = stablePow.sigmaLo;      sigma_hpow = stablePow.sigmaHi;
  turnon_pow = stablePow.turnon;        turnon_lpow = stablePow.turnonLo;    turnon_hpow = stablePow.turnonHi;

  if (mass_ALP == 14) {
    par1_pow1 = -8.8;      par1_lpow1 = -12.0;   par1_hpow1 = -6.0;
    par1_pow3 = -8.5;      par1_lpow3 = -12.0;   par1_hpow3 = -6.0;
    par3_pow3 = -5.5;      par3_lpow3 = -9.0;    par3_hpow3 = -2.0;
    par1_pow5 = -8.5;      par1_lpow5 = -12.0;   par1_hpow5 = -6.0;
    par3_pow5 = -5.5;      par3_lpow5 = -9.0;    par3_hpow5 = -2.0;
    par5_pow5 = -6.0;      par5_lpow5 = -9.0;    par5_hpow5 = -2.0;
    sigma_pow  = 4.1;      sigma_lpow  = 1.5;    sigma_hpow  = 8.0;
    turnon_pow = 114.4;    turnon_lpow = 111.0;  turnon_hpow = 118.0;
    width_pow  = 5.0;      width_lpow  = 1.0;    width_hpow  = 12.0;
  } else if (mass_ALP == 21) {
    par1_pow1 = -7.6;      par1_lpow1 = -10.0;   par1_hpow1 = -6.0;
    par1_pow3 = -7.5;      par1_lpow3 = -10.0;   par1_hpow3 = -6.0;
    par3_pow3 = -5.2;      par3_lpow3 = -8.0;    par3_hpow3 = -2.0;
    par1_pow5 = -7.5;      par1_lpow5 = -10.0;   par1_hpow5 = -6.0;
    par3_pow5 = -5.2;      par3_lpow5 = -8.0;    par3_hpow5 = -2.0;
    par5_pow5 = -6.2;      par5_lpow5 = -8.0;    par5_hpow5 = -2.0;
    sigma_pow  = 2.8;      sigma_lpow  = 1.0;    sigma_hpow  = 6.0;
    turnon_pow = 111.1;    turnon_lpow = 108.0;  turnon_hpow = 114.0;
    width_pow  = 3.0;      width_lpow  = 0.8;    width_hpow  = 8.0;
  } else if (mass_ALP == 23) {
    par1_pow1 = -7.4;      par1_lpow1 = -10.0;   par1_hpow1 = -5.5;
    par1_pow3 = -7.2;      par1_lpow3 = -10.0;   par1_hpow3 = -5.5;
    par3_pow3 = -4.8;      par3_lpow3 = -8.0;    par3_hpow3 = -2.0;
    par1_pow5 = -7.2;      par1_lpow5 = -10.0;   par1_hpow5 = -5.5;
    par3_pow5 = -4.8;      par3_lpow5 = -8.0;    par3_hpow5 = -2.0;
    par5_pow5 = -5.8;      par5_lpow5 = -8.0;    par5_hpow5 = -2.0;
    sigma_pow  = 5.5;      sigma_lpow  = 2.0;    sigma_hpow  = 10.0;
    turnon_pow = 110.5;    turnon_lpow = 108.0;  turnon_hpow = 113.0;
    width_pow  = 1.0;      width_lpow  = 0.5;    width_hpow  = 6.0;
  } else if (mass_ALP == 29) {
    par1_pow1 = -8.6;      par1_lpow1 = -11.0;   par1_hpow1 = -6.5;
    par1_pow3 = -8.3;      par1_lpow3 = -11.0;   par1_hpow3 = -6.0;
    par3_pow3 = -5.5;      par3_lpow3 = -8.5;    par3_hpow3 = -2.0;
    par1_pow5 = -8.3;      par1_lpow5 = -11.0;   par1_hpow5 = -6.0;
    par3_pow5 = -5.5;      par3_lpow5 = -8.5;    par3_hpow5 = -2.0;
    par5_pow5 = -6.4;      par5_lpow5 = -8.5;    par5_hpow5 = -2.0;
    sigma_pow  = 8.0;      sigma_lpow  = 4.0;    sigma_hpow  = 12.0;
    turnon_pow = 115.2;    turnon_lpow = 112.0;  turnon_hpow = 118.0;
    width_pow  = 1.0;      width_lpow  = 0.5;    width_hpow  = 6.0;
  }

  
  RooRealVar *sigma = new RooRealVar(Form("%s_sigma_p%d",prefix.c_str(),order),Form("%s_sigma_p%d",prefix.c_str(),order),sigma_pow,sigma_lpow,sigma_hpow);
  RooRealVar *turnon = new RooRealVar(Form("%s_turnon_p%d",prefix.c_str(),order),Form("%s_turnon_p%d",prefix.c_str(),order),turnon_pow,turnon_lpow,turnon_hpow);
  RooRealVar *width  = new RooRealVar(Form("%s_width_p%d",prefix.c_str(),order),Form("%s_width_p%d",prefix.c_str(),order),width_pow,width_lpow,width_hpow);
  
    if (order==1) {
      RooRealVar *p1          = new RooRealVar(Form("%s_p1_pow1",prefix.c_str()),Form("%s_p1_pow1",prefix.c_str()),par1_pow1,par1_lpow1,par1_hpow1);
      RooRealVar *cp1         = new RooRealVar(Form("%s_cp1_pow1",prefix.c_str()),Form("%s_cp1_pow1",prefix.c_str()),coeff1_pow1,coeff1_lpow1,coeff1_hpow1);
      // [PZ-FIX] cp1 is an overall scale factor for an internally-normalized pdf -> unconstrained.
      // Fix it to 1 to avoid Minuit2 "2nd derivative zero"/invalid Hessian warnings.
      fixUnitShapeCoeff(cp1);
      RooGenericPdf *soft_step = new RooGenericPdf(
        Form("%s_soft_pow1",prefix.c_str()),Form("%s_soft_pow1",prefix.c_str()),
        "1e-20+0.5*(1.0+TMath::Erf((@0-@1)/(@2*sqrt(2.))))*(@4*(@0)^(@3))",
        RooArgList(*obs_var,*turnon,*width,*p1,*cp1)
      );
      RooGaussModel *gau      = new RooGaussModel(Form("%s_gau_pow1",prefix.c_str()),Form("%s_gau_pow1",prefix.c_str()),*obs_var,*mean,*sigma);
      setFftObsBinning(obs_var);
      RooFFTConvPdf *gauxpow  = new RooFFTConvPdf(Form("%s1",prefix.c_str()),Form("%s_gauxpow1",prefix.c_str()),*obs_var,*soft_step,*gau);
      setCacheObsBinning(obs_var);
      gauxpow->setBufferFraction(kBkgBufferFraction);
      return gauxpow;
  } else if (order==3) {
      RooRealVar *p1          = new RooRealVar(Form("%s_p1_pow3",prefix.c_str()),Form("%s_p1_pow3",prefix.c_str()),par1_pow3,par1_lpow3,par1_hpow3);
      RooRealVar *cp1         = new RooRealVar(Form("%s_cp1_pow3",prefix.c_str()),Form("%s_cp1_pow3",prefix.c_str()),coeff1_pow3,coeff1_lpow3,coeff1_hpow3);
      RooRealVar *p3          = new RooRealVar(Form("%s_p3_pow3",prefix.c_str()),Form("%s_p3_pow3",prefix.c_str()),par3_pow3,par3_lpow3,par3_hpow3);
      RooRealVar *cp3         = new RooRealVar(Form("%s_cp3_pow3",prefix.c_str()),Form("%s_cp3_pow3",prefix.c_str()),coeff3_pow3,coeff3_lpow3,coeff3_hpow3);
      fixUnitShapeCoeff(cp1);
      RooGenericPdf *soft_step = new RooGenericPdf(
        Form("%s_soft_pow3",prefix.c_str()),Form("%s_soft_pow3",prefix.c_str()),
        "1e-20+0.5*(1.0+TMath::Erf((@0-@1)/(@2*sqrt(2.))))*(@4*(@0)^(@3)+@6*(@0)^(@5))",
        RooArgList(*obs_var,*turnon,*width,*p1,*cp1,*p3,*cp3)
      );
      RooGaussModel *gau      = new RooGaussModel(Form("%s_gau_pow3",prefix.c_str()),Form("%s_gau_pow3",prefix.c_str()),*obs_var,*mean,*sigma);
      setFftObsBinning(obs_var);
      RooFFTConvPdf *gauxpow  = new RooFFTConvPdf(Form("%s3",prefix.c_str()),Form("%s_gauxpow3",prefix.c_str()),*obs_var,*soft_step,*gau);
      setCacheObsBinning(obs_var);
      gauxpow->setBufferFraction(kBkgBufferFraction);
      return gauxpow;
  } else if (order==5) {
      RooRealVar *p1          = new RooRealVar(Form("%s_p1_pow5",prefix.c_str()),Form("%s_p1_pow5",prefix.c_str()),par1_pow5,par1_lpow5,par1_hpow5);
      RooRealVar *cp1         = new RooRealVar(Form("%s_cp1_pow5",prefix.c_str()),Form("%s_cp1_pow5",prefix.c_str()),coeff1_pow5,coeff1_lpow5,coeff1_hpow5);
      RooRealVar *p3          = new RooRealVar(Form("%s_p3_pow5",prefix.c_str()),Form("%s_p3_pow5",prefix.c_str()),par3_pow5,par3_lpow5, par3_hpow5);
      RooRealVar *cp3         = new RooRealVar(Form("%s_cp3_pow5",prefix.c_str()),Form("%s_cp3_pow5",prefix.c_str()),coeff3_pow5,coeff3_lpow5,coeff3_hpow5);
      RooRealVar *p5          = new RooRealVar(Form("%s_p5_pow5",prefix.c_str()),Form("%s_p5_pow5",prefix.c_str()),par5_pow5,par5_lpow5,par5_hpow5);
      RooRealVar *cp5         = new RooRealVar(Form("%s_cp5_pow5",prefix.c_str()),Form("%s_cp5_pow5",prefix.c_str()),coeff5_pow5,coeff5_lpow5,coeff5_hpow5);
      fixUnitShapeCoeff(cp1);
    	//RooGenericPdf *step     = new RooGenericPdf(Form("%s_step_pow5",prefix.c_str()),Form("%s_step_pow5",prefix.c_str()), "1e-20+(@0 > @1)*(@3*(@0)^(@2)+@5*(@0)^(@4)+@7*(@0)^(@6))", RooArgList(*obs_var,*turnon,*p1,*cp1,*p3,*cp3,*p5,*cp5));
      RooGenericPdf *soft_step = new RooGenericPdf(
        Form("%s_soft_pow5",prefix.c_str()),Form("%s_soft_pow5",prefix.c_str()),
        "1e-20+0.5*(1.0+TMath::Erf((@0-@1)/(@2*sqrt(2.))))*(@4*(@0)^(@3)+@6*(@0)^(@5)+@8*(@0)^(@7))",
        RooArgList(*obs_var,*turnon,*width,*p1,*cp1,*p3,*cp3,*p5,*cp5)
      );
      RooGaussModel *gau      = new RooGaussModel(Form("%s_gau_pow5",prefix.c_str()),Form("%s_gau_pow5",prefix.c_str()),*obs_var,*mean,*sigma);
      setFftObsBinning(obs_var);
      RooFFTConvPdf *gauxpow  = new RooFFTConvPdf(Form("%s5",prefix.c_str()),Form("%s_gauxpow5",prefix.c_str()),*obs_var,*soft_step,*gau);
      setCacheObsBinning(obs_var);
      gauxpow->setBufferFraction(kBkgBufferFraction);
      return gauxpow;
  } 
   else {
	return NULL;
  }
}

RooAbsPdf* PdfModelBuilder::getPowerLawGeneric(string prefix, int order){
  
  if (order%2==0){
    cerr << "ERROR -- addPowerLaw -- only odd number of params allowed" << endl;
    return NULL;
  }
  else {
    int nfracs=(order-1)/2;
    int npows=order-nfracs;
    assert(nfracs==npows-1);
    string formula="";
    RooArgList *dependents = new RooArgList();
    dependents->add(*obs_var);
    if (order>1) {
      formula += "(1.-";
      for (int i=1; i<=nfracs; i++){
        if (i<nfracs) formula += Form("@%d-",i);
        else formula += Form("@%d)*",i);
        string name =  Form("%s_f%d",prefix.c_str(),i);
        params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),0.1,0.,1.)));
        dependents->add(*params[name]);
      }
    }
    for (int i=1; i<=npows; i++){
      string pname =  Form("%s_p%d",prefix.c_str(),i);
      string fname =  Form("%s_f%d",prefix.c_str(),i-1);
      params.insert(pair<string,RooRealVar*>(pname, new RooRealVar(pname.c_str(),pname.c_str(),TMath::Max(-10.,-2.*(i+1)),-10.,0.)));
      if (i==1) {
        formula += Form("TMath::Power(@0,@%d)",nfracs+i);
        dependents->add(*params[pname]);
      }
      else {
        formula += Form(" + @%d*TMath::Power(@0,@%d)",i-1,nfracs+i);
        dependents->add(*params[pname]);
      }
    }
    cout << "FORM -- " << formula << endl;
    dependents->Print("v");
    RooGenericPdf *pow = new RooGenericPdf(prefix.c_str(),prefix.c_str(),formula.c_str(),*dependents);
    pow->Print("v");
    return pow;
  }
}

RooAbsPdf* PdfModelBuilder::getPowerLaw(string prefix, int order){
  
  RooArgList coefList;
  for (int i=0; i<order; i++){
    double start=-2.;
    double low=-10.;
    double high=0.;
    if (order>0){
      start=-0.001/double(i);
      low=-0.01;
      high=0.01;
    }
    RooRealVar *var = new RooRealVar(Form("%s_p%d",prefix.c_str(),i),Form("%s_p%d",prefix.c_str(),i),start,low,high);
    coefList.add(*var);
  }
  RooPowerLawSum *pow = new RooPowerLawSum(prefix.c_str(),prefix.c_str(),*obs_var,coefList);
  return pow;
}

RooAbsPdf* PdfModelBuilder::getExponential(string prefix, int order){
  
  RooArgList coefList;
  for (int i=0; i<order; i++){
    double start=-1.;
    double low=-2.;
    double high=0.;
    if (order>0){
      start=-0.001/double(i);
      low=-0.01;
      high=0.01;
    }
    RooRealVar *var = new RooRealVar(Form("%s_p%d",prefix.c_str(),i),Form("%s_p%d",prefix.c_str(),i),start,low,high);
    coefList.add(*var);
  }
  RooPowerLawSum *exp = new RooPowerLawSum(prefix.c_str(),prefix.c_str(),*obs_var,coefList);
  return exp;
}

RooAbsPdf* PdfModelBuilder::getExponentialStepxGau(string prefix, int order, int cat, int mass_ALP){
  if(order%2==0) return NULL;

  RooRealVar *mean = new RooRealVar(Form("%s_mean",prefix.c_str()),Form("%s_mean",prefix.c_str()),0.);
  mean->setConstant(true);
  double sigma_exp,sigma_lexp,sigma_hexp;
  double turnon_exp,turnon_lexp,turnon_hexp;
  double width_exp,width_lexp,width_hexp;
  double par1_exp1, par1_exp3, par3_exp3, par1_exp5, par3_exp5, par5_exp5;
  double par1_hexp1, par1_hexp3, par3_hexp3, par1_hexp5, par3_hexp5, par5_hexp5;
  double par1_lexp1, par1_lexp3, par3_lexp3, par1_lexp5, par3_lexp5, par5_lexp5;
  double coeff1_exp1, coeff1_exp3, coeff3_exp3, coeff1_exp5, coeff3_exp5, coeff5_exp5;
  double coeff1_hexp1, coeff1_hexp3, coeff3_hexp3, coeff1_hexp5, coeff3_hexp5, coeff5_hexp5;
  double coeff1_lexp1, coeff1_lexp3, coeff3_lexp3, coeff1_lexp5, coeff3_lexp5, coeff5_lexp5;
  
  coeff1_exp1 = 1.0;      coeff1_lexp1 = 0.;      coeff1_hexp1 = 10.;
  auto stableExp = makeStableStepGausWindow(mass_ALP, 0.4, 0.90, 1.15);
  par1_exp1 = -0.055;     par1_lexp1 = -0.10;     par1_hexp1 = -0.02;
  sigma_exp = stableExp.sigma;        sigma_lexp = stableExp.sigmaLo;         sigma_hexp = stableExp.sigmaHi;
  turnon_exp = stableExp.turnon;      turnon_lexp = stableExp.turnonLo;       turnon_hexp = stableExp.turnonHi;
  width_exp = stableExp.width;        width_lexp = stableExp.widthLo;         width_hexp = stableExp.widthHi;

  coeff1_exp3 = 1.0;        coeff1_lexp3 = 0.;    coeff1_hexp3 = 10.;
  coeff3_exp3 = 0.25;       coeff3_lexp3 = 1.e-4; coeff3_hexp3 = 1.;
  par1_exp3 = -0.055;       par1_lexp3 = -0.12;   par1_hexp3 = -0.02;
  par3_exp3 = -0.040;       par3_lexp3 = -0.10;   par3_hexp3 = -0.015;

  coeff1_exp5 = 1.0;        coeff1_lexp5 = 0.;    coeff1_hexp5 = 10.;
  coeff3_exp5 = 0.20;       coeff3_lexp5 = 1.e-4; coeff3_hexp5 = 1.;
  coeff5_exp5 = 0.05;       coeff5_lexp5 = 1.e-4; coeff5_hexp5 = 1.;
  par1_exp5 = -0.055;       par1_lexp5 = -0.12;   par1_hexp5 = -0.015;
  par3_exp5 = -0.038;       par3_lexp5 = -0.10;   par3_hexp5 = -0.015;
  par5_exp5 = -0.018;       par5_lexp5 = -0.06;   par5_hexp5 = -0.008;

  if (mass_ALP == 14) {
    par1_exp1 = -0.0566;    par1_lexp1 = -0.09;   par1_hexp1 = -0.03;
    par1_exp3 = -0.056;     par1_lexp3 = -0.15;   par1_hexp3 = -0.03;
    par3_exp3 = -0.040;     par3_lexp3 = -0.12;   par3_hexp3 = -0.02;
    par1_exp5 = -0.060;     par1_lexp5 = -0.15;   par1_hexp5 = -0.02;
    par3_exp5 = -0.040;     par3_lexp5 = -0.12;   par3_hexp5 = -0.02;
    par5_exp5 = -0.020;     par5_lexp5 = -0.08;   par5_hexp5 = -0.01;
    sigma_exp  = 1.1;       sigma_lexp  = 0.6;    sigma_hexp  = 4.0;
    turnon_exp = 113.5;     turnon_lexp = 111.0;  turnon_hexp = 116.0;
    width_exp  = 5.0;       width_lexp  = 1.0;    width_hexp  = 10.0;
  } else if (mass_ALP == 21) {
    par1_exp1 = -0.0544;    par1_lexp1 = -0.09;   par1_hexp1 = -0.03;
    par1_exp3 = -0.054;     par1_lexp3 = -0.15;   par1_hexp3 = -0.03;
    par3_exp3 = -0.040;     par3_lexp3 = -0.12;   par3_hexp3 = -0.02;
    par1_exp5 = -0.055;     par1_lexp5 = -0.15;   par1_hexp5 = -0.02;
    par3_exp5 = -0.040;     par3_lexp5 = -0.12;   par3_hexp5 = -0.02;
    par5_exp5 = -0.020;     par5_lexp5 = -0.08;   par5_hexp5 = -0.01;
    sigma_exp  = 1.1;       sigma_lexp  = 0.6;    sigma_hexp  = 4.0;
    turnon_exp = 110.8;     turnon_lexp = 109.0;  turnon_hexp = 113.0;
    width_exp  = 3.7;       width_lexp  = 0.8;    width_hexp  = 8.0;
  } else if (mass_ALP == 23) {
    par1_exp1 = -0.0566;    par1_lexp1 = -0.08;   par1_hexp1 = -0.03;
    par1_exp3 = -0.056;     par1_lexp3 = -0.12;   par1_hexp3 = -0.03;
    par3_exp3 = -0.040;     par3_lexp3 = -0.10;   par3_hexp3 = -0.02;
    par1_exp5 = -0.056;     par1_lexp5 = -0.12;   par1_hexp5 = -0.02;
    par3_exp5 = -0.040;     par3_lexp5 = -0.10;   par3_hexp5 = -0.02;
    par5_exp5 = -0.020;     par5_lexp5 = -0.08;   par5_hexp5 = -0.01;
    sigma_exp  = 6.4;       sigma_lexp  = 3.0;    sigma_hexp  = 10.0;
    turnon_exp = 110.2;     turnon_lexp = 108.5;  turnon_hexp = 112.5;
    width_exp  = 1.0;       width_lexp  = 0.5;    width_hexp  = 6.0;
  } else if (mass_ALP == 29) {
    par1_exp1 = -0.0598;    par1_lexp1 = -0.09;   par1_hexp1 = -0.03;
    par1_exp3 = -0.060;     par1_lexp3 = -0.12;   par1_hexp3 = -0.03;
    par3_exp3 = -0.040;     par3_lexp3 = -0.10;   par3_hexp3 = -0.02;
    par1_exp5 = -0.060;     par1_lexp5 = -0.12;   par1_hexp5 = -0.02;
    par3_exp5 = -0.040;     par3_lexp5 = -0.10;   par3_hexp5 = -0.02;
    par5_exp5 = -0.020;     par5_lexp5 = -0.08;   par5_hexp5 = -0.01;
    sigma_exp  = 6.8;       sigma_lexp  = 3.0;    sigma_hexp  = 10.0;
    turnon_exp = 115.5;     turnon_lexp = 113.0;  turnon_hexp = 118.0;
    width_exp  = 4.8;       width_lexp  = 1.0;    width_hexp  = 9.0;
  }

  
  RooRealVar *sigma = new RooRealVar(Form("%s_sigma_p%d",prefix.c_str(),order),Form("%s_sigma_p%d",prefix.c_str(),order),sigma_exp,sigma_lexp,sigma_hexp);
  RooRealVar *turnon = new RooRealVar(Form("%s_turnon_p%d",prefix.c_str(),order),Form("%s_turnon_p%d",prefix.c_str(),order),turnon_exp,turnon_lexp,turnon_hexp);
  RooRealVar *width  = new RooRealVar(Form("%s_width_p%d",prefix.c_str(),order),Form("%s_width_p%d",prefix.c_str(),order),width_exp,width_lexp,width_hexp);
  
    if (order==1) {
      RooRealVar *p1 = new RooRealVar(Form("%s_p1_exp1",prefix.c_str()),Form("%s_p1_exp1",prefix.c_str()),par1_exp1,par1_lexp1,par1_hexp1);
      RooRealVar *cp1 = new RooRealVar(Form("%s_cp1_exp1",prefix.c_str()),Form("%s_cp1_exp1",prefix.c_str()),coeff1_exp1,coeff1_lexp1,coeff1_hexp1);
      // [PZ-FIX] cp1 is an overall scale factor for an internally-normalized pdf -> unconstrained.
      // Fix it to 1 to avoid Minuit2 "2nd derivative zero"/invalid Hessian warnings.
      fixUnitShapeCoeff(cp1);
      RooGenericPdf *soft_step = new RooGenericPdf(
        Form("%s_soft_exp1",prefix.c_str()),Form("%s_soft_exp1",prefix.c_str()),
        "1e-20+0.5*(1.0+TMath::Erf((@0-@1)/(@2*sqrt(2.))))*(@4*TMath::Exp(@0*@3))",
        RooArgList(*obs_var,*turnon,*width,*p1,*cp1)
      );
      RooGaussModel *gau = new RooGaussModel(Form("%s_gau_exp1",prefix.c_str()),Form("%s_gau_exp1",prefix.c_str()),*obs_var,*mean,*sigma);
      setFftObsBinning(obs_var);
      RooFFTConvPdf *gauxexp = new RooFFTConvPdf(Form("%s1",prefix.c_str()),Form("%s_gauxexp1",prefix.c_str()),*obs_var,*soft_step,*gau);
      setCacheObsBinning(obs_var);
      gauxexp->setBufferFraction(kBkgBufferFraction);
      return gauxexp;
  } else if (order==3) {
      RooRealVar *p1 = new RooRealVar(Form("%s_p1_exp3",prefix.c_str()),Form("%s_p1_exp3",prefix.c_str()),par1_exp3,par1_lexp3, par1_hexp3);
      RooRealVar *cp1 = new RooRealVar(Form("%s_cp1_exp3",prefix.c_str()),Form("%s_cp1_exp3",prefix.c_str()),coeff1_exp3,coeff1_lexp3,coeff1_hexp3);
      RooRealVar *p3 = new RooRealVar(Form("%s_p3_exp3",prefix.c_str()),Form("%s_p3_exp3",prefix.c_str()),par3_exp3,par3_lexp3, par3_hexp3);
      RooRealVar *cp3 = new RooRealVar(Form("%s_cp3_exp3",prefix.c_str()),Form("%s_cp3_exp3",prefix.c_str()),coeff3_exp3,coeff3_lexp3,coeff3_hexp3);
      fixUnitShapeCoeff(cp1);
      RooGenericPdf *soft_step = new RooGenericPdf(
        Form("%s_soft_exp3",prefix.c_str()),Form("%s_soft_exp3",prefix.c_str()),
        "1e-20+0.5*(1.0+TMath::Erf((@0-@1)/(@2*sqrt(2.))))*(@4*TMath::Exp(@0*@3)+@6*TMath::Exp(@0*@5))",
        RooArgList(*obs_var,*turnon,*width,*p1,*cp1,*p3,*cp3)
      );
      RooGaussModel *gau = new RooGaussModel(Form("%s_gau_exp3",prefix.c_str()),Form("%s_gau_exp3",prefix.c_str()),*obs_var,*mean,*sigma);
      setFftObsBinning(obs_var);
      RooFFTConvPdf *gauxexp = new RooFFTConvPdf(Form("%s3",prefix.c_str()),Form("%s_gauxexp3",prefix.c_str()),*obs_var,*soft_step,*gau);
      setCacheObsBinning(obs_var);
      gauxexp->setBufferFraction(kBkgBufferFraction);
      return gauxexp;
  } else if (order==5) {
      RooRealVar *p1 = new RooRealVar(Form("%s_p1_exp5",prefix.c_str()),Form("%s_p1_exp5",prefix.c_str()),par1_exp5,par1_lexp5, par1_hexp5);
      RooRealVar *cp1 = new RooRealVar(Form("%s_cp1_exp5",prefix.c_str()),Form("%s_cp1_exp5",prefix.c_str()),coeff1_exp5,coeff1_lexp5,coeff1_hexp5);
      RooRealVar *p3 = new RooRealVar(Form("%s_p3_exp5",prefix.c_str()),Form("%s_p3_exp5",prefix.c_str()),par3_exp5,par3_lexp5, par3_hexp5);
      RooRealVar *cp3 = new RooRealVar(Form("%s_cp3_exp5",prefix.c_str()),Form("%s_cp3_exp5",prefix.c_str()),coeff3_exp5,coeff3_lexp5,coeff3_hexp5);
      RooRealVar *p5 = new RooRealVar(Form("%s_p5_exp5",prefix.c_str()),Form("%s_p5_exp5",prefix.c_str()),par5_exp5,par5_lexp5, par5_hexp5);
      RooRealVar *cp5 = new RooRealVar(Form("%s_cp5_exp5",prefix.c_str()),Form("%s_cp5_exp5",prefix.c_str()),coeff5_exp5,coeff5_lexp5,coeff5_hexp5);
      fixUnitShapeCoeff(cp1);
      RooGenericPdf *soft_step = new RooGenericPdf(
        Form("%s_soft_exp5",prefix.c_str()),Form("%s_soft_exp5",prefix.c_str()),
        "1e-20+0.5*(1.0+TMath::Erf((@0-@1)/(@2*sqrt(2.))))*(@4*TMath::Exp(@0*@3)+@6*TMath::Exp(@0*@5)+@8*TMath::Exp(@0*@7))",
        RooArgList(*obs_var,*turnon,*width,*p1,*cp1,*p3,*cp3,*p5,*cp5)
      );
      RooGaussModel *gau = new RooGaussModel(Form("%s_gau_exp5",prefix.c_str()),Form("%s_gau_exp5",prefix.c_str()),*obs_var,*mean,*sigma);
      setFftObsBinning(obs_var);
      RooFFTConvPdf *gauxexp = new RooFFTConvPdf(Form("%s5",prefix.c_str()),Form("%s_gauxexp5",prefix.c_str()),*obs_var,*soft_step,*gau);
      setCacheObsBinning(obs_var);
      gauxexp->setBufferFraction(kBkgBufferFraction);
      return gauxexp;
  } 
   else {
	return NULL;
  }
}

RooAbsPdf* PdfModelBuilder::getPowerLawSingle(string prefix, int order){
  
  if (order%2==0){
    cerr << "ERROR -- addPowerLaw -- only odd number of params allowed" << endl;
    return NULL;
  }
  else {
    int nfracs=(order-1)/2;
    int npows=order-nfracs;
    assert(nfracs==npows-1);
    RooArgList *fracs = new RooArgList();
    RooArgList *pows = new RooArgList();
    for (int i=1; i<=nfracs; i++){
      string name =  Form("%s_f%d",prefix.c_str(),i);
      params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),0.9-float(i-1)*1./nfracs,0.,1.)));
      fracs->add(*params[name]);
    }
    for (int i=1; i<=npows; i++){
      string name =  Form("%s_p%d",prefix.c_str(),i);
      string ename =  Form("%s_e%d",prefix.c_str(),i);
      params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),TMath::Max(-9.,-1.*(i+1)),-9.,1.)));
      utilities.insert(pair<string,RooAbsPdf*>(ename, new RooPower(ename.c_str(),ename.c_str(),*obs_var,*params[name])));
      pows->add(*utilities[ename]);
    }
    RooAbsPdf *pow = new RooAddPdf(prefix.c_str(),prefix.c_str(),*pows,*fracs,true); 
    return pow;

    //bing
    RooRealVar *fgaus = new RooRealVar("fgaus", "gaus fraction",0.5,0.,1.) ;
    RooRealVar *mean = new RooRealVar("mean","mean",0,-10.0,10.0) ;
    mean->setConstant(true);
    RooRealVar *sigma = new RooRealVar("sigma","sigma",1,0.,10.) ;
    RooGaussian *gaus = new RooGaussian("gaus","gaus",*obs_var,*mean,*sigma) ;

    RooRealVar *step_value = new RooRealVar("step_value", "step value", 115., 110., 130.);
    RooGenericPdf *step_func = new RooGenericPdf("step_func", "step_func", "1e-20+( @0 > @1) * @2", RooArgSet(*obs_var, *step_value, *pow));
    obs_var->setRange(-400.0,500.0);
    setFftObsBinning(obs_var);
    RooFFTConvPdf *pow_gaus = new RooFFTConvPdf("pow_gaus", "pow_gaus", *obs_var, *gaus, *step_func);
    setCacheObsBinning(obs_var);
    pow_gaus->setBufferFraction(kBkgBufferFraction);
    obs_var->setRange(110.0,180.0);
  }
}

RooAbsPdf* PdfModelBuilder::getLaurentSeries(string prefix, int order){
 
  int nlower=int(ceil(order/2.));
  int nhigher=order-nlower;
  RooArgList *pows = new RooArgList();
  RooArgList *plist = new RooArgList();
  string pname =  Form("%s_pow0",prefix.c_str());
  utilities.insert(pair<string,RooAbsPdf*>(pname, new RooPower(pname.c_str(),pname.c_str(),*obs_var,RooConst(-4.))));
  pows->add(*utilities[pname]);

  // even terms
  for (int i=1; i<=nlower; i++){
    string name = Form("%s_l%d",prefix.c_str(),i);
    params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),0.25/order,0.000001,0.999999)));
    plist->add(*params[name]);
    string pname =  Form("%s_powl%d",prefix.c_str(),i);
    utilities.insert(pair<string,RooAbsPdf*>(pname, new RooPower(pname.c_str(),pname.c_str(),*obs_var,RooConst(-4.-i))));
    pows->add(*utilities[pname]);
  }
  // odd terms
  for (int i=1; i<=nhigher; i++){
    string name = Form("%s_h%d",prefix.c_str(),i);
    params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),0.25/order,0.000001,0.999999)));
    plist->add(*params[name]);
    string pname =  Form("%s_powh%d",prefix.c_str(),i);
    utilities.insert(pair<string,RooAbsPdf*>(pname, new RooPower(pname.c_str(),pname.c_str(),*obs_var,RooConst(-4.+i))));
    pows->add(*utilities[pname]);
  }
  RooAddPdf *pdf = new RooAddPdf(prefix.c_str(),prefix.c_str(),*pows,*plist,true);
  return pdf;

    RooRealVar *mean1 = new RooRealVar("mean1","mean1",0.0) ;
    RooRealVar *sigma1 = new RooRealVar("sigma1","sigma1",5,-10.,20.) ;
    RooGaussian *gaus1 = new RooGaussian("gaus1","gaus1",*obs_var,*mean1,*sigma1) ;
    RooRealVar *step_value1 = new RooRealVar("step_value1", "step_value1",115.,110.,130.) ;
    RooGenericPdf *step_func1 = new RooGenericPdf("step_func1","step_func1","(1e-20+( @0 > @1)) * @2",RooArgSet(*obs_var,*step_value1,*pdf));
    obs_var->setRange(-400.0,500.0);
    setFftObsBinning(obs_var);
    RooFFTConvPdf *pdf_gaus = new RooFFTConvPdf("pdf_gaus","pdf_gaus", *obs_var, *gaus1, *step_func1);
    setCacheObsBinning(obs_var);
    pdf_gaus->setBufferFraction(kBkgBufferFraction);
    obs_var->setRange(110.0,180.0);
}

RooAbsPdf* PdfModelBuilder::getLaurentStepxGau(string prefix, int order, int cat, int mass_ALP){
  
  if(order>5) return NULL;

  RooRealVar *mean = new RooRealVar(Form("%s_mean",prefix.c_str()),Form("%s_mean",prefix.c_str()),0.);
  mean->setConstant(true);
  double sigma_lau,sigma_llau,sigma_hlau;
  double turnon_lau,turnon_llau,turnon_hlau;
  double width_lau,width_llau,width_hlau;
  double coeff1_lau1,   coeff1_lau2,  coeff2_lau2,    coeff1_lau3,  coeff2_lau3,  coeff3_lau3,    coeff1_lau4,  coeff2_lau4,  coeff3_lau4,  coeff4_lau4;
  double coeff1_hlau1,  coeff1_hlau2, coeff2_hlau2,   coeff1_hlau3, coeff2_hlau3, coeff3_hlau3,   coeff1_hlau4, coeff2_hlau4, coeff3_hlau4, coeff4_hlau4;
  double coeff1_llau1,  coeff1_llau2, coeff2_llau2,   coeff1_llau3, coeff2_llau3, coeff3_llau3,   coeff1_llau4, coeff2_llau4, coeff3_llau4, coeff4_llau4;
 
  coeff1_lau1 = 1.0;        coeff1_llau1 = 0.;      coeff1_hlau1 = 10.;
  auto stableLau = makeStableStepGausWindow(mass_ALP, 0.2, 1.00, 0.95);
  sigma_lau = stableLau.sigma;          sigma_llau = stableLau.sigmaLo;       sigma_hlau = stableLau.sigmaHi;
  turnon_lau = stableLau.turnon;        turnon_llau = stableLau.turnonLo;     turnon_hlau = stableLau.turnonHi;
  width_lau = stableLau.width;          width_llau = stableLau.widthLo;       width_hlau = stableLau.widthHi;

  coeff1_lau2 = 1.0;        coeff1_llau2 = 0.;      coeff1_hlau2 = 10.;
  coeff2_lau2 = 0.10;       coeff2_llau2 = 1.e-4;   coeff2_hlau2 = 0.5;
  sigma_lau = stableLau.sigma;          sigma_llau = stableLau.sigmaLo;      sigma_hlau = stableLau.sigmaHi;
  turnon_lau = stableLau.turnon;        turnon_llau = stableLau.turnonLo;    turnon_hlau = stableLau.turnonHi;
 
  coeff1_lau3 = 1.0;        coeff1_llau3 = 0.;      coeff1_hlau3 = 10.;
  coeff2_lau3 = 0.10;       coeff2_llau3 = 1.e-4;   coeff2_hlau3 = 0.5;
  coeff3_lau3 = 0.05;       coeff3_llau3 = 1.e-4;   coeff3_hlau3 = 2.0;
  sigma_lau = stableLau.sigma;         sigma_llau = stableLau.sigmaLo;      sigma_hlau = stableLau.sigmaHi;
  turnon_lau = stableLau.turnon;       turnon_llau = stableLau.turnonLo;   turnon_hlau = stableLau.turnonHi;

  coeff1_lau4 = 1.0;        coeff1_llau4 = 0.;      coeff1_hlau4 = 10.;
  coeff2_lau4 = 0.10;       coeff2_llau4 = 1.e-4;   coeff2_hlau4 = 0.5;
  coeff3_lau4 = 0.05;       coeff3_llau4 = 1.e-4;   coeff3_hlau4 = 2.0;
  coeff4_lau4 = 0.02;       coeff4_llau4 = 1.e-4;   coeff4_hlau4 = 2.0;
  sigma_lau = stableLau.sigma;         sigma_llau = stableLau.sigmaLo;      sigma_hlau = stableLau.sigmaHi;
  turnon_lau = stableLau.turnon;       turnon_llau = stableLau.turnonLo;     turnon_hlau = stableLau.turnonHi;

  if (mass_ALP == 14) {
    sigma_lau  = 4.8;       sigma_llau  = 1.5;     sigma_hlau  = 8.0;
    turnon_lau = 109.0;     turnon_llau = 106.0;   turnon_hlau = 112.0;
    width_lau  = 0.8;       width_llau  = 0.3;     width_hlau  = 5.0;
  } else if (mass_ALP == 21) {
    sigma_lau  = 3.0;       sigma_llau  = 1.0;     sigma_hlau  = 6.0;
    turnon_lau = 110.5;     turnon_llau = 108.0;   turnon_hlau = 113.0;
    width_lau  = 3.0;       width_llau  = 0.8;     width_hlau  = 8.0;
  } else if (mass_ALP == 23) {
    sigma_lau  = 5.0;       sigma_llau  = 2.0;     sigma_hlau  = 9.0;
    turnon_lau = 110.8;     turnon_llau = 108.0;   turnon_hlau = 113.0;
    width_lau  = 1.0;       width_llau  = 0.5;     width_hlau  = 6.0;
  } else if (mass_ALP == 29) {
    sigma_lau  = 7.4;       sigma_llau  = 3.0;     sigma_hlau  = 11.0;
    turnon_lau = 113.9;     turnon_llau = 111.0;   turnon_hlau = 117.0;
    width_lau  = 1.0;       width_llau  = 0.5;     width_hlau  = 6.0;
  }


  RooRealVar *sigma = new RooRealVar(Form("%s_sigma_p%d",prefix.c_str(),order),Form("%s_sigma_p%d",prefix.c_str(),order),sigma_lau,sigma_llau,sigma_hlau);
  RooRealVar *turnon = new RooRealVar(Form("%s_turnon_p%d",prefix.c_str(),order),Form("%s_turnon_p%d",prefix.c_str(),order),turnon_lau,turnon_llau,turnon_hlau);
  RooRealVar *width  = new RooRealVar(Form("%s_width_p%d",prefix.c_str(),order),Form("%s_width_p%d",prefix.c_str(),order),width_lau,width_llau,width_hlau);
  
  if (order==1) {
      RooRealVar *cp1 = new RooRealVar(Form("%s_cp1_lau1",prefix.c_str()),Form("%s_cp1_lau1",prefix.c_str()),coeff1_lau1,coeff1_llau1,coeff1_hlau1);
      // [PZ-FIX] cp1 is an overall scale factor for an internally-normalized pdf -> unconstrained.
      // Fix it to 1 to avoid Minuit2 "2nd derivative zero"/invalid Hessian warnings.
      fixUnitShapeCoeff(cp1);
      RooGenericPdf *soft_step = new RooGenericPdf(
        Form("%s_soft_lau1",prefix.c_str()),Form("%s_soft_lau1",prefix.c_str()),
        "1e-20+0.5*(1.0+TMath::Erf((@0-@1)/(@2*sqrt(2.))))*(@3*(@0)^(-4))",
        RooArgList(*obs_var,*turnon,*width,*cp1)
      );
      RooGaussModel *gau = new RooGaussModel(Form("%s_gau_lau1",prefix.c_str()),Form("%s_gau_lau1",prefix.c_str()),*obs_var,*mean,*sigma);
      setFftObsBinning(obs_var);
      RooFFTConvPdf *gauxlau = new RooFFTConvPdf(Form("%s1",prefix.c_str()),Form("%s_gauxlau1",prefix.c_str()),*obs_var,*soft_step,*gau);
      setCacheObsBinning(obs_var);
      gauxlau->setBufferFraction(kBkgBufferFraction);
      return gauxlau;
  } else if (order==2) {
      RooRealVar *cp1 = new RooRealVar(Form("%s_cp1_lau2",prefix.c_str()),Form("%s_cp1_lau2",prefix.c_str()),coeff1_lau2,coeff1_llau2,coeff1_hlau2);
      RooRealVar *cp2 = new RooRealVar(Form("%s_cp2_lau2",prefix.c_str()),Form("%s_cp2_lau2",prefix.c_str()),coeff2_lau2,coeff2_llau2,coeff2_hlau2);
      fixUnitShapeCoeff(cp1);
      RooGenericPdf *soft_step = new RooGenericPdf(
        Form("%s_soft_lau2",prefix.c_str()),Form("%s_soft_lau2",prefix.c_str()),
        "1e-20+0.5*(1.0+TMath::Erf((@0-@1)/(@2*sqrt(2.))))*(@3*(@0)^(-4)+@4*(@0)^(-5))",
        RooArgList(*obs_var,*turnon,*width,*cp1,*cp2)
      );
      RooGaussModel *gau = new RooGaussModel(Form("%s_gau_lau2",prefix.c_str()),Form("%s_gau_lau2",prefix.c_str()),*obs_var,*mean,*sigma);
      setFftObsBinning(obs_var);
      RooFFTConvPdf *gauxlau = new RooFFTConvPdf(Form("%s2",prefix.c_str()),Form("%s_gauxlau2",prefix.c_str()),*obs_var,*soft_step,*gau);
      setCacheObsBinning(obs_var);
      gauxlau->setBufferFraction(kBkgBufferFraction);
      return gauxlau;
  } else if (order==3) {
      RooRealVar *cp1 = new RooRealVar(Form("%s_cp1_lau3",prefix.c_str()),Form("%s_cp1_lau3",prefix.c_str()),coeff1_lau3,coeff1_llau3,coeff1_hlau3);
      RooRealVar *cp2 = new RooRealVar(Form("%s_cp2_lau3",prefix.c_str()),Form("%s_cp2_lau3",prefix.c_str()),coeff2_lau3,coeff2_llau3,coeff2_hlau3);
      RooRealVar *cp3 = new RooRealVar(Form("%s_cp3_lau3",prefix.c_str()),Form("%s_cp3_lau3",prefix.c_str()),coeff3_lau3,coeff3_llau3,coeff3_hlau3);
      fixUnitShapeCoeff(cp1);
      RooGenericPdf *soft_step = new RooGenericPdf(
        Form("%s_soft_lau3",prefix.c_str()),Form("%s_soft_lau3",prefix.c_str()),
        "1e-20+0.5*(1.0+TMath::Erf((@0-@1)/(@2*sqrt(2.))))*(@3*(@0)^(-4)+@4*(@0)^(-5)+@5*(@0)^(-6))",
        RooArgList(*obs_var,*turnon,*width,*cp1,*cp2,*cp3)
      );
      RooGaussModel *gau = new RooGaussModel(Form("%s_gau_lau3",prefix.c_str()),Form("%s_gau_lau3",prefix.c_str()),*obs_var,*mean,*sigma);
      setFftObsBinning(obs_var);
      RooFFTConvPdf *gauxlau = new RooFFTConvPdf(Form("%s3",prefix.c_str()),Form("%s_gauxlau3",prefix.c_str()),*obs_var,*soft_step,*gau);
      setCacheObsBinning(obs_var);
      gauxlau->setBufferFraction(kBkgBufferFraction);
      return gauxlau;
  } 
  else if (order==4) {
      RooRealVar *cp1 = new RooRealVar(Form("%s_cp1_lau4",prefix.c_str()),Form("%s_cp1_lau4",prefix.c_str()),coeff1_lau4,coeff1_llau4,coeff1_hlau4);
      RooRealVar *cp2 = new RooRealVar(Form("%s_cp2_lau4",prefix.c_str()),Form("%s_cp2_lau4",prefix.c_str()),coeff2_lau4,coeff2_llau4,coeff2_hlau4);
      RooRealVar *cp3 = new RooRealVar(Form("%s_cp3_lau4",prefix.c_str()),Form("%s_cp3_lau4",prefix.c_str()),coeff3_lau4,coeff3_llau4,coeff3_hlau4);
      RooRealVar *cp4 = new RooRealVar(Form("%s_cp4_lau4",prefix.c_str()),Form("%s_cp4_lau4",prefix.c_str()),coeff4_lau4,coeff4_llau4,coeff4_hlau4);
      fixUnitShapeCoeff(cp1);
      RooGenericPdf *soft_step = new RooGenericPdf(
        Form("%s_soft_lau4",prefix.c_str()),Form("%s_soft_lau4",prefix.c_str()),
        "1e-20+0.5*(1.0+TMath::Erf((@0-@1)/(@2*sqrt(2.))))*(@3*(@0)^(-4)+@4*(@0)^(-5)+@5*(@0)^(-6)+@6*(@0)^(-7))",
        RooArgList(*obs_var,*turnon,*width,*cp1,*cp2,*cp3,*cp4)
      );
      RooGaussModel *gau = new RooGaussModel(Form("%s_gau_lau4",prefix.c_str()),Form("%s_gau_lau4",prefix.c_str()),*obs_var,*mean,*sigma);
      setFftObsBinning(obs_var);
      RooFFTConvPdf *gauxlau = new RooFFTConvPdf(Form("%s4",prefix.c_str()),Form("%s_gauxlau4",prefix.c_str()),*obs_var,*soft_step,*gau);
      setCacheObsBinning(obs_var);
      gauxlau->setBufferFraction(kBkgBufferFraction);
      return gauxlau;
  } 
   else {
	return NULL;
  }
}

RooAbsPdf* PdfModelBuilder::getKeysPdf(string prefix){
  if (!keysPdfAttributesSet){
    cerr << "ERROR -- keysPdf attributes not set" << endl;
    exit(1);
  }
  return new RooKeysPdf(prefix.c_str(),prefix.c_str(),*obs_var,*keysPdfData,RooKeysPdf::MirrorBoth,keysPdfRho);
}

RooAbsPdf* PdfModelBuilder::getPdfFromFile(string &prefix){
  vector<string> details;
  split(details,prefix,boost::is_any_of(","));

  string fname = details[2];
  string wsname = details[1];
  string pdfname = details[0];

  TFile *tempFile = TFile::Open(fname.c_str());
  if (!tempFile){
    cerr << "PdfModelBuilder::getPdfFromFile -- file not found " << fname << endl;
    assert(0);
  }
  RooWorkspace *tempWS = (RooWorkspace*)tempFile->Get(wsname.c_str());
  if (!tempWS){
    cerr << "PdfModelBuilder::getPdfFromFile -- workspace not found " << wsname << endl;
    assert(0);
  }
  RooAbsPdf *tempPdf = (RooAbsPdf*)tempWS->pdf(pdfname.c_str());
  if (!tempPdf){
    cerr << "PdfModelBuilder::getPdfFromFile -- pdf not found " << pdfname << endl;
    assert(0);
  }
  prefix = pdfname;
  RooAbsPdf *pdf = (RooAbsPdf*)tempPdf->Clone(prefix.c_str());
  tempFile->Close();
  delete tempFile;
  return pdf;
}

RooAbsPdf* PdfModelBuilder::getExponentialSingle(string prefix, int order){
  
  if (order%2==0){
    cerr << "ERROR -- addExponential -- only odd number of params allowed" << endl;
    return NULL;
  }
  else {
    int nfracs=(order-1)/2;
    int nexps=order-nfracs;
    assert(nfracs==nexps-1);
    RooArgList *fracs = new RooArgList();
    RooArgList *exps = new RooArgList();
    for (int i=1; i<=nfracs; i++){
      string name =  Form("%s_f%d",prefix.c_str(),i);
      params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),0.9-float(i-1)*1./nfracs,0.0001,0.9999)));
      fracs->add(*params[name]);
    }
    for (int i=1; i<=nexps; i++){
      string name =  Form("%s_p%d",prefix.c_str(),i);
      string ename =  Form("%s_e%d",prefix.c_str(),i);
      params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),TMath::Max(-1.,-0.04*(i+1)),-1.,0.)));
      utilities.insert(pair<string,RooAbsPdf*>(ename, new RooExponential(ename.c_str(),ename.c_str(),*obs_var,*params[name])));
      exps->add(*utilities[ename]);
    }
    RooAbsPdf *exp = new RooAddPdf(prefix.c_str(),prefix.c_str(),*exps,*fracs,true);

    RooRealVar *mean = new RooRealVar("mean","mean",0.0) ;
    mean->setConstant(true);
    RooRealVar *sigma = new RooRealVar("sigma","sigma",5,0.2,30.0) ;
    RooGaussian *gaus = new RooGaussian("gaus","gaus",*obs_var,*mean,*sigma) ;
    RooRealVar *step_value = new RooRealVar("step_value", "step value",115.,100.,130.) ;
    RooGenericPdf *step_func = new RooGenericPdf("step_func","step_func","(1e-20+( @0 > @1)) * @2",RooArgSet(*obs_var,*step_value,*exp));
    obs_var->setRange(-200.0,200.0);
    setFftObsBinning(obs_var);
    RooFFTConvPdf *exp_gaus = new RooFFTConvPdf("exp_gaus","exp_gaus", *obs_var, *gaus, *step_func);
    setCacheObsBinning(obs_var);
    exp_gaus->setBufferFraction(kBkgBufferFraction);
    obs_var->setRange(110.0,180.0);
    return exp;

  }
}


void PdfModelBuilder::addBkgPdf(string type, int nParams, string name, bool cache){
 
  if (!obs_var_set){
    cerr << "ERROR -- obs Var has not been set!" << endl;
    exit(1);
  }
  bool found=false;
  for (vector<string>::iterator it=recognisedPdfTypes.begin(); it!=recognisedPdfTypes.end(); it++){
    if (*it==type) found=true;
  }
  if (!found){
    cerr << "Pdf of type " << type << " is not recognised!" << endl;
    exit(1);
  }
  RooAbsPdf *pdf=0;

  if (type=="Bernstein") pdf = getBernstein(name,nParams);
  if (type=="Exponential") pdf = getExponentialSingle(name,nParams);
  if (type=="PowerLaw") pdf = getPowerLawSingle(name,nParams);
  if (type=="Laurent") pdf = getLaurentSeries(name,nParams);
  if (type=="KeysPdf") pdf = getKeysPdf(name);
  if (type=="File") pdf = getPdfFromFile(name);

  // 新增：讓 "BernsteinStepxGau" 類型可被使用
  if (type=="BernsteinStepxGau") pdf = getBernsteinStepxGau(name,nParams,0);

  if (cache) {
    wsCache->import(*pdf);
    RooAbsPdf *cachePdf = wsCache->pdf(pdf->GetName());
    bkgPdfs.insert(pair<string,RooAbsPdf*>(cachePdf->GetName(),cachePdf));
  }
  else {
    bkgPdfs.insert(pair<string,RooAbsPdf*>(pdf->GetName(),pdf));
  }

}

void PdfModelBuilder::setKeysPdfAttributes(RooDataSet *data, double rho){
  keysPdfData = data;
  keysPdfRho = rho;
  keysPdfAttributesSet=true;
}

void PdfModelBuilder::setSignalPdf(RooAbsPdf *pdf, RooRealVar *norm){
  sigPdf=pdf;
  sigNorm=norm;
  signal_set=true;
}

void PdfModelBuilder::setSignalPdfFromMC(RooDataSet *data){
  
  RooDataHist *sigMCBinned = new RooDataHist(Form("roohist_%s",data->GetName()),Form("roohist_%s",data->GetName()),RooArgSet(*obs_var),*data);
  sigPdf = new RooHistPdf(Form("pdf_%s",data->GetName()),Form("pdf_%s",data->GetName()),RooArgSet(*obs_var),*sigMCBinned);
  sigNorm = new RooConstVar(Form("sig_events_%s",data->GetName()),Form("sig_events_%s",data->GetName()),data->sumEntries());
  signal_set=true;
}

void PdfModelBuilder::makeSBPdfs(bool cache){
  
  if (!signal_set){
    cerr << "ERROR - no signal model set!" << endl;
    exit(1);
  }
  if (!signal_modifier_set){
    cerr << "ERROR - no signal modifier set!" << endl;
    exit(1);
  }
 
  if (sigNorm) {
    sigYield = new RooProduct("sig_yield","sig_yield",RooArgSet(*signalModifier,*sigNorm));
  }
  else {
    sigYield = signalModifier;
  }
  bkgYield = new RooRealVar("bkg_yield","bkg_yield",1000.,0.,1.e6);

  for (map<string,RooAbsPdf*>::iterator bkg=bkgPdfs.begin(); bkg!=bkgPdfs.end(); bkg++){
    RooAbsPdf *sbMod = new RooAddPdf(Form("sb_%s",bkg->first.c_str()),Form("sb_%s",bkg->first.c_str()),RooArgList(*(bkg->second),*sigPdf),RooArgList(*bkgYield,*sigYield));
    if (cache) {
      wsCache->import(*sbMod,RecycleConflictNodes());
      RooAbsPdf *cachePdf = (RooAbsPdf*)wsCache->pdf(sbMod->GetName());
      signalModifier = (RooRealVar*)wsCache->var(signalModifier->GetName());
      sbPdfs.insert(pair<string,RooAbsPdf*>(cachePdf->GetName(),cachePdf));
    }
    else {
      sbPdfs.insert(pair<string,RooAbsPdf*>(sbMod->GetName(),sbMod));
    }
  }
}

map<string,RooAbsPdf*> PdfModelBuilder::getBkgPdfs(){
  return bkgPdfs;
}

map<string,RooAbsPdf*> PdfModelBuilder::getSBPdfs(){
  return sbPdfs;
}

RooAbsPdf* PdfModelBuilder::getSigPdf(){
  return sigPdf;
}

void PdfModelBuilder::plotPdfsToData(RooAbsData *data, int binning, string name, bool bkgOnly,string specificPdfName){
  
  TCanvas *canv = new TCanvas();
  bool specPdf=false;
  if (specificPdfName!="") specPdf=true;

  map<string,RooAbsPdf*> pdfSet;
  if (bkgOnly) pdfSet = bkgPdfs;
  else pdfSet = sbPdfs;
  
  for (map<string,RooAbsPdf*>::iterator it=pdfSet.begin(); it!=pdfSet.end(); it++){
    if (specPdf && it->first!=specificPdfName && specificPdfName!="NONE") continue;
    RooPlot *plot = obs_var->frame();
    data->plotOn(plot, Binning(binning), DataError(RooAbsData::SumW2));
    if (specificPdfName!="NONE") {
	 it->second->plotOn(plot);
	 it->second->paramOn(plot,RooFit::Layout(0.34,0.96,0.89),RooFit::Format("NEA",AutoPrecision(1)));
    }	
    plot->Draw();
    canv->Print(Form("%s_%s.pdf",name.c_str(),it->first.c_str()));
    canv->Print(Form("%s_%s.png",name.c_str(),it->first.c_str()));
  }
  delete canv;
}

void PdfModelBuilder::fitToData(RooAbsData *data, bool bkgOnly, bool cache, bool print){
  
  map<string,RooAbsPdf*> pdfSet;
  if (bkgOnly) pdfSet = bkgPdfs;
  else pdfSet = sbPdfs;

  for (map<string,RooAbsPdf*>::iterator it=pdfSet.begin(); it!=pdfSet.end(); it++){
    RooFitResult *fit = (RooFitResult*)it->second->fitTo(*data,Save(true));
    if (print){
      cout << "Fit Res Before: " << endl;
      fit->floatParsInit().Print("v");
      cout << "Fit Res After: " << endl;
      fit->floatParsFinal().Print("v");
    }
    if (cache) {
      RooArgSet *fitargs = (RooArgSet*)it->second->getParameters(*obs_var);
      fitargs->remove(*signalModifier); 
      wsCache->defineSet(Form("%s_params",it->first.c_str()),*fitargs);
      wsCache->defineSet(Form("%s_observs",it->first.c_str()),*obs_var);
      wsCache->saveSnapshot(it->first.c_str(),*fitargs,true);
      if (print) {
        cout << "Cached values: " << endl;
        fitargs->Print("v");
      }
    }
  }
  if (bkgOnly) bkgHasFit=true;
  else sbHasFit=true;
}

void PdfModelBuilder::setSeed(int seed){
  RooRandom::randomGenerator()->SetSeed(seed);
}

RooDataSet* PdfModelBuilder::makeHybridDataset(vector<float> switchOverMasses, vector<RooDataSet*> dataForHybrid){
  
  assert(switchOverMasses.size()==dataForHybrid.size()-1);

  vector<string> cut_strings;
  cut_strings.push_back("cutstring0");
  obs_var->setRange("cutstring0",obs_var->getMin(),switchOverMasses[0]);
  for (unsigned int i=1; i<switchOverMasses.size(); i++){
    cut_strings.push_back(Form("cutstring%d",i));
    obs_var->setRange(Form("cutstring%d",i),switchOverMasses[i-1],switchOverMasses[i]);
  }
  cut_strings.push_back(Form("cutstring%d",int(switchOverMasses.size())));
  obs_var->setRange(Form("cutstring%d",int(switchOverMasses.size())),switchOverMasses[switchOverMasses.size()-1],obs_var->getMax());
  
  obs_var->Print("v");
  assert(cut_strings.size()==dataForHybrid.size());
  
	RooDataSet *data=0;
  for (unsigned int i=0; i<dataForHybrid.size(); i++){
    RooDataSet *cutData = (RooDataSet*)dataForHybrid[i]->reduce(Name("hybridToy"),Title("hybridToy"),CutRange(cut_strings[i].c_str()));
    if (i==0) data=cutData;
    else data->append(*cutData);
  }
  return data;
}

void PdfModelBuilder::throwHybridToy(string postfix, int nEvents, vector<float> switchOverMasses, vector<string> functions, bool bkgOnly, bool binned, bool poisson, bool cache){
  
  assert(switchOverMasses.size()==functions.size()-1);
  toyHybridData.clear();

  throwToy(postfix,nEvents,bkgOnly,false,poisson,cache);

  vector<RooDataSet*> dataForHybrid;
  string hybridName = "hybrid";
  for (vector<string>::iterator func=functions.begin(); func!=functions.end(); func++){
    hybridName += "_"+*func;
    for (map<string,RooDataSet*>::iterator it=toyDataSet.begin(); it!=toyDataSet.end(); it++){
      if (it->first.find(*func)!=string::npos){
        dataForHybrid.push_back(it->second);
      }
    }
  }
  if (dataForHybrid.size()!=functions.size()){
    cerr << "One of the requested hybrid functions has not been found" << endl;
    exit(1);
  }

  RooDataSet *hybridData = makeHybridDataset(switchOverMasses,dataForHybrid);
  toyHybridData.clear();
  if (binned) {
    RooDataHist *hybridDataHist = hybridData->binnedClone();
    hybridDataHist->SetName(Form("%s_%s",hybridName.c_str(),postfix.c_str()));
    toyHybridData.insert(pair<string,RooAbsData*>(hybridDataHist->GetName(),hybridDataHist));
  }
  else {
    hybridData->SetName(Form("%s_%s",hybridName.c_str(),postfix.c_str()));
    toyHybridData.insert(pair<string,RooAbsData*>(hybridData->GetName(),hybridData));
  }
}

void PdfModelBuilder::throwToy(string postfix, int nEvents, bool bkgOnly, bool binned, bool poisson, bool cache){

  toyData.clear();
  toyDataSet.clear();
  toyDataHist.clear();
  map<string,RooAbsPdf*> pdfSet;
  if (bkgOnly) {
    pdfSet = bkgPdfs;
    if (!bkgHasFit) cerr << "WARNING -- bkg has not been fit to data. Are you sure this is wise?" << endl; 
  }
  else {
    pdfSet = sbPdfs;
    if (!sbHasFit) cerr << "WARNING -- sb has not been fit to data. Are you sure this is wise?" << endl;
  }
  
  for (map<string,RooAbsPdf*>::iterator it=pdfSet.begin(); it!=pdfSet.end(); it++){
    if (cache) {
      wsCache->loadSnapshot(it->first.c_str());
      cout << "Loaded snapshot, params at.." << endl;
      it->second->getVariables()->Print("v");
    }
    RooAbsData *toy;
    if (binned){
      RooDataHist *toyHist;
      if (poisson) toyHist = it->second->generateBinned(RooArgSet(*obs_var),nEvents,Extended(),Name(Form("%s_%s",it->first.c_str(),postfix.c_str())));
      else toyHist = it->second->generateBinned(RooArgSet(*obs_var),nEvents,Name(Form("%s_%s",it->first.c_str(),postfix.c_str())));
      toyDataHist.insert(pair<string,RooDataHist*>(toyHist->GetName(),toyHist));
      toy=toyHist;
    }
    else {
      RooDataSet *toySet;
      if (poisson) toySet = it->second->generate(RooArgSet(*obs_var),nEvents,Extended(),Name(Form("%s_%s",it->first.c_str(),postfix.c_str())));
      else toySet = it->second->generate(RooArgSet(*obs_var),nEvents,Name(Form("%s_%s",it->first.c_str(),postfix.c_str())));
      toyDataSet.insert(pair<string,RooDataSet*>(toySet->GetName(),toySet));
      toy=toySet;
    }
    toyData.insert(pair<string,RooAbsData*>(toy->GetName(),toy));
  }
  
}

map<string,RooAbsData*> PdfModelBuilder::getToyData(){
  return toyData;
}

map<string,RooAbsData*> PdfModelBuilder::getHybridToyData(){
  return toyHybridData;
}

void PdfModelBuilder::plotHybridToy(string prefix, int binning, vector<float> switchOverMasses, vector<string> functions, bool bkgOnly){

  map<string,RooAbsPdf*> pdfSet;
  if (bkgOnly) {
    pdfSet = bkgPdfs;
  }
  else {
    pdfSet = sbPdfs;
  }
  
  int tempColors[4] = {kBlue,kRed,kGreen+2,kMagenta};

  vector<string> cut_strings;
  cut_strings.push_back("cutstring0");
  obs_var->setRange("cutstring0",obs_var->getMin(),switchOverMasses[0]);
  for (unsigned int i=1; i<switchOverMasses.size(); i++){
    cut_strings.push_back(Form("cutstring%d",i));
    obs_var->setRange(Form("cutstring%d",i),switchOverMasses[i-1],switchOverMasses[i]);
  }
  cut_strings.push_back(Form("cutstring%d",int(switchOverMasses.size())));
  obs_var->setRange(Form("cutstring%d",int(switchOverMasses.size())),switchOverMasses[switchOverMasses.size()-1],obs_var->getMax());
  
  RooPlot *plot = obs_var->frame();
  TCanvas *canv = new TCanvas();
  int i=0;
  for (vector<string>::iterator func=functions.begin(); func!=functions.end(); func++){
    for (map<string,RooAbsPdf*>::iterator pdfIt = pdfSet.begin(); pdfIt != pdfSet.end(); pdfIt++){
      if (pdfIt->first.find(*func)!=string::npos) {
        for (map<string,RooAbsData*>::iterator toyIt = toyData.begin(); toyIt != toyData.end(); toyIt++){
          if (toyIt->first.find(pdfIt->first)!=string::npos){
            RooAbsData *data = toyIt->second->reduce(CutRange(cut_strings[i].c_str()));
            data->plotOn(plot, Binning(binning), DataError(RooAbsData::SumW2), MarkerColor(tempColors[i]), LineColor(tempColors[i]), CutRange(cut_strings[i].c_str()));
            pdfIt->second->plotOn(plot,LineColor(tempColors[i]),Range(cut_strings[i].c_str()));
            i++;
          }
        }
      }
    }
  }
  for (map<string,RooAbsData*>::iterator hybrid=toyHybridData.begin(); hybrid!=toyHybridData.end(); hybrid++){
    hybrid->second->plotOn(plot, Binning(binning), DataError(RooAbsData::SumW2), MarkerSize(0.8),MarkerStyle(kFullSquare));
    plot->SetMinimum(0.0001);
    plot->Draw();
    canv->Print(Form("%s_%s.pdf",prefix.c_str(),hybrid->first.c_str()));
  }
  delete canv;
}

void PdfModelBuilder::plotToysWithPdfs(string prefix, int binning, bool bkgOnly){
  
  map<string,RooAbsPdf*> pdfSet;
  if (bkgOnly) {
    pdfSet = bkgPdfs;
  }
  else {
    pdfSet = sbPdfs;
  }
  TCanvas *canv = new TCanvas();
  for (map<string,RooAbsPdf*>::iterator pdfIt = pdfSet.begin(); pdfIt != pdfSet.end(); pdfIt++){
    for (map<string,RooAbsData*>::iterator toyIt = toyData.begin(); toyIt != toyData.end(); toyIt++){
      if (toyIt->first.find(pdfIt->first)!=string::npos){
        RooPlot *plot = obs_var->frame();
        toyIt->second->plotOn(plot, Binning(binning), DataError(RooAbsData::SumW2));
        pdfIt->second->plotOn(plot,LineColor(kRed));
        pdfIt->second->paramOn(plot,LineColor(kRed),RooFit::Layout(0.34,0.96,0.89),RooFit::Format("NEA",AutoPrecision(1)));
        plot->Draw();
        canv->Print(Form("%s_%s.pdf",prefix.c_str(),pdfIt->first.c_str()));
        canv->Print(Form("%s_%s.png",prefix.c_str(),pdfIt->first.c_str()));
      }
    }
  }
  delete canv;

}

void PdfModelBuilder::saveWorkspace(string filename){
  TFile *outFile = new TFile(filename.c_str(),"RECREATE");
  outFile->cd();
  wsCache->Write();
  outFile->Close();
  delete outFile;
}

void PdfModelBuilder::saveWorkspace(TFile *file){
  file->cd();
  wsCache->Write();
}
