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
#include "HiggsAnalysis/CombinedLimit/interface/RooMultiPdf.h"

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
// 新增：集中式樣式設定
#include "../interface/PlotStyle.h"
#include "TStyle.h" // 新增：使用 gStyle 設定軸樣式
#include "TColor.h" // 新增：確保可設定 PS 顏色模型
// 新增：判斷字元是否為數字與取字串長度
#include <cctype>
#include <cstring>
#include <cmath> // 新增：std::isfinite

using namespace std;
using namespace RooFit;
using namespace boost;

namespace po = program_options;

// 新增：前置宣告，避免使用在定義之前
static void transferMacros(TFile* in, TFile* out);
static void ensureOrderSuffix(RooAbsPdf* pdf, int order);
static int getBestFitFunction(RooMultiPdf *mpdf, RooAbsData *data, RooCategory *catIndex, bool silent);

bool BLIND = true;
bool runFtestCheckWithToys=false;
// 新增：只畫圖不做任何 fit 的快速模式
bool PLOT_ONLY = false;

float mgglow_ =110.;//FIXME
float mgghigh_ =180;//FIXME
float mggblindlow_ =115;//FIXME
float mggblindhigh_ =135;//FIXME

float mgg_low =110.;//FIXME
float mgg_high =180.;//FIXME
float nBinsForMass = 1.*(mgg_high-mgg_low);
float mgg_blind_low =115;//FIXME
float mgg_blind_high =135;//FIXME

RooRealVar *intLumi_ = new RooRealVar("IntLumi","hacked int lumi", 1000.);

TRandom3 *RandomGen = new TRandom3();

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

  int tries = 0;
  int status = 1;
  double bestNll = 1e12;

  {
    RooArgSet* params = pdf->getParameters(*data);
    std::unique_ptr<TIterator> it(params->createIterator());
    for (RooRealVar* v = (RooRealVar*)it->Next(); v; v = (RooRealVar*)it->Next()) {
      if (v->isConstant()) continue;
      if (v->getError()<=0) {
        double step = (v->hasMax() && v->hasMin()) ? 0.1*fabs(v->getMax()-v->getMin()) : std::max(1e-2, 0.1*fabs(v->getVal()));
        v->setError(step);
      }
      if (v->hasMin() && fabs(v->getVal()-v->getMin())<1e-6) v->setVal(v->getMin()+1e-3*(v->hasMax()? (v->getMax()-v->getMin()) : 1.0));
      if (v->hasMax() && fabs(v->getVal()-v->getMax())<1e-6) v->setVal(v->getMax()-1e-3*(v->hasMin()? (v->getMax()-v->getMin()) : 1.0));
    }
  }

  while (status!=0 && tries<MaxTries) {
    std::unique_ptr<RooAbsReal> nll(pdf->createNLL(
      *data,
      RooFit::Offset(true),
      RooFit::Optimize(true),
      RooFit::SumW2Error(kFALSE)
    ));
    RooMinimizer minim(*nll);
    minim.setPrintLevel(-1);
    minim.setStrategy(1);
    minim.setOffsetting(true);
    minim.optimizeConst(1);
    minim.setEps(1e-6);

    status = minim.minimize("Minuit2","minimize");

    if (status!=0) {
      minim.setStrategy(2);
      status = minim.minimize("Minuit2","minimize");
    }
    if (status!=0) {
      status = minim.minimize("Minuit2","simplex");
      if (status==0) {
        status = minim.minimize("Minuit2","minimize");
      }
    }

    std::unique_ptr<RooFitResult> res(minim.save());
    double nllVal = nll->getVal();
    if (nllVal < bestNll) bestNll = nllVal;

    if (status!=0 && res) {
      RooArgSet* pars = pdf->getParameters((const RooArgSet*)nullptr);
      pars->assignValueOnly(res->randomizePars());
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

  if (!runFtestCheckWithToys) return prob_asym;

  int ndata = data->sumEntries();

  RooFitResult *fitNullData = pdfNull->fitTo(*data,RooFit::Save(1),RooFit::Strategy(1)
		,RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kFALSE),RooFit::Hesse(kFALSE),RooFit::PrintLevel(-1));
  RooFitResult *fitTestData = pdfTest->fitTo(*data,RooFit::Save(1),RooFit::Strategy(1)
		,RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kFALSE),RooFit::Hesse(kFALSE),RooFit::PrintLevel(-1));

  RooArgSet *params_null = pdfNull->getParameters((const RooArgSet*)(0));
  RooArgSet preParams_null;
  params_null->snapshot(preParams_null);
  RooArgSet *params_test = pdfTest->getParameters((const RooArgSet*)(0));
  RooArgSet preParams_test;
  params_test->snapshot(preParams_test);

  int ntoys =5000;
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

	int stat_n=1;
        int stat_t=1;
	int ntries = 0;
	double nllNull,nllTest;
	int MaxTries = 2;
	while (stat_n!=0){
	  if (ntries>=MaxTries) break;
	  RooFitResult *fitNull = pdfNull->fitTo(*binnedtoy,RooFit::Save(1),RooFit::Strategy(1),RooFit::SumW2Error(kFALSE)
		,RooFit::Minimizer("Minuit2","minimize"),RooFit::Minos(0),RooFit::Hesse(0),RooFit::PrintLevel(-1));

	  nllNull = fitNull->minNll();
          stat_n = fitNull->status();
	  if (stat_n!=0) params_null->assignValueOnly(fitNullData->randomizePars());
	  ntries++;
	}

	ntries = 0;
	while (stat_t!=0){
	  if (ntries>=MaxTries) break;
	  RooFitResult *fitTest = pdfTest->fitTo(*binnedtoy,RooFit::Save(1),RooFit::Strategy(1),RooFit::SumW2Error(kFALSE)
		,RooFit::Minimizer("Minuit2","minimize"),RooFit::Minos(0),RooFit::Hesse(0),RooFit::PrintLevel(-1));
	  nllTest = fitTest->minNll();
          stat_t = fitTest->status();
	  if (stat_t!=0) params_test->assignValueOnly(fitTestData->randomizePars());
	  ntries++;
	}

	toyhistStatN.Fill(stat_n);
	toyhistStatT.Fill(stat_t);

  if (stat_t !=0 || stat_n !=0) continue;
	nsuccesst++;
	double chi2_t = 2*(nllNull-nllTest);
	if (chi2_t >= chi2) npass++;
        toyhist.Fill(chi2_t);
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

  double prob;
  int ntoys = 500;
  name+="_gofTest.pdf";
  RooRealVar norm("norm","norm",data->sumEntries(),0,10E6);

  RooExtendPdf *pdf = new RooExtendPdf("ext","ext",*mpdf,norm);

  RooPlot *plot_chi2 = mass->frame();
  data->plotOn(plot_chi2,Binning(nBinsForMass),RooFit::DataError(RooAbsData::SumW2));
  pdf->plotOn(plot_chi2,Name("pdf"));
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
      pdf->fitTo(*binnedtoy,RooFit::Minimizer("Minuit2","minimize"),RooFit::Minos(0),RooFit::Hesse(0),RooFit::PrintLevel(-1),RooFit::Strategy(0),RooFit::SumW2Error(kFALSE));

      RooPlot *plot_t = mass->frame();
      binnedtoy->plotOn(plot_t,RooFit::DataError(RooAbsData::SumW2));
      pdf->plotOn(plot_t);

      double chi2_t = plot_t->chiSquare(np);
      if( chi2_t>=chi2) npass++;
      toy_chi2.push_back(chi2_t*(nBinsForMass-np));
      delete plot_t;
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
  data->plotOn(plot_chi2,Binning(nBinsForMass),RooFit::DataError(RooAbsData::SumW2));
  pdf->plotOn(plot_chi2);

  int np = pdf->getParameters(*data)->getSize()+1;
  double chi2 = plot_chi2->chiSquare(np);

  *prob = getGoodnessOfFit(mass,pdf,data,name);
  RooPlot *plot = mass->frame();
  mass->setRange("unblindReg_1",mgg_low,mgg_blind_low);
  mass->setRange("unblindReg_2",mgg_blind_high,mgg_high);

  plot->GetXaxis()->SetTitle("m_{ll#gamma#gamma} (GeV)");
  plot->GetXaxis()->SetTitleSize(0.05);
  plot->GetXaxis()->SetLabelSize(0.04);
  plot->GetYaxis()->SetTitle("Events / 1 GeV");
  plot->GetYaxis()->SetTitleSize(0.05);
  plot->GetYaxis()->SetLabelSize(0.04);

  if (BLIND) {
    data->plotOn(plot,Binning(mgg_high-mgg_low),CutRange("unblindReg_1"),RooFit::DataError(RooAbsData::SumW2));
    data->plotOn(plot,Binning(mgg_high-mgg_low),CutRange("unblindReg_2"),RooFit::DataError(RooAbsData::SumW2));
  }
  else data->plotOn(plot,Binning(nBinsForMass),RooFit::DataError(RooAbsData::SumW2));

  TCanvas *canv = new TCanvas("","",800,800);
  canv->SetLeftMargin(0.14);
  canv->SetRightMargin(0.05);
  canv->SetTopMargin(0.07);
  canv->SetBottomMargin(0.12);

  pdf->plotOn(plot);
  pdf->paramOn(plot,RooFit::Layout(0.14,0.96,0.89),RooFit::Format("NEA",AutoPrecision(1)));
  if (BLIND) plot->SetMinimum(0.0001);
  plot->SetTitle("");
  plot->GetYaxis()->SetTitleOffset(1.4);  // y 轴标题一般要更远一点，避免和数字重叠
  plot->GetXaxis()->SetTitleOffset(1.2);  // 默认大约 1，可以稍微调大

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

  TCanvas *canv = new TCanvas("","",800,800);
  canv->SetLeftMargin(PlotStyleCfg::canvasLeftMargin);
  canv->SetRightMargin(PlotStyleCfg::canvasRightMargin);
  canv->SetTopMargin(PlotStyleCfg::canvasTopMargin);
  canv->SetBottomMargin(PlotStyleCfg::canvasBottomMargin);
  gStyle->SetOptStat(PlotStyleCfg::showStatBox ? 1 : 0);
  gStyle->SetPadTickX(PlotStyleCfg::axisTickX);
  gStyle->SetPadTickY(PlotStyleCfg::axisTickY);

  // 讓 PDF 向量輸出線更粗，避免消失
  gStyle->SetLineScalePS(4.0);              // 調大：PDF 向量線寬縮放
  gStyle->SetEndErrorSize(0);               // 移除端點過大影響
  // 定義較清楚的虛線樣式（避免 PDF 端虛線太稀疏）
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
  leg->SetTextSize(0.035); // 縮小 legend 文字，避免 PNG 顯得過大
  leg->SetTextFont(42);

  RooPlot *plot = mass->frame();

  mass->setRange("unblindReg_1",mgg_low,mgg_blind_low);
  mass->setRange("unblindReg_2",mgg_blind_high,mgg_high);
  if (BLIND) {
    data->plotOn(plot,Binning(mgg_high-mgg_low),CutRange("unblindReg_1"),RooFit::DataError(RooAbsData::SumW2));
    data->plotOn(plot,Binning(mgg_high-mgg_low),CutRange("unblindReg_2"),RooFit::DataError(RooAbsData::SumW2));
    data->plotOn(plot,Binning(mgg_high-mgg_low),Invisible());
  }
  else data->plotOn(plot,Binning(nBinsForMass),RooFit::DataError(RooAbsData::SumW2));
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
    // 避免過淡顏色在 PDF 幾乎不可見
    if (col==kWhite || col==kYellow) col = kOrange+7;

    if (icat>6) { col=kBlack; style++; }
    catIndex->setIndex(icat);
    if (!PLOT_ONLY) {
      pdfs->getCurrentPdf()->fitTo(*data,RooFit::Minos(0),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kFALSE));
    }

    // 唯一命名 + 加粗線寬，取消填色，避免 PDF 遮蓋/消失
    std::string curveName = Form("bkg_curve_%d",icat);
    pdfs->getCurrentPdf()->plotOn(
      plot,
      RooFit::Binning(nBinsForMass),
      LineColor(col),
      LineStyle(style),
      LineWidth(4),                 // 調粗線寬
      Name(curveName.c_str())
    );

    RooCurve *thisCurve = dynamic_cast<RooCurve*>(plot->findObject(curveName.c_str()));
    if (thisCurve) {
      thisCurve->SetFillStyle(0);   // 取消任何填色
      thisCurve->SetLineWidth(4);   // 再次確保線寬
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

  // 刪除 pad1 上多餘的 hbplottmp 疊圖，避免 PDF 遮蔽
  // TH1D *hbplottmp = nullptr;
  // if (canDoRatio) {
  //   hbplottmp = (TH1D*) pdf->createHistogram("hbplottmp",*mass,Binning(mgg_high-mgg_low,mgg_low,mgg_high));
  //   hbplottmp->SetDirectory(nullptr);
  //   hbplottmp->Scale(plotdata->Integral());
  //   hbplottmp->SetFillStyle(0);
  //   hbplottmp->SetFillColor(0);
  //   hbplottmp->SetLineColor(bestcol);
  //   hbplottmp->SetLineWidth(1);
  //   hbplottmp->Draw("HISTSAME");
  // }

  leg->Draw("same");

  // 確保所有背景曲線最後疊在最上層（避免被其他物件覆蓋）
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
   if ((xtmp > mgg_blind_low ) && ( xtmp < mgg_blind_high) ) continue;
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
  TH1 *hdummy = new TH1D("hdummyweight","",mgg_high-mgg_low,mgg_low,mgg_high);
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
    TLine *line3 = new TLine(mgg_low,0.,mgg_high,0.);
    line3->SetLineColor(bestcol);
    line3->SetLineWidth(PlotStyleCfg::zeroLineWidth);
    line3->Draw();
    hdatasub->Draw("PESAME");
  }

  // 在儲存前明確標記更新，避免 PDF 重繪遺漏物件
  pad1->Modified(); pad1->Update();
  pad2->Modified(); pad2->Update();
  canv->Modified(); canv->Update();

  canv->SaveAs(Form("%s.pdf",name.c_str()));
  canv->SaveAs(Form("%s.png",name.c_str()));
  catIndex->setIndex(currentIndex);
  delete canv;
}

// Truth Plot 
void truth_plot(RooRealVar *mass, map<string,RooAbsPdf*> pdfs, RooAbsData *data, string name, vector<string> flashggCats_, int cat, int bestFitPdf=-1){
  if (!mass) return;
  if (!data) {
    std::cerr << "[WARN] truth_plot(map) called with null data. Skip plotting." << std::endl;
    return;
  }

  // 與 RooMultiPdf 版本相同的畫布與樣式設定
  TCanvas *canv = new TCanvas("","",800,800);
  canv->SetLeftMargin(PlotStyleCfg::canvasLeftMargin);
  canv->SetRightMargin(PlotStyleCfg::canvasRightMargin);
  canv->SetTopMargin(0.08);
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

  // Legend（保留 truths 位置設定）
  TLegend *leg = new TLegend(PlotStyleCfg::truthLegendX1,
                             PlotStyleCfg::truthLegendY1,
                             PlotStyleCfg::truthLegendX2,
                             PlotStyleCfg::truthLegendY2);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  leg->SetTextFont(42);

  RooPlot *plot = mass->frame();

  // 資料與盲區處理，與 RooMultiPdf 版本一致
  mass->setRange("unblindReg_1",mgg_low,mgg_blind_low);
  mass->setRange("unblindReg_2",mgg_blind_high,mgg_high);
  if (BLIND) {
    data->plotOn(plot,Binning(mgg_high-mgg_low),CutRange("unblindReg_1"),RooFit::DataError(RooAbsData::SumW2));
    data->plotOn(plot,Binning(mgg_high-mgg_low),CutRange("unblindReg_2"),RooFit::DataError(RooAbsData::SumW2));
    data->plotOn(plot,Binning(mgg_high-mgg_low),Invisible());
  } else {
    data->plotOn(plot,Binning(mgg_high-mgg_low),RooFit::DataError(RooAbsData::SumW2));
  }

  TObject *datLeg = plot->getObject(int(plot->numItems()-1));
  leg->AddEntry(datLeg,"Data","LEP");

  RooHist *plotdata = (RooHist*)plot->getObject(plot->numItems()-1);

  // 建立上下兩個 pad
  bool doRatioPlot_ = PlotStyleCfg::enableRatio;
  if (name.find("truths_cat0") != std::string::npos) {
    doRatioPlot_ = false;
  }

  if (!doRatioPlot_) {
    canv->cd();
    int i=0, style=1;
    for (auto it = pdfs.begin(); it != pdfs.end(); ++it, ++i) {
      RooAbsPdf* p = it->second;
      if (!p) continue;
      int col = PlotStyleCfg::colorForIndex(i);
      if (col==kWhite || col==kYellow) col = kOrange+7;
      if (i>6) { col=kBlack; style++; }
      p->plotOn(plot, RooFit::Binning(nBinsForMass), LineColor(col), LineStyle(style), LineWidth(4));
      TObject *pdfLeg = plot->getObject(int(plot->numItems()-1));
      leg->AddEntry(pdfLeg, it->first.c_str(), "L");
    }

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

    plot->SetTitle(Form("Category %s",flashggCats_[cat].c_str()));
    if (BLIND) plot->SetMinimum(0.0001);
    plot->Draw();
    CMS_lumi(canv, 22, 0);
    leg->Draw("same");
    canv->RedrawAxis();

    canv->SaveAs(Form("%s.pdf",name.c_str()));
    canv->SaveAs(Form("%s.png",name.c_str()));
    delete canv;
    return;
  }

  TPad *pad1 = new TPad("pad1","pad1",0, doRatioPlot_ ? PlotStyleCfg::pad2Height : 0.0, 1, 1);
  TPad *pad2 = doRatioPlot_ ? new TPad("pad2","pad2",0,0,1,PlotStyleCfg::pad2Height) : nullptr;
  pad1->SetTopMargin(PlotStyleCfg::pad1TopMargin);
  pad1->SetBottomMargin(doRatioPlot_ ? PlotStyleCfg::pad1BottomMargin : PlotStyleCfg::canvasBottomMargin);
  if (doRatioPlot_) {
    pad2->SetTopMargin(PlotStyleCfg::pad2TopMargin);
    pad2->SetBottomMargin(PlotStyleCfg::pad2BottomMargin);
  }
  pad1->SetLeftMargin(PlotStyleCfg::canvasLeftMargin);
  pad1->SetRightMargin(PlotStyleCfg::canvasRightMargin);
  if (doRatioPlot_) {
    pad2->SetLeftMargin(PlotStyleCfg::canvasLeftMargin);
    pad2->SetRightMargin(PlotStyleCfg::canvasRightMargin);
  }
  pad1->SetTicks(1,1);
  if (doRatioPlot_) {
    pad2->SetTicks(1,1);
    pad2->Draw();
  }
  pad1->Draw();
  pad1->cd();

  // 繪製每個 pdf，樣式與 RooMultiPdf 版本一致
  int i=0, style=1;
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
    if (col==kWhite || col==kYellow) col = kOrange+7;
    if (i>6) { col=kBlack; style++; }

    // 唯一命名，線寬、無填色
    std::string curveName = Form("bkg_curve_%d",i);
    p->plotOn(
      plot,
      RooFit::Binning(nBinsForMass),
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

    std::string ext = "";
    if (bestFitPdf==i) {
      ext=" (Best Fit) ";
      nomBkgCurve = thisCurve ? thisCurve : (RooCurve*)plot->getObject(plot->numItems()-1);
      bestcol = col;
    }
    TObject *pdfLeg = thisCurve ? (TObject*)thisCurve : plot->getObject(int(plot->numItems()-1));
    leg->AddEntry(pdfLeg,Form("%s%s",it->first.c_str(),ext.c_str()),"L");
  }
  bool canDoRatio = (nomBkgCurve!=nullptr);

  // pad1 軸樣式與標題
  plot->GetYaxis()->SetTitle("Events / 1 GeV");
  plot->GetYaxis()->SetTitleFont(PlotStyleCfg::axisTitleFont);
  plot->GetYaxis()->SetLabelFont(PlotStyleCfg::axisLabelFont);
  plot->GetYaxis()->SetTitleSize(PlotStyleCfg::pad1axisTitleSizeY);
  plot->GetYaxis()->SetLabelSize(PlotStyleCfg::pad1axisLabelSizeY);
  plot->GetYaxis()->SetTitleOffset(PlotStyleCfg::pad1axisTitleOffsetY);
  plot->GetXaxis()->SetTitle("m_{ll#gamma#gamma} GeV");
  plot->GetXaxis()->SetLabelSize(doRatioPlot_ ? 0.0 : PlotStyleCfg::pad2axisLabelSizeX);
  if (!doRatioPlot_) {
    plot->GetXaxis()->SetTitleSize(PlotStyleCfg::pad2axisTitleSizeX);
    plot->GetXaxis()->SetTitleOffset(PlotStyleCfg::pad2axisTitleOffsetX);
  }

  plot->SetTitle(Form("Category %s",flashggCats_[cat].c_str()));
  if (BLIND) plot->SetMinimum(0.0001);
  plot->Draw();
  CMS_lumi(canv, 22, 0);

  // 確保背景曲線在最上層
  for (auto* c : pdfCurves) {
    if (!c) continue;
    c->SetFillStyle(0);
    c->SetLineWidth(3);
    c->Draw("L same");
  }

  leg->Draw("same");
  pad1->RedrawAxis();

  // 計算殘差圖 Data - Best Fit（與 RooMultiPdf 版本一致）
  if (doRatioPlot_) {
    int npoints = plotdata->GetN();
    double xtmp, ytmp;
    int point = 0;
    TGraphAsymmErrors *hdatasub = new TGraphAsymmErrors(npoints);
    for (int ipoint=0; ipoint<npoints; ++ipoint) {
      plotdata->GetPoint(ipoint, xtmp, ytmp);
      double bkgval = 0.;
      if (canDoRatio) {
        bkgval = nomBkgCurve->interpolate(xtmp);
      } else {
        continue;
      }
      if (BLIND) {
        if ((xtmp > mgg_blind_low) && (xtmp < mgg_blind_high)) continue;
      }
      double errhi = plotdata->GetErrorYhigh(ipoint);
      double errlow = plotdata->GetErrorYlow(ipoint);
      bool drawZeroBins_ = 1;
      if (!drawZeroBins_) if (fabs(ytmp)<1e-5) continue;
      hdatasub->SetPoint(point, xtmp, ytmp - bkgval);
      hdatasub->SetPointError(point, 0., 0., errlow, errhi);
      point++;
    }

    pad2->cd();
    TH1 *hdummy = new TH1D("hdummyweight","",mgg_high-mgg_low,mgg_low,mgg_high);
    hdummy->SetMaximum(hdatasub->GetHistogram()->GetMaximum()+1);
    hdummy->SetMinimum(hdatasub->GetHistogram()->GetMinimum()-1);

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
      TLine *line3 = new TLine(mgg_low,0.,mgg_high,0.);
      line3->SetLineColor(bestcol);
      line3->SetLineWidth(PlotStyleCfg::zeroLineWidth);
      line3->Draw();
      hdatasub->Draw("PESAME");
    }
  }

  pad1->Modified(); pad1->Update();
  if (doRatioPlot_) {
    pad2->Modified(); pad2->Update();
  }
  canv->Modified(); canv->Update();

  canv->SaveAs(Form("%s.pdf",name.c_str()));
  canv->SaveAs(Form("%s.png",name.c_str()));
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
    RooMsgService::instance().setSilentMode(true);
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
			inWS = (RooWorkspace*)inFile->Get("cms_hgg_workspace");
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

	vector<string> functionClasses;
	functionClasses.push_back("Bernstein");
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
	double upperEnvThreshold = 0.1;

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

		mass->setBins(nBinsForMass);
		RooAbsData *data = dataFull;

		RooArgList storedPdfs("store");

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
          // 確保 legend 名稱包含 order 尾碼
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
        RooRealVar nBackground(Form("CMS_hgg_%s_%s_bkgshape_norm",catname.c_str(),ext.c_str()),"nbkg",data->sumEntries(),0,3*data->sumEntries());
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
			std::vector<int> pdforders;

			int counter =0;
			while (prob<0.05 && order < 7){
				
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


				while (prob<upperEnvThreshold){
					RooAbsPdf *bkgPdf = getPdf(pdfsModel,*funcType,order,"", mass_ALP);

          if (!bkgPdf ){
						if (order >6) { std::cout << " [WARNING] could not add ] " << std::endl; break ;}
						order++;
					}
					else {
						// 在任何擬合與加入 storedPdfs 之前，先修正名稱尾碼
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
            
						if ((prob < upperEnvThreshold) ) {

							if (gofProb > 0.01 || order == truthOrder ) {
								std::cout << "[INFO] Adding to Envelope " << bkgPdf->GetName() << " "<< gofProb
									<< " 2xNLL + c is " << myNll + bkgPdf->getVariables()->getSize() <<  std::endl;

                if (logFile) {
                  int isTruth = (order==truthOrder) ? 1 : 0;
                  fprintf(logFile,"category : %d , pdf : %s , gof : %f, isTruth : %d \n ",cat, bkgPdf->GetName(), gofProb, isTruth);
                }

								storedPdfs.add(*bkgPdf);
								pdforders.push_back(order);
								if ((myNll + bkgPdf->getVariables()->getSize()) < MinimimNLLSoFar) {
									simplebestFitPdfIndex = storedPdfs.getSize()-1;
									MinimimNLLSoFar = myNll + bkgPdf->getVariables()->getSize();
								}
							}
						}

						prev_order=order;
						prev_pdf=bkgPdf;
						order++;
					}
				}

				fprintf(resFile,"%15s & %d & %5.2f & %5.2f \\\\\n",funcType->c_str(),cache_order+1,chi2,prob);
				choices_envelope.insert(pair<string,std::vector<int> >(*funcType,pdforders));
			}
		}

		fprintf(resFile,"\\hline\n");
		choices_vec.push_back(choices);
		choices_envelope_vec.push_back(choices_envelope);
		pdfs_vec.push_back(pdfs);

		truth_plot(mass,pdfs,data,Form("%s/truths_cat%d",outDir.c_str(),cat),flashggCats_,cat);

		if (saveMultiPdf){

      if (storedPdfs.getSize() == 0) {
        std::cerr << "[WARN] storedPdfs is empty for " << catname << ". Skip MultiPdf for this category." << std::endl;
        continue;
      }

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
			RooRealVar nBackground(Form("CMS_hgg_%s_%s_bkgshape_norm",catname.c_str(),ext.c_str()),"nbkg",data->sumEntries(),0,3*data->sumEntries());
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
// 新增：將輸入 ROOT 檔中的 TMacro 複製到輸出檔
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
// 新增：若 RooAbsPdf 名稱沒有以數字結尾，補上 order 尾碼
static void ensureOrderSuffix(RooAbsPdf* pdf, int order) {
  if (!pdf) return;
  std::string name(pdf->GetName() ? pdf->GetName() : "");
  if (!name.empty() && std::isdigit(name.back())) return;
  pdf->SetName(Form("%s%d", name.c_str(), order));
}
// 新增：對 RooMultiPdf 逐一嘗試，選擇 2*NLL + #params 最小者為最佳
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