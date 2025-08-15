#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TColor.h>
#include <TSystem.h>
#include <TString.h>
#include <iostream>
#include <vector>
#include <cmath>

// 如果你已有 CMS_lumi.cpp，请保持这行；没有也可注释掉
#include "src/CMS_lumi.cpp"

using namespace std;

// 需要画的质量点（单位 GeV）
const int massPoints[] = {5, 15, 30};
const int Nmass = sizeof(massPoints)/sizeof(int);

// 读一个 combine 输出 ROOT 文件中的四分位数和可选观测值
// 返回 true 表示至少取到了中位数
bool readLimitsFromFile(const TString& fname,
                        double& q025, double& q16, double& q50, double& q84, double& q975,
                        double& obs,  bool& hasObs)
{
    q025 = q16 = q50 = q84 = q975 = -1.0;
    obs   = -1.0;
    hasObs = false;

    TFile* f = TFile::Open(fname, "READ");
    if (!f || f->IsZombie()) {
        cerr << "Cannot open file: " << fname << endl;
        return false;
    }
    TTree* t = (TTree*) f->Get("limit");
    if (!t) {
        cerr << "No TTree 'limit' in file: " << fname << endl;
        f->Close();
        return false;
    }

    double limit=0, limitErr=0, mh=0;
    float quant=0.f;
    int syst=0, iToy=0, iSeed=0, iChannel=0;
    float t_cpu=0.f, t_real=0.f;

    t->SetBranchAddress("limit", &limit);
    t->SetBranchAddress("limitErr", &limitErr);
    if (t->GetBranch("mh")) t->SetBranchAddress("mh", &mh);
    if (t->GetBranch("syst")) t->SetBranchAddress("syst", &syst);
    if (t->GetBranch("iToy")) t->SetBranchAddress("iToy", &iToy);
    if (t->GetBranch("iSeed")) t->SetBranchAddress("iSeed", &iSeed);
    if (t->GetBranch("iChannel")) t->SetBranchAddress("iChannel", &iChannel);
    if (t->GetBranch("t_cpu")) t->SetBranchAddress("t_cpu", &t_cpu);
    if (t->GetBranch("t_real")) t->SetBranchAddress("t_real", &t_real);
    if (t->GetBranch("quantileExpected")) t->SetBranchAddress("quantileExpected", &quant);

    const double eps = 1e-4;
    Long64_t n = t->GetEntries();
    for (Long64_t i=0; i<n; ++i) {
        t->GetEntry(i);
        if (t->GetBranch("quantileExpected")) {
            if (quant < 0) { hasObs = true; obs = limit; }
            else if (fabs(quant - 0.500) < eps) q50  = limit;
            else if (fabs(quant - 0.160) < eps) q16  = limit;
            else if (fabs(quant - 0.840) < eps) q84  = limit;
            else if (fabs(quant - 0.025) < eps) q025 = limit;
            else if (fabs(quant - 0.975) < eps) q975 = limit;
        } else {
            // 没有 quantileExpected 分支（极少见），则按条目顺序尝试
            // 仅作兜底，不推荐依赖
            if (n == 5) {
                // 0.025, 0.16, 0.5, 0.84, 0.975
                if (i==0) q025 = limit;
                if (i==1) q16  = limit;
                if (i==2) q50  = limit;
                if (i==3) q84  = limit;
                if (i==4) q975 = limit;
            } else if (n >= 6) {
                // 可能包含观测：-1, then 0.025, 0.16, 0.5, 0.84, 0.975
                if (i==0) { hasObs = true; obs = limit; }
                if (i==1) q025 = limit;
                if (i==2) q16  = limit;
                if (i==3) q50  = limit;
                if (i==4) q84  = limit;
                if (i==5) q975 = limit;
            }
        }
    }
    f->Close();
    return (q50 > 0);
}

// 画一个生产机制（sample: -1=pp? 0=gg? 1=qq? ——如果无区分，可忽略）
void BrazilianPlots(int sample=0, bool isInt=true, int year=0, bool APV=false, bool drawObs=false, bool setLimitsOnBR=false)
{
    // 输入文件名模板（当前目录）
    auto makeFileName = [](int m) {
        return TString::Format("/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine/output_combine_results/higgsCombine%d.AsymptoticLimits.mH125.root", m);
    };

    // 构建图形对象
    TGraph *g_exp = new TGraph();
    TGraphAsymmErrors *g_exp_1sigma = new TGraphAsymmErrors();
    TGraphAsymmErrors *g_exp_2sigma = new TGraphAsymmErrors();
    TGraph *g_obs = new TGraph();

    int idx = 0;
    for (int im=0; im<Nmass; ++im) {
        int m = massPoints[im];
        TString fname = makeFileName(m);

        double q025, q16, q50, q84, q975, obs;
        bool hasObs=false;
        if (!readLimitsFromFile(fname, q025, q16, q50, q84, q975, obs, hasObs)) {
            cerr << "Skip mass " << m << " GeV due to missing median." << endl;
            continue;
        }

        // 假设乘以某个 xs（若不需要，设为 1.0 pb）
        // 1.0 pb = 1000.0 fb
        double assume_xs = 1000.0;

        // https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHWG136TeVxsec_extrap?utm_source=chatgpt.com
        // mH=125.09, 52.17 pb, 52170 fb
        double ggF_xs = 52170; // ggF 的假设交叉截面（单位 pb）

        double exp  = q50  * assume_xs;
        double p1s  = (q84  - q50)  * assume_xs;
        double m1s  = (q50  - q16)  * assume_xs;
        double p2s  = (q975 - q50)  * assume_xs;
        double m2s  = (q50  - q025) * assume_xs;

        if (setLimitsOnBR==false)
        {
            exp  = q50  * assume_xs;
            p1s  = (q84  - q50)  * assume_xs;
            m1s  = (q50  - q16)  * assume_xs;
            p2s  = (q975 - q50)  * assume_xs;
            m2s  = (q50  - q025) * assume_xs;
        }

        if (setLimitsOnBR==true)
        {
            exp  = q50  * assume_xs / ggF_xs;
            p1s  = (q84  - q50)  * assume_xs / ggF_xs;
            m1s  = (q50  - q16)  * assume_xs / ggF_xs;
            p2s  = (q975 - q50)  * assume_xs / ggF_xs;
            m2s  = (q50  - q025) * assume_xs / ggF_xs;
        }

        g_exp->SetPoint(idx, m, exp);
        g_exp_1sigma->SetPoint(idx, m, exp);
        g_exp_2sigma->SetPoint(idx, m, exp);
        g_exp_1sigma->SetPointError(idx, 0, 0, m1s, p1s);
        g_exp_2sigma->SetPointError(idx, 0, 0, m2s, p2s);

        if (drawObs && hasObs) g_obs->SetPoint(idx, m, obs * assume_xs);

        ++idx;
    }

    // 画布与坐标框
    double lumi_fb = 62.5; // 可按需改
    TCanvas *c = new TCanvas("cLimits","",800,600);
    c->SetBottomMargin(0.12);
    c->SetRightMargin(0.05);
    c->SetLeftMargin(0.14);

    // x 轴范围用质量点的最小/最大稍微放宽
    double xmin = massPoints[0] - 1 ;
    double xmax = massPoints[Nmass-1] + 1;

    TH1F *frame = new TH1F("frame",
        Form(";m_{a} (GeV);Upper limit on #sigma(H #rightarrow Za #rightarrow 2l + 2#gamma) [fb] at 95%% CL"),
        100, xmin, xmax);
    frame->SetStats(0);
    if (setLimitsOnBR==false)
    {
        frame->GetYaxis()->SetTitle("Upper limit on #sigma(H #rightarrow Za #rightarrow 2l + 2#gamma) [fb] at 95% CL");
        frame->SetMaximum(100);    // 依据你的数值可调整
        frame->SetMinimum(8e-1);   // 依据你的数值可调整
    }

    if (setLimitsOnBR==true)
    {
        frame->GetYaxis()->SetTitle("Upper limit on Br(H #rightarrow Za #rightarrow 2l + 2#gamma) at 95% CL");
        frame->SetMaximum(3e-2);    // 依据你的数值可调整
        frame->SetMinimum(1e-6); // 依据你的数值可调整
    }

    frame->GetXaxis()->SetTitleSize(0.055);
    frame->GetXaxis()->SetLabelSize(0.05);
    frame->GetYaxis()->SetTitleSize(0.041);
    frame->GetYaxis()->SetLabelSize(0.05);
    frame->GetYaxis()->SetTitleOffset(1.5);
    frame->GetYaxis()->CenterTitle(true);
    // frame->GetXaxis()->SetMoreLogLabels();
    frame->GetXaxis()->SetNoExponent();

    // 样式（你也可以改为 kGreen/kYellow 等传统色）
    g_exp->SetMarkerStyle(24);
    g_exp->SetMarkerColor(kBlack);
    g_exp->SetMarkerSize(0.8);
    g_exp->SetLineColor(kBlack);
    g_exp->SetLineWidth(3);
    g_exp->SetLineStyle(2);
    g_exp->SetFillColor(0);

    g_exp_1sigma->SetMarkerStyle(0);
    g_exp_1sigma->SetMarkerColor(3);
    g_exp_1sigma->SetFillColor(kGreen);
    g_exp_1sigma->SetLineColor(kGreen);
    g_exp_1sigma->SetFillStyle(1001);

    g_exp_2sigma->SetMarkerStyle(0);
    g_exp_2sigma->SetMarkerColor(5);
    g_exp_2sigma->SetFillColor(kYellow);
    g_exp_2sigma->SetLineColor(kYellow);
    g_exp_2sigma->SetFillStyle(1001);

    if (drawObs) {
        g_obs->SetMarkerStyle(8);
        g_obs->SetMarkerColor(kBlack);
        g_obs->SetLineColor(kBlack);
        g_obs->SetLineWidth(3);
    }

    // 图例
    double x1, y1, x2, y2;

    if (setLimitsOnBR==false) {
        x1 = 0.59; y1 = 0.68;
        x2 = 0.97; y2 = 0.87;
    } else if (setLimitsOnBR==true) {
        x1 = 0.58; y1 = 0.65;
        x2 = 0.96; y2 = 0.87;
    }

    TLegend *leg = new TLegend(x1,y1,x2,y2);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetFillColor(0);
    leg->SetTextSize(0.045);
    leg->AddEntry(g_exp, "Median expected", "l");
    leg->AddEntry(g_exp_1sigma, "68% expected", "f");
    leg->AddEntry(g_exp_2sigma, "95% expected", "f");
    // if (drawObs) leg->AddEntry(g_obs, "Observed", "lp");

    // 绘制
    frame->Draw();
    g_exp_2sigma->Draw("3 same");
    g_exp_1sigma->Draw("3 same");
    g_exp->Draw("L same");
    if (drawObs) g_obs->Draw("LP same");
    leg->Draw("same");

    // CMS 样式（需要 src/CMS_lumi.cpp）
    CMS_lumi *lumi = new CMS_lumi;
    lumi->set_lumi(c, lumi_fb, 0);

    // c->SetLogx();
    c->SetLogy();
    // c->SetGrid();
    gPad->SetTicks(1,1);
    gPad->RedrawAxis();

    // 输出
    gSystem->mkdir("output_plots", true);
    if (setLimitsOnBR==false)
    {
        c->SaveAs("output_plots/Limits_XS.png");
        c->SaveAs("output_plots/Limits_XS.pdf");    
    }

    if (setLimitsOnBR==true)
    {
        c->SaveAs("output_plots/Limits_BR.png");
        c->SaveAs("output_plots/Limits_BR.pdf");    
    }
}

// 对三个“生产机制”各画一张（如果不区分，可以只保留一个）
void LimitPlots(int year=2018, bool APV=false)
{
    BrazilianPlots(/*sample=*/-1, /*isInt=*/true, year, APV, /*drawObs=*/true, /*setLimitsOnBR=*/false);
    BrazilianPlots(/*sample=*/-1, /*isInt=*/true, year, APV, /*drawObs=*/true, /*setLimitsOnBR=*/true);
    // BrazilianPlots(W, /*sample=*/0, /*isInt=*/true, year, APV, /*drawObs=*/true);
    // BrazilianPlots(W, /*sample=*/1, /*isInt=*/true, year, APV, /*drawObs=*/true);
    // if (W == "NWA") {
    //     BrazilianPlots(W, /*sample=*/-1, /*isInt=*/true, year, APV, /*drawObs=*/true);
    // }
}
