#include "src/CMS_lumi.cpp"

using namespace std;

const TString cats[]={"ggH"};

const int massPoints[] = {5,15,30};
const int Nmass = 3;

std::vector<double> extractNumbers(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<double> numbers;

    if (file.is_open()) {
        std::string line;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            std::string word;
            while (iss >> word) {
                if (word == "r") {
                    double num;
                    if (iss >> word && word == "<" && (iss >> num)) {
                        numbers.push_back(num);
                    }
                }
            }
        }
        file.close();
    } else {
        // std::cout << "Unable to open file" << endl;
    }

    return numbers;
}

void BrazilianPlots(TString W, int sample, bool isInt=true, int year = 0, bool APV = 0, bool isObs=false) {

    double _lumi=138;

    TGraph *g_exp=new TGraph();
    TGraphAsymmErrors *g_exp_1sigma=new TGraphAsymmErrors();
    TGraphAsymmErrors *g_exp_2sigma=new TGraphAsymmErrors();
    // TGraph *g_obs=new TGraph();

    int idx=0;
    for (int im=0;im<Nmass;im++) {
        vector<double> limit;
        TString width;
        if (W=="NWA") width="W0";
        else width=W;

        double exp, p2s, p1s, m1s, m2s, obs;

        limit=extractNumbers(Form("/your/input/log/file_M%d.log",massPoints[im]));

        if (limit.size()!=5) {
            cout<<"Something is wrong!!!"<<endl;
            continue;
        }
        
        float assume_xs = 1.;
        // obs=limit.at(0)*0.1;
        exp=limit.at(2)*assume_xs;
        p2s=limit.at(4)*assume_xs-exp;
        p1s=limit.at(3)*assume_xs-exp;
        m1s=exp-limit.at(1)*assume_xs;
        m2s=exp-limit.at(0)*assume_xs;

        // g_obs->SetPoint(idx,massPointsNew[im],obs);

        g_exp->SetPoint(idx,massPointsNew[im],exp);
        g_exp_1sigma->SetPoint(idx,massPointsNew[im],exp);
        g_exp_2sigma->SetPoint(idx,massPointsNew[im],exp);
        g_exp_1sigma->SetPointError(idx,0,0,m1s,p1s);
        g_exp_2sigma->SetPointError(idx,0,0,m2s,p2s);
        
        idx++;
    }

    CMS_lumi *lumi = new CMS_lumi;
    TCanvas *c=new TCanvas("","",0,0,1000,700);
    c->SetBottomMargin(0.12);
    c->SetRightMargin(0.05);

    TString prod[3]={"pp","gg","qq"};
    TH1F *frame=new TH1F("frame",Form(";m_{X} (GeV);Upper limit on #sigma(%s#rightarrow X#rightarrow ZZ) [pb] at 95%% CL",prod[sample+1].Data()),100,130,3001);
    frame->SetStats(0);
    frame->SetMaximum(0.5);
    frame->SetMinimum(0.0009);
    frame->GetXaxis()->SetTitleSize(0.05);
    frame->GetXaxis()->SetLabelSize(0.045);
    frame->GetYaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetTitleOffset(1);
    frame->GetYaxis()->SetLabelSize(0.045);
    frame->GetXaxis()->SetMoreLogLabels();
    frame->GetXaxis()->SetNoExponent();

    g_exp->SetMarkerStyle(24);
    g_exp->SetMarkerColor(TColor::GetColor("#F0240Bff"));
    g_exp->SetMarkerSize(0.8);
    g_exp->SetLineColor(TColor::GetColor("#F0240Bff"));
    g_exp->SetLineWidth(3);
    g_exp->SetLineStyle(2);
    g_exp->SetFillColor(0);
    
    g_exp_1sigma->SetMarkerStyle(0);
    g_exp_1sigma->SetMarkerColor(3);
    g_exp_1sigma->SetFillColor(TColor::GetColor("#FFDF7Fff"));
    g_exp_1sigma->SetLineColor(TColor::GetColor("#FFDF7Fff"));
    g_exp_1sigma->SetFillStyle(1001);

    g_exp_2sigma->SetMarkerStyle(0);
    g_exp_2sigma->SetMarkerColor(5);
    g_exp_2sigma->SetFillColor(TColor::GetColor("#85D1FBff"));
    g_exp_2sigma->SetLineColor(TColor::GetColor("#85D1FBff"));
    g_exp_2sigma->SetFillStyle(1001);

    // if (isObs) {
    //     g_obs->SetMarkerStyle(8);
    //     g_obs->SetMarkerColor(kBlack);
    //     // g_obs->SetMarkerSize(1);
    //     g_obs->SetLineColor(kBlack);
    //     g_obs->SetLineWidth(3);
    //     // g_obs->SetLineStyle(2);
    //     g_obs->SetFillColor(0);
    // }


    TString prod_[3]={"f_{VBF} floating","ggF","VBF"};
    TLegend *l=new TLegend(0.53,0.57,0.92,0.87,Form("%s, %s",prod_[sample+1].Data(),Wmap[W].Contains(" ")?("#Gamma_{X}="+Wmap[W]).Data():Wmap[W].Data()));

    l->SetFillColor(0);
    l->SetTextSize(0.043);
    l->AddEntry(g_exp,"Median expected","l");
    l->AddEntry(g_exp_1sigma,"68% expected","f");
    l->AddEntry(g_exp_2sigma,"95% expected","f");
    if (isObs) l->AddEntry(g_obs,"Observed","lp");

    frame->Draw();
    g_exp_2sigma->Draw("3same");
    g_exp_1sigma->Draw("3same");
    g_exp->Draw("Lsame");
    if (isObs) g_obs->Draw("Lpsame");
    l->Draw("same");
    lumi->set_lumi(c, _lumi, 0);

    TLatex *text = new TLatex(130, 0.000685, "130");
    text->SetTextFont(42);
    text->SetTextAlign(22);  // Center align the text
    text->SetTextSize(0.045);
    text->Draw();

    c->SetLogy();
    c->SetLogx();
    c->SetGrid();
    gPad->SetTicks(1,1);
    gPad->RedrawAxis();
    c->SaveAs("/your/output/file.png");
    c->SaveAs("/your/output/file.pdf","pdf");
}

void LimitPlots(TString W, int year = 2018, bool APV = 0) {

    BrazilianPlots(W, 0, 1, year, APV,1);
    BrazilianPlots(W, 1, 1, year, APV,1);
    if (W=="NWA") {
        BrazilianPlots(W, -1, 1, year, APV,1);
    }
}
