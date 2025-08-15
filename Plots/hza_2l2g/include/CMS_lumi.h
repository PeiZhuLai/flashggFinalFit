#ifndef CMS_lumi_h
#define CMS_lumi_h

#include "TPad.h"
#include "TLatex.h"
#include "TLine.h"
#include "TBox.h"

using namespace std;

class CMS_lumi
{

public:
      
   CMS_lumi();
   ~CMS_lumi();
   void set_lumi(TPad* pad, float lumi, int iPosX = 10 , TString text="Preliminary", bool drawLumi=true);
      
private:
      
   TString lumiText;
      
   TString cmsText   = "CMS";
   float cmsTextFont = 61; // default is helvetic-bold

   bool writeExtraText = true;
   TString extraText   = "Simulation Preliminary";
   float extraTextFont = 52;  // default is helvetica-italics

   float lumiTextSize   = 0.5; //0.5
   float lumiTextOffset = 0.07;
   float cmsTextSize    = 0.55; //0.5
   float cmsTextOffset  = 0.1;  // only used in outOfFrame version
   
   float relPosX    = 0.045;
   float relPosY    = 0.035;
   float relExtraDY = 1.2;
      
   // ratio of "CMS" and extra text size
   float extraOverCmsTextSize  = 0.95; //0.76
   
   TString lumi_sqrtS = " fb^{-1} (13.6 TeV)";
      
   bool drawLogo = false;    
};

#endif
