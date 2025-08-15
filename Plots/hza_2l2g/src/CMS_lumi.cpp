// Include classes
#include "../include/CMS_lumi.h"
#include <iostream>


// Constructor
//====================
CMS_lumi::CMS_lumi(){}
//====================


// Destructor
//=====================
CMS_lumi::~CMS_lumi(){}
//=====================


void CMS_lumi::set_lumi( TPad* pad, float lumi, int iPosX , TString text, bool drawLumi)
{           
   bool outOfFrame = false;
  
   if( iPosX / 10 == 0 )
   {
      outOfFrame = true;
   }
   
   int alignY_ = 3;
   int alignX_ = 2;
  
   if( iPosX / 10 == 0 )
   {
      alignX_ = 1;
   }
   
   if( iPosX == 0 )
   {
      alignX_ = 1;
      alignY_ = 1;     
   }
   
   if( iPosX/10 == 1 )
   {
      alignX_ = 1;
   }
   
    if( iPosX/10 == 2 )
   {
      alignX_ = 2;
   }
   
    if( iPosX/10 == 3 )
   {
      alignX_ = 3;
   }
   
   if( iPosX == 0  ) relPosX = 0.14;
  
   int align_ = 10*alignX_ + alignY_;

//   float H = pad->GetWh();
//   float W = pad->GetWw();
   float l = pad->GetLeftMargin();
   float t = pad->GetTopMargin();
   float r = pad->GetRightMargin();
   float b = pad->GetBottomMargin();
   //float e = 0.025;

   pad->cd();
   
   if (lumi>100)
      lumiText = Form("%.0f",lumi) + lumi_sqrtS;
   else
      lumiText = Form("%.1f",lumi) + lumi_sqrtS;
   
   
   TLatex latex;
   latex.SetNDC();
   latex.SetTextAngle(0);
   latex.SetTextColor(kBlack);    

   float extraTextSize = extraOverCmsTextSize*cmsTextSize;

   latex.SetTextFont(42);
   latex.SetTextAlign(31); 
   latex.SetTextSize(lumiTextSize*t);    
   if (drawLumi) {
      latex.DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText);
   }
   else {
      latex.DrawLatex(1-r,1-t+lumiTextOffset*t,"13.6 TeV");
   }

   if( outOfFrame )
   {
      latex.SetTextFont(cmsTextFont);
      latex.SetTextAlign(11); 
      latex.SetTextSize(cmsTextSize*t);    
      latex.DrawLatex(l,1-t+lumiTextOffset*t,cmsText);
   }
  
   pad->cd();

   float posX_ = 0;
   if( iPosX % 10 <= 1 )
   {
      posX_ = l + relPosX*(1-l-r);
   }
   else if( iPosX % 10 == 2 )
   {
      posX_ =  l+0.5*(1-l-r);
   }
   else if( iPosX % 10 == 3 )
   {
      posX_ =  1-r-relPosX*(1-l-r);
   }
  
   float posY_ = 1-t-relPosY*(1-t-b);
      
   if( !outOfFrame )
   {
      if( drawLogo )
      {
         /*posX_ =   l + 0.045*(1-l-r)*W/H;
         posY_ = 1-t - 0.045*(1-t-b);
         float xl_0 = posX_;
         float yl_0 = posY_ - 0.15;
         float xl_1 = posX_ + 0.15*H/W;
         float yl_1 = posY_;
         TImage* CMS_logo = new TImage("CMS-BW-label.png");
         TPad* pad_logo = new TPad("logo","logo", xl_0, yl_0, xl_1, yl_1 );
         pad_logo->Draw();
         pad_logo->cd();
         CMS_logo->Draw("X");
         pad_logo->Modified();
         pad->cd();*/
      }
      else
      {
         latex.SetTextFont(cmsTextFont);
         latex.SetTextSize(cmsTextSize*t);
         latex.SetTextAlign(align_);
         latex.DrawLatex(posX_, posY_, cmsText);
         
         if( writeExtraText ) 
         {
            latex.SetTextFont(extraTextFont);
            latex.SetTextAlign(align_);
            latex.SetTextSize(extraTextSize*t);
            const float extraDX = -0.02; 
            latex.DrawLatex(posX_ + extraDX*(1 - l - r),  // X 方向偏移
                            posY_ - relExtraDY*cmsTextSize*t,
                            text);

         }
      }
   }
   else if( writeExtraText )
   {
      if( iPosX == 0 ) 
      {
         posX_ = l+relPosX*(1-l-r);
         posY_ = 1-t+lumiTextOffset*t;
      }
      
      latex.SetTextFont(extraTextFont);
      latex.SetTextSize(extraTextSize*t);
      latex.SetTextAlign(align_);
      // 负数 = 往左挪；范围 -0.015 ~ -0.04 先试 -0.02
      const float extraDX = -0.03;

      latex.DrawLatex(posX_ + extraDX*(1 - l - r),  // 这里也加上水平偏移
                     posY_,
                     text);   }
   return;
}
