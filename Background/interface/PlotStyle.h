#ifndef FLASHGGFINALFIT_BACKGROUND_INTERFACE_PLOTSTYLE_H
#define FLASHGGFINALFIT_BACKGROUND_INTERFACE_PLOTSTYLE_H

#include "TColor.h"

namespace PlotStyleCfg {
  // Legend 位置：multipdf 視圖
  static const double multipdfLegendX1 = 0.63;
  static const double multipdfLegendY1 = 0.50;
  static const double multipdfLegendX2 = 0.95;
  static const double multipdfLegendY2 = 0.83;

  // Legend 位置：truths 視圖
  static const double truthLegendX1 = 0.60;
  static const double truthLegendY1 = 0.65;
  static const double truthLegendX2 = 0.88;
  static const double truthLegendY2 = 0.88;

  // Ratio pad 幾何與樣式
  static const bool   enableRatio      = true;   // 目前程式仍會建立 pad，這裡先作為統一開關配置用途
  static const double pad2Height       = 0.35;   // 下方 pad 高度 (0~1)
  // 縮小上下 pad 的縫隙：上 pad 底邊界小一點、下 pad 頂邊界給一點空間
  static const double pad1TopMargin    = 0.12;     // 原 0.10
  static const double pad1BottomMargin = 0.01;    // 原 0.01
  static const double pad2TopMargin    = 0.05;    // 原 0.01
  static const double pad2BottomMargin = 0.35;    // 原 0.25
  static const double zeroLineWidth    = 5.0;     // ratio = 0 之參考線線寬

  // 顏色表
  inline int colorForIndex(int i) {
    static int colors[7] = { kBlue, kRed, kMagenta, kGreen+1, kOrange+7, kAzure+10, kBlack };
    return (i >= 0 && i < 7) ? colors[i] : kBlack;
  }

  // 新增：對齊 Signal/scripts/plot_effisigma.py 的畫圖樣式
  // 來源：OFFSET=0.01；axis title/label size、offset、邊界、tick、legend、TLatex、線寬與點大小
  static const double offset = 0.01;

  // Canvas 邊界（c.SetMargin(左, 右, 下, 上)）
  static const double canvasLeftMargin   = 0.12 + offset; // = 0.13
  static const double canvasRightMargin  = 0.035;
  static const double canvasBottomMargin = 0.14;
  static const double canvasTopMargin    = 0.09;

  // 座標軸樣式（放大 title/label）
  static const int    axisTitleFont    = 42;
  static const int    axisLabelFont    = 42;
  static const double axisTitleSize    = 0.065;   // 原 0.055
  static const double axisLabelSize    = 0.055;   // 原 0.05
  static const double axisTitleOffsetX = 1.05;    // 原 1.15
  static const double axisTitleOffsetY = 1.25;    // 原 ~1.2

  static const double pad1axisTitleSizeY    = 0.075;
  static const double pad1axisLabelSizeY    = 0.065;  
  static const double pad1axisTitleOffsetY  = 0.85;   

  static const double pad2axisTitleSizeY    = 0.1;  
  static const double pad2axisLabelSizeY    = 0.1;   
  static const double pad2axisTitleOffsetY  = 0.55;    

  static const double pad2axisTitleSizeX    = 0.14;   
  static const double pad2axisLabelSizeX    = 0.12;   
  static const double pad2axisTitleOffsetX  = 1.1;    


  static const double axisLabelOffsetX = 0.009;
  static const bool   axisCenterTitle  = false;
  static const bool   axisTickX        = true;
  static const bool   axisTickY        = true;

  // 線條與標記（LINE_WIDTH, MARKER_SIZE）
  static const double defaultLineWidth  = 3.0;
  static const double defaultMarkerSize = 1.3;

  // TLatex（CMS 頂欄）
  static const int    cmsTextFont   = 42;
  static const int    cmsTextAlign  = 11;
  static const double cmsTextSize   = 0.04;
  static const double cmsTextXLeft  = 0.12 + offset; // = 0.13
  static const double cmsTextXRight = 0.82 + offset; // = 0.83
  static const double cmsTextY      = 0.92;

  // gStyle
  static const bool   showStatBox = false;
}

#endif // FLASHGGFINALFIT_BACKGROUND_INTERFACE_PLOTSTYLE_H
