"""
Python 版本的 CMS_lumi，從 CMS_lumi.cpp 簡化移植。

使用方式：
from cmsLumi import CMS_lumi
cms = CMS_lumi()
cms.set_lumi(pad, lumi_fb, iPosX=0, text="Preliminary", drawLumi=True)
"""

from ROOT import TLatex, kBlack

class CMS_lumi:
    def __init__(self):
        # 與原 C++ 版對應的屬性
        self.cmsText = "CMS"
        self.writeExtraText = True
        self.extraTextFont = 52
        self.cmsTextFont = 61
        self.extraOverCmsTextSize = 0.95
        self.lumi_sqrtS = " fb^{-1} (13.6 TeV)"
        self.relPosX = 0.045
        self.relPosY = 0.035
        self.relExtraDY = 1.2
        self.lumiTextSize = 0.5
        self.lumiTextOffset = 0.15
        self.cmsTextSize = 0.55
        self.drawLogo = False  # 保留旗標但未實作 Logo
        # 使用者可調外部文字
        self.defaultExtraText = "Preliminary"

    def set_lumi(self, pad, lumi, iPosX=10, text="Preliminary", drawLumi=True):
        pad.cd()
        outOfFrame = (iPosX // 10 == 0)

        alignY = 3
        alignX = 2
        if iPosX // 10 == 0:
            alignX = 1
        if iPosX == 0:
            alignX = 1
            alignY = 1
        if iPosX // 10 == 1:
            alignX = 1
        if iPosX // 10 == 2:
            alignX = 2
        if iPosX // 10 == 3:
            alignX = 3
        if iPosX == 0:
            self.relPosX = 0.14

        align = 10 * alignX + alignY

        l = pad.GetLeftMargin()
        t = pad.GetTopMargin()
        r = pad.GetRightMargin()
        b = pad.GetBottomMargin()

        if lumi > 100:
            lumiText = f"{lumi:.0f}{self.lumi_sqrtS}"
        else:
            lumiText = f"{lumi:.1f}{self.lumi_sqrtS}"

        latex = TLatex()
        latex.SetNDC()
        latex.SetTextColor(kBlack)

        extraTextSize = self.extraOverCmsTextSize * self.cmsTextSize

        # 繪製 luminosity
        latex.SetTextFont(42)
        latex.SetTextAlign(31)
        latex.SetTextSize(self.lumiTextSize * t)
        if drawLumi:
            latex.DrawLatex(1 - r, 1 - t + self.lumiTextOffset * t, lumiText)
        else:
            latex.DrawLatex(1 - r, 1 - t + self.lumiTextOffset * t, "13.6 TeV")

        # outOfFrame 模式
        if outOfFrame:
            latex.SetTextFont(self.cmsTextFont)
            latex.SetTextAlign(11)
            latex.SetTextSize(self.cmsTextSize * t)
            latex.DrawLatex(l, 1 - t + self.lumiTextOffset * t, self.cmsText)

        pad.cd()

        # 位置計算
        if iPosX % 10 <= 1:
            posX = l + self.relPosX * (1 - l - r)
        elif iPosX % 10 == 2:
            posX = l + 0.5 * (1 - l - r)
        else:
            posX = 1 - r - self.relPosX * (1 - l - r)

        posY = 1 - t - self.relPosY * (1 - t - b)

        # 主要文字 + Extra
        if not outOfFrame:
            latex.SetTextFont(self.cmsTextFont)
            latex.SetTextSize(self.cmsTextSize * t)
            latex.SetTextAlign(align)
            latex.DrawLatex(posX, posY, self.cmsText)

            if self.writeExtraText:
                latex.SetTextFont(self.extraTextFont)
                latex.SetTextAlign(align)
                latex.SetTextSize(extraTextSize * t)
                extraDX = -0.02
                latex.DrawLatex(posX + extraDX * (1 - l - r),
                                posY - self.relExtraDY * self.cmsTextSize * t,
                                text or self.defaultExtraText)
        else:
            if self.writeExtraText:
                if iPosX == 0:
                    posX = l + self.relPosX * (1 - l - r)
                    posY = 1 - t + self.lumiTextOffset * t
                latex.SetTextFont(self.extraTextFont)
                latex.SetTextSize(extraTextSize * t)
                latex.SetTextAlign(align)
                extraDX = -0.03
                latex.DrawLatex(posX + extraDX * (1 - l - r),
                                posY,
                                text or self.defaultExtraText)
