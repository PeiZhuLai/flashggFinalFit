#!/bin/bash

cmsenv

lable='run3'
version='ReReco'
Lumi_run3='62.5'

Lumis=( 7.98 27.01 17.61 9.53 )

# massList=( 1 2 3 4 5 6 7 8 9 10 15 20 25 30 )
massList=( 5 15 30 )
nMass=${#massList[@]}

###### background fit ######

# 設定基礎工作目錄
BASE_DIR="/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit"

cd $BASE_DIR/Background/
make 

dir_out_bkg="./ALP_BkgModel_${version}"

if [ ! -d $dir_out_bkg ];then
  mkdir $dir_out_bkg
  echo "creat $dir_out_bkg"
else
  echo "$dir_out_bkg already exist"
fi

path_in_bkg="../MVAcut/${lable}_${version}/output/data"

path_out_bkg="$dir_out_bkg/fit_results_${lable}"
mkdir -p $path_out_bkg
mkdir -p $path_out_bkg/AllFitResults
total_OutDir="$path_out_bkg/AllFitResults"


# mkdir -p "$path_out_bkg/30"
# path_bkg="$path_out_bkg/30"
# ./bin/fTest_ALP_turnOn -i $path_in_bkg/ALP_data_bkg_Am30_workspace.root --saveMultiPdf $path_bkg/CMS-HGG_mva_13TeV_multipdf.root -D $path_bkg/HZAmassInde_fTest --mass_ALP 30 -c 1 --isFlashgg 0 --isData 0 -f data, --mhLow 95 --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135 > $path_bkg/ftest.log
# exit 


for ((iBin=0; iBin<$nMass; iBin++))
# for ((iBin=0; iBin<1; iBin++))
    do
    mkdir -p "$path_out_bkg/${massList[$iBin]}"
    path_bkg="$path_out_bkg/${massList[$iBin]}"

    # 1
    # ./bin/fTest_ALP_turnOn -i $path_in_bkg/ALP_data_bkg_Am${massList[$iBin]}_workspace.root --saveMultiPdf $path_bkg/CMS-HGG_mva_13TeV_multipdf.root -D $path_bkg/HZAmassInde_fTest --mass_ALP ${massList[$iBin]} -c 1 --isFlashgg 0 --isData 0 -f data, --mhLow 95 --mhHigh 180  --mhLowBlind 115 --mhHighBlind 135 > $path_bkg/ftest.log
    # 2
    ./bin/makeBkgPlots_ALP -b $path_bkg/CMS-HGG_mva_13TeV_multipdf.root -d $path_bkg/BkgPlots --total_OutDir $total_OutDir -o $path_bkg/BkgPlots.root -S 13 --isMultiPdf --useBinnedData --massStep 2.5 --mhVal 125.0 --maVal ${massList[$iBin]} --mhLow 95 --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135 --intLumi $Lumi_run3 -c 0 --isFlashgg 0

    # Blind (Adding --mhLowBlind 115 --mhHighBlind 135)
    # careful --unblind whether added
    # ./bin/makeBkgPlots_ALP -b $path_bkg/CMS-HGG_mva_13TeV_multipdf.root -d $path_bkg/BkgPlots --total_OutDir $total_OutDir -o $path_bkg/BkgPlots.root -S 13 --isMultiPdf --useBinnedData --massStep 2.5 --mhVal 125.0 --maVal ${massList[$iBin]} --mhLow 95 --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135 --intLumi $Lumi_run3 -c 0 --isFlashgg 0
    # /afs/cern.ch/work/p/pelai/HZa/CMSSW_14_1_0_pre4/src/flashggFinalFit/Background/bin/makeBkgPlots_ALP -b $path_bkg/CMS-HGG_mva_13TeV_multipdf.root -d $path_bkg/BkgPlots --total_OutDir $total_OutDir -o $path_bkg/BkgPlots.root -S 13 --isMultiPdf --useBinnedData --unblind --massStep 2.5 --mhVal 125.0 --maVal ${massList[$iBin]} --mhLow 95 --mhHigh 180 --intLumi $Lumi_run3 -c 0 --isFlashgg 0

    done