#!/bin/bash

cmsenv

lable='run2'
version='UL'
Lumi_run2='138'

Lumis=( 16.81 19.52 41.48 59.83 )

massList=( 1 2 3 4 5 6 7 8 9 10 15 20 25 30 )
nMass=${#massList[@]}

###### background fit ######

cd ./Background/
make 

dir_out_bkg="./ALP_BkgModel_param_${version}"

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

# mkdir -p "$path_out_bkg/9"
# path_bkg="$path_out_bkg/9"
# ./bin/fTest_ALP_turnOn -i $path_in_bkg/ALP_data_bkg_Am1_workspace.root --saveMultiPdf $path_bkg/CMS-HGG_mva_13TeV_multipdf.root -D $path_bkg/HZAmassInde_fTest --mass_ALP 9 -c 1 --isFlashgg 0 --isData 0 -f data, --mhLow 95 --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135 > $path_bkg/ftest.log
# exit 

for ((iBin=0; iBin<$nMass; iBin++))
# for ((iBin=0; iBin<1; iBin++))
    do
    mkdir -p "$path_out_bkg/${massList[$iBin]}"
    path_bkg="$path_out_bkg/${massList[$iBin]}"
    
    # ./bin/makeBkgPlots_ALP -b $path_bkg/CMS-HGG_mva_13TeV_multipdf.root -d $path_bkg/BkgPlots --total_OutDir $total_OutDir -o $path_bkg/BkgPlots.root -S 13 --isMultiPdf --useBinnedData --unblind --massStep 2.5 --mhVal 125.0 --maVal ${massList[$iBin]} --mhLow 95 --mhHigh 180 --intLumi $Lumi_run2 -c 0 --isFlashgg 0
    # ./bin/fTest_ALP_turnOn -i $path_in_bkg/ALP_data_bkg_Am${massList[$iBin]}_workspace.root --saveMultiPdf $path_bkg/CMS-HGG_mva_13TeV_multipdf.root -D $path_bkg/HZAmassInde_fTest --mass_ALP ${massList[$iBin]} -c 1 --isFlashgg 0 --isData 0 -f data, --mhLow 95 --mhHigh 180  --mhLowBlind 115 --mhHighBlind 135 > $path_bkg/ftest.log

    /afs/cern.ch/work/p/pelai/HZa/CMSSW_14_1_0_pre4/src/flashggFinalFit/Background/bin/makeBkgPlots_ALP -b $path_bkg/CMS-HGG_mva_13TeV_multipdf.root -d $path_bkg/BkgPlots --total_OutDir $total_OutDir -o $path_bkg/BkgPlots.root -S 13 --isMultiPdf --useBinnedData --unblind --massStep 2.5 --mhVal 125.0 --maVal ${massList[$iBin]} --mhLow 95 --mhHigh 180 --intLumi $Lumi_run2 -c 0 --isFlashgg 0

    done