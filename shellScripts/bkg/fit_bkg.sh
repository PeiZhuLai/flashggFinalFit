#!/bin/bash
# 避免以 sh 執行（需要 Bash）
if [ -z "${BASH_VERSION:-}" ]; then
  echo "[ERROR] Please run with bash: bash $0"
  exit 2
fi
# set -euo pipefail

cmsenv

source /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/setup.sh

lable='run3'
version='ReReco'
Lumi_run3='170.84'

Lumis=( 7.98 27.01 17.61 9.53 108.95)

# massList=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 )
massList=( 14 21 23 29 )
# massList=( 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 )
# massList=( 16 17 18 19 )
# massList=( 14 )
nMass=${#massList[@]}

###### background fit ######

# 設定基礎工作目錄
BASE_DIR="/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit"

cd $BASE_DIR/Background/
# 明確指定使用小寫 makefile
make -f makefile clean && make -f makefile -B
# make clean && make vars && make -j
# make

# 檢查目標二進位是否存在
if [[ ! -x ./bin/fTest_ALP_turnOn || ! -x ./bin/makeBkgPlots_ALP ]]; then
  echo "[ERROR] build failed or binaries missing in ./bin"
  ls -l ./bin || true
  exit 1
fi

# dir_input="../MVAcut/${lable}_${version}/output/data"
dir_input="/eos/home-p/pelai/HZa/root_MVAcut/data"
dir_output="./ALP_BkgModel_${version}"

if [ ! -d $dir_output ];then
  mkdir $dir_output
  echo "creat $dir_output"
else
  echo "$dir_output already exist"
fi

path_out_bkg="$dir_output/fit_results_${lable}"
mkdir -p $path_out_bkg
mkdir -p $path_out_bkg/AllFitResults
total_OutDir="$path_out_bkg/AllFitResults"


# mkdir -p "$path_out_bkg/5"
# path_bkg="$path_out_bkg/5"
# ./bin/fTest_ALP_turnOn -i $dir_input/mA_M5/ws/run3.root --saveMultiPdf $path_bkg/CMS-HGG_mva_13TeV_multipdf.root -D $path_bkg/HZAmassInde_fTest --mass_ALP 5 -c 1 --isFlashgg 0 --isData 0 -f data, --mhLow 95 --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135 > $path_bkg/ftest.log
# ./bin/makeBkgPlots_ALP -b $path_bkg/CMS-HGG_mva_13TeV_multipdf.root -d $path_bkg/BkgPlots --total_OutDir $total_OutDir -o $path_bkg/BkgPlots.root --sqrts 13p6TeV --isMultiPdf --useBinnedData --massStep 2.5 --mhVal 125.0 --maVal 30 --mhLow 95 --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135 --intLumi $Lumi_run3 -c 0 --isFlashgg 0 --doBands
# exit 

failed_log="failed_mass_points.log"
: > "$failed_log"   # 清空舊 log

summary_log="$path_out_bkg/background_fit_summary.csv"
echo "mass,ftest_status,best_fit_pdf,n_pdfs,ftest_results,envelope_results,ftest_stdout,bkgplots_status,bkgplots_stdout" > "$summary_log"

for ((iBin=0; iBin<$nMass; iBin++))
# for ((iBin=0; iBin<1; iBin++))
    do
    echo ">>> Processing mass point: ${massList[$iBin]}"

    mkdir -p "$path_out_bkg/${massList[$iBin]}"
    path_bkg="$path_out_bkg/${massList[$iBin]}"
    ftest_outdir="$path_bkg/HZAmassInde_fTest"
    ftest_stdout="$path_bkg/ftest.stdout.log"
    bkgplots_stdout="$path_bkg/makeBkgPlots.stdout.log"
    ftest_results="$ftest_outdir/fTestResults.txt"
    envelope_results="$ftest_outdir/EnvelopeResults.txt"
    best_fit_pdf="NA"
    n_pdfs="NA"
    ftest_status="OK"
    bkgplots_status="OK"

    # Syst
    ######################################
    # 1. fTest
    ######################################
    # ./bin/fTest_ALP_turnOn -i $dir_input/mA_M${massList[$iBin]}/ws/run3.root --saveMultiPdf $path_bkg/CMS-HGG_mva_13p6TeV_multipdf.root -D $path_bkg/HZAmassInde_fTest --mass_ALP ${massList[$iBin]} -c 1 --isFlashgg 0 --isData 0 -f data, --mhLow 95 --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135 > $path_bkg/ftest.log
    ./bin/fTest_ALP_turnOn -i $dir_input/mA_M${massList[$iBin]}/ws/run3.root --saveMultiPdf $path_bkg/CMS-HGG_mva_13p6TeV_multipdf.root -D $ftest_outdir --mass_ALP ${massList[$iBin]} -c 1 --isFlashgg 0 --isData 0 -f data, --mhLow 95 --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135 > "$ftest_stdout" 2>&1
    
    if [[ $? -ne 0 ]]; then
        ftest_status="FAIL"
        echo "[FAIL][fTest] mA=${massList[$iBin]}" | tee -a "$failed_log"
        printf "%s,%s,%s,%s,%s,%s,%s,%s,%s\n" \
          "${massList[$iBin]}" "$ftest_status" "$best_fit_pdf" "$n_pdfs" "$ftest_results" "$envelope_results" "$ftest_stdout" "$bkgplots_status" "$bkgplots_stdout" >> "$summary_log"
        continue
    fi

    if [[ -f "$ftest_stdout" ]]; then
        best_fit_pdf=$(grep -E "Best Fit Pdf =" "$ftest_stdout" | tail -1 | sed 's/.*Best Fit Pdf = //' | sed 's/, /:/' || true)
        n_pdfs=$(grep -E "with a total of [0-9]+ pdfs" "$ftest_stdout" | tail -1 | sed -E 's/.*with a total of ([0-9]+) pdfs.*/\1/' || true)
        [[ -n "$best_fit_pdf" ]] || best_fit_pdf="NA"
        [[ -n "$n_pdfs" ]] || n_pdfs="NA"
    fi
    
    ######################################
    # 2. makeBkgPlots
    ######################################
    # makeBkgPlots 只負責產圖與 band，可用來診斷；datacard 真正讀的是上面的 CMS-HGG_mva_13p6TeV_multipdf.root
    ./bin/makeBkgPlots_ALP -b $path_bkg/CMS-HGG_mva_13p6TeV_multipdf.root -d $path_bkg/BkgPlots --total_OutDir $total_OutDir -o $path_bkg/BkgPlots.root --sqrts 13p6TeV --isMultiPdf --useBinnedData --massStep 2.5 --mhVal 125.0 --maVal ${massList[$iBin]} --mhLow 95 --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135 --intLumi $Lumi_run3 -c 0 --isFlashgg 0 --doBands > "$bkgplots_stdout" 2>&1

    if [[ $? -ne 0 ]]; then
        bkgplots_status="FAIL"
        echo "[FAIL][BkgPlots] mA=${massList[$iBin]}" | tee -a "$failed_log"
        printf "%s,%s,%s,%s,%s,%s,%s,%s,%s\n" \
          "${massList[$iBin]}" "$ftest_status" "$best_fit_pdf" "$n_pdfs" "$ftest_results" "$envelope_results" "$ftest_stdout" "$bkgplots_status" "$bkgplots_stdout" >> "$summary_log"
        continue
    fi

    printf "%s,%s,%s,%s,%s,%s,%s,%s,%s\n" \
      "${massList[$iBin]}" "$ftest_status" "$best_fit_pdf" "$n_pdfs" "$ftest_results" "$envelope_results" "$ftest_stdout" "$bkgplots_status" "$bkgplots_stdout" >> "$summary_log"

    echo "[OK] mA=${massList[$iBin]}"
    # Nominal
    # 1
    # ./bin/fTest_ALP_turnOn -i $dir_input/ALP_data_bkg_Am${massList[$iBin]}_workspace.root --saveMultiPdf $path_bkg/CMS-HGG_mva_13TeV_multipdf.root -D $path_bkg/HZAmassInde_fTest --mass_ALP ${massList[$iBin]} -c 1 --isFlashgg 0 --isData 0 -f data, --mhLow 95 --mhHigh 180  --mhLowBlind 115 --mhHighBlind 135 > $path_bkg/ftest.log
    # 2
    # ./bin/makeBkgPlots_ALP -b $path_bkg/CMS-HGG_mva_13TeV_multipdf.root -d $path_bkg/BkgPlots --total_OutDir $total_OutDir -o $path_bkg/BkgPlots.root -S 13 --isMultiPdf --useBinnedData --massStep 2.5 --mhVal 125.0 --maVal ${massList[$iBin]} --mhLow 95 --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135 --intLumi $Lumi_run3 -c 0 --isFlashgg 0

    # Blind (Adding --mhLowBlind 115 --mhHighBlind 135)
    # careful --unblind whether added
    # ./bin/makeBkgPlots_ALP -b $path_bkg/CMS-HGG_mva_13TeV_multipdf.root -d $path_bkg/BkgPlots --total_OutDir $total_OutDir -o $path_bkg/BkgPlots.root -S 13 --isMultiPdf --useBinnedData --unblind --massStep 2.5 --mhVal 125.0 --maVal ${massList[$iBin]} --mhLow 95 --mhHigh 180 --intLumi $Lumi_run3 -c 0 --isFlashgg 0

    done
