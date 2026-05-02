#!/bin/bash
# /bin/hostname
# gcc -v
# pwd
# For IHEP 
# export PATH=$PATH:/afs/ihep.ac.cn/soft/common/sysgroup/hep_job/bin/

source /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/setup.sh
export PYTHONPATH=$PYTHONPATH:/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/tools:/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Signal/tools

### ---- Signal With Uncertainty ----
mAs_sig=(1 2 3 4 5 6 7 8 9 10 15 20 25 30)
sig_samples=()
for m in "${mAs_sig[@]}"; do
    sig_samples+=("mA_M${m}")
done
years_sig=(2022preEE 2022postEE 2023preBPix 2023postBPix 2024)  # 信號
DO_SYSTEMATICS=1  # 設為 1 時會加上 --doSystematics
leps=(ele mu)     # 逐個 lepton 輸出各自的 ws 檔

for year in "${years_sig[@]}"; do
    for ma in "${sig_samples[@]}"; do
        for lep in "${leps[@]}"; do
            path="/eos/home-p/pelai/HZa/root_MVAcut/sig/${ma}"
            extra_args=()
            if [ "${DO_SYSTEMATICS}" = "1" ]; then
                extra_args+=(--doSystematics)
            fi
            python3 trees2ws.py --inputConfig config.py --inputTreeFile "${path}/output_${year}.root" --inputMass 125 --productionMode ggh --year "${year}" --lepton "${lep}" "${extra_args[@]}"
        done
    done
done

# # ---- Data ----
mAs_data=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30)
data_samples=()
for m in "${mAs_data[@]}"; do
    data_samples+=("mA_M${m}")
done
for ma in "${data_samples[@]}"; do
    path="/eos/home-p/pelai/HZa/root_MVAcut/data/${ma}"
    python3 trees2ws_data.py --inputConfig config.py --inputTreeFile "${path}/run3.root"
done


# # For CERN
# source /cvmfs/cms.cern.ch/cmsset_default.sh

# export PATH=/eos/home-p/pelai/App/Anaconda/Anaconda/Install/anaconda3/bin:$PATH
# source /eos/home-p/pelai/App/Anaconda/Anaconda/Install/anaconda3/etc/profile.d/conda.sh

# # >>> conda initialize >>>
# # !! Contents within this block are managed by 'conda init' !!
# __conda_setup="$('/eos/home-p/pelai/App/Anaconda/Anaconda/Install/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
# if [ $? -eq 0 ]; then
#     eval "$__conda_setup"
# else
#     if [ -f "/eos/home-p/pelai/App/Anaconda/Anaconda/Install/anaconda3/etc/profile.d/conda.sh" ]; then
#         . "/eos/home-p/pelai/App/Anaconda/Anaconda/Install/anaconda3/etc/profile.d/conda.sh"
#     else
#         export PATH="/eos/home-p/pelai/App/Anaconda/Anaconda/Install/anaconda3/bin:$PATH"
#     fi
# fi
# unset __conda_setup
# # <<< conda initialize <<<

# conda activate higgs-alp-ana

# path_code='/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/MVAcut/run3_ReReco'

# python $path_code/makeWorkspace_data.py -m $1

# ---- Signal Without Uncertainty ----
# python $path_code/makeWorkspace_sig.py -m $1
# python $path_code/makeWorkspace_sig.py -m $1 --ele
# python $path_code/makeWorkspace_sig.py -m $1 --mu

# Commission Test
# python3 $path_code/makeWorkspace_data.py -m 1
# python $path_code/makeWorkspace_sig.py -m 5
# python $path_code/makeWorkspace_sig.py -m 5 --ele
# python $path_code/makeWorkspace_sig.py -m 5 --mu