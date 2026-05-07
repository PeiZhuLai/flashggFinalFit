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
max_parallel=6
logDir="/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Trees2WS/logs"
mkdir -p "${logDir}"

pids=()
labels=()
batch_id=1
job_in_batch=0
status=0

wait_batch() {
    if [ "${#pids[@]}" -eq 0 ]; then
        return
    fi

    echo "[INFO] Waiting for Tree2WS batch ${batch_id} (${#pids[@]} jobs)"
    for i in "${!pids[@]}"; do
        if wait "${pids[$i]}"; then
            echo "[INFO] Finished ${labels[$i]}"
        else
            echo "[ERROR] Failed ${labels[$i]}; see ${logDir}/${labels[$i]}.log" >&2
            status=1
        fi
    done

    pids=()
    labels=()
    job_in_batch=0
    batch_id=$((batch_id + 1))
}

submit_job() {
    local label="$1"
    shift

    job_in_batch=$((job_in_batch + 1))
    echo "[INFO] Start Tree2WS batch ${batch_id} job ${job_in_batch}/${max_parallel}: ${label}"
    "$@" > "${logDir}/${label}.log" 2>&1 &
    pids+=("$!")
    labels+=("${label}")

    if [ "${job_in_batch}" -eq "${max_parallel}" ]; then
        wait_batch
    fi
}

for ma in "${sig_samples[@]}"; do
    for year in "${years_sig[@]}"; do
        for lep in "${leps[@]}"; do
            path="/eos/home-p/pelai/HZa/root_MVAcut/sig/${ma}"
            extra_args=()
            if [ "${DO_SYSTEMATICS}" = "1" ]; then
                extra_args+=(--doSystematics)
            fi
            submit_job "sig_${ma}_${year}_${lep}" \
                python3 trees2ws.py --inputConfig config.py --inputTreeFile "${path}/output_${year}.root" --inputMass 125 --productionMode ggh --year "${year}" --lepton "${lep}" "${extra_args[@]}"
        done
    done
done

wait_batch

# # ---- Data ----
mAs_data=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30)
data_samples=()
for m in "${mAs_data[@]}"; do
    data_samples+=("mA_M${m}")
done
for ma in "${data_samples[@]}"; do
    path="/eos/home-p/pelai/HZa/root_MVAcut/data/${ma}"
    submit_job "data_${ma}" \
        python3 trees2ws_data.py --inputConfig config.py --inputTreeFile "${path}/run3.root"
done

wait_batch
exit "${status}"


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
