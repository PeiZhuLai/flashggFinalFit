#!/bin/bash
/bin/hostname
gcc -v
pwd
# export PATH=$PATH:/afs/ihep.ac.cn/soft/common/sysgroup/hep_job/bin/
source /cvmfs/cms.cern.ch/cmsset_default.sh

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/eos/home-p/pelai/anaconda/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/eos/home-p/pelai/anaconda/etc/profile.d/conda.sh" ]; then
        . "/eos/home-p/pelai/anaconda/etc/profile.d/conda.sh"
    else
        export PATH="/eos/home-p/pelai/anaconda/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<
conda init

conda activate python_Three

path_code='/afs/cern.ch/work/p/pelai/HZa/CMSSW_14_1_0_pre4/src/flashggFinalFit/MVAcut/run2_UL'

python $path_code/makeWorkspace_data.py -m $1

python $path_code/makeWorkspace_sig.py -m $1

python $path_code/makeWorkspace_sig.py -m $1 --ele

python $path_code/makeWorkspace_sig.py -m $1 --mu

# python3 $path_code/makeWorkspace_data.py -m 3
