ALP Flashgg Final Fit
========================

Please follow the entire instruction from offical flashggfit 

```
export SCRAM_ARCH=el9_amd64_gcc12
cmsrel CMSSW_14_1_0_pre4
cd CMSSW_14_1_0_pre4/src
cmsenv

COMBINE_TAG=07b56c67ba6e4304b42c3a6cdba710d59c719192
COMBINEHARVESTER_TAG=94017ba5a3a657f7b88669b1a525b19d34ea41a2
FINALFIT_TAG=higgsdnafinalfit

# Install Combine with the latest EL9 compatible branch
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit && git fetch origin ${COMBINE_TAG} && git checkout ${COMBINE_TAG}

# Install CombineTools in CombineHarvester
cd ${CMSSW_BASE}/src
bash <(curl -s https://raw.githubusercontent.com/cms-analysis/CombineHarvester/${COMBINEHARVESTER_TAG}/CombineTools/scripts/sparse-checkout-https.sh)
cd CombineHarvester && git fetch origin ${COMBINEHARVESTER_TAG} && git checkout ${COMBINEHARVESTER_TAG}

# Compile libraries
cd ${CMSSW_BASE}/src
cmsenv
scram b clean
scram b -j 8

# Install Final Fit package
git clone -b $FINALFIT_TAG https://github.com/cms-analysis/flashggFinalFit.git
cd flashggFinalFit/
source setup.sh
```

Add the modified file into counterpart of fiels

```
cd Background

cmsenv

make clean

make

```

========================
Cloning the Repository
---------------
```
cmssw-el7

at IHEP
export PATH=/cvmfs/container.ihep.ac.cn/bin:$PATH
hep_container shell CentOS7

export SCRAM_ARCH=slc7_amd64_gcc700

cmsrel CMSSW_10_2_13

cd CMSSW_10_2_13/src

cmsenv

git cms-init

git clone git@github.com:PeiZhuLai/cmssw.git
```

Install the GBRLikelihood package which contains the RooDoubleCBFast implementation
```
git clone git@github.com:jonathon-langford/HiggsAnalysis.git
```
Install Combine as per the documentation here: cms-analysis.github.io/HiggsAnalysis-CombinedLimit/
```
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v8.2.0
```
Compile external libraries
-----------------------
```
cd ../HiggsAnalysis

cmsenv

scramv1 b clean; scramv1 b
```
Install Flashgg Final Fit packages
-----------------------
```
cd ../..

git clone https://github.com/zebingwang/flashggFinalFit.git
git clone git@github.com:PeiZhuLai/flashggFinalFit.git

install flashggFinalFit parallel to HiggsAnalysis
```

Compile Signal Model Pakage
--------------

```
cd flashggFinalFit

cd Signal

cmsenv

make clean

make

```

Compile Background Model Pakage
--------------

```
cd Background

cmsenv

make clean

make

```

Prepare dataset
--------------

Copy all the input datasets into dir: ```InputData```


Run
-----

```
./fit_HZGamma_model_UL.sh
```
