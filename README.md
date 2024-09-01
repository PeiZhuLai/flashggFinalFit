# ALP Flashgg Final Fit

run all the readme in official flashggFinalFit Github

https://github.com/cms-analysis/flashggFinalFit?tab=readme-ov-file##download-and-setup-instructions

Cloning the Repository
-----------------------
```
export SCRAM_ARCH=slc7_amd64_gcc700

cmssw-el7

cmsrel CMSSW_10_2_13

cd CMSSW_10_2_13/src

cmsenv

git cms-init

git clone git@github.com:PeiZhuLai/cmssw.git
```

Install the GBRLikelihood package which contains the RooDoubleCBFast implementation
```
git clone git@github.com:jonathon-langford/HiggsAnalysis.git
git clone https://github.com/jonathon-langford/HiggsAnalysis.git
git clone git@github.com:PeiZhuLai/HiggsAnalysis.git
```
Install Combine as per the documentation here: cms-analysis.github.io/HiggsAnalysis-CombinedLimit/
```
git clone git@github.com:cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
git clone git@github.com:PeiZhuLai/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
```
Compile external libraries
-----------------------
```
cd HiggsAnalysis

cmsenv

scram b -j

cd $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit

git fetch origin

git checkout v10.0.2

scramv1 b clean; scramv1 b # always make a clean build

```
Install Flashgg Final Fit packages
-----------------------
```
cd ..

git clone https://github.com/PeiZhuLai/flashggFinalFit.git

git clone git@github.com:PeiZhuLai/flashggFinalFit.git

cd flashggFinalFit/
```

```
cd Signal

cmsenv

make clean

make

```

Background Model
--------------

```
cd Background

cmsenv

make clean

make

```
Run
-----

```
./fit_ALP.sh
```
