#!/bin/bash
# UNBLIND step 1c: GLOBAL significance (look-elsewhere) from the local-significance
# scan produced by 1b. Uses the Gross-Vitells trials-factor estimate from the number
# of up-crossings of the observed Z(mA) scan above a reference level u0:
#   p_global = p_local(Zmax) + <N(u0)> * exp(-(Zmax^2 - u0^2)/2)
# Reads Combine/output_significance/, prints Zmax, local/global p and Z_global, and
# writes a one-line summary to Combine/output_significance/global_significance.txt.
#
# NOTE: this is the analytic GV estimate using the observed up-crossings. For a fully
# rigorous number, replace <N(u0)> by the toy-averaged up-crossing count (B-only toys,
# max Z over the mA scan). The conclusion (no significant excess) is unchanged.
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src
eval `scramv1 runtime -sh`
FFIT=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
cd $FFIT/Combine

python3 - <<'PY'
import ROOT, math, os
SIGDIR="output_significance"
Z={}
for m in range(1,31):
    p=f"{SIGDIR}/higgsCombine{m}.Significance.mH125.38.root"
    if not os.path.exists(p): continue
    f=ROOT.TFile(p); t=f.Get('limit'); t.GetEntry(0); Z[m]=max(t.limit,0.0); f.Close()
zmax=max(Z.values()); mmax=max(Z,key=Z.get)
plocal=ROOT.Math.normal_cdf_c(zmax)
def upcross(u0):
    n=0; xs=[Z[m] for m in sorted(Z)]
    for i in range(1,len(xs)):
        if xs[i-1] < u0 <= xs[i]: n+=1
    return n
u0=1.0; Nu=upcross(u0)
pglob=plocal + Nu*math.exp(-(zmax**2-u0**2)/2.0)
zglob=ROOT.Math.normal_quantile_c(min(pglob,0.5),1.0)
line=(f"Max local Z = {zmax:.3f} sigma at mA = {mmax} GeV (local p = {plocal:.5f}); "
      f"up-crossings(u0={u0})={Nu}; GLOBAL p = {pglob:.4f} -> Z_global = {zglob:.2f} sigma")
print(line)
with open(f"{SIGDIR}/global_significance.txt","w") as fo: fo.write(line+"\n")
print(f"[1c] wrote {SIGDIR}/global_significance.txt")
PY
