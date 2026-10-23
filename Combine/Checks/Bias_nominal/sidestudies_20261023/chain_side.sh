#!/usr/bin/env bash
# Side studies of the background bias, run next to the production envelope revert
# (HiggsZaAna/studies/l3_review_20260910/deploy_envrevert_20261023/chain_prod.sh).
#
#   truthout  Truth function outside the fitted envelope: toys generated from Bern6 (removed by the Bern5 cap at
#             m_a = 20 and 25; workspaces of envcheck_20261022/t2w) and fitted with the production envelope
#             (Combine/root_t2w). m_a = 25 waits for the reverted production workspace (chain_prod stage 4).
#   ma29      m_a = 29 minimizer settings (Bern5 loses ~13% of its toys to failed fits in production):
#             fallback  = production options + --cminFallbackAlgo Minuit2,1:1.0 + Minuit2,Simplex,0:1.0 (10 jobs)
#             fbrobust  = fallback + --robustFit 1, one job per (chunk, truth function) as it is ~5x slower (40 jobs)
#             Same seeds as production, so every variant fits the same toys.
#   r2        Recheck of the two envelope restrictions not covered by envcheck_20261022, with the same method:
#             v11full = m_a = 11 without the low-order restriction; v11b3 = low-order with Bernstein up to 3;
#             v27b6 = Bernstein back at m_a = 27 (cap 6). fTest with the envcheck build (env overrides), one process.
#             A variant whose envelope equals production is not submitted.
#
# All condor jobs: production runjob_bias_highstat.sh, 10 chunks x 1000 toys, stdout to /dev/null.
# Markers: logs/side.DONE / logs/side.FAILED; heartbeat logs/heartbeat.txt; registry ../PENDING_RESULTS.txt.
set -uo pipefail
S=$(cd "$(dirname "$0")" && pwd)
CMSSW_TOP=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4
F=$CMSSW_TOP/src/flashggFinalFit
C=$F/Combine
BN=$C/Checks/Bias_nominal
EC=$BN/envcheck_20261022
FR=$F/Background/ALP_BkgModel_ReReco/fit_results_run3
DEP=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/deploy_envrevert_20261023
MV=/eos/home-p/pelai/HZa/root_MVAcut/data
REG=$BN/PENDING_RESULTS.txt
EXE=$F/shellScripts/bias/Condor/runjob_bias_highstat.sh
MERGE=$F/shellScripts/bias/Condor/merge_bias_highstat.sh
L=$S/logs
mkdir -p $L/condor
hb() { echo "$(date '+%F %T') $*" | tee -a $L/heartbeat.txt; }
die() { hb "STOPPING: $*"; echo "$(date '+%F %T') $*" > $L/side.FAILED; echo "$(date '+%F %H:%M') DIED  HZa bias side studies: $*  -> $L/heartbeat.txt" >> $REG; exit 1; }
envelope() { grep -oE "pdf : [A-Za-z0-9]+ , gof : [0-9.e-]+" $1 | sed 's/pdf : //; s/ , gof : /:/' | tr '\n' ' '; }
names() { envelope $1 | sed -E 's/:[^ ]*//g'; }
FALLBACK="--cminFallbackAlgo Minuit2,1:1.0 --cminFallbackAlgo Minuit2,Simplex,0:1.0"
rm -f $L/side.DONE $L/side.FAILED

set +u
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd $CMSSW_TOP/src && eval $(scramv1 runtime -sh)
set -u

# submit NAME "ENV" "MASS[:FUNC] ..."   -> prints the cluster id; one job per chunk and entry
declare -A CLUSTER MENV MMASS
submit() {
  local name=$1 envs=$2 entries=$3 sub=$S/$1.submit e m fn c
  cat > $sub <<EOF
universe              = vanilla
executable            = $EXE
getenv                = True
request_memory        = 4000
request_cpus          = 1
+JobFlavour           = "tomorrow"
batch_name            = Bias_side_$name
environment           = "NTOYS_PER_CHUNK=1000 $envs"
on_exit_hold          = (ExitBySignal == True) || (ExitCode != 0)
periodic_release      = (NumJobStarts < 3) && ((CurrentTime - EnteredCurrentStatus) > 600)

EOF
  for e in $entries; do
    m=${e%%:*}; fn=""; [[ $e == *:* ]] && fn=${e#*:}
    for c in $(seq 0 9); do
      cat >> $sub <<EOF
log                   = $L/condor/${name}_mA${m}${fn:+_$fn}_c$c.\$(ClusterId).log
output                = /dev/null
error                 = $L/condor/${name}_mA${m}${fn:+_$fn}_c$c.\$(ClusterId).err
arguments             = $m $c $F $CMSSW_TOP
EOF
      # per-job environment, written outside the heredoc: inside ${:+} the shell drops the double quotes and
      # condor then reads the old unquoted syntax, which leaves BIAS_HS_ROOT unset
      [ -n "$fn" ] && echo "environment           = \"NTOYS_PER_CHUNK=1000 $envs SELECT_FUNCTION=$fn\"" >> $sub
      printf 'queue\n\n' >> $sub
    done
  done
  local out; out=$(condor_submit $sub 2>&1); echo "$out" > $L/submit_$name.log
  grep -oP 'submitted to cluster \K[0-9]+' <<< "$out" | head -1
}
launch() {   # launch NAME "ENV" "ENTRIES" "MERGE_ENV"; a cluster already listed in clusters.txt is reused (restart)
  local cl; cl=$(awk -v n=$1 '$1==n{print $2}' $L/clusters.txt 2>/dev/null | tail -1)
  if [ -n "$cl" ]; then
    CLUSTER[$1]=$cl; MENV[$1]=$4; MMASS[$1]=$(for e in $3; do echo ${e%%:*}; done | sort -un | tr '\n' ' ')
    hb "    $1: reusing cluster $cl"; return
  fi
  cl=$(submit "$1" "$2" "$3")
  [ -n "$cl" ] || die "submission $1 failed ($L/submit_$1.log)"
  # every job must carry the study's output root, otherwise it would write into the production bias outputs
  local nbad; nbad=$(condor_q $cl -af Environment 2>/dev/null | grep -vcE "(^| )BIAS_HS_ROOT=$S/")
  if [ "$nbad" -ne 0 ]; then condor_rm $cl > /dev/null 2>&1; die "submission $1: $nbad jobs without the study BIAS_HS_ROOT (removed $cl)"; fi
  CLUSTER[$1]=$cl; MENV[$1]=$4; MMASS[$1]=$(for e in $3; do echo ${e%%:*}; done | sort -un | tr '\n' ' ')
  echo "$1 $cl" >> $L/clusters.txt
  hb "    $1: cluster $cl, $(( $(wc -w <<< "$3") * 10 )) jobs, masses ${MMASS[$1]}"
}

# ---------------------------------------------------------------- truthout m_a = 20, ma29
hb "[1] truthout m_a = 20 and m_a = 29 minimizer variants"
TO=$S/truthout/outputs; mkdir -p $TO
launch truthout20 "BIAS_HS_ROOT=$TO ROOT_DATACARD_PATH=$EC/t2w FIT_DATACARD_PATH=$C/root_t2w SELECT_FUNCTION=Bern6" "20" \
  "BIAS_HS_ROOT=$TO ROOT_DATACARD_PATH=$EC/t2w SELECT_FUNCTION=Bern6"
mkdir -p $S/ma29/fallback $S/ma29/fbrobust
launch ma29fallback "BIAS_HS_ROOT=$S/ma29/fallback BIAS_FIT_EXTRA_OPTIONS='$FALLBACK'" "29" "BIAS_HS_ROOT=$S/ma29/fallback"
launch ma29fbrobust "BIAS_HS_ROOT=$S/ma29/fbrobust BIAS_FIT_EXTRA_OPTIONS='$FALLBACK --robustFit 1'" "29:Bern5 29:Exp3 29:Lau1 29:Pow1" \
  "BIAS_HS_ROOT=$S/ma29/fbrobust"

# ---------------------------------------------------------------- r2: envelopes, cards, workspaces
hb "[2] r2 envelopes (one fTest at a time)"
declare -A VMASS=( [v11full]=11 [v11b3]=11 [v27b6]=27 )
declare -A VENV=( [v11full]="FULL_ORDER_MASSES=11" [v11b3]="LOWORDER_BERN_EXCL=11:4" [v27b6]="BERN_ORDER_OVERRIDE=27:6" )
for v in v11full v11b3 v27b6; do
  m=${VMASS[$v]}; P=$S/r2/$v/bkg/$m
  rm -rf $P; mkdir -p $P/dat
  ( cd $P && env ${VENV[$v]} $EC/build/bin/fTest_ALP_turnOn -i $MV/mA_M$m/ws/run3.root \
      --saveMultiPdf $P/CMS-HGG_mva_13p6TeV_multipdf.root -D $P/HZAmassInde_fTest --mass_ALP $m -c 1 --isFlashgg 0 --isData 0 \
      -f data, --mhLow 95 --mhHigh 180 --mhLowBlind 115 --mhHighBlind 135 --gtoys 100 -d $P/dat/fTest.dat ) > $P/ftest.stdout.log 2>&1
  [ -s $P/CMS-HGG_mva_13p6TeV_multipdf.root ] || die "r2 $v: no multipdf"
  a=$(envelope $FR/$m/HZAmassInde_fTest/EnvelopeResults.txt); b=$(envelope $P/HZAmassInde_fTest/EnvelopeResults.txt)
  hb "    $v mA$m production: $a| variant: $b"
  if [ "$(names $FR/$m/HZAmassInde_fTest/EnvelopeResults.txt)" = "$(names $P/HZAmassInde_fTest/EnvelopeResults.txt)" ]; then
    hb "    $v: same members as production, not submitted"; continue
  fi
  mkdir -p $S/r2/$v/datacards $S/r2/$v/t2w $S/r2/$v/outputs
  prod_card=$F/Datacard/output_Datacard_leptons/${m}_pruned_datacard_leptons_xsec.txt
  card=$S/r2/$v/datacards/${m}_pruned_datacard_leptons_xsec.txt
  sed "s#$FR/$m/CMS-HGG_mva_13p6TeV_multipdf.root#$P/CMS-HGG_mva_13p6TeV_multipdf.root#g" $prod_card > $card
  [ "$(diff $prod_card $card | grep -c '^>')" -eq 2 ] || die "r2 $v: card differs in other than the two shape lines"
  text2workspace.py $card -o $S/r2/$v/t2w/${m}_Datacard_leptons.root -m 125.38 higgsMassRange=115,135 > $L/t2w_$v.log 2>&1 \
    || die "r2 $v: text2workspace"
  launch r2_$v "BIAS_HS_ROOT=$S/r2/$v/outputs ROOT_DATACARD_PATH=$S/r2/$v/t2w" "$m" "BIAS_HS_ROOT=$S/r2/$v/outputs ROOT_DATACARD_PATH=$S/r2/$v/t2w"
done

# ---------------------------------------------------------------- wait, truthout m_a = 25 once its card exists, merge
hb "[3] waiting (truthout m_a = 25 after chain_prod stage 4)"
declare -A MERGED
for _ in $(seq 1 400); do
  if [ -z "${CLUSTER[truthout25]:-}" ] && [ -f $DEP/logs/stage4.DONE ]; then
    launch truthout25 "BIAS_HS_ROOT=$TO ROOT_DATACARD_PATH=$EC/t2w FIT_DATACARD_PATH=$C/root_t2w SELECT_FUNCTION=Bern6" "25" \
      "BIAS_HS_ROOT=$TO ROOT_DATACARD_PATH=$EC/t2w SELECT_FUNCTION=Bern6"
  fi
  [ -f $DEP/logs/chain.FAILED ] && [ -z "${CLUSTER[truthout25]:-}" ] && { hb "    chain_prod failed before stage 4: truthout25 skipped"; CLUSTER[truthout25]=none; MERGED[truthout25]=1; }
  q=$(condor_q ${CLUSTER[@]/none/} -af ClusterId JobStatus 2>/dev/null) || { sleep 300; continue; }
  line=""
  for n in "${!CLUSTER[@]}"; do
    cl=${CLUSTER[$n]}; [ "$cl" = none ] && continue
    nq=$(awk -v c=$cl '$1==c' <<< "$q" | wc -l); nh=$(awk -v c=$cl '$1==c && $2==5' <<< "$q" | wc -l)
    line+=" $n=$nq${nh:+/h$nh}"
    if [ $nq -eq 0 ] && [ -z "${MERGED[$n]:-}" ]; then
      bad=$(condor_history $cl -af ExitCode 2>/dev/null | awk '$1!=0' | wc -l)
      hb "    $n ($cl) left the queue, non-zero exit codes: $bad; merging ${MMASS[$n]}"
      for m in ${MMASS[$n]}; do env ${MENV[$n]} bash $MERGE $m >> $L/merge_$n.log 2>&1 || hb "    merge $n m_a = $m rc!=0"; done
      MERGED[$n]=1
    fi
  done
  hb "   in queue:$line"
  done_all=1
  for n in "${!CLUSTER[@]}"; do [ -n "${MERGED[$n]:-}" ] || done_all=0; done
  [ -n "${CLUSTER[truthout25]:-}" ] && [ $done_all -eq 1 ] && break
  sleep 1800
done

python3 $S/summarize.py > $S/summary.txt 2>&1 || hb "    summarize.py rc!=0"
cat $S/summary.txt >> $L/heartbeat.txt
hb "SIDE DONE"
date '+%F %T' > $L/side.DONE
echo "$(date '+%F %H:%M') DONE  HZa bias side studies (truthout, ma29, r2)  -> $S/summary.txt" >> $REG
