# Envelope recheck with independent toys

The envelope changes of AN Sec 7.4 that were triggered by a bias just above the
0.2 threshold were all decided on high-statistics runs whose two splits per
chunk shared one seed (5000 distinct toys, uncertainty ~0.014 rather than the
quoted 0.010), and with earlier datacards. This study restores the pre-change
envelope at those mass points, with the current data and datacards, and reruns
the 10000-toy bias study with distinct split seeds.

| mA | restored | value that triggered the change |
|----|----------|---------------------------------|
| 2  | Bernstein (Bern2) | Bern2 -0.222 |
| 3  | Bern2 instead of Bern1 | Bern2 -0.211 |
| 4  | Bernstein (Bern4) | Bern4 -0.221 |
| 20 | Bern6 instead of Bern5 | Bern6 +0.208 |
| 24 | Laurent (Lau1) | Lau1 -0.232 |
| 25 | Bern6 instead of Bern5, and Laurent (Lau1) | Bern6 +0.212, Lau1 -0.215 |

Not rechecked: m_a = 11 (Bern3 -0.350) and 27 (Bern5 -0.299), far above threshold.

## Steps

1. `1_run_ftest.sh` — fTest with a patched copy of `fTest_ALP_turnOn` in `build/`
   (`BERN_ORDER_OVERRIDE`, `KEEP_LAURENT_MASSES`; unset = nominal). Output `bkg/<m>/`.
   `compare_envelopes.py` checks that the rebuilt envelope is the production one
   plus the restored member: all shared members have the same data 2NLL as in
   production (mA3 Pow3 and mA24 Exp1 differ in degenerate parameters only).
2. `2_make_cards.sh` — production xsec card with only the bkg multipdf path
   changed, then text2workspace with the production options. Output `t2w/`.
3. `3_submit.sh` — condor, 10 chunks x 1000 toys per mass, production job script
   and seeds (`BIAS_HS_ROOT`, `ROOT_DATACARD_PATH` overrides). Cluster 13052752.
   Shared truth functions are fitted on the same toys as the production run.
4. `watch.sh` — watchdog; merges with `merge_bias_highstat.sh` into
   `outputs/mA_<m>/merged/` and writes to `../PENDING_RESULTS.txt`.

The injected signal strength is the production q50 of each mass point.
