# Bias side studies

Run by `chain_side.sh` next to the production envelope revert
(`HiggsZaAna/studies/l3_review_20260910/deploy_envrevert_20261023/chain_prod.sh`).
All condor jobs use the production `runjob_bias_highstat.sh` (10 chunks x 1000 toys, production seeds, so the
truth functions shared with production are fitted on the same toys). Numbers: `summary.txt` (`summarize.py`).

| study | question | setup |
|---|---|---|
| `truthout/` | bias when the truth is a function the envelope no longer contains | toys from Bern6 (`envcheck_20261022/t2w`), fitted with the production workspace (`FIT_DATACARD_PATH`), m_a = 20 and 25 |
| `ma29/fallback` | does recovering the failed Bern5 fits at m_a = 29 change the bias | production options + `--cminFallbackAlgo Minuit2,1:1.0` + `Minuit2,Simplex,0:1.0` (`BIAS_FIT_EXTRA_OPTIONS`) |
| `ma29/fbrobust` | same, with the interval from `--robustFit 1` | fallback + `--robustFit 1`, one job per (chunk, truth function) |
| `r2/v11full` | is the m_a = 11 low-order restriction still needed | fTest without it (`FULL_ORDER_MASSES=11`, envcheck build) |
| `r2/v11b3` | is the m_a = 11 Bern2 cap still needed | low-order envelope with Bernstein up to 3 (`LOWORDER_BERN_EXCL=11:4`) |
| `r2/v27b6` | is the removal of Bernstein at m_a = 27 still needed | `BERN_ORDER_OVERRIDE=27:6` |

A r2 variant whose envelope has the same members as production is not submitted.
The job switches (`SELECT_FUNCTION`, `FIT_DATACARD_PATH`, `BIAS_FIT_EXTRA_OPTIONS`, `BIAS_HS_ROOT`, `ROOT_DATACARD_PATH`)
are unset in production, where the scripts behave as before.
