#!/bin/bash

default_bias_fit_combine_options() {
  cat <<'EOF'
--cminDefaultMinimizerType Minuit2 --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.01 --cminPreScan --cminFallbackAlgo Minuit2,0:0.1 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2 --freezeParameters MH --rMin -5 --rMax 5 --setParameterRanges 'r=-5,5:Exp_turnon_p1=95,125:Exp_turnon_p3=95,125:Exp_turnon_p5=95,125:Exp_width_p1=0.1,50:Exp_width_p3=0.1,50:Exp_width_p5=0.1,50:Exp_sigma_p1=0.05,20:Exp_sigma_p3=0.05,20:Exp_sigma_p5=0.05,20:Pow_turnon_p1=95,125:Pow_turnon_p3=95,125:Pow_turnon_p5=95,125:Pow_width_p1=0.1,50:Pow_width_p3=0.1,50:Pow_width_p5=0.1,50:Pow_sigma_p1=0.05,20:Pow_sigma_p3=0.05,20:Pow_sigma_p5=0.05,20:Lau_turnon_p1=95,125:Lau_turnon_p2=95,125:Lau_turnon_p3=95,125:Lau_turnon_p4=95,125:Lau_width_p1=0.1,50:Lau_width_p2=0.1,50:Lau_width_p3=0.1,50:Lau_width_p4=0.1,50:Lau_sigma_p1=0.05,20:Lau_sigma_p2=0.05,20:Lau_sigma_p3=0.05,20:Lau_sigma_p4=0.05,20:Bern_gsigma=0.05,20:Bern_step=90,130:Bern_stepWidth=0.1,50'
EOF
}

resolve_bias_fit_combine_options() {
  if [[ -n "${BIAS_FIT_COMBINE_OPTIONS:-}" ]]; then
    printf '%s' "${BIAS_FIT_COMBINE_OPTIONS}"
  else
    default_bias_fit_combine_options
  fi
}
