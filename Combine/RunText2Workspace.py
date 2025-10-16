import os, glob, sys
from optparse import OptionParser
from models import models

print(" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ HGG T2W RUN II ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ")

def get_options():
  parser = OptionParser()

  parser.add_option('--mA', dest='mA', default=5, type='int', help="ALP mass") # PZ
  parser.add_option('--year', dest='year', default='2022preEE', help="Comma separated list of years")

  parser.add_option('--mode', dest='mode', default='mu_inclusive', help="Physics Model (specified in models.py)")
  parser.add_option('--ext',dest='ext', default="", help='In case running over datacard with extension')
  parser.add_option('--common_opts',dest='common_opts', default="-m 125.38 higgsMassRange=115,135", help='Common options')
  parser.add_option('--batch', dest='batch', default='local', help="Batch system [SGE,IC,condor]")
  parser.add_option('--queue', dest='queue', default='workday', help="Condor queue")
  parser.add_option('--ncpus', dest='ncpus', default=4, type='int', help="Number of cpus")
  parser.add_option('--dryRun', dest='dryRun', action="store_true", default=False, help="Only create submission files")
  return parser.parse_args()
(opt,args) = get_options()

def leave():
  print(" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ HGG T2W RUN II (END) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ")
  exit(1)

def run(cmd):
  print("%s\n\n"%cmd)
  os.system(cmd)

if opt.mode not in models: 
  print(" --> [ERROR] opt.mode (%s) is not specified in models.py. Leaving..."%opt.mode)
  leave()

extStr = ""
print(" --> Running text2workspace for model: %s"%opt.mode)
print(" --> Input: /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Datacard/output_Datacard_leptons/%s_pruned_datacard_leptons.txt --> Output: ./root_t2w/%s_Datacard_leptons.root"%(opt.mA,opt.mA))

if not os.path.isdir(f"./t2w_jobs"): os.system(f"mkdir ./t2w_jobs")
if not os.path.isdir(f"./root_t2w"): os.system(f"mkdir ./root_t2w")

# Open submission file to write to
fsub = open("./t2w_jobs/%s_t2w.sh"%(opt.mA),"w")
fsub.write("#!/bin/bash\n\n")
fsub.write("cd %s\n\n"%os.environ['PWD'])
fsub.write("eval `scramv1 runtime -sh`\n\n")
# if not os.path.isdir("./output_Datacard%s"%extStr): os.system("mkdir ./output_Datacard%s"%extStr)
# fdataName = "./output_Datacard%s/%s_pruned_datacard_%s_%s.txt"%(extStr,opt.mA,opt.year,opt.channel)
# fsub.write("text2workspace.py Datacard%s.txt -o Datacard%s_%s.root %s %s"%(opt.ext,opt.ext,opt.mode,opt.common_opts,models[opt.mode]))
fsub.write("text2workspace.py /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Datacard/output_Datacard_leptons/%s_pruned_datacard_leptons.txt -o ./root_t2w/%s_Datacard_leptons.root %s %s"%(opt.mA,opt.mA,  opt.common_opts,models[opt.mode]))
fsub.close()

# /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Datacard/output_Datacard_leptons/15_pruned_datacard_leptons.txt

# Change permission for file
os.system("chmod 775 ./t2w_jobs/%s_t2w.sh"%(opt.mA))

# If using condor then also write submission file
if opt.batch == 'condor':
  f_cdr = open("./t2w_jobs/%s_t2w.sub"%(opt.mA),"w")
  f_cdr.write("executable          = %s/src/flashggFinalFit/Combine/t2w_jobs/%s_t2w.sh\n"%(os.environ['CMSSW_BASE'],opt.mA))
  f_cdr.write("output              = %s/src/flashggFinalFit/Combine/t2w_jobs/%s_t2w.sh.out\n"%(os.environ['CMSSW_BASE'],opt.mA))
  f_cdr.write("error               = %s/src/flashggFinalFit/Combine/t2w_jobs/%s_t2w.sh.err\n"%(os.environ['CMSSW_BASE'],opt.mA))
  f_cdr.write("log                 = %s/src/flashggFinalFit/Combine/t2w_jobs/%s_t2w.sh.log\n"%(os.environ['CMSSW_BASE'],opt.mA))
  f_cdr.write("+JobFlavour         = \"%s\"\n"%opt.queue)
  f_cdr.write("RequestCpus         = %g\n"%opt.ncpus)
  f_cdr.write("queue\n")
  f_cdr.close()

# Submit
if opt.batch == "condor": 
  if os.environ['PWD'].startswith("/eos"):
    subcmd = "condor_submit -spool ./t2w_jobs/%s_t2w.sub"%(opt.mA) 
  else:
    subcmd = "condor_submit ./t2w_jobs/%s_t2w.sub"%(opt.mA)
elif opt.batch == 'local': subcmd = "bash ./t2w_jobs/%s_t2w.sh"%(opt.mA)
else: subcmd = "qsub -q hep.q -l h_rt=6:0:0 -l h_vmem=24G ./t2w_jobs/%s_t2w.sh"%(opt.mA)
if opt.dryRun: print("[DRY RUN] %s"%subcmd)
else: run(subcmd)
