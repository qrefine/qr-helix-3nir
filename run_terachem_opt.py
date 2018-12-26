# At this moment, it is more convenitent to have independent job script for every method.
from __future__ import division
import math, os
import time
import iotbx.pdb
from scitbx.array_family import flex
import mmtbx.model
from libtbx.utils import null_out
from libtbx import easy_run

cmd_terachem=" ".join([
    "qr.refine",
    "%s",
    "restraints=qm",
    "engine_name=terachem",
    "mode=opt",
    "stpmax=0.2",
    "gradient_only=true", # no problem with running up to upper bound.
    "clustering=false",
    "use_convergence_test=true",
    "rmsd_tolerance=0.001",
    "> %s.log"
    ])

def run(cmd = cmd_terachem):
  results_prefix = "./terachem_opt"
  perturbed_prefix = "./perturbed/"
  rmsd_dirs = ["0.3/","0.6/","0.9/","1.2/","1.5/"]
  base_names = ["0","1","2","3","4","5","6","7","8","9"]
  easy_run.call("rm -rf terachem_opt")
  os.makedirs(results_prefix)
  os.chdir(results_prefix)
  for rmsd_dir in rmsd_dirs:
    print rmsd_dir
    os.makedirs(rmsd_dir)
    for bn in base_names:
      file_name = "../"+perturbed_prefix+rmsd_dir+bn+".pdb"
      if(not os.path.exists(file_name)): assert 0
      cmd_full = cmd%(file_name,bn)
      print "running command:\n%s"%(cmd_full)
      easy_run.call(cmd_full)
      easy_run.call("mv pdb/%s_refined.pdb %s/%s.pdb"%(bn,rmsd_dir,bn))
      easy_run.call("mv %s.log %s"%(bn,rmsd_dir))
      easy_run.call("rm -rf pdb")

if __name__ == "__main__":
  t0 = time.time()
  run()
  print "Time: %6.4f s" % (time.time() - t0) 
