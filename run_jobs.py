from __future__ import division
import math, os
import time
import iotbx.pdb
from scitbx.array_family import flex
import mmtbx.model
from libtbx.utils import null_out
from libtbx import easy_run

def dist(site1, site2):
  return math.sqrt(
    (site1[0]-site2[0])**2 +
    (site1[1]-site2[1])**2 +
    (site1[2]-site2[2])**2)

def get_model(file_name):
  params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  params.pdb_interpretation.use_neutron_distances = True
  model = mmtbx.model.manager(
    model_input               = iotbx.pdb.input(file_name=file_name),
    build_grm                 = True,
    log                       = null_out(),
    pdb_interpretation_params = params)
  return model
  
cmd_cctbx=" ".join([
    "qr.refine",
    "%s",
    "restraints=cctbx",
    "mode=opt",
    "stpmax=0.2",
    "gradient_only=true", # no problem with running up to upper bound.
    "clustering=false",
    "use_convergence_test=true",
    "rmsd_tolerance=0.001",
    "> %s.log"
    ])

def run(cmd = cmd_cctbx):
  results_prefix = "./cctbx_opt"
  perturbed_prefix = "./perturbed/"
  rmsd_dirs = ["0.3/","0.6/","0.9/","1.2/","1.5/"]
  base_names = ["0","1","2","3","4","5","6","7","8","9"]
  easy_run.call("rm -rf cctbx_opt")
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
