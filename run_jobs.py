from __future__ import division
import math, os
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
  #
  root="/Users/pafonine/tmp17/qr-helix-3nir/perturbed/"
  rmsd_dirs = ["0.3/","0.6/","0.9/","1.2/","1.5/"]
  file_names = ["0.pdb","1.pdb","2.pdb","3.pdb","4.pdb","5.pdb","6.pdb","7.pdb",
                "8.pdb","9.pdb"]
  base_names = ["0","1","2","3","4","5","6","7",
                "8","9"]
  for rmsd_dir in rmsd_dirs:
    print rmsd_dir
    os.makedirs(rmsd_dir)
    for bn in base_names:
      file_name = root+rmsd_dir+bn+".pdb"
      if(not os.path.exists(file_name)): assert 0
      print file_name, bn
      cmd_full = cmd%(file_name,bn)
      print cmd_full
      easy_run.call(cmd_full)
      easy_run.call("mv pdb/%s_refined.pdb %s/%s.pdb"%(bn,rmsd_dir,bn))
      easy_run.call("mv %s.log %s"%(bn,rmsd_dir))
      easy_run.call("rm -rf pdb")

if __name__ == "__main__":
  run()
