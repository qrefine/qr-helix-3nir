from __future__ import division
import math, os
import iotbx.pdb
from scitbx.array_family import flex
import mmtbx.model
from libtbx.utils import null_out

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

def run():
  rmsd_dirs = ["0.3/","0.6/","0.9/","1.2/","1.5/"]
  base_names = ["0","1","2","3","4","5","6","7","8","9"]
  # Identify i_seqs of O-H pairs involved into H bonds
  h_bonds_i_seqs = []
  h = iotbx.pdb.input(file_name="helix_3nir_6_19.pdb").construct_hierarchy()
  atoms = h.atoms()
  sites_cart_orig = atoms.extract_xyz()
  ref_dist = flex.double()
  for i, a_i in enumerate(list(atoms)):
    for j, a_j in enumerate(list(atoms)):
      if(i<j):
        if(a_i.name.strip() in ["H","O"] and a_j.name.strip() in ["O","H"] and
          a_i.name.strip() != a_j.name.strip()):
          d = dist(a_i.xyz, a_j.xyz)
          if(d<2.4 and d>1.5): # 2.2, 1.7
            #print a_i.name.strip(), a_j.name.strip(), a_j.i_seq, a_i.i_seq, \
            #  a_j.i_seq-a_i.i_seq, d
            h_bonds_i_seqs.append([a_i.i_seq, a_j.i_seq])
            ref_dist.append(d)
  print "Total number of H bonds in one file:", len(h_bonds_i_seqs)
  print "Reference H-bond distances range (min,max,mean): %8.3f %8.3f %8.3f"%\
    ref_dist.min_max_mean().as_tuple()
  print "Hbond analysis"
  print "model     min     max     mean   %conserved"
  for sub_root in ["./perturbed/","./cctbx_opt/"]:
    print sub_root
    for rmsd_dir in rmsd_dirs:
      h_bonds = flex.double()
      cntr = 0
      for fn in base_names:
        file_name = sub_root+rmsd_dir+fn+".pdb"
        if(not os.path.exists(file_name)): assert 0
        model = get_model(file_name)
        g = model.geometry_statistics().result()
        sites_cart = model.get_sites_cart()
        #print file_name, bond_rmsd_mean, flex.mean(flex.sqrt((sites_cart_orig - sites_cart).dot()))
        for pair in h_bonds_i_seqs:
          d = dist(sites_cart[pair[0]], sites_cart[pair[1]])
          h_bonds.append(d)
        cntr+=1
      sel  = h_bonds < 2.3
      sel &= h_bonds > 1.6
      if(h_bonds.size()>0):
        print rmsd_dir, "%8.3f %8.3f %8.3f"%h_bonds.min_max_mean().as_tuple(), \
          "%7.2f"%(sel.count(True)*100./(len(h_bonds_i_seqs)*cntr))
      else:
        print rmsd_dir, "N/A"

if __name__ == "__main__":
  run()
