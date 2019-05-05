# qr-helix-3nir
This is a helix extracted from Crystal structure of small protein crambin at 0.48 A resolution. Except for its high resolution feature, its sidechain is complete and diverse.

## Analysis of optimization results

### H bond analysis
   Total number of H bonds: 9.

   Reference H-bond distances range (min,max,mean):    2.005    2.325    2.139
   
   "use_hydrogens=False" is used for bond rmsn analysis.
#### perturbed:
model   | min   |  max   |  mean  |recovered |bond_rmsd| rama_favored |clashscore
:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:
0.3     |1.921  |2.639   |2.237    |64.44    |0.0171    |100.00 |0.00
0.6     |1.937  |3.279   |2.584    |22.22    |0.0186    |95.00  |4.93
0.9     |1.909  |4.178   |2.712    |28.89    |0.0190    |92.50  |0.00
1.2     |2.189  |5.793   |3.702    |3.33     |0.0185    |76.67  |0.45
1.5     |1.630  |6.461   |3.880    |5.56     |0.0189    |83.33  |3.59
#### cctbx_opt:
model   | min   |  max   |  mean  |recovered |bond_rmsd| rama_favored |clashscore
:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:
0.3    |1.988   |2.616   |2.257    |52.22     |0.0035    |100.00 |0.00
0.6    |2.094   |2.939   |2.441    |35.56     |0.0032    | 91.67 |0.00
0.9    |2.057   |4.825   |2.832    |23.33     |0.0033    | 85.83 |0.00
1.2    |2.079   |5.752   |3.472    |12.22     |0.0034    | 83.33 |0.90
1.5    |1.695   |6.724   |3.880    |11.11     |0.0035    | 82.50 |0.90

#### xtb_opt_water:
model   | min   |  max   |  mean  |recovered |bond_rmsd| rama_favored |clashscore
:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:
0.3    |1.811   |2.219   |1.936   |100.00     |0.0112     |91.67 |4.04
0.6    |1.795   |2.601   |2.006   | 88.89     |0.0113     |91.67 |4.48
0.9    |1.805   |2.722   |2.028   | 88.89     |0.0113     |91.67 |4.48
1.2    |1.769   |4.554   |2.420   | 63.33     |0.0116     |75.83 |1.79
1.5    |1.802   |4.397   |2.784   | 42.22     |0.0115     |71.67 |2.69
### xtb_opt_stpmax_0.3:
model   | min   |  max   |  mean  |recovered |bond_rmsd| rama_favored |clashscore
:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:
0.3/   |1.810   |2.186   |1.931   |100.00     |0.0094     |91.67  |4.42
0.6/   |1.798   |2.577   |2.005   | 88.89     |0.0094     |91.67  |4.87
0.9/   |1.803   |2.716   |2.030   | 88.89     |0.0095     |91.67  |6.19
1.2/   |1.766   |4.550   |2.364   | 67.78     |0.0096     |82.50  |6.64
1.5/   |1.774   |4.344   |2.711   | 44.44     |0.0096     |70.00  |0.88
### xtb_opt_stpmax_0.4:
model   | min   |  max   |  mean  |recovered |bond_rmsd| rama_favored |clashscore
:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:
0.3/   |1.811   |2.209   |1.934   |100.00     |0.0094     |91.67    |4.42
0.6/   |1.798   |2.574   |2.002   | 88.89     |0.0094     |91.67    |5.31
0.9/   |1.802   |2.716   |2.025   | 88.89     |0.0095     |91.67    |7.52
1.2/   |1.766   |4.318   |2.253   | 72.22     |0.0096     |88.33    |7.08
1.5/   |1.776   |4.133   |2.662   | 46.67     |0.0095     |72.50    |1.33
### xtb_opt_stpmax_0.5:
model   | min   |  max   |  mean  |recovered |bond_rmsd| rama_favored |clashscore
:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:
0.3/   |1.813   |2.193   |1.936   |100.00     |0.0094     |91.67    |4.42
0.6/   |1.795   |2.564   |1.995   | 88.89     |0.0094     |91.67    |8.85
0.9/   |1.797   |2.689   |2.013   | 88.89     |0.0095     |91.67    |8.41
1.2/   |1.767   |3.391   |2.181   | 75.56     |0.0095     |91.67    |7.52
1.5/   |1.753   |3.909   |2.503   | 56.67     |0.0095     |78.33    |1.77
### terachem_opt_water
model   | min   |  max   |  mean  |recovered |bond_rmsd| rama_favored |clashscore
:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:
0.3    |1.783    |2.024    |1.867    |100.00   |0.0150    |100.00 |0.00
0.6    |1.791    |2.309    |1.897    |97.78    |0.0150    |100.00 |0.00
0.9    |1.789    |2.318    |1.903    |93.33    |0.0150    |100.00 |0.00
1.2    |1.781    |2.314    |1.902    |98.89    |0.0149    |100.00 |0.90
1.5    |1.735    |3.604    |2.250    |67.78    |0.0149    |93.33  |0.00
### terachem_opt_final
model   | min   |  max   |  mean  |recovered |bond_rmsd| rama_favored |clashscore
:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:
1.5    |1.786    |2.665   |1.921    |95.56    |0.0149     |100.00  |1.35
### simulated map are generated by commands:
1. phenix.fmodel helix_3nir_6_19.pdb high_res=4 low_res=6
2. phenix.mtz2map helix_3nir_6_19.pdb.mtz include_fmodel=True
3. mv helix_3nir_6_19.pdb_fmodel.ccp4 map.ccp4
