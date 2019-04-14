run gradient 
$multibasis 
Se lanl2dz_ecp
Cl lanl2dz_ecp
Cd lanl2dz_ecp
Zn lanl2dz_ecp
Mg lanl2dz_ecp 
$end 
basis 6-31g 
epsilon 78.39 
scf diis 
coordinates ./ase/tmp_ase.xyz 
gpumem 512 
gpus 4 
charge 2 
seed 1351351 
maxit 200 
threall 1e-12 
pcm cosmo 
watcheindiis no 
method rhf 
dftd yes 
scrdir /tmp/ase/scr 
end