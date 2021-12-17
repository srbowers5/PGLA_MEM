# PGLA_MEM
Initial data for PGLA in membrane

Download NAMD version NAMD_2.14_Linux-x86_64-multicore-CUDA   from site https://www.ks.uiuc.edu/Research/namd/

The namd configuration files are generated every 2 ps for each replica.
  The input to generate the files is pgl_mem_template.namd, RepTemp.dat, and ScaleFactor.dat.

pgl_mem_template.namd - template file to generate NAMD files.
RepTemp.dat - Temperatures of replicas
ScaleFactor.dat - Scaling factors used.
spt.pdb - Used to determine which atoms are part of solute (peptide) or solvent (membrane + water)
constraints.tcl - tcl script to set constraints for Phospate atoms.
boundaries.tcl - forces to keep peptide from crossing into opposite leaflet.
mem22_pep.psf - PSF file
par_all22_prot_cmap.inp - NAMD force field file.
par_all36_lipid.prm - NAMD force field file.

initSTR1 - initial structures, velosity, and xsc files for 12 replicas for trajectory 1.
initSTR2 - initial structures, velosity, and xsc files for 12 replicas for trajectory 2.
initSTR3 - initial structures, velosity, and xsc files for 12 replicas for trajectory 3.
initNAMD2 - first namd configuration files generated for trajectory 2.
        The files are the same for trajectory 1 and 3, except the input and
        output names are changed to reflect the trajectory.
