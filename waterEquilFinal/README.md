# POLIR 128-water example

## BUILD
pgf90 -fast ls4equil.f90 r4.f ps_intra_polir_b.f90 -o water.equil

## Run
./water.equil < inpEquil

## Outputs
pos_out123.xyz -> final coordinates
vel123.out -> final velocities
temp.out -> temperature for each step
pe.out -> potential energy for each step
ke.out -> kinetic energy for each step
tot.oiut -> total energy for each step
md_traj.xyz -> XYZ trajectory file for full simulation
md_vel.out -> velocity trajectory

