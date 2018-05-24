waterEquil
=========

You will need the [PGI community compiler](https://www.pgroup.com/products/community.htm)
to compile this code. Make sure to add
`/opt/pgi/<your_os_and_architecture>/<pgi_version>/bin/`
to your $PATH. Once done, a compiled `water.equil` is a simple `make polir` away!

This directory also contains `get_TCF.cpp`, which calculates the 
time correlation function (TCFs) for the system.

Lastly, `get_ave.cpp` will calculate the average energy of the system.

`make` will generate optimized binaries for both `get_TCF` and `get_ave`

The file `inpEquil` is an example input for this program. 

Using Water.Equil
-----------------
```
water.equil < inpEquil
```

Description of the input file (`inpEquil`) from Gungor:
 * 1st row: INPUT FILES---coordinates, velocities, dimensions of the box
 * 2nd row: OUTPUT FILES---final coordinates, final velocities
                        temperature, potential energy, kinetic energy, total energy
                        dipoles, md_coordinates, md_velocities
 * 3rd row: PARAMETERS---number of steps, first time step, polarizability, c_16, c_14, c_12, c_6, mass

The file `cvec.dat` has the coordinates of the simulation box.
`x.xyz` has the starting geometry of your system.

Some quirks:
  * note that you will need to update the number of molecules by hand in `parameter.i` (the `nmol` line)


Using TCF Function
------------------
I apologize in advance for the way this is hacked together. (see the [License](https://github.com/LucianoLaratelli/POLIR/blob/master/waterEquil/CRAPL-LICENSE.txt))
```
python process_data.py dipole.out
```
Ten thousand steps produces a nice correlation function for a water dimer.
This program will also output a plot of the Fourier Transform of the TCF--hope you find it handy!


Known Errors
------------
1. Running the code on a water dimer does not move the hydrogens on the first water molecule.

TODO
----
1. Reduce known errors (above) to 0.
2. Allow for nmol to be determined dynamically from input file instead of being set as a parameter
3. Confirm that energy calculations are correct for the water dimer, namely:
  * charge-charge
  * charge-dipole
  * dipole-dipole
  * van Der Waals energy
