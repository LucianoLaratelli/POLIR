#waterEquil

You will need the [pgf community compiler](https://www.pgroup.com/products/community.htm)
to compile this code. Make sure to add `/opt/pgi/<your_os_and_architecture>/<pgi_version>/bin/`
to your $PATH. Once done, a compiled `water.equil` is a simple `make` away!

The file `inpEquil` is an example input for this program. 
Some quirks:
  * note that you will need to update the number of molecules by hand in `parameter.i` (the `nmol` line)

Description of the input file from Gungor:
 * 1st row: INPUT FILES---coordinates, velocities, dimensions of the box
 * 2nd row: OUTPUT FILES---final coordinates, final velocities
                        temperature, potential energy, kinetic energy, total energy
                        dipoles, md_coordinates, md_velocities
 * 3rd row: PARAMETERS---number of steps, first time step, polarizability, c_16, c_14, c_12, c_6, mass

The file `cvec.dat` has the coordinates of the simulation box.
`x.xyz` has the starting geometry of your system.
