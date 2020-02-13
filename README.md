# Molecular Dynamics 

This program simulates molecular interactions of monoatomic molecules and particles with many molecules in a 2D region of calculus. 
The development of this code was part of a graduante reseach project of the author's Master Degree at the Faculty of Mechanical Engineering at University of Campinas. 

## Compilation

The Fortran code can be compiled using mpiifort (Intel Fortran Compiler) or mpifort (GCC Fortran). You may change the make file to choose the compiler that you wish. 
The main file is "lennard.f90" 

## Requirements
GNU Fortran 
```
sudo apt-get install gfortran # Ubuntu Based distros
sudo zypper in gcc-fortran # OpenSUSE
```
Or [Intel Fortran Compiler](https://software.intel.com/en-us/fortran-compilers)

If using GNU Fortran, install also the OpenMPI package
```
sudo apt-get install openmpi-bin libopenmpi-dev # Ubuntu 
sudo zypper in openmpi openmpi-devel # OpenSUSE
```
For the Python code it is necessary the pyevtk package to convert the csv files to vtk. 
The package progressbar2 is optional. 

For the visualization, you can use Paraview or some program that reads VTK files. 


## Input and Output

The input needs a configuration file named settings.ini , which contains the physical properties and numerial parameters for the simulation, and files containing the positions of the molecules. Optionally a file containing the initial velocity of the particles can also be provided. 

The output contains the position and velocities of the molecules. Optionally the angle of rotation of particles, and forces related to intermolecular interactions can also be extracted.

## Postprocessing

The program automatically converts the .csv results to VTK, so one may use this result to postprocess on paraview. 
Some python scripts are provided to calculate concentration, distribution entropy, viscosity and  other things. The theory behind these program is summarized at the Master Thesis that will be published soon. 

## References

Some of the code on this repository was based on code published by 
Jason R. Blevins ACM Fortran Forum 28(3), 2–7, 2009.
https://github.com/jannisteunissen/config_fortran
Numerical Simulation in Molecular Dynamics: Numerics, Algorithms, Parallelization, Applications, by Michael Griebel (Author), Stephan Knapek (Contributor), Gerhard Zumbusch (Contributor) 

## Acknowledgments 

The author to acknowledge the support of this work by FAPESP-Brazil under the grant number 2017/05643-8 and FAEPEX Unicamp. 


