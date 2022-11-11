### Purpose and working flow with codes under this path

This path contains the user files for the potential and nonlinear force-free field extrapolation in the Cartesian coordinates with a uniform grid. 

The files have to be used in combination with MPI-AMRVAC version 2.2+ (commit fa03a48), which can be found via [this link](https://github.com/amrvac/amrvac/commit/fa03a48cfb53cee18b4f5db0c1ab1f28e122a6c3). Or, download the latest version and use the following Git command in a Git Repository:
> cd ~/codes    
> git clone https://github.com/amrvac/amrvac.git    
> cd ~/codes/amrvac    
> git checkout fa03a48

The user files and setups for both the potential field and nonlinear force-free field are provided in the files:
> ./potential    
> ./   

respectively. The extrapolations are performed in the following steps:

(1) Go to     
> ./potential/potential_boundary

Run the read_boundary.pro, which is an IDL procedure. It prepares the boundary conditions 
for the potential field. Note that allboundaries.dat contains the vector magnetic field 
on the photosphere observed by SDO/HMI. It has been processed following the steps detailed 
in [this paper](https://ui.adsabs.harvard.edu/abs/2017ScChD..60.1408G), and has been processed 
by the codes under [this path](https://github.com/njuguoyang/magnetic_modeling_codes/tree/main/example/vmf_analysis).

(2) Go to     
> ./potential

Setup, make, and run amrvac to get the potential field.

(3) Go to     
> ./nlfff_boundary

Run read_boundary.pro to prepare the boundary condition for the nonlinear force-free 
field extrapolation.

(4) Setup, make, and run amrvac to get the nonlinear force-free field in 
> ./

Note: refer to [amrvac.org](http://amrvac.org) for more information on how to setup, make, and run amrvac.
