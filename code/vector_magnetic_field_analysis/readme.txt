These IDL procedures are provided to analyze vector magnetic field. 
Please refer to the following papers for the description of the 
analyzing steps:
https://ui.adsabs.harvard.edu/abs/2013A%26A...555A..19G
https://ui.adsabs.harvard.edu/abs/2013ApJ...779..157G
https://ui.adsabs.harvard.edu/abs/2017ScChD..60.1408G

The procedures under the directory "./vector_magnetic_field_analysis_example", 
are written by Yang Guo (Nanjing Univeristy). They are examples for reading HMI data, removing 
the 180 degree ambiguity of the transverse components of the vector 
magnetic field, correcting the projection effect, preprocessing for 
nonlinear force-free field extrapolation, coaligning two magnetic field maps, 
computing the velocity and vertical electric currents. The precedures should be rewritten 
for you own purpose. These procedures may call the procedures in the following 
directories:
./ambiguity
./DAVE-DAVE4VM
./fieldline
./NPFC
./null_point
./prepro_wie
./QSL
./twist
which should be put under the search paths of IDL.

Fow example, use .cshrc in Linux to set the search paths of IDL.
add the following command in .cshrc
    setenv IDL_STARTUP /home/my_path/.idl_startup.pro
create a new file, .idl_starup.pro under the path /home/my_path
add the following command in .idl_startup.pro
    !path=expand_path('+/home/my_path/vector_magnetic_field_analysis')+':'+!path
after sourcing the .cshrc file, the path containing the vector magnetic field analysis
code will be included in the search paths of IDL.

In "./vector_magnetic_field_analysis_example", only procedures and working directories 
are provided, while data and executables are omitted. An exmaple with SDO/HMI data
is provided in a separate but larger directory "../../example/example_vector_magnetic_field_20150827.zip".

An example of the working flow is as follows:
(1) Download HMI vector magnetic field data in "./vector_magnetic_field_analysis/01data".

(2) If only the field, inclination, and azimuth data are available, one should go to
"./vector_magnetic_field_analysis/02ambiguity". In 02ambiguity, read HMI data by 
01plot_hmi. Prepare a file, "./02removing_amb/field000.dat", for the Minimum Energy (ME) disambiguiation. 
Complile and link an executable file, ambig, with the FORTRAN codes in ME0-1.1.
Set the parameters in "./02removing_amb/par" for ambig. Then use "./02removing_amb/mybat.sh"
to removing the ambiguity for multiple files. The output of ambig is azimuth***.data.
Use "./vector_magnetic_field_analysis_example/02ambiguity/03amb_removed.pro" to prepare the ambiguity 
removed data, which are stored in "./vector_magnetic_field_analysis_example/02ambiguity/03amb_removed".
Use "./vector_magnetic_field_analysis_example/02ambiguity/05projection_modified_version.pro" to correct
the projection effect.

Otherwise, if additional disambig data are available, one should go to 
"./vector_magnetic_field_analysis_example/02ambiguity_hmi".

(3) Use "./vector_magnetic_field_analysis_example/03preprocess/creb_lv3_projection_m.pro" to preprocess 
the vector magnetic field for a nonlinear force-free field extrapolation.

(4) Coalign a time series of vector magnetic field using procedures in 
"./vector_magnetic_field_analysis_example/04coalign".

(5) Compute the velocity using procedures in "./vector_magnetic_field_analysis_example/05velocity".

(6) Compute the electric currents using procedures in "./vector_magnetic_field_analysis_example/06currents".

(7) Write vtk for visulization with Paraview using procedures in 
"./vector_magnetic_field_analysis_example/07_write_vtk".

(8) Extrapolate the potential and nonlinear force-free field using codes in a separate zip file "MPI-AMRVAC_1.0.zip"
or "MPI-AMRVAC_2.0.zip", which are version 1.0 and 2.0 of MPI-AMRVAC, respectively. Some examples of the user files 
have been provided in the zip files initiated with "example".

(9) Compute magnetic null points with procedures in "./null_point".

(10) Compute quasi-separatrix layers (QSLs) with FORTRAN codes in the directory "./QSL".

(11) Compute twist with procedures in "./twist".


=====================================
More information for the codes
=====================================
The procedures under the directory "./DAVE-DAVE4VM" are written by P. W. Schuck. 
Please refer to the following papers for more information:
https://ui.adsabs.harvard.edu/abs/2005ApJ...632L..53S
https://ui.adsabs.harvard.edu/abs/2006ApJ...646.1358S
https://ui.adsabs.harvard.edu/abs/2008ApJ...683.1134S

The procedures under the directory "./fieldline" are written by Yang Guo.
It is used to compute, store, and visulize magnetic field lines in IDL.
 
The procedures under the directory "./NPFC" are written by Manolis K. Georgoulis. 
Please refer to the following paper for more information:
https://ui.adsabs.harvard.edu/abs/2005ApJ...629L..69G

The procedures under the directory "./null_point" are written 
by Kai Yang, Ze Zhong, and Yang Guo. It is used to locate magnetic null points and 
to visualize magnetic field lines in the vicinity of a magnetic null point.

The prodedures in "./prepro_wie" are writeen by T. Wiegelmann. Refer to this paper 
for more information:
https://ui.adsabs.harvard.edu/abs/2006SoPh..233..215W

The procedures under the directory "./twist" are written by Yang Guo
It is used to compute the twist of two field lines. Please refer to the following papers 
for more information:
https://ui.adsabs.harvard.edu/abs/2010ApJ...725L..38G
https://ui.adsabs.harvard.edu/abs/2017ApJ...840...40G

The procedures in the ME0-1.1 directory are written by Thomas R. Metcalf and 
further improved by K. D. Leka. Please refer to the following papers for more 
information:
https://ui.adsabs.harvard.edu/abs/1994SoPh..155..235M
https://ui.adsabs.harvard.edu/abs/2006SoPh..237..267M
https://ui.adsabs.harvard.edu/abs/2009SoPh..260...83L

The procedures in QSL are written by Kai Yang and Yang Guo.
Read these papers for more information:
https://ui.adsabs.harvard.edu/abs/2013ApJ...779..157G
https://ui.adsabs.harvard.edu/abs/2015ApJ...806..171Y
A updated version written by Kai Yang can be found in this link:
https://github.com/Kai-E-Yang/QSL
