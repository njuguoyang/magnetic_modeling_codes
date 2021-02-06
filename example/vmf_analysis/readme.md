### Purpose and working flow

Load and process SDO/HMI vector magnetic field on the photosphere. Remove the 180 degree ambiguity. Correct the projection effect. Preprocess the vector magnetic field for nonlinear force-free field extrapolation. Prepare VTK files for visualizaiton.

(1) Download HMI vector magnetic field data in "./01data".

(2) If only the field, inclination, and azimuth data are available, one should go to "./02ambiguity". In 02ambiguity, read HMI data by 01plot_hmi. Prepare a file, "./02removing_amb/field000.dat", for the Minimum Energy (ME) disambiguiation. Complile and link an executable file, ambig, with the FORTRAN codes in ME0-1.1. Set the parameters in "./02removing_amb/par" for ambig. Then use "./02removing_amb/mybat.sh" to removing the ambiguity for multiple files. The output of ambig is azimuth***.data. Use "./02ambiguity/03amb_removed.pro" to prepare the ambiguity removed data, which are stored in "./02ambiguity/03amb_removed". Use "./02ambiguity/05projection.pro" to correct the projection effect.

Otherwise, if additional disambig data are available, one should go to "./02ambiguity_hmi".

(3) Use "./03preprocess/creb_lv3.pro" to preprocess the vector magnetic field for a nonlinear force-free field extrapolation.

Read [the readme file here](https://github.com/njuguoyang/magnetic_modeling_codes/blob/main/code/vector_magnetic_field_analysis) for more details on how to use these codes.
