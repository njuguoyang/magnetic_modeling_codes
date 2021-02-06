
Set paramters, read 3D magnetic field data, and save VTK files 
with compute_q.f90. 
The codes might need to be modified to read
 your own 3D magnetic field data.


Complie the FORTRAN codes using gfortran in Linux:


gfortran -ffree-form -ffree-line-length-none compute_q.f90 -o qsl



Run the code in Linux:


./qsl
