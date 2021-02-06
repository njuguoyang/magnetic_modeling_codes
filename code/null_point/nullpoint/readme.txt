The program is used to calculated the location of the magnetic null and designed 
to use only one thread.
The main program is null.pro.


Program cube_index.pro use the Poincare theory to scan each of cube with 8 point 
in the computed domain for the possible location of a null.
 You can use this program 
as the example 'cal.pro' or just change the input parameter in it and use it as 
tapping '@cal.pro' in SSWIDL.


The time costs by the testfield.sav is about 30 min.
 The testfield.pro is used to 
output the testfield.sav and testfield.vtk. The result of testfield.sav is seved as 
an example. One could visulize it by testfield.vtk in Paraview.



The result contains:

null_position[i,*,*]:       the pixel location of the null.

null_position_error[i,*,*]: the associated error of the location of the null.

null_magnitude[i,*]:        the magnitude of the magnetic field at the null point which can be used to confirm the null.

null_matrix[i,*,*]:         the Jacobian matrix of the magnetic field at the null.

null_eigen_value[i,*]:      the eigen value of the Jacobian matrix.

null_eigen_vector[i,*,*]:   the eigen vector of the Jacobian matrix.
