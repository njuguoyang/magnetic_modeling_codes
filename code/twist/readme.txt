To calculate the twist between a field line and the axis.

1. Use "findaxis.pro" to determine the axis. The point, through which the axis
field line passes, is saved in "sample_start_point.sav". This procedure is
not universal. Different procedures for different cases should be developed
for specific purposes.

2. The procedure "twist.pro" calculates the twist between one field line and the 
axis (see Berger, A. and Prior, C. 2006, J. Phys. A: Math. Gen. 39, 8321).

3. The procedure "faxis.pro" calcualtes selected field lines, the twists between
those field lines and the axis, and plots the figures.

In order to keep the speed and increase the accuracy in calculating the twists,
the step size for the axis integration is one order larger than that for other 
field lines. The curve data with larger step size is in the directory "curve_step_
size0.02". So, the commands should be as follows when calculating the twist:

IDL> axis_filename='curve_step_size0.02/curve_info0007.sav'
IDL> line_filename='curve_info0002.sav'
IDL> ntwist=twist(axis_filename,line_filename)
