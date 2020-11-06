### Purpose and working flow with codes under this path

01plot_hmi.pro is used to read SDO/HMI data and to remove the 180 degree ambiguity using 
"*disambig.fits" in "../01data". The most important variables in 01plot_hmi.pro, which 
needs to been determined by the user, are dx0 and dy0. They are used to determine the 
center coordinate of the solar disk. We ususlly measure the four extremes of the solar 
disk to determine these parameters. One may also use other methods to do that. Setting 
dx0=0 and dy0=0 to adopt the default coordinates provided by SDO/HMI.

05projection.pro is used to correct the projection effect. Important issues include the
file:
./05projection/cmd_lat.sav
which contains the central meridian and latitudes of the interesetd data. If this file
is not there, 05projection.pro will compute a new one. If it is there, the procedure 
will read it. So, remove this file if one wants to compute a new one. Another important
variable is Lcc, which is the longitude and latitude of the point at which the heliographic 
plane is tangent to. It is set to be the longitude and latitude of the center of the field 
of view. But sometimes, one may want to control Lcc to be a specific location. Then, try
to tune this variable.

Once we get the output file from 05projection.pro, go to "./05projection", run plot_submap_m.pro.
Here, one needs to set the x- and y-ranges for the field of view. One may also use the irange
keywords in sub_map.pro to control the pixel numbers more precisely. To fulfill the requirements
of procedures in "../03preprocess" and AMRVAC, the pixel numbers of the cropped region of
interest must obey:
nx = (2^n1 * 3^n2 * 5^n3 ... + 4)*(2^(level-1))
ny = similar to nx but unnecessary to be the same
n1 is integer geter than 0, n2, n3 ... are integer number that greater than or equivalent to 0,
and level=1, 2, 3 ...
(1) nx and ny needs to be multiples of 2^(level-1). This is the requirement of [creo_lv3.pro](https://github.com/njuguoyang/magnetic_modeling_codes/blob/main/example/example_vector_magnetic_field_20150827/03preprocess/creb_lv3.pro),
    level=1 means no rebin, keeping the original spatial resoltuion. 
    level=2 means 2*2 pixel rebin. 
    level=3 means 4*4 pixel rebin.
(2) 4 pixels are added to 2^n1 * 3^n2 * 5^n3 ... This is used to set values in the ghost cells
    for AMRVAC
(3) The remaining are multiples of small primes. This is used to decompose the physical domain
    to many blocks, and use multi CPUs to parallelize the code.
For example, if we select level=3, nx could be (2^8 + 4)*4 =1040, or, nx = (2^6 * 3 + 4)*4 = 784.
