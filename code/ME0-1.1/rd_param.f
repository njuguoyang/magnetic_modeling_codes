c***********************************************************************
      subroutine rd_param(filename,bthresh,npad,nap)
c***********************************************************************
c     Read in some parameters which set options for the code.          *
c***********************************************************************
      include "include_file.sizes.f"
      include "include_file.flags.f"
      include "include_file.fits.f"
      include "include_file.lambda.f"
      include "include_file.verb.f"
      include "include_file.seed.f"
      include "include_file.anneal.f"
c-----------------------------------------------------------------------
      integer npad,nap,iseed
      real bthresh
      character*8 radkey
      character*90 filename
c-----------------------------------------------------------------------
      open(1,file='par',status='old')

c --> filename for the input file.
      read(1,'(a90)') filename

c --> irflag determines format of input: 1: FITS, 0: arrays
c --> ibflag determines whether field is specified as
c -->    0: two components and an angle
c -->    1: magnitude and two angles
c --> iaflag determines the direction of zero azimuth as
c -->    0: west
c -->    1: north (UHIM standard)
      read(1,*) irflag,ibflag,iaflag

c --> npad determines whether to surround with zeros.
c --> nap sets the number of pixels over which to apodize the edge.
      read(1,*) npad,nap

      if(nap.gt.npad) nap=npad

c --> bthresh determines the level below which to use acute angle.
      read(1,*) bthresh

c --> iaunit specifies whether azimuth is in radians or degrees. 
c --> ipunit specifies whether pointing is in radians or degrees. 
c --> incflag specifies whether inclination of 0 is vertical or
c -->         horizontal.
      read(1,*) iaunit,ipunit,incflag

c --> The following control the annealing schedule and similar.

      read(1,*) iseed,iverb,neq
      read(1,*) lambda,tfac0,tfactr

c --> seed for random number generator

      seed=iseed
      if (seed.le.0) seed=123456

c --> check other input parameters

      if (lambda.lt.0) then
         print*,'read_param: lambda<0'
         stop
      endif
      if (tfac0.lt.0) then
         print*,'read_param: tfac0<0'
         stop
      endif
      if (neq.le.0) then
         print*,'read_param: neq.le.0'
         stop
      endif
      if ((tfactr.lt.0).or.(tfactr.gt.1)) then
         print*,'read_param: tfactr problem'
         stop
      endif

c --> If the data are read from a FITS file, need some additional 
c --> information about how to extract the necessary data.

      if(irflag.eq.1) then

c --> Array elements in the FITS image containing the longitudinal and
c --> transverse field components, and the azimuth.
         read(1,*) iblong,ibtrans,ibazim
 
c --> Keyword names for longitude, latitude, p and b angles, and pixel
c --> size.
         read(1,'(a8)') cmdkey
         read(1,'(a8)') latkey
         read(1,'(a8)') pkey
         read(1,'(a8)') bkey
         read(1,'(a8)') radkey
         read(1,'(a8)') pxkey
         read(1,'(a8)') pykey

      else if(irflag.eq.2) then
c --> If the FITS file comes from the Hinode SOT/SP pipeline: 
c --> number of extension to FITS file containing field strength,
c --> inclination and azimuth, and the continuum intensity.
         read(1,*) jfs,jfi,jfa,jsf,jci
      endif

      return
      end
