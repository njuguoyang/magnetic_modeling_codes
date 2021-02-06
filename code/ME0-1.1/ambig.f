      program ambig
c***********************************************************************
c   Resolve the 180 degree ambiguity in the transverse component of the*
c   magnetic field using the Minimum Energy approach first described by*
c   Metcalf (1994, Solar Phys., 155, 235). The code reads in data from * 
c   either a FITS or text file, calculates derivatives for a potential *
c   field, then uses simulated annealing to minimize the sum of the    *
c   absolute value of the divergence plus the absolute value of the    *
c   vertical current density.                                          *
c   For pixels with a transverse field strength below a given threshold*
c   the ambiguity is resolved using an acute angle to nearest neighbors*
c   methods, based on that described in Canfield et al. (1993, ApJ,    *
c   411, 362).                                                         *
c   The potential field is calculated directly from the line of sight  *
c   component of the field on a uniform grid in image coordinates. This*
c   avoids having to interpolate to and from a uniform grid in         *
c   helioplanar coordinates in order to use FFTs. In order to minimize *
c   edge effects, the code surrounds the observed field with an area of*
c   zero field, and uses a bi-cubic spline to smoothly interpolate from*
c   observed field to zero field.                                      *
c***********************************************************************
      include "include_file.sizes.f"
      include "include_file.flags.f"
      include "include_file.field.f"
      include "include_file.potential.f"
      include "include_file.verb.f"
c-----------------------------------------------------------------------
      integer ixmin,ixmax,iymin,iymax,nxp,nyp,npad,nap
      real bthresh
      real,dimension(2*nxmax,2*nymax) :: blpad
      character*90 filename
c-----------------------------------------------------------------------
c --> Read in some parameters which set options for the code. 

      call rd_param(filename,bthresh,npad,nap)
      if(iverb.eq.1) write(*,*) 'parameters read'

c --> Read in magnetic field and pointing information. 
c --> irflag determines whether it is from a FITS file or not, and 
c --> whether the FITS file is from the Hinode SOT/SP pipeline. 

      if(irflag.eq.0) then
         call rd_field(filename,npad)
      else if(irflag.eq.1) then
         call rd_fits(filename,npad)
      else if(irflag.eq.2) then
         call rd_sotsp(filename,npad)
      endif
      if(iverb.eq.1) write(*,*) 'field read'

c --> Construct coordinate transformation matrices. 

      call transform()
      if(iverb.eq.1) write(*,*) 'transform calculated'

c --> Pad the input array with zeros.
c --> This reduces the impact of the periodic boundary conditions. 

      call buffer(Bz,blpad,ixmin,ixmax,iymin,iymax,npad,nxp,nyp)
      if(iverb.eq.1) write(*,*) 'buffering done'

c --> Smooth the transition to zero field. 
c --> This reduces ringing if the field is non-zero at the boundaries.

      if(nap.ge.1) call smooth(blpad,ixmin,ixmax,iymin,iymax,nap)
      if(iverb.eq.1) write(*,*) 'apodizing done'

c --> Calculate potential field (and derivatives) using FFTs. 

      call potential(blpad,nxp,nyp,ixmin,ixmax,iymin,iymax)
      if(iverb.eq.1) write(*,*) 'potential field done'

c --> Determine box enclosing all above-threshold pixels. 

      call boxit(bthresh)

c --> Global minimization of "energy" to resolve ambiguity. 

      call setup
      call minimise_energy
      if(iverb.eq.1) write(*,*) 'optimization done'

c --> Perform acute angle to nearest neighbor ambiguity resolution on
c --> pixels below a given threshold. 

      call nacute5(bthresh)
      if(iverb.eq.1) write(*,*) 'smoothing done'

c --> If the data came from a FITS file, update it; otherwise, write
c --> out the azimuth array. 

      if(irflag.eq.0) then
         call write_field()
      else if(irflag.eq.1) then
         call update_fits(filename)
      else if(irflag.eq.2) then
         call write_field()
      endif
      if(iverb.eq.1) write(*,*) 'results written'

      stop '***END***'
      end
