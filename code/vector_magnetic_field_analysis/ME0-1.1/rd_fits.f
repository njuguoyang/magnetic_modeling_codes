c***********************************************************************
      subroutine rd_fits(filename,npad)
c***********************************************************************
c     Reads in magnetic field data plus pointing information from a    *
c     FITS file. Based on routines given in cookbook.f supplied with   *
c     FITSIO.                                                          *
c***********************************************************************
      include "include_file.sizes.f"
      include "include_file.point.f"
      include "include_file.field.f"
      include "include_file.fits.f"
      include "include_file.flags.f"
      include "include_file.pix_size.f"
c-----------------------------------------------------------------------
      integer i,j,i1d
      integer status,unit,readwrite,blocksize
      integer firstpix,group,nfound,npixels,narray
      real nullval,pi,phase
      real,dimension(nxmax,nymax) :: bl,bt,ba
      real,dimension(nxmax*nymax) :: bl1d,bt1d,ba1d
      logical anynull
      character*90 filename
      character*80 record,comment
      common/bobs/ bl,bt,ba
c-----------------------------------------------------------------------
      data pi/3.1415926535897932/
c-----------------------------------------------------------------------
c --> The status parameter must always be initialized.
      status=0

c --> Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)

c --> Open the FITS file, with read-only access.  The returned blocksize
c --> parameter is obsolete and should be ignored. 
      readwrite=0
      call ftopen(unit,filename,readwrite,blocksize,status)

c --> Determine the size of the image.
      call ftgknj(unit,'NAXIS',1,3,naxes,nfound,status)
      if(nfound.ne.3) stop '*** Unexpected number of axes found ***'

      nx=naxes(1)
      ny=naxes(2)

c --> Check that various array dimensions are consistent.
      if(nx.gt.nxmax.or.ny.gt.nymax) then
         write(*,*) 'Increase array dimensions in include_file.size.f'
         write(*,*) 'nxmax=',nxmax,' must be greater than nx=',nx
         write(*,*) 'nymax=',nymax,' must be greater than ny=',ny
         stop
      endif

      if(nx+2*npad.gt.2*nxmax.or.ny+2*npad.gt.2*nymax) then
         write(*,*) 'Increase array dimensions in include_file.size.f'
         write(*,*) 'or decrease padding in file par'
         write(*,*) '2 nxmax=',2*nxmax,' must be greater than nx+2npad='
     +               ,nx+2*npad
         write(*,*) '2 nymax=',2*nymax,' must be greater than ny+2npad='
     +               ,ny+2*npad
         stop
      endif

c --> Get the pointing information from the keywords. 
c --> b=solar b-angle
c --> p=solar p-angle

      call ftgkye(unit,cmdkey,theta,comment,status)
      call ftgkye(unit,latkey,phi,comment,status)
      call ftgkye(unit,pkey,p,comment,status)
      call ftgkye(unit,bkey,b,comment,status)

c --> Convert to radians, if necessary.
      if(ipunit.eq.1) then
         theta=theta*pi/180.
         phi=phi*pi/180.
         p=p*pi/180.
         b=b*pi/180.
      endif

      if(status.gt.0) call printerror(status)

c --> Get the pixel size from the keywords.

      call ftgkye(unit,pxkey,dxi,comment,status)
      call ftgkye(unit,pykey,dyi,comment,status)

c --> Initialize variables for reading field. 
      narray=naxes(1)*naxes(2)
      npixels=narray*naxes(3)
      group=1
      nullval=-999.

c --> Different steps required if field is specified as a magnitude plus
c --> two angles versus two components plus one angle. 

      phase=0.5*float(iaflag)*pi
      if(ibflag.eq.0) then

c --> Read in line of sight component.
         firstpix=1+narray*(iblong-1)
         call ftgpve(unit,group,firstpix,narray,nullval,
     &               bl1d,anynull,status)

         do 2 i=1,nx
            do 1 j=1,ny
               bl(i,j)=bl1d(i+(j-1)*nx)
    1       continue
    2    continue

c --> Read in transverse component.
         firstpix=1+narray*(ibtrans-1)
         call ftgpve(unit,group,firstpix,narray,nullval,
     &               bt1d,anynull,status)

         do 4 i=1,nx
            do 3 j=1,ny
               bt(i,j)=bt1d(i+(j-1)*nx)
    3       continue
    4    continue

c --> Read in azimuthal angle.
         firstpix=1+narray*(ibazim-1)
         call ftgpve(unit,group,firstpix,narray,nullval,
     &               ba1d,anynull,status)

c --> Convert to radians if necessary. 
         if(iaunit.eq.1) then
            dtor=pi/180.
            do 5 i=1,narray
               ba1d(i)=ba1d(i)*dtor
    5       continue
         endif

         do 7 i=1,nx
            do 6 j=1,ny
               ba(i,j)=ba1d(i+(j-1)*nx)
c --> Calculate Cartesian (image) components.
               Bx(i,j)=bt(i,j)*cos(ba(i,j)+phase)
               By(i,j)=bt(i,j)*sin(ba(i,j)+phase)
               Bz(i,j)=bl(i,j)
    6       continue
    7    continue

      else

c --> Read in magnitude.
         firstpix=1+narray*(iblong-1)
         call ftgpve(unit,group,firstpix,narray,nullval,
     &               bl1d,anynull,status)

c --> Read in inclination angle.
         firstpix=1+narray*(ibtrans-1)
         call ftgpve(unit,group,firstpix,narray,nullval,
     &               bt1d,anynull,status)

c --> Read in azimuthal angle.
         firstpix=1+narray*(ibazim-1)
         call ftgpve(unit,group,firstpix,narray,nullval,
     &               ba1d,anynull,status)

c --> Convert to radians if necessary. 
         if(iaunit.eq.1) then
            dtor=pi/180.
            do 8 i=1,narray
               bt1d(i)=bt1d(i)*dtor
               ba1d(i)=ba1d(i)*dtor
    8       continue
         endif

c --> Transform to two components. 
         if(incflag.eq.0) then
            do 10 i=1,nx
               do 9 j=1,ny
                  i1d=i+(j-1)*nx
                  bl(i,j)=bl1d(i1d)*cos(bt1d(i1d))
                  bt(i,j)=bl1d(i1d)*sin(bt1d(i1d))
                  ba(i,j)=ba1d(i1d)
c --> Calculate Cartesian (image) components.
                  Bx(i,j)=bt(i,j)*cos(ba(i,j)+phase)
                  By(i,j)=bt(i,j)*sin(ba(i,j)+phase)
                  Bz(i,j)=bl(i,j)
    9          continue
   10       continue
         else
            do 12 i=1,nx
               do 11 j=1,ny
                  i1d=i+(j-1)*nx
                  bl(i,j)=bl1d(i1d)*sin(bt1d(i1d))
                  bt(i,j)=bl1d(i1d)*cos(bt1d(i1d))
                  ba(i,j)=ba1d(i1d)
c --> Calculate Cartesian (image) components.
                  Bx(i,j)=bt(i,j)*cos(ba(i,j)+phase)
                  By(i,j)=bt(i,j)*sin(ba(i,j)+phase)
                  Bz(i,j)=bl(i,j)
   11          continue
   12       continue
         endif

      endif

c --> The FITS file must always be closed before exiting the program. 
c --> Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      call ftclos(unit, status)
      call ftfiou(unit, status)

c --> Check for any error, and if so print out error messages.
      if(status.gt.0) call printerror(status)

      return
      end
