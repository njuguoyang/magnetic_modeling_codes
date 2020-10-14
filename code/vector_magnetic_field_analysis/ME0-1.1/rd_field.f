c***********************************************************************
      subroutine rd_field(filename,npad)
c***********************************************************************
c     Reads in magnetic field data plus pointing information.          *
c***********************************************************************
      include "include_file.sizes.f"
      include "include_file.point.f"
      include "include_file.flags.f"
      include "include_file.field.f"
      include "include_file.pix_size.f"
c-----------------------------------------------------------------------
      integer i,j
      real phase,pi,dtor,radius
      real,dimension(nxmax,nymax) :: bl,bt,ba
      common/bobs/ bl,bt,ba
      character*90 filename
      data pi/3.1415926535897932/
c-----------------------------------------------------------------------
      open(2,file=filename,status='old')

c --> Read in array dimensions.

      read(2,*) nx,ny

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

c --> Read in pixel size.

      read(2,*) dxi,dyi

c --> Read in ephemeris information:
c --> b=solar b-angle
c --> p=solar p-angle
c --> radius=solar radius

      read(2,*) b,p,radius

c --> Read in pointing information:
c --> theta=central meridian angle
c --> phi=latitude

      read(2,*) theta,phi

c --> Check whether conversion to degrees is needed for pointing and
c --> ephemeris.
      if(ipunit.eq.1) then
         dtor=pi/180.
         b=b*dtor
         p=p*dtor
         theta=theta*dtor
         phi=phi*dtor
      endif

c --> Read in line of sight field.
      do 1 j=1,ny
         read(2,*) (bl(i,j),i=1,nx)
    1 continue

c --> Read in transverse field.
      do 2 j=1,ny
         read(2,*) (bt(i,j),i=1,nx)
    2 continue

c --> Read in azimuthal angle.
c --> Check whether conversion to degrees is needed for azimuth.
      if(iaunit.eq.1) then
         dtor=pi/180.
         do 4 j=1,ny
            read(2,*) (ba(i,j),i=1,nx)
            do 3 i=1,nx
               ba(i,j)=ba(i,j)*dtor
    3       continue
    4    continue
      else
         do 5 j=1,ny
            read(2,*) (ba(i,j),i=1,nx)
    5    continue
      endif

c --> Calculate Cartesian (image) components.
      phase=0.5*float(iaflag)*pi
      do i=1,nx
         do j=1,ny
            Bx(i,j)=bt(i,j)*cos(ba(i,j)+phase)
            By(i,j)=bt(i,j)*sin(ba(i,j)+phase)
            Bz(i,j)=bl(i,j)
         enddo
      enddo

      return
      end
