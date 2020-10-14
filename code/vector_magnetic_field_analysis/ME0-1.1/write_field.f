c***********************************************************************
      subroutine write_field()
c***********************************************************************
c     Writes out the new azimuth angle in the same units as the input. *
c***********************************************************************
      include "include_file.sizes.f"
      include "include_file.field.f"
      include "include_file.flags.f"
c-----------------------------------------------------------------------
      integer i,j
      real pi,rtod,phase,ba(nxmax,nymax)
c-----------------------------------------------------------------------
      data pi/3.1415926535897932/
c-----------------------------------------------------------------------
      open(5,file='azimuth.dat')

c --> Recompute azimuth from ambiguity-resolved Cartesian components. 

      phase=0.5*float(iaflag)*pi
      do i=1,nx
         do j=1,ny
           ba(i,j)=atan2(By(i,j),Bx(i,j))-phase
         enddo
      enddo

c --> Convert back to degrees, if necessary.
      if(iaunit.eq.1) then
         rtod=180./pi
         do 2 i=1,nx
            do 1 j=1,ny
               ba(i,j)=ba(i,j)*rtod
    1       continue
    2    continue
      endif

c --> write out updated azimuth angle
      do 100 j=1,ny
         write(5,*) (ba(i,j),i=1,nx)
  100 continue

      return
      end
