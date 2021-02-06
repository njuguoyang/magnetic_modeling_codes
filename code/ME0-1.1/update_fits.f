c***********************************************************************
      subroutine update_fits(filename)
c***********************************************************************
c     Overwrites the azimuth angle in a FITS file. Adds a keyword      *
c     comment that the ambiguity has been resolved using ME0.          *
c***********************************************************************
      include "include_file.sizes.f"
      include "include_file.field.f"
      include "include_file.fits.f"
      include "include_file.flags.f"
      include "include_file.point.f"
c-----------------------------------------------------------------------
      integer status,unit,readwrite,blocksize
      integer fpixels(3),lpixels(3),group,i,j
      real rtod,phase
      real pi,ba1d(nxmax*nymax)
      real,dimension(nxmax,nymax) :: bl,bt,ba
      character*90 filename
      character*80 comment
      common/bobs/ bl,bt,ba
      data pi/3.1415926535897932/
c-----------------------------------------------------------------------
c --> The status parameter must always be initialized.
      status=0

c --> Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)

c --> Open the FITS file, with read-write access.  The returned blocksize
c --> parameter is obsolete and should be ignored. 
      readwrite=1
      call ftopen(unit,filename,readwrite,blocksize,status)

c --> Recompute azimuth from ambiguity-resolved Cartesian components
c --> unless transverse field is zero, in which case use original angle.

      phase=0.5*float(iaflag)*pi
      do i=1,nx
         do j=1,ny
            if(bt(i,j).ne.0.) then
               ba1d(i+(j-1)*nxmax)=atan2(By(i,j),Bx(i,j))-phase
            else
               ba1d(i+(j-1)*nxmax)=ba(i,j)
            endif
         enddo
      enddo

c --> Convert back to degrees, if necessary.
      if(iaunit.eq.1) then
         rtod=180./pi
         do i=1,nxmax*nymax
            ba1d(i)=ba1d(i)*rtod
         enddo
      endif

c --> Initialize variables for writing azimuth. 
      group=1

c --> First and last pixels in each array dimension.
      fpixels(1)=1
      fpixels(2)=1
      fpixels(3)=ibazim
      lpixels(1)=naxes(1)
      lpixels(2)=naxes(2)
      lpixels(3)=ibazim

c --> Write new azimuthal angle.

      call ftpsse(unit,group,3,naxes,fpixels,lpixels,ba1d,status)

c --> Append a comment that the ambiguity has been resolved using ME0.

      comment='Ambiguity has been resolved using ME0 code.'

      call ftpcom(unit,comment,status)

c --> The FITS file must always be closed before exiting the program. 
c --> Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      call ftclos(unit, status)
      call ftfiou(unit, status)

c --> Check for any error, and if so print out error messages.
c --> The PRINTERROR subroutine is listed near the end of this file.
      if(status.gt.0) call printerror(status)

      return
      end
