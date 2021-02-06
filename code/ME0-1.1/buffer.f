c***********************************************************************
      subroutine buffer(bl,blpad,ixmin,ixmax,iymin,iymax,npad,nxp,nyp)
c***********************************************************************
c     Pad the input array with zeros.                                  *
c***********************************************************************
      include "include_file.sizes.f"
c-----------------------------------------------------------------------
      integer i,j
      integer ixmin,ixmax,iymin,iymax,npad,nxp,nyp
      real bl(nxmax,nymax),blpad(2*nxmax,2*nymax)
c-----------------------------------------------------------------------
c --> non-zero data will be centered

      ixmin=npad+1
      ixmax=ixmin+nx-1
      iymin=npad+1
      iymax=iymin+ny-1

      nxp=2*npad+nx
      nyp=2*npad+ny

      do 2 i=1,ixmin-1
         do 1 j=iymin,iymax
            blpad(i,j)=0.
    1    continue
    2 continue

      do 4 i=ixmax+1,nxp
         do 3 j=iymin,iymax
            blpad(i,j)=0.
    3    continue
    4 continue

      do 6 i=1,nxp
         do 5 j=1,iymin-1
            blpad(i,j)=0.
    5    continue
    6 continue

      do 8 i=1,nxp
         do 7 j=iymax+1,nyp
            blpad(i,j)=0.
    7    continue
    8 continue

      do 10 i=1,nx
         do 9 j=1,ny
            blpad(i+ixmin-1,j+iymin-1)=bl(i,j)
    9    continue
   10 continue

      return
      end
