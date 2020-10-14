c***********************************************************************
      subroutine smooth(blpad,ixmin,ixmax,iymin,iymax,nap)
c***********************************************************************
c     Use a bi-cubic spline to make a smooth transition to zero field. *
c***********************************************************************
      include "include_file.sizes.f"
c-----------------------------------------------------------------------
      integer i,j,nmax
      integer ixmin,ixmax,iymin,iymax,nap
      real x1,y1,b,blpad(2*nxmax,2*nymax)
      parameter(nmax=2*max(nxmax,nymax))
      real,dimension(nmax) :: x,y
      real,dimension(nmax,nmax) :: bdum,adum
      real,dimension(nmax) :: xa,bi
c-----------------------------------------------------------------------
c --> Do interpolation for each boundary separately. 

      do 1 j=1,ny
         y(j)=float(j-1)/float(ny-1)
    1 continue

c --> x=0 boundary
      x(1)=-1./float(nx-1)
      x(2)=0.
      x(3)=float(nap)/float(nx-1)
   
      do j=1,ny
         bdum(1,j)=blpad(ixmin+1,j+iymin-1)
         bdum(2,j)=blpad(ixmin,j+iymin-1)
         bdum(3,j)=blpad(ixmin-1,j+iymin-1)
      enddo

      call cuspl2d(y,bdum,adum,3,ny)

      do i=1,nap-1
         xa(i)=float(i)/float(nx-1)
      enddo

      do j=iymin,iymax
         y1=float(j-iymin)/float(ny-1)
         call cuspleval2d(x,y,xa,y1,bdum,adum,bi,3,ny,nap-1)
         do i=1,nap-1
            blpad(ixmin-i,j)=bi(i)
         enddo
      enddo

c --> x=L boundary
      x(1)=1.-1./float(nx-1)
      x(2)=1.
      x(3)=1.+float(nap)/float(nx-1)

      do j=1,ny
         bdum(1,j)=blpad(ixmax-1,j+iymin-1)
         bdum(2,j)=blpad(ixmax,j+iymin-1)
         bdum(3,j)=blpad(ixmax+1,j+iymin-1)
      enddo

      call cuspl2d(y,bdum,adum,3,ny)

      do i=1,nap-1
         xa(i)=1.+float(i)/float(nx-1)
      enddo

      do j=iymin,iymax
         y1=float(j-iymin)/float(ny-1)
         call cuspleval2d(x,y,xa,y1,bdum,adum,bi,3,ny,nap-1)
         do i=1,nap-1
            blpad(ixmax+i,j)=bi(i)
         enddo
      enddo

c --> y=0 boundary
      do i=2,nx+1
         y(i)=float(i-2)/float(nx-1)
      enddo
      y(1)=-float(nap)/float(nx-1)
      y(nx+2)=1.+float(nap)/float(nx-1)

      x(1)=-1./float(ny-1)
      x(2)=0.
      x(3)=float(nap)/float(ny-1)

      do i=2,nx+1
         bdum(1,i)=blpad(i+ixmin-2,iymin+1)
         bdum(2,i)=blpad(i+ixmin-2,iymin)
         bdum(3,i)=blpad(i+ixmin-2,iymin-1)
      enddo
      bdum(1,1)=blpad(ixmin-nap,iymin+1)
      bdum(2,1)=blpad(ixmin-nap,iymin)
      bdum(3,1)=blpad(ixmin-nap,iymin-1)
      bdum(1,nx+2)=blpad(ixmax+nap,iymin+1)
      bdum(2,nx+2)=blpad(ixmax+nap,iymin)
      bdum(3,nx+2)=blpad(ixmax+nap,iymin-1)

      call cuspl2d(y,bdum,adum,3,nx+2)

      do i=1,nap-1
         xa(i)=float(i)/float(ny-1)
      enddo

      do i=ixmin-nap,ixmax+nap
         y1=float(i-ixmin)/float(nx-1)
         call cuspleval2d(x,y,xa,y1,bdum,adum,bi,3,nx+2,nap-1)
         do j=1,nap-1
            blpad(i,iymin-j)=bi(j)
         enddo
      enddo

c --> y=L boundary
      x(1)=1.-1./float(ny-1)
      x(2)=1.
      x(3)=1.+float(nap)/float(ny-1)

      do 13 i=2,nx+1
         bdum(1,i)=blpad(i+ixmin-2,iymax-1)
         bdum(2,i)=blpad(i+ixmin-2,iymax)
         bdum(3,i)=blpad(i+ixmin-2,iymax+1)
   13 continue
      bdum(1,1)=blpad(ixmin-nap,iymax-1)
      bdum(2,1)=blpad(ixmin-nap,iymax)
      bdum(3,1)=blpad(ixmin-nap,iymax+1)
      bdum(1,nx+2)=blpad(ixmax+nap,iymax-1)
      bdum(2,nx+2)=blpad(ixmax+nap,iymax)
      bdum(3,nx+2)=blpad(ixmax+nap,iymax+1)

      call cuspl2d(y,bdum,adum,3,nx+2)

      do i=1,nap-1
         xa(i)=1.+float(i)/float(ny-1)
      enddo

      do i=ixmin-nap,ixmax+nap
         y1=float(i-ixmin)/float(nx-1)
         call cuspleval2d(x,y,xa,y1,bdum,adum,bi,3,nx+2,nap-1)
         do j=1,nap-1
            blpad(i,iymax+j)=bi(j)
         enddo
      enddo

      return
      end
