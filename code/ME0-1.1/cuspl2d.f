c***********************************************************************
      subroutine cuspl2d(y,f,a,n1,n2)
c***********************************************************************
c     Computes 2-D spline coefficients.                                *
c***********************************************************************
      include "include_file.sizes.f"
c-----------------------------------------------------------------------
      integer i,j,n1,n2,nmax
      parameter(nmax=2*max(nxmax,nymax))
      real,dimension(nmax,nmax) :: f,a
      real,dimension(nmax) :: y,fdum,adum
c-----------------------------------------------------------------------
      do i=1,n1
         do j=1,n2
            fdum(j)=f(i,j)
         enddo

         call cuspl(y,fdum,adum,n2)

         do j=1,n2
            a(i,j)=adum(j)
         enddo
      enddo

      return
      end
