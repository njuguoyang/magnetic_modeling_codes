c***********************************************************************
      subroutine cuspleval2d(x,y,x1,y1,f,a,fi,n1,n2,n3)
c***********************************************************************
c     Evaluates 2-D cubic spline interpolation with coefficients       *
c     computed in cuspl2d.f                                            *
c***********************************************************************
      include "include_file.sizes.f"
c-----------------------------------------------------------------------
      integer i,j,n1,n2,nmax
      real y1,fdum
      parameter(nmax=2*max(nxmax,nymax))
      real,dimension(nmax,nmax) :: f,a
      real,dimension(nmax) :: x1,fi,x,y,f1d,ffdum,a1d
c-----------------------------------------------------------------------
c --> Evaluate the row splines.
      do i=1,n1
         do j=1,n2
            f1d(j)=f(i,j)
            a1d(j)=a(i,j)
         enddo
         call cuspleval(y,y1,f1d,a1d,ffdum(i),n2)
      enddo

c --> Construct and evaluate the column spline.
      call cuspl(x,ffdum,a1d,n1)

c --> Loop over values of x1, evaluating the interpolation.
      do i=1,n3
         call cuspleval(x,x1(i),ffdum,a1d,fdum,n1)
         fi(i)=fdum
      enddo

      return
      end
