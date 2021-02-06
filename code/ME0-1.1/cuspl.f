c***********************************************************************
      subroutine cuspl(x,f,a,n)
c***********************************************************************
c     Computes one dimensional cubic spline coefficients.              *
c***********************************************************************
      include "include_file.sizes.f"
c-----------------------------------------------------------------------
      integer i,n
      real,dimension(2*max(nxmax,nymax)) :: x,f,a,g
c-----------------------------------------------------------------------
c --> Natural boundary condition. 
      a(1)=0.
      g(1)=0.
c --> Require first derivative to be 0. at the first boundary. 
c      a(1)=-0.5
c      g(1)=3.*(f(2)-f(1))/(x(2)-x(1))**2
      do i=2,n-1
         dxp=x(i+1)-x(i)
         dxm=x(i)-x(i-1)
         a(i)=-dxp/(a(i-1)*dxm+2.*(dxm+dxp))
         g(i)=(6.*((f(i+1)-f(i))/dxp-(f(i)-f(i-1))/dxm)
     +        -dxm*g(i-1))/(a(i-1)*dxm+2.*(dxm+dxp))
      enddo

c --> Require first derivative to be 0. at the first boundary. 
      g(n)=-3.*(f(n)-f(n-1))/(x(n)-x(1))**2
      a(n)=(2.*g(n)-g(n-1))/(a(n-1)+2.)
c --> Natural boundary condition.
c      a(n)=0.
      do i=n-1,1,-1
         a(i)=a(i)*a(i+1)+g(i)
      enddo

      return
      end
