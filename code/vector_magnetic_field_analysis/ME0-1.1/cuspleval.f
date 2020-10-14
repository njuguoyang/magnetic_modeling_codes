c***********************************************************************
      subroutine cuspleval(x,x1,f,a,fi,n)
c***********************************************************************
c     Evaluates cubic spline interpolation with coefficients computed  *
c     in cuspl.f                                                       *
c***********************************************************************
      include "include_file.sizes.f"
c-----------------------------------------------------------------------
      integer i,im,ip,n
      real g1,g2,dx,x1,fi
      real,dimension(2*max(nxmax,nymax)) :: x,f,a
c-----------------------------------------------------------------------
c --> Use bisection to find the indices of the points bracketing the 
c --> interpolation point. 
      im=1
      ip=n
      do while(ip-im.gt.1)
         i=(ip+im)/2
         if(x(i).gt.x1)then
            ip=i
         else
            im=i
         endif
      enddo
      dx=x(ip)-x(im)

c --> Check bracketing points are distinct.
      if(dx.le.0.) then
         write(*,*) 'cuspleval',x(ip),x1,x(im)
         stop
      endif

c --> Evaluate the interpolation. 
      g1=(x(ip)-x1)/dx
      g2=(x1-x(im))/dx
      fi=g1*f(im)+g2*f(ip)
     +  +((g1**3-g1)*a(im)+(g2**3-g2)*a(ip))*(dx**2)/6.

      return
      end
