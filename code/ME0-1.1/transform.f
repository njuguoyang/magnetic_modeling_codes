
c***********************************************************************
      subroutine transform()
c***********************************************************************
c     Construct coordinate transformation matrices from Gary and       *
c     Hagyard (1990).                                                  *
c     Note that the matrix C is the transpose of A, which is also the  *
c     inverse of A, since A is orthogonal.                             *
c***********************************************************************
      include "include_file.point.f"
c-----------------------------------------------------------------------

      c11=-sin(b)*sin(p)*sin(theta) + cos(p)*cos(theta)
      c12=-sin(phi)*(sin(b)*sin(p)*cos(theta) + cos(p)*sin(theta))
     +    -cos(phi)*cos(b)*sin(p)
      c21= sin(b)*cos(p)*sin(theta) + sin(p)*cos(theta)
      c22= sin(phi)*(sin(b)*cos(p)*cos(theta) - sin(p)*sin(theta))
     +    +cos(phi)*cos(b)*cos(p)

      a13=-cos(b)*sin(theta)
      a23=-cos(b)*sin(phi)*cos(theta) + sin(b)*cos(phi)
      a31= cos(phi)*(sin(b)*sin(p)*cos(theta) + cos(p)*sin(theta))
     +    -sin(phi)*cos(b)*sin(p)
      a32=-cos(phi)*(sin(b)*cos(p)*cos(theta) - sin(p)*sin(theta))
     +    +sin(phi)*cos(b)*cos(p)
      a33= cos(phi)*cos(b)*cos(theta) + sin(phi)*sin(b)

      return
      end
