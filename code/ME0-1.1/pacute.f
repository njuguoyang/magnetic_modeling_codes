
c***********************************************************************
      subroutine pacute(bthresh)
c***********************************************************************
c     In pixels where the magnitude of the transverse field is below a *
c     given threshold, resolve the 180 degree ambiguity by picking the *
c     direction closest to the potential field.                        *
c***********************************************************************
      include "include_file.sizes.f"
      include "include_file.potential.f"
      include "include_file.point.f"
      include "include_file.verb.f"
c-----------------------------------------------------------------------
      real,dimension(nxmax,nymax) :: bl,bt,ba
      common/bobs/ bl,bt,ba
      data pi/3.1415926535897932/
c-----------------------------------------------------------------------
c --> This is a local test, so loop over pixels and apply at each one. 

      nflip=0
      do 2 i=1,nx
         do 1 j=1,ny

c --> Only apply to points below threshold.
            if(bt(i,j).le.bthresh) then

c --> Potential field is given in helioplanar components, so need to 
c --> transform to image plane components
               bipx=c11*Bpx(i,j)+c12*Bpy(i,j)+a31*Bpz(i,j)
               bipy=c21*Bpx(i,j)+c22*Bpy(i,j)+a32*Bpz(i,j)
               bipz=a13*Bpx(i,j)+a23*Bpy(i,j)+a33*Bpz(i,j)

c --> Compute dot product of observed transverse field with potential field.
               bdotb=bipy*cos(ba(i,j))-bipx*sin(ba(i,j))

c --> Flip direction of transverse field if dot product is negative. 
               if(bdotb.lt.0.) then
                  nflip=nflip+1
                  if(ba(i,j).lt.0.) then
                     ba(i,j)=ba(i,j)+pi
                  else
                     ba(i,j)=ba(i,j)-pi
                  endif
               endif

            endif

    1    continue
    2 continue
      if(iverb.ge.1) then
         write(*,*) 'number of flips',nflip
         write(*,*) 'fraction flipped',float(nflip)/float(nx*ny)
      endif

      return
      end
