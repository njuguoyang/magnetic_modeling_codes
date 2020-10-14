c***********************************************************************
      subroutine nacute5(bthresh)
c***********************************************************************
c     In pixels where the magnitude of the transverse field is below a *
c     given threshold, resolve the 180 degree ambiguity by picking the *
c     direction closest to the average direction of the field at       *
c     neighboring pixels. The solution is iterated as the average may  *
c     change. Based on the approach in UHIM code.                      *
c     This version starts with the pixel with the strongest field and  *
c     works radially out from there.                                   *
c***********************************************************************
      include "include_file.sizes.f"
      include "include_file.field.f"
      include "include_file.verb.f"
c-----------------------------------------------------------------------
      integer i,j,i1d,j1d,iter,niter,nflip,nflip1,nvisit
      integer imax,jmax
      integer indx(nxmax*nymax)
      integer,dimension(nxmax,nymax) :: nnn
      real Bxs,Bys,bdotb,bthresh,pi,bmax,bmag,dmax
      real,dimension(nxmax,nymax) :: bl,bt,ba,wt,wtemp
      real,dimension(nxmax*nymax) :: bt1d,ba1d,dist1d
      common/bobs/ bl,bt,ba
      equivalence(bt,bt1d)
      equivalence(ba,ba1d)
      data pi/3.1415926535897932/
c-----------------------------------------------------------------------
c --> Determine the pixel with the strongest field. 

      bmax=-1.
      do 12 i=1,nx
         do 11 j=1,ny
            bmag=Bx(i,j)**2+By(i,j)**2+Bz(i,j)**2
            if(bmag.gt.bmax) then
              bmax=bmag
              imax=i
              jmax=j
            endif
   11    continue
   12 continue
      if(iverb.eq.2) write(*,*) 'max field at',imax,jmax,bmax

c --> Calculate distance from each pixel to strongest field pixel.
      dmax=float(nxmax)*float(nymax)
      do 14 i=1,nxmax
         do 13 j=1,nymax
            i1d=(j-1)*nxmax+i
            if(i.gt.nx.or.j.gt.ny) then
               dist1d(i1d)=dmax
            else
               dist1d(i1d)=sqrt(float(i-imax)**2+float(j-jmax)**2)
            endif
   13    continue
   14 continue

c --> Create an index array for pixels sorted on distance.

      call sortrx(nxmax*nymax,dist1d,indx)

c --> Maximum number of iterations. 
      niter=50

      nflip=0
      nflip1=1
      iter=1
      dowhile(nflip1.gt.0.and.iter.lt.niter)

c --> This is a local test, so loop over pixels and apply at each one. 

         nflip1=0
         do 2 i1d=1,nx*ny
            j1d=indx(i1d)

c --> 2-d indices for this pixel.
            j=1+j1d/nxmax
            i=j1d-(j-1)*nxmax

c --> Only apply to points below threshold.

            if(bt1d(j1d).le.bthresh) then

c --> Compute dot product of transverse field at this pixel with
c --> transverse field at neighboring pixels. 

c --> If edge pixel, has fewer neighbors.
               if(i.eq.1) then
                  if(j.eq.1) then
                     Bxs=Bx(i+1,j+1)+Bx(i,j+1)+Bx(i+1,j)
                     Bys=By(i+1,j+1)+By(i,j+1)+By(i+1,j)
                  elseif(j.eq.ny) then
                     Bxs=Bx(i+1,j-1)+Bx(i,j-1)+Bx(i+1,j)
                     Bys=By(i+1,j-1)+By(i,j-1)+By(i+1,j)
                  else
                     Bxs=Bx(i,j-1)+Bx(i+1,j-1)+Bx(i+1,j)
     +                  +Bx(i,j+1)+Bx(i+1,j+1)
                     Bys=By(i,j-1)+By(i+1,j-1)+By(i+1,j)
     +                  +By(i,j+1)+By(i+1,j+1)
                  endif
               elseif(i.eq.nx) then
                  if(j.eq.1) then
                     Bxs=Bx(i-1,j+1)+Bx(i,j+1)+Bx(i-1,j)
                     Bys=By(i-1,j+1)+By(i,j+1)+By(i-1,j)
                  elseif(j.eq.ny) then
                     Bxs=Bx(i-1,j-1)+Bx(i,j-1)+Bx(i-1,j)
                     Bys=By(i-1,j-1)+By(i,j-1)+By(i-1,j)
                  else
                     Bxs=Bx(i,j-1)+Bx(i-1,j-1)+Bx(i-1,j)
     +                  +Bx(i,j+1)+Bx(i-1,j+1)
                     Bys=By(i,j-1)+By(i-1,j-1)+By(i-1,j)
     +                  +By(i,j+1)+By(i-1,j+1)
                  endif
               else
                  if(j.eq.1) then
                     Bxs=Bx(i-1,j+1)+Bx(i-1,j)+Bx(i,j+1)
     +                  +Bx(i+1,j+1)+Bx(i+1,j)
                     Bys=By(i-1,j+1)+By(i-1,j)+By(i,j+1)
     +                  +By(i+1,j+1)+By(i+1,j)
                  elseif(j.eq.ny) then
                     Bxs=Bx(i-1,j-1)+Bx(i-1,j)+Bx(i,j-1)
     +                  +Bx(i+1,j-1)+Bx(i+1,j)
                     Bys=By(i-1,j-1)+By(i-1,j)+By(i,j-1)
     +                  +By(i+1,j-1)+By(i+1,j)
                  else
                     Bxs=Bx(i-1,j-1)+Bx(i,j-1)+Bx(i+1,j-1)
     +                  +Bx(i-1,j+1)+Bx(i,j+1)+Bx(i+1,j+1)
     +                  +Bx(i-1,j)+Bx(i+1,j)
                     Bys=By(i-1,j-1)+By(i,j-1)+By(i+1,j-1)
     +                  +By(i-1,j+1)+By(i,j+1)+By(i+1,j+1)
     +                  +By(i-1,j)+By(i+1,j)
                  endif
               endif

               bdotb=Bxs*Bx(i,j)+Bys*By(i,j)

c --> Flip direction of transverse field if dot product is negative.
               if(bdotb.lt.0.) then
                  nflip1=nflip1+1
                  Bx(i,j)=-Bx(i,j)
                  By(i,j)=-By(i,j)
               endif

            endif

    2    continue
         iter=iter+1
         nflip=nflip+nflip1
         if(iverb.eq.2) write(*,*) 'iteration',iter,nflip1,nflip
      enddo
      if(iverb.eq.2) write(*,*) 'total number of flips',nflip
      if(iverb.eq.2) write(*,*) 'fraction flipped',float(nflip)/float(nx*ny)

      return
      end
