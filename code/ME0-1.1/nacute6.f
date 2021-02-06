c***********************************************************************
      subroutine nacute6(bthresh)
c***********************************************************************
c     In pixels where the magnitude of the transverse field is below a *
c     given threshold, resolve the 180 degree ambiguity by picking the *
c     direction closest to the average direction of the field at       *
c     neighboring pixels. The solution is iterated as the average may  *
c     change. Based on the approach in UHIM code.                      *
c     This version starts with the pixels which have the most above-   *
c     theshold neightbors and works down from there, and allows all    *
c     pixels not yet visited to change.                                *
c***********************************************************************
      include "include_file.sizes.f"
      include "include_file.field.f"
      include "include_file.verb.f"
c-----------------------------------------------------------------------
      integer i,j,i1d,j1d,iter,niter,nflip,nflip1,nvisit
      integer i1,i2,i3,i4,i5,i6,i7,i8,itot,nmax,iflipt
      integer inn,jnn
      integer indx(nxmax*nymax)
      integer wt(nxmax,nymax),iflip(3,3)
      integer, dimension(9) :: iexp,ishift,jshift
      integer,dimension(nxmax,nymax) :: nnn,nnsav
      real Bxs,Bys,bdot,bdotb,bthresh,pi,bdotm
      real,dimension(nxmax,nymax) :: bl,bt,ba,bdotmax
      real,dimension(nxmax*nymax) :: bt1d,ba1d
      common/bobs/ bl,bt,ba
      equivalence(bt,bt1d)
      equivalence(ba,ba1d)
      data pi/3.1415926535897932/
c-----------------------------------------------------------------------
c --> Determine the number of above threshold neighbors for each pixel. 
c --> Also flag which pixels are themselves above threshold. 

      do 12 i=2,nx-1
         do 11 j=2,ny-1
            if(bt(i,j).gt.bthresh) then
               wt(i,j)=0
            else
               wt(i,j)=1
            endif
            nnn(i,j)=0
            if(bt(i+1,j-1).gt.bthresh) nnn(i,j)=nnn(i,j)+1
            if(bt(i,j-1).gt.bthresh) nnn(i,j)=nnn(i,j)+1
            if(bt(i-1,j-1).gt.bthresh) nnn(i,j)=nnn(i,j)+1
            if(bt(i-1,j).gt.bthresh) nnn(i,j)=nnn(i,j)+1
            if(bt(i+1,j).gt.bthresh) nnn(i,j)=nnn(i,j)+1
            if(bt(i-1,j+1).gt.bthresh) nnn(i,j)=nnn(i,j)+1
            if(bt(i,j+1).gt.bthresh) nnn(i,j)=nnn(i,j)+1
            if(bt(i+1,j+1).gt.bthresh) nnn(i,j)=nnn(i,j)+1
   11    continue
   12 continue

c --> Edge pixels need to be treated separately. 
      do 13 i=2,nx-1
         if(bt(i,1).gt.bthresh) then
            wt(i,1)=0
         else
            wt(i,1)=1
         endif
         nnn(i,1)=0
         if(bt(i-1,1).gt.bthresh) nnn(i,1)=nnn(i,1)+1
         if(bt(i+1,1).gt.bthresh) nnn(i,1)=nnn(i,1)+1
         if(bt(i-1,2).gt.bthresh) nnn(i,1)=nnn(i,1)+1
         if(bt(i,2).gt.bthresh) nnn(i,1)=nnn(i,1)+1
         if(bt(i+1,2).gt.bthresh) nnn(i,1)=nnn(i,1)+1

         if(bt(i,ny).gt.bthresh) then
            wt(i,ny)=0
         else
            wt(i,ny)=1
         endif
         nnn(i,ny)=0
         if(bt(i-1,ny-1).gt.bthresh) nnn(i,ny)=nnn(i,ny)+1
         if(bt(i,ny-1).gt.bthresh) nnn(i,ny)=nnn(i,ny)+1
         if(bt(i+1,ny-1).gt.bthresh) nnn(i,ny)=nnn(i,ny)+1
         if(bt(i-1,ny).gt.bthresh) nnn(i,ny)=nnn(i,ny)+1
         if(bt(i+1,ny).gt.bthresh) nnn(i,ny)=nnn(i,ny)+1
   13 continue

      do 14 j=2,ny-1
         if(bt(1,j).gt.bthresh) then
            wt(1,j)=0
         else
            wt(1,j)=1
         endif
         nnn(1,j)=0
         if(bt(1,j-1).gt.bthresh) nnn(1,j)=nnn(1,j)+1
         if(bt(1,j+1).gt.bthresh) nnn(1,j)=nnn(1,j)+1
         if(bt(2,j-1).gt.bthresh) nnn(1,j)=nnn(1,j)+1
         if(bt(2,j).gt.bthresh) nnn(1,j)=nnn(1,j)+1
         if(bt(2,j+1).gt.bthresh) nnn(1,j)=nnn(1,j)+1

         if(bt(nx,j).gt.bthresh) then
            wt(nx,j)=0
         else
            wt(nx,j)=1
         endif
         nnn(nx,j)=0
         if(bt(nx-1,j-1).gt.bthresh) nnn(nx,j)=nnn(nx,j)+1
         if(bt(nx-1,j).gt.bthresh) nnn(nx,j)=nnn(nx,j)+1
         if(bt(nx-1,j+1).gt.bthresh) nnn(nx,j)=nnn(nx,j)+1
         if(bt(nx,j-1).gt.bthresh) nnn(nx,j)=nnn(nx,j)+1
         if(bt(nx,j+1).gt.bthresh) nnn(nx,j)=nnn(nx,j)+1
   14 continue

c --> Finally corner pixels. 
      if(bt(1,1).gt.bthresh) then
         wt(1,1)=0
      else
         wt(1,1)=1
      endif
      nnn(1,1)=0
      if(bt(1,2).gt.bthresh) nnn(1,1)=nnn(1,1)+1
      if(bt(2,1).gt.bthresh) nnn(1,1)=nnn(1,1)+1
      if(bt(2,2).gt.bthresh) nnn(1,1)=nnn(1,1)+1

      if(bt(1,ny).gt.bthresh) then
         wt(1,ny)=0
      else
         wt(1,ny)=1
      endif
      nnn(1,ny)=0
      if(bt(1,ny-1).gt.bthresh) nnn(1,ny)=nnn(1,ny)+1
      if(bt(2,ny).gt.bthresh) nnn(1,ny)=nnn(1,ny)+1
      if(bt(2,ny-1).gt.bthresh) nnn(1,ny)=nnn(1,ny)+1

      if(bt(nx,1).gt.bthresh) then
         wt(nx,1)=0
      else
         wt(nx,1)=1
      endif
      nnn(nx,1)=0
      if(bt(nx-1,1).gt.bthresh) nnn(nx,1)=nnn(nx,1)+1
      if(bt(nx,2).gt.bthresh) nnn(nx,1)=nnn(nx,1)+1
      if(bt(nx-1,2).gt.bthresh) nnn(nx,1)=nnn(nx,1)+1

      if(bt(nx,ny).gt.bthresh) then
         wt(nx,ny)=0
      else
         wt(nx,ny)=1
      endif
      nnn(nx,ny)=0
      if(bt(nx-1,ny).gt.bthresh) nnn(nx,ny)=nnn(nx,ny)+1
      if(bt(nx,ny-1).gt.bthresh) nnn(nx,ny)=nnn(nx,ny)+1
      if(bt(nx-1,ny-1).gt.bthresh) nnn(nx,ny)=nnn(nx,ny)+1

      do i=1,nx
         do j=1,ny
            nnsav(i,j)=nnn(i,j)
         enddo
      enddo

c --> Arrays containing shifts about present pixel.
      ishift(1)=-1
      ishift(2)=0
      ishift(3)=1
      ishift(4)=-1
      ishift(5)=0
      ishift(6)=1
      ishift(7)=-1
      ishift(8)=0
      ishift(9)=1
      jshift(1)=-1
      jshift(2)=-1
      jshift(3)=-1
      jshift(4)=0
      jshift(5)=0
      jshift(6)=0
      jshift(7)=1
      jshift(8)=1
      jshift(9)=1

c --> Create an index array for pixels sorted on transverse field.

      call sortrx(nxmax*nymax,-bt1d,indx)

      jnn=8
      nflip=0
      dowhile(jnn.gt.0)

c --> Start with pixels with most neighbors above threshold. 
         inn=8

         nvisit=1

c --> Keep iterating as long as any pixels remain with a given value of
c --> nnn.
         dowhile(inn.ge.jnn)
            nvisit=0

c --> This is a local test, so loop over pixels and apply at each one. 

            do 2 i1d=1,nx*ny
               j1d=indx(i1d)

c --> 2-d indices for this pixel.
               j=1+j1d/nxmax
               i=j1d-(j-1)*nxmax

c --> Only apply to points below threshold and not already visited.

               if(wt(i,j).gt.0.5.and.nnsav(i,j).ge.inn) then

c --> Compute dot product of transverse field at this pixel with
c --> transverse field at neighboring pixels. 
c --> Also update neighbors to indicate this pixel has been visited.

c --> If edge pixel, has fewer neighbors.
                  if(i.eq.1) then
                     if(j.eq.1) then
                        wt(1,1)=0
                        nnn(1,2)=nnn(1,2)+1
                        nnn(2,1)=nnn(2,1)+1
                        nnn(2,2)=nnn(2,2)+1
                        Bxs=Bx(i+1,j+1)+Bx(i,j+1)+Bx(i+1,j)
                        Bys=By(i+1,j+1)+By(i,j+1)+By(i+1,j)
                     elseif(j.eq.ny) then
                        wt(1,ny)=0
                        nnn(1,ny-1)=nnn(1,ny-1)+1
                        nnn(2,ny)=nnn(2,ny)+1
                        nnn(2,ny-1)=nnn(2,ny-1)+1
                        Bxs=Bx(i+1,j-1)+Bx(i,j-1)+Bx(i+1,j)
                        Bys=By(i+1,j-1)+By(i,j-1)+By(i+1,j)
                     else
                        wt(1,j)=0
                        nnn(1,j-1)=nnn(1,j-1)+1
                        nnn(2,j-1)=nnn(2,j-1)+1
                        nnn(2,j)=nnn(2,j)+1
                        nnn(2,j+1)=nnn(2,j+1)+1
                        nnn(1,j+1)=nnn(1,j+1)+1
                        Bxs=Bx(i,j-1)+Bx(i+1,j-1)+Bx(i+1,j)
     +                     +Bx(i,j+1)+Bx(i+1,j+1)
                        Bys=By(i,j-1)+By(i+1,j-1)+By(i+1,j)
     +                     +By(i,j+1)+By(i+1,j+1)
                     endif
               bdotb=Bxs*Bx(i,j)+Bys*By(i,j)
               iflipt=0
               if(bdotb.lt.0.) iflipt=1
                  elseif(i.eq.nx) then
                     if(j.eq.1) then
                        wt(nx,1)=0
                        nnn(nx,2)=nnn(nx,2)+1
                        nnn(nx-1,1)=nnn(nx-1,1)+1
                        nnn(nx-1,2)=nnn(nx-2,1)+1
                        Bxs=Bx(i-1,j+1)+Bx(i,j+1)+Bx(i-1,j)
                        Bys=By(i-1,j+1)+By(i,j+1)+By(i-1,j)
                     elseif(j.eq.ny) then
                        wt(nx,ny)=0
                        nnn(nx,ny-1)=nnn(nx,ny-1)+1
                        nnn(nx-1,ny)=nnn(nx-1,ny)+1
                        nnn(nx-1,ny-1)=nnn(nx-1,ny-1)+1
                        Bxs=Bx(i-1,j-1)+Bx(i,j-1)+Bx(i-1,j)
                        Bys=By(i-1,j-1)+By(i,j-1)+By(i-1,j)
                     else
                        wt(nx,j)=0
                        nnn(nx,j-1)=nnn(nx,j-1)+1
                        nnn(nx-1,j-1)=nnn(nx-1,j-1)+1
                        nnn(nx-1,j)=nnn(nx-1,j)+1
                        nnn(nx-1,j+1)=nnn(nx-1,j+1)+1
                        nnn(nx,j+1)=nnn(nx,j+1)+1
                        Bxs=Bx(i,j-1)+Bx(i-1,j-1)+Bx(i-1,j)
     +                     +Bx(i,j+1)+Bx(i-1,j+1)
                        Bys=By(i,j-1)+By(i-1,j-1)+By(i-1,j)
     +                     +By(i,j+1)+By(i-1,j+1)
                     endif
               bdotb=Bxs*Bx(i,j)+Bys*By(i,j)
               iflipt=0
               if(bdotb.lt.0.) iflipt=1
                  else
                     if(j.eq.1) then
                        wt(i,1)=0
                        nnn(i-1,1)=nnn(i-1,1)+1
                        nnn(i-1,2)=nnn(i,2-1)+1
                        nnn(i,2)=nnn(i,2)+1
                        nnn(i+1,2)=nnn(i,2+1)+1
                        nnn(i+1,1)=nnn(i+1,1)+1
                        Bxs=Bx(i-1,j+1)+Bx(i-1,j)+Bx(i,j+1)
     +                     +Bx(i+1,j+1)+Bx(i+1,j)
                        Bys=By(i-1,j+1)+By(i-1,j)+By(i,j+1)
     +                     +By(i+1,j+1)+By(i+1,j)
               bdotb=Bxs*Bx(i,j)+Bys*By(i,j)
               iflipt=0
               if(bdotb.lt.0.) iflipt=1
                     elseif(j.eq.ny) then
                        wt(i,ny)=0
                        nnn(i-1,ny)=nnn(i-1,ny)+1
                        nnn(i-1,ny-1)=nnn(i-1,ny-1)+1
                        nnn(i,ny-1)=nnn(i,ny-1)+1
                        nnn(i+1,ny-1)=nnn(i+1,ny-1)+1
                        nnn(i+1,ny)=nnn(i+1,ny)+1
                        Bxs=Bx(i-1,j-1)+Bx(i-1,j)+Bx(i,j-1)
     +                     +Bx(i+1,j-1)+Bx(i+1,j)
                        Bys=By(i-1,j-1)+By(i-1,j)+By(i,j-1)
     +                     +By(i+1,j-1)+By(i+1,j)
               bdotb=Bxs*Bx(i,j)+Bys*By(i,j)
               iflipt=0
               if(bdotb.lt.0.) iflipt=1
                     else
                        nnn(i-1,j-1)=nnn(i-1,j-1)+1
                        nnn(i-1,j)=nnn(i-1,j)+1
                        nnn(i-1,j+1)=nnn(i-1,j+1)+1
                        nnn(i,j-1)=nnn(i,j-1)+1
                        nnn(i,j+1)=nnn(i,j+1)+1
                        nnn(i+1,j-1)=nnn(i+1,j-1)+1
                        nnn(i+1,j)=nnn(i+1,j)+1
                        nnn(i+1,j+1)=nnn(i+1,j+1)+1
c --> Loop over possible directions of transverse field at each pixel. 
c --> If pixel has already been visited or is above threshold, do not
c --> allow it to flip.
            bdotm=-1.e9
            do 21 i1=0,wt(i-1,j-1)
             iexp(1)=i1
             do 22 i2=0,wt(i,j-1)
              iexp(2)=i2
              do 23 i3=0,wt(i+1,j-1)
               iexp(3)=i3
               do 24 i4=0,wt(i-1,j)
                iexp(4)=i4
                do 25 i5=0,wt(i,j)
                 iexp(5)=i5
                 do 26 i6=0,wt(i+1,j)
                  iexp(6)=i6
                  do 27 i7=0,wt(i-1,j+1)
                   iexp(7)=i7
                   do 28 i8=0,wt(i,j+1)
                    iexp(8)=i8
                    do 29 i9=0,wt(i+1,j+1)
                     iexp(9)=i9
                     bdot=0.
                     do 20 j1=1,9
                      Bxs=Bx(i+ishift(j1),j+jshift(j1))*(-1.)**iexp(j1)
                      Bys=By(i+ishift(j1),j+jshift(j1))*(-1.)**iexp(j1)
                      if(j1.ne.5) then
                      bdot=bdot+(Bxs*Bx(i,j)+Bys*By(i,j))*(-1.)**iexp(5)
                      endif
   20                continue


                     if(bdot.gt.bdotm) then
                       bdotm=bdot
                       iflipt=i5
                     endif
   29               continue
   28              continue
   27             continue
   26            continue
   25           continue
   24          continue
   23         continue
   22        continue
   21       continue
                        wt(i,j)=0

                     endif
                  endif

                  nvisit=nvisit+1

c --> Flip direction of transverse field if it improves dot product.
                  if(iflipt.gt.0.5) then
                     nflip=nflip+1
                     Bx(i,j)=-Bx(i,j)
                     By(i,j)=-By(i,j)
                  endif

               endif
    2       continue
         if(iverb.eq.2) write(*,*) 'iteration',nflip,inn,nvisit
            if(nvisit.eq.0) then
               inn=inn-1
            else
               inn=8
               do i=1,nx
                  do j=1,ny
                     nnsav(i,j)=nnn(i,j)
                  enddo
               enddo
            endif
         enddo

         jnn=jnn-1
      enddo
      if(iverb.eq.2) then
         write(*,*) 'total number of flips',nflip
         write(*,*) 'fraction flipped',float(nflip)/float(nx*ny)
      endif

      return
      end
