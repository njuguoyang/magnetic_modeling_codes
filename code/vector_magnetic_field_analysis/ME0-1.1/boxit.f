
c***********************************************************************
      subroutine boxit(bthresh)
c***********************************************************************
c     Determine the corners of a single box which encloses all of the  *
c     pixels with transverse field strength above a threshold. If      *
c     there are no pixels above threshold, print a warning and do a    *
c     potential field acute angle ambiguity resolution.                *
c***********************************************************************
      include "include_file.sizes.f"
      include "include_file.bounds.f"
      include "include_file.field.f"
      include "include_file.verb.f"
c-----------------------------------------------------------------------
      integer i,j,nthresh
      real,dimension(nxmax,nymax) :: bl,bt,ba
      common/bobs/ bl,bt,ba
c-----------------------------------------------------------------------
c --> Initialize.
      nxa=nx+1
      nxb=-1
      nya=ny+1
      nyb=-1
      nthresh=0
    
c --> Loop over pixels, checking which are above threshold. 

      do 2 i=1,nx
         do 1 j=1,ny

            if(bt(i,j).ge.bthresh) then
               nthresh=nthresh+1

               if(nxb.lt.i) nxb=i
               if(nxa.gt.i) nxa=i
               if(nyb.lt.j) nyb=j
               if(nya.gt.j) nya=j

            endif

    1    continue
    2 continue

      if(nthresh.eq.0) then
         write(*,*) '*** WARNING: NO PIXELS ABOVE THRESHOLD ***'
         write(*,*) '    check value of bthresh in par file'
         call pacute(bthresh)
      endif
      if(iverb.eq.2) write(*,*) 'box:',nxa,nxb,nya,nyb

      return
      end
