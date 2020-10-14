c***********************************************************************************************************************************
      subroutine get_ran_pix
      implicit none
c===================================================================================================================================
c This subroutine determines the coordinates for a sequence of pixels with a random start point
c
c 2009, Ashley Crouch, ash@cora.nwra.com
c===================================================================================================================================
      include "include_file.sizes.f"
      include "include_file.bounds.f"
      include "include_file.seed.f"
      include "include_file.ran_pix.f"

      integer max_tries
      parameter (max_tries=50)

      integer i,j,ii,jj,ntries,get_pix
      real ran3

      if ((nxjump.eq.1).and.(nyjump.eq.1)) then
c
c Select a single random pixel
c
         i=int(nxny*ran3(seed))
         if (i.ge.nxny) i=nxny-1
         j=i/dnx
         i=i-j*dnx+nxa
         j=j+nya

         nxg=1
         ivec(nxg)=i
         nyg=1
         jvec(nyg)=j
      else
c
c Select a random start point and avoid the grid from the previous sequence
c
         ntries=0
         get_pix=1
         do while ((get_pix.eq.1).and.(ntries.le.max_tries))
            ia=int(ran3(seed)*jumpsq)
            if (ia.eq.jumpsq) ia=jumpsq-1
            ja=ia/jump
            ia=ia-ja*jump+nxa
            ja=ja+nya
            if ((ia.ne.ia_prev).and.(ja.ne.ja_prev)) get_pix=0
            ntries=ntries+1
         enddo
         if (ntries.gt.max_tries) then
            print*,'too many attempts'
            stop
         endif
         ia_prev=ia
         ja_prev=ja
c
c Construct the sequence of pixels with random start point and separation: jump
c
         nxg=0
         ii=1
         get_pix=1
         do while (get_pix.eq.1)
            i=ia+(ii-1)*jump
            if (i.le.nxb) then
               nxg=nxg+1
               ivec(nxg)=i
            else
               get_pix=0
            endif
            ii=ii+1
         enddo
         nyg=0
         jj=1
         get_pix=1
         do while (get_pix.eq.1)
            j=ja+(jj-1)*jump
            if (j.le.nyb) then
               nyg=nyg+1
               jvec(nyg)=j
            else
               get_pix=0
            endif
            jj=jj+1
         enddo
      endif
      return
      end
c***********************************************************************************************************************************
