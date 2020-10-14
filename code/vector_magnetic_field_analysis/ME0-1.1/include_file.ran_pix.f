c***********************************************************************************************************************************
c
c Parameter jump controls the size of the jump between pixels in the sequence of attempted reconfigurations
c
      integer jump,jumpsq
      parameter (jump=5)
      parameter (jumpsq=jump*jump)
      integer nxmax_jump,nymax_jump
      parameter (nxmax_jump=nxmax/jump)
      parameter (nymax_jump=nymax/jump)
c
c Common stuff for the sequence of random pixels
c
      integer nxjump,nyjump,nxg,nyg,ia,ia_prev,ja,ja_prev
      common /ran_pix_com/ nxjump,nyjump,nxg,nyg,ia,ia_prev,ja,ja_prev
      integer ivec(nxmax_jump+1),jvec(nymax_jump+1)
      common /ran_vec_com/ ivec,jvec
c
c The number of pixels
c
      integer nxny
      common /ijkcom/ nxny
c***********************************************************************************************************************************
