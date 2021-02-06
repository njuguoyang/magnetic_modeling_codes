c***********************************************************************************************************************************
      subroutine CalcDE_reconfig(i,j)
c===================================================================================================================================
c For a reconfiguration of the azimuthal angle at pixel (i,j) this subroutine computes Div(B) and Jz at each point with a finite 
c differencing stencil that references pixel (i,j). For most points this involves four neighbouring pixels. At the corners and edges 
c it involves less than four pixels.
c
c 2009, Ashley Crouch, ash@cora.nwra.com
c===================================================================================================================================
      implicit none

      include "include_file.sizes.f"
      include "include_file.field.f"
      include "include_file.energy.f"
      include "include_file.reconfig.f"
      include "include_file.derivs.f"
      include "include_file.bounds.f"
      include "include_file.point.f"

      integer i,j
      real Bx_ambig,By_ambig

      Bx_ambig=c11*Bx(i,j)+c21*By(i,j)
      By_ambig=c12*Bx(i,j)+c22*By(i,j)

      if ((i.eq.nxa).and.(j.eq.nya)) then
         n_reconfig=1

         i_reconfig(1)=i
         j_reconfig(1)=j
         DivB_reconfig(1)=DivB(i,j)-ddxc(1)*Bx_ambig-ddyc(1)*By_ambig
         Jz_reconfig(1)=Jz(i,j)-ddxc(1)*By_ambig+ddyc(1)*Bx_ambig

         return
      endif

      if ((i.eq.nxb).and.(j.eq.nyb)) then
         n_reconfig=1
 
         i_reconfig(1)=i-1
         j_reconfig(1)=j-1
         DivB_reconfig(1)=DivB(i-1,j-1)-ddxc(4)*Bx_ambig-ddyc(4)*By_ambig
         Jz_reconfig(1)=Jz(i-1,j-1)-ddxc(4)*By_ambig+ddyc(4)*Bx_ambig

         return
      endif

      if ((i.eq.nxa).and.(j.eq.nyb)) then
         n_reconfig=1

         i_reconfig(1)=i
         j_reconfig(1)=j-1
         DivB_reconfig(1)=DivB(i,j-1)-ddxc(3)*Bx_ambig-ddyc(3)*By_ambig
         Jz_reconfig(1)=Jz(i,j-1)-ddxc(3)*By_ambig+ddyc(3)*Bx_ambig

         return
      endif

      if ((i.eq.nxb).and.(j.eq.nya)) then
         n_reconfig=1

         i_reconfig(1)=i-1
         j_reconfig(1)=j
         DivB_reconfig(1)=DivB(i-1,j)-ddxc(2)*Bx_ambig-ddyc(2)*By_ambig
         Jz_reconfig(1)=Jz(i-1,j)-ddxc(2)*By_ambig+ddyc(2)*Bx_ambig

         return
      endif

      if (j.eq.nya) then
         n_reconfig=2

         i_reconfig(1)=i
         j_reconfig(1)=j
         DivB_reconfig(1)=DivB(i,j)-ddxc(1)*Bx_ambig-ddyc(1)*By_ambig
         Jz_reconfig(1)=Jz(i,j)-ddxc(1)*By_ambig+ddyc(1)*Bx_ambig

         i_reconfig(2)=i-1
         j_reconfig(2)=j
         DivB_reconfig(2)=DivB(i-1,j)-ddxc(2)*Bx_ambig-ddyc(2)*By_ambig
         Jz_reconfig(2)=Jz(i-1,j)-ddxc(2)*By_ambig+ddyc(2)*Bx_ambig

         return
      endif

      if (i.eq.nxa) then
         n_reconfig=2

         i_reconfig(1)=i
         j_reconfig(1)=j
         DivB_reconfig(1)=DivB(i,j)-ddxc(1)*Bx_ambig-ddyc(1)*By_ambig
         Jz_reconfig(1)=Jz(i,j)-ddxc(1)*By_ambig+ddyc(1)*Bx_ambig

         i_reconfig(2)=i
         j_reconfig(2)=j-1
         DivB_reconfig(2)=DivB(i,j-1)-ddxc(3)*Bx_ambig-ddyc(3)*By_ambig
         Jz_reconfig(2)=Jz(i,j-1)-ddxc(3)*By_ambig+ddyc(3)*Bx_ambig

         return
      endif

      if (j.eq.nyb) then
         n_reconfig=2

         i_reconfig(1)=i
         j_reconfig(1)=j-1
         DivB_reconfig(1)=DivB(i,j-1)-ddxc(3)*Bx_ambig-ddyc(3)*By_ambig
         Jz_reconfig(1)=Jz(i,j-1)-ddxc(3)*By_ambig+ddyc(3)*Bx_ambig

         i_reconfig(2)=i-1
         j_reconfig(2)=j-1
         DivB_reconfig(2)=DivB(i-1,j-1)-ddxc(4)*Bx_ambig-ddyc(4)*By_ambig
         Jz_reconfig(2)=Jz(i-1,j-1)-ddxc(4)*By_ambig+ddyc(4)*Bx_ambig

         return
      endif

      if (i.eq.nxb) then
         n_reconfig=2

         i_reconfig(1)=i-1
         j_reconfig(1)=j
         DivB_reconfig(1)=DivB(i-1,j)-ddxc(2)*Bx_ambig-ddyc(2)*By_ambig
         Jz_reconfig(1)=Jz(i-1,j)-ddxc(2)*By_ambig+ddyc(2)*Bx_ambig

         i_reconfig(2)=i-1
         j_reconfig(2)=j-1
         DivB_reconfig(2)=DivB(i-1,j-1)-ddxc(4)*Bx_ambig-ddyc(4)*By_ambig
         Jz_reconfig(2)=Jz(i-1,j-1)-ddxc(4)*By_ambig+ddyc(4)*Bx_ambig

         return
      endif

      n_reconfig=4

      i_reconfig(1)=i
      j_reconfig(1)=j
      DivB_reconfig(1)=DivB(i,j)-ddxc(1)*Bx_ambig-ddyc(1)*By_ambig
      Jz_reconfig(1)=Jz(i,j)-ddxc(1)*By_ambig+ddyc(1)*Bx_ambig

      i_reconfig(2)=i
      j_reconfig(2)=j-1
      DivB_reconfig(2)=DivB(i,j-1)-ddxc(3)*Bx_ambig-ddyc(3)*By_ambig
      Jz_reconfig(2)=Jz(i,j-1)-ddxc(3)*By_ambig+ddyc(3)*Bx_ambig

      i_reconfig(3)=i-1
      j_reconfig(3)=j
      DivB_reconfig(3)=DivB(i-1,j)-ddxc(2)*Bx_ambig-ddyc(2)*By_ambig
      Jz_reconfig(3)=Jz(i-1,j)-ddxc(2)*By_ambig+ddyc(2)*Bx_ambig

      i_reconfig(4)=i-1
      j_reconfig(4)=j-1
      DivB_reconfig(4)=DivB(i-1,j-1)-ddxc(4)*Bx_ambig-ddyc(4)*By_ambig
      Jz_reconfig(4)=Jz(i-1,j-1)-ddxc(4)*By_ambig+ddyc(4)*Bx_ambig

      return
      end
c***********************************************************************************************************************************
