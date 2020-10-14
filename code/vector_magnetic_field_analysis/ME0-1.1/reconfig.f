c***********************************************************************************************************************************
      subroutine reconfig(i,j)
c===================================================================================================================================
c This subroutine implements the reconfiguration caused by changing the azimuthal angle at pixel (i,j). This involves:
c (1) reversing the sign of the transverse components of the field at pixel (i,j); and 
c (2) updating div(B) and Jz for each finite differencing stencil that is affected by the reconfiguration.
c
c 2009, Ashley Crouch, ash@cora.nwra.com
c===================================================================================================================================
      implicit none

      include "include_file.sizes.f"
      include "include_file.field.f"
      include "include_file.energy.f"
      include "include_file.reconfig.f"
      
      integer i,j,ir
c
c Reverse the sign of the transvserse components of the field
c
      Bx(i,j)=-Bx(i,j)
      By(i,j)=-By(i,j)
c
c Update div(B) and Jz
c
      do ir=1,n_reconfig
         DivB(i_reconfig(ir),j_reconfig(ir))=DivB_reconfig(ir)
         Jz(i_reconfig(ir),j_reconfig(ir))=Jz_reconfig(ir)
      enddo
      return
      end
c***********************************************************************************************************************************
