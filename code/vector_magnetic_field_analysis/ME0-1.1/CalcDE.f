c***********************************************************************************************************************************
      subroutine CalcDE(i,j,de)
c===================================================================================================================================
c This subroutine calculates the change in the "energy" caused by reconfiguring the azimuthal angle at pixel (i,j). Given Div(B) and
c Jz at the relevant finite differencing stencils for the reconfigured state, the change in the energy is determined by examining
c the appropriate sums over the affected stencils
c
c 2009, Ashley Crouch, ash@cora.nwra.com
c===================================================================================================================================
      implicit none

      include "include_file.sizes.f"
      include "include_file.energy.f"
      include "include_file.reconfig.f"
      include "include_file.lambda.f"

      integer i,j,ir
      real E,Eprev,de
c
c Compute Div(B) and Jz after the reconfiguration
c
      call CalcDE_reconfig(i,j)
c
c Calculate the sum of the contributions to the energy before and after the reconfiguration
c
      Eprev=0.
      E=0.
      do ir=1,n_reconfig
         Eprev=Eprev+abs(DivB(i_reconfig(ir),j_reconfig(ir)))+lambda*abs(Jz(i_reconfig(ir),j_reconfig(ir)))
         E=E+abs(DivB_reconfig(ir))+lambda*abs(Jz_reconfig(ir))
      enddo
c
c Calculate the change in energy due to the reconfiguration
c
      de=E-Eprev
      return
      end
c***********************************************************************************************************************************
