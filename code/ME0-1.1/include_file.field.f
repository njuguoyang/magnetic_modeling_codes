c***********************************************************************************************************************************
c must be preceeded by include_file.sizes.f
c***********************************************************************************************************************************
c
c The arrays containing the field related data
c
      real Bx(nxmax,nymax),By(nxmax,nymax),Bz(nxmax,nymax)
      real dBxdz(nxmax,nymax),dBydz(nxmax,nymax),dBzdz(nxmax,nymax)
      common /mgram_data/ Bx,By,Bz,dBxdz,dBydz,dBzdz
c***********************************************************************************************************************************
