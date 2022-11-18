module mod_usr   ! computing potential field 
  use mod_mhd
  use mod_lfff
  implicit none
  ! some global parameters
  double precision :: k_B,miu0,mass_H,usr_grav,radius_sun,rhob,Tiso,lalpha,llift
  logical, save :: firstusrglobaldata=.true.
  double precision, save :: x00,y00,phi00,lat00,B00

contains

  subroutine usr_init()
    use mod_global_parameters
    use mod_usr_methods

    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => initonegrid_usr
    usr_special_bc      => specialbound_usr
    usr_gravity         => gravity
    usr_refine_grid     => special_refine_grid
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output
    usr_special_convert => unstructuredvtkB_usr

    call set_coordinate_system("Cartesian_3D")
    call mhd_activate()
  end subroutine usr_init

  subroutine initglobaldata_usr()
    use mod_global_parameters

    ! normalization unit in CGS Unit
    k_B = 1.3806d-16          ! erg*K^-1
    miu0 = 4.d0*dpi           ! Gauss^2 cm^2 dyne^-1
    mass_H = 1.67262d-24      ! g
    unit_length        = 1.d9 ! cm
    unit_temperature   = 1.d6 ! K
    unit_numberdensity = 1.d9 ! cm^-3
    unit_density       = 1.4d0*mass_H*unit_numberdensity               ! 2.341668000000000E-015 g*cm^-3
    unit_pressure      = 2.3d0*unit_numberdensity*k_B*unit_temperature ! 0.317538000000000 erg*cm^-3
    unit_magneticfield = dsqrt(miu0*unit_pressure)                     ! 1.99757357615242 Gauss
    unit_velocity      = unit_magneticfield/dsqrt(miu0*unit_density)   ! 1.16448846777562E007 cm/s = 116.45 km/s
    unit_time          = unit_length/unit_velocity                     ! 85.8746159942810 s 

    usr_grav=-2.74d4*unit_length/unit_velocity**2  ! solar gravity
    radius_sun=6.955d10/unit_length                   ! Solar radius
    rhob=1.d0                                      ! bottom density
    Tiso=1.d6/unit_temperature                     ! uniform temperature

    lalpha=0.d0    ! alpha coefficient for linear force-free field
    llift=0.d0     ! lift up domain above bottom magnetogram
                   ! Usually, we need to lift up domain above bottom magnetogram by 0.5*dz, since the potential
                   ! field assumes the bottom boundary lies on the cell surface between the
                   ! ghost layers and the physcial domain, while we want the bottom boundary
                   ! lies on the cell center of the inner ghost layer.

    ! We want to convert the local Cartesian coordinates to a physical Cartesian
    ! coordiantes, whose x-axis is towards the observer, y-axis to the right, and
    ! z-axis to the upward. At the middle point, the bottom of the local Cartesian 
    ! box is tagent to the point (radius_sun, theta0, phi0) in the physical one.
    ! Local x, y, z are towards the westward, northward, and radial directions.
    x00 = half*(xprobmax1+xprobmin1)+0.0  ! to coalign with HMI line-of-sight magnetic field
    y00 = half*(xprobmax2+xprobmin2)+0.0
    ! Positions of SDO, STEREO A and B for 2011-06-21 01:11 UT
    !                                 STEREO-B           Earth        STEREO-A
    ! Heliographic (HEEQ) longitude    -92.472          -0.000          96.481
    ! Heliographic (HEEQ) latitude      -7.224           1.682           6.887

    ! viewed by SDO
    phi00 =  0.0973727d0              ! in the unit of radian. It is W5.57905
    lat00 =  0.2480280d0              ! lat is N14.2110 for this case
    B00   =  1.682d0*dpi/180.0d0      ! B0 angle of the soalr disk center
    ! viewed by STEREO-A
 !   phi00 = 0.0973727d0 - 96.481d0*dpi/180.0d0
 !   lat00 = 0.2480280d0
 !   B00   = 6.887d0*dpi/180.0d0
    ! viewed by STEREO-B 
 !   phi00 = 0.0973727d0 + 92.472d0*dpi/180.0d0
 !   lat00 = 0.2480280d0
 !   B00   =-7.224d0*dpi/180.0d0
 
    ! prepare magnetogram at bottom
    if(firstusrglobaldata) then
      call init_b_fff_data('potential_boundary/potential_boundary.dat',unit_length,unit_magneticfield)
      firstusrglobaldata=.false.
    endif
  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    ! initialize one grid
    use mod_global_parameters
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: Bf(ixI^S,1:ndir)
    logical, save :: first=.true.
    if(mype==0 .and. first) then
      write(*,*)'initializing grids ...'
      first=.false.
    endif
    ! use green function method to extrapolate B
    call calc_lin_fff(ixI^L,ixO^L,Bf,x,lalpha,llift)
    w(ixO^S,mag(:))=Bf(ixO^S,1:3)
    w(ixO^S,mom(:))=zero
    w(ixO^S,rho_)=rhob*dexp(usr_grav*radius_sun**2/Tiso*&
                  (1.d0/radius_sun-1.d0/(x(ixO^S,3)+radius_sun)))
  end subroutine initonegrid_usr

  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    ! special boundary types, user defined
    use mod_global_parameters
    use mod_physics
    integer, intent(in) :: ixO^L, iB, ixI^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: ft,tfstop,tramp1,tramp2,coeffrho,vlimit,vsign
    double precision :: delxdely,delxdelz,delydelx,delydelz,delzdelx,delzdely
    double precision :: xlen^D,dxa^D,startpos^D
    integer :: ix^D,ixIM^L,ixbc^D,af

    select case(iB)
    case(1)

    case(2)

    case(3)

    case(4)

    case(5)

    case(6)

    case default
      call mpistop("Special boundary is not defined for this region")
    end select
  end subroutine specialbound_usr

  subroutine getggrav(ggrid,ixI^L,ixO^L,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: ggrid(ixI^S)

    ggrid(ixO^S)=usr_grav*(radius_sun/(radius_sun+x(ixO^S,3)))**2
  end subroutine

  subroutine gravity(ixI^L,ixO^L,wCT,x,gravity_field)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(out)   :: gravity_field(ixI^S,ndim)
    double precision                :: ggrid(ixI^S)

    gravity_field=0.d0
    call getggrav(ggrid,ixI^L,ixO^L,x)
    gravity_field(ixO^S,3)=ggrid(ixO^S)
  end subroutine gravity

  subroutine special_refine_grid(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
  ! Enforce additional refinement or coarsening
  ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
    use mod_global_parameters

    integer, intent(in) :: igrid, level, ixI^L, ixO^L
    double precision, intent(in) :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
    integer, intent(inout) :: refine, coarsen

    ! fix the bottom layer to the highest level
    if (any(x(ixO^S,3)<=xprobmin3+0.3d0)) then
      refine=1
      coarsen=-1
    endif
  end subroutine special_refine_grid

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
    use mod_global_parameters
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)
    double precision                   :: tmp(ixI^S),dip(ixI^S),divb(ixI^S),B2(ixI^S)
    double precision, dimension(ixI^S,1:ndir) :: Btotal,qvec,curlvec
    integer                            :: ix^D,idirmin,idims,idir,jdir,kdir

    ! Btotal & B^2
    Btotal(ixI^S,1:ndir)=w(ixI^S,mag(1:ndir))
    B2(ixO^S)=sum((Btotal(ixO^S,:))**2,dim=ndim+1)
    ! output Alfven wave speed B/sqrt(rho)
    w(ixO^S,nw+1)=dsqrt(B2(ixO^S)/w(ixO^S,rho_))
    ! output divB1
    call divvector(Btotal,ixI^L,ixO^L,divb)
    w(ixO^S,nw+2)=0.5d0*divb(ixO^S)/dsqrt(B2(ixO^S))/(^D&1.0d0/dxlevel(^D)+)
    ! output the plasma beta p*2/B**2
    w(ixO^S,nw+3)=w(ixO^S,rho_)*Tiso*two/B2(ixO^S)
    ! store current
    call curlvector(Btotal,ixI^L,ixO^L,curlvec,idirmin,1,ndir)
    do idir=1,ndir
      w(ixO^S,nw+3+idir)=curlvec(ixO^S,idir)
    end do
    ! find magnetic dips
!    dip=0.d0
!    do idir=1,ndir
!      call gradient(w(ixI^S,mag(3)),ixI^L,ixO^L,idir,tmp)
!      dip(ixO^S)=dip(ixO^S)+w(ixO^S,b0_+idir)*tmp(ixO^S)
!    end do
!    where(dabs(w(ixO^S,mag(3)))<0.08d0 .and. dip(ixO^S)>=0.d0)
!      w(ixO^S,nw+8)=1.d0
!    elsewhere
!      w(ixO^S,nw+8)=0.d0
!    end where
  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    use mod_global_parameters
    character(len=*) :: varnames

    varnames='Alfv divB beta j1 j2 j3'
  end subroutine specialvarnames_output

  !==============================================================================
  ! Purpose: user output for vtu format to paraview
  !==============================================================================
subroutine unstructuredvtkB_usr(qunit)

! user output for vtu format to paraview, binary version output, cell_corner
! not parallel, uses calc_grid_user to compute nwauxio variables
! allows renormalizing using convert factors
use mod_forest, only: Morton_start, Morton_stop, sfc_to_igrid
use mod_global_parameters
use mod_calculate_xw

integer, intent(in) ::    qunit

double precision ::  x_VTK(1:3)

double precision, dimension(ixMlo^D-1:ixMhi^D,ndim) :: xC_TMP
double precision, dimension(ixMlo^D:ixMhi^D,ndim)   :: xCC_TMP
double precision, dimension(ixMlo^D-1:ixMhi^D,ndim) :: xC
double precision, dimension(ixMlo^D:ixMhi^D,ndim)   :: xCC

double precision, dimension(ixMlo^D-1:ixMhi^D,nw+nwauxio):: wC_TMP
double precision, dimension(ixMlo^D:ixMhi^D,nw+nwauxio)  :: wCC_TMP
double precision :: normconv(0:nw+nwauxio)

integer, allocatable :: intstatus(:,:)
integer :: itag,ipe,igrid,level,icel,ixC^L,ixCC^L,Morton_no,Morton_length
integer :: nx^D,nxC^D,nc,np,VTK_type,ix^D,filenr
integer*8 :: offset

integer::  k,iw
integer::  length,lengthcc,length_coords,length_conn,length_offsets
character::  buf
character(len=80)::  filename
character(len=name_len) :: wnamei(1:nw+nwauxio),xandwnamei(1:ndim+nw+nwauxio)
character(len=1024) :: outfilehead

logical ::   fileopen,cell_corner=.true.
logical, allocatable :: Morton_aim(:),Morton_aim_p(:)
!-----------------------------------------------------------------------------

normconv=one
Morton_length=Morton_stop(npe-1)-Morton_start(0)+1
allocate(Morton_aim(Morton_start(0):Morton_stop(npe-1)))
allocate(Morton_aim_p(Morton_start(0):Morton_stop(npe-1)))
Morton_aim=.false.
Morton_aim_p=.false.
do Morton_no=Morton_start(mype),Morton_stop(mype)
  igrid=sfc_to_igrid(Morton_no)
  level=node(plevel_,igrid)
  ! we can clip parts of the grid away, select variables, levels etc.
  if(writelevel(level)) then
   ! only output a grid when fully within clipped region selected
   ! by writespshift array
   if(({rnode(rpxmin^D_,igrid)>=xprobmin^D+(xprobmax^D-xprobmin^D)&
         *writespshift(^D,1)|.and.}).and.({rnode(rpxmax^D_,igrid)&
        <=xprobmax^D-(xprobmax^D-xprobmin^D)*writespshift(^D,2)|.and.})) then
     Morton_aim_p(Morton_no)=.true.
   end if
  end if
end do
call MPI_ALLREDUCE(Morton_aim_p,Morton_aim,Morton_length,MPI_LOGICAL,MPI_LOR,icomm,ierrmpi)
select case(convert_type)
 case('vtuB','vtuBmpi')
   cell_corner=.true.
 case('vtuBCC','vtuBCCmpi','user')
   cell_corner=.false.
end select
if (mype /= 0) then
 do Morton_no=Morton_start(mype),Morton_stop(mype)
   if(.not. Morton_aim(Morton_no)) cycle
   igrid=sfc_to_igrid(Morton_no)
   call calc_x_usr(igrid,xC,xCC)
   call calc_grid_usr(qunit,igrid,xC,xCC,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
                     ixC^L,ixCC^L,.true.)
   itag=Morton_no
   call MPI_SEND(xC_TMP,1,type_block_xc_io, 0,itag,icomm,ierrmpi)
   if(cell_corner) then
     call MPI_SEND(wC_TMP,1,type_block_wc_io, 0,itag,icomm,ierrmpi)
   else
     call MPI_SEND(wCC_TMP,1,type_block_wcc_io, 0,itag,icomm,ierrmpi)
   endif
 end do

else
 ! mype==0
 offset=0
 
 inquire(qunit,opened=fileopen)
 if(.not.fileopen)then
   ! generate filename 
    filenr=snapshotini
    if (autoconvert) filenr=snapshotnext
   write(filename,'(a,i4.4,a)') TRIM(base_filename),filenr,".vtu"
   ! Open the file for the header part
   open(qunit,file=filename,status='replace')
 endif
 
 call getheadernames(wnamei,xandwnamei,outfilehead)
 
 ! generate xml header
 write(qunit,'(a)')'<?xml version="1.0"?>'
 write(qunit,'(a)',advance='no') '<VTKFile type="UnstructuredGrid"'
 write(qunit,'(a)')' version="0.1" byte_order="LittleEndian">'
 write(qunit,'(a)')'<UnstructuredGrid>'
 write(qunit,'(a)')'<FieldData>'
 write(qunit,'(2a)')'<DataArray type="Float32" Name="TIME" ',&
                    'NumberOfTuples="1" format="ascii">'
 write(qunit,*) real(global_time*time_convert_factor)
 write(qunit,'(a)')'</DataArray>'
 write(qunit,'(a)')'</FieldData>'
 
 ! number of cells, number of corner points, per grid.
 nx^D=ixMhi^D-ixMlo^D+1;
 nxC^D=nx^D+1;
 nc={nx^D*}
 np={nxC^D*}
 
 length=np*size_real
 lengthcc=nc*size_real
 
 length_coords=3*length
 length_conn=2**^ND*size_int*nc
 length_offsets=nc*size_int

 ! Note: using the w_write, writelevel, writespshift
 do Morton_no=Morton_start(0),Morton_stop(0)
    if(.not. Morton_aim(Morton_no)) cycle
    if(cell_corner) then
      ! we write out every grid as one VTK PIECE
      write(qunit,'(a,i7,a,i7,a)') &
         '<Piece NumberOfPoints="',np,'" NumberOfCells="',nc,'">'
      write(qunit,'(a)')'<PointData>'
      do iw=1,nw+nwauxio
         if(iw<=nw) then 
            if(.not.w_write(iw)) cycle
         endif

         write(qunit,'(a,a,a,i16,a)')&
             '<DataArray type="Float32" Name="',TRIM(wnamei(iw)), &
             '" format="appended" offset="',offset,'">'
         write(qunit,'(a)')'</DataArray>'
         offset=offset+length+size_int
      enddo
      write(qunit,'(a)')'</PointData>'

      write(qunit,'(a)')'<Points>'
      write(qunit,'(a,i16,a)') &
   '<DataArray type="Float32" NumberOfComponents="3" format="appended" offset="',offset,'"/>'
      ! write cell corner coordinates in a backward dimensional loop, always 3D output
      offset=offset+length_coords+size_int
      write(qunit,'(a)')'</Points>'
    else
      ! we write out every grid as one VTK PIECE
      write(qunit,'(a,i7,a,i7,a)') &
         '<Piece NumberOfPoints="',np,'" NumberOfCells="',nc,'">'
      write(qunit,'(a)')'<CellData>'
      do iw=1,nw+nwauxio
         if(iw<=nw) then 
            if(.not.w_write(iw)) cycle
         endif

         write(qunit,'(a,a,a,i16,a)')&
             '<DataArray type="Float32" Name="',TRIM(wnamei(iw)), &
             '" format="appended" offset="',offset,'">'
         write(qunit,'(a)')'</DataArray>'
         offset=offset+lengthcc+size_int
      enddo
      write(qunit,'(a)')'</CellData>'

      write(qunit,'(a)')'<Points>'
      write(qunit,'(a,i16,a)') &
   '<DataArray type="Float32" NumberOfComponents="3" format="appended" offset="',offset,'"/>'
      ! write cell corner coordinates in a backward dimensional loop, always 3D output
      offset=offset+length_coords+size_int
      write(qunit,'(a)')'</Points>'
    end if
   
    write(qunit,'(a)')'<Cells>'

    ! connectivity part
    write(qunit,'(a,i16,a)')&
      '<DataArray type="Int32" Name="connectivity" format="appended" offset="',offset,'"/>'
    offset=offset+length_conn+size_int    

    ! offsets data array
    write(qunit,'(a,i16,a)') &
      '<DataArray type="Int32" Name="offsets" format="appended" offset="',offset,'"/>'
    offset=offset+length_offsets+size_int    

    ! VTK cell type data array
    write(qunit,'(a,i16,a)') &
      '<DataArray type="Int32" Name="types" format="appended" offset="',offset,'"/>' 
    offset=offset+size_int+nc*size_int

    write(qunit,'(a)')'</Cells>'

    write(qunit,'(a)')'</Piece>'
 end do
 ! write metadata communicated from other processors
 if(npe>1)then
  do ipe=1, npe-1
    do Morton_no=Morton_start(ipe),Morton_stop(ipe)
      if(.not. Morton_aim(Morton_no)) cycle
      if(cell_corner) then
        ! we write out every grid as one VTK PIECE
        write(qunit,'(a,i7,a,i7,a)') &
           '<Piece NumberOfPoints="',np,'" NumberOfCells="',nc,'">'
        write(qunit,'(a)')'<PointData>'
        do iw=1,nw+nwauxio
         if(iw<=nw) then 
            if(.not.w_write(iw)) cycle
         endif

           write(qunit,'(a,a,a,i16,a)')&
               '<DataArray type="Float32" Name="',TRIM(wnamei(iw)), &
               '" format="appended" offset="',offset,'">'
           write(qunit,'(a)')'</DataArray>'
           offset=offset+length+size_int
        enddo
        write(qunit,'(a)')'</PointData>'

        write(qunit,'(a)')'<Points>'
        write(qunit,'(a,i16,a)') &
     '<DataArray type="Float32" NumberOfComponents="3" format="appended" offset="',offset,'"/>'
        ! write cell corner coordinates in a backward dimensional loop, always 3D output
        offset=offset+length_coords+size_int
        write(qunit,'(a)')'</Points>'
      else
        ! we write out every grid as one VTK PIECE
        write(qunit,'(a,i7,a,i7,a)') &
           '<Piece NumberOfPoints="',np,'" NumberOfCells="',nc,'">'
        write(qunit,'(a)')'<CellData>'
        do iw=1,nw+nwauxio
         if(iw<=nw) then 
            if(.not.w_write(iw)) cycle
         endif

           write(qunit,'(a,a,a,i16,a)')&
               '<DataArray type="Float32" Name="',TRIM(wnamei(iw)), &
               '" format="appended" offset="',offset,'">'
           write(qunit,'(a)')'</DataArray>'
           offset=offset+lengthcc+size_int
        enddo
        write(qunit,'(a)')'</CellData>'

        write(qunit,'(a)')'<Points>'
        write(qunit,'(a,i16,a)') &
     '<DataArray type="Float32" NumberOfComponents="3" format="appended" offset="',offset,'"/>'
        ! write cell corner coordinates in a backward dimensional loop, always 3D output
        offset=offset+length_coords+size_int
        write(qunit,'(a)')'</Points>'
      end if
     
      write(qunit,'(a)')'<Cells>'

      ! connectivity part
      write(qunit,'(a,i16,a)')&
        '<DataArray type="Int32" Name="connectivity" format="appended" offset="',offset,'"/>'
      offset=offset+length_conn+size_int    

      ! offsets data array
      write(qunit,'(a,i16,a)') &
        '<DataArray type="Int32" Name="offsets" format="appended" offset="',offset,'"/>'
      offset=offset+length_offsets+size_int    

      ! VTK cell type data array
      write(qunit,'(a,i16,a)') &
        '<DataArray type="Int32" Name="types" format="appended" offset="',offset,'"/>' 
      offset=offset+size_int+nc*size_int

      write(qunit,'(a)')'</Cells>'

      write(qunit,'(a)')'</Piece>'
    end do
  end do
 end if

 write(qunit,'(a)')'</UnstructuredGrid>'
 write(qunit,'(a)')'<AppendedData encoding="raw">'
 close(qunit)
 open(qunit,file=filename,access='stream',form='unformatted',position='append')
 buf='_'
 write(qunit) TRIM(buf)

 do Morton_no=Morton_start(0),Morton_stop(0)
   if(.not. Morton_aim(Morton_no)) cycle
   igrid=sfc_to_igrid(Morton_no)
   call calc_x_usr(igrid,xC,xCC)
   call calc_grid_usr(qunit,igrid,xC,xCC,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
                  ixC^L,ixCC^L,.true.)
   do iw=1,nw+nwauxio
         if(iw<=nw) then 
            if(.not.w_write(iw)) cycle
         endif
     if(cell_corner) then
       write(qunit) length
       write(qunit) {(|}real(wC_TMP(ix^D,iw)*normconv(iw)),{ix^D=ixCmin^D,ixCmax^D)}
     else
       write(qunit) lengthcc
       write(qunit) {(|}real(wCC_TMP(ix^D,iw)*normconv(iw)),{ix^D=ixCCmin^D,ixCCmax^D)}
     end if
   enddo

   write(qunit) length_coords
   {do ix^DB=ixCmin^DB,ixCmax^DB \}
     x_VTK(1:3)=zero;
     x_VTK(1:ndim)=xC_TMP(ix^D,1:ndim)*normconv(0);
     do k=1,3
      write(qunit) real(x_VTK(k))
     end do
   {end do \}

   write(qunit) length_conn
   {do ix^DB=1,nx^DB\}
   {^IFONED write(qunit)ix1-1,ix1 \}
   {^IFTWOD
   write(qunit)(ix2-1)*nxC1+ix1-1, &
   (ix2-1)*nxC1+ix1,ix2*nxC1+ix1-1,ix2*nxC1+ix1
    \}
   {^IFTHREED
   write(qunit)&
   (ix3-1)*nxC2*nxC1+(ix2-1)*nxC1+ix1-1, &
   (ix3-1)*nxC2*nxC1+(ix2-1)*nxC1+ix1,&
   (ix3-1)*nxC2*nxC1+    ix2*nxC1+ix1-1,&
   (ix3-1)*nxC2*nxC1+    ix2*nxC1+ix1,&
    ix3*nxC2*nxC1+(ix2-1)*nxC1+ix1-1,&
    ix3*nxC2*nxC1+(ix2-1)*nxC1+ix1,&
    ix3*nxC2*nxC1+    ix2*nxC1+ix1-1,&
    ix3*nxC2*nxC1+    ix2*nxC1+ix1
    \}
   {end do\}

   write(qunit) length_offsets
   do icel=1,nc
     write(qunit) icel*(2**^ND)
   end do


  {^IFONED VTK_type=3 \}
  {^IFTWOD VTK_type=8 \}
  {^IFTHREED VTK_type=11 \}
   write(qunit) size_int*nc
   do icel=1,nc
     write(qunit) VTK_type
   end do
 end do
 allocate(intstatus(MPI_STATUS_SIZE,1))
 if(npe>1)then
  ixCCmin^D=ixMlo^D; ixCCmax^D=ixMhi^D;
  ixCmin^D=ixMlo^D-1; ixCmax^D=ixMhi^D;
  do ipe=1, npe-1
    do Morton_no=Morton_start(ipe),Morton_stop(ipe)
      if(.not. Morton_aim(Morton_no)) cycle
      itag=Morton_no
      call MPI_RECV(xC_TMP,1,type_block_xc_io, ipe,itag,icomm,intstatus(:,1),ierrmpi)
      if(cell_corner) then
        call MPI_RECV(wC_TMP,1,type_block_wc_io, ipe,itag,icomm,intstatus(:,1),ierrmpi)
      else
        call MPI_RECV(wCC_TMP,1,type_block_wcc_io, ipe,itag,icomm,intstatus(:,1),ierrmpi)
      end if
      do iw=1,nw+nwauxio
         if(iw<=nw) then 
            if(.not.w_write(iw)) cycle
         endif
        if(cell_corner) then
          write(qunit) length
          write(qunit) {(|}real(wC_TMP(ix^D,iw)*normconv(iw)),{ix^D=ixCmin^D,ixCmax^D)}
        else
          write(qunit) lengthcc
          write(qunit) {(|}real(wCC_TMP(ix^D,iw)*normconv(iw)),{ix^D=ixCCmin^D,ixCCmax^D)}
        end if
      enddo

      write(qunit) length_coords
      {do ix^DB=ixCmin^DB,ixCmax^DB \}
        x_VTK(1:3)=zero;
        x_VTK(1:ndim)=xC_TMP(ix^D,1:ndim)*normconv(0);
        do k=1,3
         write(qunit) real(x_VTK(k))
        end do
      {end do \}

      write(qunit) length_conn
      {do ix^DB=1,nx^DB\}
      {^IFONED write(qunit)ix1-1,ix1 \}
      {^IFTWOD
      write(qunit)(ix2-1)*nxC1+ix1-1, &
      (ix2-1)*nxC1+ix1,ix2*nxC1+ix1-1,ix2*nxC1+ix1
       \}
      {^IFTHREED
      write(qunit)&
      (ix3-1)*nxC2*nxC1+(ix2-1)*nxC1+ix1-1, &
      (ix3-1)*nxC2*nxC1+(ix2-1)*nxC1+ix1,&
      (ix3-1)*nxC2*nxC1+    ix2*nxC1+ix1-1,&
      (ix3-1)*nxC2*nxC1+    ix2*nxC1+ix1,&
       ix3*nxC2*nxC1+(ix2-1)*nxC1+ix1-1,&
       ix3*nxC2*nxC1+(ix2-1)*nxC1+ix1,&
       ix3*nxC2*nxC1+    ix2*nxC1+ix1-1,&
       ix3*nxC2*nxC1+    ix2*nxC1+ix1
       \}
      {end do\}

      write(qunit) length_offsets
      do icel=1,nc
        write(qunit) icel*(2**^ND)
      end do
      {^IFONED VTK_type=3 \}
      {^IFTWOD VTK_type=8 \}
      {^IFTHREED VTK_type=11 \}
      write(qunit) size_int*nc
      do icel=1,nc
        write(qunit) VTK_type
      end do
    end do
  end do
 end if
 close(qunit)
 open(qunit,file=filename,status='unknown',form='formatted',position='append')
 write(qunit,'(a)')'</AppendedData>'
 write(qunit,'(a)')'</VTKFile>'
 close(qunit)
 deallocate(intstatus)
end if

deallocate(Morton_aim,Morton_aim_p)
if (npe>1) then
  call MPI_BARRIER(icomm,ierrmpi)
endif

end subroutine unstructuredvtkB_usr

  !==============================================================================
  ! Purpose: 
  !==============================================================================
subroutine calc_grid_usr(qunit,igrid,xC,xCC,xC_TMP,xCC_TMP,wC_TMP,wCC_TMP,normconv,&
                     ixC^L,ixCC^L,first)

! this subroutine computes both corner as well as cell-centered values
! it handles how we do the center to corner averaging, as well as 
! whether we switch to cartesian or want primitive or conservative output,
! handling the addition of B0 in B0+B1 cases, ...
!
! the normconv is passed on to usr_aux_output for extending with
! possible normalization values for the nw+1:nw+nwauxio entries
use mod_usr_methods, only: usr_aux_output
use mod_global_parameters
use mod_limiter
use mod_physics, only: phys_energy, physics_type, phys_to_primitive
use mod_calculate_xw

integer, intent(in) :: qunit, igrid
double precision, intent(in), dimension(ixMlo^D-1:ixMhi^D,ndim) :: xC
double precision, intent(in), dimension(ixMlo^D:ixMhi^D,ndim)   :: xCC
integer :: ixC^L,ixCC^L
logical, intent(in) :: first

double precision, dimension(ixMlo^D-1:ixMhi^D,ndim) :: xC_TMP
double precision, dimension(ixMlo^D:ixMhi^D,ndim)   :: xCC_TMP
double precision, dimension(ixMlo^D-1:ixMhi^D,nw+nwauxio)   :: wC_TMP
double precision, dimension(ixMlo^D:ixMhi^D,nw+nwauxio)     :: wCC_TMP
double precision,dimension(0:nw+nwauxio),intent(out)       :: normconv 

double precision :: ldw(ixG^T), dwC(ixG^T)
double precision, dimension(ixMlo^D-1:ixMhi^D,nw+nwauxio)   :: wC
double precision, dimension(ixMlo^D:ixMhi^D,nw+nwauxio)     :: wCC
double precision, dimension(ixG^T,1:nw+nwauxio)   :: w
double precision :: dx^D
integer :: nxCC^D,idims,jxC^L,iwe
integer :: nx^D, nxC^D, ix^D, ix, iw, level, idir
logical, save :: subfirst=.true.
double precision   :: bxx,byy,bzz,bx1,by1,bz1,bx2,by2,bz2

!-----------------------------------------------------------------------------
ixCmin^D=ixMlo^D-1; ixCmax^D=ixMhi^D; ! Corner indices
ixCCmin^D=ixMlo^D; ixCCmax^D=ixMhi^D; ! Center indices

nx^D=ixMhi^D-ixMlo^D+1;
level=node(plevel_,igrid)
dx^D=dx(^D,level);

normconv(0) = length_convert_factor
normconv(1:nw) = w_convert_factor

w(ixG^T,1:nw)=ps(igrid)%w(ixG^T,1:nw)

if (nwextra>0) then
 ! here we actually fill the ghost layers for the nwextra variables using 
 ! continuous extrapolation (as these values do not exist normally in ghost cells)
 do idims=1,ndim
  select case(idims)
   {case(^D)
     jxCmin^DD=ixGhi^D+1-nghostcells^D%jxCmin^DD=ixGlo^DD;
     jxCmax^DD=ixGhi^DD;
     do ix^D=jxCmin^D,jxCmax^D
         w(ix^D^D%jxC^S,nw-nwextra+1:nw) = w(jxCmin^D-1^D%jxC^S,nw-nwextra+1:nw)
     end do 
     jxCmin^DD=ixGlo^DD;
     jxCmax^DD=ixGlo^D-1+nghostcells^D%jxCmax^DD=ixGhi^DD;
     do ix^D=jxCmin^D,jxCmax^D
         w(ix^D^D%jxC^S,nw-nwextra+1:nw) = w(jxCmax^D+1^D%jxC^S,nw-nwextra+1:nw)
     end do \}
  end select
 end do
end if

! next lines needed when usr_aux_output uses gradients
! and later on when dwlimiter2 is used 
^D&dxlevel(^D)=rnode(rpdx^D_,igrid);
if(nwauxio>0)then
  ! auxiliary io variables can be computed and added by user
  ! next few lines ensure correct usage of routines like divvector etc
  ! default (no) normalization for auxiliary variables
  normconv(nw+1:nw+nwauxio)=one

  if (.not. associated(usr_aux_output)) then
     call mpistop("usr_aux_output not defined")
  else
     call usr_aux_output(ixG^LL,ixM^LL^LADD1,w,ps(igrid)%x,normconv)
  end if
endif

! In case primitives to be saved: use primitive subroutine
!  extra layer around mesh only needed when storing corner values and averaging
if(saveprim.and.first) call phys_to_primitive(ixG^LL,ixM^LL^LADD1,w(ixG^T,1:nw),ps(igrid)%x)

if(B0field) then
! B0+B1 split handled here
  if(.not.saveprim.and.phys_energy) then
    w(ixG^T,iw_e)=w(ixG^T,iw_e)+0.5d0*sum(ps(igrid)%B0(ixG^T,:,0)**2,dim=ndim+1) &
          + sum(w(ixG^T,iw_mag(:))*ps(igrid)%B0(ixG^T,:,0),dim=ndim+1)
  end if
  w(ixG^T,iw_mag(:))=w(ixG^T,iw_mag(:))+ps(igrid)%B0(ixG^T,:,0)
end if

! rotate the magentic field vectors from the local Cartesian coordinates to a physical coordinate system
  do ix3=ixCCmin3,ixCCmax3
    do ix2=ixCCmin2,ixCCmax2
      do ix1=ixCCmin1,ixCCmax1
          bxx = w(ix1,ix2,ix3,mag(3))
          byy = w(ix1,ix2,ix3,mag(1))
          bzz = w(ix1,ix2,ix3,mag(2))
          bx1 = bxx*cos(lat00)-bzz*sin(lat00)
          by1 = byy
          bz1 = bxx*sin(lat00)+bzz*cos(lat00)

          bx2 = bx1*cos(phi00) - by1*sin(phi00)
          by2 = bx1*sin(phi00) + by1*cos(phi00)
          bz2 = bz1

          bx1 = bx2*cos(B00)+bz2*sin(B00)
          by1 = by2
          bz1 =-bx2*sin(B00)+bz2*cos(B00)
          
          w(ix1,ix2,ix3,mag(1))=bx1
          w(ix1,ix2,ix3,mag(2))=by1
          w(ix1,ix2,ix3,mag(3))=bz1
      end do
    end do
  end do

! compute the cell-center values for w first
! cell center values obtained from mere copy
wCC(ixCC^S,:)=w(ixCC^S,:)

! compute the corner values for w now by averaging

if(slab) then
   ! for slab symmetry: no geometrical info required
   do iw=1,nw+nwauxio
     {do ix^DB=ixCmin^DB,ixCmax^DB\}
        wC(ix^D,iw)=sum(w(ix^D:ix^D+1,iw))/dble(2**ndim)
     {end do\}
   end do
else
   do iw=1,nw+nwauxio
     {do ix^DB=ixCmin^DB,ixCmax^DB\}
       wC(ix^D,iw)=sum(w(ix^D:ix^D+1,iw)*ps(igrid)%dvolume(ix^D:ix^D+1)) &
                /sum(ps(igrid)%dvolume(ix^D:ix^D+1))
     {end do\}
   end do
endif

if(nocartesian) then
  ! keep the coordinate and vector components
  xC_TMP(ixC^S,1:ndim)          = xC(ixC^S,1:ndim)
  wC_TMP(ixC^S,1:nw+nwauxio)    = wC(ixC^S,1:nw+nwauxio)
  xCC_TMP(ixCC^S,1:ndim)        = xCC(ixCC^S,1:ndim)
  wCC_TMP(ixCC^S,1:nw+nwauxio)  = wCC(ixCC^S,1:nw+nwauxio)
else
  ! do all conversions to cartesian coordinates and vector components
  ! start for the corner values
  call to_cartesian(xC_TMP,wC_TMP,ixC^L,xC,wC)
  ! then cell center values
  call to_cartesian(xCC_TMP,wCC_TMP,ixCC^L,xCC,wCC)
endif

! Warning: differentiate between idl/idlCC/tecplot/tecplotCC/vtu(B)/vtu(B)CC
if(nwaux>0 .and. mype==0 .and. first.and.subfirst) then
  ! when corner values are computed and auxiliaries present: warn!
  if(convert_type=='idl'.or.convert_type=='tecplot' &
   .or.convert_type=='vtu'.or.convert_type=='vtuB') &
      write(*,*) 'Warning: also averaged auxiliaries within calc_grid'
  subfirst=.false.
endif

end subroutine calc_grid_usr

  !==============================================================================
  ! Purpose: 
  !==============================================================================
  subroutine calc_x_usr(igrid,xC,xCC)
  use mod_global_parameters

  integer, intent(in)               :: igrid
  double precision, intent(out)     :: xC(ixMlo^D-1:ixMhi^D,ndim)
  double precision, intent(out)     :: xCC(ixMlo^D:ixMhi^D,ndim)
  ! .. local ..
  integer                           :: ixC^L, ixCC^L, idims, level, ix^D
  double precision   :: phi11,theta11,r11,rr,x11,y11,z11,x22,y22,z22

  level=node(plevel_,igrid)

  ! coordinates of cell centers
  xCC(ixM^T,:)=ps(igrid)%x(ixM^T,:)

  ! We want to convert the local Cartesian coordinates to a physical Cartesian
  ! coordiantes, whose x-axis is towards the observer, y-axis to the right, and
  ! z-axis to the upward. At the middle point, the bottom of the local Cartesian 
  ! box is tagent to the point (radius_sun, theta0, phi0) in the physical one.
  ! Local x, y, z are towards the westward, northward, and radial directions.

  ixCCmin^D=ixMlo^D; ixCCmax^D=ixMhi^D; ! coordinates of cell centers
  do ix3=ixCCmin3,ixCCmax3
    do ix2=ixCCmin2,ixCCmax2
      do ix1=ixCCmin1,ixCCmax1
       x11 = xCC(ix1,ix2,ix3,1)
       y11 = xCC(ix1,ix2,ix3,2)
       z11 = xCC(ix1,ix2,ix3,3)
       phi11 = atan((x11-x00)/(radius_sun+z11))
       r11 = sqrt((radius_sun+z11)**2 + (x11-x00)**2)
       theta11 = atan((y11-y00)/r11)
       rr = sqrt(r11**2 + (y11-y00)**2)
       x22 = rr*cos(theta11)*cos(phi11)
       y22 = rr*cos(theta11)*sin(phi11)
       z22 = rr*sin(theta11)

       x11 = cos(lat00)*x22 - sin(lat00)*z22
       y11 = y22
       z11 = sin(lat00)*x22 + cos(lat00)*z22       

       x22 = cos(phi00)*x11 - sin(phi00)*y11
       y22 = sin(phi00)*x11 + cos(phi00)*y11
       z22 = z11

       x11 = cos(B00)*x22 + sin(B00)*z22
       y11 = y22
       z11 =-sin(B00)*x22 + cos(B00)*z22 
       xCC(ix1,ix2,ix3,1) = x11
       xCC(ix1,ix2,ix3,2) = y11
       xCC(ix1,ix2,ix3,3) = z11
      end do
    end do
  end do

  ! coordinates of cell corners
  ixCmin^D=ixMlo^D-1; ixCmax^D=ixMhi^D;
  if(slab)then
     !do idims=1,ndim
     !  xC(ixC^S,idims)=ps(igrid)%x(ixC^S,idims)+0.5d0*dx(idims,level)
     !end do
    do ix3=ixCmin3,ixCmax3
     do ix2=ixCmin2,ixCmax2
      do ix1=ixCmin1,ixCmax1
       x11 = ps(igrid)%x(ix^D,1)+0.5d0*dx(1,level)
       y11 = ps(igrid)%x(ix^D,2)+0.5d0*dx(2,level)
       z11 = ps(igrid)%x(ix^D,3)+0.5d0*dx(3,level)
       phi11 = atan((x11-x00)/(radius_sun+z11))
       r11 = sqrt((radius_sun+z11)**2 + (x11-x00)**2)
       theta11 = atan((y11-y00)/r11)
       rr = sqrt(r11**2 + (y11-y00)**2)
       x22 = rr*cos(theta11)*cos(phi11)
       y22 = rr*cos(theta11)*sin(phi11)
       z22 = rr*sin(theta11)

       x11 = cos(lat00)*x22 - sin(lat00)*z22
       y11 = y22
       z11 = sin(lat00)*x22 + cos(lat00)*z22       

       x22 = cos(phi00)*x11 - sin(phi00)*y11
       y22 = sin(phi00)*x11 + cos(phi00)*y11
       z22 = z11

       x11 = cos(B00)*x22 + sin(B00)*z22
       y11 = y22
       z11 =-sin(B00)*x22 + cos(B00)*z22 
       xC(ix1,ix2,ix3,1) = x11
       xC(ix1,ix2,ix3,2) = y11
       xC(ix1,ix2,ix3,3) = z11
      end do
     end do
    end do
  else
     ! for any non-cartesian or stretched coordinate (allow multiple stretched directions)
    do ix3=ixCmin3,ixCmax3
     do ix2=ixCmin2,ixCmax2
      do ix1=ixCmin1,ixCmax1
       x11 = ps(igrid)%x(ix^D,1)+0.5d0*ps(igrid)%dx(ix^D,1)
       y11 = ps(igrid)%x(ix^D,2)+0.5d0*ps(igrid)%dx(ix^D,2)
       z11 = ps(igrid)%x(ix^D,3)+0.5d0*ps(igrid)%dx(ix^D,3)
       phi11 = atan((x11-x00)/(radius_sun+z11))
       r11 = sqrt((radius_sun+z11)**2 + (x11-x00)**2)
       theta11 = atan((y11-y00)/r11)
       rr = sqrt(r11**2 + (y11-y00)**2)
       x22 = rr*cos(theta11)*cos(phi11)
       y22 = rr*cos(theta11)*sin(phi11)
       z22 = rr*sin(theta11)

       x11 = cos(lat00)*x22 - sin(lat00)*z22
       y11 = y22
       z11 = sin(lat00)*x22 + cos(lat00)*z22       

       x22 = cos(phi00)*x11 - sin(phi00)*y11
       y22 = sin(phi00)*x11 + cos(phi00)*y11
       z22 = z11

       x11 = cos(B00)*x22 + sin(B00)*z22
       y11 = y22
       z11 =-sin(B00)*x22 + cos(B00)*z22 
       xC(ix1,ix2,ix3,1) = x11
       xC(ix1,ix2,ix3,2) = y11
       xC(ix1,ix2,ix3,3) = z11
      end do
     end do
    end do
  endif

  end subroutine calc_x_usr


end module mod_usr





