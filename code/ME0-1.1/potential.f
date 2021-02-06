c***********************************************************************
      subroutine potential(blpad,nxp,nyp,ixmin,ixmax,iymin,iymax)
c***********************************************************************
c     Use FFTW routines called through the Math Kernel Library to      *
c     compute the potential field and some of its gradients.           *
c***********************************************************************
      Use MKL_DFTI
      include "include_file.sizes.f"
      include "include_file.field.f"
      include "include_file.potential.f"
      include "include_file.point.f"
      include "include_file.pix_size.f"
c-----------------------------------------------------------------------
      integer i,j,i1d,im,jm,im1d
      integer ixmin,ixmax,iymin,iymax,nxp,nyp
      type(DFTI_DESCRIPTOR), POINTER :: Desc_Handle
      type(DFTI_DESCRIPTOR), POINTER :: Desc_Handle2
      integer   Status
      integer   lengths(2)
      integer   strides_in(3)
      integer   strides_out(3)
      real bz0,bl0,pi
      real bipx,bipy,db1dl,db2dl,db3dl
      real blpad(2*nxmax,2*nymax)
      real,dimension(:),allocatable :: xi,yi
      complex gradlos
      complex,dimension(:),allocatable :: blpad1d,fftbl1d
      complex,dimension(:),allocatable :: db1dl1d,db2dl1d,db3dl1d
      complex,dimension(:),allocatable :: bipx1d,bipy1d,bipz1d
      complex,dimension(:),allocatable :: fftbix1d,fftbiy1d,fftbiz1d
      complex,dimension(:),allocatable :: fftdb11d,fftdb21d,fftdb31d
c-----------------------------------------------------------------------
      data pi/3.1415926535897932/
c------------------------------------------------------------------------
c --> Subtract uniform vertical field to achieve flux balance.

      bl0=0.
      do 2 i=1,nxp
         do 1 j=1,nyp
            bl0=bl0+blpad(i,j)
    1    continue
    2 continue
      bl0=bl0/float(nxp)/float(nyp)
      bz0=bl0/a33

c --> Convert 2D real array of line of sight field to 1D complex array. 

      allocate(blpad1d(nxp*nyp),fftbl1d(nxp*nyp))
      do 4 i=1,nxp
         do 3 j=1,nyp
            i1d=i+(j-1)*nxp
            blpad1d(i1d)=cmplx(blpad(i,j)-bl0,0.)
    3    continue
    4 continue

c --> Fourier transform line of sight component of the field.

      lengths(1) = nxp
      lengths(2) = nyp

      strides_in(1) = 0
      strides_in(2) = 1
      strides_in(3) = nxp

      strides_out(1) = 0
      strides_out(2) = 1
      strides_out(3) = nxp

c --> Create MKL descriptor for 2D real to complex transform
      Status = DftiCreateDescriptor(Desc_Handle, DFTI_SINGLE,
     +                              DFTI_COMPLEX, 2, lengths)
      if (.not. DftiErrorClass(Status, DFTI_NO_ERROR)) then
          call dfti_example_status_print(Status)
          goto 999
      end if

      Status = DftiSetValue(Desc_Handle, DFTI_PLACEMENT,
     +                      DFTI_NOT_INPLACE)
      if (.not. DftiErrorClass(Status, DFTI_NO_ERROR)) then
          call dfti_example_status_print(Status)
          goto 999
      end if

      Status = DftiSetValue(Desc_Handle,DFTI_INPUT_STRIDES,strides_in)
      if (.not. DftiErrorClass(Status, DFTI_NO_ERROR)) then
          call dfti_example_status_print(Status)
          goto 999
      end if

      Status = DftiSetValue(Desc_Handle,DFTI_OUTPUT_STRIDES,strides_out)
      if (.not. DftiErrorClass(Status, DFTI_NO_ERROR)) then
          call dfti_example_status_print(Status)
          goto 999
      end if

c --> Set the scale
      sigma=1./(float(nxp)*float(nyp))
      Status = DftiSetValue(Desc_Handle, DFTI_FORWARD_SCALE, sigma)
      if (.not. DftiErrorClass(Status, DFTI_NO_ERROR)) then
          call dfti_example_status_print(Status)
          goto 999
      end if

      Status = DftiCommitDescriptor( Desc_Handle )
      if (.not. DftiErrorClass(Status, DFTI_NO_ERROR)) then
          call dfti_example_status_print(Status)
          goto 999
      end if

c --> Compute 2D complex to complex transform
      Status = DftiComputeForward(Desc_Handle,blpad1d,fftbl1d)
      if (.not. DftiErrorClass(Status, DFTI_NO_ERROR)) then
          call dfti_example_status_print(Status)
          goto 999
      end if

c --> Wavenumbers (recall the order in which FFTW returns the frequencies)
      allocate(xi(nxp),yi(nyp))

      do 5 i=1,nxp/2+1
         xi(i)=float(i-1)/float(nxp)/dxi
    5 continue
      do 6 i=1,(nxp-1)/2
         xi(nxp-i+1)=-xi(i+1)
    6 continue

      do 7 j=1,nyp/2+1
         yi(j)=float(j-1)/float(nyp)/dyi
    7 continue
      do 8 j=1,(nyp-1)/2
         yi(nyp-j+1)=-yi(j+1)
    8 continue

c --> Multiply by appropriate functions of wavenumber and coordinate
c --> transform. 
      allocate(fftbix1d(nxp*nyp),fftbiy1d(nxp*nyp),fftbiz1d(nxp*nyp))
      allocate(fftdb11d(nxp*nyp),fftdb21d(nxp*nyp),fftdb31d(nxp*nyp))

      do 10 i=1,nxp
         do 9 j=1,nyp
            i1d=i+(j-1)*nxp
            xkappa=2.*pi*sqrt(((c11**2+c12**2)*xi(i)**2+
     +                      2.*(c11*c21+c12*c22)*xi(i)*yi(j)+
     +                         (c21**2+c22**2)*yi(j)**2))
            gradlos=cmplx(-a33*xkappa,2.*pi*((c11*a13+c12*a23)*xi(i)+
     +                                       (c21*a13+c22*a23)*yi(j)))

            fftbix1d(i1d)=fftbl1d(i1d)*xi(i)/gradlos
            fftbiy1d(i1d)=fftbl1d(i1d)*yi(j)/gradlos
            fftbiz1d(i1d)=fftbl1d(i1d)*xkappa/gradlos

c --> Vertical derivative of helioplanar components.
            fftdb11d(i1d)=fftbl1d(i1d)*xi(i)*xkappa/gradlos
            fftdb21d(i1d)=fftbl1d(i1d)*yi(j)*xkappa/gradlos
            fftdb31d(i1d)=fftbl1d(i1d)*xkappa**2/gradlos

    9    continue
   10 continue
      deallocate(xi,yi)

c --> Force constant term to be 0.
      fftbix1d(1)=cmplx(0.,0.)
      fftbiy1d(1)=cmplx(0.,0.)
      fftbiz1d(1)=cmplx(0.,0.)

c --> Vertical derivatives also need constant term to be 0.
      fftdb11d(1)=cmplx(0.,0.)
      fftdb21d(1)=cmplx(0.,0.)
      fftdb31d(1)=cmplx(0.,0.)

c --> Inverse transform to get potential field components
      allocate(bipx1d(nxp*nyp),bipy1d(nxp*nyp),bipz1d(nxp*nyp))

c --> Compute 2D backward transform for x-component
      Status = DftiComputeBackward(Desc_Handle,fftbix1d,bipx1d)
      if (.not. DftiErrorClass(Status, DFTI_NO_ERROR)) then
          call dfti_example_status_print(Status)
          goto 999
      end if

c --> Compute 2D backward transform for y-component
      Status = DftiComputeBackward(Desc_Handle,fftbiy1d,bipy1d)
      if (.not. DftiErrorClass(Status, DFTI_NO_ERROR)) then
          call dfti_example_status_print(Status)
          goto 999
      end if

c --> Compute 2D backward transform for z-component
      Status = DftiComputeBackward(Desc_Handle,fftbiz1d,bipz1d)
      if (.not. DftiErrorClass(Status, DFTI_NO_ERROR)) then
          call dfti_example_status_print(Status)
          goto 999
      end if
      deallocate(fftbix1d,fftbiy1d,fftbiz1d)

c --> Inverse transform to get line of sight gradients of the
c --> potential field
      allocate(db1dl1d(nxp*nyp),db2dl1d(nxp*nyp),db3dl1d(nxp*nyp))

c --> Compute 2D backward transform for x-component
      Status = DftiComputeBackward(Desc_Handle,fftdb11d,db1dl1d)
      if (.not. DftiErrorClass(Status, DFTI_NO_ERROR)) then
          call dfti_example_status_print(Status)
          goto 999
      end if

c --> Compute 2D backward transform for y-component
      Status = DftiComputeBackward(Desc_Handle,fftdb21d,db2dl1d)
      if (.not. DftiErrorClass(Status, DFTI_NO_ERROR)) then
          call dfti_example_status_print(Status)
          goto 999
      end if

c --> Compute 2D backward transform for z-component
      Status = DftiComputeBackward(Desc_Handle,fftdb31d,db3dl1d)
      if (.not. DftiErrorClass(Status, DFTI_NO_ERROR)) then
          call dfti_example_status_print(Status)
          goto 999
      end if
      deallocate(fftdb11d,fftdb21d,fftdb31d)

c --> Extract field and line of sight gradients

      do 12 i=1,nx
         im=i+ixmin-1
         do 11 j=1,ny
            jm=j+iymin-1
            im1d=im+(jm-1)*nxp

            bipx=-2.*pi*aimag(bipx1d(im1d))
            bipy=-2.*pi*aimag(bipy1d(im1d))

            Bpx(i,j)=c11*bipx+c21*bipy
            Bpy(i,j)=c12*bipx+c22*bipy
            Bpz(i,j)=bz0-real(bipz1d(im1d))

            db1dl=2.*pi*aimag(db1dl1d(im1d))
            db2dl=2.*pi*aimag(db2dl1d(im1d))
            db3dl=real(db3dl1d(im1d))

c --> Vertical derivative of helioplanar components.
            dBxdz(i,j)=c11*db1dl+c21*db2dl
            dBydz(i,j)=c12*db1dl+c22*db2dl
            dBzdz(i,j)=db3dl

   11    continue
   12 continue
      deallocate(bipx1d,bipy1d,bipz1d)
      deallocate(db1dl1d,db2dl1d,db3dl1d)

  999 return
      end
