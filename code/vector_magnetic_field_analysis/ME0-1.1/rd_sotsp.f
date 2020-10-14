c***********************************************************************
      subroutine rd_sotsp(filename,npad)
c***********************************************************************
c     Reads in magnetic field data plus pointing information from a    *
c     FITS file, in the format used by Hinode SOT/SP from NAOJ. Based  *
c     on routines given in cookbook.f supplied with FITSIO.            *
c***********************************************************************
      include "include_file.sizes.f"
      include "include_file.point.f"
      include "include_file.field.f"
      include "include_file.fits.f"
      include "include_file.flags.f"
      include "include_file.pix_size.f"
c-----------------------------------------------------------------------
      integer i,j,i1d
      integer status,unit,readwrite,blocksize
      integer firstpix,group,nfound,narray
      integer colnum,irow,felem,nelems,naxis,maxdim
      integer datacode,repeat,width
      integer nkeys,nspace,hdutype,k
      integer nlimb
      real nullval,pi,dtor,phase,radius,xcen,ycen
      real cifrac,cimean
      real,dimension(nxmax,nymax) :: bl,bt,ba
      real,dimension(nxmax*nymax) :: bs1d,bi1d,ba1d,sl1d,ci1d
      logical anynull
      character*90 filename
      character*80 record,comment
      common/bobs/ bl,bt,ba
c-----------------------------------------------------------------------
      data pi,radius,cifrac/3.1415926535897932,959.689,0.05/
c-----------------------------------------------------------------------
c --> The status parameter must always be initialized.
      status=0

c --> Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)

c --> Open the FITS file, with read-only access.  The returned blocksize
c --> parameter is obsolete and should be ignored. 
      readwrite=0
      call ftopen(unit,filename,readwrite,blocksize,status)

c --> Get the pixel size from the keywords.

      call ftgkye(unit,'XSCALE',dxi,comment,status)
      call ftgkye(unit,'YSCALE',dyi,comment,status)

c --> Get the location of the FOV relative to disk center. 

      call ftgkye(unit,'XCEN',xcen,comment,status)
      call ftgkye(unit,'YCEN',ycen,comment,status)

c --> Determine latitude and longitude from this. 
c --> theta=central meridian angle
c --> phi=latitude
      theta=atan(xcen/sqrt(radius**2-xcen**2-ycen**2))
      phi=asin(ycen/radius)

c --> Zero p and b angles. 
      p=0.
      b=0.

c --> Jump straight to the HDU containing the field strength.

      call ftmahd(unit,jfs,hdutype,status)

c --> Check this HDU exists.

      if(status.ne.0) then
        write(*,*) 'error encountered while looking for field strength'
        call printerror(status)
        stop
      endif

c --> Have to treat this differently depending on whether it's an image 
c --> extension (hdutype=0) or a binary table extension (hdutype=2).

      if(hdutype.eq.0) then

c --> Check this is the field strength.
         call ftgkys(unit,'PAR_TYPE',record,comment,status)
         if(record.ne.'Field strength') then
            write(*,*) 'unexpected PAR_TYPE for field strength',record
         endif

c --> Determine the size of the extension.
         call ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)
         if(nfound.ne.2) then
            write(*,*) '*** Unexpected number of axes found while 
     +                  reading field strength ***',nfound
            write(*,*) naxes(1),naxes(2)
         endif

c --> Initialize variables for reading field strength. 
         narray=naxes(1)*naxes(2)
         group=1
         nullval=-999.
         firstpix=1

c --> Read in field strength.
         call ftgpve(unit,group,firstpix,narray,nullval,bs1d,anynull,
     &               status)

      elseif(hdutype.eq.2) then

c --> Check this is the field strength.
         call ftgkys(unit,'EXTNAME',record,comment,status)
         if(record.ne.'FieldStrength') then
            write(*,*) 'unexpected EXTNAME for field strength',record
         endif

c --> Determine the size of the extension.
         maxdim=2
         colnum=1
         call ftgtcl(unit,colnum,datacode,repeat,width,status)
         call printerror(status)

         call ftgtdm(unit,colnum,maxdim,naxis,naxes,status)
         call printerror(status)

         if(naxis.ne.2) then
            write(*,*) '*** Unexpected number of axes found while
     + reading field strength ***',naxis
            write(*,*) naxes(1),naxes(2)
         endif

         if(repeat.ne.naxes(1)*naxes(2)) then
            write(*,*) '*** Inconsistent row length found while
     + reading field strength ***'
            write(*,*) repeat,naxes(1),naxes(2)
         endif

c --> Initialize variables for reading field strength. 
         irow=1
         felem=1
         nelems=repeat

c --> Read in field strength.
         if(datacode.eq.42) then
            call ftgcve(unit,colnum,irow,felem,nelems,nullval,bs1d,
     &                  anynull,status)
         else
            write(*,*) 'unexpected data type for field strength',
     +                 datacode
            stop
         endif

      else
         write(*,*) 'unrecognized HDU type for field strength',hdutype
         stop
      endif

      nx=naxes(1)
      ny=naxes(2)

c --> Repeat for the HDU containing the field inclination.

      call ftmahd(unit,jfi,hdutype,status)

c --> Check this HDU exists.

      if(status.ne.0) then
        write(*,*) 'error encountered while looking for field
     +  inclination'
        call printerror(status)
        stop
      endif

c --> Have to treat this differently depending on whether it's an image 
c --> extension (hdutype=0) or a binary table extension (hdutype=2).

      if(hdutype.eq.0) then

c --> Check this is the field inclination.
         call ftgkys(unit,'PAR_TYPE',record,comment,status)
         if(record.ne.'Field inclination') then
           write(*,*) 'unexpected PAR_TYPE for field inclination',record
         endif

c --> Determine the size of the image.
         call ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)
         if(nfound.ne.2) stop '*** Unexpected number of axes found ***'

c --> Initialize variables for reading field inclination. 
         narray=naxes(1)*naxes(2)
         group=1
         nullval=-999.
         firstpix=1

c --> Read in field inclination.
         call ftgpve(unit,group,firstpix,narray,nullval,bi1d,anynull,
     &               status)

      elseif(hdutype.eq.2) then

c --> Check this is the field inclination.
         call ftgkys(unit,'EXTNAME',record,comment,status)
         if(record.ne.'FieldInclination') then
            write(*,*) 'unexpected EXTNAME for field inclination',record
         endif

c --> Determine the size of the extension.
         maxdim=2
         colnum=1
         call ftgtcl(unit,colnum,datacode,repeat,width,status)
         call printerror(status)

         call ftgtdm(unit,colnum,maxdim,naxis,naxes,status)
         call printerror(status)

         if(naxis.ne.2) then
            write(*,*) '*** Unexpected number of axes found while
     + reading field inclination ***',naxis
            write(*,*) naxes(1),naxes(2)
         endif

         if(repeat.ne.naxes(1)*naxes(2)) then
            write(*,*) '*** Inconsistent row length found while
     + reading field inclination ***'
            write(*,*) repeat,naxes(1),naxes(2)
         endif

c --> Initialize variables for reading field inclination. 
         irow=1
         felem=1
         nelems=repeat

c --> Read in field inclination.
         if(datacode.eq.42) then
            call ftgcve(unit,colnum,irow,felem,nelems,nullval,bi1d,
     &                  anynull,status)
         else
            write(*,*) 'unexpected data type for field inclination',
     +                 datacode
            stop
         endif

      else
        write(*,*) 'unrecognized HDU type for field inclination',hdutype
        stop
      endif

c --> Check the dimensions of this array match the previous dimensions.
      if(nx.ne.naxes(1).or.ny.ne.naxes(2)) then
        write(*,*) 'array dimensions do not match for field inclination'
         write(*,*) 'nx,ny',ny,ny
         write(*,*) 'naxes',naxes(1),naxes(2)
         stop
      endif

c --> Repeat for the HDU containing the field azimuth.

      call ftmahd(unit,jfa,hdutype,status)

c --> Check this HDU exists.

      if(status.ne.0) then
        write(*,*) 'error encountered while looking for field azimuth'
        call printerror(status)
        stop
      endif

c --> Have to treat this differently depending on whether it's an image 
c --> extension (hdutype=0) or a binary table extension (hdutype=2).

      if(hdutype.eq.0) then

c --> Check this is the field azimuth.
         call ftgkys(unit,'PAR_TYPE',record,comment,status)
         if(record.ne.'Field azimuth') then
           write(*,*) 'unexpected PAR_TYPE for field azimuth',record
         endif

c --> Determine the size of the image.
         call ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)
         if(nfound.ne.2) stop '*** Unexpected number of axes found while
     + reading field azimuth ***'

c --> Initialize variables for reading field azimuth. 
         narray=naxes(1)*naxes(2)
         group=1
         nullval=-999.
         firstpix=1

c --> Read in field azimuth.
         call ftgpve(unit,group,firstpix,narray,nullval,ba1d,anynull,
     &               status)

      elseif(hdutype.eq.2) then

c --> Check this is the field azimuth.
         call ftgkys(unit,'EXTNAME',record,comment,status)
         if(record.ne.'FieldAzimuth') then
            write(*,*) 'unexpected EXTNAME for field azimuth',record
         endif

c --> Determine the size of the extension.
         maxdim=2
         colnum=1
         call ftgtcl(unit,colnum,datacode,repeat,width,status)
         call printerror(status)

         call ftgtdm(unit,colnum,maxdim,naxis,naxes,status)
         call printerror(status)

         if(naxis.ne.2) then
            write(*,*) '*** Unexpected number of axes found while
     + reading field azimuth ***',naxis
            write(*,*) naxes(1),naxes(2)
         endif

         if(repeat.ne.naxes(1)*naxes(2)) then
            write(*,*) '*** Inconsistent row length found while
     + reading field azimuth ***'
            write(*,*) repeat,naxes(1),naxes(2)
         endif

c --> Initialize variables for reading field azimuth. 
         irow=1
         felem=1
         nelems=repeat

c --> Read in field azimuth.
         if(datacode.eq.42) then
            call ftgcve(unit,colnum,irow,felem,nelems,nullval,ba1d,
     &                  anynull,status)
         else
            write(*,*) 'unexpected data type for field azimuth',
     +                 datacode
            stop
         endif

      else
        write(*,*) 'unrecognized HDU type for field azimuth',hdutype
        stop
      endif

c --> Check the dimensions of this array match the previous dimensions.
      if(nx.ne.naxes(1).or.ny.ne.naxes(2)) then
        write(*,*) 'array dimensions do not match for field azimuth'
         write(*,*) 'nx,ny',ny,ny
         write(*,*) 'naxes',naxes(1),naxes(2)
         stop
      endif

c --> Repeat for the HDU containing the straylight fraction.

      call ftmahd(unit,jsf,hdutype,status)

c --> Check this HDU exists.

      if(status.ne.0) then
        write(*,*) 'error encountered while looking for straylight
     +  fraction'
        call printerror(status)
        stop
      endif

c --> Have to treat this differently depending on whether it's an image 
c --> extension (hdutype=0) or a binary table extension (hdutype=2).

      if(hdutype.eq.0) then

c --> Check this is the straylight fraction.
         call ftgkys(unit,'PAR_TYPE',record,comment,status)
         if(record.ne.'Straylight fraction') then
            write(*,*) 'unexpected PAR_TYPE for Straylight fraction',
     +                 record
         endif

c --> Determine the size of the image.
         call ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)
         if(nfound.ne.2) stop '*** Unexpected number of axes found while
     + reading Staylight fraction ***'

c --> Initialize variables for reading straylight fraction. 
         narray=naxes(1)*naxes(2)
         group=1
         nullval=-999.
         firstpix=1

c --> Read in straylight fraction.
         call ftgpve(unit,group,firstpix,narray,nullval,sl1d,anynull,
     &               status)

      elseif(hdutype.eq.2) then

c --> Check this is the stray light fill factor.
         call ftgkys(unit,'EXTNAME',record,comment,status)
         if(record.ne.'Alpha') then
            write(*,*) 'unexpected EXTNAME for Stray Light Fill Factor',
     +                 record
         endif

c --> Determine the size of the extension.
         maxdim=2
         colnum=1
         call ftgtcl(unit,colnum,datacode,repeat,width,status)
         call printerror(status)

         call ftgtdm(unit,colnum,maxdim,naxis,naxes,status)
         call printerror(status)

         if(naxis.ne.2) then
            write(*,*) '*** Unexpected number of axes found while
     + reading stray light fill factor ***',naxis
            write(*,*) naxes(1),naxes(2)
         endif

         if(repeat.ne.naxes(1)*naxes(2)) then
            write(*,*) '*** Inconsistent row length found while
     + reading stray light fill factor ***'
            write(*,*) repeat,naxes(1),naxes(2)
         endif

c --> Initialize variables for reading stray light fill factor. 
         irow=1
         felem=1
         nelems=repeat

c --> Read in stray light fill factor.
         if(datacode.eq.42) then
            call ftgcve(unit,colnum,irow,felem,nelems,nullval,sl1d,
     &                  anynull,status)
         else
            write(*,*) 'unexpected data type for stray light fill
     + factor',datacode
            stop
         endif

      else
        write(*,*) 'unrecognized HDU type for Straylight fraction',
     +             hdutype
        stop
      endif

c --> Check the dimensions of this array match the previous dimensions.
      if(nx.ne.naxes(1).or.ny.ne.naxes(2)) then
         write(*,*) 'array dimensions do not match for straylight
     + fraction'
         write(*,*) 'nx,ny',ny,ny
         write(*,*) 'naxes',naxes(1),naxes(2)
         stop
      endif

c --> Repeat for the HDU containing the continuum intensity.

      call ftmahd(unit,jci,hdutype,status)

c --> Check this HDU exists.

      if(status.ne.0) then
        write(*,*) 'error encountered while looking for continuum
     +  intensity'
        call printerror(status)
        stop
      endif

c --> Have to treat this differently depending on whether it's an image 
c --> extension (hdutype=0) or a binary table extension (hdutype=2).

      if(hdutype.eq.0) then

c --> Check this is the continuum intensity.
         call ftgkys(unit,'PAR_TYPE',record,comment,status)
         if(record.ne.'Continuum intensity') then
            write(*,*) 'unexpected PAR_TYPE for Continuum intensity',
     +                 record
         endif

c --> Determine the size of the image.
         call ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)
         if(nfound.ne.2) stop '*** Unexpected number of axes found while
     + reading continuum intensity ***'

c --> Initialize variables for reading continuum intensity. 
         narray=naxes(1)*naxes(2)
         group=1
         nullval=-999.
         firstpix=1

c --> Read in continuum intensity. 
         call ftgpve(unit,group,firstpix,narray,nullval,ci1d,anynull,
     &               status)

      elseif(hdutype.eq.2) then

c --> Check this is the line strength.
         call ftgkys(unit,'EXTNAME',record,comment,status)
c --> Note typo: this is what appears to be in the level 2 data.
         if(record.ne.'LineStength') then
            write(*,*) 'unexpected EXTNAME for Line Strength',
     +                 record
         endif

c --> Determine the size of the extension.
         maxdim=2
         colnum=1
         call ftgtcl(unit,colnum,datacode,repeat,width,status)
         call printerror(status)

         call ftgtdm(unit,colnum,maxdim,naxis,naxes,status)
         call printerror(status)

         if(naxis.ne.2) then
            write(*,*) '*** Unexpected number of axes found while
     + reading line strength ***',naxis
            write(*,*) naxes(1),naxes(2)
         endif

         if(repeat.ne.naxes(1)*naxes(2)) then
            write(*,*) '*** Inconsistent row length found while
     + reading line strength ***'
            write(*,*) repeat,naxes(1),naxes(2)
         endif

c --> Initialize variables for reading line strength. 
         irow=1
         felem=1
         nelems=repeat

c --> Read in line strength.
         if(datacode.eq.42) then
            call ftgcve(unit,colnum,irow,felem,nelems,nullval,ci1d,
     &                  anynull,status)
         else
            write(*,*) 'unexpected data type for line strength',datacode
            stop
         endif

      else
        write(*,*) 'unrecognized HDU type for Continuum intensity',
     +             hdutype
        stop
      endif

c --> Check the dimensions of this array match the previous dimensions.
      if(nx.ne.naxes(1).or.ny.ne.naxes(2)) then
         write(*,*) 'array dimensions do not match for continuum
     + intensity'
         write(*,*) 'nx,ny',ny,ny
         write(*,*) 'naxes',naxes(1),naxes(2)
         stop
      endif

c --> Set a continuum intensity to define off-limb observations. 
c --> Default is less than 5% of the mean continuum intensity. 
      cimean=0.
      do j=1,nxmax*nymax
         cimean=cimean+ci1d(j)
      enddo
      cimean=cifrac*cimean/float(nx*ny)

c --> Convert field strength to mean field ("flux") and zero the field
c --> strength off the limb.
      nlimb=nx*ny-nxmax*nymax
      do j=1,nxmax*nymax
         if(ci1d(j).gt.cimean) then
            bs1d(j)=bs1d(j)*(1.0-sl1d(j))
         else
            bs1d(j)=0.
            nlimb=nlimb+1
         endif
      enddo

c --> Check that various array dimensions are consistent.
      if(nx.gt.nxmax.or.ny.gt.nymax) then
         write(*,*) 'Increase array dimensions in include_file.sizes.f'
         write(*,*) 'nxmax=',nxmax,' must be greater than nx=',nx
         write(*,*) 'nymax=',nymax,' must be greater than ny=',ny
         stop
      endif

      if(nx+2*npad.gt.2*nxmax.or.ny+2*npad.gt.2*nymax) then
         write(*,*) 'Increase array dimensions in include_file.sizes.f'
         write(*,*) 'or decrease padding in file par'
         write(*,*) '2 nxmax=',2*nxmax,' must be greater than nx+2npad='
     +               ,nx+2*npad
         write(*,*) '2 nymax=',2*nymax,' must be greater than ny+2npad='
     +               ,ny+2*npad
         stop
      endif

      if(status.gt.0) call printerror(status)

      phase=0.5*float(iaflag)*pi

c --> Convert to radians if necessary. 
      if(iaunit.eq.1) then
         dtor=pi/180.
         do 8 i=1,nx*ny
            bi1d(i)=bi1d(i)*dtor
            ba1d(i)=ba1d(i)*dtor
    8    continue
      endif

c --> Transform to two components. 
      if(incflag.eq.0) then
         do 10 i=1,nx
            do 9 j=1,ny
               i1d=i+(j-1)*nx
               bl(i,j)=bs1d(i1d)*cos(bi1d(i1d))
               bt(i,j)=bs1d(i1d)*sin(bi1d(i1d))
               ba(i,j)=ba1d(i1d)
c --> Calculate Cartesian (image) components.
               Bx(i,j)=bt(i,j)*cos(ba(i,j)+phase)
               By(i,j)=bt(i,j)*sin(ba(i,j)+phase)
               Bz(i,j)=bl(i,j)
    9       continue
   10    continue
      else
         do 12 i=1,nx
            do 11 j=1,ny
               i1d=i+(j-1)*nx
               bl(i,j)=bs1d(i1d)*sin(bi1d(i1d))
               bt(i,j)=bs1d(i1d)*cos(bi1d(i1d))
               ba(i,j)=ba1d(i1d)
c --> Calculate Cartesian (image) components.
               Bx(i,j)=bt(i,j)*cos(ba(i,j)+phase)
               By(i,j)=bt(i,j)*sin(ba(i,j)+phase)
               Bz(i,j)=bl(i,j)
   11       continue
   12    continue
      endif

c --> The FITS file must always be closed before exiting the program. 
c --> Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      call ftclos(unit, status)
      call ftfiou(unit, status)

c --> Check for any error, and if so print out error messages.
      if(status.gt.0) call printerror(status)

      return
      end
