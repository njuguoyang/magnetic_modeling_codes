!======================================================================================================
! Purpose :
!         share data for the program which calculate the squashing degree Q
!         Data is used by qsls3d.f90 and squashingdegreeq
!         virable should be set here : nx ny nz h delta nlevel0
!         nx ny nz is the dimension of the magntic data
!         h is integation step
!         delta is partial difference at the middle of the field line 
!         nlevel0 is resolution level 
!======================================================================================================
MODULE magnetic_data
       
IMPLICIT NONE

INTEGER :: nx0,ny0,nz0,nlevel0
INTEGER :: nx,ny,nz
INTEGER :: nx_start,ny_start,nz_start
INTEGER :: nx_end,ny_end,nz_end     
REAL :: h0,delta0
REAL :: h,delta 
REAL,DIMENSION(:,:,:),ALLOCATABLE :: Bx,By,Bz
REAL,DIMENSION(:,:,:),ALLOCATABLE :: qvalue

END MODULE magnetic_data

!======================================================================================================
! Program:
!          QSL3D
! Purpose: 
!          calculate squashing degree Q of a given magnetic field
! Functions needed:
!          squashingdegree,fieldline,integral,corner,differential3,xitp,rk8,friendposition,partialdifference,acrossposition,
! Data input:
!          Bx,By,Bz should be 3D data
! Parameters:
!          delta is the parameter during the calculation of the squashing degree Q, which means the difference distance of the origion field line at origion point
!          h is the parameter during the calculation of the field line, which is the step of the intergration
!          Nlevel is the level of refinement
! Algorithm:
!          The Q factor is proposed by Titov et al. 2002 (JGR)
!          This program is based on the method proposed by Pariat & Demoulin 2012 (A&A)
! History:
!          2009.11 GUO Yang (NJU), magnetic field line integration codes in IDL 
!          2014.09 YANG Kai (NJU), squashing degree computation in FORTRAN 
!======================================================================================================
SUBROUTINE qsl3d

USE magnetic_data
IMPLICIT NONE
      
INTEGER :: xn,yn,zn,i,j,k,totle
REAL :: sum0=0.0
REAL :: nlevel1,q
REAL :: flog1,error_flog,temp_q,temp_qvalue   
REAL :: xindex,yindex,zindex   
REAL :: time_begin,time_end
 CHARACTER*20 :: filename_b,outFileName

write(*,*)'QSL3D: reading parameters'

OPEN(UNIT=3,FILE='par',STATUS='OLD',ACTION='READ',POSITION='REWIND')
READ(3,'(a20)') filename_b
READ(3,'(a20)') outFileName
READ(3,'(I6.0)') nx0
READ(3,'(I6.0)') ny0
READ(3,'(I6.0)') nz0
READ(3,'(I6.0)') nx_start
READ(3,'(I6.0)') ny_start
READ(3,'(I6.0)') nz_start
READ(3,'(I6.0)') nx_end
READ(3,'(I6.0)') ny_end
READ(3,'(I6.0)') nz_end
READ(3,'(I6.0)') nlevel0
READ(3,'(e15.8)') h0
READ(3,'(e15.8)') delta0
 CLOSE(3)

delta=delta0
h=h0 
nx=nx0-2
ny=ny0-2
nz=nz0-2

ALLOCATE(Bx(0:(nx0-1),0:(ny0-1),0:(nz0-1)))
ALLOCATE(By(0:(nx0-1),0:(ny0-1),0:(nz0-1)))
ALLOCATE(Bz(0:(nx0-1),0:(ny0-1),0:(nz0-1)))
ALLOCATE(qvalue(0:(nlevel0*(nx_end-nx_start)),0:(nlevel0*(ny_end-ny_start)),0:(nlevel0*(nz_end-nz_start))))

 CALL CPU_TIME(time_begin)
nlevel1=REAL(nlevel0)
xn=nlevel0*(nx_end-nx_start)+1
yn=nlevel0*(ny_end-ny_start)+1
zn=nlevel0*(nz_end-nz_start)+1
totle=xn*yn*zn 

!Read magnetic field data
 
write(*,*)'QSL3D: loading magnetic data'
OPEN(UNIT=3,STATUS='OLD',ACTION='READ',FILE=filename_b,POSITION='REWIND',FORM='UNFORMATTED')
DO k=0,(nz0-1)
  DO j=0,(ny0-1)
    DO i=0,(nx0-1)
      READ(3) Bx(i,j,k),By(i,j,k),Bz(i,j,k)
    END DO
  END DO
END DO 
 CLOSE(3)

write(*,*)'QSL3D: calculation begin !'

301    FORMAT(A,F6.2,A)
200    FORMAT(A,I5,A,I5)
DO k=0,zn-1
  WRITE(*,200) 'Total step:',zn,'. Now step:',k
  DO j=0,yn-1
    DO i=0,xn-1    
      xindex = nx_start+i/nlevel1 
      yindex = ny_start+j/nlevel1           
      zindex = nz_start+k/nlevel1 
      CALL squashingdegree(xindex,yindex,zindex,h,q)
      qvalue(i,j,k)=q
    END DO
  END DO
  WRITE(*,301)'Percent complete:',100.0*(k+1.0)/(zn*1.0),'%'
END DO 

 CALL CPU_TIME(time_end)
write(*,*)'calculation takes time (hours):'
write(*,*)(time_end-time_begin)/3600.0
WRITE(*,*) 'Q computation complete!'

OPEN(UNIT=4,FILE=outFileName,ACTION='WRITE',STATUS='REPLACE',POSITION='REWIND',FORM='UNFORMATTED')
DO k=0,zn-1
  DO j=0,yn-1
    DO i=0,xn-1
      WRITE(4) qvalue(i,j,k)
    END DO
  END DO
END DO
 CLOSE(4)

DEALLOCATE(Bx)
DEALLOCATE(By)
DEALLOCATE(Bz)
DEALLOCATE(qvalue)    
!WRITE(*,*)'mission complete!'

END SUBROUTINE qsl3d

!=============================================================================================================
! Purpose:
!	calculate the factor Q. This subroutine is call by the main program qsls3d.f90
! Subroutines called here:
!	fieldline.f90 that calculate the position of the footpoints of a field line by RK8 integral
!	friendposition.f90 that gives the other 4 points' position which is needed during the partial difference calculation	
!	pritial.f90 that gives the partial elements	
! Syntax:
!	inputs
!		xindex, yindex, zindex is the position of the point for the Q value.
!		hh is the integral step usually 0.0625
!	outputs
!		squashing is the returned Q value
!=============================================================================================================
SUBROUTINE squashingdegree(xindex,yindex,zindex,hh,squashing)
     
USE magnetic_data
IMPLICIT NONE

REAL,INTENT(IN) :: xindex,yindex,zindex,hh
REAL :: squashing
REAL,DIMENSION(0:2) :: normal_vactor
REAL,DIMENSION(0:2) :: p1position0,p1position1,p1position2,p1position3,p1position4
REAL,DIMENSION(0:2) :: p2position0,p2position1,p2position2,p2position3,p2position4
REAL :: company(0:3,0:2)
REAL :: f0,factor,bz1,bz2,bnorm
REAL,DIMENSION(0:3) :: partiald1,partiald2
REAL :: dx1xc,dy1xc,dx1yc,dy1yc,dx2xc,dy2xc,dx2yc,dy2yc
REAL,DIMENSION(0:7,0:2) :: line0,line1,line2,line3,line4
INTEGER :: i,j,k

 CALL fieldline(xindex,yindex,zindex,line0,hh)
normal_vactor=line0(3:5,1)
 CALL friendposition(line0(0:2,1),normal_vactor,company)
bnorm=line0(6,1)
inner: IF(bnorm == 0.0 .or. line0(7,0) == 0 .or. line0(7,2) == 0) THEN
         squashing=1.0d0
       ELSE
         CALL fieldline(company(0,0),company(0,1),company(0,2),line1,hh)
         CALL fieldline(company(1,0),company(1,1),company(1,2),line2,hh)
         CALL fieldline(company(2,0),company(2,1),company(2,2),line3,hh)
         CALL fieldline(company(3,0),company(3,1),company(3,2),line4,hh)
         if (line0(7,0) == 1) bz1=line0(6,0)*line0(3,0)
         if (line0(7,0) == 2) bz1=line0(6,0)*line0(4,0)
         if (line0(7,0) == 3) bz1=line0(6,0)*line0(5,0)
         if (line0(7,2) == 1) bz2=line0(6,2)*line0(3,2)
         if (line0(7,2) == 2) bz2=line0(6,2)*line0(4,2)
         if (line0(7,2) == 3) bz2=line0(6,2)*line0(5,2)
         CALL partial(line0(0:5,0),line1(0:5,0),line2(0:5,0),&
         line3(0:5,0),line4(0:5,0),partiald1)
         CALL partial(line0(0:5,2),line1(0:5,2),line2(0:5,2),&
         line3(0:5,2),line4(0:5,2),partiald2)
         !write(*,*) 'partiald1,partiald2:',partiald1,partiald2
         dx1xc = partiald1(0)
         dy1xc = partiald1(1)
         dx1yc = partiald1(2)
         dy1yc = partiald1(3)
         dx2xc = partiald2(0)
         dy2xc = partiald2(1)
         dx2yc = partiald2(2)
         dy2yc = partiald2(3)
         f0=abs(  (bz1*bz2)/(bnorm**2*delta**4)  )
         factor=(dx2xc*dy1yc-dx2yc*dy1xc)**2+&
                (dx2yc*dx1xc-dx2xc*dx1yc)**2+&
                (dy2xc*dy1yc-dy2yc*dy1xc)**2+&
                (dy2yc*dx1xc-dy2xc*dx1yc)**2
         squashing = f0*factor
         !write(*,*) 'bnorm,bz1,bz2,delta:',bnorm,bz1,bz2,delta
         !write(*,*) 'factor,f0:',factor,f0
         !write(*,*) 'Q:',squashing
       END IF inner
END SUBROUTINE squashingdegree

!===============================================================================================================
! Purpose :
!	calculate the two footpoints of a field line.
!	The subroutine integral.f90 is called
! Syntax :
!	inputs
!		xindex,yindex,zindex they are the initial position of the integral
!		hh is the integral step
!	outputs
!		line, this varible returns the position coordinates and the the vector of direction of the field line.
!===============================================================================================================
SUBROUTINE fieldline(xindex,yindex,zindex,line,hh)

USE magnetic_data
IMPLICIT NONE

REAL,DIMENSION(0:7,0:2),INTENT(OUT) :: line 
REAL,INTENT(IN) :: xindex,yindex,zindex,hh
REAL :: sig
REAL,DIMENSION(0:7) :: line0,line1,line2
 
line0(0)=xindex
line0(1)=yindex
line0(2)=zindex
line0(3)=0.0
line0(4)=0.0
line0(5)=0.0
line0(6)=0.0
line0(7)=0.0  

sig = 1.0
 CALL integral(line0,line1,sig,hh)
sig = -1.0
 CALL integral(line0,line2,sig,hh)
line(0:7,0) = line1
line(0:7,1) = line0
line(0:7,2) = line2

END SUBROUTINE fieldline

!===============================================================================================================
! Purpose:
!      return the footpoints of the mangetic field lines
!===============================================================================================================
SUBROUTINE integral(posi1,posi2,sig,hh)

USE magnetic_data
IMPLICIT NONE

REAL,INTENT(IN) :: sig,hh
REAL,DIMENSION(0:7),INTENT(INOUT) :: posi1,posi2
REAL :: xyz(0:2),xv(0:2,0:7),dx1,dy1,dz1
REAL :: xm1,ym1,zm1,xm2,ym2,zm2,coords(0:2)
REAL :: tang0(0:3),tang1(0:3),tang2(0:3),dxdy(0:3),result0(0:2)
REAL :: hhh
INTEGER :: n1n,n2n,n3n
INTEGER :: ln,isn

ln  = 0
isn = 300000      
xm1 = posi1(0)
ym1 = posi1(1)
zm1 = posi1(2)
xm2 = xm1
ym2 = ym1
zm2 = zm1
hhh=hh

! Note that nx=nx0-2, and we adopt an indexing convention starting from 0. Therefore, subscript >=1.0 and <= nx 
! mean that one layer is truncated at the left and right sides. It is similar for other sides. The treatment is
! due to the fact that Q is not defined on the boundaries with the present formula. 
DO  WHILE(xm2 >= 1.0 .and. xm2 <= nx .and. ym2 >= 1.0 .and. ym2 <= ny .and. zm2 >= 1.0 .and. zm2 <= nz .and. ln <= isn )       
  xm1   = xm2
  ym1   = ym2
  zm1   = zm2
  tang1 = tang2
  n1n   = floor(xm1)
  n2n   = floor(ym1)
  n3n   = floor(zm1)
  CALL corner(n1n,n2n,n3n,xv)
  dx1    = xm1-n1n
  dy1    = ym1-n2n
  dz1    = zm1-n3n
  coords = (/xm1,ym1,zm1/)
  CALL differential3(coords,dxdy,dx1,dy1,dz1,sig,xv)
  tang2  = dxdy
  IF (ln == 0) THEN
    tang0 = tang2
  END IF
  IF(xm1 < 2 .or. ym1 < 2 .or. zm1 < 2 .or. xm1 > nx-1 .or. ym1 > ny-1 .or. zm1 > nz-1) THEN
    hhh  = 0.5*hh
  ELSE
    hhh  = hh
  END IF
  CALL rk8(coords,result0,sig,hhh)
  xm2 = result0(0)
  ym2 = result0(1)
  zm2 = result0(2)
  ln  = ln+1
END DO

posi1(3:6) = tang0
posi1(7)   = 1
      
 CALL acrossposition(xm1,ym1,zm1,xm2,ym2,zm2,posi2(0:2))
xm1=posi2(0)
ym1=posi2(1)
zm1=posi2(2)
if ((isNaN(xm1) .eqv. .TRUE.) .or. (isNaN(ym1) .eqv. .TRUE.) .or. (isNaN(zm1) .eqv. .TRUE.)) then
  n1n   = floor(xm2)
  n2n   = floor(ym2)
  n3n   = floor(zm2)
else
  n1n   = floor(xm1)
  n2n   = floor(ym1)
  n3n   = floor(zm1)
end if
if (xm1 < 0.0) n1n = 0
if (ym1 < 0.0) n2n = 0
if (zm1 < 0.0) n3n = 0
if (xm1 > nx) n1n = nx
if (ym1 > ny) n2n = ny
if (zm1 > nz) n3n = nz
 CALL corner(n1n,n2n,n3n,xv)
dx1    = xm1-n1n
dy1    = ym1-n2n
dz1    = zm1-n3n
coords = (/xm1,ym1,zm1/)
 CALL differential3(coords,dxdy,dx1,dy1,dz1,sig,xv)
posi2(3:6)= dxdy
posi2(7)=0
IF (xm2<1.0 .or. xm2>nx) THEN
  posi2(7) = 1
END IF
IF (ym2<1.0 .or. ym2>ny) THEN
  posi2(7) = 2
END IF
IF (zm2<1.0 .or. zm2>nz) THEN
  posi2(7) = 3
END IF
IF ((isNaN(xm1) .eqv. .TRUE.) .or. (isNaN(ym1) .eqv. .TRUE.) .or. (isNaN(zm1) .eqv. .TRUE.)) THEN
  posi2(7) = 4
END IF

END SUBROUTINE integral

!===============================================================================================================
! Purpose: 
!      return the mangnetic field on the 8 vertexes
!===============================================================================================================
SUBROUTINE corner (n1,n2,n3,xv)

USE magnetic_data
IMPLICIT NONE

INTEGER,INTENT(IN) :: n1,n2,n3
REAL,DIMENSION(0:2,0:7),INTENT(INOUT) :: xv
INTEGER :: n1p,n2p,n3p

n1p = n1+1
n2p = n2+1
n3p = n3+1

xv(0,0) = Bx(n1,n2,n3)
xv(0,1) = Bx(n1p,n2,n3)
xv(1,0) = By(n1,n2,n3)
xv(1,1) = By(n1p,n2,n3)
xv(2,0) = Bz(n1,n2,n3)
xv(2,1) = Bz(n1p,n2,n3)
xv(0,2) = Bx(n1,n2p,n3)
xv(0,3) = Bx(n1p,n2p,n3)
xv(1,2) = By(n1,n2p,n3)   
xv(1,3) = By(n1p,n2p,n3)
xv(2,2) = Bz(n1,n2p,n3)
xv(2,3) = Bz(n1p,n2p,n3)
xv(0,4) = Bx(n1,n2,n3p)
xv(0,5) = Bx(n1p,n2,n3p)
xv(1,4) = By(n1,n2,n3p)
xv(1,5) = By(n1p,n2,n3p)
xv(2,4) = Bz(n1,n2,n3p)
xv(2,5) = Bz(n1p,n2,n3p)
xv(0,6) = Bx(n1,n2p,n3p)
xv(0,7) = Bx(n1p,n2p,n3p)
xv(1,6) = By(n1,n2p,n3p)
xv(1,7) = By(n1p,n2p,n3p)
xv(2,6) = Bz(n1,n2p,n3p)
xv(2,7) = Bz(n1p,n2p,n3p)

END SUBROUTINE corner

!===============================================================================================================
! Purpose:
!      compute the right hand side of the equation dx/ds = Bx/B, dy/ds = By/B, dz/ds = Bz/B
!===============================================================================================================
SUBROUTINE differential3 (coords,dxdy,dx,dy,dz,sig,xv)
    
USE magnetic_data
IMPLICIT NONE

REAL,DIMENSION(0:2),INTENT(IN) :: coords
REAL,INTENT(INOUT),DIMENSION(0:3) :: dxdy
REAL,INTENT(IN) :: dx,dy,dz,sig
REAL,DIMENSION(0:2,0:7),INTENT(IN) :: xv
REAL :: eps,bxinter,byinter,bzinter,binter
eps = epsilon(1.0d0)
 CALL xitp(0,bxinter,xv,dx,dy,dz)
 CALL xitp(1,byinter,xv,dx,dy,dz)
 CALL xitp(2,bzinter,xv,dx,dy,dz)
 binter = SQRT(bxinter**2+byinter**2+bzinter**2)
IF (binter < eps) THEN
  dxdy = (/0.0d0,0.0d0,0.0d0,0.0d0/)
ELSE
  dxdy(0) = sig*bxinter/binter
  dxdy(1) = sig*byinter/binter
  dxdy(2) = sig*bzinter/binter
  dxdy(3) = binter
END IF

END SUBROUTINE differential3

!===============================================================================================================
! Purpose:
!      return the magnetic field at the given location (not necessarily on the vertexes)
! Inputs
!      magnetic field at the 8 vertexes and the coordinates refering to the corner with the 
!      lowest indices
!===============================================================================================================
SUBROUTINE xitp (nv,interp,xv,dx1,dy1,dz1)

USE magnetic_data
IMPLICIT NONE

INTEGER,INTENT(IN) :: nv
REAL,INTENT(OUT) :: interp
REAL,DIMENSION(0:2,0:7),INTENT(IN) :: xv
REAL,INTENT(IN)::dx1,dy1,dz1
REAL :: dx2,dy2,dz2

dx2 = 1.0d0 - dx1
dy2 = 1.0d0 - dy1
dz2 = 1.0d0 - dz1

interp=xv(nv,0)*dx2*dy2*dz2+xv(nv,1)*dx1*dy2*dz2+ &
       xv(nv,2)*dx2*dy1*dz2+xv(nv,3)*dx1*dy1*dz2+ &
       xv(nv,4)*dx2*dy2*dz1+xv(nv,5)*dx1*dy2*dz1+ &
       xv(nv,6)*dx2*dy1*dz1+xv(nv,7)*dx1*dy1*dz1

END SUBROUTINE xitp

!===============================================================================================================
! Purpose:
!       Runge-Kutta method to integrate the ordinary differential equation for the magnetic field line.
!===============================================================================================================
SUBROUTINE rk8 (x,result0,sig,h3)

USE magnetic_data
IMPLICIT NONE

INTEGER :: i
REAL,INTENT(IN) :: h3,x(0:2),sig
REAL,INTENT(OUT) :: result0(0:2)
REAL,DIMENSION(0:2) :: x0,x1
REAL,DIMENSION(0:3) :: y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12
REAL :: dx1,dy1,dz1
REAL,DIMENSION(0:2,0:7) :: xv
INTEGER :: m1m,m2m,m3m

DO i=0,2
  x0(i)=x(i)
END DO
DO i=0,2
  x1(i)=x(i)
END DO

! ------- 1 --------
m1m=FLOOR(x1(0))
m2m=FLOOR(x1(1))
m3m=FLOOR(x1(2))
if (x1(0) < 0.0) m1m = 0
if (x1(1) < 0.0) m2m = 0
if (x1(2) < 0.0) m3m = 0
if (x1(0) > nx)  m1m = nx
if (x1(1) > ny)  m2m = ny
if (x1(2) > nz)  m3m = nz
 CALL corner(m1m,m2m,m3m,xv)
dx1=x1(0)-m1m
dy1=x1(1)-m2m
dz1=x1(2)-m3m
 CALL differential3(x1,y0,dx1,dy1,dz1,sig,xv)
DO i=0,2
  x1(i)=x(i)+h3*2.0D0/27.0D0*y0(i)
END DO

! ------- 2 --------
m1m=FLOOR(x1(0))
m2m=FLOOR(x1(1))
m3m=FLOOR(x1(2))
if (x1(0) < 0.0) m1m = 0
if (x1(1) < 0.0) m2m = 0
if (x1(2) < 0.0) m3m = 0
if (x1(0) > nx)  m1m = nx
if (x1(1) > ny)  m2m = ny
if (x1(2) > nz)  m3m = nz
 CALL corner(m1m,m2m,m3m,xv)
dx1=x1(0)-m1m
dy1=x1(1)-m2m
dz1=x1(2)-m3m
 CALL differential3(x1,y1,dx1,dy1,dz1,sig,xv)
DO i=0,2
  x1(i)=x(i)+h3*(y0(i)+3.0D0*y1(i))/36.0D0
END DO

! ------- 3 --------
m1m=FLOOR(x1(0))
m2m=FLOOR(x1(1))
m3m=FLOOR(x1(2))
if (x1(0) < 0.0) m1m = 0
if (x1(1) < 0.0) m2m = 0
if (x1(2) < 0.0) m3m = 0
if (x1(0) > nx)  m1m = nx
if (x1(1) > ny)  m2m = ny
if (x1(2) > nz)  m3m = nz
 CALL corner(m1m,m2m,m3m,xv)
dx1=x1(0)-m1m
dy1=x1(1)-m2m
dz1=x1(2)-m3m
 CALL differential3(x1,y2,dx1,dy1,dz1,sig,xv)
DO i=0,2
  x1(i)=x(i)+h3*(y0(i)+3.0D0*y2(i))/24.0D0
END DO 

! ------- 4 --------
m1m=FLOOR(x1(0))
m2m=FLOOR(x1(1))
m3m=FLOOR(x1(2))
if (x1(0) < 0.0) m1m = 0
if (x1(1) < 0.0) m2m = 0
if (x1(2) < 0.0) m3m = 0
if (x1(0) > nx)  m1m = nx
if (x1(1) > ny)  m2m = ny
if (x1(2) > nz)  m3m = nz
 CALL corner(m1m,m2m,m3m,xv)
dx1=x1(0)-m1m
dy1=x1(1)-m2m
dz1=x1(2)-m3m
 CALL differential3(x1,y3,dx1,dy1,dz1,sig,xv)
DO i=0,2
  x1(i)=x(i)+h3*(y0(i)*20.0D0+(-y2(i)+y3(i))*75.0D0)/48.0D0
END DO

! ------- 5 --------
m1m=FLOOR(x1(0))
m2m=FLOOR(x1(1))
m3m=FLOOR(x1(2))
if (x1(0) < 0.0) m1m = 0
if (x1(1) < 0.0) m2m = 0
if (x1(2) < 0.0) m3m = 0
if (x1(0) > nx)  m1m = nx
if (x1(1) > ny)  m2m = ny
if (x1(2) > nz)  m3m = nz
 CALL corner(m1m,m2m,m3m,xv)
dx1=x1(0)-m1m
dy1=x1(1)-m2m
dz1=x1(2)-m3m
 CALL differential3(x1,y4,dx1,dy1,dz1,sig,xv)
DO i=0,2
  x1(i)=x(i)+h3*(y0(i)+y3(i)*5.0D0+y4(i)*4.0D0)/20.0D0
END DO

! ------- 6 --------
m1m=FLOOR(x1(0))
m2m=FLOOR(x1(1))
m3m=FLOOR(x1(2))
if (x1(0) < 0.0) m1m = 0
if (x1(1) < 0.0) m2m = 0
if (x1(2) < 0.0) m3m = 0
if (x1(0) > nx)  m1m = nx
if (x1(1) > ny)  m2m = ny
if (x1(2) > nz)  m3m = nz
 CALL corner(m1m,m2m,m3m,xv)
dx1=x1(0)-m1m
dy1=x1(1)-m2m
dz1=x1(2)-m3m
 CALL differential3(x1,y5,dx1,dy1,dz1,sig,xv)
DO i=0,2
  x1(i)=x(i)+h3*(-y0(i)*25.0D0+y3(i)*125.0D0-y4(i)*260.0D0+y5(i)*250.0D0)/108.0D0
END DO

! ------- 7 --------
m1m=FLOOR(x1(0))
m2m=FLOOR(x1(1))
m3m=FLOOR(x1(2))
if (x1(0) < 0.0) m1m = 0
if (x1(1) < 0.0) m2m = 0
if (x1(2) < 0.0) m3m = 0
if (x1(0) > nx)  m1m = nx
if (x1(1) > ny)  m2m = ny
if (x1(2) > nz)  m3m = nz
 CALL corner(m1m,m2m,m3m,xv)
dx1=x1(0)-m1m
dy1=x1(1)-m2m
dz1=x1(2)-m3m
 CALL differential3(x1,y6,dx1,dy1,dz1,sig,xv)
DO i=0,2
  x1(i)=x(i)+h3*(y0(i)*93.0D0+y4(i)*244.0D0-y5(i)*200.0D0+y6(i)*13.0D0)/900.0D0
END DO

! ------- 8 --------
m1m=FLOOR(x1(0))
m2m=FLOOR(x1(1))
m3m=FLOOR(x1(2))
if (x1(0) < 0.0) m1m = 0
if (x1(1) < 0.0) m2m = 0
if (x1(2) < 0.0) m3m = 0
if (x1(0) > nx)  m1m = nx
if (x1(1) > ny)  m2m = ny
if (x1(2) > nz)  m3m = nz
 CALL corner(m1m,m2m,m3m,xv)
dx1=x1(0)-m1m
dy1=x1(1)-m2m
dz1=x1(2)-m3m
 CALL differential3(x1,y7,dx1,dy1,dz1,sig,xv)
DO i=0,2
  x1(i)=x(i)+h3*(y0(i)*180.0D0-y3(i)*795.0D0+y4(i)*1408.0D0-y5(i)*1070.0D0+y6(i)*67.0D0+y7(i)*270.0D0)/90.0D0
END DO

! ------- 9 --------
m1m=FLOOR(x1(0))
m2m=FLOOR(x1(1))
m3m=FLOOR(x1(2))
if (x1(0) < 0.0) m1m = 0
if (x1(1) < 0.0) m2m = 0
if (x1(2) < 0.0) m3m = 0
if (x1(0) > nx)  m1m = nx
if (x1(1) > ny)  m2m = ny
if (x1(2) > nz)  m3m = nz
 CALL corner(m1m,m2m,m3m,xv)
dx1=x1(0)-m1m
dy1=x1(1)-m2m
dz1=x1(2)-m3m
 CALL differential3(x1,y8,dx1,dy1,dz1,sig,xv)
DO i=0,2
  x1(i)=x(i)+h3*(-y0(i)*455.0D0+y3(i)*115.0D0-y4(i)*3904.0D0+y5(i)*3110.0D0-y6(i)*171.0D0+y7(i)*1530.0D0-y8(i)*45.0D0)/540.0D0
END DO

! ------- 10 --------
m1m=FLOOR(x1(0))
m2m=FLOOR(x1(1))
m3m=FLOOR(x1(2))
if (x1(0) < 0.0) m1m = 0
if (x1(1) < 0.0) m2m = 0
if (x1(2) < 0.0) m3m = 0
if (x1(0) > nx)  m1m = nx
if (x1(1) > ny)  m2m = ny
if (x1(2) > nz)  m3m = nz
 CALL corner(m1m,m2m,m3m,xv)
dx1=x1(0)-m1m
dy1=x1(1)-m2m
dz1=x1(2)-m3m
 CALL differential3(x1,y9,dx1,dy1,dz1,sig,xv)
DO i=0,2
  x1(i)=x(i)+h3*(y0(i)*2383.0D0-y3(i)*8525.0D0+y4(i)*17984.0D0-y5(i)*15050.0D0+y6(i)*2133.0D0+y7(i)*2250.0D0+y8(i)*1125.0D0+y9(i)*1800.0D0)/4100.0D0
END DO

! ------- 11 --------
m1m=FLOOR(x1(0))
m2m=FLOOR(x1(1))
m3m=FLOOR(x1(2))
if (x1(0) < 0.0) m1m = 0
if (x1(1) < 0.0) m2m = 0
if (x1(2) < 0.0) m3m = 0
if (x1(0) > nx)  m1m = nx
if (x1(1) > ny)  m2m = ny
if (x1(2) > nz)  m3m = nz
 CALL corner(m1m,m2m,m3m,xv)
dx1=x1(0)-m1m
dy1=x1(1)-m2m
dz1=x1(2)-m3m
 CALL differential3(x1,y10,dx1,dy1,dz1,sig,xv)
DO i=0,2
  x1(i)=x(i)+h3*(y0(i)*60.0D0-y5(i)*600.0D0-y6(i)*60.0D0+(y8(i)-y7(i)+2.0D0*y9(i))*300.0D0)/4100.0D0
END DO

! ------- 12 --------
m1m=FLOOR(x1(0))
m2m=FLOOR(x1(1))
m3m=FLOOR(x1(2))
if (x1(0) < 0.0) m1m = 0
if (x1(1) < 0.0) m2m = 0
if (x1(2) < 0.0) m3m = 0
if (x1(0) > nx)  m1m = nx
if (x1(1) > ny)  m2m = ny
if (x1(2) > nz)  m3m = nz
 CALL corner(m1m,m2m,m3m,xv)
dx1=x1(0)-m1m
dy1=x1(1)-m2m
dz1=x1(2)-m3m
 CALL differential3(x1,y11,dx1,dy1,dz1,sig,xv)
DO i=0,2
  x1(i)=x(i)+h3*(-y0(i)*1777.0D0-y3(i)*8525.0D0+y4(i)*17984.0D0-y5(i)*14450.0D0+y6(i)*2193.0D0+y7(i)*2550.0D0+y8(i)*825.0D0+y9(i)*1200.0D0+y11(i)*4100.0D0)/4100.0D0
END DO

! ------- 13 --------
m1m=FLOOR(x1(0))
m2m=FLOOR(x1(1))
m3m=FLOOR(x1(2))
if (x1(0) < 0.0) m1m = 0
if (x1(1) < 0.0) m2m = 0
if (x1(2) < 0.0) m3m = 0
if (x1(0) > nx)  m1m = nx
if (x1(1) > ny)  m2m = ny
if (x1(2) > nz)  m3m = nz
 CALL corner(m1m,m2m,m3m,xv)
dx1=x1(0)-m1m
dy1=x1(1)-m2m
dz1=x1(2)-m3m
 CALL differential3(x1,y12,dx1,dy1,dz1,sig,xv)
DO i=0,2
  result0(i)=x(i)+h3*(y5(i)*272.0D0+(y6(i)+y7(i))*216.0D0+(y8(i)+y9(i))*27.0D0+(y11(i)+y12(i))*41.0D0)/840.0D0
END DO

END SUBROUTINE rk8

!===============================================================================================================
! Purpose:
!      fix the mangetic field line position in last step when it crosses the computation boundaries
!===============================================================================================================
SUBROUTINE acrossposition (x0,y0,z0,x1,y1,z1,acrossp)
           
USE magnetic_data
IMPLICIT NONE

REAL,INTENT(IN) :: x0,y0,z0,x1,y1,z1
REAL,INTENT(out) :: acrossp(0:2)
REAL,DIMENSION(0:2) :: position0,position1,position2,position3,position4,position5,position6,vactor
REAL :: k1,k2,k3,k4,k5,k6

position0(0)=x1
position0(1)=y1
position0(2)=z1
vactor(0)=x0-x1
vactor(1)=y0-y1
vactor(2)=z0-z1

! foot point on plane 1 x=1
k1=(1-x1)/vactor(0)
position1=position0+k1*vactor
position1(0)=1.0d0
! foot point on plane 2 x=nx, where nx=nx0-2. Since we use the index convention such that the first position starts from 0, nx indicates the position one layer away from the boundary.
k2=(nx-x1)/vactor(0)
position2=position0+k2*vactor
position2(0)=nx

! foot point on plane 3 y=1
k3=(1-y1)/vactor(1)
position3=position0+k3*vactor
position3(1)=1.0d0
! foot point on plane 4 y=ny
k4=(ny-y1)/vactor(1)
position4=position0+k4*vactor
position4(1)=ny

! foot point on plane 5 z=1
k5=(1-z1)/vactor(2)
position5=position0+k5*vactor
position5(2)=1.0d0
! foot point on plane 6 z=zn-1
k6=(nz-z1)/vactor(2)
position6=position0+k6*vactor
position6(2)=nz

IF (position1(0) >= 1. .AND. position1(0) <= nx .AND. position1(1) >= 1. .AND. position1(1) <= ny .AND. position1(2) >= 1. .AND. position1(2) <= nz .AND. k1 >= 0. .AND. k1 <= 1.) THEN
  acrossp=position1
END IF
IF (position2(0) >= 1. .AND. position2(0) <= nx .AND. position2(1) >= 1. .AND. position2(1) <= ny .AND. position2(2) >= 1. .AND. position2(2) <= nz .AND. k2 >= 0. .AND. k2 <= 1.) THEN
  acrossp=position2
END IF
IF (position3(0) >= 1. .AND. position3(0) <= nx .AND. position3(1) >= 1. .AND. position3(1) <= ny .AND. position3(2) >= 1. .AND. position3(2) <= nz .AND. k3 >= 0. .AND. k3 <= 1.) THEN
  acrossp=position3
END IF
IF (position4(0) >= 1. .AND. position4(0) <= nx .AND. position4(1) >= 1. .AND. position4(1) <= ny .AND. position4(2) >= 1. .AND. position4(2) <= nz .AND. k4 >= 0. .AND. k4 <= 1.) THEN
  acrossp=position4
END IF
IF (position5(0) >= 1. .AND. position5(0) <= nx .AND. position5(1) >= 1. .AND. position5(1) <= ny .AND. position5(2) >= 1. .AND. position5(2) <= nz .AND. k5 >= 0. .AND. k5 <= 1.) THEN
  acrossp=position5
END IF
IF (position6(0) >= 1. .AND. position6(0) <= nx .AND. position6(1) >= 1. .AND. position6(1) <= ny .AND. position6(2) >= 1. .AND. position6(2) <= nz .AND. k6 >= 0. .AND. k6 <= 1.) THEN
  acrossp=position6
END IF
END SUBROUTINE acrossposition

!==============================================================================================
! Purpose:
!       Calculate the position of 4 points company the origion point which is the midule of the field line
!       Used by file squashingdegree.f90
!==============================================================================================
SUBROUTINE friendposition (coordination,nvactor,company)

USE magnetic_data
IMPLICIT NONE

REAL,INTENT(IN),DIMENSION(0:2)::coordination,nvactor
REAL,DIMENSION(0:2)::vactor1,vactor2,position1,position2,position3,position4,norm0,xyz1
REAL,INTENT(OUT),DIMENSION(0:3,0:2)::company
REAL::n1,n2,n3,flag,const

n1=1.0d0
n2=1.0d0
n3=1.0d0
flag=1.0d0
const=1.0d0

xyz1=coordination
n1=nvactor(0)
n2=nvactor(1)
n3=nvactor(2)

vactor1(0)=-n2/(SQRT(1.0d0-n3**2))
vactor1(1)=n1/(SQRT(1.0d0-n3**2))
vactor1(2)=0.0d0

norm0(0)=-n3+(n2*n2*n3)/(1.0d0-n3*n3)
norm0(1)=-n1*n2*n3/(1.0d0-n3*n3)
norm0(2)=n1

const=SQRT(norm0(0)**2+norm0(1)**2+norm0(2)**2) 
vactor2(0)=norm0(0)/const
vactor2(1)=norm0(1)/const
vactor2(2)=norm0(2)/const

flag=n1*(vactor1(1)*vactor2(2)-vactor1(2)*vactor2(1))-n2*(vactor1(0)*vactor2(2)-vactor1(2)*vactor2(0))+n3*(vactor1(0)*vactor2(1)-vactor1(1)*vactor2(0))

IF (flag > 0) THEN 
  position1=xyz1+delta*vactor1
  position2=xyz1+delta*vactor2
  position3=xyz1-delta*vactor1
  position4=xyz1-delta*vactor2
ELSE
  position1=xyz1+delta*vactor2
  position2=xyz1+delta*vactor1
  position3=xyz1-delta*vactor2
  position4=xyz1-delta*vactor1
END IF
company(0,:)=position1
company(1,:)=position2
company(2,:)=position3
company(3,:)=position4

END SUBROUTINE friendposition  

!==============================================================================================
! Purpose:
!       Calculate the partial difference at the end of the field line
!       used by file squashingdegree.f90
!==============================================================================================
SUBROUTINE partial (line0,line1,line2,line3,line4,partiald)

USE magnetic_data
IMPLICIT NONE

REAL,INTENT(IN),DIMENSION(0:5)::line0,line1,line2,line3,line4
REAL,INTENT(OUT),DIMENSION(0:3)::partiald
REAL,DIMENSION(0:2)::position0,position1,position2,position3,position4,position5,position6,vactor1,vactor2,vactor3,vactor4,long15,long26,long46,long35,long06,long05,normal
REAL::sig1,sig2,sig3,sig4,flag0,flag31,flag32,flag41,flag42
REAL::dxnxc,dynxc,dxnyc,dynyc,fenzi1,fenzi2,fenmu1,fenmu2,k1,k2
 
sig1   = 1.0d0
sig2   = 1.0d0
sig3   = 1.0d0
sig4   = 1.0d0
flag0  = 1.0d0
flag31 = 1.0d0
flag32 = 1.0d0
flag41 = 1.0d0
flag42 = 1.0d0

long15 = (/0.0d0,0.0d0,0.0d0/)
long26 = (/0.0d0,0.0d0,0.0d0/)
long46 = (/0.0d0,0.0d0,0.0d0/)
long35 = (/0.0d0,0.0d0,0.0d0/)
long06 = (/0.0d0,0.0d0,0.0d0/) 
long05 = (/0.0d0,0.0d0,0.0d0/)

position0 = line0(0:2)
position1 = line1(0:2)
position2 = line2(0:2)
position3 = line3(0:2)
position4 = line4(0:2)
normal    = line0(3:5)

vactor1   = position0-position1
vactor2   = position0-position2
vactor3   = position3-position0
vactor4   = position4-position0
!write(*,*) 'vactor1:',vactor1

!fenzi1    = (position3(0)-position0(0))*vactor1(0)+(position3(1)-position0(1))*vactor1(1)+(position3(2)-position0(2))*vactor1(2)
fenzi1    = vactor3(0)*vactor1(0)+vactor3(1)*vactor1(1)+vactor3(2)*vactor1(2)
fenmu1    = vactor1(0)**2+vactor1(1)**2+vactor1(2)**2
k1        = fenzi1/fenmu1
!fenzi2    = (position4(0)-position0(0))*vactor2(0)+(position4(1)-position0(1))*vactor2(1)+(position4(2)-position0(2))*vactor2(2)
fenzi2    = vactor4(0)*vactor2(0)+vactor4(1)*vactor2(1)+vactor4(2)*vactor2(2)
fenmu2    = vactor2(0)**2+vactor2(1)**2+vactor2(2)**2
k2        = fenzi2/fenmu2

position5 = position0+k1*vactor1
position6 = position0+k2*vactor2
!write(*,*) 'k1,k2:',k1,k2

long15    = k1*vactor1
long26    = k2*vactor2  
long46    = position6-position4
long35    = position5-position3 
long05    = position0-position5
long06    = position0-position6

flag0     = normal(0)*(vactor1(1)*vactor2(2)-vactor1(2)*vactor2(1))-normal(1)*(vactor1(0)*vactor2(2)-vactor1(2)*vactor2(0))+normal(2)*(vactor1(0)*vactor2(1)-vactor1(1)*vactor2(0))
flag31    = normal(0)*(vactor3(1)*vactor1(2)-vactor3(2)*vactor1(1))-normal(1)*(vactor3(0)*vactor1(2)-vactor3(2)*vactor1(0))+normal(2)*(vactor3(0)*vactor1(1)-vactor3(1)*vactor1(0))          
flag32    = normal(0)*(vactor3(1)*vactor2(2)-vactor3(2)*vactor2(1))-normal(1)*(vactor3(0)*vactor2(2)-vactor3(2)*vactor2(0))+normal(2)*(vactor3(0)*vactor2(1)-vactor3(1)*vactor2(0))     
flag41    = normal(0)*(vactor4(1)*vactor1(2)-vactor4(2)*vactor1(1))-normal(1)*(vactor4(0)*vactor1(2)-vactor4(2)*vactor1(0))+normal(2)*(vactor4(0)*vactor1(1)-vactor4(1)*vactor1(0))          
flag42    = normal(0)*(vactor4(1)*vactor2(2)-vactor4(2)*vactor2(1))-normal(1)*(vactor4(0)*vactor2(2)-vactor4(2)*vactor2(0))+normal(2)*(vactor4(0)*vactor2(1)-vactor4(1)*vactor2(0))    

outer:IF(flag0 > 0) THEN
 
        IF(flag31 > 0 .AND. flag32 > 0) THEN
          sig2=1.0d0
          sig4=1.0d0
        END IF
  
        IF(flag31 < 0 .AND. flag32 > 0) THEN
          sig2=1.0d0
          sig4=-1.0d0
        END IF

        IF(flag31 < 0 .AND. flag32 < 0) THEN
          sig2=-1.0d0
          sig4=-1.0d0
        END IF

        IF(flag31 > 0 .AND. flag32 < 0) THEN
          sig2=-1.0d0
          sig4=1.0d0
        END IF       

        IF(flag41 > 0 .AND. flag42 > 0) THEN
          sig1=-1.0d0
          sig3=-1.0d0
        END IF

        IF(flag41 < 0 .AND. flag42 > 0) THEN
          sig1=1.0d0
          sig3=-1.0d0
        END IF

        IF(flag41 < 0 .AND. flag42 < 0) THEN
          sig1=1.0d0
          sig3=1.0d0
        END IF

        IF(flag41 > 0 .AND. flag42 < 0) THEN
          sig1=-1.0d0
          sig3=1.0d0
        END IF      

      ELSE

        IF(flag31 > 0 .AND. flag32 < 0) THEN
          sig2=1.0d0
          sig4=-1.0d0
        END IF

        IF(flag31 < 0 .AND. flag32 < 0) THEN
          sig2=1.0d0
          sig4=1.0d0
        END IF

        IF(flag31 < 0 .AND. flag32 > 0) THEN
          sig2=-1.0d0
          sig4=1.0d0
        END IF

        IF(flag31 > 0 .AND. flag32 > 0) THEN
          sig2=-1.0d0
          sig4=-1.0d0
        END IF

        IF(flag41 > 0 .AND. flag42 < 0) THEN
          sig1=-1.0d0
          sig3=-1.0d0
        END IF

        IF(flag41 < 0 .AND. flag42 < 0) THEN
          sig1=1.0d0
          sig3=-1.0d0
        END IF

        IF(flag41 < 0 .AND. flag42 > 0) THEN
          sig1=1.0d0
          sig3=1.0d0
        END IF

        IF(flag41 > 0 .AND. flag42 > 0) THEN
          sig1=-1.0d0
          sig3=1.0d0
        END IF  

      END IF outer  

Dynyc = sig1*SQRT(long06(0)**2+long06(1)**2+long06(2)**2)
Dxnxc = sig2*SQRT(long05(0)**2+long05(1)**2+long05(2)**2)
Dxnyc = sig3*SQRT(long46(0)**2+long46(1)**2+long46(2)**2)
Dynxc = sig4*SQRT(long35(0)**2+long35(1)**2+long35(2)**2)

partiald(0) = dxnxc
partiald(1) = dynxc
partiald(2) = dxnyc
partiald(3) = dynyc

END SUBROUTINE partial
