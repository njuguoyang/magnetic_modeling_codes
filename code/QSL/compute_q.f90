include 'qsl3d.f90'

PROGRAM compute_q

!USE ieee_arithmetic
IMPLICIT NONE

INTEGER :: i,j,k,i1,nf
INTEGER :: nx0,ny0,nz0
INTEGER :: x1,x2,y1,y2,z1,z2
INTEGER :: nlevel,nx,ny,nz
REAL,DIMENSION(:,:,:),ALLOCATABLE :: qvalue,logQ
REAL    :: h0,delta0
REAL    :: res
REAL    :: time_begin,time_end
CHARACTER*1024 :: fmt,fn_b,fn_q,fn_q_vtk,fn_b_vtk,fn_sub,fn_bxyz
LOGICAL tf,aexist

nf = 9
res = 1.0
fmt = '(I4.4)'
tf = .TRUE.
200    FORMAT(A,I5,A,I5)

CALL CPU_TIME(time_begin)

do i1=nf,nf 
  write(*,'(A)') '========================================='
  write(*,'(A)') 'COMPUTE_Q:'
  WRITE(*,200) 'Total snapshots:',nf,'. Now snapshot:',i1
  write(*,'(A)') '========================================='
  write(fn_sub, fmt) i1
  fn_b = 'tiltresmhdthreedCforq'//trim(fn_sub)//'.blk'     ! Input file for output_b
  fn_q =                  'qmap'//trim(fn_sub)//'.dat'     ! Output file for qsl3d
  fn_q_vtk   =            'qmap'//trim(fn_sub)//'.vtk'     ! VTK file name for the Q value
  fn_b_vtk   =           'field'//trim(fn_sub)//'.vtk'     ! Output file for output_b 
  fn_bxyz =                                 'bxyz.dat'     ! Output file for output_b
  
  inquire(file='par',exist=aexist)
  call output_b(fn_b,fn_q,fn_b_vtk,fn_bxyz,nx0,ny0,nz0)

  ! To save the compuation time, only a subdomain is selected
  x1 = 1 + 149                !start point in x (it should be greater than or equal to 1)
  y1 = 1 + 20                 !start point in y (it should be greater than or equal to 1)
  z1 = 1                      !start point in z (it should be greater than or equal to 1)
  x2 = nx0-2 -148             !end point in x (it should be less than or equal to the x dimension-2)
  y2 = ny0-2 -20              !end point in y (it should be less than or equal to the y dimension-2)
  z2 = nz0-2                  !end point in z (it should be less than or equal to the z dimension-2)
  nlevel = 1                  !grid level. The Q value will be computed nlevel times finer than the original spatial resolution.
  h0 = 1.0/16.0               !integral step in computing magnetic field lines
  delta0 = 9.53674e-05        !scale for the interp

  OPEN(UNIT=3,FILE='par',ACTION='WRITE',STATUS='REPLACE',POSITION='REWIND',FORM='FORMATTED')
  write(3,'(a20)') fn_bxyz            !file name for the x,y,z components of the magnetic field
  write(3,'(a20)') fn_q               !file name for the output of qsl3d
  write(3,'(I6.0)') nx0               !dimension for x
  write(3,'(I6.0)') ny0               !dimension for y
  write(3,'(I6.0)') nz0               !dimension for z
  write(3,'(I6.0)') x1    
  write(3,'(I6.0)') y1
  write(3,'(I6.0)') z1
  write(3,'(I6.0)') x2
  write(3,'(I6.0)') y2
  write(3,'(I6.0)') z2
  write(3,'(I6.0)') nlevel
  write(3,'(e15.8)') h0
  write(3,'(e15.8)') delta0
  close(3)

  call qsl3d

  nx=nlevel*(x2-x1)+1
  ny=nlevel*(y2-y1)+1
  nz=nlevel*(z2-z1)+1
  ALLOCATE(qvalue(0:(nx-1),0:(ny-1),0:(nz-1)))
  ALLOCATE(logQ(0:(nx-1),0:(ny-1),0:(nz-1)))

  OPEN(UNIT=3,FILE=fn_q,STATUS='OLD',ACTION='READ',POSITION='REWIND',FORM='UNFORMATTED')
  DO k=0,nz-1
    DO j=0,ny-1
      DO i=0,nx-1
        read(3) qvalue(i,j,k)
      END DO
    END DO
  END DO
  CLOSE(3)
  
  !WHERE (ieee_is_nan(qvalue) == tf)
  !  qvalue = 2.0
  !ELSEWHERE
  !
  !END WHERE
  where (isNaN(qvalue) .eqv. tf)
     qvalue=2.0
  elsewhere
    !do nothing
  end where
  WHERE (qvalue < 2.0)
    qvalue = 2.0
  ELSEWHERE
    !do nothing
  END WHERE
  logQ = alog10(qvalue)

  OPEN(UNIT=3,FILE=fn_q_vtk,ACTION='WRITE',STATUS='REPLACE',POSITION='REWIND',FORM='FORMATTED')
  write(3,'(a)') '# vtk DataFile Version 2.0'
  write(3,'(a)') 'Volume example'
  write(3,'(a)') 'ASCII'
  write(3,'(a)') 'DATASET STRUCTURED_POINTS'
  write(3,'(a,3I6.0)') 'DIMENSIONS',nx,ny,nz
  write(3,'(a,3F20.8)') 'ASPECT_RATIO',res/nlevel,res/nlevel,res/nlevel
  write(3,'(a,3F20.8)') 'ORIGIN',x1*res,y1*res,z1*res
  write(3,'(a,I16.0)') 'POINT_DATA',nx*ny*nz
  write(3,'(a,I6.0)') 'SCALARS volume_scalars FLOAT',1
  write(3,'(a)') 'LOOKUP_TABLE default'

  102    FORMAT(1F20.8)
  DO k=0,nz-1
    DO j=0,ny-1
      DO i=0,nx-1
        write(3,102) logQ(i,j,k)
      END DO
    END DO
  END DO 
  close(3)

  DEALLOCATE(qvalue)
  DEALLOCATE(logQ)
end do

CALL CPU_TIME(time_end)
write(*,'(A)') '======================================================'
write(*,*)'calculation of all the snapshots takes time (hours):'
write(*,*)(time_end-time_begin)/3600.0
write(*,'(A)') '======================================================'

END PROGRAM compute_q

!==========================================================
!Purpose: loading data and preparing input files for QSL3D
!==========================================================
SUBROUTINE output_b(fn_b,fn_q,fn_b_vtk,fn_bxyz,nx,ny,nz)

IMPLICIT NONE

INTEGER :: ntotal,i,j,k
DOUBLE PRECISION :: t,res
REAL,DIMENSION(:,:,:),ALLOCATABLE :: x,y,z
REAL,DIMENSION(:,:,:),ALLOCATABLE :: rho
REAL,DIMENSION(:,:,:),ALLOCATABLE :: Bx,By,Bz
 CHARACTER*1024 :: header
 CHARACTER*1024, INTENT(IN) :: fn_b, fn_q, fn_b_vtk, fn_bxyz
INTEGER,INTENT(OUT) :: nx,ny,nz

res = 1.0

write(*,*)'OUTPUT_B: loading data and preparing input files for QSL3D'
OPEN(UNIT=3,STATUS='OLD',ACTION='READ',FILE=fn_b,POSITION='REWIND',FORM='UNFORMATTED')
read(3) header
read(3) ntotal,nx,ny,nz
read(3) t
ALLOCATE(x(0:(nx-1),0:(ny-1),0:(nz-1)))
ALLOCATE(y(0:(nx-1),0:(ny-1),0:(nz-1)))
ALLOCATE(z(0:(nx-1),0:(ny-1),0:(nz-1)))
!ALLOCATE(rho(0:(nx-1),0:(ny-1),0:(nz-1)))
ALLOCATE(Bx(0:(nx-1),0:(ny-1),0:(nz-1)))
ALLOCATE(By(0:(nx-1),0:(ny-1),0:(nz-1)))
ALLOCATE(Bz(0:(nx-1),0:(ny-1),0:(nz-1)))
DO k=0,nz-1
  DO j=0,ny-1
    DO i=0,nx-1
      READ(3) x(i,j,k),y(i,j,k),z(i,j,k),Bx(i,j,k),By(i,j,k),Bz(i,j,k)
    END DO
  END DO
END DO 
 CLOSE(3)

!DO k=0,0
!  DO j=0,0
!    DO i=0,2
!      write(*,*) x(i,j,k),y(i,j,k),z(i,j,k),rho(i,j,k),Bx(i,j,k),By(i,j,k),Bz(i,j,k)
!    END DO
!  END DO
!END DO 

OPEN(UNIT=3,FILE=fn_bxyz,ACTION='WRITE',STATUS='REPLACE',POSITION='REWIND',FORM='UNFORMATTED')
DO k=0,(nz-1)
  DO j=0,(ny-1)
    DO i=0,(nx-1)
      write(3) Bx(i,j,k),By(i,j,k),Bz(i,j,k)
    END DO
  END DO
END DO 
 CLOSE(3)

OPEN(UNIT=3,FILE=fn_b_vtk,ACTION='WRITE',STATUS='REPLACE',POSITION='REWIND',FORM='FORMATTED')
write(3,'(a)') '# vtk DataFile Version 2.0'
write(3,'(a)') 'Volume example'
write(3,'(a)') 'ASCII'
write(3,'(a)') 'DATASET STRUCTURED_POINTS'
write(3,'(a,3I6.0)') 'DIMENSIONS',nx,ny,nz
write(3,'(a,3F20.8)') 'ASPECT_RATIO',res,res,res
write(3,'(a,3F20.8)') 'ORIGIN',0.0,0.0,0.0
write(3,'(a,I16.0)') 'POINT_DATA',ntotal
write(3,'(a)') 'VECTORS volume_vectors FLOAT'

102    FORMAT(3F20.8)
DO k=0,nz-1
  DO j=0,ny-1
    DO i=0,nx-1
      write(3,102) Bx(i,j,k),By(i,j,k),Bz(i,j,k)
    END DO
  END DO
END DO 
 close(3)

DEALLOCATE(x)
DEALLOCATE(y)
DEALLOCATE(z)
!DEALLOCATE(rho) 
DEALLOCATE(Bx)
DEALLOCATE(By)
DEALLOCATE(Bz)
 
write(*,*)'OUTPUT_B: input files have been prepared! Magnetic field VTK file has been written!'
END SUBROUTINE output_b
