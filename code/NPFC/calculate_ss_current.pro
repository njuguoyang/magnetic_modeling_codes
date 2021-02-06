PRO CALCULATE_SS_CURRENT,Btr,phi,lamda,Jz

; PURPOSE: Calculate the absolute value of the minimum vertical
; electric current density Jz as specified by Semel & Skumanich
; (1998), where Btr is the transverse field and phi is its azimuth
; PROGRAMMER: Manolis K. Georgoulis (JHU/APL) 10/06/05

Bx=Btr*cos(phi)
By=Btr*sin(phi)

dBxBy_x=(1./(2.*lamda))*(shift(Bx*By,-1,0)-shift(Bx*By,1,0))
dBxBy_y=(1./(2.*lamda))*(shift(Bx*By,0,-1)-shift(Bx*By,0,1))
dBx_By2_x=(1./(2.*lamda))*(shift((Bx^2.-By^2.),-1,0)-shift((Bx^2.-By^2.),1,0))
dBy_Bx2_y=(1./(2.*lamda))*(shift((By^2.-Bx^2.),0,-1)-shift((By^2.-Bx^2.),0,1))
dBx2_y=(1./(2.*lamda))*(shift(Bx^2.,0,-1)-shift(Bx^2.,0,1))
dBy2_x=(1./(2.*lamda))*(shift(By^2.,-1,0)-shift(By^2.,1,0))

gy=(1./Btr^4.)*(Bx^2. * dBxBy_x - 0.5*Bx*By*dBx_By2_x - 0.5*Btr^2.*dBx2_y)
gx=(1./Btr^4.)*(By^2. * dBxBy_y - 0.5*Bx*By*dBy_Bx2_y - 0.5*Btr^2.*dBy2_x)

Jz=sqrt(Bx^2.*gy^2. + By^2.*gx^2. -2.*Bx*By*gx*gy)
r1=where(Jz*0. ne 0. or Btr lt 1.,ico) & if ico gt 0. then Jz(r1)=0.

RETURN
END
