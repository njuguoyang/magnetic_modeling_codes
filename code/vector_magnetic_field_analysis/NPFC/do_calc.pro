PRO DO_CALC,Bx,By,Bz,dB_zp,Fz,Jhp,lamda

; PURPOSE:
; Calculation of the minimum-structure vertical magnetic field
; gradient (dB_zp), the vertical Lorentz force (Fz), and the minimum
; cross-field electric current density (Jhp)
; NOTES:
; (1) The calculation is done with the heliographic field components
;     on the image plane
; (2) The Gaussian unit system is used
; (3) lamda is the pixel size in cm
; PROGRAMMER: Manolis K. Georgoulis (JHU/APL, 10/13/05)

c=3.d10
B=sqrt(Bx^2.+By^2.+Bz^2.)
Bnx=Bx/B & Bny=By/B & Bnz=Bz/B 
r1=where(B eq 0.,ico) & if ico gt 0. then begin &$
Bnx(r1)=0. & Bny(r1)=0. & Bnz(r1)=0. &endif
dBz_x=(1./(2*lamda))*(shift(Bz,-1,0)-shift(Bz,1,0))
dBz_y=(1./(2*lamda))*(shift(Bz,0,-1)-shift(Bz,0,1))
dBx_x=(1./(2*lamda))*(shift(Bx,-1,0)-shift(Bx,1,0))
dBy_y=(1./(2*lamda))*(shift(By,0,-1)-shift(By,0,1))
dBz_z=-(dBx_x + dBy_y)

dB_x=(1./(2*lamda))*(shift(B,-1,0)-shift(B,1,0))
dB_y=(1./(2*lamda))*(shift(B,0,-1)-shift(B,0,1))
dB_zp=(Bnz/(Bnx^2.+Bny^2.))*(Bnx*dB_x + Bny*dB_y)
r1=where(dB_zp*0. ne 0.,ico) & if ico gt 0. then dB_zp(r1)=0.
RETRIEVE,dB_zp,Bnx^2.+Bny^2.,0.05

dBnx_x=(1./(2*lamda))*(shift(Bnx,-1,0)-shift(Bnx,1,0))
dBny_y=(1./(2*lamda))*(shift(Bny,0,-1)-shift(Bny,0,1))
dBnz_z=(-1./B)*((1./(Bnx^2.+Bny^2.))*(Bnx*dB_x + Bny*dB_y) + $
                B*(dBnx_x + dBny_y))

Fz=(1./(4.*!dpi))*(Bx*dBz_x + By*dBz_y + Bz*dBz_z - B*dB_zp)
Bh=sqrt(Bx^2.+By^2.)
Jhp=c*abs(Fz)/Bh
RETRIEVE,Jhp,Bh,10.

RETURN
END
