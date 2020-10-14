PRO CALCULATE_VERTICAL_CURRENT,Bx,By,lamda,Jz

; PURPOSE:
; The vertical component Jz of the electric current density is calculated. 
; NOTES:
; (1) The Gaussian unit system is used
; (2) lamda is the pixel size in cm
; PROGRAMMER: Manolis K. Georgoulis (JHU/APL, 10/13/05)

t1=(1./(2.*lamda))*(shift(By,-1,0)-shift(By,1,0))
t2=(1./(2.*lamda))*(shift(Bx,0,-1)-shift(Bx,0,1))
c=3.d10 & Jz=(c/(4.*!dpi))*(t1-t2)


RETURN
END
