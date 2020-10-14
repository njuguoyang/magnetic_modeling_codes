PRO RESTRUCTURE_ANGLE,sinx,cosx,x

; PURPOSE: Find an angle x in radians from its sine and cosine sinx
; and cosx, respectively
; PROGRAMMER: Manolis K. Georgoulis

x=acos(cosx)
r1=where(x*0. ne 0.,ico) & if ico gt 0. then x(r1)=0.
r1=where(sinx lt 0.,ico) & if ico gt 0. then x(r1)=2.*!dpi - x(r1)

RETURN
END
