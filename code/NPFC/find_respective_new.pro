PRO FIND_RESPECTIVE_NEW,Bz,Bz1,Bz2,phi,phi1,phi2

; PURPOSE: Choose a solution Bz between two possible solutions Bz1 and
; Bz2 if a solution phi between two possible solutions phi1 and phi2
; has been reached
; PROGRAMMER: Manolis K. Georgoulis (JHU/APL, 10/13/05)

res=size(phi) & idim1=res(1) & idim2=res(2)
Bz=dblarr(idim1,idim2)
r1=where(phi eq phi1,ic) & if ic gt 0. then Bz(r1)=Bz1(r1)
r1=where(phi eq phi2,ic) & if ic gt 0. then Bz(r1)=Bz2(r1)

RETURN
END
