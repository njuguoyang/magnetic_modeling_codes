PRO CLOSER_ANGLE,phi,phi1,phi2,phif

; PURPOSE: Find the closest observed angle between phi1 and phi2 when a certain
; angle phi is available
; PROGRAMMER: Manolis K. Georgoulis (JHU/APL, 10/13/05)

res=size(phi) & id1=res(1) & id2=res(2)
sin_phi=sin(phi)
cos_phi=cos(phi)
d1=abs(cos_phi-cos(phi1))+abs(sin_phi-sin(phi1))
d2=abs(cos_phi-cos(phi2))+abs(sin_phi-sin(phi2))
phif=dblarr(id1,id2)
r1=where(d1 le d2,ico) & if ico gt 0. then phif(r1)=phi1(r1)
r1=where(d1 gt d2,ico) & if ico gt 0. then phif(r1)=phi2(r1)

RETURN
END
