PRO HELIOGRAPHIC_IMAGE_PLANE,Bz,Btr,phi,Bz1,Bz2,Bh1,Bh2,phi1,phi2,$
idim1,idim2,B0,L0,P,B,L

; PURPOSE: From the LOS components (Bz,Btr,phi) calculate the two
; possible heliographic solutions (Bz1,Bh1,phi1) and (Bz2,Bh2,phi2) on
; the image plane, where: 
; B0, L0 --> Heliographic coordinates of the solar disk center
; P --> The solar P-angle
; B,L --> Heliographic coordinates of the image plane locations
; idim1,idim2 --> Linear size of the image plane
; PROGRAMMER: Manolis K. Georgoulis (JHU/APL, 10/13/05)

Tx1=-sin(B0)*sin(P)*sin(L-L0) + cos(P)*cos(L-L0)
Tx2=sin(B0)*cos(P)*sin(L-L0) + sin(P)*cos(L-L0)
Tx3=-cos(B0)*sin(L-L0)
Ty1=-sin(B)*(sin(B0)*sin(P)*cos(L-L0) + cos(P)*sin(L-L0)) $
             -cos(B)*cos(B0)*sin(P)
Ty2=sin(B)*(sin(B0)*cos(P)*cos(L-L0) - sin(P)*sin(L-L0)) $
                  +cos(B)*cos(B0)*cos(P)
Ty3=-cos(B0)*sin(B)*cos(L-L0) + sin(B0)*cos(B)
Tz1=cos(B)*(sin(B0)*sin(P)*cos(L-L0) + cos(P)*sin(L-L0)) $ 
            -sin(B)*cos(B0)*sin(P)
Tz2=-cos(B)*(sin(B0)*cos(P)*cos(L-L0) - sin(P)*sin(L-L0)) $ 
                   +sin(B)*cos(B0)*cos(P)
Tz3=cos(B)*cos(B0)*cos(L-L0) + sin(B)*sin(B0)
;
Bt=sqrt(Btr^2.+Bz^2.)
iflag=0.
recalc:
if iflag eq 0. then phi_los=phi else phi_los=phi+!dpi
Bx=Btr*cos(phi_los)
By=Btr*sin(phi_los)
Bxn=Bx/Bt
Byn=By/Bt
Bzn=Bz/Bt
r1=where(Bt eq 0.,ico) & if ico gt 0. then begin & Bxn(r1)=0. & Byn(r1)=0. &$
Bzn(r1)=0. &endif

Bx_loc=Bt*(Bxn*Tx1+ Byn*Tx2+ Bzn*Tx3)
By_loc=Bt*(Bxn*Ty1+ Byn*Ty2+ Bzn*Ty3)
Bz_loc=Bt*(Bxn*Tz1+ Byn*Tz2+ Bzn*Tz3)
Bh=sqrt(Bx_loc^2.+By_loc^2.)
cos_phi=Bx_loc/Bh
sin_phi=By_loc/Bh
RESTRUCTURE_ANGLE,sin_phi,cos_phi,phi_loc
r1=where(phi_loc*0. ne 0.,ico) & if ico gt 0. then phi_loc(r1)=0. 
if iflag eq 0. then begin 
  Bz1=Bz_loc
  Bh1=Bh
  phi1=phi_loc
  iflag=1
  goto, recalc
endif else begin 
  Bz2=Bz_loc
  Bh2=Bh
  phi2=phi_loc
endelse 

RETURN
END

