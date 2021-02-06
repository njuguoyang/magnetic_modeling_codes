PRO HELIOGRAPHIC_PLANE,Bzi,Bhi,phii,Bzh,Bhh,phih,idim1,idim2,$
B0,L0,P,Bc,Lc,limb,iapp_tot

; PURPOSE: Define the heliographic plane and interpolate the
; heliographic components (Bzi,Bhi,phii) of the image plane into the
; heliographic components (Bzh,Bhh,phih) of the heliographic plane, where 
; B0, L0 --> Heliographic coordinates of the solar disk center
; P --> The solar P-angle
; Bc,Lc --> Heliographic coordinates of the image plane locations
; idim1,idim2 --> Linear size of the image plane
; limb --> An indicator of whether the target is close (limb=1) or far
; (limb=0) from the solar limb
; iapp_tot --> The mask of the nonzero field locations
; PROGRAMMER: Manolis K. Georgoulis (JHU/APL, 10/13/05)

COMMON PAR1,id1,id2,xi,yi,ksi,eta,mask
COMMON PAR2,iapp_toth

limb=0.
c11=cos(P)*cos(Lc-L0)- sin(P)*sin(B0)*sin(Lc-L0)
c12=-cos(P)*sin(Bc)*sin(Lc-L0) - sin(P)*(cos(B0)*cos(Bc) + sin(B0)*sin(Bc)*$
                                         cos(Lc-L0))
c21=sin(P)*cos(Lc-L0) + cos(P)*sin(B0)*sin(Lc-L0)
c22=-sin(P)*sin(Bc)*sin(Lc-L0) + cos(P)*(cos(B0)*cos(Bc) + sin(B0)*sin(Bc)*$
                                         cos(Lc-L0))

denom=1./(c11*c22 - c12*c21)
r1=where(denom gt 3.5 or denom lt 0.,ico) & if ico gt 0. then begin & $ 
denom(r1)=0. & limb=1. & endif
alpha=c22*denom
beta=-c21*denom
gamm=-c12*denom
delta=c11*denom
xi=dblarr(idim1,idim2)
yi=dblarr(idim1,idim2)
if n_elements(alpha) gt 1. then begin 
;  for i=0,idim1-1 do for j=0,idim2-1 do xi(i,j)=alpha(i,j)*i+gamm(i,j)*j    ;Original version by Manolis K. Georgoulis (JHU/APL, 10/13/05)
;  for i=0,idim1-1 do for j=0,idim2-1 do yi(i,j)=beta(i,j)*i+delta(i,j)*j
  for i=0,idim1-1 do for j=0,idim2-1 do xi(i,j)=alpha(i,j)*(i-0.5*(idim1-1)) + gamm(i,j)*(j-0.5*(idim2-1))    ;Modified by Yang Guo (NJU, 2011 Sep. 30)
  for i=0,idim1-1 do for j=0,idim2-1 do yi(i,j)=beta(i,j)*(i-0.5*(idim1-1)) + delta(i,j)*(j-0.5*(idim2-1))
endif else begin 
  for i=0,idim1-1 do for j=0,idim2-1 do xi(i,j)=alpha*i+gamm*j
  for i=0,idim1-1 do for j=0,idim2-1 do yi(i,j)=beta*i+delta*j
endelse
ksi=c11*xi + c12*yi
eta=c21*xi + c22*yi
ksi=ksi + 0.5*(idim1-1)    ;Modified by Yang Guo (NJU, 2011 Sep. 30) to be consistent with line 39
eta=eta + 0.5*(idim2-1)    ;Modified by Yang Guo (NJU, 2011 Sep. 30) to be consistent with line 40

;id1=fix(max(xi)-min(xi)+0.5)   ;Original version by Manolis K. Georgoulis (JHU/APL, 10/13/05)
;id2=fix(max(yi)-min(yi)+0.5)
id1=fix(max(xi)-min(xi)+1.0)   ;Modified by Yang Guo (NJU, 2011 Jun. 5)
id2=fix(max(yi)-min(yi)+1.0)
xi=xi-min(xi)
yi=yi-min(yi)

INTERPOLATE_NPFC,iapp_tot,iapp_toth,idim1,idim2,0,/fit_mask
r1=where(iapp_toth lt 0.98,ico) & if ico gt 0. then iapp_toth(r1)=0.
r1=where(iapp_toth ne 0.,ico) & if ico gt 0. then iapp_toth(r1)=1.
INTERPOLATE_NPFC,phii,phih,idim1,idim2,1
;INTERPOLATE_NPFC,sin(phii),sin_phih,idim1,idim2,1
;INTERPOLATE_NPFC,cos(phii),cos_phih,idim1,idim2,1
;RESTRUCTURE_ANGLE,sin_phih,cos_phih,phih
INTERPOLATE_NPFC,Bhi,Bhh,idim1,idim2,1
INTERPOLATE_NPFC,Bzi,Bzh,idim1,idim2,1

RETURN
END
