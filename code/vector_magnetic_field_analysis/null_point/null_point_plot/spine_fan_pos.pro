;+
; NAME :
;   spine_fan_pos
; PURPOSE :
;   determine the start points of magnetic field lines in the vicinity of a null point
; CATEGORY :
; CALLING SEQUENCE :
; INPUTS :
; OPTIONAL INPUT PARAMETERS (KEYWORD PARAMETERS):
; OUTPUTS :
; COMMON BLOCKS :
; SIDE EFFECTS :
; RESTRICTIONS :
; PROCEDURE :
; MODIFICATION HISTORY :
;   2010.05 Guo Yang@Nanjing University
;
;   This software is provided without any warranty. Permission to use,
;   copy, modify. Distributing modified or unmodified copies is granted,
;   provided this disclaimer and information are included unchanged.
;-

pro spine_fan_pos, r, theta, nphi, null_pos, evecx, evecy, evecz, xindex=xindex, yindex=yindex, zindex=zindex

phi=2.0d*!dpi*findgen(nphi)/double(nphi)
xl=r*sin(theta)*cos(phi)         ;Cartesian coordinates in the local spherical coordinate system
yl=r*sin(theta)*sin(phi)
zl=r*cos(theta)*(dblarr(nphi)+1.0d)
thetaz=acos(evecz[2]/sqrt(evecz[0]*evecz[0]+evecz[1]*evecz[1]+evecz[2]*evecz[2])) 
if (abs(evecz[0]) ge 0.000001 or abs(evecz[1]) ge 0.000001) then begin
  phiz=atan(evecz[1],evecz[0])
  if (abs(evecz[2]) ge 0.000001) then begin
    pvzl=[evecz[0],evecz[1],-(evecz[0]*evecz[0]+evecz[1]*evecz[1])/evecz[2]]  ;Projection vector of eigen vector z on xy-plane in the local spherical coordinate system
  endif else begin
    pvzl=[0.0,0.0,1.0]
  endelse
  phixz=acos((pvzl[0]*evecx[0]+pvzl[1]*evecx[1]+pvzl[2]*evecx[2])/ $
        sqrt(pvzl[0]*pvzl[0]+pvzl[1]*pvzl[1]+pvzl[2]*pvzl[2])/       $
        sqrt(evecx[0]*evecx[0]+evecx[1]*evecx[1]+evecx[2]*evecx[2]))        ;Angle between x_local and pvzl
  phiyz=acos((pvzl[0]*evecy[0]+pvzl[1]*evecy[1]+pvzl[2]*evecy[2])/ $
        sqrt(pvzl[0]*pvzl[0]+pvzl[1]*pvzl[1]+pvzl[2]*pvzl[2])/       $
        sqrt(evecy[0]*evecy[0]+evecy[1]*evecy[1]+evecy[2]*evecy[2]))        ;Angle between y_local and pvzl
  ;Rotate about z-axis
  if (phiyz le !dpi/2.0d) then begin
    x2= xl*cos(phixz)+yl*sin(phixz)
    y2=-xl*sin(phixz)+yl*cos(phixz)
    z2= zl
  endif else begin
    x2= xl*cos(phixz)-yl*sin(phixz)
    y2= xl*sin(phixz)+yl*cos(phixz)
    z2=zl
  endelse
  if (evecz[2] ge 0.000001) then begin
    ;Rotate about y-axis
    xl=x2*cos(thetaz)+z2*sin(thetaz)
    yl=y2
    zl=-x2*sin(thetaz)+z2*cos(thetaz)
    ;Rotate about z-axis
    x2= xl*cos(phiz)-yl*sin(phiz)
    y2= xl*sin(phiz)+yl*cos(phiz)
    z2= zl
  endif else begin
    ;Rotate about y-axis
    xl=x2*cos(thetaz)-z2*sin(thetaz)
    yl=y2
    zl=x2*sin(thetaz)+z2*cos(thetaz)
    ;Rotate about z-axis
    x2= xl*cos(phiz+!dpi)-yl*sin(phiz+!dpi)
    y2= xl*sin(phiz+!dpi)+yl*cos(phiz+!dpi)
    z2= zl
  endelse 
endif else begin
  phix=atan(evecx[1],evecx[0])
  if (evecz gt 0.0) then begin
    ;Rotate about z-axis
    x2= xl*cos(phix)-yl*sin(phix)
    y2= xl*sin(phix)+yl*cos(phix)
    z2= zl
  endif else begin
    ;Rotate about x-axis to flip z-axis
    x2= xl
    y2=-yl
    z2=-zl
    ;Rotate about z-axis
    xl= x2*cos(phix)-y2*sin(phix)
    yl= x2*sin(phix)+y2*cos(phix)
    zl= z2
    x2=xl
    y2=yl
    z2=zl
  endelse
endelse

xindex = null_pos[0]+x2
yindex = null_pos[1]+y2
zindex = null_pos[2]+z2

end
