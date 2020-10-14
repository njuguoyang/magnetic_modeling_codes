;+
; NAME :
;   xitp
; PURPOSE :
;   Interpolation
; CATEGORY :
;
; CALLING SEQUENCE :
;   Called by fieldline3d.pro, which is used for drawing 3 dimensional field lines
; INPUTS :
;
; OUTPUTS :
;
; COMMON BLOCKS : 
; MODIFICATION HISTORY :
;   originally coded by M.T.Song by FORTRAN
;   2007.02 Guo Yang translated these codes from FORTRAN to IDL
;-

function xitp,nv,dx1,dy1,dz1,xv
;for interpolation
dx2=1.0-dx1
dy2=1.0-dy1
dz2=1.0-dz1
XITP=xv[nv,0]*dx2*dy2*dz2+xv[nv,1]*dx1*dy2*dz2+xv[nv,2]*dx2*dy1*dz2+xv[nv,3]*dx1*dy1*dz2+xv[nv,4]*dx2*dy2*dz1+xv[nv,5]*dx1*dy2*dz1+xv[nv,6]*dx2*dy1*dz1+xv[nv,7]*dx1*dy1*dz1
return,XITP
end
