;+
; NAME :
;   differential3
; PURPOSE :
;   return rhs of ordinary differential equations
; CATEGORY :
;
; CALLING SEQUENCE :
;   Called by fieldline3d.pro, which is used for drawing 3 dimensional field lines
; INPUTS :
;
; OUTPUTS :
;
; COMMON BLOCKS :  common prm1,dx1,dy1,dz1,sig
;                  common bvertex,xv,DELX1,DELY1,DELZ1
; MODIFICATION HISTORY :
;   2008.10 Guo Yang Refer to SSW/packages/nlfff/idl/differential2.pro written by
;   jmm, jimm@ssl.berkeley.edu. The interpolation scheme is referred to FORTRAN code written by
;   Song M.T.
;-

function field_interp,bx,by,bz,xindex,yindex,zindex

msp=machar()   ;Machine-Specific Parameters for floating-point airthmetic. MACHAR(/DOUBLE) for double type 
eps=msp.eps

floor_x=FIX(xindex)
floor_y=FIX(yindex)
floor_z=FIX(zindex)

dx1=xindex-floor_x
dy1=yindex-floor_y
dz1=zindex-floor_z

xv=dblarr(3,8)

corner,floor_x,floor_y,floor_z,bx,by,bz,xv

Bxinterp=XITP(0,dx1,dy1,dz1,xv)
Byinterp=XITP(1,dx1,dy1,dz1,xv)
Bzinterp=XITP(2,dx1,dy1,dz1,xv)
Binterp=SQRT(Bxinterp*Bxinterp+Byinterp*Byinterp+Bzinterp*Bzinterp)

if (Binterp le eps) then begin
  result=[0.0,0.0,0.0]
;  print,'Null Points here !!!'
endif else begin
  result=[bxinterp,byinterp,bzinterp]
endelse
return,result

end
