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
; COMMON BLOCKS :  
;                  common bvertex,xv,Bx0,By0,Bz0
;                  common prm1,orig,DELX1,DELY1,DELZ1,sig
; MODIFICATION HISTORY :
;   2008.10 Guo Yang Refer to SSW/packages/nlfff/idl/differential2.pro written by
;   jmm, jimm@ssl.berkeley.edu. The interpolation scheme is referred to FORTRAN code written by
;   Song M.T.
;-

function differential3,s,coords

common bvertex,xv,Bx0,By0,Bz0,nx,ny,nz
common prm1,orig,DELX1,DELY1,DELZ1,sig
common prm2,Binterp

xd=coords[0]
yd=coords[1]
zd=coords[2]
xm=(xd-orig[0])/DELX1     ;convert to pixel position
ym=(yd-orig[1])/DELY1     ;convert to pixel position
zm=(zd-orig[2])/DELZ1     ;convert to pixel position
n1n=floor(xm)
n2n=floor(ym)
n3n=floor(zm)
if (n1n le -1 or n1n ge nx-1 or $
    n2n le -1 or n2n ge ny-1 or $
    n3n le -1 or n3n ge nz-1) then begin
  result=[0.0,0.0,0.0]
  return,result
endif else begin
  CORNER,n1n,n2n,n3n,Bx0,By0,Bz0,xv
  dx1=xm-n1n
  dy1=ym-n2n
  dz1=zm-n3n
  Bxinterp=XITP(0,dx1,dy1,dz1,xv)
  Byinterp=XITP(1,dx1,dy1,dz1,xv)
  Bzinterp=XITP(2,dx1,dy1,dz1,xv)
  Binterp=SQRT(Bxinterp*Bxinterp+Byinterp*Byinterp+Bzinterp*Bzinterp)

  msp=machar()   ;Machine-Specific Parameters for floating-point airthmetic. MACHAR(/DOUBLE) for double type 
  eps=msp.eps
  if (Binterp le eps) then begin
    result=[0.0,0.0,0.0]
    print,'Null Points Encountered!'
  endif else begin
    result=sig*[bxinterp/binterp,byinterp/binterp,bzinterp/binterp]
  endelse
  return,result
endelse
end
