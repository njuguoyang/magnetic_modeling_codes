;+
; NAME:
;   null_newton_raphson
;
; PURPOSE:calculate the null position in the 3D computational domain 
;     	  by the Newton-Raphson scheme
;
; CALLING SEQUENCE:
;	Result = null_newton_raphson(cube_index,bx,by,bz,eps,matr,note)
;
; INPUTS:
;	cube_index: the box where the null locate with the poincare principle
;	bx,by,bz  : magnetic field
;	eps       : threshold
;
; KEYWORD PARAMETERS:
;         matr: the matrix of the null-point
;         note: the label of the null point position
;
; OUTPUTS:
;	return the position of the null-point.
;         
; MODIFICATION HISTORY:
; written by K. Yang 2014 Dec 3 NJU
; Modified by Z. Zhong (NJU, 2017 June 1) 
;-

FUNCTION null_newton_raphson,cube_index,bx,by,bz,eps,matr,note

xx=cube_index[0]
yy=cube_index[1]
zz=cube_index[2]
xindex = double(cube_index[0]);+0.5
yindex = double(cube_index[1]);+0.5
zindex = double(cube_index[2]);+0.5; give the initial position for the calculation in the computational box.
fieldb = DBLARR(3)
fieldb = field_interp(bx,by,bz,xindex,yindex,zindex)
matr = DBLARR(3,3)
step = 0
note = 0
ss=size(bz)

;WHILE ((abs(fieldb[0]) GT eps OR abs(fieldb[1]) GT eps OR abs(fieldb[2]) GT eps) AND step LT 10000 AND xindex GT xx AND xindex LT xx+1 AND yindex GT yy AND yindex LT yy+1 AND zindex GT zz ANd zindex LT zz+1) DO BEGIN

WHILE ((abs(fieldb[0]) GT eps OR abs(fieldb[1]) GT eps OR abs(fieldb[2]) GT eps) AND step LT 10000 AND xindex GE 1 AND xindex LT ss[1]-2 AND yindex GE 1 AND yindex LT ss[2]-2 AND zindex GE 1 ANd zindex LT ss[3]-2) DO BEGIN

fieldb = field_interp(bx,by,bz,xindex,yindex,zindex)
if (abs(fieldb[0]) le (machar()).eps and abs(fieldb[1]) le (machar()).eps and abs(fieldb[2]) le (machar()).eps) then break
matr = matrix_interp(bx,by,bz,xindex,yindex,zindex)
inverse_matr = LA_INVERT(matr)
xindex = xindex-(inverse_matr[0,0]*fieldb[0]+inverse_matr[0,1]*fieldb[1]+inverse_matr[0,2]*fieldb[2])
yindex = yindex-(inverse_matr[1,0]*fieldb[0]+inverse_matr[1,1]*fieldb[1]+inverse_matr[1,2]*fieldb[2])
zindex = zindex-(inverse_matr[2,0]*fieldb[0]+inverse_matr[2,1]*fieldb[1]+inverse_matr[2,2]*fieldb[2])
;fieldb = field_interp(bx,by,bz,xindex,yindex,zindex)

;print,xindex,yindex,zindex,'             ',fieldb

step = step+1
ENDWHILE

;print,'step',step

if (xindex LT 1 or xindex GT ss[1]-2 or yindex LT 1 or yindex GT ss[2]-2 or zindex LT 0 or zindex GT ss[3]-2) then begin
	note=0
endif else begin
	if (abs(fieldb[0]) gt eps or abs(fieldb[1]) gt eps or abs(fieldb[2]) gt eps) then begin
		note = 0
	endif else begin
		note=1
	endelse
endelse

;matr = matrix_interp(bx,by,bz,xindex,yindex,zindex)
posi = DBLARR(3)
posi[0] = xindex
posi[1] = yindex
posi[2] = zindex
return,posi

END
