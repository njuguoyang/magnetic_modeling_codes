;+
; NAME:
;   matrix_interp
;
; PURPOSE:calculate the differencial matrix at the given 
;		  point by a linear interpolation.
;
; CALLING SEQUENCE:
;	Result = matrix_interp(bx,by,bz,xindex,yindex,zindex)  	  

; INPUTS:
;	bx,by,bz  : magnetic field
;	xindex    : The X location of null point
;	yindex    : The Y location of null point
;	zindex    : The Z location of null point
;
; OUTPUTS:
;	return the the matrix of the at the given point 
;         
; MODIFICATION HISTORY:
; written by Z. Zhong (NJU, 2017 June 1) 
;-

function matrix_interp,bx,by,bz,xindex,yindex,zindex

Matrix_B=dblarr(3,3)
delta=double(1e-6)

x_000=xindex-delta
y_000=yindex
z_000=zindex
field_000 = field_interp(bx,by,bz,x_000,y_000,z_000)

x_001=xindex+delta
y_001=yindex
z_001=zindex
field_001 = field_interp(bx,by,bz,x_001,y_001,z_001)

x_010=xindex
y_010=yindex-delta
z_010=zindex
field_010 = field_interp(bx,by,bz,x_010,y_010,z_010)

x_011=xindex
y_011=yindex+delta
z_011=zindex
field_011 = field_interp(bx,by,bz,x_011,y_011,z_011)

x_100=xindex
y_100=yindex
z_100=zindex-delta
field_100 = field_interp(bx,by,bz,x_100,y_100,z_100)

x_101=xindex
y_101=yindex
z_101=zindex+delta
field_101 = field_interp(bx,by,bz,x_101,y_101,z_101)

if (zindex eq 0) then begin
  
Matrix_B[0,0]=(field_001[0]-field_000[0])/(delta*2.0d)
Matrix_B[0,1]=(field_011[0]-field_010[0])/(delta*2.0d)
Matrix_B[0,2]=(field_101[0]-field_100[0])/(delta)

Matrix_B[1,0]=(field_001[1]-field_000[1])/(delta*2.0d)
Matrix_B[1,1]=(field_011[1]-field_010[1])/(delta*2.0d)
Matrix_B[1,2]=(field_101[1]-field_100[1])/(delta)

Matrix_B[2,0]=(field_001[2]-field_000[2])/(delta*2.0d)
Matrix_B[2,1]=(field_011[2]-field_010[2])/(delta*2.0d)
Matrix_B[2,2]=(field_101[2]-field_100[2])/(delta)

endif else begin
  
Matrix_B[0,0]=(field_001[0]-field_000[0])/(delta*2.0d)
Matrix_B[0,1]=(field_011[0]-field_010[0])/(delta*2.0d)
Matrix_B[0,2]=(field_101[0]-field_100[0])/(delta*2.0d)

Matrix_B[1,0]=(field_001[1]-field_000[1])/(delta*2.0d)
Matrix_B[1,1]=(field_011[1]-field_010[1])/(delta*2.0d)
Matrix_B[1,2]=(field_101[1]-field_100[1])/(delta*2.0d)

Matrix_B[2,0]=(field_001[2]-field_000[2])/(delta*2.0d)
Matrix_B[2,1]=(field_011[2]-field_010[2])/(delta*2.0d)
Matrix_B[2,2]=(field_101[2]-field_100[2])/(delta*2.0d)

endelse
return,Matrix_B
end
