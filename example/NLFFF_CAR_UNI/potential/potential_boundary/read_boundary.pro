;pro read_boundary

;=============================================
; Parameter needed to be set for your own case
;=============================================
nx = 244                ; x-size of the preprocessed magnetic field
ny = 132                ; y-size
arcsec2cm=7.32825e7     ; 1 arcsec in cm
xc=759.86469*arcsec2cm  ; x-coordinate of the center of the field of view
yc=-345.34813*arcsec2cm ; y-coordinate of the center
dx=4.0*0.50d*arcsec2cm  ; x spatial resolution of the preprocessed magnetic field
dy=dx                   ; y spatial resolution
;==============================================

x = 0.0
bx=fltarr(nx,ny)
by=fltarr(nx,ny)
bz=fltarr(nx,ny)
get_lun,u
openr,u,'./allboundaries.dat'
for j=0,ny-1 do begin
  for i=0,nx-1 do begin
    readf,u,x
    bx[i,j]=x
    readf,u,x
    by[i,j]=x
    readf,u,x
    bz[i,j]=x
  endfor
endfor
close,u
free_lun,u
tvscl,bz
;v=vector(bx,by)

tmp = Bz[2:nx-3,2:ny-3]
Bz = tmp
sizebz=size(Bz)
nx1=sizebz[1]
nx2=sizebz[2]
filename='potential_boundary.dat'
openw,lun,filename,/get_lun
writeu,lun,nx1
writeu,lun,nx2
writeu,lun,double(xc)
writeu,lun,double(yc)
writeu,lun,double(dx)
writeu,lun,double(dy)
writeu,lun,double(Bz)
free_lun,lun
print,'Bz range (Gauss):', min(Bz),max(Bz)
;wdef,2,800,800
;tvscl,Bz
print,'Computation domain for potential field:'
print,'nx1,nx2',nx1,nx2
print,'xc,yc (cm)',xc,yc
print,'dx,dy (cm)',dx,dy
x1=xc-nx1*dx/2
x2=xc+nx1*dx/2
y1=yc-nx2*dy/2
y2=yc+nx2*dy/2
; output 
print,'x,y, and z range (10 Mm):'
print,'        xprobmin1=',strtrim(string(x1*1.e-9),2),'d0'
print,'        xprobmax1=',strtrim(string(x2*1.e-9),2),'d0'
print,'        xprobmin2=',strtrim(string(y1*1.e-9),2),'d0'
print,'        xprobmax2=',strtrim(string(y2*1.e-9),2),'d0'
print,'        xprobmin3=',strtrim(string(0.0*1.e-9+0.1),2),'d0'   ; to lift the domain 1 Mm above 0
print,'        xprobmax3=',strtrim(string((y2-y1)*1.e-9),2),'d0'

end
