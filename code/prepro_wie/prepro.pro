FUNCTION prepro,mu3,mu4,nx,ny,nz,NOPRINT=noprint
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Written by Thomas Wiegelmann 04.04.2006;Aim: Preprocessing of vectormagnetograms for
;a nonlinear force-free reconstruction
;
; The method is described in:
; T. Wiegelmann, B. Inhester, T.Sakurai:
; Preprocessing of vector magnetograph data
; for a nonlinear force-free magnetic field reconstruction.
; (Solar Physics, Vol. 233, 215-232, 2006.)
;
; See there for the meaning of the parameters
; mu3 (data) and mu4 (smoothing).
; often used values are in the range of mu3, mu4 <~0.1
; too large values of mu3 and mu4 might cause problems.
;
; Please note that the vector magnetogram must exist as
; global variables "Bx, By, Bz" in the common block "vecB"
; This program overwrites Bx,By,Bz with the preprocessed
; data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
common vecB,bx,by,bz
if (n_elements(bx) eq 0) or  (n_elements(by) eq 0) or $
 (n_elements(bz) eq 0) then begin
print, 'Global 2D-fields Bx,By,Bz must exist and have same size.'
stop
endif
if n_elements(noprint) eq 0 then noprint=0
x1=findgen(nx)/(nx-1)
y1=findgen(ny)/(ny-1)
z1=findgen(nz)/(nz-1)
;
dx=1.0/(nx-1) & dy=1.0/(ny-1) & dxdy=dx*dy
x=fltarr(nx,ny)
y=fltarr(nx,ny)
for iy=0,ny-1 do x[*,iy]=x1
for ix=0,nx-1 do y[ix,*]=y1
; correct_magnetogram
bave2d=total(sqrt(bx^2+by^2+bz^2))/N_ELEMENTS(bx)
bx=bx/bave2d & by=by/bave2d & bz=bz/bave2d

bxo=bx & byo=by &bzo=bz
bx1=bx & by1=by &bz1=bz
;
mu1=1.0  ; Force
mu2=mu1 ; Torque
;mu3=0.1  ; Observations
;mu4=0.1  ; Smooth
;
fac1=10.0
mu1=mu1/fac1 & mu2=mu2/fac1
mu3=mu3/fac1 & mu4=mu4/fac1
if noprint eq 0 then $
print,'mu1,mu2,mu3,mu4:',mu1,mu2,mu3,mu4
emag=total(bx^2+by^2+bz^2)
emag2=emag*emag
;
r=sqrt(x^2+y^2)
ihelp=total(r*(bx^2+by^2+bz^2))
ihelp2=ihelp*ihelp
;
it=-1
L=10 & dL=1.0
print,'dL,eps_force,eps_torque,eps_smooth:'
print,'it,L1,L2,L3,L4:'

while it lt 2000 and dL gt 1.0e-4 do begin
it=it+1
;
;force
;
term1a=total(bx*bz)/emag
term1b=total(by*bz)/emag
term1c=(total(bz^2)-total(bx^2+by^2))/emag
;
;torque
;
term2a=(total(x*bz^2)-total(x*(bx^2+by^2)))/ihelp
term2b=(total(y*bz^2)-total(y*(bx^2+by^2)))/ihelp
term2c=(total(y*bx*bz)-total(x*by*bz))/ihelp
;
;observation
;
term3a=(bx1-bxo)^2/emag
term3b=(by1-byo)^2/emag
term3c=(bz1-bzo)^2/emag
;
;smooth
;
term4a=2.0*(-4.0*laplace(bx)+laplace(shift(bx,1,0))+ $
laplace(shift(bx,-1,0))+laplace(shift(bx,0,1))+laplace(shift(bx,0,-1)))
term4b=2.0*(-4.0*laplace(by)+laplace(shift(by,1,0))+ $
laplace(shift(by,-1,0))+laplace(shift(by,0,1))+laplace(shift(by,0,-1)))
term4c=2.0*(-4.0*laplace(bz)+laplace(shift(bz,1,0))+ $
laplace(shift(bz,-1,0))+laplace(shift(bz,0,1))+laplace(shift(bz,0,-1)))

eps_smooth = total(abs(laplace(bx))) + total(abs(laplace(by))) + total(abs(laplace(bz)))
eps_smooth = eps_smooth*dxdy/sqrt(emag)
;
L1=(term1a^2+term1b^2+term1c^2);^2
eps_force = abs(term1a) + abs(term1b) + abs(term1c)
L2=(term2a^2+term2b^2+term2c^2);^2
eps_torque= abs(term2a) + abs(term2b) + abs(term2c)
L12=L1+L2
L3=total((term3a+term3b+term3c))/emag      
; YG20160214: it seems that L3 does not need to be divided by emag again. But this operation
; does not change the preprocessing of the magnetic field. It only affects L3 itself.
L4=dxdy*total((laplace(bx))^2+(laplace(by))^2+(laplace(bz))^2)/emag
L=L1+L2+L3+L4
;
if it gt 0 then $
dL=abs(L12-oldL12)/L12+abs(L3-oldL3)/L3+abs(L4-oldL4)/L4;

if (it mod 50) eq 0 then $
if noprint eq 0 then print,dL,eps_force,eps_torque,eps_smooth
oldL12=L12 & oldL3=L3 & oldL4=L4
if (it mod 50) eq 0 then $
if noprint eq 0 then  print, it,L1,L2,L3,L4
;if L lt Lold then begin
bx1=bx+mu1*(-2.*term1a*bz+4.*term1c*bx)-mu3*2.*(bx-bxo) $
- mu4*term4a $
+mu2*(4.*term2a*x*bx+4.*term2b*y*bx-2.*term2c*y*bz)
by1=by+mu1*(-2.*term1b*bz+4.*term1c*by)-mu3*2.*(by-byo) $
- mu4*term4b $
+mu2*(4.*term2a*x*by+4.*term2b*y*by+2.*term2c*x*bz)
bz1=bz-mu3*2.*(bz-bzo)- mu4*term4c
bx=bx1
by=by1
bz=bz1
 ;mu=mu*1.1
;Lold=L
;endif else mu=mu/2.
endwhile
;endfor
bx=bx*bave2d & by=by*bave2d & bz=bz*bave2d
if noprint eq 0 then $
print,'Correct magnetogram finished'
subL=[L1+L2,L3,L4]
return,subL
END
;
;
