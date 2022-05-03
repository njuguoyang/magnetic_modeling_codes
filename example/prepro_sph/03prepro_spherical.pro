;+
;
; PURPOSE : pre-process bottom boundaries in spherical geometry
;
; CATEGORY :
; CALLING SEQUENCE :
; INPUTS :
;   vector magnetic field projected in the spherical coordinates
; OPTIONAL INPUT PARAMETERS (KEYWORD PARAMETERS):
;   mu3 : coefficient of L3
;   mu4 : coefficient of L4
;   LL  : save L12, L3, L4
;   itr : iteration steps
;   dt  : steps to print information
;   scl : scale multiplied to mu1, mu2, mu3, and mu4 to control the iteration step size
;
; OUTPUTS :
;    preprocessed vector magentic field
; COMMON BLOCKS :
; SIDE EFFECTS :
; RESTRICTIONS :
; PROCEDURE :
; MODIFICATION HISTORY :
;   2012.03 Guo Yang@Nanjing University
;
;   This software is provided without any warranty. Permission to use,
;   copy, modify. Distributing modified or unmodified copies is granted,
;   provided this disclaimer and information are included unchanged.
;-


FUNCTION laplace2,field
return,-4*field+shift(field,1,0)+shift(field,-1,0)+shift(field,0,1)+shift(field,0,-1)
END

pro 03prepro_spherical,mu3=mu3,mu4=mu4,LL=LL,itr=itr,dt=dt,scl=scl

;  read noisy boundary data
;filename='.//field0.dat'                       ; high resolution with 0.03 degrees per cell
filename='./field0_midresolution.dat'           ; middle resolution with 0.06 degrees per cell
nr   = 1
nlat = 1
nlon = 1
openr, unit, filename, /get_lun
readu, unit, nr, nlat, nlon
print,nr,nlat,nlon
;nlon=2400
nlon=1201
Br0 = fltarr(nlat,nlon)
Bt0 = fltarr(nlat,nlon)
Bp0 = fltarr(nlat,nlon)
rix   = fltarr(nr)
theta = fltarr(nlat)
phi   = fltarr(nlon)
readu, unit, Br0
readu, unit, Bt0
readu, unit, Bp0
;print,'minmax Brtp:',minmax(Br0),minmax(Bt0),minmax(Bp0)
readu, unit, rix
readu, unit, theta
readu, unit, phi
free_lun, unit

if not (keyword_set(scl)) then begin
  scl=1.0d
endif else begin
  scl=scl
endelse
mu1=1.0d  ; Force
mu2=1.0d  ; Torque
if not (keyword_set(mu3)) then begin
  mu3=0.01d  ; Observations
endif else begin
  mu3=mu3
endelse
if not (keyword_set(mu4)) then begin 
  mu4=0.1d  ; Smooth
endif else begin
  mu4=mu4
endelse
mu1 = mu1*scl
mu2 = mu2*scl
mu3 = mu3*scl
mu4 = mu4*scl
print,'    mu1, mu2, mu3, mu4:'
print,mu1, mu2, mu3, mu4
if not (keyword_set(itr)) then begin
  itr=5000 
endif else begin
  itr=itr
endelse
if not (keyword_set(dt)) then begin
  dt=200
endif else begin
  dt=dt
endelse

;  Correct the singular point at the north and south poles
;theta[0] = 0.5*theta[1]
;theta[nlat-1] = theta[nlat-2] + 0.5*theta[1]
;Br0[0,*] = 0.5*(Br0[0,*] + Br0[1,*])
;Bt0[0,*] = 0.5*(Bt0[0,*] + Bt0[1,*])
;Bp0[0,*] = 0.5*(Bp0[0,*] + Bp0[1,*])
;Br0[nlat-1,*] = 0.5*(Br0[nlat-1,*] + Br0[nlat-2,*])
;Bt0[nlat-1,*] = 0.5*(Bt0[nlat-1,*] + Bt0[nlat-2,*])
;Bp0[nlat-1,*] = 0.5*(Bp0[nlat-1,*] + Bp0[nlat-2,*])

;nb  = 3
;Br0 = Br0[nb:nlat-nb-1,*]
;Bt0 = Bt0[nb:nlat-nb-1,*]
;Bp0 = Bp0[nb:nlat-nb-1,*]
;theta = theta[nb:nlat-nb-1]
print,'    theta0:',theta[0]*!radeg
;nlat = nlat - 2*nb

bave2d = total(sqrt(Br0^2+Bt0^2+Bp0^2))/N_ELEMENTS(Br0)
Br0 = Br0/bave2d
Bt0 = Bt0/bave2d
Bp0 = Bp0/bave2d
Bro = Br0
Bto = Bt0
Bpo = Bp0
print,'    the flux balance coefficient eps_flux:',total(Br0)/total(abs(Br0))
dth = mean(theta[1:nlat-1] - theta[0:nlat-2])
dph = mean(phi[1:nlon-1] - phi[0:nlon-2])
print,'    dth*dph:',dth*dph

;  Construct the meshes
theta2 = fltarr(nlat,nlon)
for i=0,nlon-1 do begin
  theta2[*,i] = theta
endfor
phi2 = fltarr(nlat,nlon)
for i=0,nlat-1 do begin
  phi2[i,*] = phi
endfor

;  Integration of the magnetic field energy on the bottom surface. The constant Delta_theta*Delta_phi is omitted
Eb = total(sin(theta2)*(Bt0*Bt0 + Bp0*Bp0 + Br0*Br0))
print,'    Eb:',Eb

it=-1
L=9000.0d
oldL=9999.0d
dL=1.0d
while ((it lt itr) and (dL gt 1.0e-4) and (L lt oldL)) do begin
  it = it + 1
  if (it gt 0) then begin
    oldL = L
  endif

  ;  Construct some variables. Refer to T. Tadesse et al. (2009) A&A, 508, 421
  Ebm = 0.5*(Bt0*Bt0 + Bp0*Bp0 - Br0*Br0)
  B1  = Bt0*cos(theta2)*cos(phi2) - Bp0*sin(phi2)
  B2  = Bt0*cos(theta2)*sin(phi2) + Bp0*cos(phi2)
  B3  = Bp0*cos(theta2)*cos(phi2) + Bt0*sin(phi2)
  B4  = Bp0*cos(theta2)*sin(phi2) - Bt0*cos(phi2)

  ;  force
  term1a = total(sin(theta2)*(Ebm*sin(theta2)*cos(phi2) - Br0*B1))/Eb
  term1b = total(sin(theta2)*(Ebm*sin(theta2)*sin(phi2) - Br0*B2))/Eb
  term1c = total(sin(theta2)*(Ebm*cos(theta2) + Br0*Bt0*sin(theta2)))/Eb

  ;  torque
  term2a = total(sin(theta2)*Br0*B3)/Eb
  term2b = total(sin(theta2)*Br0*B4)/Eb
  term2c = total(sin(theta2)*sin(theta2)*Br0*Bp0)/Eb

  ;  Observations
  term3a = (Br0 - Bro)/Eb
  term3b = (Bt0 - Bto)/Eb
  term3c = (Bp0 - Bpo)/Eb

  ;  smooth
  term4a = laplace2(Br0)/Eb
  term4b = laplace2(Bt0)/Eb
  term4c = laplace2(Bp0)/Eb

  L1 = (term1a^2+term1b^2+term1c^2)
  eps_force = abs(term1a) + abs(term1b) + abs(term1c)
  L2 = (term2a^2+term2b^2+term2c^2)
  eps_torque= abs(term2a) + abs(term2b) + abs(term2c)
  L12= L1+L2
  L3 = total((term3a^2+term3b^2+term3c^2))
  L4 = total((term4a^2+term4b^2+term4c^2))
  L=mu1*L1+mu2*L2+mu3*L3+mu4*L4

  if (it gt 0) then begin
    dL = abs(L12-oldL12)/L12+abs(L3-oldL3)/L3+abs(L4-oldL4)/L4
  endif
  oldL12=L12
  oldL3=L3
  oldL4=L4
  if (it mod dt) eq 0 then begin
    if (it eq 0) then print,'     it,             dL,         L12,          L3,          L4,      L=mu1*L12+mu3*L3+mu4*L4,    eps_force,    eps_torque:'
    print,it, dL, L12, L3, L4, L, eps_force, eps_torque
  endif

  Br1 = Br0 - 2.0*mu3*term3a $
            - 2.0*mu4*laplace2(term4a)

  Bt1 = Bt0 - 2.0*mu1*term1a*(Bt0*sin(theta2)*sin(theta2)*cos(phi2) - Br0*sin(theta2)*cos(theta2)*cos(phi2)) $
            - 2.0*mu1*term1b*(Bt0*sin(theta2)*sin(theta2)*sin(phi2) - Br0*sin(theta2)*cos(theta2)*sin(phi2)) $
            - 2.0*mu1*term1c*(Bt0*sin(theta2)*cos(theta2) + Br0*sin(theta2)*sin(theta2)) $
            - 2.0*mu2*(term2a*Br0*sin(theta2)*sin(phi2) - term2b*Br0*sin(theta2)*cos(phi2)) $
            - 2.0*mu3*term3b $
            - 2.0*mu4*laplace2(term4b)

  Bp1 = Bp0 - 2.0*mu1*term1a*(Bp0*sin(theta2)*sin(theta2)*cos(phi2) + Br0*sin(theta2)*sin(phi2)) $
            - 2.0*mu1*term1b*(Bp0*sin(theta2)*sin(theta2)*sin(phi2) - Br0*sin(theta2)*cos(phi2)) $
            - 2.0*mu1*term1c*(Bp0*sin(theta2)*cos(theta2)) $
            - 2.0*mu2*(term2a*Br0*cos(theta2)*cos(phi2)*sin(theta2) + term2b*Br0*cos(theta2)*sin(phi2)*sin(theta2) + term2c*Br0*sin(theta2)*sin(theta2)) $
            - 2.0*mu3*term3c $
            - 2.0*mu4*laplace2(term4c)
  Br0 = Br1
  Bt0 = Bt1
  Bp0 = Bp1

endwhile

LL=[0.0,0.0,0.0]
LL=[oldL12,oldL3,oldL4]

Br0 = Br0*bave2d
Bt0 = Bt0*bave2d
Bp0 = Bp0*bave2d
print,'Correct magnetogram finished'

device,decomposed=0
loadct,4
;wdef,1,800,800
window,xs=600,ys=600
tmp = congrid(Br0,600,600)
;contour, rotate(congrid(Br0,600,600),4), congrid(phi,600)*!radeg, 90-congrid(theta,600)*!radeg, nlevels=55,xstyle=1,ystyle=1,/fill,zrange=[-430,430]
print,min(tmp),max(tmp)
ind = where(tmp le -1500)
tmp[ind] = -1500
ind = where(tmp ge 1500)
tmp[ind] = 1500
tvscl,rotate(tmp,3)
write_png,'Brad.png',tvrd(true=1)

tmp = congrid(Bt0,600,600)
;contour, rotate(congrid(Bt0,600,600),4), congrid(phi,600)*!radeg, 90-congrid(theta,600)*!radeg, nlevels=55,xstyle=1,ystyle=1,/fill,zrange=[-300,300]
print,min(tmp),max(tmp)
ind = where(tmp le -1000)
tmp[ind] = -1000
ind = where(tmp ge 1000)
tmp[ind] = 1000
tvscl,rotate(tmp,3)
write_png,'Bthe.png',tvrd(true=1)

tmp = congrid(Bp0,600,600)
;contour, rotate(congrid(Bp0,600,600),4), congrid(phi,600)*!radeg, 90-congrid(theta,600)*!radeg, nlevels=55,xstyle=1,ystyle=1,/fill,zrange=[-300,300]
print,min(tmp),max(tmp)
ind = where(tmp le -1300)
tmp[ind] = -1300
ind = where(tmp ge 1300)
tmp[ind] = 1300
tvscl,rotate(tmp,3)
write_png,'Bphi.png',tvrd(true=1)

;  save in ASCII file
filename='field0_preprocessed.dat'
openw, unit, filename, /get_lun
writeu, unit, nr, nlat, nlon
writeu, unit, Br0
writeu, unit, Bt0
writeu, unit, Bp0
writeu, unit, rix
writeu, unit, theta
writeu, unit, phi
free_lun, unit

end

