PRO MATCH_BNP,Bzin,Bz1,Bz2,Bcx,Bcy,lamda,fmask,B0,P,Bcc,Lcc,px,Bx1,Bx2,$
              By1,By2,Bz,mirror=mirror,quiet=quiet

; PURPOSE: Perform and check the convergence of the final distribution
; Bz of the vertical magnetic field. Convergence is checked by means
; of the number of vector flips or, equivalently, changed values of Bz
; in each iterations. The arguments are
; Bzin --> The initial Bz-distribution
; Bz1,Bz2 --> The two possible solutions of Bz on the image plane
; Bcx, Bcy --> The nonpotential magnetic field components (updated in
;              each iteration)
; lamda --> The linear dimension of each pixel (assumed a rectangle)
; fmask --> The mask with the strong-field locations where convergence
;           is checked
; B0 --> Heliographic latitude of the solar disk center
; P --> The solar P_angle
; Bcc,Lcc --> Heliographic coordinates of the center of the image
;             plane
; px --> Pixel size in arcsec
; Bx1,Bx2 --> The two possible solutions of Bx on the image plane
; By1,By2 --> The two possible solutions of By on the image plane
; Bz--> The final (converged) configuration of Bz
; 
; CURRENT SETTINGS:
; Maximum allowed number of iterations (itmax): 300
; Allowed maximum number of flips in each iteration for convergence to
; be assumed complete (min_number): 5
; Number of consecutive iterations in which the number of flips must
; be kept < min_number to assume convergence (iter_num): 10
; THE CONVERGENCE STOPS WHEN THE NUMBER OF FLIPS IS KEPT =< min_number FOR iter_num
; CONSECUTIVE ITERATIONS
;
; NOTE: Locations where the vector flipps for iter_num or more consecutive
; iterations are not counted in the convergence process (by using the
; fitting masks smask and pmask)
;
; PROGRAMMER: Manolis K. Georgoulis (JHU/APL, 10/13/05)

itmax=300.
ic=0.
iterm=0
min_number=5
iter_num=10

res=size(Bz1) & id1=res(1) & id2=res(2)
Bz=dblarr(id1,id2) & Bzr=Bzin 
smask=fltarr(id1,id2) & pmask=fltarr(id1,id2)

tJz_min=1.d20
redo: 
ic=ic+1.  
if ic gt itmax then begin & Bz=Bzr & return &endif

if keyword_set(mirror) then begin 
  GET_POTENTIAL_FIELD,Bzr,B0,P,Bcc,Lcc,px,Bpx,Bpy,/mirror
  ASSIGN_VALUE,Bpx,Bpy,Bcx,Bcy,Bx,Bx1,Bx2,By,By1,By2,Bzr,Bz,Bz1,Bz2,/mirror
  FIND_PROXY,Bx,By,lamda,Jzr,Bcx,Bcy,/mirror
endif else begin 
  GET_POTENTIAL_FIELD,Bzr,B0,P,Bcc,Lcc,px,Bpx,Bpy
  ASSIGN_VALUE,Bpx,Bpy,Bcx,Bcy,Bx,Bx1,Bx2,By,By1,By2,Bzr,Bz,Bz1,Bz2
  FIND_PROXY,Bx,By,lamda,Jzr,Bcx,Bcy
endelse 
r0=where(float(Bz-Bzr)*fmask ne 0.,idiff)
if idiff gt 0 then begin & pmask(r0)=1. & smask(r0)=smask(r0)+1. &endif 
r1=where(smask ge float(iter_num) and pmask eq 1.,ics) & idiff=idiff-ics

if keyword_set(quiet) eq 0. then begin 
  CALCULATE_VERTICAL_CURRENT,Bx,By,lamda,Jz
  tJz=(total(abs(Jz)*fmask)/3.d5)*(lamda/100.)^2. & if tJz lt tJz_min then tJz_min=tJz
  print,'  '
  print,strcompress('Iteration '+string(long(ic)))
  print,strcompress('Total current in this iteration: '+string(tJz)+' Ampere')
  print,strcompress('Minimum achieved total current: '+string(tJz_min)+' Ampere')
  print,strcompress(string(idiff)+' vector flips performed in this iteration')
endif  

Bzr=Bz
if float(total(r0)) ne -1. then pmask(r0)=0.
if idiff gt min_number then begin 
  iterm=0
  goto,redo
endif else begin 
  iterm=iterm+1 
  if iterm lt iter_num then goto,redo 
endelse

RETURN
END




