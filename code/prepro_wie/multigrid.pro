;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;written by Thomas Wiegelmann 03.05.2006
;AIMS: The program reads vectormagnetogram
;      data (as example for Aads model) and prepares
;      all input file needed for a nonlinear force-free
;      extrapolation with a multigrid-like version of
;      the Optimization code:
;      Wiegelmann,Solar Physics, Vol. 219: 87-108, 2004.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; this program reads the vectormagnetogram for
; Aads model, as prepared by Tom Metcalf
; written by Thomas Wiegelmann
; The files grid1-3.ini and alboundaries1-3.dat
; are used for the nonlinear force-free code
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO multigrid,folder=folder,file=file,level=level,prep=prep
common vecB,bx,by,bz

if n_elements(folder) eq 0 then folder=''
if n_elements(file) eq 0 then file=''
restore,folder+file,/verbose
nx=(size(bz))[1] & ny=(size(bz))[2] & nz=Min([nx,ny]) & nd=Min([nx,ny])/8
if n_elements(level) eq 0 then level=3
;3 stands for 3 level multigrid
;But this file has been modified to output only one grid, say, level=3 then 
;nx=nx/4, level=1 then nx=nx. Modified by GUO Yang, 2007-2008
if n_elements(prep) eq 0 then prep=1
;prep=1 means that we preprocess the magnetogram on each level
;prep=0 meens that the vectormagnetogram is
;NOT preprocessed (for test cases, like Aads model)
;set the second parameter to 1 to preprocess the data
;See Wiegelmann et al. (2006): Sol. Phys., 233, p.215

reducer=2^(level-1)
if ((nd mod reducer) ne 0) then nd=nd-(nd mod reducer)
help, folder,nx,ny,nz,nd

; test if nx,ny,nz,nd are multiples of reducer
dummy=(nx mod reducer) +(ny mod reducer) + (nz mod reducer) +  (nd mod reducer)

if (dummy ne 0) then print,'nx,ny,nz,nd must be multiples of',reducer
if (dummy eq 0) then begin
print,'Multigrid on', level, ' grids.'

get_lun,u
openw,u,folder+'nd.ini'
printf,u,'nd'
printf,u,nd
printf,u,'nxmax'
printf,u,nx
printf,u,'nymax'
printf,u,ny
printf,u,'nzmax'
printf,u,nz
close,u
free_lun,u

;Rebin grid
nx=nx/reducer & ny=ny/reducer & nz=nz/reducer & nd=nd/reducer
bx=rebin(bx,nx,ny)
by=rebin(by,nx,ny)
bz=rebin(bz,nx,ny)
Lx=1.0*(nx-1) & Ly=1.0*(ny-1) & Lz=1.0*(nz-1)

print,'Processing grid',nx,ny,nz,nd
; Preprocessing yes or no? (prep=1 means yes)
if (prep eq 1) then begin
  bxorig=bx & byorig=by & bzorig=bz
  dummy=prepro(0.001,0.01,nx,ny,nz)
  ; here mu3=0.001 and mu4=0.01
  ; These parameters are for preprocessing of vectormagnetograms
  ; The method is described in:
  ; T. Wiegelmann, B. Inhester, T.Sakurai:
  ; Preprocessing of vector magnetograph data
  ; for a nonlinear force-free magnetic field reconstruction.
  ; (Solar Physics, Vol. 233, 215-232, 2006.)
  ; Please note:
  ; Test cases like Low&Lou, Titov&Demoulin or Aads model
  ; do NOT need preprocessing
endif else begin
  print,'Magnetogram is not preprocessed'
endelse
get_lun,u
gridstring='grid'+'.ini'
openw,u,folder+gridstring
printf,u,'nx'
printf,u,nx
printf,u,'ny'
printf,u,ny
printf,u,'nz'
printf,u,nz
printf,u,'mu'
printf,u,0.1
printf,u,'nd'
printf,u,nd
close,u
allstring='allboundaries'+'.dat'
print,'Saving: ',gridstring
print,'Saving: ', allstring
openw, u, folder+allstring
for iy=0,ny-1 do begin
for ix=0,nx-1 do begin
printf,u, bx[ix,iy]
printf,u, by[ix,iy]
printf,u, bz[ix,iy]
endfor
endfor
close,u
free_lun,u
if (prep eq 1) then begin
bx=bxorig & by=byorig & bz=bzorig
endif

ENDIF
END
