PRO CALCULATE_FZ_ANY,Bl,Btr,phi_los,Bz,Bh,phi,B0,L0,P,Bc,Lc,pixel_size,sig_Bl,$
sig_Btr,sig_phil,Jz,dB_zp,Fz,Jhp,err_Jz,err_dBz,err_Fz,err_Jhp,quiet=quiet,help=help

; PURPOSE :   Calculate the vertical Lorentz force and the minimum
;             cross-field electric current density from any vector magnetogram
;             (Calculation is made with the heliographic components on
;              the image plane)
; PROGRAMMER: Manolis K. Georgoulis (JHU/APL) - May/June 2004
;
; CALLING SEQUENCE: 
; CALCULATE_FZ_ANY,Bl,Btr,phi_los,Bz,Bh,phi,B0,L0,P,Bc,Lc,pixel_size,$
;                  sig_Bl,sig_Btr,sig_phil,Jz,dB_zp,Fz,Jhp,err_Jz,err_dBz,$
;                  err_Fz,err_Jhp,/quiet,/help
;
; INPUT:
; Bl --> Longitudinal magnetic field (G)
; Btr --> Transverse magnetic field (G)
; phi_los --> Line-of-sight azimuth angle (radians)
; Bz --> Vertical magnetic field (G)
; Bh --> Horizontal magnetic field (G)
; phi --> Azimuth angle (radians)
; B0 --> Latitude of solar disk center (radians)
; L0 --> Longitude of solar disk center (radians)
; P --> Position of the northern extremity of solar rotation axis
;       (radians)
; Bc --> Heliographic latitude at every pixel of the magnetogram
;        (radians)
; Lc --> Heliographic longitude at every pixel of the magnetogram
;        (radians)
; pixel_size --> Pixel size in arcsec
; sig_Bl -->  sigma of Bl - can be a single value or an array (Gauss)
; sig_Btr --> sigma of Btr - can be a single value or an array (Gauss)
; sig_phil --> sigma of phi_los - can be a single value or an array (degrees)
;
; OUTPUT:
; Jz --> Vertical current density (statAmpere / cm^2)
; dB_zp --> Minimum-structure vertical magnetic field gradient (G/cm)
; Fz --> Vertical Lorentz force (dyn/cm^2)
; Jhp --> Minimum cross-field electric current density (statAmpere /
;         cm^2)
; err_Jz --> Relative error of Jz
; err_dBz --> Relative error of dB_zp
; err_Fz --> Relative error of Fz
; err_Jhp --> Relative error of Jhp
;
; KEYWORDS:
; /quiet -->    If set, it suppresses informational messages regarding
;               the mean and median values and the mean and median
;               errors of Jz,Fz, and Jhp
; /help -->     Use this keyword by itself to obtain the calling
;               sequence and the above information 
;
if keyword_set(help) then begin 
  print,'  '
  print,'SYNTAX:'
  print,'CALCULATE_FZ_ANY,Bl,Btr,phi_los,Bz,Bh,phi,B0,L0,P,Bc,Lc,pixel_size,$'
  print,'sig_Bl,sig_Btr,sig_phil,Jz,dB_zp,Fz,Jhp,err_Jz,err_dBz,err_Fz,err_Jhp,/help'
  print,'  '
  print,'INPUT:'
  print,'Bl --> Longitudinal magnetic field (G)'
  print,'Btr --> Transverse magnetic field (G)'
  print,'phi_los --> Line-of-sight azimuth angle (radians)'
  print,'Bz --> Vertical magnetic field (G)'
  print,'Bh --> Horizontal magnetic field (G)'
  print,'phi --> Azimuth angle (radians)'
  print,'B0 --> Latitude of solar disk center (radians)'
  print,'L0 --> Longitude of solar disk center (radians)'
  print,'P --> Position of the northern extremity of solar rotation axis (radians)'
  print,'Bc --> Heliographic latitude at every pixel of the magnetogram (radians)'
  print,'Lc --> Heliographic longitude at every pixel of the magnetogram (radians)'
  print,'pixel_size --> Pixel size in arcsec'
  print,'sig_Bl -->  sigma of Bl (Gauss)'
  print,'sig_Btr --> sigma of Btr (Gauss)'
  print,'sig_phil --> sigma of phi_los (degrees)'
  print,'  '
  print,'OUTPUT:'
  print,'Jz --> Vertical current density (statAmpere / cm^2)'
  print,'dB_zp --> Minimum-structure vertical magnetic field gradient (G/cm)'
  print,'Fz --> Vertical Lorentz force (dyn/cm^2)'
  print,'Jhp --> Minimum cross-field electric current density (statAmpere / cm^2)'
  print,'err_Jz --> Relative error of Jz'
  print,'err_dBz --> Relative error of dB_zp'
  print,'err_Fz --> Relative error of Fz'
  print,'err_Jhp --> Relative error of Jhp'
  print,'  '
  print,'KEYWORDS:'
  print,'/QUIET --> If set, it suppresses informational messages about the results'
  print,'/HELP --> Displays this message'
  print,'  '
  return
endif

sig_phil=sig_phil*!dtor
r1=where(abs(Bz) gt 5.d3,ico) & if ico gt 0. then Bz(r1)=0.
r1=where(abs(Bh) gt 5.d3,ico) & if ico gt 0. then Bh(r1)=0.
res=size(Bz) & id1=res(1) & id2=res(2)
lamda=pixel_size*7.25d7
Bx=Bh*cos(phi) & By=Bh*sin(phi)

CALCULATE_HELIO_ERRORS_ANY,Bl,Btr,phi_los,Bz,Bh,phi,B0,L0,P,Bc,Lc,sig_Bl,$
sig_Btr,sig_phil,sig_Bh,sig_Bz,sig_B
DO_CALC,Bx,By,Bz,dB_zp,Fz,Jhp,lamda
CALCULATE_VERTICAL_CURRENT,Bx,By,lamda,Jz
err_Jz=fltarr(id1,id2) 
r1=where(Bh ne 0.) & err_Jz(r1)=sig_Bh(r1)/Bh(r1)
CALCULATE_ERRORS,sig_Bz,sig_Bh,sig_B,Bz,Bh,err_dBz,err_Fz,err_Jhp,lamda

if keyword_set(quiet) eq 0. then begin 

per_thres=0.3
r1=where(err_Fz lt per_thres and err_Jhp lt per_thres and err_Jz lt per_thres $
and err_Fz*err_Jhp*err_Jz ne 0.)
Jhp_mean=mean(Jhp(r1))*1.d3/3.d5
Jhp_med=median(Jhp(r1))*1.d3/3.d5
mean_error_Jhp=mean(err_Jhp(r1))*Jhp_mean
median_error_Jhp=median(err_Jhp(r1))*Jhp_mean
Fz_mean=mean(abs(Fz(r1)))
Fz_med=median(abs(Fz(r1)))
mean_error_Fz=mean(err_Fz(r1))*Fz_mean
median_error_Fz=median(err_Fz(r1))*Fz_mean
Jz_mean=mean(abs(Jz(r1)))*1.d3/3.d5
Jz_med=median(abs(Jz(r1)))*1.d3/3.d5
mean_error_Jz=mean(err_Jz(r1))*Jz_mean
median_error_Jz=median(err_Jz(r1))*Jz_mean

print,'  '
print,strcompress('  Max Tolerable error in cross-field currents:'+$
                  string(fix(per_thres*100))+'%')
print,'  '
print,'Cross-field electric current density'
print,strcompress('  Mean Jhp '+string(Jhp_mean)+' mA/m2')
print,strcompress('  Median Jhp '+string(Jhp_med)+' mA/m2')
print,strcompress('  Mean error '+string(mean_error_Jhp)+' mA/m2')
print,strcompress('  Median error '+string(median_error_Jhp)+' mA/m2')
print,'  '
print,'Vertical Lorentz force:'
print,strcompress('  Mean Fz '+string(Fz_mean)+' dyn/cm2')
print,strcompress('  Median Fz '+string(Fz_med)+' dyn/cm2')
print,strcompress('  Mean error '+string(mean_error_Fz)+' dyn/cm2')
print,strcompress('  Median error '+string(median_error_Fz)+' dyn/cm2')
print,'  '
print,'Vertical electric current density'
print,strcompress('  Mean Jz '+string(Jz_mean)+' mA/m2')
print,strcompress('  Median Jz '+string(Jz_med)+' mA/m2')
print,strcompress('  Mean error '+string(mean_error_Jz)+' mA/m2')
print,strcompress('  Median error '+string(median_error_Jz)+' mA/m2')
endif

END
