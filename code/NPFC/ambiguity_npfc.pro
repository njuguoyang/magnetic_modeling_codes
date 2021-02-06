PRO AMBIGUITY_NPFC,raw,amb_res,sig_Bl,sig_Btr,hms,mirror=mirror,ff=ff,quiet=quiet,help=help

; PURPOSE :   Azimuth disambiguation in solar vector magnetograms
; FEATURES :  Default imput --> data in a IVM-type structure (using routine IVM_GET_NEW)
;             Similar routine must be used if instrument other than IVM
; PROGRAMMER: Manolis K. Georgoulis (JHU/APL) - March 2005

if keyword_set(help) then begin 
  print,'  '
  print,'CALLING SEQUENCE:' 
  print,'AMBIGUITY_NPFC,raw,amb_res,sig_Bl,sig_Btr,hms,/mirror,/ff,/quiet'
  print,'  '
  print,'ARGUMENTS:'
  print,'raw -->     IVM data structure containing the LOS magnetogram information' 
  print,'amb_res --> output structure containing the results' 
  print,'            (default is a format giving information on both the image and the' 
  print,'             heliographic planes - if the format of the ambiguity resolution'
  print,'             workshop is preferred, simply uncomment the relative section at'
  print,'             the end of the AMBIGUITY_NPFC.PRO routine)'
  print,'sig_Bl -->  uncertainty in longitudinal B above which convergence is checked (Gauss)'
  print,'sig_Btr --> uncertainty in transverse B above which convergence is checked (Gauss)' 
  print,'hms -->     sigma significance; how many sigmas (1,2,3,etc.)'
  print,'  '
  print,'KEYWORDS:'
  print,'/mirror --> The differential and Fourier transform calculations are executed by'
  print,'            padding the original images with 0s. Setting this keyword, the image'
  print,'            extensions are nonzero and contain a mirror of the original image'     
  print,'            HINT: USE MOSTLY WHEN STRONG FLUX RESIDES ON THE IMAGE"S BOUNDARY'
  print,'/ff -->     If set, the IVM filling factor is included in the calculations'
  print,'/quiet -->  Set to suppress informational messages on the convergence process'
  print,'/help -->   Use this keyword by itself to obtain the above information'
  return
endif
;
; MODIFICATION HISTORY:
; * Added keyword /potential for a simple potential-field disambiguation (05/25/05)
; * Added keyword /no_smooth to avoid smoothing of the initial
;   disambiguation solution (05/25/05)
; * Improved the initial ambiguity solution to make them as smooth as
;   possible - routine ARRANGE_PAIRS added (06/03/05) 
; * Improved the methodology - non-potential B-field fixed from proxy
;   current Jz_ref; looking for the optimal distribution of Bz - Old
;   routine can be found in obsolete/ambiguity_npfc_old.pro (06/16/05)
; * Introduced slight smoothing of Bz in MATCH_BNP for off-disk cases
;   and introduced a buffer of 0.1 in the determination of the
;   closer solution for the horizontal field in ASSIGN_VALUE (09/09/05) 
; * (1) Work on the boundary conditions (padding with zeroes or
;       introducing a mirror image (/mirror keyword added)
; * (2) Work on the calculation of the proxy current density. The
;       initial proxy is now given by the formula
;       Jz_pr=Jz(Bl)/|Jz(Bl)|* [sin(Lc)^2. * |Jz(Bl)| + cos(Lc)^2. * |Jz(Btr)|]
;       where Jz(Bl) is the mean vertical current density due to Bl
;       and Jz(Btr) is the Semel & Skumanich (1998) absolute value of
;       the ambiguity-free vertical current density (routine
;       CALCULATE_SS_CURRENT ADDED)
;   (3) The proxy current is now updated in each iteration based on the
;       actual calculated vertical current. Moreover, large current
;       values are ignored in order not to affect the calculation of
;       the nonpotential magnetic field components.
;   (4) The convergence process has been changed to reflect number of
;       vector flips or, equivalently, number of changed Bz-values
;       from iteration to iteration. Process stops when less than 5
;       Bz-values are changed for 10 consecutive iterations
;   (5) A consistency check has been added in the ASSIGN_VALUE routine. The
;       program applies mean neighborhood statistics to determine
;       whether a pair of (Bx,By)-values need to be changed in a given
;       iteration or not. 
;    (Modifications applied between 10/01 and 10/12 2005).

COMMON PAR1,id1,id2,xi,yi,ksi,eta,mask
COMMON PAR2,iapp_toth

; Extracting the necessary info by the IVM structure files - use 
; appropriate routine if instrument is other than IVM
; * IVM_GET can now be applied to BBSO structures created by BBSO_MAG_STRUCTURE
; * 

if keyword_set(ff) then $
IVM_GET_NEW,raw,Bl,Btr,phi_unres,iapp_ref,iapp_tot,pixel_size,Bc,$
Lc,Bcc,Lcc,B0,L0,P,sig_Bl,sig_Btr,hms,/ff else $
IVM_GET_NEW,raw,Bl,Btr,phi_unres,iapp_ref,iapp_tot,pixel_size,Bc,$
Lc,Bcc,Lcc,B0,L0,P,sig_Bl,sig_Btr,hms
lamda=pixel_size*7.25d7
lntd=raw.point.cmd*!radeg
; IVM or BBSO Information extracted

res=size(Bl) & idim1=res(1) & idim2=res(2)
HELIOGRAPHIC_IMAGE_PLANE,Bl,Btr,phi_unres,Bz1_im,Bz2_im,Bh1_im,Bh2_im,$
phi1_im,phi2_im,idim1,idim2,B0,L0,P,Bc,Lc
HELIOGRAPHIC_PLANE,Bz1_im,Bh1_im,phi1_im,Bz1,Bh1,phi1,idim1,idim2,$
B0,L0,P,Bc,Lc,limb,iapp_tot
HELIOGRAPHIC_PLANE,Bz2_im,Bh2_im,phi2_im,Bz2,Bh2,phi2,idim1,idim2,$
B0,L0,P,Bc,Lc,limb,iapp_tot
MAKE_FRAME_NEW,Bz1,iapp,frame,iapp_ref,idim1,idim2
iapp=iapp*iapp_toth
frame_tot=fltarr(idim1,idim2)
r1=where(iapp_tot eq 1. and shift(iapp_tot,-1,0)*shift(iapp_tot,1,0)* $
         shift(iapp_tot,0,-1)*shift(iapp_tot,0,1) eq 0.,ico) 
if ico gt 0. then frame_tot(r1)=1.
Bx1=Bh1_im*cos(phi1_im) & By1=Bh1_im*sin(phi1_im)
Bx2=Bh2_im*cos(phi2_im) & By2=Bh2_im*sin(phi2_im)

Bzin=0.5*(Bz1_im+Bz2_im)
if keyword_set(mirror) then begin 
  FIND_PROXY_INITIAL,Lc,Btr,phi_unres,0.5*(Bx1+Bx2),0.5*(By1+By2),lamda,$
                     Jzpr,Bcx,Bcy,/mirror
  if keyword_set(quiet) then $
  MATCH_BNP,Bzin,Bz1_im,Bz2_im,Bcx,Bcy,lamda,iapp_ref,B0,P,Bcc,Lcc,pixel_size,$
            Bx1,Bx2,By1,By2,Bz_best,/mirror,/quiet else $
  MATCH_BNP,Bzin,Bz1_im,Bz2_im,Bcx,Bcy,lamda,iapp_ref,B0,P,Bcc,Lcc,pixel_size,$
            Bx1,Bx2,By1,By2,Bz_best,/mirror 
endif else begin 
  FIND_PROXY_INITIAL,Lc,Btr,phi_unres,0.5*(Bx1+Bx2),0.5*(By1+By2),lamda,Jzpr,Bcx,Bcy
  if keyword_set(quiet) then $
  MATCH_BNP,Bzin,Bz1_im,Bz2_im,Bcx,Bcy,lamda,iapp_ref,B0,P,Bcc,Lcc,pixel_size,$
            Bx1,Bx2,By1,By2,Bz_best,/quiet else $
  MATCH_BNP,Bzin,Bz1_im,Bz2_im,Bcx,Bcy,lamda,iapp_ref,B0,P,Bcc,Lcc,pixel_size,$
            Bx1,Bx2,By1,By2,Bz_best
endelse
Bz_res_im=Bz_best

FIND_RESPECTIVE_NEW,phi_res_im,phi1_im,phi2_im,Bz_res_im,Bz1_im,Bz2_im
FIND_RESPECTIVE_NEW,Bh_res_im,Bh1_im,Bh2_im,Bz_res_im,Bz1_im,Bz2_im

INTERPOLATE_NPFC,phi_res_im,phi_res_hel,idim1,idim2,1
CLOSER_ANGLE,phi_res_hel,phi1,phi2,phi_res_hel
FIND_RESPECTIVE_NEW,Bz_res_hel,Bz1,Bz2,phi_res_hel,phi1,phi2
FIND_RESPECTIVE_NEW,Bh_res_hel,Bh1,Bh2,phi_res_hel,phi1,phi2
LOS_IMAGE_SOLUTION,phi_res_im,phi1_im,phi2_im,phi_unres,phi_res_los,idim1,idim2
CALCULATE_FZ_ANY,Bl,Btr,phi_res_los,Bz_res_im,Bh_res_im,phi_res_im,B0,L0,P,$
Bc,Lc,pixel_size,sig_Bl,sig_Btr,hms,Jzc,dB_zim,Fz,Jhp,err_Jz,err_dBz,err_Fz,err_Jhp,/quiet
INTERPOLATE_NPFC,Jzc,Jzc_hel,idim1,idim2,1
INTERPOLATE_NPFC,dB_zim,dB_z,idim1,idim2,1
INTERPOLATE_NPFC,Fz,Fz_hel,idim1,idim2,1
INTERPOLATE_NPFC,Jhp,Jhp_hel,idim1,idim2,1

;##################################################################
;                       AMBIGUITY WORKSHOP FORMAT
;##################################################################
;
;Bx=Bh_res_im*cos(phi_res_im)
;By=Bh_res_im*sin(phi_res_im)
;Bz=Bz_res_im
;phi_l=phi_res_los*!radeg - 90.
;r1=where(phi_l lt 0.,ico) & if ico gt 0. then phi_l(r1)=float(phi_l(r1)+360.)
;r1=where(phi_l gt 180.,ico) & if ico gt 0. then phi_l(r1)=phi_l(r1)-360.
;amb_res=create_struct(['I_CONT','LATITUDE','CMD','B_LONG','B_TRANS',$
;'B_AZIM','POINT','BX','BY','BZ'],raw.I_cont,raw.latitude,raw.cmd,$
;float(Bl),float(Btr),float(phi_l),raw.point,float(Bx),float(By),float(Bz))
;return
;####################################################################

if keyword_set(ff) then lab='Filling factor taken into account' else $
lab='Filling factor ignored'
header=strarr(23)
header(0)='Bl --> Longitudinal field'
header(1)='Btr --> Transverse field'
header(2)='phi_los --> LOS solution for the azimuth'
header(3)='Bz --> Solution for the vertical field on the image plane'
header(4)='Bh --> Solution for the horizontal field on the image plane'
header(5)='phi --> Solution for the azimuth on the image plane'
header(6)='Bz_hel --> Solution for the vertical field on the heliographic plane'
header(7)='Bh_hel --> Solution for the horizontal field on the heliographic plane'
header(8)='phi_hel --> Solution for the azimuth on the heliographic plane'
header(9)='Jz --> Vertical current density on the image plane'
header(10)='Jz_hel --> Vertical current density on the heliographic plane'
header(11)='Fz --> Vertical Lorentz force on the image plane'
header(12)='Jhp --> Cross-field horizontal current on the image plane'
header(13)='dB_z --> Partial vertical field gradient on the image plane'
header(14)='dB_z_hel --> Partial vertical field gradient on the heliographic plane'
header(15)='Fz_hel --> Vertical Lorentz force on the heliographic plane'
header(16)='Jhp_hel --> Cross-field horizontal current on the heliographic plane' 
header(17)='err_Jz --> Relative error of Jz on the image plane'
header(18)='err_dBz --> Relative error of dB_z on the image plane'
header(19)='err_Fz --> Relative error of Fz on the image plane'
header(20)='err_Jhp --> Relative error of Jhp on the image plane'
header(21)='fmask_image --> Fitting mask on the image plane'
header(22)='fmask_helio --> Fitting mask on the heliographic plane'
amb_res=create_struct(['Bl','Btr','phi_los','Bz',$
                       'Bh','phi','Bz_hel','Bh_hel',$ 
                       'phi_hel','Jz','Jz_hel','Fz','Jhp','dB_z','dB_z_hel',$;
                       'Fz_hel','Jhp_hel','err_Jz','err_dBz','err_Fz',$
                       'err_Jhp','fmask_image','fmask_helio','ff_impact',$
                       'header'], $
                       float(Bl),float(Btr),float(phi_res_los),float(Bz_res_im),$ 
                       float(Bh_res_im),float(phi_res_im),float(Bz_res_hel),$
                       float(Bh_res_hel),float(phi_res_hel),float(Jzc),$
                       float(Jzc_hel),float(Fz),float(Jhp),float(dB_zim),$
                       float(dB_z),float(Fz_hel),float(Jhp_hel),float(err_Jz),$
                       float(err_dBz),float(err_Fz),float(err_Jhp),$
                       float(iapp_ref),float(iapp),lab,header)


END








