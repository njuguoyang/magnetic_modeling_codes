PRO IVM_GET_NEW,raw,Bz,Btr,phi_unres,iapp_ref,iapp_tot,pixel_size,Bc,$
Lc,Bcc,Lcc,B0,L0,P,sig_Bl,sig_Btr,hms,ff=ff

; PURPOSE: From a IVM-type structure RAW, extract all the neseccary
; information for the calculations, namely
; Bz --> The longitudinal magnetic field
; Btr --> The transverse magnetic field 
; phi_unres --> The azimuth of Btr, 0-POINT is NORTHWARD (UPWARD), counterclockwise, in degree or radian
; iapp_ref --> A mask of the strong-field locations, i.e. the
;              locations where Bz > hms*sig_Bl and Btr > hms*sig_Btr
;              or where the total field strength is > 500 G
; iapp_tot --> A mask of the nonzero field locations
; pizel_size --> Pixel size in arcsec
; Bc,Lc --> Heliographic coordinates of the image plane
; Bcc,Lcc --> Heliographic coordinates of the image plane center
; B0, L0 --> Heliographic coordinates of the solar disk center
; P --> The solar P-angle
; sig_Bl --> Uncertainty in the longitudinal field Bz
; sig_Btr --> Uncertainty in the transverse field Btr
; hms --> Sigma threshold

; NOTE: FOR AN INPUT STRUCTURE WITH DIFFERENT FORMAT, A DIFFERENT
; ROUTINE HAS TO BE CREATED!
; PROGRAMMER: Manolis K. Georgoulis (JHU/APL, 10/13/05)

Bl=double(raw.B_long)
inf=size(Bl) & idim1=inf(1) & idim2=inf(2)
Btr=double(raw.B_trans)
Bz=Bl
; INCLUDE FILLING FACTOR
if keyword_set(ff) then begin 
  Bfill=double(raw.B_fill) ;  IVM Format that I use
;  Bfill=ff                  ;  Format for the ambiguity workshop
  Bz=Bl*Bfill
  Btr=Btr*sqrt(Bfill)
endif
phi_unres=double(raw.B_azim)
if max(abs(phi_unres)) gt 2.*!dpi then phi_unres=phi_unres*!dtor
q1=where(phi_unres lt 0.,ico) & if ico gt 0. then $ 
   phi_unres(q1)=2.*!dpi + phi_unres(q1)
;
; IVM 0-POINT IN AZIMUTH: NORTHWARD (UPWARD) - CORRECTION BELOW
phi_unres=!dpi/2.+phi_unres
r1=where(phi_unres eq 0. or phi_unres gt 2.*!dpi,ico) & if ico gt 0. then $ 
   phi_unres(r1)=phi_unres(r1)-2.*!dpi
;
pixel_size1=double(raw.point.pix_size(0))
pixel_size2=double(raw.point.pix_size(1))
if pixel_size1 eq pixel_size2 then pixel_size=pixel_size1 else $ 
message,'DO SOMETHING WITH THE PIXEL SIZE!'

B=sqrt(Btr^2.+Bz^2.)
Bref=B
iapp=fltarr(idim1,idim2)
r1=where(abs(Bz) gt hms*sig_Bl and Btr gt hms*sig_Btr,ico) 
if ico gt 0. then iapp(r1)=1.
r1=where(iapp eq 0. and B gt 500.,ico) & if ico gt 0. then iapp(r1)=1.
r1=where(B gt 5.d3,ico) 
if ico gt 0. then begin 
  Bz(r1)=0. & Btr(r1)=0. & phi_unres(r1)=0.
  iapp(r1)=0. & iapp(r1-1)=0. & iapp(r1-1)=0. & iapp(r1+idim1)=0.
  iapp(r1-idim1)=0.
endif
r1=where(iapp eq 1. and shift(B,-1,0)*shift(B,1,0)*shift(B,0,-1)*shift(B,0,1) $
         eq 0.,ico) & if ico gt 0. then iapp(r1)=0.
iapp(0,*)=0. & iapp(idim1-1,*)=0. & iapp(*,idim2-1)=0. & iapp(*,0)=0.
B=Bref & iapp_ref=iapp

str_res=tag_names(raw)
r1=where(str_res eq 'LATITUDE',ico)
if ico gt 0. then begin 
  Bc=double(raw.latitude)*!dtor
  Lc=double(raw.cmd)*!dtor
endif else begin 
  Bc=double(raw.point.lat)
  Lc=double(raw.point.cmd)
endelse
Bcc=double(raw.point.lat)
Lcc=double(raw.point.cmd)
B0=double(raw.point.B0)
L0=0.
P=double(raw.point.P)

fitmask=fltarr(idim1,idim2) 
r1=where(B ne 0.) & fitmask(r1)=1.

iapp_tot=fltarr(idim1,idim2) & iapp_tot(*,*)=1.
iapp_tot=iapp_tot*fitmask
iapp_tot(0,*)=0. & iapp_tot(idim1-1,*)=0. & iapp_tot(*,0)=0. & iapp_tot(*,idim2-1)=0.

RETURN
END




