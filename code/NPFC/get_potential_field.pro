PRO GET_POTENTIAL_FIELD,Bz,B0,P,Bcc,Lcc,px,Bpx,Bpy,mirror=mirror

; PURPOSE: For a given distribution Bz of the vertical magnetic field
; calculate the horizontal potential field (Bpx,Bpy) where:
; B0 --> Heliographic latitude of the solar disk center
; P --> The solar P-angle
; Bcc,Lcc --> Heliographic coordinates of the image center
; px --> Pixel size in arcsec
; PROGRAMMER: Manolis K. Georgoulis (JHU/APL, 10/13/05)

res=size(Bz) & id1=res(1) & id2=res(2)

if keyword_set(mirror) then EXPAND_IMAGE,Bz,Bzg,id1,id2,stx,sty,/mirror else $
                            EXPAND_IMAGE,Bz,Bzg,id1,id2,stx,sty 
Bpx=dblarr(id1,id2) & Bpy=dblarr(id1,id2)
b1=lff(Bzg,b0=B0,pangle=P,lat=Bcc,cmd=Lcc,pixel=px,$
       alpha=0.,/normal,/quiet,z=0.)
if n_elements(size(b1)) eq 6 then begin 
  tmp=b1(stx:stx+id1-1,sty:sty+id2-1,*) & b1=tmp
endif else begin 
  tmp=b1(stx:stx+id1-1,sty:sty+id2-1,*,*) & b1=tmp
endelse
Bp1=dblarr(3,id1,id2) & GET_EXTRAP_COMPONENTS,Bp1,b1
Bpx(*,*)=Bp1(0,*,*) & Bpy(*,*)=Bp1(1,*,*)

RETURN
END
