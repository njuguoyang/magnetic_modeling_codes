PRO INTERPOLATE_NPFC,im,hel,idim1,idim2,iflag,fit_mask=fit_mask

; PURPOSE: From an image im on the image plane interpolate into an
; image hel on the heliographic plane, where idim1,idim2 are the
; linear dimensions of the image plane
; PROGRAMMER: Manolis K. Georgoulis (JHU/APL, 10/13/05)

COMMON PAR1,id1,id2,xi,yi,ksi,eta,mask
COMMON PAR2,iapp_toth

hel=dblarr(id1,id2)
hel(fix(xi+0.5),fix(yi+0.5))=im(fix(ksi+0.5),fix(eta+0.5))

if iflag eq 0 then begin 
  mask=fltarr(id1,id2)
  q1=where(hel eq 0.,count) & if count gt 0. then mask(q1)=1.
  tmp=mask
  ms=shift(mask,-1,0)+shift(mask,1,0)+shift(mask,0,-1)+shift(mask,0,1)+ $
     shift(mask,-1,-1)+shift(mask,-1,1)+shift(mask,1,-1)+shift(mask,1,1)
  q1=where(mask eq 1. and ms eq 8,count) & if count gt 0. then tmp(q1)=0.
  mask=tmp
endif 
hel_in=hel

r1=where(mask eq 1.,ico)
if ico gt 0. then for m=0,49 do $
hel(r1)=(1./4.)*(hel(r1-1)+hel(r1+1)+hel(r1+id1)+hel(r1-id1))
if keyword_set(fit_mask) eq 0. then hel=hel*iapp_toth

RETURN
END
