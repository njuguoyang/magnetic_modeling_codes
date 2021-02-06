PRO MAKE_FRAME_NEW,hel,iapp,frame,iapp_ref,idim1,idim2 

; PURPOSE: Use any mask iapp_ref of the image plane to interpolate
; into a mask iapp on the heliographic plane and find its boundary
; frame, where idim1,idim2 are the linear dimensions of the image
; plane 
; PROGRAMMER: Manolis K. Georgoulis (JHU/APL, 10/13/05)

COMMON PAR1,id1,id2,xi,yi,ksi,eta,mask

iapp=fltarr(id1,id2)
for i=0,idim1-1 do for j=0,idim2-1 do if iapp_ref(i,j) eq 1. then $
    iapp(fix(xi(i,j)+0.5),fix(yi(i,j)+0.5))=1.
r1=where(mask eq 1. and shift(iapp,-1,0)+shift(iapp,1,0)+shift(iapp,0,-1)+ $
         shift(iapp,0,1) ne 0.,ic) & if ic gt 0. then iapp(r1)=1. 
iapp(0,*)=0.
iapp(id1-1,*)=0.
iapp(*,0)=0.
iapp(*,id2-1)=0.
tmp=iapp
r1=where(iapp eq 1. and shift(iapp,-1,0)+shift(iapp,1,0)+shift(iapp,0,-1)+ $
         shift(iapp,0,1) lt 3.,ic) & if ic gt 0. then tmp(r1)=0.
iapp=tmp
frame=fltarr(id1,id2)
r1=where(iapp eq 0. and shift(iapp,-1,0)+shift(iapp,1,0)+shift(iapp,0,-1)+ $
         shift(iapp,0,1) ne 0.,ic)
if ic gt 0. then frame(r1)=1.

RETURN
END
