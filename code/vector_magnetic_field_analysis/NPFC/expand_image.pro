PRO EXPAND_IMAGE,im,mim,ext_x,ext_y,stx,sty,mirror=mirror

; PURPOSE: Enlarge the linear dimensions of an image by a (ext_x,ext_y)
; and put the initial image at the center. If the keyword /mirror is
; set, the additional space corresponds to a mirror image of the
; initial image.
; PROGRAMMER: Manolis K. Georgoulis (JHU/APL, 09/30/05)

res=size(im) & id1=res(1) & id2=res(2)

mim=dblarr(id1+ext_x,id2+ext_y)
stx=fix(float(ext_x)/2.+0.5) 
sty=fix(float(ext_y)/2.+0.5) 
mim(stx:stx+id1-1,sty:sty+id2-1)=im(*,*)

if keyword_set(mirror) then begin   
  if stx le id1 then xmr=stx else xmr=id1
  mr1=im(0:xmr-1,*) & mr1=reverse(mr1)
  mr2=im(id1-xmr:id1-1,*) & mr2=reverse(mr2)
  PLACE_MIRROR,mim,0,stx-1,sty,sty+id2-1,mr1
  PLACE_MIRROR,mim,stx+id1,id1+ext_x-1,sty,sty+id2-1,mr2

  if sty le id2 then ymr=sty else ymr=id2
  mr1=mim(*,ymr:2*ymr-1) & mr1=reverse(mr1,2)
  mr2=mim(*,id2:ymr+id2-1) & mr2=reverse(mr2,2)
  PLACE_MIRROR,mim,0,id1+ext_x-1,0,sty-1,mr1
  PLACE_MIRROR,mim,0,id1+ext_x-1,sty+id2,id2+ext_y-1,mr2
endif

RETURN
END
