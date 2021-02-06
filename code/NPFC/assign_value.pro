PRO ASSIGN_VALUE,Bpx,Bpy,Bcx,Bcy,Bx,Bx1,Bx2,By,By1,By2,Bzr,Bz,Bz1,Bz2,mirror=mirror

; PURPOSE: For a configuration (Bpx,Bpy,Bzr) of the potential magnetic
; field and a configuration (Bcx,Bcy,0) of the nonpontential field,
; choose the configuration (Bx,By,Bz) of the actual magnetic field from
; the two possible solutions (Bx1,By1,Bz1) and (Bx2,By2,Bz2)
; FEATURE: Apply a mean-neighborhood consistency check to determine
; whether a given vector flip should be allowed or not
; PROGRAMMER: Manolis K. Georgoulis (JHU/APL, 10/13/05)

res=size(Bcx) & id1=res(1) & id2=res(2)

Bxr=Bpx+Bcx & Byr=Bpy+Bcy
len1=(Bxr-Bx1)^2.+(Byr-By1)^2.+(Bzr-Bz1)^2.
len2=(Bxr-Bx2)^2.+(Byr-By2)^2.+(Bzr-Bz2)^2.
r1=where(len1 le len2,ico) & if ico gt 0. then Bz(r1)=Bz1(r1)
r1=where(len1 gt len2,ico) & if ico gt 0. then Bz(r1)=Bz2(r1)
FIND_RESPECTIVE_NEW,Bx,Bx1,Bx2,Bz,Bz1,Bz2
FIND_RESPECTIVE_NEW,By,By1,By2,Bz,Bz1,Bz2

expf=20
tx=Bx & ty=By 
tx1=Bx1 & ty1=By1
tx2=Bx2 & ty2=By2
if keyword_set(mirror) then begin 
  EXPAND_IMAGE,tx,Bx,expf,expf,stx,sty,/mirror
  EXPAND_IMAGE,ty,By,expf,expf,sty,sty,/mirror
  EXPAND_IMAGE,tx1,Bx1,expf,expf,stx,sty,/mirror
  EXPAND_IMAGE,ty1,By1,expf,expf,stx,sty,/mirror
  EXPAND_IMAGE,tx2,Bx2,expf,expf,stx,sty,/mirror
  EXPAND_IMAGE,ty2,By2,expf,expf,stx,sty,/mirror
endif else begin 
  EXPAND_IMAGE,tx,Bx,expf,expf,stx,sty
  EXPAND_IMAGE,ty,By,expf,expf,stx,sty
  EXPAND_IMAGE,tx1,Bx1,expf,expf,stx,sty
  EXPAND_IMAGE,ty1,By1,expf,expf,stx,sty
  EXPAND_IMAGE,tx2,Bx2,expf,expf,stx,sty
  EXPAND_IMAGE,ty2,By2,expf,expf,stx,sty
endelse
MEAN_STATISTICS,Bx,4
MEAN_STATISTICS,By,4
lx1=abs(Bx-Bx1) & lx2=abs(Bx-Bx2) & dx=float((lx1-lx2)/sqrt(Bx^2.+By^2.))
ly1=abs(By-By1) & ly2=abs(By-By2) & dy=float((ly1-ly2)/sqrt(Bx^2.+By^2.))

r1=where(dx*dy gt 0.,ico)
if ico gt 0. then begin 
  r2=where(dx(r1) gt 0.,ic1) 
  if ic1 gt 0. then begin & Bx(r1(r2))=Bx2(r1(r2)) & By(r1(r2))=By2(r1(r2)) &endif
  r2=where(dx(r1) lt 0.,ic1)
  if ic1 gt 0. then begin & Bx(r1(r2))=Bx1(r1(r2)) & By(r1(r2))=By1(r1(r2)) &endif
endif
r1=where(dx*dy le 0.,ico) 
if ico gt 0. then begin 
  r2=where(abs(dx(r1)) ge abs(dy(r1)),ic1)
  if ic1 gt 0. then begin 
    r3=where(dx(r1(r2)) ge 0.,ic2) 
    if ic2 gt 0. then begin & Bx(r1(r2(r3)))=Bx2(r1(r2(r3))) & By(r1(r2(r3)))=By2(r1(r2(r3))) &endif
    r3=where(dx(r1(r2)) lt 0.,ic2) 
    if ic2 gt 0. then begin & Bx(r1(r2(r3)))=Bx1(r1(r2(r3))) & By(r1(r2(r3)))=By1(r1(r2(r3))) &endif    
  endif
  r2=where(abs(dx(r1)) lt abs(dy(r1)),ic1)
  if ic1 gt 0. then begin 
    r3=where(dy(r1(r2)) gt 0.,ic2) 
    if ic2 gt 0. then begin & Bx(r1(r2(r3)))=Bx2(r1(r2(r3))) & By(r1(r2(r3)))=By2(r1(r2(r3))) &endif
    r3=where(dy(r1(r2)) lt 0.,ic2)
    if ic2 gt 0. then begin & Bx(r1(r2(r3)))=Bx1(r1(r2(r3))) & By(r1(r2(r3)))=By1(r1(r2(r3))) &endif
  endif
endif
Bx=Bx(stx:stx+id1-1,sty:sty+id2-1) 
By=By(stx:stx+id1-1,sty:sty+id2-1)
Bx1=tx1 & By1=ty1
Bx2=tx2 & By2=ty2

FIND_RESPECTIVE_NEW,Bz,Bz1,Bz2,Bx,Bx1,Bx2

RETURN
END





