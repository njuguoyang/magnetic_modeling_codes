PRO FIND_PROXY,Bx,By,lamda,Jz,Bcx,Bcy,mirror=mirror

; PURPOSE: Calculate the vertical current density Jz and the
; respective nonpotential magnetic field (Bcx,Bcy) for a horizontal
; magnetic field (Bx,By)
; PROGRAMMER: Manolis K. Georgoulis (JHU/APL, 10/13/05)

res=size(Bx) & id1=res(1) & id2=res(2)

tx=Bx & ty=By
if keyword_set(mirror) then begin 
  EXPAND_IMAGE,tx,Bx,10,10,stx,sty,/mirror
  EXPAND_IMAGE,ty,By,10,10,stx,sty,/mirror
endif else begin 
  EXPAND_IMAGE,tx,Bx,10,10,stx,sty
  EXPAND_IMAGE,ty,By,10,10,stx,sty
endelse
CALCULATE_VERTICAL_CURRENT,Bx,By,lamda,Jz
mJz=max(abs(Jz))
r1=where(abs(Jz) gt 0.4*mJz,ico) & if ico gt 0. then Jz(r1)=0.

Jz=(4.*!dpi/3.d10)*Jz(stx:stx+id1-1,sty:sty+id2-1)
if keyword_set(mirror) eq 0. then begin 
  Jz(*,0)=0. & Jz(*,id2-1)=0. & Jz(0,*)=0. & Jz(id1-1,*)=0.
endif
 
if keyword_set(mirror) then DIV_FREE_SOLUTION_RING,Jz,lamda,Bnp,/mirror else $ 
                            DIV_FREE_SOLUTION_RING,Jz,lamda,Bnp 
Bcx=dblarr(id1,id2) & Bcy=dblarr(id1,id2)
Bcx(*,*)=Bnp(*,*,0) & Bcy(*,*)=Bnp(*,*,1)

Jz=Jz*(3.d10/(4.*!dpi))
Bx=tx & By=ty

RETURN
END
