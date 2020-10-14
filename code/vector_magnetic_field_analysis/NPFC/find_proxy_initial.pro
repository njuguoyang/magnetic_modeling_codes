PRO FIND_PROXY_INITIAL,Lc,Btr,phi,Bx,By,lamda,Jz,Bcx,Bcy,mirror=mirror

; PURPOSE: Calculates an initial proxy vertical current density and
; the respective nonpotential field for the calculations. This proxy
; current is given by 
;          Jz=s(Jzl)*[sin(Lc)^2.*|Jzl| + cos(Lc)^2.*Jzss]
; where Jzl is the vertical current density obtained from the mean
; heliographic magnetic field vector (Georgoulis 2005), s(Jzl) is the
; sign of Jzl, Jzss is the absolute value of the minimum current
; derived by Semel & Skumanich (1998), and Lc is the heliographic
; longitude. 
; PROGRAMMER: Manolis K. Georgoulis (JHU/APL, 10/13/2005)

res=size(Bx) & id1=res(1) & id2=res(2)

tx=Bx & ty=By & ttr=Btr & tphi=phi & tLc=Lc
if n_elements(Lc) eq 1 then begin 
  Lc=fltarr(2*id1,2*id2) & Lc(*,*)=tLc
endif else EXPAND_IMAGE,tLc,Lc,id1,id2,stx,sty,/mirror
if keyword_set(mirror) then begin   
  EXPAND_IMAGE,tx,Bx,id1,id2,stx,sty,/mirror
  EXPAND_IMAGE,ty,By,id1,id2,stx,sty,/mirror
  EXPAND_IMAGE,ttr,Btr,id1,id2,stx,sty,/mirror
  EXPAND_IMAGE,tphi,phi,id1,id2,stx,sty,/mirror
endif else begin 
  EXPAND_IMAGE,tx,Bx,id1,id2,stx,sty
  EXPAND_IMAGE,ty,By,id1,id2,stx,sty
  EXPAND_IMAGE,ttr,Btr,id1,id2,stx,sty
  EXPAND_IMAGE,tphi,phi,id1,id2,stx,sty
endelse
CALCULATE_VERTICAL_CURRENT,Bx,By,lamda,Jz
CALCULATE_SS_CURRENT,Btr,phi,lamda,Jzm & Jzm=Jzm*(3.d10/(4.*!dpi))
r1=where(Jz ne 0.) 
Jz(r1)=Jz(r1)/abs(Jz(r1))*(sin(Lc(r1))^2.*abs(Jz(r1)) + cos(Lc(r1))^2.*Jzm(r1))
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
Bx=tx & By=ty & Btr=ttr & phi=tphi & Lc=tLc

RETURN
END
