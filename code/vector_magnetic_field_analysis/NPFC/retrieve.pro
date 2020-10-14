PRO RETRIEVE,arg,crit,thres

; PURPOSE:
; Smoothes the artificially large values of arg when crit is very
; small. arg is a function of crit which goes to infinity when crit
; goes to sero
; PROGRAMMER: Manolis K. Georgoulis (JHU/APL, 10/13/05)

res=size(arg) & id1=res(1) & id2=res(2)
r1=where(crit lt thres,ico)
if ico gt 0. then begin 
  arg(r1)=0.
  for i=0,19 do $ 
  arg(r1)=(1./4.)*(arg(r1+1)+arg(r1-1)+arg(r1+id1)+arg(r1-id1))
endif
  
RETURN
END

