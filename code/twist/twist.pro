;+
; NAME :
;   twist
; PURPOSE :
;   Calculate the twist between one field line and the axis (see Berger, A.
;   and Prior, C. 2006, J. Phys. A: Math. Gen. 39, 8321).
; CATEGORY :
;
; CALLING SEQUENCE :
;   twist,axis_filename,line_filename
; INPUTS :
;   The names of the files that contain the Cartesian coordinates of the 
;   axis, and the Cartesian coordinates of another field line which twists 
;   around the axis.
; OUTPUTS :
;   The twist number in the unit of radian or turns (1 turn= 2 pi).
; COMMON BLOCKS :
;   None
; MODIFICATION HISTORY :
;   2010.03 Guo Yang @ Nanjing University, 
;
;   This software is provided without any warranty. Permission to use,
;   copy, modify. Distributing modified or unmodified copies is granted,
;   provided this disclaimer and information are included unchanged.
;-

function twist,axis_filename,line_filename,rad=rad

If(n_params() Lt 2) Then Begin
  print,'================================='
  print,'ERROR!'
  print,'CALLING SEQUENCE :'
  print,'twist,axis_filename,line_filename'
  print,'================================='
  Return,0
Endif

If((size(axis_filename,/type) ne 7) OR (size(line_filename,/type) ne 7)) Then Begin
  print,'============================================================'
  print,'ERROR!'
  print,'Please specify the file name of the axis and the field line!'
  print,'============================================================'
  Return,0
Endif

restore,axis_filename    ;,/verbose      ;restore coordinates of the axis curve
axis_tangentx=tangentx
axis_tangenty=tangenty
axis_tangentz=tangentz
axis_linex=linex
axis_liney=liney
axis_linez=linez
restore,line_filename    ;,/verbose      ;restore coordinates of the secondary curve

epsvdb=0.01    ;Errors allowed for VdotB, which is a the inner product of the vector, pointing from the axis curve to the secondary curve, 
               ;and the field direction B. If VdotB=0, the vector V maps a point on the axis curve to the secondary curve, where the 
               ;infinitesimal twist is defined at this point on the secondary curve with respect to the point on the axis curve.
;window,1
;window,2
ss=size(axis_linex,/n_elements)
print,'Points of the axis curve:',ss
index0=lonarr(ss)
vdbmin=dblarr(ss)
ssl=size(linex,/n_elements)
print,'Points of the secondary curve:',ssl
for i=0,ss-1 do begin
  veca2s=[[linex-axis_linex[i]],[liney-axis_liney[i]],[linez-axis_linez[i]]]     ;the vector points from the i_th point on the axis curve to the secondary curve
;  print,'================='
;  print,'Step:',i
;  help,veca2s
  vnorm=SQRT(veca2s[*,0]*veca2s[*,0]+veca2s[*,1]*veca2s[*,1]+veca2s[*,2]*veca2s[*,2])
  vectn=[[veca2s[*,0]/vnorm],[veca2s[*,1]/vnorm],[veca2s[*,2]/vnorm]]
;  help,vectn
  VdotB=axis_tangentx[i]*vectn[*,0]+axis_tangenty[i]*vectn[*,1]+axis_tangentz[i]*vectn[*,2]
;  wset,1
;  plot,VdotB,xstyle=1+2,ystyle=1,yrange=[-1,2],title='V dot B distribution along the secondary curve'
;  oplot,vectn[*,0]*vectn[*,0]+vectn[*,1]*vectn[*,1]+vectn[*,2]*vectn[*,2]
  vdb=0.8*(shift(vdotB,-1)-shift(vdotB,1))
  vdb[0]=0.0 & vdb[ssl-1]=0.0   ; Set the first and last points to 0.0
  index1 = where(abs(vdb) le 0.5*max(abs(vdb)))    ; It stores the index where VdotB changes slow. We expect that VdotB changes the fastest close to VdotB=0.0
  VdotB[index1] = max(abs(VdotB))
  index0[i]=where(abs(VdotB) eq min(abs(VdotB)))   ;it contains the index along the secondary curve, where the vector V poiting from the axis to this point has a minmum VdotB.
  VdBmin[i]=VdotB[index0[i]]
;  print,'VdBmin:',VdBmin[i]
;  oplot,[index0[i],index0[i]],[-1,2]
;  wset,2
;  plot,vdb,xstyle=1+2,title='Derivative of V dot B distribution along the secondary curve'
;  oplot,[index0[i],index0[i]],[-1,2]
;  wait,0.05
endfor

;window,3
;plot,index0,xstyle=1+2,title='Index of the secondary curve to map a ribbon with the axis curve'
;window,4
;plot,VdBmin,xstyle=1+2,title='min V dot B distribution along the axis curve',charsize=2.0
;print,'max, min of VdBmin:',max(vdbmin),min(vdbmin)
;for i=0,ss-1 do begin
;  if (abs(VdBmin[i]) le epsvdb) then begin
;    istart=i
;    break
;  endif
;endfor
;for i=istart+1,ss-1 do begin
;  if ((abs(VdBmin[i]) gt epsvdb) or (index0[i] lt index0[i-1]) or (i eq ss-1)) then begin
;    iend=i-1
;    break
;  endif
;endfor
;print,'istart, iend along the axis:',istart,iend
;;oplot,[istart,istart],[-1,1]
;;oplot,[iend,iend],[-1,1]
;;xyouts,0.23,0.85,'Max of VdBmin between istart and iend:  '+strtrim(string(max(vdbmin[istart:iend]),format='(f6.3)'),2),/normal,charsize=2.0

index2 = where(abs(VdBmin) le epsvdb)
vect1=[[linex[index0[index2]]-axis_linex[index2]],[liney[index0[index2]]-axis_liney[index2]],[linez[index0[index2]]-axis_linez[index2]]]
vnorm=SQRT(vect1[*,0]*vect1[*,0]+vect1[*,1]*vect1[*,1]+vect1[*,2]*vect1[*,2])
vectn=[[vect1[*,0]/vnorm],[vect1[*,1]/vnorm],[vect1[*,2]/vnorm]]
;help,vectn
tangn=[[axis_tangentx[index2]],[axis_tangenty[index2]],[axis_tangentz[index2]]]
;help,tangn
ssv=size(vectn[*,0],/n_elements)
if (ssv ge 3) then begin
  Deltv=dblarr(ssv,3)
  Deltv[*,0]=deriv(vectn[*,0])
  Deltv[*,1]=deriv(vectn[*,1])
  Deltv[*,2]=deriv(vectn[*,2])
  ;Deltwist=tangn[*,0]*(vectn[*,1]*Deltv[*,2]-vectn[*,2]*Deltv[*,1])+tangn[*,1]*(vectn[*,2]*Deltv[*,0]-vectn[*,0]*Deltv[*,2])+tangn[*,2]*(vectn[*,0]*Deltv[*,1]-vectn[*,1]*Deltv[*,0])
  Deltwist=(tangn[*,1]*vectn[*,2]-tangn[*,2]*vectn[*,1])*Deltv[*,0]+(tangn[*,2]*vectn[*,0]-tangn[*,0]*vectn[*,2])*Deltv[*,1]+(tangn[*,0]*vectn[*,1]-tangn[*,1]*vectn[*,0])*Deltv[*,2]
;  window,5
;  plot,Deltwist,xstyle=1+2,title='Twist density distribution along the axis curve',charsize=2.0
;  print,'Mean and Max deltwist:',mean(Deltwist),max(Deltwist)
  ;window,6
  ;plot,Deltwist[100:200],xstyle=1+2
  index3=where(Deltwist lt 0 and abs(Deltwist) le 10.0*abs(mean(Deltwist)))   ; We use these consditions to exclude spurious data points, which might be caused by the highly curved 
                                                                              ; axis curve. The vector V could point to a further second place where V is still perpendicular to B. 
                                                                              ; Such a case induces very large fake twsit densities or twist with an opposite sign. But the parameter 
                                                                              ; 10 might not be valid for case with real extremely large twist density.
  n3 = n_elements(index3)
  index4=where(Deltwist gt 0 and abs(Deltwist) le 10.0*abs(mean(Deltwist)))
  n4 = n_elements(index4)
  if (n3 ge n4) then begin
    Totwist=total(Deltwist[index3])
    ;window,7
    ;plot,Deltwist[index3],xstyle=1+2,title='Negative twist density distribution along the axis curve',charsize=2.0
    ;print,'Only negative twist is counted!'
    print,'Points with negative twist along the axis:',n3
  endif else begin
    Totwist=total(Deltwist[index4]) 
    ;window,7
    ;plot,Deltwist[index4],xstyle=1+2,title='Positive twist density distribution along the axis curve',charsize=2.0
    ;print,'Only positive twist is counted!'
    print,'Points with positive twist along the axis:',n4    
  endelse 
  print,'Total twist (Radian):',Totwist
  print,'Total twist (Turn):',Totwist/2.0/!pi
  ;wdelete,1
  ;wdelete,2
  if (keyword_set (rad)) then begin
    return,Totwist
  endif else begin
    return,Totwist/2.0/!pi
  endelse
endif else begin
  print,'Total twist (Radian):',0.0
  print,'Total twist (Turn):',0.0
  return,0.0
endelse

end
