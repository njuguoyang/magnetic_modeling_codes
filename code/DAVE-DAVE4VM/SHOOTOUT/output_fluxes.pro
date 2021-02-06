pro output_fluxes,adjusted,dave,truth,mag,xtitle,ytitle,filename,good=good,pazo=pazo
;
    if (N_ELEMENTS(PAZO) eq 0) then PAZO=0
;
    scl0=1.d5/1.d10
    scl1=1.d10/1.d19
;
    FTX=mag.Bz*truth.vx-mag.Bx*truth.vz
    FTY=mag.Bz*truth.vy-mag.By*truth.vz
    HT=-2*(FTX*mag.AX+FTY*mag.AY)*scl1
    PT=-(FTX*mag.BX+FTY*mag.BY)/(4*!dpi)*scl0   
;
    HF=-2*(adjusted.X*mag.AX+adjusted.Y*mag.AY)*scl1
    PF=-(adjusted.X*mag.BX+adjusted.Y*mag.BY)/(4*!dpi)*scl0   
;
    DHF=-2*(dave.X*mag.AX+dave.Y*mag.AY)*scl1
    DPF=-(dave.X*mag.BX+dave.Y*mag.BY)/(4*!dpi)*scl0   
;
    R_HF=R_CORRELATE(HF[good],DHF[good])
    C_HF=CORRELATE(HF[good],DHF[good])
    L_HF=ladfit(HF[good],DHF[good])
    THF=total(HF[good],/double)
    RAT_HF=total(DHF[good],/double)/THF
;
    R_PF=R_CORRELATE(PF[good],DPF[good])
    C_PF=CORRELATE(PF[good],DPF[good])
    L_PF=ladfit(PF[good],DPF[good])
    TPF=total(PF[good],/double)
    RAT_PF=total(DPF[good],/double)/TPF
;
    print,'helicity flux'
    print,'  R=',R_HF
    print,'  C=',C_HF
    print,'  L=',L_HF
    print,'RAT=',RAT_HF
    print,'PIECES=',moment(DHF[good],/double)
    print,'PIECES=',moment(HF[good],/double)
    print,'estimate total=',total(DHF[good],/double)
    print,'truth total=',THF
    print,'total helicity=',total(HT[good],/double)
;
    print,'Poynting flux'
    print,'  R=',R_PF
    print,'  C=',C_PF
    print,'  L=',L_PF
    print,'RAT=',RAT_PF
    print,'PIECES=',moment(DPF[good],/double)
    print,'PIECES=',moment(PF[good],/double)
    print,'estimate total=',total(DPF[good],/double)
    print,'truth total=',TPF
    print,'total Poynting=',total(PT[good],/double)
;
;     Make a vector of 33 points, A[i] = 2pi/32:  
    A = FINDGEN(33) * (!dPI*2/32.)  
;     Define the symbol to be a unit circle with 33 points,   
;     and set the filled flag:  
    USERSYM, COS(A), SIN(A), /FILL  
;      setup
    th=4
    cth=2
    cs=1.25
    symsize=.1
    extra={ticklen:-0.02,xthick:th,ythick:th,charsize:cs,thick:th,$
           symsize:symsize,psym:8,charthick:cth}
    !X.MARGIN=[10,3]
    !Y.MARGIN=[4,2]
    XM0=!X.margin
    YM0=!Y.margin
    !P.font=0

    if (PAZO eq 1) then begin
       DEVICE, set_font='PazoMath-BoldItalic', FONT_INDEX=20 
       DEVICE, set_font='PazoMath-Italic', FONT_INDEX=10
    endif else DEVICE, set_font='Symbol', FONT_INDEX=20 
    d=0.5

    loadct,39
    dummy=fsc_color(COLORSTRUCTURE=psc,/allcolors)

    device,filename=filename,/encapsulate,/color,bits_per_pixel=8
    device,xsize=5,ysize=5,/inches
    pos=aspect(1.0,margin=0.08)+.06
    label=[' [10!U19!N Mx!U2!N/cm!U2!N/s]',' [10!U10!N ergs/cm!U2!N/s]']
    labels=['','']
    plot,HF[good],DHF[good],psym=8,charsize=cs,xrange=[-10,6],yrange=[-10,6],$
         xstyle=1,ystyle=1,_extra=extra,color=psc.black,/nodata,$
         xtitle=xtitle[0]+'    (red) '+xtitle[1]+labels[1]+'    (blue) '+xtitle[2]+labels[0],$
         ytitle=ytitle[0]+'    (red) '+ytitle[1]+labels[1]+'     (blue) '+ytitle[2]+labels[0],$
         position=aspect(1.0,margin=0.05)
    oplot,HF[good],DHF[good],psym=8,color=psc.blue,_extra=extra
    oplot,PF[good],DPF[good],psym=8,color=psc.red,_extra=extra
    oplot,[-10,6],[-10,6],color=psc.black,thick=extra.thick,linestyle=2
;
    xyouts,0.35,.85,string(R_PF[0],format='("!9r=!3",F4.2)'),$
                         align=0.0,color=psc.red,/normal,_extra=extra
    xyouts,0.35,.8,string(R_HF[0],format='("!9r=!3",F4.2)'),$
                         align=0.0,color=psc.blue,/normal,_extra=extra
   xyouts,0.3,0.825,'Spearman:',align=1,color=psc.black,_extra=extra,/normal
;
    xyouts,0.35,.725,string(C_PF,format='("C!3=",F4.2)'),$
                         align=0.0,color=psc.red,/normal,_extra=extra
    xyouts,0.35,.675,string(C_HF,format='("C!3=",F4.2)'),$
                         align=0.0,color=psc.blue,/normal,_extra=extra
   xyouts,0.3,0.7,'Pearson:',align=1,color=psc.black,_extra=extra,/normal
;
   xyouts,0.7,0.325,'Slopes:',align=1,color=psc.black,_extra=extra,/normal
   xyouts,0.75,.35,string(L_PF[1],format='("S!3=",F4.2)'),$
                         align=0.0,color=psc.red,/normal,_extra=extra
   xyouts,0.75,.3,string(L_HF[1],format='("S!3=",F4.2)'),$
                         align=0.0,color=psc.blue,/normal,_extra=extra

   xyouts,0.7,0.2,'Ratio of Totals:',align=1,color=psc.black,_extra=extra,/normal
   xyouts,0.75,.225,string(Rat_PF[0],format='("R!3=",F5.2)'),$
                         align=0.0,color=psc.red,/normal,_extra=extra
   xyouts,0.75,.175,string(Rat_HF[0],format='("R!3=",F5.2)'),$
                         align=0.0,color=psc.blue,/normal,_extra=extra
;
;
    device,/close_file
;
    print & print & print
;
end
