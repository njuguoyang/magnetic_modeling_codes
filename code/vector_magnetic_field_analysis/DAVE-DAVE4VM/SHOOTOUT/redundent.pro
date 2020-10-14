pro redundent,no_elliptic2d=no_elliptic2d,pazo=pazo
;
;     This code was written in response to Pascal Demonlin's questions
;     about the original manuscript, do the vertical velocities v_z in
;     DAVE4VM contain 'redundent' information
;
;     makes figures 7a,8a & b, 13, and 14a & b from Schuck 2008
;
if (N_ELEMENTS(PAZO) eq 0) then PAZO=0
;
   DIR='SCHUCK/'                                        
;
   scl=1.d5/1.d8
   set_plot,'x'
   mag=shootout('./old_shootout/data_for_analysis.sav',0,truth=truth,no_elliptic2d=no_elliptic2d)
;
   angle=20*!dpi/180.d0                     ; angle between LOS and z-axis
   angle=0.d0

   @optimum_windows
   WINDOW_SIZE1=DAVE4VM_WINDOW_SIZE
   WINDOW_SIZE2=DAVE_WINDOW_SIZE
;   
   mag=compute_los_velocity(mag,truth,angle)
;
   INDEX=where(abs(mag.BZ) gt mag.thr,complement=bad,ND)
   print,ND,mag.thr,format='("ND=",I4," data with |Bz|>",F6.2," Gauss")'
   plot,mag.bzt[index],(truth.ftvXX+truth.ftvYY)[index],psym=2
   print,'ANMHD=',correlate(mag.bzt[index],(truth.ftvXX+truth.ftvYY)[index])
;
;     set up comparison region
   xr=[95,205]
   yr=[90,200]
   mask=fix(mag.bz)
   mask[*]=0
   mask[xr[0]:xr[1],yr[0]:yr[1]]=1
;
   tm0=systime(1,/seconds)
   dave4vm=dave4vm(mag,window_size1,SV=SV,/double)
   tm1=systime(1,/seconds)
   dtm1=tm1-tm0
   dave=dave(mag,window_size2,/double)
   tm2=systime(1,/seconds)
   dtm2=tm2-tm1
   print,'*** Time to execute DAVE4VM: ',dtm1,' WINDOW_SIZE=',WINDOW_SIZE1
   print,'*** Time to execute    DAVE: ',dtm2,' WINDOW_SIZE=',WINDOW_SIZE2
;
   dave4vm=compute_fluxes(mag,dave4vm)
;
   sz=size(dave.U0)
   dot=(dave.U0*mag.Bx+dave.V0*mag.By)
   dave=add_tag(dave,dblarr(sz[1],sz[2]),'VPX')
   dave=add_tag(dave,dblarr(sz[1],sz[2]),'VPY')
   dave=add_tag(dave,dblarr(sz[1],sz[2]),'VPZ')
   B2=mag.bx^2+mag.by^2+mag.bz^2
   dave.VPX=(dave.U0-dot*mag.bx/B2)
   dave.VPY=(dave.V0-dot*mag.by/B2)
   dave.VPZ=-(dot*mag.bz/B2)
;
   results=compare_velocities(mag,truth,dave4vm,mask=mask,good=good,welsch=welsch)
;
   A = FINDGEN(33) * (!dPI*2/32.)  
;    Define the symbol to be a unit circle with 32 points,   
;    and set the filled flag:  
   USERSYM, COS(A), SIN(A), /FILL  
;     setup
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
   perp='!9^!3'
;
   set_plot,'ps'
;
   if (PAZO eq 1) then begin
      DEVICE, set_font='PazoMath-BoldItalic', FONT_INDEX=20 
      DEVICE, set_font='PazoMath-Italic', FONT_INDEX=10
   endif else DEVICE, set_font='Symbol', FONT_INDEX=20 
;
   DEVICE, /helvetica,/oblique, FONT_INDEX=5 
   DEVICE, /helvetica,/bold,/oblique, FONT_INDEX=6 
;
   loadct,39
   dummy=fsc_color(COLORSTRUCTURE=psc,/allcolors)
;
   dd=0.02                      ; get this number from make_plots.pro
   spearman='Spearman:'
   
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;  make plots of vz*Bx versus vx*bz
;
   xr=[-.5,2]-1
   yr=[-.8,0.8]+.3
   rat=(yr[1]-yr[0])/(xr[1]-xr[0])
   XM=!X.margin
   !X.margin[1]=0
;
   device,filename=DIR+'f8a.eps',/color,bits_per_pixel=8,xsize=5,ysize=5,$
          /inches,/encapsulate
   pos=aspect(rat,margin=0.05)
   plot,(mag.Bz[good]*dave4vm.VPX[good])*scl,(mag.Bx[good]*dave4vm.VPZ[good])*scl,$
        _extra=extra,color=psc.black,position=pos,$
        xrange=xr,yrange=yr,xstyle=1,ystyle=1,$
        xtitle='v!D'+perp+'x!N!5B!3!Dz!N (red)  !5V!3!D'+perp+'x!N!5B!3!Dz!N (blue) [10!U8!N G cm/s]',$
        ytitle='v!D'+perp+'z!N!5B!3!Dx!N (red)  !5V!3!D'+perp+'z!N!5B!3!Dx!N (blue) [10!U8!N G cm/s]',/nodata
;
   oplot,(mag.Bz[good]*dave4vm.VPX[good])*scl,(mag.Bx[good]*dave4vm.VPZ[good])*scl,$
         color=psc.red,_extra=extra
   oplot,(mag.Bz[good]*truth.VPX[good])*scl,(mag.Bx[good]*truth.VPZ[good])*scl,$
         color=psc.blue,_extra=extra
   oplot,[-2,2],[-2,2],linestyle=2,thick=extra.thick,color=psc.black
;
   U=R_Correlate(mag.Bz[good]*dave4vm.VPX[good],mag.Bx[good]*dave4vm.VPZ[good])
   V=R_Correlate(mag.Bz[good]*truth.VPX[good],mag.Bx[good]*truth.VPZ[good])
   dy=0.2
   xyouts,0.35,.925-dy,string(U[0],format='("!9r!3=",F5.2)'),color=psc.red,$
                         align=0.0,/normal,_extra=extra
   xyouts,0.35,.875-dy,string(V[0],format='("!9r!3=",F5.2)'),color=psc.blue,$
                         align=0.0,/normal,_extra=extra
   xyouts,0.3,0.9-dy,spearman,align=1,color=psc.black,_extra=extra,/normal
;
   U=Correlate(mag.Bz[good]*dave4vm.VPX[good],mag.Bx[good]*dave4vm.VPZ[good])
   V=Correlate(mag.Bz[good]*truth.VPX[good],mag.Bx[good]*truth.VPZ[good])
   xyouts,0.35,.5-dy,string(U[0],format='("!3C!N=",F5.2)'),color=psc.red,$
                         align=0.0,/normal,_extra=extra
   xyouts,0.35,.45-dy,string(V[0],format='("!3C!N=",F5.2)'),color=psc.blue,$
                         align=0.0,/normal,_extra=extra
   xyouts,0.3,0.475-dy,'Pearson:',align=1,color=psc.black,_extra=extra,/normal

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   device,filename=DIR+'f8b.eps',/color,bits_per_pixel=8,/encapsulate
   xr=[-.5,2]
   yr=[-.8,0.8]
   rat=(yr[1]-yr[0])/(xr[1]-xr[0])
;
   plot,(mag.Bz[good]*dave4vm.VPY[good])*scl,(mag.By[good]*dave4vm.VPZ[good])*scl,$
        color=psc.black,_extra=extra,position=aspect(rat,margin=0.05),$
        xrange=xr,yrange=yr,xstyle=1,ystyle=1,$
        xtitle='v!D'+perp+'y!N!5B!3!Dz!N (red) !5V!3!D'+perp+'y!N!5B!3!Dz!N (blue) [10!U8!N G cm/s]',ytitle='v!D'+perp+'z!N!5B!3!Dy!N (red) v!D'+perp+'z!N!5B!3!Dy!N (blue) [10!U8!N G cm/s]',/nodata
;
   oplot,(mag.Bz[good]*dave4vm.VPY[good])*scl,(mag.By[good]*dave4vm.VPZ[good])*scl,$
         color=psc.red,_extra=extra
   oplot,(mag.Bz[good]*truth.VPY[good])*scl,(mag.By[good]*truth.VPZ[good])*scl,$
         color=psc.blue,_extra=extra
   oplot,[-2,2],[-2,2],linestyle=2,thick=extra.thick,color=psc.black

   U=R_Correlate(mag.Bz[good]*dave4vm.VPY[good],mag.BY[good]*dave4vm.VPZ[good])
   V=R_Correlate(mag.Bz[good]*truth.VPY[good],mag.BY[good]*truth.VPZ[good])
   xyouts,0.7,.925-dy,string(U[0],format='("!9r!3=",F5.2)'),color=psc.red,$
                         align=0.0,/normal,_extra=extra
   xyouts,0.7,.875-dy,string(V[0],format='("!9r!3=",F5.2)'),color=psc.blue,$
                         align=0.0,/normal,_extra=extra
   xyouts,0.65,0.9-dy,spearman,align=1,color=psc.black,_extra=extra,/normal
;
   U=Correlate(mag.Bz[good]*dave4vm.VPY[good],mag.BY[good]*dave4vm.VPZ[good])
   V=Correlate(mag.Bz[good]*truth.VPY[good],mag.BY[good]*truth.VPZ[good])
   xyouts,0.7,.5-dy,string(U[0],format='("!3C!N=",F5.2)'),color=psc.red,$
                         align=0.0,/normal,_extra=extra
   xyouts,0.7,.45-dy,string(V[0],format='("!3C!N=",F5.2)'),color=psc.blue,$
                         align=0.0,/normal,_extra=extra
   xyouts,0.65,0.475-dy,'Pearson:',align=1,color=psc.black,_extra=extra,/normal

   device,/close_file
   !X.margin=XM
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
   d4vmx=$
    [(R_correlate((mag.Bz*dave4vm.VPX)[good],(mag.Bx*dave4vm.VPZ)[good]))[0],$
     (correlate((mag.Bz*dave4vm.VPX)[good],(mag.Bx*dave4vm.VPZ)[good]))[0]]
   trthx=$
     [(R_correlate((mag.Bz*truth.VPX)[good],(mag.Bx*truth.VPZ)[good]))[0],$
     (correlate((mag.Bz*truth.VPX)[good],(mag.Bx*truth.VPZ)[good]))[0]]
;
   d4vmy=[$
   (R_correlate((mag.Bz*dave4vm.VPY)[good],(mag.By*dave4vm.VPZ)[good]))[0],$
     (correlate((mag.Bz*dave4vm.VPY)[good],(mag.By*dave4vm.VPZ)[good]))[0]]
   trthy=[$
   (R_correlate((mag.Bz*truth.VPY)[good],(mag.By*truth.VPZ)[good]))[0],$
     (correlate((mag.Bz*truth.VPY)[good],(mag.By*truth.VPZ)[good]))[0]]
;
   openw,out,DIR+'redundent.tex',/get_lun
   printf,out,'%auto-ignore'
   printf,out,'$v_{\perp x}\,B_z$&$v_{\perp z}\,B_x$&'+string(d4vmx,format='(F5.2,"&",F5.2)')+'&'+string(trthx,format='(F5.2,"&",F5.2,"\\")')

   printf,out,'$v_{\perp y}\,B_z$&$v_{\perp z}\,B_y$&'+string(d4vmy,format='(F5.2,"&",F5.2)')+'&'+string(trthy,format='(F5.2,"&",F5.2,"\\")')

   d4vmx=[(R_correlate((dave4vm.VPX)[good],(dave4vm.VPZ)[good]))[0],$
     (correlate((dave4vm.VPX)[good],(dave4vm.VPZ)[good]))[0]]
   trthx=[(R_correlate((truth.VPX)[good],(truth.VPZ)[good]))[0],$
     (correlate((truth.VPX)[good],(truth.VPZ)[good]))[0]]
;
   d4vmy=[(R_correlate((dave4vm.VPY)[good],(dave4vm.VPZ)[good]))[0],$
     (correlate((dave4vm.VPY)[good],(dave4vm.VPZ)[good]))[0]]
   trthy=[(R_correlate((truth.VPY)[good],(truth.VPZ)[good]))[0],$
     (correlate((truth.VPY)[good],(truth.VPZ)[good]))[0]]
;
   printf,out,'$v_{\perp x}$&$v_{\perp z}$&'+string(d4vmx,format='(F5.2,"&",F5.2)')+'&'+string(trthx,format='(F5.2,"&",F5.2,"\\")')

   printf,out,'$v_{\perp y}$&$v_{\perp z}$&'+string(d4vmy,format='(F5.2,"&",F5.2)')+'&'+string(trthy,format='(F5.2,"&",F5.2,"\\")')
   close,out
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;  Plot DAVE velocities versus ANMHD perpendicular horizontal plasma velocities
;
   device,filename=DIR+'f21a.eps',/color,bits_per_pixel=8,xsize=5,ysize=5,/inches,/encapsulate
   pos=aspect(1.0,margin=0.08)+.06
;
   plot,[(mag.Bz*truth.VPX)[good],(mag.Bz*truth.VPY)[good]]*scl,[(mag.Bz*dave.U0)[good],(mag.Bz*dave.V0)[good]]*scl,psym=8,color=psc.black,position=pos,$
   xtitle=  'ANMHD  !5V!3!D'+perp+'x!N!5B!3!Dz!N (red) & !5V!3!D'+perp+'y!N!5B!3!Dz!N (blue) [10!U8!N G cm/s]',$
   ytitle='DAVE  !10J!3!Dx!N!5B!3!Dz!N (red) & !10J!3!Dy!N!5B!3!Dz!N (blue) [10!U8!N G cm/s]',$
          _extra=extra,/nodata
   oplot,(mag.Bz*truth.VPX)[good]*scl,(mag.Bz*dave.U0)[good]*scl,psym=8,color=psc.red,$
        _extra=extra
   oplot,(mag.Bz*truth.VPY)[good]*scl,(mag.Bz*dave.V0)[good]*scl,psym=8,color=psc.blue,$
        _extra=extra
   oplot,[-1.5,1.5],[-1.5,1.5],linestyle=2,thick=extra.thick,color=psc.black
;
   ux=[correlate((mag.Bz*truth.VPX)[good],(mag.Bz*dave.U0)[good]),(r_correlate((mag.Bz*truth.VPX)[good],(mag.Bz*dave.U0)[good]))[0]]
   uy=[correlate((mag.Bz*truth.VPY)[good],(mag.Bz*dave.V0)[good]),(r_correlate((mag.Bz*truth.VPY)[good],(mag.Bz*dave.V0)[good]))[0]]
;
   Lx=ladfit((mag.Bz[good]*truth.VPX[good]),(mag.Bz[good]*dave.U0[good]))
   Ly=ladfit((mag.Bz[good]*truth.VPY[good]),(mag.Bz[good]*dave.V0[good]))
;
   xyouts,0.5,.85,string(ux[1],format='("!9r!3!Dx!N=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
   xyouts,0.5,.8,string(uy[1],format='("!9r!3!Dy!N=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
   xyouts,0.45,0.825,spearman,align=1,color=psc.black,_extra=extra,/normal
;
   xyouts,0.65,.30+dd,string(ux[0],format='("C!Dx!N=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
   xyouts,0.65,.25+dd,string(uy[0],format='("C!Dy!N=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
   xyouts,0.6,0.275+dd,'Pearson:',align=1,color=psc.black,_extra=extra,/normal
;
   xyouts,0.65,0.20+dd,string(Lx[1],format='("S!Dx!N=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
   xyouts,0.65,0.15+dd,string(Ly[1],format='("S!Dy!N=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
   xyouts,0.6,0.175+dd,'Slopes:',align=1,color=psc.black,_extra=extra,/normal


   print,ux,uy
;
   device,/close_file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;  Plot DAVE velocities versus ANMHD horizontal plasma velocities
;
;
   device,filename=DIR+'f13.eps',/color,bits_per_pixel=8,/encapsulate
   pos=aspect(1.0,margin=0.08)+.06
   scl=1.d5/1.d8
   plot,[(mag.Bz*truth.VX)[good],(mag.Bz*truth.VY)[good]]*scl,[(mag.Bz*dave.U0)[good],(mag.Bz*dave.V0)[good]]*scl,psym=8,color=psc.black,position=pos,$
   xtitle=  'ANMHD  !5V!3!Dx!N!5B!3!Dz!N (red) & !5V!3!Dy!N!5B!3!Dz!N (blue) [10!U8!N G cm/s]',$
   ytitle='DAVE  !10J!3!Dx!N!5B!3!Dz!N (red) & !10J!3!Dy!N!5B!3!Dz!N (blue) [10!U8!N G cm/s]',$
          _extra=extra,xrange=[-1,1],yrange=[-1,1],/nodata
   oplot,(mag.Bz*truth.VX)[good]*scl,(mag.Bz*dave.U0)[good]*scl,psym=8,color=psc.red,$
        _extra=extra
   oplot,(mag.Bz*truth.VY)[good]*scl,(mag.Bz*dave.V0)[good]*scl,psym=8,color=psc.blue,$
        _extra=extra
   oplot,[-1.5,1.5],[-1.5,1.5],linestyle=2,thick=extra.thick,color=psc.black
;
   ux=[correlate((mag.Bz*truth.VX)[good],(mag.Bz*dave.U0)[good]),(r_correlate((mag.Bz*truth.VX)[good],(mag.Bz*dave.U0)[good]))[0]]
   uy=[correlate((mag.Bz*truth.VY)[good],(mag.Bz*dave.V0)[good]),(r_correlate((mag.Bz*truth.VY)[good],(mag.Bz*dave.V0)[good]))[0]]
;
   spearman='Spearman:'
   Lx=ladfit((mag.Bz*truth.VX)[good],(mag.Bz*dave.U0)[good])
   Ly=ladfit((mag.Bz*truth.VY)[good],(mag.Bz*dave.V0)[good])
;
   xyouts,0.45,.85,string(ux[1],format='("!9r!3!Dx!N=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
   xyouts,0.45,.8,string(uy[1],format='("!9r!3!Dy!N=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
   xyouts,0.4,0.825,spearman,align=1,color=psc.black,_extra=extra,/normal
;
   xyouts,0.75,.325+dd,string(ux[0],format='("C!Dx!N=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
   xyouts,0.75,.275+dd,string(uy[0],format='("C!Dy!N=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
   xyouts,0.7,0.3+dd,'Pearson:',align=1,color=psc.black,_extra=extra,/normal
;
   xyouts,0.75,0.20+dd,string(Lx[1],format='("S!Dx!N=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
   xyouts,0.75,0.15+dd,string(Ly[1],format='("S!Dy!N=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
   xyouts,0.7,0.175+dd,'Slopes:',align=1,color=psc.black,_extra=extra,/normal
;
   device,/close_file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  parallel velocity
;
   device,filename=DIR+'f7a.eps',color=0,/encapsulate
;
   R=r_correlate(truth.vb[good],dave4vm.vb[good])
   C=correlate(truth.vb[good],dave4vm.vb[good])
   S=ladfit(truth.vb[good],dave4vm.vb[good])
   plot,truth.vb[good],dave4vm.vb[good],xstyle=1,ystyle=1,$
         xrange=[-.4,.4],yrange=[-.4,.4],$
        _extra=extra,position=pos,$
                      xtitle='ANMHD  !5V!3!D||!N [km/s]',$
                     ytitle='DAVE4VM  !5v!3!D||!N [km/s]'
   oplot,[-.4,.4],[-.4,.4],thick=extra.thick,linestyle=2
   xyouts,0.5,.9,string(R[0],format='("!9r!3=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
   xyouts,0.5,.85,string(R[0],format='("C=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
   xyouts,0.5,.8,string(S[1],format='("S=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
;
   device,/close_file
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   X=truth.vb[good]
   Y=dave4vm.vb[good]
   N_TRIALS=10000L
   R=R_correlate(X,Y,probd=probd,zd=zd)
   print,'********************************'
   print,'Parallel Velocity'
   print,'Computing spearman permutation....'
   dist=spearman_permutation(X,Y,N_TRIALS)
   mom=moment(dist,/double)
   print,'Moments of the test distribution'
   print,'    MEAN=',mom[0],sqrt(mom[1])
   print,'    SKEW=',mom[2]
   print,'KURTOSIS=',mom[3]
   print
   Z=fisher_z_transform(R[0])
   print,'       R=',R[0]
   print,' Z_SCORE=',Z
   print,'',Z/sqrt(mom[1])
   print,'********************************'
;
;     compute chi2
   FPX=truth.vpx*mag.bz
   FPY=truth.vpy*mag.bz
   FX=truth.vx*mag.bz
   FY=truth.vy*mag.bz
;
   vbh=(truth.vx*mag.Bx+truth.vy*mag.By)/(mag.bx^2+mag.by^2)
   ftvpX=(truth.vx-vbh*mag.Bx)*mag.bz
   ftvpY=(truth.vy-vbh*mag.By)*mag.bz
;
   
;
   FAUX_TRANSPORT=replicate({NAME:'Horizontal Velocities',X:FX,Y:FY},4)
   FAUX_TRANSPORT[1].NAME='Perpendicular Horizontal Velocities'
   FAUX_TRANSPORT[1].X=FPX
   FAUX_TRANSPORT[1].Y=FPY
   FAUX_TRANSPORT[2].NAME='Flux Transport Velocities'
   FAUX_TRANSPORT[2].X=truth.ftvx
   FAUX_TRANSPORT[2].Y=truth.ftvy
;
   FAUX_TRANSPORT[3].NAME='Flux Transport Velocities'
   FAUX_TRANSPORT[3].X=ftvpx
   FAUX_TRANSPORT[3].Y=ftvpy
;
   DT=replicate(FAUX_TRANSPORT[0],3)
;
   DT[0].NAME='RAW DAVE'
   DT[0].X=dave.U0*mag.bz
   DT[0].Y=dave.V0*mag.bz
;
   DT[1].NAME='DAVE PERPENDICULARIZED'
   DT[1].X=dave.vpx*mag.bz
   DT[1].Y=dave.vpy*mag.bz
;
   DT[2].NAME='SCHUCK HYPOTHESIS'

   DT[2].X=dave.U0*mag.bz-truth.vz*mag.bx
   DT[2].Y=dave.V0*mag.bz-truth.vz*mag.by
;
   set_plot,'ps'
   cdot='!9'+string(215B)+'!3'
;
   print,"Pascal's HYPOTHESIS: U from DAVE is the footpoint velocity"
   filename=DIR+'f14a.eps'
   xtitle=['ANMHD','-!6B!3!dh!N'+cdot+'!6UB!3!Dz!N/(4!9p!3)','-2!6A!3!dp!N'+cdot+'!6U!5B!3!Dz!N']
   ytitle=['DAVE','-!6B!3!dh!N'+cdot+'!20J!5B!3!Dz!N/(4!9p!3)','-2!6A!3!dp!N'+cdot+'!20J!5B!3!Dz!N']
   output_fluxes,FAUX_TRANSPORT[2],DT[0],truth,mag,xtitle,ytitle,filename,good=good,pazo=pazo
;
   print & print,"Chae's/Pete's HYPOTHESIS: U from DAVE is V_H"
   filename=DIR+'f14b.eps'
   xtitle=['ANMHD','-!6B!3!dh!N'+cdot+'!6U!5B!3!Dz!N/(4!9p!3)','-2!6A!3!dp!N'+cdot+'!6U!5B!3!Dz!N']
   ytitle=['DAVE+ANMHD','-!6B!3!dh!N'+cdot+'(!20J!5B!3!Dz!N-!5V!3!dz!n!6B!3!Dh!N)/(4!9p!3)','-2!6A!3!dp!N'+cdot+'(!20J!5B!3!Dz!N-!5V!3!dz!n!6B!3!Dh!N)']
   output_fluxes,FAUX_TRANSPORT[2],DT[2],truth,mag,xtitle,ytitle,filename,good=good,pazo=pazo

end
