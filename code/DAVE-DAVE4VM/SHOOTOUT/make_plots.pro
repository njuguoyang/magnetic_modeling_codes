;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro make_plots,dir,mag,truth,vel,good,tag,name,results,pazo=pazo
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;   The routine makes most of the figures for the manuscript
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (N_ELEMENTS(PAZO) eq 0) then PAZO=0
;
   eps='.eps'
   files={$
          ftv:'f2',$
          ftv_corr:'f3',$
          vp_corr:'f4',$
          ftv_theta:'f5',$
          induction_corr:'f9',$
          E_corr:'f10',$
          poynting_corr:'f11',$
          helicity_corr:'f12',$
          v_corr:'f7'$
          }
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   common constants,speed_of_light
;
   spearman='!3Spearman:'
   RO='!9r!3'
;
   dd=0.02
;
   xr=results.xr
   yr=results.yr
;
;    Make a vector of 16 points, A[i] = 2pi/16:  
   A = FINDGEN(33) * (!dPI*2/32.)  
;    Define the symbol to be a unit circle with 33 points,   
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
   d=0.5
   set_plot,'ps'
   !P.font=0

   if (PAZO eq 1) then begin
       DEVICE, set_font='PazoMath-BoldItalic', FONT_INDEX=20 
       DEVICE, set_font='PazoMath-Italic', FONT_INDEX=10
   endif else DEVICE, set_font='Symbol', FONT_INDEX=20 

   DEVICE, /helvetica,/oblique, FONT_INDEX=5 
   DEVICE, /helvetica,/bold,/oblique, FONT_INDEX=6 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;     Induction equation
   device,filename=dir+files.induction_corr+tag+eps,/encapsulate,color=0
   device,xsize=5,ysize=5,/inches

   pos=aspect(1.0,margin=0.08)+.06
;
   scl=1.d0
   cdot=string(215B)
;  !10J is vartheta
   del=string(209B)
   cross=string(180B)
   if (name eq 'DAVE') then $
     ytitle='!3'+name+'  !9'+del+'!3!Dh!N!9'+cdot+'!3(!20J!5B!3!Dz!N) [G/s]' $
   else $
     ytitle='!3'+name+'  !9'+del+'!3!Dh!N!9'+cdot+'!3(u!5B!3!Dz!N) [G/s]'
   print,':',name,':' 

   plot,mag.bzt[good],Vel.ftvxx[good]+vel.ftvyy[good],$
        _extra=extra,position=pos,$
        /nodata,color=0,xrange=[-2,2],yrange=[-2,2],$
        xstyle=1,ystyle=1,$
        xtitle='!3ANMHD  !9D!5B!3!Dz!N/!9D!5t!3 [G/s]',$
;         changes here to make plots 'in terms of electric field'
;        ytitle='!3'+name+'  [!9'+del+cross+'!3(v!9'+cross+'!3B)]!Dz!N [G/s]'
        ytitle=ytitle

   oplot,mag.bzt[good],vel.ftvxx[good]+vel.ftvyy[good],$
         _extra=extra,color=0
   oplot,[-2,2],[2,-2],linestyle=2,thick=extra.thick,color=0
;
   xyouts,0.5,.25+dd,string(results.induction_R,format='("!9r!3=",F5.2)'),$
                         align=0.0,color=0,/normal,_extra=extra
;
   xyouts,0.5,.2+dd,string(results.induction_C,format='("!3C=",F5.2)'),$
                         align=0.0,color=0,/normal,_extra=extra
   xyouts,0.5,.15+dd,string(results.induction_L[1],format='("!3S=",F5.2)'),$
                         align=0.0,color=0,/normal,_extra=extra
   xyouts,0.45,0.25+dd,spearman,align=1,color=0,_extra=extra,/normal
   xyouts,0.45,0.2+dd,'Pearson:',align=1,color=0,_extra=extra,/normal
   xyouts,0.45,0.15+dd,'Slope:',align=1,color=0,_extra=extra,/normal
;
   device,/close_file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;     flux transport velocity field
   sz=size(mag.bz)
   X=lindgen(sz[1])
   Y=lindgen(sz[2])
;     set cutoff at 2.5% of max(Bz)
;   mxbz=max(mag.Bz)
;   kill=where(abs(mag.bz) lt 0.025*mxbz)
;
   mx=0.5*max(abs(mag.Bz[xr[0]:xr[1],yr[0]:yr[1]]))
   mx=3000
   loadct,0
   black=0
;
;     color bar
   device,filename=dir+files.ftv+'c'+eps,/encapsulate,/color,bits_per_pixel=8
   device,xsize=5,ysize=1,/inches
   hextra=extra
   hextra.ticklen=-0.15
   Colorbar,_extra=hextra,font=1,range=[-mx,mx],$
            position=[0.05,0.6,.95,.9],minor=10
   xyouts,0.5,0.1,'!5B!3!Dz!N [G]',_extra=extra,align=0.5,/normal
   device,/close_file
;
;     magnetic field
   device,filename=dir+files.ftv+'d'+eps,/encapsulate,/color,bits_per_pixel=8
   device,xsize=5,ysize=5,/inches
   pos=aspect(1.0,margin=0.08)
;
   loadct,0
   plotimage,bytscl(mag.Bz[xr[0]:xr[1],yr[0]:yr[1]],min=-mx,max=mx),$
             _extra=extra,color=black,$
             xtitle='!3[Pixels]',ytitle='!3[Pixels]',position=pos+.06,$
             imgxrange=[xr[0]-d,xr[1]+d],imgyrange=[yr[0]-d,yr[1]+d]
  
;
   loadct,39
   dummy=fsc_color(COLORSTRUCTURE=psc,/allcolors)
   BZ=smooth(mag.bz,7)
   contour,Bz[xr[0]:xr[1],yr[0]:yr[1]],X[xr[0]:xr[1]],Y[yr[0]:yr[1]],$
           color=psc.blue,levels=[0],/overplot

   SKIP=3L
   IX=where(X mod skip ne 0,KX)
   IY=where(Y mod skip ne 1,KY)
;
   BX=mag.BX
   BY=mag.BY
;
   for i=0L,kx-1L do BX[IX[i],*]=!values.d_nan
   for j=0L,ky-1L do BY[*,IY[j]]=!values.d_nan
;   
   velovect,bX[xr[0]:xr[1],yr[0]:yr[1]],bY[xr[0]:xr[1],yr[0]:yr[1]],$
            X[xr[0]:xr[1]],Y[yr[0]:yr[1]],/overplot,$
            length=skip,color=psc.aquamarine,thick=3

   device,/close_file
;
;     Bh angles
   device,filename=dir+files.ftv+'e'+eps,/encapsulate,color=0
   device,xsize=5,ysize=5,/inches
;
   bh=atan(mag.By[good],mag.bx[good])
   plothist,bh*180/!dpi,xrange=[-180,180],$
            xtitle='ArcTan(!5B!3!Dy!N,!5B!3!Dx!N) [Degrees]',$
            ytitle='!3Counts',$
            bins=2,ticklen=extra.ticklen,thick=extra.thick,$
            xthick=extra.thick,ythick=extra.thick,charsize=extra.charsize,$
            charthick=extra.charthick,bin=2,xstyle=1,color=psc.black,$
                       xtickv=[-180,-90,0,90,180],xticks=4,xminor=9
   
   device,/close_file
;
   loadct,0
;  
   device,filename=dir+files.ftv+tag+eps,/encapsulate,/color,bits_per_pixel=8
   device,xsize=5,ysize=5,/inches
   pos=aspect(1.0,margin=0.08)
;
;
   plotimage,bytscl(mag.Bz[xr[0]:xr[1],yr[0]:yr[1]],min=-mx,max=mx),$
             _extra=extra,$
             xtitle='!3[Pixels]',ytitle='!3[Pixels]',position=pos+.06,$
             imgxrange=[xr[0]-d,xr[1]+d],imgyrange=[yr[0]-d,yr[1]+d]
;
   loadct,39
   dummy=fsc_color(COLORSTRUCTURE=psc,/allcolors)
   VZ=smooth(truth.VPZ,3)
;   contour,vz[xr[0]:xr[1],yr[0]:yr[1]],X[xr[0]:xr[1]],Y[yr[0]:yr[1]],$
;            levels=[0.04,0.08,0.12],color=psc.blue,/follow,thick=2,/overplot
;   contour,vz[xr[0]:xr[1],yr[0]:yr[1]],X[xr[0]:xr[1]],Y[yr[0]:yr[1]],$
;            levels=-0.04,color=psc.red,/follow,thick=2,/overplot
   BZ=smooth(mag.bz,7)
   contour,Bz[xr[0]:xr[1],yr[0]:yr[1]],X[xr[0]:xr[1]],Y[yr[0]:yr[1]],$
           color=psc.blue,levels=[0],/overplot
;
   FX=truth.FTVX
   FY=truth.FTVY
   SKIP=3L
   IX=where(X mod skip ne 0,KX)
   IY=where(Y mod skip ne 1,KY)
;;;;;;;;   FX[kill]=!values.d_nan
;;;;;;;;   FY[kill]=!values.d_nan
   for i=0L,kx-1L do FX[IX[i],*]=!values.d_nan
   for j=0L,ky-1L do FX[*,IY[j]]=!values.d_nan
;
   FTVX=Vel.FTVX
   FTVY=Vel.FTVY
   for i=0L,kx-1L do FTVX[IX[i],*]=!values.d_nan
   for j=0L,ky-1L do FTVX[*,IY[j]]=!values.d_nan
;;;;;;;;   FTVX[kill]=!values.d_nan
;;;;;;;;   FTVY[kill]=!values.d_nan
;
   
   velovect,FX[xr[0]:xr[1],yr[0]:yr[1]],FY[xr[0]:xr[1],yr[0]:yr[1]],$
            X[xr[0]:xr[1]],Y[yr[0]:yr[1]],/overplot,$
            length=skip*sqrt(2),color=psc.green,thick=3
;
   mxF=max(sqrt(FX[xr[0]:xr[1],yr[0]:yr[1]]^2+FY[xr[0]:xr[1],yr[0]:yr[1]]^2),/nan)
;
   FTVX0=FTVX
   FTVY0=FTVY
;
   neg=where(mag.bz gt 0,complement=pos)
;  FTVX[pos]=!values.D_NAN
;   FTVY[pos]=!values.D_NAN
;
   mxBU=max(sqrt(FTVX[xr[0]:xr[1],yr[0]:yr[1]]^2+FTVY[xr[0]:xr[1],yr[0]:yr[1]]^2),/nan)
   velovect,FTVX[xr[0]:xr[1],yr[0]:yr[1]],FTVY[xr[0]:xr[1],yr[0]:yr[1]],$
            X[xr[0]:xr[1]],Y[yr[0]:yr[1]],length=skip*sqrt(2)*mxBU/mxF,$
            color=psc.red,thick=2,/overplot
;
;   FTVX=FTVX0
;   FTVY=FTVY0
;
;   FTVX[neg]=!values.D_NAN
;   FTVY[neg]=!values.D_NAN
;
   mxBU=max(sqrt(FTVX[xr[0]:xr[1],yr[0]:yr[1]]^2+FTVY[xr[0]:xr[1],yr[0]:yr[1]]^2),/nan)
;      this one should be white eventially
;   velovect,FTVX[xr[0]:xr[1],yr[0]:yr[1]],FTVY[xr[0]:xr[1],yr[0]:yr[1]],$
;            X[xr[0]:xr[1]],Y[yr[0]:yr[1]],length=skip*sqrt(2)*mxBU/mxF,$
;            color=psc.white,thick=2,/overplot
;
   device,/close_file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;     flux transport velocity correlation plot
   device,filename=dir+files.ftv_corr+tag+eps,/encapsulate,/color,bits_per_pixel=8
   device,xsize=5,ysize=5,/inches
   pos=aspect(1.0,margin=0.08)+.06
;
   scl=1.d5/1.d8
;
   if (NAME eq 'DAVE') then $
      ytitle=name+'  !10J!3!Dx!NB!Dz!N (red) & !10J!3!Dy!N!5B!3!Dz!N (blue) [10!U8!N G cm/s]' $
   else $
      ytitle=name+'  !5u!3!Dx!NB!Dz!N (red) & !5u!3!Dy!N!5B!3!Dz!N (blue) [10!U8!N G cm/s]' 
;
   plot,truth.ftvx[good]*scl,Vel.FTVX[good]*scl,$
        _extra=extra,position=pos,$
        /nodata,color=psc.black,xrange=[-1.5,1.5],yrange=[-1.5,1.5],$
        xstyle=1,ystyle=1,$
    xtitle=  'ANMHD  !5U!3!Dx!N!5B!3!Dz!N (red) & !5U!3!Dy!N!5B!3!Dz!N (blue) [10!U8!N G cm/s]',$
    ytitle=ytitle
   oplot,truth.ftvx[good]*scl,Vel.FTVX[good]*scl,color=psc.red,_extra=extra
   oplot,truth.ftvy[good]*scl,Vel.FTVY[good]*scl,color=psc.blue,_extra=extra
   oplot,[-1.5,1.5],[-1.5,1.5],linestyle=2,thick=extra.thick,color=psc.black
;
   xyouts,0.45,.85,string(results.FTVX_R,format='("!9r!3!Dx!N=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
   xyouts,0.45,.8,string(results.FTVY_R,format='("!9r!3!Dy!N=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
   xyouts,0.4,0.825,spearman,align=1,color=psc.black,_extra=extra,/normal
;
   xyouts,0.75,.40+dd,string(results.FTVX_C,format='("!3C!Dx!N=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
   xyouts,0.75,.35+dd,string(results.FTVY_C,format='("!3C!Dy!N=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
    xyouts,0.7,0.375+dd,'Pearson:',align=1,color=psc.black,_extra=extra,/normal
;
   xyouts,0.75,0.25+dd,string(results.FTVX_L[1],format='("!3S!Dx!N=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
   xyouts,0.75,0.2+dd,string(results.FTVY_L[1],format='("!3S!Dy!N=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
    xyouts,0.7,0.225+dd,'Slopes:',align=1,color=psc.black,_extra=extra,/normal
;
   device,/close_file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;     flux transport angle histogram
   device,filename=dir+files.ftv_theta+tag+eps,/encapsulate,/color,bits_per_pixel=8
   device,xsize=5,ysize=5,/inches
   pos=aspect(1.0,margin=0.08)+.06
;
   X1=truth.FTVX[good] & Y1=truth.FTVY[good]
   X2=vel.FTVX[good] & Y2=vel.FTVY[good]
   M1=sqrt(X1^2+Y1^2)
   M2=sqrt(X2^2+Y2^2)
;
   cs=(X1*X2+Y1*Y2)/(M1*M2)
   sn=(X1*Y2-Y1*X2)/(M1*M2)
   theta=180.d0*atan(sn,cs)/!dpi
   theta_stats=moment(theta,/double)
;
   !P.color=psc.black
   plothist,theta,xrange=[-180,180],xtitle='!9q!3 [Degrees]',ytitle='!3Counts',$
            color=psc.black,bins=2,ticklen=extra.ticklen,thick=extra.thick,$
            xthick=extra.thick,ythick=extra.thick,charsize=extra.charsize,$
            charthick=extra.charthick,bin=2,xstyle=1,$
                       xtickv=[-180,-90,0,90,180],xticks=4,xminor=9
   
   xyouts,0.25,.875,'!3<!9q!3>='+string(results.ftv_theta[0],format='(F4.1)')+$
                       '!9'+string(177B)+'!3'+$
                       string(sqrt(results.ftv_theta[1]),format='(F4.1)')+$
                       '!9'+string(176B)+'!3!N',/normal,_extra=extra
;
;
   device,/close_file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;     perpendicular velocity correlation plot
   device,filename=dir+files.vp_corr+tag+eps,/encapsulate,/color,bits_per_pixel=8
   device,xsize=5,ysize=5,/inches
   pos=aspect(1.0,margin=0.08)+.06
;
   scl=1.d0
   perp=string("122)
   plot,truth.vpx[good]*scl,Vel.vpx[good]*scl,$
        _extra=extra,position=pos,$
        /nodata,color=psc.black,xrange=[-0.5,0.5],yrange=[-0.5,0.5],$
        xstyle=1,ystyle=1,$
    xtitle='!3ANMHD  !5V!3!D!9^!3x!N (red) & !5V!3!D!9^!3y!N (blue) & !5V!3!D!9^!3z!N (black) [km/s]',$
    ytitle=name+'  !5v!3!D!9^!3x!N (red) & !5v!3!D!9^!3y!N (blue) & !5v!3!D!9^!3z!N (black) [km/s]'
   oplot,truth.vpx[good]*scl,Vel.vpx[good]*scl,color=psc.red,_extra=extra
   oplot,truth.vpy[good]*scl,Vel.vpy[good]*scl,color=psc.blue,_extra=extra
   oplot,truth.vpz[good]*scl,Vel.vpz[good]*scl,color=psc.black,_extra=extra
   oplot,[-1.5,1.5],[-1.5,1.5],linestyle=2,color=psc.black,thick=extra.thick
;
   xyouts,0.45,.875,string(results.VPX_R,format='("!9r!3!Dx!N=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
   xyouts,0.45,.825,string(results.VPY_R,format='("!9r!3!Dy!N=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
   xyouts,0.45,.775,string(results.VPZ_R,format='("!9r!3!Dz!N=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
   xyouts,0.4,0.825,spearman,align=1,color=psc.black,_extra=extra,/normal
;
   xyouts,0.75,.3,string(results.VPX_C,format='("!3C!Dx!N=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
   xyouts,0.75,.25,string(results.VPY_C,format='("!3C!Dy!N=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
   xyouts,0.75,.2,string(results.VPZ_C,format='("!3C!Dz!N=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
   xyouts,0.7,0.25,'Pearson:',align=1,color=psc.black,_extra=extra,/normal
;
;
   xyouts,0.8,.525,string(results.VPX_L[1],format='("!3S!Dx!N=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
   xyouts,0.8,.475,string(results.VPY_L[1],format='("!3S!Dy!N=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
   xyouts,0.8,.425,string(results.VPZ_L[1],format='("!3S!Dz!N=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
   xyouts,0.75,0.475,'Slopes:',align=1,color=psc.black,_extra=extra,/normal
;
;
   device,/close_file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;     perpendicular electric field correlation plot
   device,filename=dir+files.E_corr+tag+eps,/encapsulate,/color,bits_per_pixel=8
   device,xsize=5,ysize=5,/inches
   pos=aspect(1.0,margin=0.08)+.06
;
   q1=speed_of_light
   q2=10.d0^fix(alog10(speed_of_light))
   q=q1/q2
   scl=1.d5/(1.d-4/q) ; convert statvolts/cm to V/m
   plot,truth.vpx[good]*scl,Vel.vpx[good]*scl,$
        _extra=extra,position=pos,$
        /nodata,color=psc.black,xrange=[-200,200],yrange=[-200,200],$
        xstyle=1,ystyle=1,$
    xtitle='!3ANMHD  !5E!3!D!9^!3x!N (red) & !5E!3!D!9^!3y!N (blue) & !5E!3!D!9^!3z!N (black) [V/m]',$
    ytitle=name+'  !5e!3!D!9^!3x!N (red) & !5e!3!D!9^!3y!N (blue) & !5e!3!D!9^!3z!N (black) [V/m]'
   oplot,truth.Ex[good]*scl,Vel.Ex[good]*scl,color=psc.red,_extra=extra
   oplot,truth.Ey[good]*scl,Vel.Ey[good]*scl,color=psc.blue,_extra=extra
   oplot,truth.Ez[good]*scl,Vel.Ez[good]*scl,color=psc.black,_extra=extra
   oplot,[-200,200],[-200,200],linestyle=2,color=psc.black,thick=extra.thick
;
   xyouts,0.45,.9,string(results.EX_R,format='("!9r!3!Dx!N=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
   xyouts,0.45,.85,string(results.EY_R,format='("!9r!3!Dy!N=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
   xyouts,0.45,.8,string(results.EZ_R,format='("!9r!3!Dz!N=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
   xyouts,0.4,0.85,spearman,align=1,color=psc.black,_extra=extra,/normal
;
   xyouts,0.75,.3,string(results.EX_C,format='("!3C!Dx!N=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
   xyouts,0.75,.25,string(results.EY_C,format='("!3C!Dy!N=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
   xyouts,0.75,.2,string(results.EZ_C,format='("!3C!Dz!N=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
   xyouts,0.7,0.255,'Pearson:',align=1,color=psc.black,_extra=extra,/normal
;
;
   xyouts,0.335,.7,string(results.EX_L[1],format='("!3S!Dx!N=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
   xyouts,0.335,.65,string(results.EY_L[1],format='("!3S!Dy!N=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
   xyouts,0.335,.6,string(results.EZ_L[1],format='("!3S!Dz!N=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
   xyouts,0.285,0.65,'Slopes:',align=1,color=psc.black,_extra=extra,/normal
;
   device,/close_file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;     Poynting flux
   device,filename=dir+files.poynting_corr+tag+eps,/encapsulate,/color,bits_per_pixel=8
   device,xsize=5,ysize=5,/inches
   pos=aspect(1.0,margin=0.08)+.06
   if (NAME eq 'DAVE') then $
      ytitle=NAME+'  -!6B!3!Dh!N!9'+cdot+'!3!20J!5B!3!Dz!N/(4!9p!3) [10!U10!N ergs/cm!U2!N/s]' $
   else $
      ytitle=NAME+'  -!6B!3!Dh!N!9'+cdot+'!3!6u!5B!3!Dz!N/(4!9p!3) [10!U10!N ergs/cm!U2!N/s]' 

;
   scl=1.d5/1.d10 ; convert to ergs/cm^2/s
   plot,truth.poynting[good]*scl,Vel.poynting[good]*scl,$
        _extra=extra,position=pos,$
        /nodata,color=psc.black,xrange=[-2,6],yrange=[-2,6],$
        xstyle=1,ystyle=1,$
        xtitle='!3ANMHD  -!6B!3!Dh!N!9'+cdot+'!6U!5B!3!Dz!N/(4!9p!3) [10!U10!N ergs/cm!U2!N/s]',$
        ytitle=ytitle
   oplot,truth.poynting[good]*scl,Vel.poynting[good]*scl,$
        _extra=extra,color=psc.black
   oplot,[-2,6],[-2,6],linestyle=2,thick=extra.thick,color=psc.black
;
   xyouts,0.5,.9,string(results.poynting_R,format='("!9r!3=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
;
   xyouts,0.5,.85,string(results.poynting_C,format='("!3C=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
;
   xyouts,0.5,0.8,string(results.poynting_L[1],format='("!3S=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra

   xyouts,0.45,0.9,spearman,align=1,color=psc.black,_extra=extra,/normal
   xyouts,0.45,0.85,'Pearson:',align=1,color=psc.black,_extra=extra,/normal
   xyouts,0.45,0.8,'Slope:',align=1,color=psc.black,_extra=extra,/normal
   xyouts,0.5,.15,$
           string(results.poynting_rat,format='("!3Ratio of Totals=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
;

   device,/close_file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;     Helicity flux
   device,filename=dir+files.helicity_corr+tag+eps,/encapsulate,/color,bits_per_pixel=8
   device,xsize=5,ysize=5,/inches
   pos=aspect(1.0,margin=0.08)+.06
;
   if (NAME eq 'DAVE') then $
      ytitle=name+'  -2!6A!3!Dp!N!9'+cdot+'!3!20J!5B!3!Dz!N [10!U19!N Mx!U2!N/cm!U2!N/s]' $
   else $
      ytitle=name+'  -2!6A!3!Dp!N!9'+cdot+'!3!6u!5B!3!Dz!N [10!U19!N Mx!U2!N/cm!U2!N/s]'
;
   scl=1.d10/1.d19 ; convert to [10^19 Mx^2/cm^2/s]
   plot,truth.helicity[good]*scl,Vel.helicity[good]*scl,$
        _extra=extra,position=pos,$
        /nodata,color=psc.black,xrange=[-10,5],yrange=[-10,5],$
        xstyle=1,ystyle=1,$
        xtitle='!3ANMHD  -2!6A!3!Dp!N!9'+cdot+'!6U!5B!3!Dz!N [10!U19!N Mx!U2!N/cm!U2!N/s]',$
        ytitle=ytitle
   oplot,truth.helicity[good]*scl,Vel.helicity[good]*scl,$
        _extra=extra,color=psc.black
   oplot,[-10,5],[-10,5],linestyle=2,thick=extra.thick,color=psc.black
;
   xyouts,0.5,.9,string(results.helicity_R,format='("!9r!3=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
;
   xyouts,0.5,.85,string(results.helicity_C,format='("!3C=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
   xyouts,0.5,0.8,string(results.helicity_L[1],format='("!3S=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
;
   xyouts,0.45,0.9,spearman,align=1,color=psc.black,_extra=extra,/normal
   xyouts,0.45,0.85,'Pearson:',align=1,color=psc.black,_extra=extra,/normal
   xyouts,0.45,0.8,'Slope:',align=1,color=psc.black,_extra=extra,/normal
   xyouts,0.5,.15,$
           string(results.helicity_rat,format='("!3Ratio of Totals=",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
;

   device,/close_file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;     total velocity correlation plot
   device,filename=dir+files.v_corr+tag+eps,/encapsulate,/color,bits_per_pixel=8
   device,xsize=5,ysize=5,/inches
   pos=aspect(1.0,margin=0.08)+.06
;
   scl=1.d0
   perp=string("122)
   plot,truth.vx[good]*scl,Vel.U0[good]*scl,$
        _extra=extra,position=pos,$
        /nodata,color=psc.black,xrange=[-0.4,0.4],yrange=[-0.4,0.4],$
        xstyle=1,ystyle=1,$
    xtitle='!3ANMHD  V!Dx!N (red) & V!Dy!N (blue) & V!Dz!N (black) [km/s]',$
    ytitle=name+'  v!Dx!N (red) & v!Dy!N  (blue) & v!Dz!N  (black) [km/s]'
   oplot,truth.vx[good]*scl,Vel.U0[good]*scl,color=psc.red,_extra=extra
   oplot,truth.vy[good]*scl,Vel.V0[good]*scl,color=psc.blue,_extra=extra
   oplot,truth.vz[good]*scl,Vel.W0[good]*scl,color=psc.black,_extra=extra
   oplot,[-1.5,1.5],[-1.5,1.5],linestyle=2,color=psc.black,thick=extra.thick
;
;     test null hypothesis
   VX_R=R_Correlate(truth.vx[good],truth.vpx[good])
   VY_R=R_Correlate(truth.vy[good],truth.vpy[good])
   VZ_R=R_Correlate(truth.vz[good],truth.vpz[good])
   VX_C=Correlate(truth.vx[good],truth.vpx[good])
   VY_C=Correlate(truth.vy[good],truth.vpy[good])
   VZ_C=Correlate(truth.vz[good],truth.vpz[good])
   VX_L=ladfit(truth.vx[good],truth.vpx[good])
   VY_L=ladfit(truth.vy[good],truth.vpy[good])
   VZ_L=ladfit(truth.vz[good],truth.vpz[good])
;
   VX_R=R_Correlate(truth.vx[good],vel.vpx[good])
   VY_R=R_Correlate(truth.vy[good],vel.vpy[good])
   VZ_R=R_Correlate(truth.vz[good],vel.vpz[good])
   VX_C=Correlate(truth.vx[good],vel.vpx[good])
   VY_C=Correlate(truth.vy[good],vel.vpy[good])
   VZ_C=Correlate(truth.vz[good],vel.vpz[good])
   VX_L=ladfit(truth.vx[good],vel.vpx[good])
   VY_L=ladfit(truth.vy[good],vel.vpy[good])
   VZ_L=ladfit(truth.vz[good],vel.vpz[good])
;
   xyouts,0.4,.9,string(results.VX_R,format='("!9r!3!Dx!N=",F4.2)')+$
                   string(VX_R[0],format='(",",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
   xyouts,0.4,.85,string(results.VY_R,format='("!9r!3!Dy!N=",F4.2)')+$
                   string(VY_R[0],format='(",",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
   xyouts,0.4,.8,string(results.VZ_R,format='("!9r!3!Dz!N=",F4.2)')+$
                   string(VZ_R[0],format='(",",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
   xyouts,0.35,0.85,spearman,align=1,color=psc.black,_extra=extra,/normal
;
   xyouts,0.725,.275,string(results.VX_C,format='("!3C!Dx!N=",F4.2)')+$
                   string(VX_C,format='(",",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
   xyouts,0.725,.225,string(results.VY_C,format='("!3C!Dy!N=",F4.2)')+$
                   string(VY_C,format='(",",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
   xyouts,0.725,.175,string(results.VZ_C,format='("!3C!Dz!N=",F4.2)')+$
                   string(VZ_C,format='(",",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
;   xyouts,0.675,0.225,'Pearson:',align=1,color=psc.black,_extra=extra,/normal
   xyouts,0.75,0.325,'Pearson:',align=0,color=psc.black,_extra=extra,/normal
;
;
   xyouts,0.35,.275,string(results.VX_L[1],format='("!3S!Dx!N=",F5.2)')+$
                   string(VX_L[1],format='(",",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
   xyouts,0.35,.225,string(results.VY_L[1],format='("!3S!Dy!N=",F5.2)')+$
                   string(VY_L[1],format='(",",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
   xyouts,0.35,.175,string(results.VZ_L[1],format='("!3S!Dz!N=",F5.2)')+$
                   string(VZ_L[1],format='(",",F4.2)'),$
                         align=0.0,color=psc.black,/normal,_extra=extra
   xyouts,0.3,0.225,'Slopes:',align=1,color=psc.black,_extra=extra,/normal
;
;
   device,/close_file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
    set_plot,'x'
;   stop
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
end
