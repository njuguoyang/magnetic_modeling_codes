pro dave4vm_calibration

  XM0=[10.0000,3.00000]
  method='DAVE4VM '
  ergscl=1.d28
  helicityscl=1.d37
  DIR='SCHUCK/'
  dave4vm=mrdfits(DIR+'dave4vm_optimize.fts',1)
  opt=mrdfits(DIR+'dave4vm.fts',3)
;
  set_plot,'x'
  plot,dave4vm.window_size[0],dave4vm.helicity_d/helicityscl
  YC=!Y.CRANGE
  !P.font=1
  th=4
  extra={charsize:1.35,ticklen:-0.02,thick:th,charthick:2,xthick:th,ythick:th}
  set_plot,'ps'
  loadct,39
  device,filename=DIR+'f1d.eps',/encapsulate,/color,bits_per_pixel=8
  !X.MARGIN=[6,5]
  yr=[min([0,dave4vm.poynting_d,dave4vm.poynting_anmhd],MAX=MX),MX]/ergscl
  plot,dave4vm.window_size[0],dave4vm.poynting_d/ergscl,$
        xtitle='Window Size [Pixels]',xstyle=1,yrange=yr,$
        ytitle=method+'Poynting Flux [10!U28!N ergs/s]',$
       _extra=extra,color=0,ystyle=8,position=aspect(1.d0)
  NN=N_ELEMENTS(dave4vm)
  oplot,dave4vm.window_size[0],replicate(dave4vm[0].poynting_anmhd/ergscl,NN),$
            _extra=extra,linestyle=2
;
  XR1=!X.crange
  YR1=!Y.crange
  oplot,opt.window_size,YR1,_extra=extra,linestyle=3
  hc=254
  yr=[min([0,helicityscl,dave4vm.helicity_d,dave4vm.helicity_anmhd]/helicityscl,MAX=MX),MX]
;  print,yr
  axis,XR1[1],0,/yaxis,/data,_extra=extra,/save,color=hc,$
       ytitle=method+'Helicity Flux [10!U37!N Wb!U2!N/s]',yrange=yr
  oplot,dave4vm.window_size[0],dave4vm.helicity_d/helicityscl,_extra=extra,$
        color=hc
  oplot,dave4vm.window_size[0],replicate(dave4vm[0].helicity_anmhd/helicityscl,NN),$
            _extra=extra,linestyle=2,color=hc
  device,/close_file

  plot,dave4vm.window_size[0],dave4vm.induction_L[1]
  YC=!Y.CRANGE
  device,filename=DIR+'f1b.eps',/encapsulate
  !X.MARGIN=[7,7]
  plot,dave4vm.window_size[0],dave4vm.induction_c,$
        xtitle='Window Size [Pixels]',xstyle=1,$
        ytitle=method+'Correlation',yrange=[-1,0],$
       _extra=extra,color=0,ystyle=9,position=aspect(1.d0)
  oplot,dave4vm.window_size[0],dave4vm.induction_r,_extra=extra,linestyle=2
  XR1=!X.crange
  YR1=!Y.crange
  oplot,opt.window_size,YR1,_extra=extra,linestyle=3
  hc=254
  yr=[min([dave4vm.induction_L[1],-.92,-1.02],MAX=MX),MX]
  axis,XR1[1],0,/yaxis,/data,yrange=YR,_extra=extra,ystyle=1,/save,color=hc,$
       ytitle=method+'Slope'
  oplot,dave4vm.window_size[0],dave4vm.induction_L[1],_extra=extra,$
        color=hc
  corners=[0.823780   ,  0.704488  ,   0.379254  ,   0.808425]
;  print,corners
  corners[2*indgen(2)+1]=corners[2*indgen(2)+1]-0.05
  corners[2*indgen(2)]=corners[2*indgen(2)]+0.05
  pos=corners[2:3]
;  print,pos
  legend,['Pearson C','Spearman !9r!3'],linestyle=[0,2],color=0,box=0,charsize=1.2,thick=extra.thick,charthick=extra.charthick,corners=corners,position=pos,/normal
;  print,corners
  device,/close_file
;
;
;
  set_plot,'x'
  !X.MARGIN=XM0
 
end
