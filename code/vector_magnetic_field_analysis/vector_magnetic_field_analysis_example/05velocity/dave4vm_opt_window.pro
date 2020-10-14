;+
;PURPOSE: To find the optimized window size for DAVE4VM method.
;-

pro dave4vm_opt_window,thr=thr,verbose=verbose,_extra=extra
;
DIR='GUO2/'
common constants,speed_of_light
if (N_ELEMENTS(speed_of_light) eq 0) then speed_of_light=2.998d10 ; cm/s
if (N_ELEMENTS(thr) eq 0) then thr=370.d0 else thr=double(thr)
solar_r=6.955d5 ;solar radius in km
;
file=find_file('GUO2/mapbxyz*')
;
nn=n_elements(file)-1
nn=5    ;for test
nw = 13
wins = lindgen(nw)*3+13
pearson = fltarr(nw)
slope = fltarr(nw)
for i=nn-1,nn-1 do begin
   ;print,'*** ',i,100*i/(1.*nn),'%'
   file1=file[i]
   file2=file[i+1]
   print,'*** file1:',file1
   print,'*** file2:',file2
   restore,file1,verbose=verbose,_extra=extra
   temp=pb0r(smapbz.time)
   DX=double(smapbz.dx*solar_r/temp[2]/60.0)  ;in kilometers
   DY=DX
   bz_start=smapbz.data
   bx_start=smapbx.data
   by_start=smapby.data
   internal_time=anytim2utc(smapbz.time)      ;we do not have to consider about the leap seconde here. If the 
                                              ;time ranges span over one day, the leap second maybe have to be 
                                              ;considered in some cases.
   t_start=double(internal_time.time/1000.0)  ;in seconds
   restore,file2,verbose=verbose,_extra=extra
   bz_stop=smapbz.data
   bx_stop=smapbx.data
   by_stop=smapby.data
   internal_time=anytim2utc(smapbz.time)  
   t_stop=double(internal_time.time/1000.0)
;
;
   DT=double(T_stop-T_start)
   BZT=(BZ_stop-BZ_start)/DT
   sz=size(BZT)
;
;     use time_centered spatial variables
   BX=(BX_stop+BX_start)/2.d0
   BY=(BY_stop+BY_start)/2.d0
   BZ=(BZ_stop+BZ_start)/2.d0
;   
;     compute 5-point optimized derivatives
   odiffxy5,BX,BXX,BXY
   odiffxy5,BY,BYX,BYY
   odiffxy5,BZ,BZX,BZY
;
   MAG={BZT:BZT,BX:BX,BXX:BXX/DX,BXY:BXY/DY,BY:BY,BYX:BYX/DX,BYY:BYY/DY,$
        BZ:BZ,BZX:BZX/DX,BZY:BZY/DY,DX:DX,DY:DY,DT:DT,THR:THR}
;     unit vector (direction of the magnetic field)
   B=sqrt(bx^2+by^2+bz^2)
   qx=bx/B
   qy=by/B
   qz=bz/B
;
   mag=add_tag(mag,qx,'QX')
   mag=add_tag(mag,qy,'QY')
   mag=add_tag(mag,qz,'QZ')
;
;   @optimum_windows
;   WINDOW_SIZE=DAVE4VM_WINDOW_SIZE
 for j=0,nw-1 do begin
 print,'*** ',j,100*j/(1.*nw),'%'
   WINDOW_SIZE = wins[j]

   tm0=systime(1,/seconds)
   vel=dave4vm(mag,window_size,/double)
   tm=systime(1,/seconds)
   dtm=tm-tm0
   print,'*** Time to execute DAVE4VM: ',dtm,' Secondes. WINDOW_SIZE=',WINDOW_SIZE
;
   vel=compute_fluxes(mag,vel)
   INDEX=where(abs(mag.BZ) gt thr,complement=bad,ND)
   print,ND,thr,format='("ND=",I6," data with |Bz|>",F6.2," Gauss")'
   plot,mag.bzt[index],(vel.ftvXX+vel.ftvYY)[index],psym=2
   pearson[j] = correlate(mag.bzt[index],(vel.ftvXX+vel.ftvYY)[index])
   temp = poly_fit(mag.bzt[index],(vel.ftvXX+vel.ftvYY)[index],1)
   slope[j] = temp[0,1] 
   print,'DAVE4VM Correlation:',pearson[j]
   print,'DAVE4VM Slope:',slope[j]
 endfor
 print,'*** ',j,100*j/(1.*nw),'%'
 ;print,'*** ',i,100*i/(1.*nn),'%'
endfor
;
save,wins,pearson,slope,filename='dave4vm_opt_window.sav'
window,/free
!p.background = 255
plot,wins,pearson,xtitle='Window Size [Pixels]',ytitle='DAVE4VM Correlation',charsize=2,ystyle=1+8,color=0,pos=[0.15,0.15,0.85,0.96],thick=2
axis,yaxis=1,yrange=[-0.1,-0.06],/sav,color=100,charsize=2.0,ystyle=1,ytitle='DAVE4VM Slope'
oplot,wins,slope,color=0
write_png,'dave4vm_opt_window.png',tvrd()
!p.background = 0

end
