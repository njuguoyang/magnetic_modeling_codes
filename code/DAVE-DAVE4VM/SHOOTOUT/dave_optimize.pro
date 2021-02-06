pro dave_optimize,no_elliptic2d=no_elliptic2d

   mag=shootout('old_shootout/data_for_analysis.sav',0,truth=truth,no_elliptic2d=no_elliptic2d)
   DIR='SCHUCK/'
;
   angle=0.d0
   Beta=0.d-3
;
   NW=25
   NB=1
;
   mag=compute_los_velocity(mag,truth,angle)
   b2=sqrt(mag.Bx^2+mag.by^2+mag.bz^2)
;
;
   INDEX=where(abs(mag.BZ) gt mag.thr,complement=bad,ND)
   print,ND,mag.thr,format='("ND=",I4," data with |Bz|>",F6.2," Gauss")'
   plot,mag.bzt[index],(truth.ftvXx+truth.ftvYY)[index],psym=2
   print,'ANMHD=',correlate(mag.bzt[index],(truth.ftvXX+truth.ftvYY)[index])
;
;     set up comparison region
   xr=[95,205]
   yr=[90,200]
   mask=fix(mag.bz)
   mask[*]=0
   mask[xr[0]:xr[1],yr[0]:yr[1]]=1
;
   wins=rebin(lindgen(NW)*2+3,NW,NB,1)
   if (NB gt 1) then begin
      db=alog10(1.d-5)/(NB-1)
      B=transpose(rebin(10.d0^(dindgen(NB)*db),NB,NW))
   endif else B=dblarr(NW,NB)

   window_size=wins[0]
   Beta=B[0]
   ftvs=dave(mag,window_size,/double)
;
   scl=1.d0
   sz=size(mag.bz)
   nan=replicate(!values.d_nan,sz[1],sz[2])
   
   vel={u0:nan,v0:nan,w0:nan,window_size:[0L,0L]}
   dot=(ftvs.U0*mag.Bx+ftvs.V0*mag.By)
   vel.U0=(ftvs.U0-dot*mag.bx/B2)/scl
   vel.V0=(ftvs.V0-dot*mag.by/B2)/scl
   vel.W0=-(dot*mag.bz/B2)/scl
   vel.window_size=window_size

;
   vel=compute_fluxes(mag,vel)
   results=compare_velocities(mag,truth,vel,mask=mask,good=good)
   results=add_tag(results,vel[0,0].window_size,'window_size')
   work=replicate(results,NW,NB)
;
   bsclx=1.d0
   bscly=1.d0
;
   for j=0L,NB-1L do begin
       for i=0L,NW-1L do begin
           if ((i eq 0) and (j eq 0)) then goto,skip
           window_size=wins[i,j]
           print,'DAVE WINDOW SIZE=',WINDOW_SIZE
;
           ftvs=dave(mag,window_size,/double)
;
;              convert 2D flux transport velocities into '3D' plasma velocities
           dot=(ftvs.U0*mag.Bx+ftvs.V0*mag.By)
           vel.U0=(ftvs.U0-dot*mag.bx/B2)/scl
           vel.V0=(ftvs.V0-dot*mag.by/B2)/scl
           vel.W0=-(dot*mag.bz/B2)/scl
           vel.window_size=ftvs.window_size
;

           vel=compute_fluxes(mag,vel)
           results=compare_velocities(mag,truth,vel,mask=mask,good=good)
           results=add_tag(results,vel[0,0].window_size,'window_size')
           work[i,j]=results
           output_results,results,mag,vel
           mwrfits,work,DIR+'dave_optimize.fts',/create
skip:
       endfor
   endfor
;

end
