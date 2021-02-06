pro dave_results,welsch=welsch,no_elliptic2d=no_elliptic2d,_extra=extra
;
   DIR='SCHUCK/'                                        ;new mask?     |B|>370
   if (N_ELEMENTS(welsch) eq 1) then begin              ;original mask |bz|>370
      if (welsch eq 1) then DIR='WELSCH/' else welsch=0
   endif else welsch=0


;
   tag='a'
   set_plot,'x'
;
   mag=shootout('./old_shootout/data_for_analysis.sav',0,truth=truth,no_elliptic2d=no_elliptic2d)
   B2=(mag.Bx^2+mag.By^2+mag.Bz^2)
   sz=size(mag.Bz)

   dum=dblarr(sz[1],sz[2])
   nvel={U0:dum,V0:dum,W0:dum,WINDOW_SIZE:[0,0]}
   ovel=nvel
   nan=0.d0 ;;; REMOVE
   fvel=replicate({U0:nan,V0:nan,W0:nan,WINDOW_SIZE:[0,0]},sz[1],sz[2]) ;;; REMOVE
;
; *****************************************************************************
; *****************************************************************************
;   DO NOT CHANGE THIS CODE!
; *****************************************************************************
; *****************************************************************************
;     original shootout results 
   odave=mrdfits('old_shootout/original_shootout.fts',1)
   scl=1.d5
;
;     best attempt a replicating the original data  
   old_dave,mag.bz,mag.bzt,mag.bzx*mag.dx,mag.bzy*mag.dy,1.d0,1.d0,odave.window_size,vel,aperture,sv=sv,sigma2=sigma2,chi2=chi2,/errors,missing_value=0.,/double,threshold=0.,advect=0 
   odave.UFX=vel[0,*,*]*mag.dx*scl
   odave.UFY=vel[1,*,*]*mag.dy*scl
;      convert 2D flux transport velocities into '3D' plasma velocities
   dot=(odave.UFX*mag.Bx+odave.UFY*mag.By)
   ovel.U0=(odave.UFX-dot*mag.bx/B2)/scl
   ovel.V0=(odave.UFY-dot*mag.by/B2)/scl
   ovel.W0=-(dot*mag.bz/B2)/scl
   ovel.window_size=odave.window_size
; *****************************************************************************
; *****************************************************************************
   goto,natural
;  completely rescaled example
   maggit=mag
   maggit.bzx=mag.bzx*mag.dx
   maggit.bzy=mag.bzy*mag.dy
   maggit.bzt=mag.bzt*mag.dt
   maggit.dx=1.d0
   maggit.dy=1.d0
   dave=dave(maggit,odave.window_size,AM=AM,/double)
   dave.U0=dave.U0*mag.dx*1.d5/mag.dt
   dave.V0=dave.V0*mag.dy*1.d5/mag.dt
   print,'RESCALED'
;
   goto,rescale
;
;  OR
;
natural:
   @optimum_windows
   window_size=DAVE_WINDOW_SIZE
   print,'NATURAL UNITS'
;
   tm0=systime(1,/seconds)
   dave=dave(mag,window_size,AM=BM,/double)
   tm=systime(1,/seconds)
   dtm=tm-tm0
   print,'*** Time to execute DAVE: ',dtm,' WINDOW_SIZE=',WINDOW_SIZE
   dave.U0=dave.U0*1.d5
   dave.V0=dave.V0*1.d5
;
;
rescale:
;     convert 2D flux transport velocities into '3D' plasma velocities
   dot=(dave.U0*mag.Bx+dave.V0*mag.By)
   nvel.U0=(dave.U0-dot*mag.bx/B2)/scl
   nvel.V0=(dave.V0-dot*mag.by/B2)/scl
   nvel.W0=-(dot*mag.bz/B2)/scl
   nvel.window_size=dave.window_size
;
   print,'OLD THRESHOLDS'
   mwrfits,dave,'compare.fts',/create
;;
   INDEX=where((abs(mag.BZ) gt mag.thr) and $
               (finite(nvel.u0,/nan) ne 1) and $
               (finite(nvel.V0,/nan) ne 1),complement=bad,ND)
;
   print,R_correlate(truth.vpx[index],nvel.u0[index])
   print,R_correlate(truth.vpy[index],nvel.v0[index])
   print,R_correlate(truth.vpz[index],nvel.w0[index])

   print,correlate(truth.vpx[index],nvel.u0[index])
   print,correlate(truth.vpy[index],nvel.v0[index])
   print,correlate(truth.vpz[index],nvel.w0[index])
;
;     set up comparison region
   xr=[95,205]
   yr=[90,200]
   mask=fix(mag.bz)
   mask[*]=0
   mask[xr[0]:xr[1],yr[0]:yr[1]]=1
;
;     output old_dave results
   ovel=compute_fluxes(mag,ovel)
   oresults=compare_velocities(mag,truth,ovel,mask=mask,good=good,welsch=welsch)
;
   oresults=add_tag(oresults,xr,"XR")
   oresults=add_tag(oresults,yr,"YR")
;
;   output_results,oresults,mag,ovel
;   output_results,oresults,mag,ovel,filename='old_shootout/old_dave_shootout_results.txt'
;
;     output dave results
   nvel=compute_fluxes(mag,nvel)
   nresults=compare_velocities(mag,truth,nvel,mask=mask,good=good,welsch=welsch)

   nresults=add_tag(nresults,xr,"XR")
   nresults=add_tag(nresults,yr,"YR")
;
   output_results,nresults,mag,nvel
   output_results,nresults,mag,nvel,filename=DIR+'dave.txt'
   make_plots,dir,mag,truth,nvel,good,'a','DAVE',nresults,_extra=extra
   mwrfits,mag,DIR+'dave.fts',/create
   mwrfits,truth,DIR+'dave.fts'
   mwrfits,nvel,DIR+'dave.fts'
   mwrfits,nresults,DIR+'dave.fts'

end
