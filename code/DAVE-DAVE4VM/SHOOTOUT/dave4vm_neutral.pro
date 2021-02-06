pro dave4vm_neutral,welsch=welsch,dave=dave,no_elliptic2d=no_elliptic2d,_extra=extra
;
   DIR='SCHUCK/'                                        ;new mask?     |B|>370
   if (N_ELEMENTS(welsch) eq 1) then begin              ;original mask |bz|>370
      if (welsch eq 1) then DIR='WELSCH/' else welsch=0
   endif else welsch=0
;

   set_plot,'x'
   mag=shootout('./old_shootout/data_for_analysis.sav',0,truth=truth,no_elliptic2d=no_elliptic2d)
;
   angle=20*!dpi/180.d0                     ; angle between LOS and z-axis
   angle=0.d0

   @optimum_windows
   WINDOW_SIZE=DAVE4VM_WINDOW_SIZE
;   
   mag=compute_los_velocity(mag,truth,angle)
;
   magg=mag
   if (N_ELEMENTS(dave) eq 1) then begin              
      if (dave eq 1) then begin
;           set up for comparison with dave=dave4vm equivalence
         WINDOW_SIZE=DAVE_WINDOW_SIZE  ; this should be set to 
                                       ;  DAVE's window size 
                                       ;  NOT DAVE4VM's window size
         magg=mag
         magg.bx=0.d0
         magg.by=0.d0
         magg.bxx=0.d0
         magg.byx=0.d0
         magg.bxy=0.d0
         magg.byy=0.d0
      endif
   endif 
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
   vel=dave4vm(magg,window_size,/double)
   tm=systime(1,/seconds)
   dtm=tm-tm0
   print,'*** Time to execute DAVE4VM: ',dtm,' WINDOW_SIZE=',WINDOW_SIZE
;
   vel=compute_fluxes(mag,vel)

;
   estimated=neutral_line(mag,vel.u0,vel.v0,vel.w0)
   ground=neutral_line(mag,truth.vx,truth.vy,truth.vz)
   JUNK=where((abs(mag.BZ) gt mag.thr),complement=bad,ND)
;
   
   plot,ground.FU_PERP[junk],estimated.FU_PERP[junk],psym=3,_extra=extra
   oplot,ground.FU_PAR[junk],estimated.FU_PAR[junk],psym=3,_extra=extra
   plot,ground.FV_PERP[junk],estimated.FV_PERP[junk],psym=3,_extra=extra
   oplot,ground.FV_PAR[junk],estimated.FV_PAR[junk],psym=3,_extra=extra
   plot,ground.FV_PERP[junk]+ground.FU_PERP[junk],estimated.FV_PERP[junk]+estimated.FU_PERP[junk],psym=2

;

end
