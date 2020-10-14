pro dave4vm_results,welsch=welsch,dave=dave,no_elliptic2d=no_elliptic2d,_extra=extra
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
   vel=dave4vm(magg,window_size,AM=AM,KERNEL=KERNEL,/double)
   tm=systime(1,/seconds)
   dtm=tm-tm0
   print,'*** Time to execute DAVE4VM: ',dtm,' WINDOW_SIZE=',WINDOW_SIZE
;
   vel=compute_fluxes(mag,vel)
   results=compare_velocities(mag,truth,vel,mask=mask,good=good,welsch=welsch)
   results=add_tag(results,xr,"XR")
   results=add_tag(results,yr,"YR")
   make_plots,dir,mag,truth,vel,good,'b','DAVE4VM',results,_extra=extra
   output_results,results,mag,vel,filename=DIR+'dave4vm.txt'
   mwrfits,mag,DIR+'dave4vm.fts',/create ; ext#1  mag data
   mwrfits,truth,DIR+'dave4vm.fts'       ; ext#2  "correct" answers from MHD
   mwrfits,vel,DIR+'dave4vm.fts'         ; ext#3  estimates from DAVE4VM
   mwrfits,results,DIR+'dave4vm.fts'     ; ext#4  comparisons beween correct and estimates
   mwrfits,kernel,DIR+'dave4vm.fts'      ; ext#5  kernels used for convolutions   
   mwrfits,AM,DIR+'dave4vm.fts'          ; ext#6  "A" matrix for Jacob's cross-comparison
;

end
