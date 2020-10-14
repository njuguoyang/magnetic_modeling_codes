pro dave4vm_optimize,no_elliptic2d=no_elliptic2d

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
;
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
   wins=rebin(lindgen(NW)*2+3,NW,NB,1)
   if (NB gt 1) then begin
      db=alog10(1.d-5)/(NB-1)
      B=transpose(rebin(10.d0^(dindgen(NB)*db),NB,NW))
   endif else B=dblarr(NW,NB)

   window_size=wins[0]
   Beta=B[0]
   vel=dave4vm(mag,window_size,/double)
   vel=compute_fluxes(mag,vel)
   results=compare_velocities(mag,truth,vel,mask=mask,good=good)
   results=add_tag(results,window_size,'window_size')
   results=add_tag(results,beta,'beta')
   work=replicate(results,NW,NB)
;
   for j=0L,NB-1L do begin
       for i=0L,NW-1L do begin
           if ((i eq 0) and (j eq 0)) then goto,skip
           print,i,j
           window_size=wins[i,j]
           print,'DAVE4VM WINDOW SIZE=',window_size
           beta=B[i,j]
           vel=dave4vm(mag,window_size,/double)
           vel=compute_fluxes(mag,vel)
           results=compare_velocities(mag,truth,vel,mask=mask,good=good)
           results=add_tag(results,window_size,'window_size')
           results=add_tag(results,beta,'beta')
           work[i,j]=results
           mwrfits,work,DIR+'dave4vm_optimize.fts',/create
           output_results,results,mag,vel
skip:
       endfor
   endfor
;

end
