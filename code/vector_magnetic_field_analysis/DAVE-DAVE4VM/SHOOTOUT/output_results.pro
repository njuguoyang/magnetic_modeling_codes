pro output_results,results,mag,vel,filename=filename
;
   sz=size(filename,/structure)
   unit=-1
   if (N_ELEMENTS(sz) eq 0) then goto,stdout
   if (sz.type ne 7) then goto,stdout
   if (strlen(filename) eq 0) then goto,stdout
   openw,unit,filename,/get_lun
;
stdout:
;
   printf,unit,'=============================================================='
   printf,unit,'WINDOW_SIZE=',vel[0].window_size
   CHK=where(TAG_NAMES(vel) eq 'BETA',NCHK)
   if (NCHK eq 1) then printf,unit,'BETA=',vel[0].beta
   printf,unit
   printf,unit,"FLUX TRANSPORT VELOCITIES"
   printf,unit,'FTVX_R=',results.FTVX_R;,R_correlate(vel[good].BZUX,truth[good].BZUX)
   printf,unit,'FTVY_R=',results.FTVY_R;,R_correlate(vel[good].BZUY,truth[good].BZUY)
   printf,unit,'FTVX_C=',results.FTVX_C;,correlate(vel[good].BZUX,truth[good].BZUX,/double)
   printf,unit,'FTVY_C=',results.FTVY_C;,correlate(vel[good].BZUY,truth[good].BZUY,/double)
   printf,unit,'FTVX_L=',results.FTVX_L[1]
   printf,unit,'FTVY_L=',results.FTVY_L[1]
   printf,unit
   printf,unit,'Fractional Vector Error=',results.DFTV[0],'+/-',$
                                                     sqrt(results.DFTV[1])
   printf,unit,' Fractional Speed Error=',results.DFTV_SPEED[0],'+/-',$
                                          sqrt(results.DFTV_SPEED[1])
   printf,unit,'     Vector Correlation=',results.C_FTV[0],'+/-',$
                                          sqrt(results.C_FTV[1])
   printf,unit,'  Direction Correlation=',results.CS_FTV[0],'+/-',$
                                          sqrt(results.CS_FTV[1])
   printf,unit
   printf,unit,"POYNTING FLUX"
   printf,unit,'POYNTING_R=',results.POYNTING_R;,R_correlate(vel[good].Poynting,truth[good].Poynting)  
   printf,unit,'POYNTING_C=',results.POYNTING_C;,correlate(vel[good].Poynting,truth[good].Poynting,/double)  
   printf,unit,'POYNTING_L=',results.POYNTING_L[1]  
   printf,unit,'POYNTING_D=',results.poynting_D
   printf,unit,'POYNTING_T=',results.poynting_T
   printf,unit,'POYNTING_RAT=',results.POYNTING_RAT
   printf,unit,'POYNTING_ANMHD=',results.POYNTING_ANMHD
   printf,unit
   printf,unit,"INDUCTION EQUATION"
   printf,unit,'INDUCTION_R=',results.induction_R;,R_correlate(-mag.bzt[good],vel[good].bzuxx+vel[good].bzuyy)
   printf,unit,'INDUCTION_C=',results.induction_C;,correlate(-mag.bzt[good],vel[good].bzuxx+vel[good].bzuyy,/double)
   printf,unit,'INDUCTION_L=',results.induction_L
   printf,unit,'Total=',results.induction_T
   printf,unit
   printf,unit,'HELICITY'   
   printf,unit,'HELICITY_R=',results.helicity_R
   printf,unit,'HELICITY_C=',results.helicity_C
   printf,unit,'HELICITY_L=',results.helicity_L[1]
   printf,unit,'HELICITY_D=',results.helicity_D
   printf,unit,'HELICITY_T=',results.helicity_T
   printf,unit,'HELICITY_RAT=',results.helicity_RAT
   printf,unit,'HELICITY_ANMHD=',results.helicity_ANMHD
   printf,unit
   printf,unit,'TOTAL VELOCITY'
   printf,unit,'VX_R=',results.vx_r
   printf,unit,'VX_C=',results.vx_c
   printf,unit,'VX_L=',results.vx_L[1]
   printf,unit
   printf,unit,'VY_R=',results.vy_r
   printf,unit,'VY_C=',results.vy_c
   printf,unit,'VY_L=',results.vy_l[1]
   printf,unit
   printf,unit,'VZ_R=',results.vz_r
   printf,unit,'VZ_C=',results.vz_c
   printf,unit,'VZ_L=',results.vz_l[1]
   printf,unit
   printf,unit,'PERPENDICULAR VELOCITY'
   printf,unit,'VPX_R=',results.VPx_r
   printf,unit,'VPX_C=',results.VPx_c
   printf,unit,'VPX_L=',results.VPx_l[1]
   printf,unit
   printf,unit,'VPY_R=',results.VPy_r
   printf,unit,'VPY_C=',results.VPy_c
   printf,unit,'VPY_L=',results.VPy_l[1]
   printf,unit
   printf,unit,'VPZ_R=',results.VPz_r
   printf,unit,'VPZ_C=',results.VPz_c
   printf,unit,'VPZ_L=',results.VPz_l[1]
   printf,unit
   printf,unit,'PARALLEL VELOCITY'
   printf,unit,'VB_R=',results.VB_r
   printf,unit,'VB_C=',results.VB_c
   printf,unit,'VB_L=',results.VB_l[1]
   printf,unit
   printf,unit,'Fractional Vector Error=',results.DVP[0],'+/-',$
                                                     sqrt(results.DVP[1])
   printf,unit,' Fractional Speed Error=',results.DVP_SPEED[0],'+/-',$
                                          sqrt(results.DVP_SPEED[1])
   printf,unit,'     Vector Correlation=',results.C_VP[0],'+/-',$
                                          sqrt(results.C_VP[1])
   printf,unit,'  Direction Correlation=',results.CS_VP[0],'+/-',$
                                          sqrt(results.CS_VP[1])
   printf,unit
   printf,unit,'ELECTRIC FIELD'
   printf,unit,'EX_R=',results.Ex_r
   printf,unit,'EX_C=',results.Ex_c
   printf,unit,'EX_L=',results.Ex_l[1]
   printf,unit
   printf,unit,'EY_R=',results.Ey_r
   printf,unit,'EY_C=',results.Ey_c
   printf,unit,'EY_L=',results.Ey_l[1]
   printf,unit
   printf,unit,'EZ_R=',results.Ez_r
   printf,unit,'EZ_C=',results.Ez_c
   printf,unit,'EZ_L=',results.Ez_l[1]
;
   if (unit ne -1) then begin
      close,unit
      free_lun,unit
   endif
;
end
