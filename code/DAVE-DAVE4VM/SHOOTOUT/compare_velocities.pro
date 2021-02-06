function compare_velocities,mag,truth,vel,mask=mask,good=good,welsch=welsch

   sz=size(mag.bz)
   if (N_ELEMENTS(mask) eq 0) then mask=intarr(sz[1],sz[2])
;
   B=sqrt(mag.Bx^2+mag.By^2+mag.Bz^2)
;
   criteria=B
   cstring='|B|'
;
   if (check_keyword(welsch) eq 1) then begin
      if (welsch ne 0) then begin     
         criteria=abs(mag.bz)
         cstring='|Bz|'
      endif
   endif
;
   good=where(criteria gt mag.thr ,complement=bad,NG)
   print,NG,cstring,mag.thr,$
         format='("ND=",I4," data with ",A4,">",F6.2," Gauss (in simulation)")'
;
   good=where( (criteria gt mag.thr) and $
              (mask eq 1)  ,complement=bad,NG)
   print,NG,cstring,mag.thr,$
         format='("ND=",I4," data with ",A4,">",F6.2," Gauss (in mask)")'
;
   junk=Vel.FTVX
   kill=where(mask ne 1)
   junk[kill]=!values.d_nan
   
   good=where( (criteria gt mag.thr) and (finite(junk,/nan) ne 1)  $
                                        and (finite(vel.ftvx,/nan) ne 1) $
                                        and (finite(vel.ftvxx,/nan) ne 1) $
                                        and (finite(vel.ftvy,/nan) ne 1)  $
                                        and (finite(vel.ftvyy,/nan) ne 1) $
                                        and (finite(mag.Ax,/nan) ne 1) $
                                        and (finite(mag.Ay,/nan) ne 1) ,NG)
; 

;
   print,NG,cstring,mag.thr,$
      format='("ND=",I4," data with ",A4,">",F6.2," Gauss (in mask and finite)")'
;
;;;;;;;;;;;;;; Flux Transport Velocities ;;;;;;;;;;;;;; 
;
   FTVX_R=R_correlate(truth.FTVX[good],vel.FTVX[good],$
                       D=DX,PROBD=PROBDX,ZD=ZDX)
   FTVY_R=R_correlate(truth.FTVY[good],vel.FTVY[good],$
                       D=DX,PROBD=PROBDY,ZD=ZDY)
   FTVX_C=correlate(truth.FTVX[good],vel.FTVX[good])
   FTVY_C=correlate(truth.FTVY[good],vel.FTVY[good])
   
   FTVX_L=ladfit(truth.FTVX[good],vel.FTVX[good],/double)
   FTVY_L=ladfit(truth.FTVY[good],vel.FTVY[good],/double)
;
;
   FTV=sqrt(truth.FTVX[good]^2+truth.FTVY[good]^2)
   FTVEL=sqrt(vel.FTVX[good]^2+vel.FTVY[good]^2)
   metric=sqrt((truth.FTVX[good]-vel.FTVX[good])^2$
                   +(truth.FTVY[good]-vel.FTVY[good])^2)
;
;
;      compute fractional error
   DFTV=moment(metric/FTV,/double)
;
;      fractional error in magnitude
   metric=(FTVEL-FTV)/FTV
   DFTV_SPEED=moment(metric,/double)
;
;      vector correlation
   DOT=truth.FTVX[good]*vel.FTVX[good]+truth.FTVY[good]*vel.FTVY[good]
   metric=DOT/sqrt(mean(FTV^2,/double)*mean(FTVEL^2,/double))
   C_FTV=moment(metric,/double)
;
;      direction correlation
   metric=DOT/(FTV*FTVEL)
   CS_FTV=moment(metric,/double)
;
   X1=truth.FTVX[good] & Y1=truth.FTVY[good]
   X2=vel.FTVX[good] & Y2=vel.FTVY[good]
   M1=sqrt(X1^2+Y1^2)
   M2=sqrt(X2^2+Y2^2)
;
   cs=(X1*X2+Y1*Y2)/(M1*M2)
   sn=(X1*Y2-Y1*X2)/(M1*M2)
   theta=180.d0*atan(sn,cs)/!dpi
   FTV_THETA=moment(theta,/double)
;
;      weighted cosine
   NORM=total(FTVEL,/double)
   WCOS_FTV=total(cos(theta*!dpi/180.d0)*FTVEL/norm,/double)
;    
;;;;;;;;;;;;;;       Poynting Flux       ;;;;;;;;;;;;;; 
;
   POYNTING_R=(R_correlate(truth.Poynting[good],vel.Poynting[good]))[0]  
   POYNTING_C=correlate(truth.Poynting[good],vel.Poynting[good])  
   POYNTING_L=ladfit(truth.Poynting[good],vel.Poynting[good])  
   POYNTING_D=total(vel.poynting[good]*1.d5,/double)*mag.dx*mag.dy*1.d10
   POYNTING_T=total(truth.poynting[good]*1.d5,/double)*mag.dx*mag.dy*1.d10
   POYNTING_ANMHD=total(truth.poynting*1.d5,/double)*mag.dx*mag.dy*1.d10
   POYNTING_RAT=POYNTING_D/POYNTING_T
;
;;;;;;;;;;;;;;       Helicity Flux       ;;;;;;;;;;;;;; 
;
   HELICITY_R=(R_correlate(vel.Helicity[good],truth.Helicity[good]))[0]  
   HELICITY_C=correlate(vel.Helicity[good],truth.Helicity[good])  
   HELICITY_L=ladfit(truth.Helicity[good],vel.Helicity[good])  
   HELICITY_D=total(vel.helicity[good]*1.d10,/double)*mag.dx*mag.dy*1.d10
   HELICITY_T=total(truth.helicity[good]*1.d10,/double)*mag.dx*mag.dy*1.d10
   HELICITY_ANMHD=total(truth.helicity*1.d10,/double,/nan)*mag.dx*mag.dy*1.d10
   HELICITY_RAT=HELICITY_D/HELICITY_T
;
;;;;;;;;;;;;;;    Induction Equation     ;;;;;;;;;;;;;; 
;
   INDUCTION_R=(R_correlate(mag.bzt[good],vel.ftvxx[good]+vel.ftvyy[good]))[0]
   INDUCTION_C=correlate(mag.bzt[good],vel.ftvxx[good]+vel.ftvyy[good])
   INDUCTION_L=ladfit(mag.bzt[good],vel.ftvxx[good]+vel.ftvyy[good])
   INDUCTION_T=mean((mag.bzt[good]+vel.ftvxx[good]+vel.ftvyy[good])^2,/double)
   if (finite(induction_T[0],/nan) eq 1) then stop
;
;;;;;;;;;;;;;;      Total Velocity       ;;;;;;;;;;;;;; 
   VX_R=(R_correlate(truth.Vx[good],vel.U0[good]))[0]
   VY_R=(R_correlate(truth.Vy[good],vel.V0[good]))[0]
   VZ_R=(R_correlate(truth.Vz[good],vel.W0[good]))[0]
   VB_R=(R_correlate(truth.VB[good],vel.VB[good]))[0]
;
   VX_C=correlate(truth.Vx[good],vel.U0[good])
   VY_C=correlate(truth.Vy[good],vel.V0[good])
   Vz_C=correlate(truth.Vz[good],vel.W0[good])
   VB_C=correlate(truth.VB[good],vel.VB[good])
;
   VX_L=ladfit(truth.Vx[good],vel.U0[good])
   VY_L=ladfit(truth.Vy[good],vel.V0[good])
   Vz_L=ladfit(truth.Vz[good],vel.W0[good])
   VB_L=ladfit(truth.VB[good],vel.VB[good])
   
;
;;;;;;;;;;;;;;  Perpendicular Velocity   ;;;;;;;;;;;;;; 
   VPX_R=(R_correlate(truth.VPX[good],vel.VPX[good]))[0]
   VPY_R=(R_correlate(truth.VPY[good],vel.VPY[good]))[0]
   VPZ_R=(R_correlate(truth.VPZ[good],vel.VPZ[good]))[0]
   VPX_C=correlate(truth.VPX[good],vel.VPX[good])
   VPY_C=correlate(truth.VPY[good],vel.VPY[good])
   VPz_C=correlate(truth.VPZ[good],vel.VPZ[good])
   VPX_L=ladfit(truth.VPX[good],vel.VPX[good])
   VPY_L=ladfit(truth.VPY[good],vel.VPY[good])
   VPz_L=ladfit(truth.VPZ[good],vel.VPZ[good])
;
   VP=sqrt(truth.VPX[good]^2+truth.VPY[good]^2+truth.VPZ[good]^2)
   VELP=sqrt(vel.VPX[good]^2+vel.VPY[good]^2+vel.VPZ[good]^2)
   metric=sqrt((truth.VPX[good]-vel.VPX[good])^2$
                   +(truth.VPY[good]-vel.VPY[good])^2$
                   +(truth.VPZ[good]-vel.VPZ[good])^2)
;
;      compute fractional error
   DVP=moment(metric/VP,/double)
;
;      fractional error in magnitude
   metric=(VELP-VP)/VP
   DVP_SPEED=moment(metric,/double)
;
;      vector correlation
   DOT=truth.VPX[good]*vel.VPX[good]+truth.VPY[good]*vel.VPY[good]+$
                                    truth.VPZ[good]*vel.VPZ[good]
   metric=DOT/sqrt(mean(VP^2,/double)*mean(VELP^2,/double))
   C_VP=moment(metric,/double)
;
;      direction correlation
   metric=DOT/(VP*VELP)
   CS_VP=moment(metric,/double)
;
;      weighted cosine
   WCOS_VP=total(metric*FTVEL/norm,/double)
;    
;;;;;;;;;;;;;;      Electric Field       ;;;;;;;;;;;;;; 
   EX_R=(R_correlate(truth.EX[good],vel.EX[good]))[0]
   EX_C=correlate(truth.EX[good],vel.EX[good])
   EX_L=ladfit(truth.EX[good],vel.EX[good])
   EY_R=(R_correlate(truth.EY[good],vel.EY[good]))[0]
   EY_C=correlate(truth.EY[good],vel.EY[good])
   EY_L=ladfit(truth.EY[good],vel.EY[good])
   EZ_R=(R_correlate(truth.EZ[good],vel.EZ[good]))[0]
   Ez_C=correlate(truth.EZ[good],vel.EZ[good])
   Ez_L=ladfit(truth.EZ[good],vel.EZ[good])
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   results={NUMBER_OF_POINTS:NG,$
            FTVX_R:FTVX_R[0],FTVY_R:FTVY_R[0],$
            FTVX_C:FTVX_C,FTVY_C:FTVY_C,$
            FTVX_L:FTVX_L,FTVY_L:FTVY_L,$
            DFTV:DFTV,DFTV_SPEED:DFTV_SPEED,$
            C_FTV:C_FTV,CS_FTV:CS_FTV,WCOS_FTV:WCOS_FTV,$
            FTV_THETA:FTV_THETA,$
            POYNTING_R:POYNTING_R,POYNTING_C:POYNTING_C,$
            POYNTING_RAT:POYNTING_RAT,$
            POYNTING_ANMHD:POYNTING_ANMHD,$
            POYNTING_L:POYNTING_L,$
            POYNTING_T:POYNTING_T,$
            POYNTING_D:POYNTING_D,$
            HELICITY_R:HELICITY_R,HELICITY_C:HELICITY_C,$
            HELICITY_RAT:HELICITY_RAT,$
            HELICITY_ANMHD:HELICITY_ANMHD,$
            HELICITY_L:HELICITY_L,$
            HELICITY_T:HELICITY_T,$
            HELICITY_D:HELICITY_D,$
            INDUCTION_R:INDUCTION_R,$
            INDUCTION_C:INDUCTION_C,$
            INDUCTION_L:INDUCTION_L,$
            INDUCTION_T:INDUCTION_T,$
            VPX_R:VPX_R,VPX_C:VPX_C,VPX_L:VPX_L,$
            VPY_R:VPY_R,VPY_C:VPY_C,VPY_L:VPY_L,$
            VPZ_R:VPZ_R,VPZ_C:VPZ_C,VPZ_L:VPZ_L,$
            DVP:DVP,DVP_SPEED:DVP_SPEED,C_VP:C_VP,$
            CS_VP:CS_VP,WCOS_VP:WCOS_VP,$
            EX_R:EX_R,EX_C:EX_C,EX_L:EX_L,$
            EY_R:EY_R,EY_C:EY_C,EY_L:EY_L,$
            EZ_R:EZ_R,EZ_C:EZ_C,EZ_L:EZ_L,$
            VX_R:VX_R,VX_C:VX_C,VX_L:VX_L,$
            VY_R:VY_R,VY_C:VY_C,VY_L:VY_L,$
            VZ_R:VZ_R,VZ_C:VZ_C,VZ_L:VZ_L,$
            VB_R:VB_R,VB_C:VB_C,VB_L:VB_L}
;
return,results

end






