pro old_dave,f,dft,dfx,dfy,dx,dy,window_size,vel,aperture,advect=advect,sv=sv,sigma2=sigma2,chi2=chi2,errors=errors,threshold=threshold,noise=noise,double=double,missing_value=missing_value
;
 on_error,2 ; Return to the caller of the program unit that established the ON_ERROR condition.
;
 nparms=N_Params()
 if (nparms ne 9) then begin
;
  Message,$
   '====================================================================='$
   ,/info  
  Message,'',/info
  Message,'Usage:',/info 
  Message,'IDL>dave,dft,dfx,dfy,dx,dy,window_size,vel,aperture,advect=advect,$',/info
  Message,'IDL>sv=sv,sigma2=sigma2,chi2=chi2,errors=errors,threshold=threshold,$',/info
  Message,'IDL>noise=noise,double=double,missing_value=missing_value',/info
  Message,'',/info
  Message,$
  'Title: DAVE - Differential Affine Velocity Estimator',/info 
  Message,'',/info
  Message,'Author: Peter W. Schuck',/info  
  Message,'schuck@ppdmail.nrl.navy.mil',/info  
  Message,'Plasma Physics Division',/info  
  Message,'United States Naval Research Laboratory',/info  
  Message,'Washington, DC, 20375',/info  
  Message,'',/info
  Message,'VERSION 1.0 written: 05-01-2005',/info  
  Message,'VERSION 1.1 written: 12-06-2005',/info  
  Message,'           released: 01-01-2006',/info  
  Message,'',/info
  Message,$
   '====================================================================='$
   ,/info  
  Message,'',/info
  Message,"Estimates the optical flow in the center image 'f' from an",/info
  Message,'ordered sequence of three images: f1, f, f3.',/info
  Message,'',/info
  Message,'Language: IDL v. 5.6',/info
  Message,'',/info
  Message,'Code tested under the following compilers/operating systems:',/info
  Message,'          SUSE 9.0, 10.0',/info
  Message,'',/info
  Message,'Description of input data:',/info
  Message,'',/info
  Message,'   F1, F, and F3 - three same-size floating-point 2D images [NX,NY]',/info
  Message,'',/info
  Message,'        DX and DY - The floating point pixel separation scales',/info
  Message,'                    (in physical units).',/info
  Message,'',/info
  Message,'               DT - The floating point time separation',/info
  Message,'                    (in physical units).',/info
  Message,'',/info
  Message,'      WINDOW_SIZE - Window size (in pixels) that defines a ',/info
  Message,'                    neighborhood for the optical  flow estimate',/info
  Message,'                  * Can be anisotropic, i.e. a two element vector.',/info
  Message,'',/info
  Message,'Optional Inputs:',/info
  Message,'',/info
  Message,'           NOISE  - Noise threshold for minimum eigenvalue of',/info
  Message,'                    the image structure tensor (default=1.e-7).',/info
  Message,'                    lambda_i/lambda_6 must be greater than NOISE or',/info
  Message,'                    lambda_i is assumed to be zero.',/info
  Message,'',/info
  Message,'         THRESHOLD - Rough threshold for resolving the aperture problem',/info
  Message,'                    (default=1.). The trace of the structure tensor',/info
  Message,'                    must be larger than THRESH or the DAVE assumes',/info
  Message,'                    that the aperture problem cannot be resolved at',/info
  Message,'                    this location.',/info
  Message,'',/info
  Message,'    MISSING_VALUE - Value to replace any missing/invalid data',/info
  Message,'                    (default=!Values.F_NAN if DOUBLE=0 or',/info 
  Message,'                    !Values.D_NAN if DOUBLE=1).',/info
  Message,'',/info
  Message,'           ERRORS - If set, then compute velocity uncertainties',/info
  Message,'                    returned in CHI2 and SIGMA2.',/info
  Message,'                *   (these calculations take extra time)',/info
  Message,'',/info
  Message,'           ADVECT - If set, use the advection equation instead of the',/info
  Message,'                    continuity equation.',/info
  Message,'',/info
  Message,'           DOUBLE - If set, perform computations in double precision.',/info
  Message,'',/info
  Message,'Description of output data:',/info
  Message,'',/info
  Message,'     VEL[6,NX,NY] - A 3D floating point array of velocities that',/info
  Message,'                    optimally satisfies the window weighted',/info
  Message,'                    continuity (advect=0) or advection (advect=1)',/info
  Message,'                    operator.',/info
  Message,'     VEL[1,NX,NY] - U (X component of velocity)',/info
  Message,'     VEL[2,NX,NY] - V (Y component of velocity)',/info
  Message,'        shears',/info
  Message,'     VEL[3,NX,NY] - U_X (X gradient in the X component of velocity).',/info
  Message,'     VEL[4,NX,NY] - V_Y (Y gradient in the Y component of velocity).',/info
  Message,'     VEL[5,NX,NY] - U_Y (Y gradient in the X component of velocity).',/info
  Message,'     VEL[6,NX,NY] - V_X (X gradient in the Y component of velocity).',/info
  Message,'',/info
  Message,'  APERTURE[NX,NY] - Describes aperture resolution: unique/ambiguous',/info
  Message,'                    (0 le APERTURE[i,j] lt 128) unique velocity',/info
  Message,'                    (128 le APERTURE[i,j] lt 256) ambiguous velocity',/info
  Message,'                    Particular values of APERTURE[i,j] within these',/info
  Message,'                    two ranges may have specific meanings in future',/info
  Message,'                    versions of DAVE',/info
  Message,'',/info
  Message,'Optional Outputs:',/info
  Message,'',/info
  Message,'      SV[6,NX,NY] - A 3D floating point array of ordered singular',/info  
  Message,'                    values of the intensity structure tensor',/info
  Message,'                    SV[0,i,j]<SV[1,i,j]<....<SV[6,i,j].',/info
  Message,'',/info
  Message,'      CHI2[NX,NY] - Chi-squared value for the neighborhood.',/info
  Message,'',/info
  Message,'SIGMA2[6,6,NX,NY] - Velocity uncertainties squared.',/info
  Message,'',/info
  Message,'',/info
  Message,'System requirements:',/info
  Message,'',/info
  Message,'Calls to Solarsoft routines: PSF_GAUSSIAN',/info
  Message,'',/info
  Message,'Additional comments:',/info
  Message,'',/info
  Message,'If you use this routine, please reference:',/info
  Message,'',/info
  Message,'Schuck, "Tracking magnetic footpoints with the magnetic',/info
  Message,'         induction equation", Ap. J., 2005',/info
  Message,'',/info
  Message,$
   '====================================================================='
;
 endif
;
;     check keywords
   if keyword_set(advect) then nu=0 else nu=1
   if not keyword_set(errors) then errors=0
   if not keyword_set(noise) then noise=0.d0
   if not keyword_set(double) then double=0 else begin
      double=1
      f=double(f)
      dft=double(dft)
      dfx=double(dfx)
      dfy=double(dfy)
   endelse
   if not(keyword_set(missing_value)) then  $  
      if (double) then missing_value=!values.D_NAN else  missing_value=!values.F_NAN
;   print,threshold,ARG_PRESENT(threshold)
   if N_ELEMENTS(threshold) eq 0 then threshold=1.
;
   print,'nu=',nu
;
;    define arrays
   sz=size(dft)
   if (sz[0] ne 2) then begin
      message,'Image arrays dft,dfx,dfy must be two dimensional'
   endif
   vel=make_array(6,sz[1],sz[2],double=double,float=float)
   vel[*]=missing_value
   A=make_array(7,7,sz[1],sz[2],double=double,float=float)
   SV=make_array(6,sz[1],sz[2],double=double,float=float)
   WW=make_array(6,6,double=double,float=float)
   chi2=make_array(sz[1],sz[2],double=double,float=float)
   sigma2=make_array(6,6,sz[1],sz[2],double=double,float=float)
;
   chi2[*]=missing_value
   sigma2[*]=missing_value
   aperture=replicate(255B,sz[1],sz[2]) 
;
   id=indgen(6)
;
;    construct weighting functions
   nw=fix(2*fix(window_size/2)+1)
   if (n_elements(nw) eq 1) then nw=[nw[0],nw[0]]
   x=rebin((lindgen(nw[0])-nw[0]/2)*dx,nw[0],nw[1])
   y=transpose(rebin((lindgen(nw[1])-nw[1]/2)*dy,nw[1],nw[0]))
   psf=make_array(nw[0],nw[1],double=double,float=float)
;
   psf[*]=psf_gaussian( NPixel=nw,2,fwhm=window_size/4,/NORM )
;   kill=where(psf lt 0.5*max(psf))
;   print,'gaussian'
   psf=(hanning(nw[0]+1,nw[1]+1,double=double))[1:*,1:*]
;   stop
;   kill=where(psf lt 0.5*max(psf))
;   psf[kill]=0.
   psf[*]=1.d0
;   psf=psf
   psf[*]=psf/total(psf,double=double)*dx*dy*nw[0]*nw[1]
;
   psfx=psf*x
   psfy=psf*y
   psfxx=psf*x^2
   psfyy=psf*y^2
   psfxy=psf*x*y
;   stop
;
;    construct matrix elements for LKA algorithm
;
   GG=f*f
   G=reform(convol(GG,psf,/center),1,sz[1],sz[2])                       ; 1
;
   GGx=f*dfx
   Gx=reform(convol(GGx,psf,/center),1,sz[1],sz[2])                     ; 2
;
   xGx=reform(convol(GGx,psfx,/center),1,sz[1],sz[2])                   ; 3
;
   yGx=reform(convol(GGx,psfy,/center),1,sz[1],sz[2])                   ; 4
;
   GGy=f*dfy
   Gy=reform(convol(GGy,psf,/center),1,sz[1],sz[2])                     ; 5
;  
   xGy=reform(convol(GGy,psfx,/center),1,sz[1],sz[2])                   ; 6
;
   yGy=reform(convol(GGy,psfy,/center),1,sz[1],sz[2])                   ; 7

   GGt=dft*f
   Ht=reform(convol(GGt,psf,/center),1,sz[1],sz[2])                     ; 8
;   Ht[*]=0
;
   GGxx=dfx*dfx
   Gxx=reform(convol(GGxx,psf,/center),1,sz[1],sz[2])                   ; 9
;
   GGyy=dfy*dfy
   Gyy=reform(convol(GGyy,psf,/center),1,sz[1],sz[2])                   ; 10
;
   GGxy=dfx*dfy
   Gxy=reform(convol(GGxy,psf,/center),1,sz[1],sz[2])                   ; 11
;
   GGtx=dft*dfx
   Gtx=reform(convol(GGtx,psf,/center),1,sz[1],sz[2])                   ; 12
;
   GGty=dft*dfy
   Gty=reform(convol(GGty,psf,/center),1,sz[1],sz[2])                   ; 13
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
   xGxx=reform(convol(GGxx,psfx,/center),1,sz[1],sz[2])                 ; 14
;
   xGyy=reform(convol(GGyy,psfx,/center),1,sz[1],sz[2])                 ; 15
;
   xGxy=reform(convol(GGxy,psfx,/center),1,sz[1],sz[2])                 ; 16
;
   xGtx=reform(convol(GGtx,psfx,/center),1,sz[1],sz[2])                 ; 17
;
   xGty=reform(convol(GGty,psfx,/center),1,sz[1],sz[2])                 ; 18
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
   yGxx=reform(convol(GGxx,psfy,/center),1,sz[1],sz[2])                 ; 19
;
   yGyy=reform(convol(GGyy,psfy,/center),1,sz[1],sz[2])                 ; 20
;
   yGxy=reform(convol(GGxy,psfy,/center),1,sz[1],sz[2])                 ; 21
;
   yGtx=reform(convol(GGtx,psfy,/center),1,sz[1],sz[2])                 ; 22
;
   yGty=reform(convol(GGty,psfy,/center),1,sz[1],sz[2])                 ; 23
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
   xxGxx=reform(convol(GGxx,psfxx,/center),1,sz[1],sz[2])               ; 24
;
   xxGxy=reform(convol(GGxy,psfxx,/center),1,sz[1],sz[2])               ; 25
;
   xxGyy=reform(convol(GGyy,psfxx,/center),1,sz[1],sz[2])               ; 26
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
   xyGxx=reform(convol(GGxx,psfxy,/center),1,sz[1],sz[2])               ; 27
;
   xyGyy=reform(convol(GGyy,psfxy,/center),1,sz[1],sz[2])               ; 28
;
   xyGxy=reform(convol(GGxy,psfxy,/center),1,sz[1],sz[2])               ; 29
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
   yyGxx=reform(convol(GGxx,psfyy,/center),1,sz[1],sz[2])               ; 30
;
   yyGxy=reform(convol(GGxy,psfyy,/center),1,sz[1],sz[2])               ; 31
;
   yyGyy=reform(convol(GGyy,psfyy,/center),1,sz[1],sz[2])               ; 32
;
   GGtt=dft*dft
   Gtt=reform(convol(GGtt,psf,/center),1,sz[1],sz[2])                   ; 33
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   
;   
      A[*]=reform([Gxx,Gxy,nu*Gx+xGxx,nu*Gx+yGxy,yGxx,xGxy,Gtx,$
          Gxy,Gyy,nu*Gy+xGxy,nu*Gy+yGyy,yGxy,xGyy,Gty,$
          nu*Gx+xGxx,nu*Gy+xGxy,nu*(G*nu+2*xGx)+xxGxx,nu*(G*nu+xGx+yGy)+xyGxy,nu*yGx+xyGxx,nu*xGy+xxGxy,nu*Ht+xGtx,$   
          nu*Gx+yGxy,nu*Gy+yGyy,nu*(G*nu+xGx+yGy)+xyGxy,nu*(G*nu+2*yGy)+yyGyy,nu*yGx+yyGxy,nu*xGy+xyGyy,nu*Ht+yGty,$
          yGxx,yGxy,nu*yGx+xyGxx,nu*yGx+yyGxy,yyGxx,xyGxy,yGtx,$
          xGxy,xGyy,nu*xGy+xxGxy,nu*xGy+xyGyy,xyGxy,xxGyy,xGty,$
          Gtx,Gty,nu*Ht+xGtx,nu*Ht+yGty,yGtx,xGty,Gtt],7,7,sz[1],sz[2])
;
;
;    estimate trace
    trc=total((reform(A,49,sz[1],sz[2]))[id*7+id,*,*],1,/double)
;    
;    find locations where the aperture problem could be resolved 
   index=where(trc gt threshold,N) 
;
   if (N eq 0) then begin
     message,'The input images f1,f,f3 are pathological or the window',/info
     message,'size SIGMA is too small. The aperture problem cannot be',/info 
;             'resolved.',level=-1
     return
   endif
;
;    loop over good pixels
   for ii=0L,N-1L do begin
     j=index[ii]/sz[1]
     i=index[ii] mod sz[1]
;
     AA=A[*,*,i,j]
     GA=AA[0:5,0:5]
     FA=-reform(AA[0:5,6],6)
     SVDC,GA,WA,UA,VA ,double=double
     is=sort(WA)
     WA=WA(is)
     UA(id,*)=UA(is,*)
     VA(id,*)=VA(is,*)
     SV[*,i,j]=WA
     mn=min(WA,max=mx,imn)
;
      xm=total(psfx,double=doublw)
      ym=total(psfy,double=double)
      xym=total(psfxy,double=double)
      x2m=total(psfxx,double=double)
      y2m=total(psfyy,double=double)
;      TT=[[1., 0, xm, 0, ym, 0, 0], $
;          [0, 1., 0, ym, 0, xm, 0], $
;          [xm, 0, 1. + x2m, 1., xym, 0, 0],$ 
;          [0, ym, 1., 1. + y2m, 0, xym, 0], $
;          [ym, 0, xym, 0, y2m, 0, 0],$ 
;          [0, xm, 0, xym, 0, x2m, 0],$ 
;          [0, 0, 0, 0, 0, 0, 1.]]

      TT=[[1, 0, 0, 0, 0, 0, 0], $
          [0, 1, 0, 0, 0, 0, 0], $
          [0, 0, x2m,0.,0, 0, 0],$ 
          [0, 0, 0.,y2m, 0, 0, 0], $
          [0, 0, 0, 0, y2m, 0, 0],$ 
          [0, 0, 0,0, 0, x2m, 0],$ 
          [0, 0, 0, 0, 0, 0, 1]]
;
    if (mx gt 0) then begin
       WA=WA*(WA/WA[5] gt noise)
       if (WA[0] gt noise) then aperture[i,j]=0B else aperture[i,j]=128B
       vel[*,i,j]=SVSOL(UA, WA, VA, FA,double=double)
;       if (abs(f[i,j]) gt 600) then stop 
;       SWA=LA_EIGENQL(AA,TT,eigenvectors=SVA,generalized=0,method=2,/double)
;       vel[*,i,j]=float(SVA[0:5,0])/SVA[6,0]
;      if (WA[imn] gt 0.) then begin
;         WW[id,id]=1/WA[id]^2
;         if errors then begin
;;            compute chisq
;           ix1=i-nw[0]/2
;           ix2=i+nw[0]/2
;           iy1=j-nw[1]/2
;           iy2=j+nw[1]/2
;           vx=vel[0,i,j]+vel[2,i,j]*x+vel[4,i,j]*y
;           vy=vel[1,i,j]+vel[3,i,j]*y+vel[5,i,j]*x
;           chi2[i,j]=total(psf*(dfx[ix1:ix2,iy1:iy2]*vx+dfy[ix1:ix2,iy1:iy2]*vy+nu*(vel[2,i,j]+vel[3,i,j])*f[ix1:ix2,iy1:iy2]+dft[ix1:ix2,iy1:iy2])^2,/double)      
;           sigma2[*,*,i,j]=chi2[i,j]*VA##WW##transpose(VA)
;
;          note: GA-UA##[[WA[0],0],[0,WA[1]]]##Transpose(VA)~=~0
;
;         endif else sigma2[*,*,i,j]=VA##WW##transpose(VA)
;       endif else sigma2[*,*,i,j]=missing_value
     endif     
 endfor
; stop
;
return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
