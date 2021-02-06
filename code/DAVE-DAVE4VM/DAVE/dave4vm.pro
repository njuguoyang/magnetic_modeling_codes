function dave4vm,mag,window_size,$
                    sv=sv,sigma2=sigma2,chi2=chi2,errors=errors,$
                    threshold=threshold,noise=noise,double=double,AM=AM,$
                    float=float,missing_value=missing_value,help=help,$
                    kernel=kernel
;
; on_error,2 ; Return to the caller of the program unit that established the ON_ERROR condition.
;
 nparms=N_Params()
;
 if check_keyword(help) then begin
info:
    Message,$
    '============================================================================'$
     ,/info  
    Message,'',/info
    Message,'Usage:',/info 
    Message,'IDL>vel=dave4vm(mag,window_size,$',/info
    Message,'IDL>sv=sv,sigma2=sigma2,chi2=chi2,errors=errors,threshold=threshold,$',/info
    Message,'IDL>noise=noise,double=double,missing_value=missing_value)',/info
    Message,'',/info
    Message,$
    'Title: DAVE4VM - Differential Affine Velocity Estimator',/info 
    Message,$
    '                                                   for Vector Magnetograms',/info
    Message,'',/info
    Message,'Computes plasma velocities from a structure of magnetic field measurements',/info
    Message,'',/info
    Message,'Please reference:',/info
    Message,'Schuck, P. W., Tracking Vector Magnetograms with the Magnetic Induction',/info
    Message,'              Equation, Submitted to ApJ, 2008.',/info
    Message,'',/info
    Message,'Author: Peter W. Schuck',/info  
    Message,'schuck@ppdmail.nrl.navy.mil',/info  
    Message,'Plasma Physics Division',/info  
    Message,'United States Naval Research Laboratory',/info  
    Message,'Washington, DC, 20375',/info  
    Message,'',/info
    message,'VERSION HISTORY',/info
    Message,'VERSION 1.0 written: 07-11-2006 "dave_vm"',/info  
    Message,'VERSION 2.0 written: 01-15-2008 "dave_vm"',/info  
    Message,'VERSION 2.1 written: 01-31-2008 "dave4vm"',/info  
    Message,'',/info
    Message,'INPUT:',/info
    Message,'        MAG - structure of vector magnetic field measurements',/info
    Message,'              MAG.DX    X spatial scale (used to compute B?X)',/info
    Message,'              MAG.DY    Y spatial scale (used to compute B?Y)',/info
    Message,'                        ? = X Y Z',/info
    Message,'                                                           ',/info
    Message,'              MAG.BZT   Array[NX,NY]  time derivative of Bz',/info
    Message,'              MAG.BX    Array[NX,NY]  X component of B (Bx)',/info
    Message,'              MAG.BXX   Array[NX,NY]  X derivative of Bx',/info
    Message,'              MAG.BXY   Array[NX,NY]  Y derivative of By',/info
    Message,'              MAG.BY    Array[NX,NY]  Y component of B (By)',/info
    Message,'              MAG.BYX   Array[NX,NY]  X derivative of By',/info
    Message,'              MAG.BYY   Array[NX,NY]  Y derivative of By',/info
    Message,'              MAG.BZ    Array[NX,NY]  Z component of B (Bz)',/info
    Message,'              MAG.BZX   Array[NX,NY]  X derivative of Bz',/info
    Message,'              MAG.BZY   Array[NX,NY]  Y derivative of Bz',/info
    Message,'',/info
    Message,'WINDOW_SIZE - A one or two element vector for the window aperture',/info
    Message,'',/info
    Message,'OUTPUT:',/info
    Message,'        VEL - Array[NX,NY] of structures of coefficients',/info
    Message,'              U0        X-Velocity ',/info
    Message,'              V0        Y-Velocity ',/info
    Message,'              W0        Z-Velocity ',/info
    Message,'              UX        Local X derivative of the X-Velocity',/info
    Message,'              VX        Local X derivative of the Y-Velocity',/info
    Message,'              WX        Local X derivative of the Z-Velocity',/info
    Message,'              UY        Local Y derivative of the X-Velocity',/info
    Message,'              VY        Local Y derivative of the Y-Velocity',/info
    Message,'              WY        Local Y derivative of the Z-Velocity',/info
    Message,'              WINDOW_SIZE Local window size',/info
    Message,'',/info
    dave_keywords
    Message,'',/info
    Message,'***************************************************************',/info
    Message,'',/info
    Message,'Important! Velocities must be orthogonalized to obtain plasma velocities ',/info
    Message,'perpendicular to the magnetic field!',/info
    Message,'',/info
    Message,'***************************************************************',/info
    Message,'',/info
    Message,'AUTHORIZATION TO USE AND DISTRIBUTE',/info
    Message,'I hereby agree to the following terms governing the use and',/info
    Message,'redistribution of the DAVE4VM software release originally',/info
    Message,'written and developed by Dr. P. W. Schuck',/info

    Message,'',/info
    Message,'Redistribution and use in source and binary forms, with or',/info 
    Message,'without modification, are permitted provided that (1) source',/info 
    Message,'code distributions retain this paragraph in its entirety, (2)',/info 
    Message,'distributions including binary code include this paragraph in',/info 
    Message,'its entirety in the documentation or other materials provided',/info 
    Message,'with the distribution, (3) improvements, additions and',/info 
    Message,'upgrades to the software will be provided to NRL Authors in',/info 
    Message,'computer readable form, with an unlimited, royalty-free',/info 
    Message,'license to use these improvements, additions and upgrades',/info 
    Message,'and the authority to grant unlimited royalty-free sublicenses',/info
    Message,'to these improvements and (4) all published research using',/info
    Message,'this software display the following acknowledgment',/info
    Message,'``This work uses the DAVE/DAVE4VM codes written and developed',/info
    Message,'by the Naval Research Laboratory.''',/info
    Message,'',/info
    Message,'Neither the name of NRL or its contributors, nor any entity',/info
    Message,'of the United States Government may be used to endorse or',/info
    Message,'promote products derived from this software, nor does the',/info 
    Message,'inclusion of the NRL written and developed software directly',/info
    Message,"or indirectly suggest NRL's or the United States Government's",/info
    Message,'endorsement of this product.',/info
    Message,'',/info
    Message,'THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS',/info
    Message,'OR IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE',/info
    Message,'IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A',/info
    Message,'PARTICULAR PURPOSE.',/info


    Message,'',/info
    Message,'***************************************************************',/info
    Message,$
   '============================================================================',/info
    return,-1
 endif
;
 if (nparms ne 2) then begin
    Message,'Incorrect number of parameters',/info
    message,'for help use:',/info
    message,'IDL> dave4vm,/help'
 endif
;
;
;     check keywords
   if check_keyword(errors,/knot) then errors=0
;
   if check_keyword(noise,/knot) then noise=0.d0
;
   if check_keyword(double,/knot) then double=0 else double=1
   if check_keyword(float,/knot) then float=0 else float=1
;
   if (float eq double) then begin
      if (float eq 0) then begin
         float=(size(mag.bz,/type) eq 4)
         double=(size(mag.bz,/type) eq 5)
      endif else begin
         Message,'both keywords "DOUBLE" and "FLOAT" cannot be set',/info
         message,'for help use:',/info
         message,'IDL> dave4vm,/help'
      endelse
   endif
;
   if check_keyword(missing_value,/knot) then  $  
      if (double) then missing_value=!values.D_NAN $
                  else missing_value=!values.F_NAN
;
   if  check_keyword(threshold,/knot) then threshold=1.d0
;
   NWS=N_ELEMENTS(window_size)
   if ((NWS lt 1) or (nws gt 2)) then begin
      Message,'WINDOW_SIZE must be one or two dimensional',/info
      message,'for help use:',/info
      message,'IDL> dave4vm,/help'
   endif
;
   MAG_LIST=['BX','BY','BZ','BXX','BYX','BZX','BXY','BYY','BZY','BZT']
   if (check_tags(mag_list,mag,/knot,/verbose)) then goto,info
;
;    define arrays
   sz=size(mag.Bz)
   if (sz[0] ne 2) then begin
      message,'Image arrays must be two dimensional',/info
      message,'for help use:',/info
      message,'IDL> dave4vm,/help'
   endif
;
;      define arrays
   vel=missing_value
   dum=replicate(vel,sz[1],sz[2])
   v={U0:dum,UX:dum,UY:dum,$
      V0:dum,VX:dum,VY:dum,$
      W0:dum,WX:dum,WY:dum}
   
   SV=make_array(9,sz[1],sz[2],double=double,float=float)
   WW=make_array(9,9,double=double,float=float)
   chi2=make_array(sz[1],sz[2],double=double,float=float)
   sigma2=make_array(9,9,sz[1],sz[2],double=double,float=float)
;
   chi2[*]=missing_value
   sigma2[*]=missing_value
   aperture=replicate(255B,sz[1],sz[2]) 
;
   id=indgen(9)
;
;  *** NOTE IMPORTANT! ********
;  it is necessary for the integer (pixel) values of X and Y to
;  be multiplied by dx and dy respectively to ensure rescaling
;  of derivatives leads to the identical answer
;
;    construct weighting functions
   nw=fix(2*fix(window_size/2)+1)
   if (n_elements(nw) eq 1) then nw=[nw[0],nw[0]]
   x=rebin((lindgen(nw[0])-nw[0]/2),nw[0],nw[1])*mag.dx
   y=transpose(rebin((lindgen(nw[1])-nw[1]/2),nw[1],nw[0]))*mag.dy
;
   psf=make_array(nw[0],nw[1],double=double,float=float)
;     default window is top-hat
   psf[*]=1.d0
;     normalize
   psf[*]=psf/total(psf,/double)
;     moments
   psfx=psf*x
   psfy=psf*y
   psfxx=psf*x^2
   psfyy=psf*y^2
   psfxy=psf*x*y
;
   v=add_tag(v,nw,'WINDOW_SIZE')
;
   AM=dave4vm_matrix(mag.Bx,mag.Bxx,mag.Bxy,$
                    mag.By,mag.Byx,mag.Byy,$
                    mag.Bz,mag.Bzx,mag.Bzy,$
                    mag.Bzt,psf,psfx,psfy,psfxx,psfyy,psfxy,$
                    double=double,float=float)
;
   kernel={psf:psf,psfx:psfx,psfy:psfy,psfxx:psfxx,psfyy:psfyy,psfxy:psfxy}
;
;    estimate trace
   trc=total((reform(AM,100,sz[1],sz[2]))[id*10+id,*,*],1,/double)
;    
;    find locations where the aperture problem could be resolved 
   index=where(trc gt threshold,N) 
;
   if (N eq 0) then begin
     message,'The input images MAG.Bx are pathological or the window',/info
     message,'size SIGMA is too small. The aperture problem cannot be',/info 
;             'resolved.',level=-1
     return,-1
   endif
;

;    loop over good pixels
   for ii=0L,N-1L do begin
       j=index[ii]/sz[1]
       i=index[ii] mod sz[1]
;
       AA=AM[*,*,i,j]
       GA=AA[0:8,0:8]
       FA=-reform(AA[0:8,9],9)
       DESIGN=GA
       SOURCE=FA
;
;       SVDC,DESIGN,W,UA,VA
;      SV[*,i,j]=W[sort(W)]
;       mn=0.d0
;       zero=where(W le mn)
;       if (zero[0] ne -1) then W[zero]=0.d0
;       VECTOR=SVSOL(UA,W,VA,SOURCE,/double)
;
;         USE divide-and-conquer algorithm 
;
;       rcondition=1.d-6
       vector=LA_LEAST_SQUARES(DESIGN,SOURCE,double=double,$
              method=3,residual=chisq,rcondition=rcondition)
;
       v.U0[i,j]=VECTOR[0]
       v.V0[i,j]=VECTOR[1]
       v.UX[i,j]=VECTOR[2]
       v.VY[i,j]=VECTOR[3]
       v.UY[i,j]=VECTOR[4]
       v.VX[i,j]=VECTOR[5]
       v.W0[i,j]=VECTOR[6]
       v.WX[i,j]=VECTOR[7]
       v.WY[i,j]=VECTOR[8]
;
    endfor
;
return,v
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
