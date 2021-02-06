 function dave,mag,window_size,advect=advect,$
                   sv=sv,sigma2=sigma2,chi2=chi2,errors=errors,$
                   threshold=threshold,noise=noise,double=double,AM=AM,$
                   float=float,missing_value=missing_value,help=help,$
                   junk=junk
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
    Message,'IDL>vel=dave(mag,window_size,$',/info
    Message,'IDL>sv=sv,sigma2=sigma2,chi2=chi2,errors=errors,threshold=threshold,$',/info
    Message,'IDL>noise=noise,double=double,missing_value=missing_value)',/info
    Message,'',/info
    Message,$
    'Title: DAVE - Differential Affine Velocity Estimator',/info 
    Message,$
    '                                                   for Vector Magnetograms',/info
    Message,'',/info
    Message,'Computes plasma velocities from a structure of magnetic field measurements',/info
    Message,'',/info
    Message,'Please reference:',/info
    Message,'Schuck, P. W., Tracking magnetic footpoints with the magnetic',/info
    Message,'         induction equation", ApJ, 646, 1385, 2006',/info
    Message,'Schuck, P. W., Local correlation tracking and the magnetic',/info
    Message,'         induction equation", ApJ, 632, 53, 2005',/info
    Message,'',/info
    Message,'Author: Peter W. Schuck',/info  
    Message,'schuck@ppdmail.nrl.navy.mil',/info  
    Message,'Plasma Physics Division',/info  
    Message,'United States Naval Research Laboratory',/info  
    Message,'Washington, DC, 20375',/info  
    Message,'',/info
    message,'VERSION HISTORY',/info
    Message,'VERSION 1.0 written: 05-01-2005 "dave"',/info  
    Message,'VERSION 1.1 written: 12-06-2005 "dave"',/info  
    Message,'        1.2 released: 02-15-2008 "dave"',/info  
    Message,'',/info
    Message,'INPUT:',/info
    Message,'        MAG - structure of vector magnetic field measurements',/info
    Message,'              MAG.DX    X spatial scale (used to compute B?X)',/info
    Message,'              MAG.DY    Y spatial scale (used to compute B?Y)',/info
    Message,'                        ',/info
    Message,'              MAG.BZT   Array[NX,NY]  time derivative of Bz',/info
    Message,'              MAG.BZ    Array[NX,NY]  Z component of B (Bz)',/info
    Message,'              MAG.BZX   Array[NX,NY]  X derivative of Bz',/info
    Message,'              MAG.BZY   Array[NX,NY]  Y derivative of Bz',/info
    Message,'',/info
    Message,'WINDOW_SIZE - A one or two element vector for the window aperture',/info
    Message,'',/info
    Message,'KEYWORDS: (input)',/info
    Message,'     ADVECT - If set, use the advection equation instead of the',/info 
    Message,'              continuity equation',/info
    Message,'',/info
    Message,'OUTPUT:',/info
    Message,'        VEL - Array[NX,NY] of structures of coefficents',/info
    Message,'              U0         X-Flux transport velocity ',/info
    Message,'              V0         Y-Flux transport velocity ',/info
    Message,'              UX         Local X derivative of the X-Flux transport velocity',/info
    Message,'              VX         Local X derivative of the Y-Flux transport velocity',/info
    Message,'              UY         Local Y derivative of the X-Flux transport velocity',/info
    Message,'              VY         Local Y derivative of the Y-Flux transport velocity',/info
    Message,'              WINDOW_SIZE Local window size',/info
    Message,'',/info
    dave_keywords
    Message,'',/info
    Message,'***************************************************************',/info
    Message,'',/info
    Message,'AUTHORIZATION TO USE AND DISTRIBUTE',/info
    Message,'I hereby agree to the following terms governing the use and',/info
    Message,'redistribution of the DAVE software release originally',/info
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
    message,'IDL> dave,/help'
 endif
;
;     check keywords
   if check_keyword(advect,/knot) then nu=1 else if (advect eq 1) then nu=0 else nu=1 
;
   print,'nu=',nu
   if check_keyword(errors,/knot) then errors=0L
;
   if check_keyword(noise,/knot) then noise=1.d-16
;
   if check_keyword(double,/knot) then double=0L else double=1L
   if check_keyword(float,/knot) then float=0L else float=1L
;
   if (float eq double) then begin
      if (float eq 0) then begin
         float=(size(mag.bz,/type) eq 4)
         double=(size(mag.bz,/type) eq 5)
      endif else begin
         Message,'both keywords "DOUBLE" and "FLOAT" cannot be set',/info
         message,'for help use:',/info
         message,'IDL> dave,/help'
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
      message,'IDL> dave,/help'
   endif
;
   MAG_LIST=['BZ','BZX','BZY','BZT']
   if (check_tags(mag_list,mag,/knot,/verbose)) then goto,info
;
;    define arrays
   sz=size(mag.Bz)
   if (sz[0] ne 2) then begin
      message,'Image arrays must be two dimensional',/info
      message,'for help use:',/info
      message,'IDL> dave,/help
   endif
;
   vel=replicate(missing_value,sz[1],sz[2])
   v={U0:VEL,V0:VEL,UX:VEL,VY:VEL,UY:VEL,VX:VEL,NU:NU}
;
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
;  *** NOTE IMPORTANT! ********
;  it is necessary of the integer (pixel) values of X and Y to
;  be multiplied by dx and dy respectively to ensure rescaleing
;  of derivatives leads to the identical answer
;
;    construct weighting functions
   nw=fix(2*fix(window_size/2)+1)
   if (n_elements(nw) eq 1) then nw=[nw[0],nw[0]]
   window_size=nw
;
   x=rebin((lindgen(nw[0])-nw[0]/2),nw[0],nw[1])*mag.dx
   y=transpose(rebin((lindgen(nw[1])-nw[1]/2),nw[1],nw[0]))*mag.dy
;
   psf=make_array(nw[0],nw[1],double=double,float=float)
;     default window is top-hat
;   psf[*]=(hanning(nw[0]+1,nw[1]+1,double=double))[1:*,1:*]
   psf[*]=1 
;     normalize
   psf[*]=psf/total(psf,double=double)
;     moments
   psfx=psf*x
   psfy=psf*y
   psfxx=psf*x^2
   psfyy=psf*y^2
   psfxy=psf*x*y
;
   v=add_tag(v,nw,'WINDOW_SIZE')
;
   AM=dave_matrix(mag.Bz,mag.Bzx,mag.Bzy,mag.Bzt,psf,psfx,psfy,psfxx,$
                                      psfyy,psfxy,nu,double=double,float=float)
;
;     estimate trace
   trc=total((reform(AM,49,sz[1],sz[2]))[id*7+id,*,*],1,/double)
;    
;    find locations where the aperture problem could be resolved 
   index=where(trc gt threshold,N) 
;
   if (N eq 0) then begin
     message,'The input images MAG.BZ are pathological or the window',/info
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
       GA=AA[0:5,0:5]
       FA=-reform(AA[0:5,6],6)
       DESIGN=GA
       SOURCE=FA
;
;         USE divide-and-conquer algorithm 
;
;       vector=LA_LEAST_SQUARES(DESIGN,SOURCE,double=double,$
;              method=2,residual=chisq)
       vector=LA_LINEAR_EQUATION(DESIGN,SOURCE,double=double,status=status)
       if (status gt 0) then goto,skip
       goto,use_la
;
       SVDC,GA,WA,UA,VA ,double=double
       is=sort(WA)
       WA=WA(is)
       UA(id,*)=UA(is,*)
       VA(id,*)=VA(is,*)
       SV[*,i,j]=WA
       mn=min(WA,max=mx,imn)
;
       if (mx gt 0) then begin
          WA=WA*(WA/WA[5] gt noise)
          if (WA[0] gt noise) then aperture[i,j]=0B else aperture[i,j]=128B
          VECTOR=SVSOL(UA, WA, VA, FA,double=double)
       endif
;
use_la: 
;
       v.U0[i,j]=VECTOR[0]
       v.V0[i,j]=VECTOR[1]
       v.UX[i,j]=VECTOR[2]
       v.VY[i,j]=VECTOR[3]
       v.UY[i,j]=VECTOR[4]
       v.VX[i,j]=VECTOR[5]
;
skip:
;
   endfor
;
return,v
;
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
