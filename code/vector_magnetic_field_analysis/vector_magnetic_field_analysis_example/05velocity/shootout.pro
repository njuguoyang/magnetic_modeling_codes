function shootout,file1,file2,thr=thr,verbose=verbose,no_elliptic2d=no_elliptic2d,_extra=extra
;
   common constants,speed_of_light
   if (N_ELEMENTS(speed_of_light) eq 0) then speed_of_light=2.998d10 ; cm/s
   if (N_ELEMENTS(no_elliptic2d) eq 0) then no_elliptic2d=0
   if (N_ELEMENTS(thr) eq 0) then thr=370.d0 else thr=double(thr)
   solar_r=6.955d5 ;solar radius in km
;
   restore,file1,verbose=verbose,_extra=extra
   temp=pb0r(smapbz.time)
   DX=double(smapbz.dx*solar_r/temp[2]/60.0)  ;in kilometers
   DY=DX
   bz_start=smapbz.data
   bx_start=smapbx.data
   by_start=smapby.data
   internal_time=anytim2utc(smapbz.time)      ;we do not have to consider about the leap seconde here. If the 
                                              ;time ranges span over one day, the leap second maybe have to be considered in some cases.
   t_start=double(internal_time.time/1000.0)  ;in seconds
   restore,file2,verbose=verbose,_extra=extra
   bz_stop=smapbz.data
   bx_stop=smapbx.data
   by_stop=smapby.data
   internal_time=anytim2utc(smapbz.time)  
   t_stop=double(internal_time.time/1000.0)
;
;
   DT=double(T_stop-T_start)
   BZT=(BZ_stop-BZ_start)/DT
   sz=size(BZT)
;
;     use time_centered spatial variables
   BX=(BX_stop+BX_start)/2.d0
   BY=(BY_stop+BY_start)/2.d0
   BZ=(BZ_stop+BZ_start)/2.d0
;   
;     compute 5-point optimized derivatives
   odiffxy5,BX,BXX,BXY
   odiffxy5,BY,BYX,BYY
   odiffxy5,BZ,BZX,BZY
;
;     setup for mudpack solution of vector potential in main region
   N1=sz[1] & N2=sz[2]
   source=fltarr(N1,N2)
   xr=float([1,n1]*dx)
   yr=float([1,n2]*dy)
;
   ixp=2
   iex=9
   jyq=2
   jey=9
   nx1=ixp*(2^(iex-1))+1
   ny2=jyq*(2^(jey-1))+1
;
   mudpack={nx:long(nx1),ny:long(ny2),ixp:long(ixp),iex:long(iex),jyq:long(jyq),jey:long(jey),iguess:long(0),maxcy:long(5),method:long(3),loud:long(0),boundary:long(1),dx:dx,dy:dy,xr:xr,yr:yr}
;
   osx=12
   osy=12
;
   if (no_elliptic2d eq 0) then begin
;
   Bz1=dblarr(nx1,ny2)
   Bz1[osx:n1+osx-1,osy:n2+osy-1]=Bz
   Bzt1=dblarr(nx1,ny2)
   Bzt1[osx:n1+osx-1,osy:n2+osy-1]=Bzt
;
      source=Bz1
      compute_vector_potential,mudpack,source,phi
      mwrfits,phi,'phi.fts',/create   ; save for use without elliptic2d
;
      source=Bzt1
      compute_vector_potential,mudpack,source,phit
      mwrfits,phit,'phit.fts',/create ; save for use without elliptic2d
;
   endif else begin
;
      phi=mrdfits('phi.fts',1)
      phit=mrdfits('phit.fts',1)
;
   endelse
;
;   chi=poisson(BZ,dx)
;   dchidx = xderiv(chi, dx)        ; and grad phi
;   dchidy = yderiv(chi, dy)
;
   MAG={BZT:BZT,BX:BX,BXX:BXX/DX,BXY:BXY/DY,BY:BY,BYX:BYX/DX,BYY:BYY/DY,$
        BZ:BZ,BZX:BZX/DX,BZY:BZY/DY,DX:DX,DY:DY,DT:DT,THR:THR}
;
;     vector potential
   phiy=phi.y[osx:n1+osx-1,osy:n2+osy-1]
   phix=phi.x[osx:n1+osx-1,osy:n2+osy-1]
   phity=phit.y[osx:n1+osx-1,osy:n2+osy-1]
   phitx=phit.x[osx:n1+osx-1,osy:n2+osy-1]
   AX=-(phiy+shift(phiy,0,1))/2
   AY=(phix+shift(phix,1,0))/2
   ATX=-(phity+shift(phity,0,1))/2
   ATY=(phitx+shift(phitx,1,0))/2
   mag=add_tag(mag,AX,"AX")     ; [G-km]
   mag=add_tag(mag,AY,"AY")     ; [G-km]
   mag=add_tag(mag,ATX,"ATX")   ; [G-km/s]
   mag=add_tag(mag,ATY,"ATY")   ; [G-km/s]
;
;     unit vector (direction of the magnetic field)
   B=sqrt(bx^2+by^2+bz^2)
   qx=bx/B
   qy=by/B
   qz=bz/B
;
   mag=add_tag(mag,qx,'QX')
   mag=add_tag(mag,qy,'QY')
   mag=add_tag(mag,qz,'QZ')
;
   return,mag

end
