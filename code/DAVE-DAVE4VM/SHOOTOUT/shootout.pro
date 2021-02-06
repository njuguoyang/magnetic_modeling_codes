function shootout,file,istep,truth=truth,verbose=verbose,no_elliptic2d=no_elliptic2d,_extra=extra
;
   common constants,speed_of_light
   if (N_ELEMENTS(speed_of_light) eq 0) then speed_of_light=2.998d10 ; cm/s
   if (N_ELEMENTS(no_elliptic2d) eq 0) then no_elliptic2d=0
;
   restore,file,verbose=verbose,_extra=extra
   DX=double(DX)  ; (in kilometers)
   DY=DX
   thr=double(thr)
;
   if (N_ELEMENTS(istep) ne 1) then istep=0L ; default to first time-step
;
   DT=double((T_stop-T_start)[istep])
   BZT=(BZ_stop-BZ_start)[*,*,istep]/DT
   sz=size(BZT)
;
;     use time_centered spatial variables
   BX=(BX_stop+BX_start)[*,*,istep]/2.d0
   BY=(BY_stop+BY_start)[*,*,istep]/2.d0
   BZ=(BZ_stop+BZ_start)[*,*,istep]/2.d0
;
   VX=(VX_stop+VX_start)[*,*,istep]/2.d5  ; [km/s]
   VY=(VY_stop+VY_start)[*,*,istep]/2.d5  ; [km/s]
   VZ=(VZ_stop+VZ_start)[*,*,istep]/2.d5  ; [km/s]
;   
;     compute 5-point optimized derivatives
   odiffxy5,BX,BXX,BXY
   odiffxy5,BY,BYX,BYY
   odiffxy5,BZ,BZX,BZY
;
;     setup for mudpack solution of vector potential in main region
   N1=257 & N2=257
   source=fltarr(N1,N2)
   xr=float([1,n1]*dx)
   yr=float([1,n2]*dy)
;
   ixp=2
   iex=8
   jyq=2
   jey=8
;
   mudpack={nx:long(n1),ny:long(n2),ixp:long(ixp),iex:long(iex),jyq:long(jyq),jey:long(jey),iguess:long(0),maxcy:long(5),method:long(3),loud:long(0),boundary:long(1),dx:dx,dy:dy,xr:xr,yr:yr}
;
   osx=15
   osy=15
;
   if (no_elliptic2d eq 0) then begin
;
      source=Bz[osx:n1+osx-1,osy:n2+osy-1]
      compute_vector_potential,mudpack,source,phi
;      mwrfits,phi,'phi.fts',/create   ; save for use without elliptic2d
;
      source=Bzt[osx:n1+osx-1,osy:n2+osy-1]
      compute_vector_potential,mudpack,source,phit
;      mwrfits,phit,'phit.fts',/create ; save for use without elliptic2d
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
   nan=replicate(!values.d_nan,sz[1],sz[2])
   
   truth={Vx:nan,Vy:nan,Vz:nan,FTVX:nan,FTVY:nan,$
         FTVXX:nan,FTVYY:nan,$
         VPX:nan,VPY:nan,VPZ:nan,VB:nan,QX:nan,QY:nan,QZ:nan,$
         POYNTING:nan,helicity:nan,$
         Ex:nan,Ey:nan,Ez:nan,$
         PHI:nan,PHIX:nan,PHIY:nan,PHIXX:nan,PHIYY:nan,PHIXY:nan,$
         PHIRES:nan,$
         PHIT:nan,PHITX:nan,PHITY:nan,PHITXX:nan,PHITYY:nan,PHITXY:nan,$
         PHITRES:nan}
;
   truth.phi[osx:n1+osx-1,osy:n2+osy-1]=phi.scalar
   truth.phix[osx:n1+osx-1,osy:n2+osy-1]=phi.x
   truth.phiy[osx:n1+osx-1,osy:n2+osy-1]=phi.y
   truth.phixx[osx:n1+osx-1,osy:n2+osy-1]=phi.xx
   truth.phiyy[osx:n1+osx-1,osy:n2+osy-1]=phi.yy
   truth.phixy[osx:n1+osx-1,osy:n2+osy-1]=phi.xy
   truth.phires[osx:n1+osx-1,osy:n2+osy-1]=phi.res
;
   truth.phit[osx:n1+osx-1,osy:n2+osy-1]=phit.scalar
   truth.phitx[osx:n1+osx-1,osy:n2+osy-1]=phit.x
   truth.phity[osx:n1+osx-1,osy:n2+osy-1]=phit.y
   truth.phitxx[osx:n1+osx-1,osy:n2+osy-1]=phit.xx
   truth.phityy[osx:n1+osx-1,osy:n2+osy-1]=phit.yy
   truth.phitxy[osx:n1+osx-1,osy:n2+osy-1]=phit.xy
   truth.phitres[osx:n1+osx-1,osy:n2+osy-1]=phit.res
;
;     vector potential
   mag=add_tag(mag,-(truth.phiy+shift(truth.phiy,0,1))/2,"AX")     ; [G-km]
   mag=add_tag(mag,(truth.phix+shift(truth.phix,1,0))/2,"AY")      ; [G-km]
   mag=add_tag(mag,-(truth.phity+shift(truth.phity,0,1))/2,"ATX")  ; [G-km/s]
   mag=add_tag(mag,(truth.phitx+shift(truth.phitx,1,0))/2,"ATY")   ; [G-km/s]
;
   iw=where((finite(mag.Aty,/nan) ne 1) and (finite(mag.Atx,/nan) ne 1) $
                                        and (abs(mag.Bz) gt mag.thr))
;
;  print,correlate(-mag.ATY[good],truth[good].bzux)
;  print,correlate(mag.ATX[good],truth[good].bzuy)  
;
;   bzux=-mag.ATY[good]
;   bzuy=mag.ATX[good]
;
;   mag=add_tag(mag,-dchidy,"AX")
;   mag=add_tag(mag,dchidx,"AY")
;
;     checks for second order differencing
;   junkxx=(truth.phix-shift(truth.phix,1,0))/mag.dx
;   junkyy=(truth.phiy-shift(truth.phiy,0,1))/mag.dy
;   surface,junkxx-truth.phixx,charsize=4
;   surface,junkyy-truth.phiyy,charsize=4
;
   truth.VX=VX                                                    ; [km/s]
   truth.VY=VY                                                    ; [km/s]
   truth.VZ=VZ                                                    ; [km/s]
   truth.FTVX=Vx*mag.Bz-Vz*mag.BX                                 ; [G-km/s]
   truth.FTVY=Vy*mag.Bz-Vz*mag.BY                                 ; [G-km/s]
   truth.poynting=-(truth.FTVX*Mag.Bx+truth.FTVY*Mag.By)/(4*!dpi) ; [G^2-km/s]
   truth.helicity=-2*(truth.FTVX*Mag.Ax+truth.FTVY*Mag.Ay)        ;[G^2-km^2/s]
;
   good=where( (abs(mag.BZ) gt mag.thr) and (finite(truth.phix,/nan) eq 0),NG)
;
;   truth.BZUXX=xderiv(truth.bzux,mag.dx)  
;   truth.BZUYY=yderiv(truth.bzuy,mag.dy)  
;
   odiffxy5,truth.FTVX,bzuxx,bzuxy
   odiffxy5,truth.FTVY,bzuyx,bzuyy
   truth.FTVXX=bzuxx/mag.dx
   truth.FTVYY=bzuyy/mag.dy
;
;     unit vector (direction of the magnetic field)
   B=sqrt(bx^2+by^2+bz^2)
   truth.qx=bx/B
   truth.qy=by/B
   truth.qz=bz/B
;
   mag=add_tag(mag,truth.qx,'QX')
   mag=add_tag(mag,truth.qy,'QY')
   mag=add_tag(mag,truth.qz,'QZ')
;
   truth.VB=truth.qx*truth.vx+truth.qy*truth.vy+truth.qz*truth.vz
   truth.vpx=truth.vx-truth.vB*truth.qx
   truth.vpy=truth.vy-truth.vb*truth.qy
   truth.vpz=truth.vz-truth.vb*truth.qz
;
   truth.Ex=-(vy*bz-vz*by)/speed_of_light
   truth.Ey=-(vz*bx-vx*bz)/speed_of_light
   truth.Ez=-(vx*by-vy*bx)/speed_of_light
;
   return,mag

end
