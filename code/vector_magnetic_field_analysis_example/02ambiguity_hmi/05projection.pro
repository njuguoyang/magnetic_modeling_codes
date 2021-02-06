;projection

restore,'./01plot_hmi/bxyz_submap.sav',/ver
mapbz=smapbz
mapbx=smapbx
mapby=smapby

ssmtr=size(mapbz.data)
dxmtr=mapbz.dx
dymtr=mapbz.dy
xcmtr=mapbz.xc
ycmtr=mapbz.yc
time_mtr=mapbz.time

pix_size=[dxmtr,dymtr]
ans = get_rb0p(time_mtr, _extra=extra)
b0in = reform(ans[1,*])
sunrin = reform(ans[0,*])
result=xy2lonlat([xcmtr,ycmtr],time_mtr)
cmd=result[0]*!dtor
lat=result[1]*!dtor
ang=pb0r(time_mtr)
B0=ang[1]*!dtor
P=0.0
xindex=findgen(ssmtr[1])-(ssmtr[1]-1)/2.0
yindex=findgen(ssmtr[2])-(ssmtr[2]-1)/2.0
xcoor=xcmtr+dxmtr*xindex
ycoor=ycmtr+dymtr*yindex
cmd1=dblarr(ssmtr[1],ssmtr[2])
latitude=dblarr(ssmtr[1],ssmtr[2])

files=find_file('./05projection/cmd_lat.sav',count=count)
if (count eq 0) then begin
  for j=0,ssmtr[2]-1 do begin
    if (j mod 50 eq 0) then print,'Total step is:',ssmtr[2],'  Current step is:',j+1
    for i=0,ssmtr[1]-1 do begin
      result=xy2lonlat([xcoor[i],ycoor[j]],time_mtr, b0=b0in, radius=sunrin)
      cmd1[i,j]=result[0]
      latitude[i,j]=result[1]
    endfor
  endfor
  save,cmd1,latitude,filename='./05projection/cmd_lat.sav'
endif else begin
  restore,'./05projection/cmd_lat.sav',/ver
endelse

Blong=mapbz.data
Btrans=(mapbx.data^2+mapby.data^2)^0.5
azim=fltarr(ssmtr[1],ssmtr[2])
azim=-atan(mapbx.data,mapby.data)*!radeg
;===NPFC codes accept the azimuth angles both in degree or radian, but require the angles measuring from the north and counterclockwise===;
raw=create_struct('B_long',Blong,'B_trans',Btrans,'B_azim',azim,'latitude',latitude,'cmd',cmd1,$
    'point',create_struct('pix_size',pix_size,'cmd',cmd,'lat',lat,'B0',B0,'P',P))
avesz=(raw.point.pix_size[0]+raw.point.pix_size[1])/2.0
raw.point.pix_size=[avesz,avesz]

sig_Bl=5.0
sig_Btr=10.0
hms=3.0
IVM_GET_NEW,raw,Bl,Btr,phi_unres,iapp_ref,iapp_tot,pixel_size,Bc,Lc,Bcc,Lcc,B0,L0,P,sig_Bl,sig_Btr,hms
;==============================================
;from LOS components to heliographic components
;==============================================
Tx1=-sin(B0)*sin(P)*sin(Lc-L0) + cos(P)*cos(Lc-L0)
Tx2=sin(B0)*cos(P)*sin(Lc-L0) + sin(P)*cos(Lc-L0)
Tx3=-cos(B0)*sin(Lc-L0)
Ty1=-sin(Bc)*(sin(B0)*sin(P)*cos(Lc-L0) + cos(P)*sin(Lc-L0)) $
             -cos(Bc)*cos(B0)*sin(P)
Ty2=sin(Bc)*(sin(B0)*cos(P)*cos(Lc-L0) - sin(P)*sin(Lc-L0)) $
                  +cos(Bc)*cos(B0)*cos(P)
Ty3=-cos(B0)*sin(Bc)*cos(Lc-L0) + sin(B0)*cos(Bc)
Tz1=cos(Bc)*(sin(B0)*sin(P)*cos(Lc-L0) + cos(P)*sin(Lc-L0)) $ 
            -sin(Bc)*cos(B0)*sin(P)
Tz2=-cos(Bc)*(sin(B0)*cos(P)*cos(Lc-L0) - sin(P)*sin(Lc-L0)) $ 
                   +sin(Bc)*cos(B0)*cos(P)
Tz3=cos(Bc)*cos(B0)*cos(Lc-L0) + sin(Bc)*sin(B0)
;
Bt=sqrt(Btrans^2.+Blong^2.)
Bx=Btr*cos((azim + 90.0)*!dtor)
By=Btr*sin((azim + 90.0)*!dtor)
Bxn=Bx/Bt
Byn=By/Bt
Bzn=Blong/Bt
r1=where(Bt eq 0.,ico) & if ico gt 0. then begin & Bxn(r1)=0. & Byn(r1)=0. &$
Bzn(r1)=0. &endif
Bx_loc=Bt*(Bxn*Tx1+ Byn*Tx2+ Bzn*Tx3)
By_loc=Bt*(Bxn*Ty1+ Byn*Ty2+ Bzn*Ty3)
Bz_loc=Bt*(Bxn*Tz1+ Byn*Tz2+ Bzn*Tz3)
Bh_loc=sqrt(Bx_loc^2+By_loc^2)
Phi_loc=-atan(Bx_loc,By_loc)   ;in the unit of radian, the angles are measured from the north (CCD+y) and counterclockwise
;======================================
;from image plane to heliographic plane
;======================================
res=size(Bl) & idim1=res(1) & idim2=res(2)
limb=0.
;  Original version by Manolis K. Georgoulis (JHU/APL, 10/13/05)
;c11=cos(P)*cos(Lc-L0)- sin(P)*sin(B0)*sin(Lc-L0)
;c12=-cos(P)*sin(Bc)*sin(Lc-L0) - sin(P)*(cos(B0)*cos(Bc) + sin(B0)*sin(Bc)*$
;                                         cos(Lc-L0))
;c21=sin(P)*cos(Lc-L0) + cos(P)*sin(B0)*sin(Lc-L0)
;c22=-sin(P)*sin(Bc)*sin(Lc-L0) + cos(P)*(cos(B0)*cos(Bc) + sin(B0)*sin(Bc)*$
;                                         cos(Lc-L0))   
     
;  Modified by Yang Guo (NJU, 2016 Nov. 22)
xyc_ll = [0.5*(max(lc) + min(lc))*!radeg,0.5*(max(bc) + min(bc))*!radeg]    ;The longitude and latitude of the field of view
Lcc = xyc_ll[0]*!dtor  &  Bcc = xyc_ll[1]*!dtor                             ;The longitude and latitude of the point at which the heliographic plane is tangent
print,'Lcc, Bcc (degree):',Lcc*!radeg,Bcc*!radeg
print,'Longitude (radian):',Lcc
print,'Latitude (radian):',Bcc
print,'B0 (radian):',B0
c11=cos(P)*cos(Lcc-L0)- sin(P)*sin(B0)*sin(Lcc-L0)
c12=-cos(P)*sin(Bcc)*sin(Lcc-L0) - sin(P)*(cos(B0)*cos(Bcc) + sin(B0)*sin(Bcc)*$
                                         cos(Lcc-L0))
c21=sin(P)*cos(Lcc-L0) + cos(P)*sin(B0)*sin(Lcc-L0)
c22=-sin(P)*sin(Bcc)*sin(Lcc-L0) + cos(P)*(cos(B0)*cos(Bcc) + sin(B0)*sin(Bcc)*$
                                         cos(Lcc-L0))

denom=1./(c11*c22 - c12*c21)
r1=where(denom gt 3.5 or denom lt 0.,ico) & if ico gt 0. then begin & $ 
denom(r1)=0. & limb=1. & endif
alpha=c22*denom
beta=-c21*denom
gamm=-c12*denom
delta=c11*denom
xi=dblarr(idim1,idim2)
yi=dblarr(idim1,idim2)
if n_elements(alpha) gt 1. then begin 
;  for i=0,idim1-1 do for j=0,idim2-1 do xi(i,j)=alpha(i,j)*i+gamm(i,j)*j   ;Original version by Manolis K. Georgoulis (JHU/APL, 10/13/05)
;  for i=0,idim1-1 do for j=0,idim2-1 do yi(i,j)=beta(i,j)*i+delta(i,j)*j   ;Original version by Manolis K. Georgoulis (JHU/APL, 10/13/05)
  for i=0,idim1-1 do for j=0,idim2-1 do xi(i,j)=alpha(i,j)*(i-0.5*(idim1-1)) + gamm(i,j)*(j-0.5*(idim2-1))    ;Modified by Yang Guo (NJU, 2011 Sep. 30)
  for i=0,idim1-1 do for j=0,idim2-1 do yi(i,j)=beta(i,j)*(i-0.5*(idim1-1)) + delta(i,j)*(j-0.5*(idim2-1))    ;Modified by Yang Guo (NJU, 2011 Sep. 30)
endif else begin 
;  for i=0,idim1-1 do for j=0,idim2-1 do xi(i,j)=alpha*i+gamm*j   ;Original version by Manolis K. Georgoulis (JHU/APL, 10/13/05)
;  for i=0,idim1-1 do for j=0,idim2-1 do yi(i,j)=beta*i+delta*j   ;Original version by Manolis K. Georgoulis (JHU/APL, 10/13/05)
  ;print,'Test OK!'
  for i=0,idim1-1 do for j=0,idim2-1 do xi(i,j)=alpha*(i-0.5*(idim1-1)) + gamm*(j-0.5*(idim2-1))    ;Modified by Yang Guo (NJU, 2016 Nov. 22)
  for i=0,idim1-1 do for j=0,idim2-1 do yi(i,j)=beta*(i-0.5*(idim1-1)) + delta*(j-0.5*(idim2-1))    ;Modified by Yang Guo (NJU, 2016 Nov. 22)
endelse
ksi=c11*xi + c12*yi
eta=c21*xi + c22*yi
ksi = ksi - min(ksi)   ;Modified by Yang Guo (NJU, 2012 Apr. 23) to be consistent with the above modifications
eta = eta - min(eta)   ;Modified by Yang Guo (NJU, 2012 Apr. 23) to be consistent with the above modifications

;id1=fix(max(xi)-min(xi)+0.5)   ;original version by Manolis K. Georgoulis (JHU/APL, 10/13/05)
;id2=fix(max(yi)-min(yi)+0.5)
id1=fix(max(xi)-min(xi)+1.0)   ;modified by Yang Guo (NJU, 2011 Jun. 5)
id2=fix(max(yi)-min(yi)+1.0)
xi=xi-min(xi)
yi=yi-min(yi)

bz_hel=dblarr(id1,id2)
bh_hel=dblarr(id1,id2)
phi_hel=dblarr(id1,id2)
bz_hel(fix(xi+0.5),fix(yi+0.5))=bz_loc(fix(ksi+0.5),fix(eta+0.5))
  mask=fltarr(id1,id2)
  q1=where(bz_hel eq 0.,count) & if count gt 0. then mask(q1)=1.
  tmp=mask
  ms=shift(mask,-1,0)+shift(mask,1,0)+shift(mask,0,-1)+shift(mask,0,1)+ $
     shift(mask,-1,-1)+shift(mask,-1,1)+shift(mask,1,-1)+shift(mask,1,1)
  q1=where(mask eq 1. and ms eq 8,count) & if count gt 0. then tmp(q1)=0.
  mask=tmp
r1=where(mask eq 1.,ico)
if ico gt 0. then for m=0,49 do $
;bz_hel(r1)=(1./4.)*(bz_hel(r1-1)+bz_hel(r1+1)+bz_hel(r1+id1)+bz_hel(r1-id1))
bz_hel(r1)=(1./8.)*(bz_hel(r1-1)+bz_hel(r1+1)+bz_hel(r1+id1)+bz_hel(r1+id1+1)+bz_hel(r1+id1-1)+bz_hel(r1-id1)+bz_hel(r1-id1+1)+bz_hel(r1-id1-1))

bh_hel(fix(xi+0.5),fix(yi+0.5))=bh_loc(fix(ksi+0.5),fix(eta+0.5))
  mask=fltarr(id1,id2)
  q1=where(bh_hel eq 0.,count) & if count gt 0. then mask(q1)=1.
  tmp=mask
  ms=shift(mask,-1,0)+shift(mask,1,0)+shift(mask,0,-1)+shift(mask,0,1)+ $
     shift(mask,-1,-1)+shift(mask,-1,1)+shift(mask,1,-1)+shift(mask,1,1)
  q1=where(mask eq 1. and ms eq 8,count) & if count gt 0. then tmp(q1)=0.
  mask=tmp
r1=where(mask eq 1.,ico)
if ico gt 0. then for m=0,49 do $
;bh_hel(r1)=(1./4.)*(bh_hel(r1-1)+bh_hel(r1+1)+bh_hel(r1+id1)+bh_hel(r1-id1))
bh_hel(r1)=(1./8.)*(bh_hel(r1-1)+bh_hel(r1+1)+bh_hel(r1+id1)+bh_hel(r1+id1+1)+bh_hel(r1+id1-1)+bh_hel(r1-id1)+bh_hel(r1-id1+1)+bh_hel(r1-id1-1))

phi_hel(fix(xi+0.5),fix(yi+0.5))=phi_loc(fix(ksi+0.5),fix(eta+0.5))
  mask=fltarr(id1,id2)
  q1=where(phi_hel eq 0.,count) & if count gt 0. then mask(q1)=1.
  tmp=mask
  ms=shift(mask,-1,0)+shift(mask,1,0)+shift(mask,0,-1)+shift(mask,0,1)+ $
     shift(mask,-1,-1)+shift(mask,-1,1)+shift(mask,1,-1)+shift(mask,1,1)
  q1=where(mask eq 1. and ms eq 8,count) & if count gt 0. then tmp(q1)=0.
  mask=tmp
r1=where(mask eq 1.,ico)
if ico gt 0. then for m=0,49 do $
phi_hel(r1)=(1./8.)*(phi_hel(r1-1)+phi_hel(r1+1)+phi_hel(r1+id1)+phi_hel(r1+id1+1)+phi_hel(r1+id1-1)+phi_hel(r1-id1)+phi_hel(r1-id1+1)+phi_hel(r1-id1-1))

ss=size(bz_hel)
date0=time_mtr
pbr=pb0r(date0)
xyc = lonlat2xy(xyc_ll,date0)                       ;The x- and y-coordinates in arcsec of the center of the field of view
;dx=2.0*tan(0.5*(max(lc)-min(lc)))*cos(min(bc))*pbr[2]*60.0/float(ss[1])
;dy=2.0*tan(0.5*(max(bc)-min(bc)))*pbr[2]*60.0/float(ss[2])
dx = (tan(max(lc) - Lcc) + tan(Lcc - min(lc)))*cos(Bcc)*pbr[2]*60.0/float(ss[1])   ;Modified by Yang Guo (NJU, 2013 Apr. 2)
dy = (tan(max(bc) - Bcc) + tan(Bcc - min(bc)))*pbr[2]*60.0/float(ss[2])            ;Modified by Yang Guo (NJU, 2013 Apr. 2)
data=bz_hel
mapbz=create_struct('data',data,'xc',xyc[0],'yc',xyc[1],'dx',dx,'dy',dy $
                   ,'time',date0,'dur',0.0,'units','arcsec','ID','Bz' $
                   ,'roll_angle',0.0,'roll_center',[0.0,0.0])

data=bh_hel*cos(phi_hel + 0.5*!dpi)   ;The azimuth angles are measured counterclockwise from the solar north
mapbx=create_struct('data',data,'xc',xyc[0],'yc',xyc[1],'dx',dx,'dy',dy $
                   ,'time',date0,'dur',0.0,'units','arcsec','ID','Bx' $
                   ,'roll_angle',0.0,'roll_center',[0.0,0.0])

data=bh_hel*sin(phi_hel + 0.5*!dpi)   ;The azimuth angles are measured counterclockwise from the solar north
mapby=create_struct('data',data,'xc',xyc[0],'yc',xyc[1],'dx',dx,'dy',dy $
                   ,'time',date0,'dur',0.0,'units','arcsec','ID','By' $
                   ,'roll_angle',0.0,'roll_center',[0.0,0.0])

save,mapbx,mapby,mapbz,filename='./05projection/mapbxyz_m.sav'
wdef,1,1300,800
!p.background=255
mapbz.id='Heliographic Bxyz'
plot_map,mapbz,color=0,bcolor=0,charsize=2.0,dmax=800,dmin=-800
write_png,'./05projection/bz_m.png',tvrd()
plot_vmap,/over,mapbx,mapby,mapbz=mapbz,limit=180,scale=0.008,iskip=11,jskip=11,$
          v_color=255,axis_color=0,/Nolabels,v_thick=2.0,/Noaxis
!p.background=0
write_png,'./05projection/bxyz_m.png',tvrd()


end
