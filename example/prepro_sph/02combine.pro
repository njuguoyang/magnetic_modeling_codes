;+
; NAME: 01combine
; PURPOSE: 1) Prepare the sub-region data, here it is the vector magnetic field by HMI
;          2) Give the sample grids of longitudes and latitudes
;          3) Get the magnetic field on the sample grids by interpolation.
;
; MODIFICATION HISTORY :
;   2012.04 Guo Yang@Nanjing University
;
;   This software is provided without any warranty. Permission to use,
;   copy, modify. Distributing modified or unmodified copies is granted,
;   provided this disclaimer and information are included unchanged.
;-

;pro 01combine

;============
;     1
;============
restore,'./01data/mapbrtp.sav',/ver
Br_loc  = mapbr.data
Bth_loc = mapbt.data
Bph_loc = mapbp.data

;============
;     2
;============
nr_phy=1500
nt_phy=3000         ;The spatial resolution is 0.06 degrees per cell, while that of HMI observations is 0.03 degrees per cell.
np_phy=6000         ;The spatial resolution is 0.06 degrees per cell

rmin_phy=1.0             ;Radius range [rmin_phy,rmax_phy] 
rmax_phy=2.5             ;in the unit of the solar radius
tmin_phy=0.0 *!dtor      ;Theta range [tmin_phy,tmax_phy], where theta 
tmax_phy=180.0*!dtor     ;is the colatitude. Theta=0 is the north pole, theta=180 is the south pole
pmin_phy=-180.0*!dtor    ;Phi range [pmin_phy,pmax_phy], where phi
pmax_phy= 180.0*!dtor    ;is the longitude. The origin is defined at 
                         ;the X-axis, which is defined as directing to an observer
radp = (rmax_phy-rmin_phy)*(findgen(nr_phy)+0.0)/(nr_phy-1)+rmin_phy   ;Note that position is defined on the cell edge
thep = (tmax_phy-tmin_phy)*(findgen(nt_phy)+0.5)/(nt_phy-0)+tmin_phy   ;Note that position is defined on the cell center
phip = (pmax_phy-pmin_phy)*(findgen(np_phy)+0.5)/(np_phy-0)+pmin_phy   ;Note that position is defined on the cell center
lat_phy   = reverse(!pi/2.0 - thep)

fn = './01data/hmi_test.B_720s_CEA.911402.20120123_030000_TAI.lat.fits'
mreadfits,fn,index,lat_loc  ;Latitude of the vector magnetic field data.
lat_loc = lat_loc*!dtor
colat_loc = reverse(!pi/2.0 - lat_loc,2)    ;Colatitude of the vector magnetic field data. 
fn = './01data/hmi_test.B_720s_CEA.911402.20120123_030000_TAI.lon.fits'
mreadfits,fn,index,lon_loc  ;Longitude of the vector magnetic field data.
carlon = tim2carr(index.date_obs)   ;The Carrington longitude of the central meridian of the sun
lon_loc = lon_loc - carlon[0]
lon_loc = lon_loc*!dtor
ss=size(lat_loc)

tmin_loc = max(lat_loc[*,0])
tmax_loc = min(lat_loc[*,ss[2]-1])
print,'Min and max theta [degree]:',tmin_loc*!radeg,tmax_loc*!radeg
pmin_loc = max(lon_loc[0,*])
pmax_loc = min(lon_loc[ss[1]-1,*])
print,'Min and max phi [degree]:',pmin_loc*!radeg,pmax_loc*!radeg

Br_phy  = fltarr(np_phy,nt_phy)
Bth_phy = fltarr(np_phy,nt_phy)
Bph_phy = fltarr(np_phy,nt_phy)
tres = !pi/nt_phy
tind1 = min(where(lat_phy ge tmin_loc-0.5*tres))
tind2 = max(where(lat_phy le tmax_loc+0.5*tres))
pres = 2.0*!pi/np_phy
pind1 = min(where(phip ge pmin_loc-0.5*tres))
pind2 = max(where(phip le pmax_loc+0.5*tres))
print,'pind1, pind2, tind1, tind2:',pind1, pind2, tind1, tind2
nt_phys = tind2 - tind1 + 1
np_phys = pind2 - pind1 + 1
print,'nt_phys, np_phys:',nt_phys,np_phys

;============
;     3
;============
dloninterp = get_interpolation_index(lon_loc[*,0],phip[pind1:pind2]) 
dlatinterp = get_interpolation_index(lat_loc[0,*],lat_phy[tind1:tind2])

Br_phy[pind1:pind2,tind1:tind2]  = interpolate(Br_loc,dloninterp,dlatinterp,/grid)
Bth_phy[pind1:pind2,tind1:tind2] = interpolate(Bth_loc,dloninterp,dlatinterp,/grid)
Bph_phy[pind1:pind2,tind1:tind2] = interpolate(Bph_loc,dloninterp,dlatinterp,/grid)

;  read the PFSS data
fn = './01data/hmi_test.Synoptic_Mr_polfil_720s.2119.synopMr_polfil.fits'
mreadfits,fn,index,Br_pfss_loc  ;Note that the latitude of these data are uniformly divided in sin(lat)
ss = size(Br_pfss_loc,/dim)
tmin_phy=-1.0      ;Sin(theta) range [tmin_phy,tmax_phy], where theta 
tmax_phy= 1.0      ;is the latitude.
pmin_phy=-180.0*!dtor    ;Phi range [pmin_phy,pmax_phy], where phi
pmax_phy= 180.0*!dtor    ;is the longitude. The origin is defined at 
                         ;the X-axis, which is defined as directing to an observer
lon_pfss_loc1 = (pmax_phy-pmin_phy)*findgen(ss[0])/(ss[0]-0)+pmin_phy   ;N points divide a circle into n parts
lat_pfss_loc1 = (tmax_phy-tmin_phy)*findgen(ss[1])/(ss[1]-1)+tmin_phy   ;N points divide a line into n-1 parts
lat_pfss_loc1 = asin(lat_pfss_loc1)

dloninterp = get_interpolation_index(lon_pfss_loc1,phip) 
dlatinterp = get_interpolation_index(lat_pfss_loc1,lat_phy)

Br_pfss_phy = interpolate(Br_pfss_loc,dloninterp,dlatinterp,/grid)

k0 = -118
nk = 10
coef = fltarr(nk)
for kk=0,nk-1 do begin
  rbr_pfss  = shift(Br_pfss_phy,k0+kk)
  subr_pfss = rbr_pfss[pind1:pind2,tind1:tind2]
  coef[kk]  = correlate(subr_pfss,Br_phy[pind1:pind2,tind1:tind2])
endfor
kk = where(coef eq max(coef))
print,'Shifting index:',k0+kk
rbr_pfss = shift(Br_pfss_phy,k0+kk) 
device,decomposed=0
loadct,4
window,xs=1600,ys=900,/free
contour,rebin(rbr_pfss,1000,500),rebin(phip,1000)*!radeg, rebin(lat_phy,500)*!radeg, nlevels=55,xstyle=1,ystyle=1,/fill,zrange=[-300,300],charsize=3
write_png,'synoptic_map_midresolution.png',tvrd(true=1)
;
print,'The normalized error of Br1 in the synotic map and Br2 of HMI snapshot in a ROI:'
print,total(abs(rbr_pfss[pind1:pind2,tind1:tind2] - Br_phy[pind1:pind2,tind1:tind2]))/total(abs(Br_phy[pind1:pind2,tind1:tind2]))
;
;Br_phy = rbr_pfss             ; only the synoptic map
rbr_pfss[pind1:pind2,tind1:tind2] = 0.0
Br_phy = rbr_pfss + Br_phy           ; combine the synoptic map and the HMI snapshot

;wdef,1,800,450
window,xs=1600,ys=900,/free
contour,rebin(Br_phy,1000,500),rebin(phip,1000)*!radeg, rebin(lat_phy,500)*!radeg, nlevels=55,xstyle=1,ystyle=1,/fill,zrange=[-300,300],charsize=3
write_png,'synoptic_frame_midresolution.png',tvrd(true=1)
loadct,0

;  save the combined data in SAV file. The location coordinates are in (phi,theta) = (longitude,latitude)
save,nr_phy,nt_phy,np_phy,Br_phy,Bth_phy,Bph_phy,radp,thep,phip,tind1,tind2,pind1,pind2,filename='Brtp_combine_midresolution.sav'  
                                                  ;Note that tind1 and tind2 is for the latitude. So the indices for
                                                  ;the colatitude should be nt_phy-1-tind2 and nt_phy-1-tind1

;  save the observed vector data in ASCII file. The location coordinates are converted to (theta,phi) = (colatitude,longitude)
;  the component of Bth_phy is also convert to colatitude coordinate system, i.e., south is positive
filename='field0_midresolution.dat'
openw, unit, filename, /get_lun
writeu, unit, nr_phy,nt_phys,np_phys
writeu, unit, rotate(Br_phy[pind1:pind2,tind1:tind2],1)
writeu, unit, (-1.0)*rotate(Bth_phy[pind1:pind2,tind1:tind2],1)
writeu, unit, rotate(Bph_phy[pind1:pind2,tind1:tind2],1)
writeu, unit, radp
writeu, unit, thep[nt_phy-1-tind2:nt_phy-1-tind1]
writeu, unit, phip[pind1:pind2]
free_lun, unit

openw,lun,'pfss_boundary_midres.dat',/get_lun
writeu,lun,long(np_phy)          ; number of points in longitude
writeu,lun,long(nt_phy)           ; number of points in latitude
writeu,lun,double(reverse(thep))       ; coordinates in the unit of radian in latitude
writeu,lun,double(phip+!dpi)           ; coordinates in the unit of radian in longitude. Note that the coordinate of the left edge is converted to be 0.
writeu,lun,double(Br_phy)              ; in the (longitude,latitude) coordinate system
free_lun,lun


end
