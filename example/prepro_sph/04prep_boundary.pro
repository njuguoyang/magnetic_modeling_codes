

;restore,'../01combined_boundary/Brtp_combine.sav',/ver
restore,'./Brtp_combine_midresolution.sav',/ver

filename='./field0_preprocessed.dat'
nr   = 1
nlat = 1
nlon = 1
openr, unit, filename, /get_lun
readu, unit, nr, nlat, nlon
print,nr, nlat, nlon
Br0 = dblarr(nlat,nlon)
Bt0 = dblarr(nlat,nlon)
Bp0 = dblarr(nlat,nlon)
readu, unit, Br0
readu, unit, Bt0
readu, unit, Bp0
free_lun, unit

rix = radp
theta = thep
phi = phip

Br0 = rotate(Br0,3)  ;We have to rotate the image back in the (longitude,latitude) coordinate system
Br_phy[pind1:pind2,tind1:tind2] = Br0
window,xs=1600,ys=900,/free
device,decomposed=0
loadct,4
contour,rebin(Br_phy,1200,600),rebin(phi,1200)*!radeg, 90-rebin(reverse(theta),600)*!radeg, nlevels=55,xstyle=1,ystyle=1,/fill,zrange=[-430,430],charsize=3
loadct,0
write_png,'04prep_boundary_pfss.png',tvrd(true=1)

openw,lun,'04pfss_boundary.dat',/get_lun
writeu,lun,long(6000)           ; number of points in longitude
writeu,lun,long(3000)           ; number of points in latitude
writeu,lun,double(reverse(thep))       ; coordinates in the unit of radian in latitude
writeu,lun,double(phip)                ; coordinates in the unit of radian in longitude
writeu,lun,double(Br_phy)
free_lun,lun

Br_phy = rotate(Br_phy,1)   ;Convert the image in the (colatitude,longitude) coordinate system
save,Br_phy,rix,theta,phi,tind1,tind2,pind1,pind2,filename='04prep_boundary_pfss.sav'

nr1  =240
nlat1=240
nlon1=240
lat00=390
lat0 = nt_phy-1-tind2+lat00
lat1 = nt_phy-1-tind2+lat00+nlat1-1
lon00=340
lon0 = pind1+lon00
lon1 = pind1+lon00+nlon1-1

dixB=2
mxnest=1
ngl=dixB*(2^(mxnest-1))

rmin_phy=1.0             ;Radius range [rmin_phy,rmax_phy] 
rmax_phy=2.5             ;in the unit of the solar radius
tmin_phy=0.0 *!dtor      ;Theta range [tmin_phy,tmax_phy], where theta 
tmax_phy=180.0*!dtor     ;is the colatitude. Theta=0 is the north pole, theta=180 is the south pole
pmin_phy=-180.0*!dtor    ;Phi range [pmin_phy,pmax_phy], where phi
pmax_phy= 180.0*!dtor    ;is the longitude. The origin is defined at 
                         ;the X-axis, which is defined as directing to an observer
r_res = double((rmax_phy-rmin_phy))/(nr_phy)
t_res = double((tmax_phy-tmin_phy))/nt_phy
p_res = double((pmax_phy-pmin_phy))/np_phy

print,'mxnest',mxnest
print,'nxlone1',nr1/2^(mxnest-1)
print,'nxlone2',nlat1/2^(mxnest-1)
print,'nxlone3',nlon1/2^(mxnest-1)
print,'xprobmin1',1.0+0.5*r_res
print,'xprobmax1',1.0+(0.5+double(nr1))*r_res
print,'xprobmin2',(thep[lat0]-0.5d*t_res)/(2.0d*!dpi)
print,'xprobmax2',(thep[lat1]+0.5d*t_res)/(2.0d*!dpi)
print,'xprobmin3',(phip[lon0]-0.5d*p_res+!dpi)/(2.0d*!dpi)
print,'xprobmax3',(phip[lon1]+0.5d*p_res+!dpi)/(2.0d*!dpi)

Bnl = dblarr(nlat1+2*ngl,nlon1+2*ngl,3)
Br0 = rotate(Br0,1)  ;We have to rotate the image back again in the (colatitude, longitude) coordinate system
Bnl[*,*,0] = double(Br0[lat00-ngl:lat00+nlat1-1+ngl,lon00-ngl:lon00+nlon1-1+ngl])
Bnl[*,*,1] = double(Bt0[lat00-ngl:lat00+nlat1-1+ngl,lon00-ngl:lon00+nlon1-1+ngl])
Bnl[*,*,2] = double(Bp0[lat00-ngl:lat00+nlat1-1+ngl,lon00-ngl:lon00+nlon1-1+ngl])
;  save another format in ASCII file
filename='04field0_cutout.dat'
openw, unit, filename, /get_lun
writeu, unit, long(nlat1+2*ngl), long(nlon1+2*ngl)
writeu, unit, Bnl
free_lun, unit

device,decomposed=0
loadct,4
window,xs=600,ys=600
tmp = congrid(reform(Bnl[*,*,0]),600,600)
print,min(tmp),max(tmp)
ind = where(tmp le -1500)
tmp[ind] = -1500
ind = where(tmp ge 1500)
tmp[ind] = 1500
tvscl,rotate(tmp,3)
write_png,'Brad_cutout.png',tvrd(true=1)

end
