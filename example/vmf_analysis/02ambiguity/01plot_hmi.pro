;pro plot_hmi

fn1 = '../01data/hmi.B_720s.20150827_052400_TAI.field.fits'
read_sdo,fn1,index1,data1    ;add /noshell keyword in Windows
index2map,index1,data1,mapb
fn2 = '../01data/hmi.B_720s.20150827_052400_TAI.inclination.fits'
read_sdo,fn2,index2,data2    ;add /noshell keyword in Windows
index2map,index2,data2,mapi
fn3='../01data/hmi.B_720s.20150827_052400_TAI.azimuth.fits'
read_sdo,fn3,index3,data3    ;add /noshell keyword in Windows
index2map,index3,data3,mapa

mapbz = mapb
mapbz.data = rotate(mapb.data*cos(mapi.data*!dtor),2)
mapbz.roll_angle = 0.0
dx0 = 0.5*(965.0 - 935.0)
dy0 = 0.5*(944.0 - 954.0)
mapbz.xc = mapbz.xc - dx0         ;Alignment by comparing the position of the solar limbs
mapbz.yc = mapbz.yc - dy0
mapbx = mapb
mapbx.data = rotate(mapb.data*sin(mapi.data*!dtor)*cos((mapa.data + 270.0)*!dtor),2)  ;The azimuth angle is measured from the CCD+y direction, which is the south, since the solar P angle is ~180 degree
mapbx.roll_angle = 0.0
index = where(mapbx.data le -10000.0)
mapbx.data[index] = 0.0
mapbx.xc = mapbz.xc
mapbx.yc = mapbz.yc
mapby = mapb
mapby.data = rotate(mapb.data*sin(mapi.data*!dtor)*sin((mapa.data + 270.0)*!dtor),2)  ;The azimuth angle is measured from the CCD+y direction, which is the south, since the solar P angle is ~180 degree
mapby.roll_angle = 0.0
mapby.data[index] = 0.0
mapby.xc = mapbz.xc
mapby.yc = mapbz.yc

intimei = utc2int(mapbz.time)
mjdi  = intimei.mjd + fix((intimei.time + 120000.0)/(24.0*3600.0*1000.0))
timei = (intimei.time + 120000.0) mod (24.0*3600.0*1000.0)
intimei.mjd = mjdi
intimei.time= timei
utctime = anytim2utc(intimei,/vms)
mapbz.time = utctime
mapbx.time = utctime
mapby.time = utctime

loadct,0
wdef,1,800,800
!p.background = 255
plot_map,mapbz,dmax=2000,dmin=-2000,color=0,charsize=1.8
write_png,'./01plot_hmi/bz.png',tvrd()
wdef,2,800,800
plot_map,mapbz,dmax=2000,dmin=-2000,color=0,charsize=1.8
plot_vmap,/over,mapbx,mapby,mapbz=mapbz,limit=50,scale=0.1,iskip=30,jskip=30,$
          v_color=255,axis_color=0,/Nolabels,v_thick=2.0,/No_arrow_head  ;,/Noaxis
write_png,'./01plot_hmi/bxyz.png',tvrd()

sub_map,mapbz,smapbz,xrange=[300,800.0],yrange=[-500,-100.0]
print,'Flux balance coefficient:',total(smapbz.data)/total(abs(smapbz.data))
sub_map,mapbx,smapbx,ref_map=smapbz                                 
sub_map,mapby,smapby,ref_map=smapbz
wdef,1,800,800
plot_map,smapbz,dmax=2000,dmin=-2000,color=0,charsize=1.8
write_png,'./01plot_hmi/sbz.png',tvrd()
wdef,2,800,800
plot_map,smapbz,dmax=2000,dmin=-2000,color=0,charsize=1.8
plot_vmap,/over,smapbx,smapby,mapbz=smapbz,limit=180,scale=0.012,iskip=15,jskip=15,$
          v_color=255,axis_color=0,/Nolabels,v_thick=2.0 ;,/No_arrow_head  ;,/Noaxis
write_png,'./01plot_hmi/sbxyz.png',tvrd()
i=0
save,smapbx,smapby,smapbz,filename='./01plot_hmi/bxyz_submap.sav'

end
