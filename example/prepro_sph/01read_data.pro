;pro read_data

fn = './01data/hmi_test.Synoptic_Mr_polfil_720s.2119.synopMr_polfil.fits'
read_sdo,fn,index,data  ;Note that the latitude of these data are uniformly divided in sin(lat)
help,index,data
wdef,1,3600,1440
device,decomposed=0
loadct,6
tvscl,data
write_png,'./01data/hmi_test.Synoptic_Mr_polfil_720s.2119.synopMr_polfil.png',tvrd(true=1)

fnlon='./01data/hmi_test.B_720s_CEA.911402.20120123_030000_TAI.lon.fits'
read_sdo,fnlon,indlon,lon
fnlat='./01data/hmi_test.B_720s_CEA.911402.20120123_030000_TAI.lat.fits'
read_sdo,fnlat,indlat,lat

fnBr = './01data/hmi_test.B_720s_CEA.911402.20120123_030000_TAI.Br.fits'
read_sdo,fnBr,indBr,Br
index2map,indBr,Br,mapbr
lon1 = mapbr.xc - tim2carr(mapbr.time)  ;The longitude difference between the map center to the central meridian of the Sun
mapbr.xc = lon1
mapbr.dx = indBr.cdelt1
mapbr.dy = indBr.cdelt2
dims=size(br,/dimensions)
print,'dims of Br:',dims
loadct,0
!p.background=255
window,xs=2000,ys=1600,/free
plot_map,mapbr,color=0,bcolor=0,charsize=4.0,charthick=2,dmin=-1000,dmax=1000
write_png,'./01data/mapbr.png',tvrd()

fnBp = './01data/hmi_test.B_720s_CEA.911402.20120123_030000_TAI.Bp.fits'
read_sdo,fnBp,indBp,Bp
index2map,indBp,Bp,mapbp
mapbp.xc = mapbr.xc
mapbp.yc = mapbr.yc
mapbp.dx = mapbr.dx
mapbp.dy = mapbr.dy
fnBt = './01data/hmi_test.B_720s_CEA.911402.20120123_030000_TAI.Bt.fits'
read_sdo,fnBt,indBt,Bt
index2map,indBt,Bt,mapbt
mapbt.data=-mapbt.data
mapbt.xc = mapbr.xc
mapbt.yc = mapbr.yc
mapbt.dx = mapbr.dx
mapbt.dy = mapbr.dy

window,xs=2000,ys=1600,/free
loadct,0
!p.background=255
plot_map,mapbr,color=0,bcolor=0,charsize=4.0,charthick=2,dmin=-1000,dmax=1000
loadct,13
plot_vmap,/over,mapbp,mapbt,mapbz=mapbr,limit=200,scale=0.002,iskip=15,jskip=15,v_thick=2.0,/sample,index_size=1000,$
          index_color=233,index_charsize=3.0,v_color=255,axis_color=0,title_charsize=0.7,/Nolabels,/Noaxis
write_png,'./01data/mapbrpt.png',tvrd(true=1)

window,xs=2000,ys=1600,/free
loadct,0
!p.background=255
plot_map,mapbr,color=0,bcolor=0,charsize=4.0,charthick=2,dmin=-1000,dmax=1000,xrange=[0,55],yrange=[5,45]
loadct,13
plot_vmap,/over,mapbp,mapbt,mapbz=mapbr,limit=200,scale=0.002,iskip=15,jskip=15,v_thick=2.0,/sample,index_size=1000,$
          index_color=233,index_charsize=3.0,v_color=255,axis_color=0,title_charsize=0.7,/Nolabels,/Noaxis
write_png,'./01data/mapbrpt_sub.png',tvrd(true=1)

save,mapbr,mapbt,mapbp,filename='./01data/mapbrtp.sav'
!p.background=0

end
