
fn1 = '../01data/hmi.B_720s.20150827_052400_TAI.field.fits'
read_sdo,fn1,index1,data1
index2map,index1,data1,mapb
fn2 = '../01data/hmi.B_720s.20150827_052400_TAI.inclination.fits'
read_sdo,fn2,index2,data2
index2map,index2,data2,mapi

mapbz = mapb
mapbz.data = rotate(mapb.data*cos(mapi.data*!dtor),2)
mapbz.roll_angle = 0.0
dx0 = 0.5*(965.0 - 935.0)
dy0 = 0.5*(944.0 - 954.0)
mapbz.xc = mapbz.xc - dx0         ;Alignment by comparing the position of the solar limbs
mapbz.yc = mapbz.yc - dy0

ss = size(mapbz.data)
nx = ss[1]
ny = ss[2]
r = 69.55
dx = mapbz.dx*r/mapbz.rsun
dy = mapbz.dy*r/mapbz.rsun
x00 = (0.5 - nx/2.0)*mapbz.dx + mapbz.xc
x00 = x00*r/mapbz.rsun
y00 = (0.5 - ny/2.0)*mapbz.dy + mapbz.yc
y00 = y00*r/mapbz.rsun

openw,1,'20150827_hmi_bz_structured.vtk'
printf,1,'# vtk DataFile Version 2.0'
printf,1,'hmi bz'
printf,1,'ASCII'
printf,1,'DATASET STRUCTURED_POINTS'
printf,1,'DIMENSIONS 1 4096 4096'
printf,1,'ORIGIN 0 '+strtrim(string(x00),2)+' '+strtrim(string(y00),2)
printf,1,'SPACING 0.0369620 0.0369620 0.0369620'

printf,1,'POINT_DATA '+strtrim(string(long(nx)*ny),2)
printf,1,'SCALARS hmi_bz int 1'
printf,1,'LOOKUP_TABLE default'
for j=0,ny-1 do begin
  for i=0,nx-1 do begin   
    printf,1,mapbz.data[i,j],format='(i8)'
  endfor
endfor

close,1

end
