
fn = './data/aia.lev1_euv_12s.2010-11-11T060013Z.171.image_lev1.fits'
read_sdo,fn,out_index,out_data 
aia_lct,wave=out_index.wavelnth,/load
index2map,out_index,out_data,map
map.data = map.data/map.dur
map.dur = 1.0
map.xc = map.xc + 10.0
ss = size(map.data)
nx = ss[1]
ny = ss[2]
r = 69.55
dx = map.dx*r/map.rsun
dy = map.dy*r/map.rsun
x00 = (0.5 - nx/2.0)*map.dx + map.xc
x00 = x00*r/map.rsun
y00 = (0.5 - ny/2.0)*map.dy + map.yc
y00 = y00*r/map.rsun

openw,1,'20101111_171_structured.vtk'
printf,1,'# vtk DataFile Version 2.0'
printf,1,'aia 171'
printf,1,'ASCII'
printf,1,'DATASET STRUCTURED_POINTS'
printf,1,'DIMENSIONS 1 4096 4096'
printf,1,'ORIGIN 0 '+strtrim(string(x00),2)+' '+strtrim(string(y00),2)
printf,1,'SPACING 0.043049437 0.043049437 0.043049437'

printf,1,'POINT_DATA '+strtrim(string(long(nx)*ny),2)
printf,1,'SCALARS aia_171 int 1'
printf,1,'LOOKUP_TABLE default'
for j=0,ny-1 do begin
  for i=0,nx-1 do begin   
    printf,1,map.data[i,j],format='(i8)'
  endfor
endfor

close,1

end