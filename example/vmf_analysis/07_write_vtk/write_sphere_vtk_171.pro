
fn = './data/aia.lev1_euv_12s.2010-11-11T060013Z.171.image_lev1.fits'
read_sdo,fn,out_index,out_data 
aia_lct,wave=out_index.wavelnth,/load
index2map,out_index,out_data,map
map.data = map.data/map.dur
map.dur = 1.0
map.xc = map.xc + 10.0

nlon = 4001l
nlat = 4001l

openw,1,'20101111_171_sphere.vtk'
printf,1,'# vtk DataFile Version 2.0'
printf,1,'solar sphere'
printf,1,'ASCII'
printf,1,'DATASET POLYDATA'
printf,1,'POINTS '+strtrim(string(nlon*nlat),2)+' float'
r = 69.55
for i=0,nlon-1 do begin
  for j=0,nlat-1 do begin   
    theta = j/float(nlat-1)*!dpi
    phi=i*!dpi/float(nlon-1)-!dpi/2.0 
    x=R*sin(theta)*cos(phi)
    y=R*sin(theta)*sin(phi)
    z=R*cos(theta)
    printf,1,x,y,z,format='(f10.5,1x,f10.5,1x,f10.5)'
  endfor
endfor

printf,1,'POLYGONS '+strtrim(string((nlon-1)*(nlat-1)),2)+' '+strtrim(string((nlon-1)*(nlat-1)*5l),2)
for i=0,nlon-2 do begin
  for j=0,nlat-2 do begin               
    printf,1,4,j+i*nlat,j+(i+1l)*nlat,j+(i+1l)*nlat+1l,j+i*nlat+1l,format='(i5,1x,i16,1x,i16,1x,i16,1x,i16)'   
  endfor
endfor

ss = size(map.data)
nx = ss[1]
ny = ss[2]
xx_loc = fltarr(nx)
yy_loc = fltarr(ny)
for i=0,nx-1 do begin
  xx_loc[i] = (i + 0.5 - nx/2.0)*map.dx + map.xc
endfor
for j=0,ny-1 do begin
  yy_loc[j] = (j + 0.5 - ny/2.0)*map.dy + map.yc
endfor
data = fltarr(nlon-1l,nlat-1l)
R0 = map.rsun
yy_phy = fltarr(nlon-1l,nlat-1l)
zz_phy = fltarr(nlon-1l,nlat-1l)
for i=0,nlon-2l do begin
  for j=0,nlat-2l do begin   
    theta = j/float(nlat-1)*!dpi
    phi=i*!dpi/float(nlon-1)-!dpi/2.0 
    ;xx_phy=R0*sin(theta)*cos(phi)  
    yy_phy[i,j]=R0*sin(theta)*sin(phi)
    zz_phy[i,j]=R0*cos(theta)
  endfor  
endfor
print,'OK!!!!!!!!!'
indx = get_interpolation_index(xx_loc,yy_phy)
indy = get_interpolation_index(yy_loc,zz_phy) 
data = interpolate(map.data,indx,indy)
data = reform(data,nlon-1,nlat-1)
data = alog10(float(data))

printf,1,'CELL_DATA '+strtrim(string((nlon-1)*(nlat-1)),2)
printf,1,'SCALARS cell_scalars float 1'
printf,1,'LOOKUP_TABLE default'
for i=0,nlon-2l do begin
  for j=0,nlat-2l do begin   
    printf,1,data[i,j],format='(f10.6)'
  endfor  
endfor

;tvlct,red,green,blue,/get
;red = float(red)/255.0
;green = float(green)/255.0
;blue = float(blue)/255.0
;min_data = min(data)
;max_data = max(data)
;printf,1,'LOOKUP_TABLE aia193 '+strtrim(string((nlon-1)*(nlat-1)),2)
;for i=0,nlon-2l do begin
;  for j=0,nlat-2l do begin  
;    ind = fix(255.0*(data[i,j]-min_data)/(max_data-min_data)) 
;    printf,1,red[ind],green[ind],blue[ind],1.0,format='(f9.7,1x,f9.7,1x,f9.7,1x,f3.1)'
;  endfor  
;endfor

close,1

end