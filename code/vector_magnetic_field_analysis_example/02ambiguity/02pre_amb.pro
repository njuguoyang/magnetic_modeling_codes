;+
; NAME :
;   01pre_amb
; PURPOSE :
;   Prepare the input files for the Minimum Energy precedure to remove the 180 degree ambiguity.
;-

fn = find_file('./01plot_hmi/bxyz_submap*sav',count=count)

for i=0,count-1 do begin
;for i=0,1 do begin
  restore,fn[i],/ver
  lonlat=xy2lonlat([smapbz.xc,smapbz.yc],smapbz.time)
  theta = lonlat[0]*!dtor    ;theta is the central meridian angle of the center of the field of view.
  phi   = lonlat[1]*!dtor    ;phi is the latitude of the center of the field of view.
  openw,lun,'./02removing_amb/field'+string(i,format='(I03)')+'.dat',/get_lun
  ss=size(smapbz.data)
  printf,lun,ss[1],ss[2]    ;nx, ny are the dimension of the input arrays.
  printf,lun,smapbz.dx,smapbz.dy    ;xpix, ypix are the pixel sizes in the x- and y-directions.
  printf,lun,smapbz.B0*!dtor,0.0,smapbz.RSUN    ;b is the solar b-angle. p is the solar p-angle. radius is the solar radius.
  printf,lun,theta,phi    ;theta is the central meridian angle of the center of the field of view. phi is the latitude of the center of the field of view. 
  printf,lun,smapbz.data   ;This is an nx by ny array of the line of sight component of the field (or of the magnitude of the field).
  btrans=sqrt(smapbx.data^2+smapby.data^2)
  printf,lun,btrans    ;This is an nx by ny array of the transverse component of the field (or of the inclination angle of the field). 
  bazim=atan(smapby.data,smapbx.data)
  printf,lun,bazim    ;This is an nx by ny array of the azimuthal angle, containing the ambiguity.
  close,lun
  free_lun,lun
endfor
end





