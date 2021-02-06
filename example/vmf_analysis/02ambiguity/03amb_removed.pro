;+
; NAME :
;   03amb_removed
; PURPOSE :
;   Prepare the 180 degree ambiguity removed maps.
;-

;pro 03amb_removed

file=find_file('./01plot_hmi/bxyz*sav')
filea=find_file('./02removing_amb/azimuth*dat')

for i=0,n_elements(file)-1 do begin
  restore,file[i]
  btrans=sqrt(smapbx.data^2+smapby.data^2)
  ss=size(smapbz.data)
  bazim=fltarr(ss[1],ss[2])
  openr,1,filea[i]
  readf,1,bazim
  close,1
  smapbx.data=btrans*cos(bazim)
  smapby.data=btrans*sin(bazim)
  sub_map,smapbz,smapbz1,xrange=[500.0,800.0],yrange=[-500,-100.0]
  sub_map,smapbx,smapbx1,xrange=[500.0,800.0],yrange=[-500,-100.0]
  sub_map,smapby,smapby1,xrange=[500.0,800.0],yrange=[-500,-100.0]
  smapbz=smapbz1
  smapbx=smapbx1
  smapby=smapby1
  save,smapbz,smapbx,smapby,filename='./03amb_removed/smapbxyz'+string(i,format='(I03)')+'.sav'
  wdef,1,800,800
  !P.background=255
  plot_map,smapbz,charsize=2.0,color=0,dmax=2000,dmin=-2000
  write_png,'./03amb_removed/smapbz'+string(i,format='(I03)')+'.png',tvrd()
  wdef,2,800,800
  plot_map,smapbz,color=0,bcolor=0,charsize=2.0,dmax=2000,dmin=-2000
  plot_vmap,/over,smapbx,smapby,mapbz=smapbz,limit=180,scale=0.004,iskip=11,jskip=11,$
          v_color=255,axis_color=0,/Nolabels,/Noaxis,v_thick=2.0
  !P.background=0
  write_png,'./03amb_removed/smapbxyz'+string(i,format='(I03)')+'.png',tvrd()
endfor

end
