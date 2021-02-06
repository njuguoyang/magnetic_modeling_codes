pro make_vmap02

DIR='./GUO2/'
DIR1='./PIC2/'
file=find_file(DIR+'mapbxyz*')
!p.background=255
wdef,1,800,800
for i=0,n_elements(file)-2 do begin
  file1=file[i]
  file2=file[i+1]
  restore,file2
  internal_time=anytim2utc(smapbz.time)  
  t_stop=double(internal_time.time)
  restore,file1
  internal_time=anytim2utc(smapbz.time)
  t_start=double(internal_time.time)
  t_middle=0.5*(t_stop+t_start)
  internal_time.time=t_middle
  utc=anytim2utc(internal_time,/vms)
  dmap=drot_map(smapbz,time=utc)
  restore,DIR+'dave4vm'+strtrim(string(i,format='(I03)'),2)+'.sav'
  sz=size(vel.u0)
  x1=vel.window_size[0]
  x2=sz[1]-1-vel.window_size[0]
  y1=vel.window_size[1]
  y2=sz[2]-1-vel.window_size[1]
  sub_map,dmap,smap,xrange=[x1,x2],yrange=[y1,y2],/pixel
  mapvx=smap
  mapvx.data=vel.u0[x1:x2,y1:y2]
  ;mapvx.data=vel.u0[x1:x2,y1:y2]-vel.w0[x1:x2,y1:y2]/mag.bz[x1:x2,y1:y2]*mag.bx[x1:x2,y1:y2]
  mapvy=smap
  mapvy.data=vel.v0[x1:x2,y1:y2]
  ;mapvy.data=vel.v0[x1:x2,y1:y2]-vel.w0[x1:x2,y1:y2]/mag.bz[x1:x2,y1:y2]*mag.by[x1:x2,y1:y2]
  mapvz=smap
  mapvz.data=vel.w0[x1:x2,y1:y2]
  mapvz.id='DAVE4VM Velocity Map'
  smap.id='DAVE4VM Velocity Map'
  plot_map,smap,color=0,bcolor=0,charsize=1.8
  plot_vmap,/over,mapvx,mapvy,mapbz=smap,limit=0.05,scale=9,iskip=9,jskip=9,$
            v_color=255,axis_color=0,/Nolabels,/Noaxis,v_thick=2.0,/sample,index_color=0,index_size=1,index_charsize=2
  write_png,DIR1+'mapvxybz'+strtrim(string(i,format='(I03)'),2)+'.png',tvrd()  
  plot_map,mapvz,color=0,bcolor=0,charsize=1.8
  plot_vmap,/over,mapvx,mapvy,mapbz=mapvz,limit=0.05,scale=9,iskip=9,jskip=9,$
            v_color=255,axis_color=0,/Nolabels,/Noaxis,v_thick=2.0,/sample,index_color=0,index_size=1,index_charsize=2
  write_png,DIR1+'mapvxyz'+strtrim(string(i,format='(I03)'),2)+'.png',tvrd()
endfor
!p.background=0
end

