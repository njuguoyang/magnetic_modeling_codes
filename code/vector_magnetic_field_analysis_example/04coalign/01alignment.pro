;+
; NAME :
;   01alignment
; PURPOSE :
;   Align one BBSO DVMG Bz image with the nearest MDI line-of-sight magnetogram.
;-

;pro 01alignment

file=find_file('/home/guoyang/study/20050115/chengxin/bbsodvmg_20050115/bbso_dmgbz*')
fits2map,file[0],map1
filemdi='fd_M_96m_01d.4397.0011.fits'
fits2map,filemdi,map2
map2=map2earth(map2,/remap)
coalign,map1,map2

mapbz=map1
mapbz.ID='BBSO 10" DMG 6103!sA!r!u!9 %!3!n'
save,mapbz,filename='mapbz20050115_165235.sav'
!p.background=255
wdef,0,800,800
plot_map,mapbz,color=0,charsize=2.0
write_png,'mapbz20050105_165235.png',tvrd()

filex=find_file('/home/guoyang/study/20050115/chengxin/bbsodvmg_20050115/bbso_dmgbx*')
fits2map,filex[0],mapbx
mapbx.xc=mapbz.xc
mapbx.yc=mapbz.yc
filey=find_file('/home/guoyang/study/20050115/chengxin/bbsodvmg_20050115/bbso_dmgby*')
fits2map,filey[0],mapby
mapby.xc=mapbz.xc
mapby.yc=mapbz.yc
btrans=sqrt(mapbx.data^2+mapby.data^2)
bazim=atan(mapby.data,mapbx.data)/!dtor
mapbx.data=btrans*cos((bazim-22.)*!dtor)
mapby.data=btrans*sin((bazim-22.)*!dtor)
wdef,1,800,800
plot_map,mapbz,color=0,bcolor=0,charsize=2.0
plot_vmap,/over,mapbx,mapby,mapbz=mapbz,limit=30,scale=0.01,iskip=9,jskip=9,$
          v_color=255,axis_color=0,/Nolabels,/Noaxis,v_thick=2.0
write_png,'mapbxyz20050115_165235.png',tvrd()

!p.background=0
end
