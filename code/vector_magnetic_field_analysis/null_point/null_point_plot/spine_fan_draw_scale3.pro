;+
; NAME :
;   spine_fan_draw_scale3
; PURPOSE :
;   plot magnetic field lines in the vicinity of a null point
; CATEGORY :
; CALLING SEQUENCE :
; INPUTS :
; OPTIONAL INPUT PARAMETERS (KEYWORD PARAMETERS):
; OUTPUTS :
; COMMON BLOCKS :
; SIDE EFFECTS :
; RESTRICTIONS :
; PROCEDURE :
; MODIFICATION HISTORY :
;   2012.11 Guo Yang@Nanjing University
;
;   This software is provided without any warranty. Permission to use,
;   copy, modify. Distributing modified or unmodified copies is granted,
;   provided this disclaimer and information are included unchanged.
;-

pro spine_fan_draw_scale3

wdef,1,1200,800
device,decomposed=0
loadct,0
!p.background=255

fits2map,'trf20020726_185826_a1.fits',maptr
maptr=shift_map(maptr,6,-4)
sub_map,maptr,smap2,xrange = [-420,-170],yrange=[-475,-315]
data = smap2.data
index = where(data le 100.0)
data[index] = 100.0
index = where(data ge 300.0)
data[index] = 300.0
smap2.data = data

restore,'../../02ambiguity/03amb_removed/mapbxyz000.sav',/ver
smap = rebin_map(mapbz,256,256)
sub_map,smap,smap3,irange=iran,ref_map=smap2


xr1=-420
xr2=-170
yr1=-475
yr2=-315
zr1=0
zr2=200
;scale3,xrange=[xr1,xr2],yrange=[yr1,yr2],zrange=[zr1,zr2],ax=40,az=-50
;t3d,trans=[-0.20,0.00,0.0],scale=[1.4,1.4,1.4],perspective=3.0
scale3,xrange=[xr1,xr2],yrange=[yr1,yr2],zrange=[zr1,zr2],ax=90,az=0
t3d,trans=[-0.15,-0.20,0.0],scale=[1.4,1.4,1.4],perspective=3.0

;底部的contour/image
plot_map,smap2,color=0,charsize=2.5,/t3d,/no_data,/normal,ticklen=-0.02
;trace_colors,1700
loadct,10
ss = size(smap2.data)
xs = ss[1]
ys = ss[2]
x0=(0.0-0.5*(float(ss[1])-1.0))*smap2.dx+smap2.xc
y0=(0.0-0.5*(float(ss[2])-1.0))*smap2.dy+smap2.yc
xcoef = 1.0005
ycoef = 1.002
x=findgen(xs)*smap2.dx*xcoef+x0+smap2.dx*0.5
y=findgen(ys)*smap2.dy*ycoef+y0-smap2.dy
anzc=50
maxi = max(smap2.data)
mini = min(smap2.data)
farbvec=findgen(anzc)*(maxi-mini)/(anzc-1)+(mini)
contour, smap2.data, x, y, lev=farbvec, /fill, nlevel=anzc, charsize=2,  zvalue=0, xstyle=1, ystyle=1, /t3d, color=0 ,/over, /noclip
loadct,0

;画三维磁力线
loadct,13
filename = find_file('./inner_fl/curve*.sav',count=count)
for i=0,count-1 do begin
  restore,filename[i]
  plots, linex, liney, linez, /t3d, color=0, thick=1
endfor
filename = find_file('./outer_fl/curve*.sav',count=count)
for i=0,count-1 do begin
  restore,filename[i]
  plots, linex, liney, linez, /t3d, color=255, thick=1
endfor
write_png,'spine_fan.png',tvrd(true=1)
!p.background=0

end
