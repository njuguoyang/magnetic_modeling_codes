;+
; NAME :
;   spine_fan_draw
; PURPOSE :
;   plot and save magnetic field lines in the vicinity of a null point
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

pro spine_fan_draw

!p.background=255
restore,'../../04extrapolation/Bout.sav',/ver
restore,'../../02ambiguity/03amb_removed/mapbxyz000.sav',/ver
smap = rebin_map(mapbz,256,256)

fits2map,'trf20020726_185826_a1.fits',maptr
maptr=shift_map(maptr,6,-4)
sub_map,maptr,smap2,xrange = [-420,-170],yrange=[-475,-315]
data = smap2.data
index = where(data le 100.0)
data[index] = 100.0
index = where(data ge 300.0)
data[index] = 300.0
smap2.data = data
;wdef,1,1000,600
;trace_colors,1700
;plot_map,smap2,color=0,charsize=2
;loadct,0

ss=size(smap.data)
x0=(0.0-0.5*(float(ss[1])-1.0))*smap.dx+smap.xc
y0=(0.0-0.5*(float(ss[2])-1.0))*smap.dy+smap.yc
sub_map,smap,smap3,irange=iran,ref_map=smap2
bx=twbox.bx[iran[0]:iran[1],iran[2]:iran[3],*]
by=twbox.by[iran[0]:iran[1],iran[2]:iran[3],*]
bz=twbox.bz[iran[0]:iran[1],iran[2]:iran[3],*]
help,bx,by,bz
print,'irange:',iran
x1=(iran[0]-0.5*(float(ss[1])-1.0))*smap.dx+smap.xc
y1=(iran[2]-0.5*(float(ss[2])-1.0))*smap.dy+smap.yc
xyz0=[x1,y1,0.0]
print,xyz0

r        = 3.0 
theta    = !dpi/4.0*3.0
nphi     = 15 
null_pos = [-10.002025138102224,       -2.0288121612976582,        8.7311354266446362]
evecx    = [0.50310834424143847,      -0.66892397229728084,       0.54719531544108146]
evecy    = [0.67321984602116669,       0.73903488690504537,      -2.45453633194421747E-002]
evecz    = [-0.35561227063902695,      -7.96773145625667423E-002,  0.93123114129363893]
xindex=0.0
yindex=0.0
zindex=0.0
spine_fan_pos, r, theta, nphi, null_pos, evecx, evecy, evecz, xindex=xindex, yindex=yindex, zindex=zindex
xindex=xindex+155.0+29.5
yindex=yindex+95.0+29.5
xindex=xindex*smap.dx+x0
yindex=yindex*smap.dy+y0
zindex=zindex*smap.dx
xfoot = (null_pos[0]+155.0+29.5)*smap.dx+x0
yfoot = (null_pos[1]+95.0+29.5)*smap.dy+y0
zfoot = null_pos[2]*smap.dx
print,'Position of the Null Point:',xfoot,yfoot,zfoot
fieldline3d,bx,by,bz,dx=smap.dx,dy=smap.dy,dz=smap.dx,xyz0=xyz0,line_thick=1,xindex=xindex,yindex=yindex,zindex=zindex,/both,/noimage,line_color=[0,255,0],/sav_fl
spawn,'mv curve*.sav ./inner_fl'

theta    = !dpi/4.0
nphi     = 11 
xindex=0.0
yindex=0.0
zindex=0.0
spine_fan_pos, r, theta, nphi, null_pos, evecx, evecy, evecz, xindex=xindex, yindex=yindex, zindex=zindex
xindex=xindex+155.0+29.5
yindex=yindex+95.0+29.5
xindex=xindex*smap.dx+x0
yindex=yindex*smap.dy+y0
zindex=zindex*smap.dx
fieldline3d,bx,by,bz,dx=smap.dx,dy=smap.dy,dz=smap.dx,xyz0=xyz0,line_thick=1,xindex=xindex,yindex=yindex,zindex=zindex,/both,line_color=[0,0,255],/noimage,/nocontour,/sav_fl
spawn,'mv curve*.sav ./outer_fl'

ss=size(smap2.data)
x0=(0.0-0.5*(float(ss[1])-1.0))*smap2.dx+smap2.xc
y0=(0.0-0.5*(float(ss[2])-1.0))*smap2.dy+smap2.yc
xyz0=[x0,y0,0.0]
trace_colors,1700
iimage,smap2.data,image_location=[xyz0[0],xyz0[1]],image_dimensions=[ss[1]*smap2.dx,ss[2]*smap2.dy],overplot=1,zvalue=xyz0[2]
loadct,0

end
