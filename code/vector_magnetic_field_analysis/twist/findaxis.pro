;+
; NAME :
;   findaxis
; PURPOSE :
;   To determine the axis of the flux rope found on 2005 May 27, 10:17 UT
; CATEGORY :
;
; CALLING SEQUENCE :
;   
; INPUTS :
;
; OUTPUTS :
;
; COMMON BLOCKS :
;
; MODIFICATION HISTORY :
;   2010.02 Guo Yang @ Nanjing University, 
;
;   This software is provided without any warranty. Permission to use,
;   copy, modify. Distributing modified or unmodified copies is granted,
;   provided this disclaimer and information are included unchanged.
;-

pro findaxis

filename_r='mapbxyz20050527_1017'
filename1='../../20050527/01data_02ambiguity/20050527MTR/ambiguity_removal_current_flip/'+filename_r+'.sav'
restore,filename1,/verbose
bz0=mapbz.data
ss=size(bz0)
x1=165
x2=292
y1=83
y2=210
restore,'../Bxyz_IDL.sav',/verbose
ssbz=size(bz)
sx=ssbz[1]
sy=ssbz[2]
contour,bz[*,*,0],levels=[0.0],position=[0.0,0.0,1.0,1.0],xstyle=1,ystyle=1,path_xy=pathxy,closed=0,/path_double
pathx=(pathxy[0,*]*sx+x1-0.5*(ss[1]-1))*mapbz.dx+mapbz.xc
pathy=(pathxy[1,*]*sy+y1-0.5*(ss[2]-1))*mapbz.dy+mapbz.yc
xrange=[(x1-0.5*(ss[1]-1))*mapbz.dx+mapbz.xc,(x2-0.5*(ss[1]-1))*mapbz.dx+mapbz.xc]
yrange=[(y1-0.5*(ss[2]-1))*mapbz.dy+mapbz.yc,(y2-0.5*(ss[2]-1))*mapbz.dy+mapbz.yc]
pathx=reform(pathx)
pathy=reform(pathy)
xp1=(45.20+x1-0.5*(ss[1]-1))*mapbz.dx+mapbz.xc
xp2=(47.20+x1-0.5*(ss[1]-1))*mapbz.dx+mapbz.xc
xind1=where(pathx ge xp1 AND pathx le xp2)
yp1=(57.16+y1-0.5*(ss[2]-1))*mapbz.dy+mapbz.yc
yp2=(59.16+y1-0.5*(ss[2]-1))*mapbz.dy+mapbz.yc
yind1=where(pathy[xind1] ge yp1 AND pathy[xind1] le yp2)
ind1=xind1[yind1[0]]
print,'Start point of the polarity inversion line:',pathx[ind1],pathy[ind1]
xp1=(57.61+x1-0.5*(ss[1]-1))*mapbz.dx+mapbz.xc
xp2=(59.61+x1-0.5*(ss[1]-1))*mapbz.dx+mapbz.xc
xind2=where(pathx ge xp1 AND pathx le xp2)
yp1=(72.90+y1-0.5*(ss[2]-1))*mapbz.dy+mapbz.yc
yp2=(74.90+y1-0.5*(ss[2]-1))*mapbz.dy+mapbz.yc
yind2=where(pathy[xind2] ge yp1 AND pathy[xind2] le yp2)
ind2=xind2[yind2[0]]
print,'Stop point of the polarity inversion line:',pathx[ind2],pathy[ind2]
coeff=POLY_FIT(pathx[ind1:ind2],pathy[ind1:ind2],3,/double)
x0=pathx[(ind1+ind2)/2]
y0=pathy[(ind1+ind2)/2]
k0=x0^2*coeff[3]*3+x0*coeff[2]*2+coeff[1]
theta=atan(k0)
theta1=theta+!pi/2.0
npoint=80
smax=1.0
smin=-1.0
xindex=x0-(findgen(npoint)-0.5*(float(npoint)-1))*(smax-smin)/(float(npoint)-1)*cos(theta1)
yindex=y0-(findgen(npoint)-0.5*(float(npoint)-1))*(smax-smin)/(float(npoint)-1)*sin(theta1)
zindex=findgen(npoint)*(smax-smin)/(float(npoint)-1)
print,'xindex:',xindex
print,'yindex:',yindex
print,'zindex:',zindex

xindex1=xindex
yindex1=yindex
zindex1=zindex
for i=0,npoint-2 do begin
  xindex1=[[xindex1],[xindex]]
  yindex1=[[yindex1],[yindex]]
  zindex1=[[zindex1],[zindex]]
endfor
zindex1=transpose(zindex1)
xindex1=reform(xindex1,npoint*npoint)
yindex1=reform(yindex1,npoint*npoint)
zindex1=reform(zindex1,npoint*npoint)
;iplot,xindex1,yindex1,zindex1,sym_index=4,/scatter

x00=(165.0-0.5*(ss[1]-1.0))*mapbz.dx+mapbz.xc   ;data coordinates of the low left corner
y00=(83.0-0.5*(ss[2]-1.0))*mapbz.dy+mapbz.yc
fieldline3d2,bx,by,bz,xyz0=[x00,y00,0.0],dx=mapbz.dx,dy=mapbz.dy,dz=mapbz.dx,xindex=xindex1,yindex=yindex1,zindex=zindex1,/both,/noplot,/nocontour,/noimage

xerr=1.0/16.0*mapbz.dx
yerr=4.0*mapbz.dy
device,decomposed=0
loadct,0
window,1,xs=1000,ys=1000
plot,pathx,pathy,color=0,background=255,xstyle=1,xrange=[-135,-75],ystyle=1,yrange=[-155,-95]
oplot,[pathx[ind1],pathx[ind2]],[pathy[ind1],pathy[ind2]],psym=2,symsize=3,color=0,thick=1

linef=[9999L]
objf=[10.d0^9]
for k=0L,long(npoint*npoint-1) do begin
  filename='curve_info'+string(k,format='(i4.4)')+'.sav'
  restore,filename
  xind1=where(linex gt pathx[ind1]-xerr AND linex lt pathx[ind1]+xerr)
  if (xind1[0] ne -1) then begin
    yind1=where(liney[xind1] gt pathy[ind1]-yerr AND liney[xind1] lt pathy[ind1]+yerr)
    if (yind1[0] ne -1) then begin
      indl1=xind1[yind1[0]]
    endif else begin
      ;print,'k:',k,'  Ymin does not arrive.'
      ;oplot,linex,liney,linestyle=1,color=0
      ;wait,3
      continue
    endelse
  endif else begin
    ;print,'k:',k,'  Xmin does not arrive.'
    ;oplot,linex,liney,linestyle=2,color=0
    ;wait,3
    continue
  endelse
  xind2=where(linex gt pathx[ind2]-xerr AND linex lt pathx[ind2]+xerr)
  if (xind2[0] ne -1) then begin
    yind2=where(liney[xind2] gt pathy[ind2]-yerr AND liney[xind2] lt pathy[ind2]+yerr)
    if (yind2[0] ne -1) then begin
      indl2=xind2[yind2[0]]
    endif else begin
      ;print,'k:',k,'  Ymax does not arrive.'
      ;oplot,linex,liney,linestyle=3,color=0
      ;wait,3
      continue
    endelse
  endif else begin
   ; print,'k:',k,'  Xmax does not arrive.'
   ; oplot,linex,liney,linestyle=4,color=0
   ; wait,3
    continue
  endelse
  ;oplot,linex,liney,linestyle=0,color=0
  linex1=linex[indl1:indl2]
  liney1=liney[indl1:indl2]
  linez1=linez[indl1:indl2]
  pil1=linex1^3*coeff[3]+linex1^2*coeff[2]+linex1*coeff[1]+coeff[0]
  diff1=abs(liney1-pil1)
  ;plot,linex1,liney1,color=0,background=255
  ;oplot,linex1,pil1,color=0,linestyle=4
  ;window,2
  ;plot,linex1,diff1,color=0,linestyle=2,background=255
  area1=INT_TABULATED(linex1,diff1,/double)
  diff2=abs(linez1-mean(linez1))
  area2=INT_TABULATED(linex1,diff2,/double)
  objf1=area1+area2
  ;oplot,linex1,liney1,psym=3,symsize=1,color=0
  ;print,'k:',k,'  area:',area1,'  standard deviation of heights:',stddev1
  ;wait,2
  linef=[linef,k]      ;The No. of the file that store the magnetic field line
  objf=[objf,objf1]    ;The objective function
endfor
indobj=where(objf eq min(objf))
axispos=linef(indobj)
print,'Axis postion:',axispos
filename='curve_info'+string(axispos,format='(i4.4)')+'.sav'
restore,filename
loadct,3
oplot,linex,liney,color=128
window,4
plot,linex,linez,color=0,linestyle=2,background=255
save,xindex1,yindex1,zindex1,axispos,filename='sample_start_point.sav'

spawn,'rm curve_info*.sav'


end
