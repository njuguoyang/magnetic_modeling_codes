;+
; NAME :
;   faxis
; PURPOSE :
;   
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

pro faxis

lm=0.07 ;left margin
rm=0.01 ;right margin
bm=0.07 ;bottom margin
tm=0.01 ;top margin
cd=0.07 ;column distance
rd=0.07 ;row distance
col=3
row=2
cw=(1.0-lm-rm-(col-1)*cd)/col ;column width
rw=(1.0-bm-tm-(row-1)*rd)/row ;row width
!P.MULTI = [0, col, row]
device,decomposed=0
loadct,0

color1=[255,100,255]
color1='ff64ff'x
x1=52.64 & y1=67.96 & z1=1.921
color2=[240,240,0]
color2='00f0f0'x
x2=54.15 & y2=68.49 & z2=2.232
color3=[255,0,0]              ;;;;;;;;;;;;;;;red
color3='0000ff'x
x3=56.26 & y3=69.48 & z3=2.031
color4=[0,150,255]
color4='ff9600'x
x4=57.86 & y4=72.17 & z4=1.661
color5=[0,240,0]              ;;;;;;;;;;;;;;;green
color5='00f000'x
x5=53.34 & y5=69.96 & z5=0.5084
color6=[105,45,194]
color6='c22d69'x
x6=53.89 & y6=68.01 & z6=3.374
color7=[0,0,255]              ;;;;;;;;;;;;;;;blue
color7='ff0000'x
x7=54.71 & y7=68.25 & z7=0.4005
xindex=[x1,x2,x3,x4,x5,x6,x7]
yindex=[y1,y2,y3,y4,y5,y6,y7]
zindex=[z1,z2,z3,z4,z5,z6,z7]

filename_r='mapbxyz20050527_1017'
filename1='../20050527/01data_02ambiguity/20050527MTR/ambiguity_removal_current_flip/'+filename_r+'.sav'
restore,filename1,/verbose
bz0=mapbz.data
ss=size(bz0)
x1=165
x2=292
y1=83
y2=210
restore,'Bxyz_IDL.sav',/verbose
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
x00=(165.0-0.5*(ss[1]-1.0))*mapbz.dx+mapbz.xc   ;data coordinates of the low left corner
y00=(83.0-0.5*(ss[2]-1.0))*mapbz.dy+mapbz.yc
xindex=(xindex+x1-0.5*(ss[1]-1))*mapbz.dx+mapbz.xc
yindex=(yindex+y1-0.5*(ss[2]-1))*mapbz.dy+mapbz.yc
zindex=zindex*mapbz.dx

;================================
;find z-axis by determing the position where Bv changes its sign
;================================
;xm=(x0-x00)/mapbz.dx 		;convert to pixel position
;ym=(y0-y00)/mapbz.dy		;convert to pixel position
;n1n=floor(xm)
;n2n=floor(ym)
;dx1=xm-n1n
;dy1=ym-n2n
;xv=fltarr(3,8)
;nz=5
;zm=findgen(nz)
;Bxinterp=fltarr(nz)
;Byinterp=fltarr(nz)
;Bzinterp=fltarr(nz)
;for i=0,nz-1 do begin
;  n3n=floor(zm[i])
;  dz1=zm[i]-n3n
;  print,'n1n,n2n,n3n,dx1,dy1,dz1:',n1n,n2n,n3n,dx1,dy1,dz1
;  CORNER,n1n,n2n,n3n,Bx,By,Bz,xv
;  Bxinterp[i]=XITP(0,dx1,dy1,dz1,xv)
;  Byinterp[i]=XITP(1,dx1,dy1,dz1,xv)
;  Bzinterp[i]=XITP(2,dx1,dy1,dz1,xv)
;endfor
;theta=atan(k0)
;print,'Theta:',theta*!radeg
;bp=byinterp*sin(theta)+bxinterp*cos(theta)
;bv=byinterp*cos(theta)-bxinterp*sin(theta)
;bzz=bzinterp
;xaxis=findgen(nz)*mapbz.dx
;coeff=POLY_FIT(xaxis,bv,1,/double)
;zaxis=-coeff[0]/coeff[1]
;print,'zaxis:',zaxis
;xindex=[xindex,x0]
;yindex=[yindex,y0]
;zindex=[zindex,zaxis]

;i=1
;j=0
;x0p=lm+(cw+cd)*j
;y0p=bm+(rw+rd)*(row-1-i)
;x1p=lm+(cw+cd)*j+cw
;y1p=bm+(rw+rd)*(row-1-i)+rw
;plot,xaxis,bp,pos=[x0p,y0p,x1p,y1p],charsize=1,psym=6,xrange=[-0.27,2.08],yrange=[-450,950],$
;             color=0,background=255,symsize=1.0,xstyle=1,ystyle=1,xminor=1,yminor=2,xtitle='Height (arcsec)'
;oplot,[-0.27,2.08],[0.0,0.0],color=0,linestyle=1
;oplot,xaxis,bv,psym=5,color=0,symsize=1.0
;oplot,xaxis,bzz,psym=2,color=0,symsize=1.0

;================================
;find z-axis by 'findaxis.pro'
;================================
restore,'./findaxis/sample_start_point.sav',/verbose
xindex=[xindex,xindex1[axispos]]        ;x-, y-, z-index at "axispos" are the coordinates where the axis passes by. 
yindex=[yindex,yindex1[axispos]]
zindex=[zindex,zindex1[axispos]]

;step_size=0.001
;fieldline3d2,bx,by,bz,xyz0=[x00,y00,0.0],dx=mapbz.dx,dy=mapbz.dy,dz=mapbz.dx,xindex=xindex,yindex=yindex,zindex=zindex,/both,/noplot,/nocontour,/noimage,step_size=step_size,isn=60000
;fieldline3d,bx,by,bz,xyz0=[x00,y00,0.0],dx=mapbz.dx,dy=mapbz.dy,dz=mapbz.dx,xindex=xindex,yindex=yindex,zindex=zindex,/both
;================================
;
;================================
restore,'curve_step_size0.02/curve_info0000.sav',/ver
;window,/free
;plot,(tangentx^2+tangenty^2+tangentz^2)^0.5,xstyle=1,ystyle=1,yrange=[0.99,1.01]
;window,/free
;linedx=linex-shift(linex,1) 
;linedx[0]=linedx[1]
;linedy=liney-shift(liney,1) 
;linedy[0]=linedy[1]
;linedz=linez-shift(linez,1) 
;linedz[0]=linedz[1]
;plot,(linedx^2+linedy^2+linedz^2)^0.5,xstyle=1,zstyle=1,yrange=[1.0/16.0-0.01,1.0/16.0+0.01]
tangentdx=deriv(tangentx)
tangentdy=deriv(tangenty)
tangentdz=deriv(tangentz)
kappa1=(tangentdx^2+tangentdy^2+tangentdz^2)^0.5

restore,'curve_step_size0.02/curve_info0001.sav'
tangentdx=deriv(tangentx)
tangentdy=deriv(tangenty)
tangentdz=deriv(tangentz)
kappa2=(tangentdx^2+tangentdy^2+tangentdz^2)^0.5

restore,'curve_step_size0.02/curve_info0002.sav'
linex3=linex & liney3=liney & linez3=linez
tangentdx=deriv(tangentx)
tangentdy=deriv(tangenty)
tangentdz=deriv(tangentz)
kappa3=(tangentdx^2+tangentdy^2+tangentdz^2)^0.5

restore,'curve_step_size0.02/curve_info0003.sav'
tangentdx=deriv(tangentx)
tangentdy=deriv(tangenty)
tangentdz=deriv(tangentz)
kappa4=(tangentdx^2+tangentdy^2+tangentdz^2)^0.5

restore,'curve_step_size0.02/curve_info0004.sav'
linex5=linex & liney5=liney & linez5=linez
tangentdx=deriv(tangentx)
tangentdy=deriv(tangenty)
tangentdz=deriv(tangentz)
kappa5=(tangentdx^2+tangentdy^2+tangentdz^2)^0.5

restore,'curve_step_size0.02/curve_info0005.sav'
tangentdx=deriv(tangentx)
tangentdy=deriv(tangenty)
tangentdz=deriv(tangentz)
kappa6=(tangentdx^2+tangentdy^2+tangentdz^2)^0.5

restore,'curve_step_size0.02/curve_info0006.sav'
linex7=linex & liney7=liney & linez7=linez
tangentdx=deriv(tangentx)
tangentdy=deriv(tangenty)
tangentdz=deriv(tangentz)
kappa7=(tangentdx^2+tangentdy^2+tangentdz^2)^0.5

restore,'curve_step_size0.02/curve_info0007.sav'
linex8=linex & liney8=liney & linez8=linez
tangentdx=deriv(tangentx)
tangentdy=deriv(tangenty)
tangentdz=deriv(tangentz)
kappa8=(tangentdx^2+tangentdy^2+tangentdz^2)^0.5

;================================
;fig1
;================================
window,1,xs=1500,ys=1000
device,decomposed=0
loadct,0
i=0
j=0
x0p=lm+(cw+cd)*j
y0p=bm+(rw+rd)*(row-1-i)
x1p=lm+(cw+cd)*j+cw
y1p=bm+(rw+rd)*(row-1-i)+rw
plot,pathx,pathy,xrange=xrange,yrange=yrange,xstyle=1,ystyle=1,linestyle=3,$
     thick=2,xcharsize=4,ycharsize=4,background=255,color=0,xtitle='x (arcsec)',$
     ytitle='y (arcsec)',pos=[x0p,y0p,x1p,y1p]

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
oplot,[pathx[ind1],pathx[ind2]],[pathy[ind1],pathy[ind2]],psym=2,symsize=3,color=0,thick=1

coeff=POLY_FIT(pathx[ind1:ind2],pathy[ind1:ind2],3,/double)
nn=100
xx=dindgen(nn)/99.0*(pathx[ind2]+8-pathx[ind1]+8)+pathx[ind1]-8
yy=xx^3*coeff[3]+xx^2*coeff[2]+xx*coeff[1]+coeff[0]
coeff3=coeff
oplot,xx[20:85],yy[20:85],color=0,thick=2

x0=pathx[(ind1+ind2)/2]
y0=pathy[(ind1+ind2)/2]
oplot,[x0],[y0],psym=2,symsize=3,color=0,thick=1
xyouts,-111,-128,'1',charsize=2.0,color=0
xyouts,-108,-124,'2',charsize=2.0,color=0
xyouts,-105,-120,'3',charsize=2.0,color=0

k0=x0^2*coeff3[3]*3+x0*coeff3[2]*2+coeff3[1]
theta=atan(k0)
theta1=theta+!pi/2.0
npoint=2
smax=5.0
smin=-5.0
xindex=x0-(findgen(npoint)-0.5*(float(npoint)-1))*(smax-smin)/(float(npoint)-1)*cos(theta1)
yindex=y0-(findgen(npoint)-0.5*(float(npoint)-1))*(smax-smin)/(float(npoint)-1)*sin(theta1)
oplot,xindex,yindex,color=0,thick=2

;================================
;fig2
;================================
i=0
j=1
x0p=lm+(cw+cd)*j
y0p=bm+(rw+rd)*(row-1-i)
x1p=lm+(cw+cd)*j+cw
y1p=bm+(rw+rd)*(row-1-i)+rw
plot,pathx,pathy,xrange=[-118,-89],yrange=[-138,-109],xstyle=1,ystyle=1,linestyle=3,$
     thick=2,xcharsize=4,ycharsize=4,color=0,xtitle='x (arcsec)',$
     ytitle='y (arcsec)',pos=[x0p,y0p,x1p,y1p]
oplot,[pathx[ind1],pathx[ind2]],[pathy[ind1],pathy[ind2]],psym=2,symsize=3,color=0,thick=1
device,decomposed=1
oplot,linex3,liney3,color=color3,thick=2
oplot,linex5,liney5,color=color5,thick=2
oplot,linex7,liney7,color=color7,thick=2
oplot,linex8,liney8,color='000000'x,thick=2

;================================
;fig3
;================================
i=0
j=2
x0p=lm+(cw+cd)*j
y0p=bm+(rw+rd)*(row-1-i)
x1p=lm+(cw+cd)*j+cw
y1p=bm+(rw+rd)*(row-1-i)+rw
plot,[0],xrange=[-118,-89],yrange=[-0.1,2.9],xstyle=1,ystyle=1,color=0,xcharsize=4,ycharsize=4,$
     xtitle='x (arcsec)',ytitle='z (arcsec)',pos=[x0p,y0p,x1p,y1p]
oplot,linex3,linez3,color=color3,thick=2
oplot,linex5,linez5,color=color5,thick=2
oplot,linex7,linez7,color=color7,thick=2
oplot,linex8,linez8,color='000000'x,thick=2

;================================
;fig4
;================================
i=1
j=0
x0p=lm+(cw+cd)*j
y0p=bm+(rw+rd)*(row-1-i)
x1p=lm+(cw+cd)*j+cw
y1p=bm+(rw+rd)*(row-1-i)+rw
plot,[0],xrange=[pathx[ind1]-0.3,pathx[ind2]+0.3],yrange=[-0.1,2.9],xstyle=1,ystyle=1,color='000000'x,$
     xcharsize=4,ycharsize=4,xtitle='x (arcsec)',ytitle='!9 ! !3y-y!ia!n!9 ! !3 (arcsec)',pos=[x0p,y0p,x1p,y1p]

xerr=1.0/16.0*mapbz.dx
yerr=8.0*mapbz.dy

xind1=where(linex3 gt pathx[ind1]-xerr AND linex3 lt pathx[ind1]+xerr)
yind1=where(liney3[xind1] gt pathy[ind1]-yerr AND liney3[xind1] lt pathy[ind1]+yerr)
indl1=xind1[yind1[0]]
xind2=where(linex3 gt pathx[ind2]-xerr AND linex3 lt pathx[ind2]+xerr)
yind2=where(liney3[xind2] gt pathy[ind2]-yerr AND liney3[xind2] lt pathy[ind2]+yerr)
indl2=xind2[yind2[0]]
print,'indl1,indl2:',indl1,indl2
linex3=linex3[indl1:indl2]
liney3=liney3[indl1:indl2]
linez3=linez3[indl1:indl2]
pil3=linex3^3*coeff[3]+linex3^2*coeff[2]+linex3*coeff[1]+coeff[0]
diff3=abs(liney3-pil3)
kappa3=kappa3[indl1:indl2]

xind1=where(linex5 gt pathx[ind1]-xerr AND linex5 lt pathx[ind1]+xerr)
yind1=where(liney5[xind1] gt pathy[ind1]-yerr AND liney5[xind1] lt pathy[ind1]+yerr)
indl1=xind1[yind1[0]]
xind2=where(linex5 gt pathx[ind2]-xerr AND linex5 lt pathx[ind2]+xerr)
yind2=where(liney5[xind2] gt pathy[ind2]-yerr AND liney5[xind2] lt pathy[ind2]+yerr)
indl2=xind2[yind2[0]]
print,'indl1,indl2:',indl1,indl2
linex5=linex5[indl1:indl2]
liney5=liney5[indl1:indl2]
linez5=linez5[indl1:indl2]
pil5=linex5^3*coeff[3]+linex5^2*coeff[2]+linex5*coeff[1]+coeff[0]
diff5=abs(liney5-pil5)
kappa5=kappa5[indl1:indl2]

xind1=where(linex7 gt pathx[ind1]-xerr AND linex7 lt pathx[ind1]+xerr)
yind1=where(liney7[xind1] gt pathy[ind1]-yerr AND liney7[xind1] lt pathy[ind1]+yerr)
indl1=xind1[yind1[0]]
xind2=where(linex7 gt pathx[ind2]-xerr AND linex7 lt pathx[ind2]+xerr)
yind2=where(liney7[xind2] gt pathy[ind2]-yerr AND liney7[xind2] lt pathy[ind2]+yerr)
indl2=xind2[yind2[0]]
print,'indl1,indl2:',indl1,indl2
linex7=linex7[indl1:indl2]
liney7=liney7[indl1:indl2]
linez7=linez7[indl1:indl2]
pil7=linex7^3*coeff[3]+linex7^2*coeff[2]+linex7*coeff[1]+coeff[0]
diff7=abs(liney7-pil7)
kappa7=kappa7[indl1:indl2]

xind1=where(linex8 gt pathx[ind1]-xerr AND linex8 lt pathx[ind1]+xerr)
yind1=where(liney8[xind1] gt pathy[ind1]-yerr AND liney8[xind1] lt pathy[ind1]+yerr)
indl1=xind1[yind1[0]]
xind2=where(linex8 gt pathx[ind2]-xerr AND linex8 lt pathx[ind2]+xerr)
yind2=where(liney8[xind2] gt pathy[ind2]-yerr AND liney8[xind2] lt pathy[ind2]+yerr)
indl2=xind2[yind2[0]]
print,'indl1,indl2:',indl1,indl2
linex8=linex8[indl1:indl2]
liney8=liney8[indl1:indl2]
linez8=linez8[indl1:indl2]
pil8=linex8^3*coeff[3]+linex8^2*coeff[2]+linex8*coeff[1]+coeff[0]
diff8=abs(liney8-pil8)
kappa8=kappa8[indl1:indl2]

oplot,linex3,diff3,color=color3,thick=2
oplot,linex5,diff5,color=color5,thick=2
oplot,linex7,diff7,color=color7,thick=2
oplot,linex8,diff8,color=0,thick=2

;================================
;fig5
;================================
i=1
j=1
x0p=lm+(cw+cd)*j
y0p=bm+(rw+rd)*(row-1-i)
x1p=lm+(cw+cd)*j+cw
y1p=bm+(rw+rd)*(row-1-i)+rw
plot,[0],xrange=[pathx[ind1]-0.3,pathx[ind2]+0.3],yrange=[-0.1,0.9],xstyle=1,ystyle=1,color='000000'x,$
     xcharsize=4,ycharsize=4,xtitle='x (arcsec)',ytitle='!9 ! !3z-<z>!9 ! !3 (arcsec)',pos=[x0p,y0p,x1p,y1p]
oplot,linex3,abs(linez3-mean(linez3)),color=color3,thick=2
oplot,linex5,abs(linez5-mean(linez5)),color=color5,thick=2
oplot,linex7,abs(linez7-mean(linez7)),color=color7,thick=2
oplot,linex8,abs(linez8-mean(linez8)),color=0,thick=2

;================================
;fig6
;================================
result=pb0r('27-May-2005 10:17:00.000')
mpa=695.5/result[2]/60.0  ;Mm per arcsec. 695.5Mn is the radius of the Sun
ds=1.0/16.0*mapbz.dx*mpa  ;the default value
;ds=step_size*mpa
kappa3=kappa3/ds
kappa5=kappa5/ds
kappa7=kappa7/ds
kappa8=kappa8/ds

linedx=linex3-shift(linex3,1) 
linedx[0]=linedx[1]
linedy=liney3-shift(liney3,1) 
linedy[0]=linedy[1]
linedz=linez3-shift(linez3,1) 
linedz[0]=linedz[1]
arcld3=(linedx^2+linedy^2+linedz^2)^0.5
arcl3=fltarr(n_elements(kappa3))
arcl3[0]=0.0
for i=1,n_elements(kappa3)-1 do arcl3[i]=arcl3[i-1]+arcld3[i]

linedx=linex5-shift(linex5,1) 
linedx[0]=linedx[1]
linedy=liney5-shift(liney5,1) 
linedy[0]=linedy[1]
linedz=linez5-shift(linez5,1) 
linedz[0]=linedz[1]
arcld5=(linedx^2+linedy^2+linedz^2)^0.5
arcl5=fltarr(n_elements(kappa5))
arcl5[0]=0.0
for i=1,n_elements(kappa5)-1 do arcl5[i]=arcl5[i-1]+arcld5[i]

linedx=linex7-shift(linex7,1) 
linedx[0]=linedx[1]
linedy=liney7-shift(liney7,1) 
linedy[0]=linedy[1]
linedz=linez7-shift(linez7,1) 
linedz[0]=linedz[1]
arcld7=(linedx^2+linedy^2+linedz^2)^0.5
arcl7=fltarr(n_elements(kappa7))
arcl7[0]=0.0
for i=1,n_elements(kappa7)-1 do arcl7[i]=arcl7[i-1]+arcld7[i]

linedx=linex8-shift(linex8,1) 
linedx[0]=linedx[1]
linedy=liney8-shift(liney8,1) 
linedy[0]=linedy[1]
linedz=linez8-shift(linez8,1) 
linedz[0]=linedz[1]
arcld8=(linedx^2+linedy^2+linedz^2)^0.5
arcl8=fltarr(n_elements(kappa8))
arcl8[0]=0.0
for i=1,n_elements(kappa8)-1 do arcl8[i]=arcl8[i-1]+arcld8[i]

ymax=max([kappa3,kappa5,kappa7,kappa8])
xmax=max([arcl3,arcl5,arcl7,arcl8])
print,'average of kappa3,5,7,8:',mean(kappa3),mean(kappa5),mean(kappa7),mean(kappa8)

i=1
j=2
x0p=lm+(cw+cd)*j
y0p=bm+(rw+rd)*(row-1-i)
x1p=lm+(cw+cd)*j+cw
y1p=bm+(rw+rd)*(row-1-i)+rw
plot,[0],xstyle=1,xrange=[0,xmax],ystyle=1,yrange=[0,ymax],color='000000'x,$
     xtitle='Arclength (arcsec)',ytitle='Curvature (Mm!e-1!n)',xcharsize=4,ycharsize=4,$
     pos=[x0p,y0p,x1p,y1p]

kappa3=smooth(kappa3,5)
kappa5=smooth(kappa5,5)
kappa7=smooth(kappa7,5)
kappa8=smooth(kappa8,5)
oplot,arcl3,kappa3,color=color3,thick=2
oplot,arcl5,kappa5,color=color5,thick=2
oplot,arcl7,kappa7,color=color7,thick=2
oplot,arcl8,kappa8,color=0,thick=2
xyouts,0.79,0.45,'Average Curvature:',/normal,charsize=2,color=0
xyouts,0.92,0.45,strtrim(string(mean(kappa5),format='(f4.2)'),2)+' Mm!e-1!n',/normal,charsize=2,color=color5
xyouts,0.92,0.425,strtrim(string(mean(kappa7),format='(f4.2)'),2),/normal,charsize=2,color=color7
xyouts,0.92,0.40,strtrim(string(mean(kappa3),format='(f4.2)'),2),/normal,charsize=2,color=color3
xyouts,0.92,0.375,strtrim(string(mean(kappa8),format='(f4.2)'),2),/normal,charsize=2,color=0

xyouts,0.86,0.95,'Twist:',/normal,charsize=2,color=0
xyouts,0.90,0.95,'-1.36 turns',/normal,charsize=2,color=color5
xyouts,0.90,0.925,'-1.61',/normal,charsize=2,color=color7
xyouts,0.90,0.90,'-1.92',/normal,charsize=2,color=color3

write_png,'axisf3.png',tvrd(true=1)
device,decomposed=0
!P.MULTI = 0
end

