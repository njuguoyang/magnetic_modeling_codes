;+
; NAME :
;   current_flip
;  PURPOSE:
;   Flip vectors in the selected subregions to remove the ambiguity
; CATEGORY :
;
; CALLING SEQUENCE :
;
; PROCEDURE :
;   1. Display the vertical current map
;   2. Select a subregion (Region 1) to analyze
;   3. If zoom in, select a subregion (Region 2) and then define an irregular region
;      and a direction by click a start and an end point
;   4. If not zoom in, define an irregular region in the subregion 1. The vectors
;      in this irregular region will be flipped
; MODIFICATION HISTORY :
;   2009.06 Guo Yang@Observatoire de Paris & Nanjing University, introduced
;   This software is provided without any warranty. Permission to use,
;   copy, modify. Distributing modified or unmodified copies is granted,
;   provided this disclaimer and information are included unchanged.
;-

pro current_flip,mapbx,mapby,mapbz

bxo=mapbx.data
byo=mapby.data
bzo=mapbz.data

;=====================================
;Calculate the line-of-sight current
;=====================================
;lamda is the pixel size in m
lamda=3.3d5

mu0=!dpi*4.d-3  ;The unit of mu0 is G*m/A
cc=0.0001       ;An arbitrary constant to control the contrast of current
t1=(1./(2.*lamda))*(shift(Byo,-1,0)-shift(Byo,1,0))
t2=(1./(2.*lamda))*(shift(Bxo,0,-1)-shift(Bxo,0,1))
Jzo=(1./mu0)*(t1-t2)
;current=-alog(abs(Jzo)/(abs(bzo)+cc))
current=max(abs(Jzo))-abs(Jzo)          ;to reverse the color
print,'max, min of currents:',max(current),min(current)
!p.background=255
ss=size(current)
xs=ss[1]
ys=ss[2]

scl0=2
rec0='Y'
wdef,0,scl0*xs,scl0*ys
while (rec0 eq 'Y' or rec0 eq 'y') do begin
  print,'The present scale factor is ',strtrim(string(scl0),2),'. Do you want to change it (Y/N)?'
  read,rec0,prompt='Input Y/N:  '
  if (rec0 eq 'Y' or rec0 eq 'y') then begin
    read,scl0,prompt='Input scale factor (should be integers):  '
    wdef,0,scl0*xs,scl0*ys
  endif
endwhile
;====================================
;To establishe Data coordinate system
;because FCCNVE.PRO needs all the plot information
;====================================
plot,[0,scl0*xs],[0,scl0*ys],xstyle=1,ystyle=1,position=[0,0,scl0*xs,scl0*ys]
c0=congrid(current,scl0*xs,scl0*ys,/center,cubic=-0.5)
bx0=congrid(mapbx.data,scl0*xs,scl0*ys,/center,cubic=-0.5)
by0=congrid(mapby.data,scl0*xs,scl0*ys,/center,cubic=-0.5)
dmax=0.88*max(c0)
dmin=0.12*max(c0)
c0[where(c0 gt dmax)]=dmax
c0[where(c0 lt dmin)]=dmin
tvscl,c0
bz1=smooth(congrid(mapbz.data,scl0*xs,scl0*ys,/center,cubic=-0.5),9)
bzmax=max(mapbz.data)
bzmin=min(mapbz.data)
levp=bzmax*[0.3,0.5,0.7,0.9,0.99]
levm=bzmin*[0.99,0.9,0.7,0.5,0.3]
loadct,0
;contour,bz1,/overplot,levels=levp,c_linestyle=0,c_thick=1.0,color=0
contour,bz1,/overplot,levels=[0.0],c_linestyle=4,c_thick=1.0,color=0
;contour,bz1,/overplot,levels=levm,c_linestyle=5,c_thick=1.0,color=0
;fccnve,bx0,by0,/nofil,/nocon,/over,limit=10,scale=0.005,iskip=3,jskip=3,$
;           v_color=255,/Nolabels,/Noaxis

print,'Select a subregion (Region 1) to analyze.'
box_cursor,x1,y1,nx1,ny1
print,'x1,y1,nx1,ny1:',x1,y1,nx1,ny1
plots,[x1,x1+nx1,x1+nx1,x1,x1],[y1,y1,y1+ny1,y1+ny1,y1],/device,thick=2.0,color=0
scl1=5
rec1='Y'
while (rec1 eq 'Y' or rec1 eq 'y') do begin
  wdef,2,scl1*nx1,scl1*ny1
  print,'The present scale factor is ',strtrim(string(scl1),2),'. Do you want to change it (Y/N)?'
  read,rec1,prompt='Input Y/N:  '
  if (rec1 eq 'Y' or rec1 eq 'y') then begin
    read,scl1,prompt='Input scale factor (should be integers):  '
  endif
endwhile

rec='Y'
while (rec eq 'Y' or rec eq 'y') do begin
  wdef,2,scl1*nx1,scl1*ny1
  plot,[0,scl1*nx1],[0,scl1*ny1],xstyle=1,ystyle=1,position=[0,0,scl1*nx1,scl1*ny1]
  c1=congrid(current[x1/scl0:(x1+nx1)/scl0,y1/scl0:(y1+ny1)/scl0],scl1*nx1,scl1*ny1,/center,cubic=-0.5)
  bx1=congrid(mapbx.data[x1/scl0:(x1+nx1)/scl0,y1/scl0:(y1+ny1)/scl0],scl1*nx1,scl1*ny1,/center,cubic=-0.5)
  by1=congrid(mapby.data[x1/scl0:(x1+nx1)/scl0,y1/scl0:(y1+ny1)/scl0],scl1*nx1,scl1*ny1,/center,cubic=-0.5)
  tvscl,c1
  fccnve,bx1,by1,/nofil,/nocon,/over,limit=10,scale=0.04,iskip=19,jskip=19,$
             v_color=0,/Nolabels,/Noaxis
  flag='N'
  print,'Do you want to zoom in (Y/N)?'
  read,flag,prompt='Input Y/N:  '
  if (flag eq 'Y' or flag eq 'y') then begin
    box_cursor,x2,y2,nx2,ny2
    plots,[x2,x2+nx2,x2+nx2,x2,x2],[y2,y2,y2+ny2,y2+ny2,y2],/device,thick=2.0,color=0
    scl2=5
    rec2='Y'
    while (rec2 eq 'Y' or rec2 eq 'y') do begin
      wdef,3,scl2*nx2,scl2*ny2
      print,'The present scale factor is ',strtrim(string(scl2),2),'. Do you want to change it (Y/N)?'
      read,rec2,prompt='Input Y/N:  '
      if (rec2 eq 'Y' or rec2 eq 'y') then begin
        read,scl2,prompt='Input scale factor (should be integers):  '
      endif
    endwhile
    wdef,3,scl2*nx2,scl2*ny2
    plot,[0,scl2*nx2],[0,scl2*ny2],xstyle=1,ystyle=1,position=[0,0,scl2*nx2,scl2*ny2]
    xind0=(x2/scl1+x1)/scl0
    xind1=((x2+nx2)/scl1+x1)/scl0
    yind0=(y2/scl1+y1)/scl0
    yind1=((y2+ny2)/scl1+y1)/scl0
    c2=congrid(current[xind0:xind1,yind0:yind1],scl2*nx2,scl2*ny2,/center,cubic=-0.5)
    bx2=congrid(mapbx.data[xind0:xind1,yind0:yind1],scl2*nx2,scl2*ny2,/center,cubic=-0.5)
    by2=congrid(mapby.data[xind0:xind1,yind0:yind1],scl2*nx2,scl2*ny2,/center,cubic=-0.5)
    tvscl,c2
    fccnve,bx2,by2,/nofil,/nocon,/over,limit=10,scale=0.04,iskip=23,jskip=23,$
               v_color=0,/Nolabels,/Noaxis
    result=defroi(scl2*nx2,scl2*ny2,zoom=scl0*scl1*scl2)

    x20=xind0
    y20=yind0
    nx20=long(scl2*nx2)
    indx=long((result mod nx20)+x20)
    indy=long(result/nx20+y20)
    result=indy*xs+indx
    index=result
    uniqind=index[uniq(index,sort(index))]
    indsize=(size(uniqind))[1]
    help,index

    print,'Click a START point for defining a reference direction in this region.'
    cursor,xsta,ysta,/device,/down
    print,'xsta,ysta:',xsta,ysta
    xsta=double(xsta)
    ysta=double(ysta)
    print,'Click an END point for defining a reference direction in this region.'
    cursor,xend,yend,/device,/down
    print,'xend,yend:',xend,yend
    xend=double(xend)
    yend=double(yend)
    xref=xend-xsta
    yref=yend-ysta
    r = 0.3                          ;len of arrow head
    angle = 22.5 * !dtor            ;Angle of arrowhead
    st = r * sin(angle)             ;sin 22.5 degs * length of head
    ct = r * cos(angle)
    plots,[xsta,xend,xend-(ct*xref+st*yref),xend,xend-(ct*xref-st*yref)],$
          [ysta,yend,yend-(ct*yref-st*xref),yend,yend-(ct*yref+st*xref)],$
          color=0,thick=2.0
    for i=0,indsize-1 do begin
      if ((mapbx.data[uniqind[i]]*xref+mapby.data[uniqind[i]]*yref) ge 0.0) then sign=1.0 else sign=-1.0
      mapbx.data[uniqind[i]]=sign*mapbx.data[uniqind[i]]
      mapby.data[uniqind[i]]=sign*mapby.data[uniqind[i]]
    endfor

    bxo=mapbx.data
    byo=mapby.data
    t1=(1./(2.*lamda))*(shift(Byo,-1,0)-shift(Byo,1,0))
    t2=(1./(2.*lamda))*(shift(Bxo,0,-1)-shift(Bxo,0,1))
    Jzo=(1./mu0)*(t1-t2)
   ; current=-alog(abs(Jzo)/(abs(bzo)+cc))
    current=max(abs(Jzo))-abs(Jzo)

    c2=congrid(current[xind0:xind1,yind0:yind1],scl2*nx2,scl2*ny2,/center,cubic=-0.5)
    bx2=congrid(mapbx.data[xind0:xind1,yind0:yind1],scl2*nx2,scl2*ny2,/center,cubic=-0.5)
    by2=congrid(mapby.data[xind0:xind1,yind0:yind1],scl2*nx2,scl2*ny2,/center,cubic=-0.5)
    wdef,4,scl2*nx2,scl2*ny2
    plot,[0,scl2*nx2],[0,scl2*ny2],xstyle=1,ystyle=1,position=[0,0,scl2*nx2,scl2*ny2]
    tvscl,c2
    fccnve,bx2,by2,/nofil,/nocon,/over,limit=10,scale=0.04,iskip=23,jskip=23,$
               v_color=0,/Nolabels,/Noaxis

  endif else begin
    wdef,1,scl1*nx1,scl1*ny1
    plot,[0,scl1*nx1],[0,scl1*ny1],xstyle=1,ystyle=1,position=[0,0,scl1*nx1,scl1*ny1]
    tvscl,c1
    fccnve,bx1,by1,/nofil,/nocon,/over,limit=10,scale=0.04,iskip=19,jskip=19,$
             v_color=0,/Nolabels,/Noaxis
    result=defroi(scl1*nx1,scl1*ny1,zoom=scl0*scl1)

    x10=x1/scl0
    y10=y1/scl0
    nx10=long(scl1*nx1)
    indx=long((result mod nx10)+x10)
    indy=long(result/nx10+y10)
    result=indy*xs+indx

    index=result
    uniqind=index[uniq(index,sort(index))]
    mapbx.data[uniqind]=-1.0*mapbx.data[uniqind]
    mapby.data[uniqind]=-1.0*mapby.data[uniqind]
    help,index

    bxo=mapbx.data
    byo=mapby.data
    t1=(1./(2.*lamda))*(shift(Byo,-1,0)-shift(Byo,1,0))
    t2=(1./(2.*lamda))*(shift(Bxo,0,-1)-shift(Bxo,0,1))
    Jzo=(1./mu0)*(t1-t2)
   ; current=-alog(abs(Jzo)/(abs(bzo)+cc))
    current=max(abs(Jzo))-abs(Jzo)

    c1=congrid(current[x1/scl0:(x1+nx1)/scl0,y1/scl0:(y1+ny1)/scl0],scl1*nx1,scl1*ny1,/center,cubic=-0.5)
    bx1=congrid(mapbx.data[x1/scl0:(x1+nx1)/scl0,y1/scl0:(y1+ny1)/scl0],scl1*nx1,scl1*ny1,/center,cubic=-0.5)
    by1=congrid(mapby.data[x1/scl0:(x1+nx1)/scl0,y1/scl0:(y1+ny1)/scl0],scl1*nx1,scl1*ny1,/center,cubic=-0.5)
    wdef,3,scl1*nx1,scl1*ny1
    plot,[0,scl1*nx1],[0,scl1*ny1],xstyle=1,ystyle=1,position=[0,0,scl1*nx1,scl1*ny1]
    tvscl,c1
    fccnve,bx1,by1,/nofil,/nocon,/over,limit=10,scale=0.04,iskip=19,jskip=19,$
               v_color=0,/Nolabels,/Noaxis
  endelse
  print,'Do you want to analyze another region (Y/N)?'
  read,rec,prompt='Input Y/N:  '
endwhile

wdef,5,scl0*xs,scl0*ys

bxo=mapbx.data
byo=mapby.data
t1=(1./(2.*lamda))*(shift(Byo,-1,0)-shift(Byo,1,0))
t2=(1./(2.*lamda))*(shift(Bxo,0,-1)-shift(Bxo,0,1))
Jzo=(1./mu0)*(t1-t2)
; current=-alog(abs(Jzo)/(abs(bzo)+cc))
current=max(abs(Jzo))-abs(Jzo)
;====================================
;To establishe Data coordinate system
;because FCCNVE.PRO needs all the plot information
;====================================
plot,[0,scl0*xs],[0,scl0*ys],xstyle=1,ystyle=1,position=[0,0,scl0*xs,scl0*ys]
c0=congrid(current,scl0*xs,scl0*ys,/center,cubic=-0.5)
bx0=congrid(mapbx.data,scl0*xs,scl0*ys,/center,cubic=-0.5)
by0=congrid(mapby.data,scl0*xs,scl0*ys,/center,cubic=-0.5)
dmax=0.88*max(c0)
dmin=0.12*max(c0)
c0[where(c0 gt dmax)]=dmax
c0[where(c0 lt dmin)]=dmin
c00=(c0-min(c0))/(max(c0)-min(c0))*60+195
tv,c00
bz1=smooth(congrid(mapbz.data,scl0*xs,scl0*ys,/center,cubic=-0.5),9)
bzmax=max(mapbz.data)
bzmin=min(mapbz.data)
levp=bzmax*[0.3,0.5,0.7,0.9,0.99]
levm=bzmin*[0.99,0.9,0.7,0.5,0.3]
loadct,0
;contour,bz1,/overplot,levels=levp,c_linestyle=0,c_thick=1.0,color=0
contour,bz1,/overplot,levels=[0.0],c_linestyle=4,c_thick=1.0,color=0
;contour,bz1,/overplot,levels=levm,c_linestyle=5,c_thick=1.0,color=0

save,mapbx,mapby,mapbz,filename='ambiguity_removed_data.sav'

!p.background=0
end
