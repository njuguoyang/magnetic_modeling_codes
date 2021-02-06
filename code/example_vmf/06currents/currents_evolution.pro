
fn = ['..\20100917_0000_hmi\02ambiguity\05projection\mapbxyz.sav', $
      '..\20100917_0100_hmi\02ambiguity\05projection\mapbxyz.sav', $
      '..\20100917_0200_hmi\02ambiguity\05projection\mapbxyz.sav', $
      '..\20100917_0300_hmi\02ambiguity\05projection\mapbxyz.sav', $
      '..\20100917_0400_hmi\02ambiguity\05projection\mapbxyz.sav', $
      '..\20100917_0500_hmi\02ambiguity\05projection\mapbxyz.sav', $
      '..\20100917_0600_hmi\02ambiguity\05projection\mapbxyz.sav', $
      '..\20100917_0800_hmi\02ambiguity\05projection\mapbxyz.sav', $
      '..\20100917_0900_hmi\02ambiguity\05projection\mapbxyz.sav', $
      '..\20100917_0936_hmi\02ambiguity\05projection\mapbxyz.sav']
      
count = n_elements(fn)

Jz_tot = fltarr(count)
Jz_tot_abs = fltarr(count)
Jz_tot_pos = fltarr(count)
Jz_tot_neg = fltarr(count)
error_jz = fltarr(count)

xd = 1.5*5
yd = 0.0*5
xrange0 = [-195,-125]
yrange0 = [-425,-370]
wdef,1,1300,800
!p.background=255

for i=0,count-1 do begin
  restore,fn[i],/ver
  if (i le 6) then begin
    xrange = xrange0 + xd*(i-3)
    yrange = yrange0 + yd*(i-3)
  endif else begin
    xrange = xrange0 + xd*(i-2)
    yrange = yrange0 + yd*(i-2)  
    if (i eq 9) then begin        ;for the time 09:36
      xrange = xrange - 0.4*xd
      yrange = yrange - 0.4*yd
    endif
  endelse
  xr1 = xrange
  xr1[0] = xr1[0]+1
  xr1[1] = xr1[1]-1
  yr1 = yrange
  yr1[0] = yr1[0]+1
  yr1[1] = yr1[1]-1
  time1 = utc2int(mapbz.time)
  time1.time = time1.time + 120000      ;120000 ms = 120 s = 2 min
  time2 = int2utc(time1,/vms)
  mapbx.time = time2
  mapby.time = time2
  mapbz.time = time2 
;  if (i eq 5 or i eq 10) then begin
;    dx = -1.0
;    dy = -7.0
;    mapbx.xc = mapbx.xc + dx
;    mapby.xc = mapby.xc + dx
;    mapbz.xc = mapbz.xc + dx    
;    mapbx.yc = mapbx.yc + dy
;    mapby.yc = mapby.yc + dy
;    mapbz.yc = mapbz.yc + dy
;  endif  
  sub_map,mapbx,smapbx,xrange=xrange,yrange=yrange
  sub_map,mapby,smapby,xrange=xrange,yrange=yrange
  sub_map,mapbz,smapbz,xrange=xrange,yrange=yrange
  smapbz.id='SDO/HMI Heliographic Bz'
  plot_map,smapbz,color=0,bcolor=0,charsize=2.0,dmax=800,dmin=-800,xrange=xr1,yrange=yr1
  write_png,'./sbz20100917_'+strtrim(string(i,format='(i4.4)'),2)+'.png',tvrd()
  smapbz.id='SDO/HMI Heliographic Bxyz'
  plot_map,smapbz,color=0,bcolor=0,charsize=2.0,dmax=800,dmin=-800,xrange=xr1,yrange=yr1
  plot_vmap,/over,smapbx,smapby,mapbz=smapbz,limit=200,scale=0.0025,iskip=3,jskip=3,$
            v_color=255,axis_color=0,/Nolabels,v_thick=2.0,/Noaxis
  write_png,'./sbxyz20100917_'+strtrim(string(i,format='(i4.4)'),2)+'.png',tvrd()

  bxo=smapbx.data
  byo=smapby.data
  bzo=smapbz.data
  temp = pb0r(mapbz.time,/arcsec)
  mpa = 695500000.0/temp[2]    ;meters per arcsec
  lamdax = mapbz.dx*mpa        ;lamda is the pixel size in m
  lamday = mapbz.dy*mpa     
  mu0=!dpi*4.d-3  ;The unit of mu0 is G*m/A
  ;cc=0.0001      ;An arbitrary constant to control the contrast of current
  t1=(1./(2.*lamdax))*(shift(Byo,-1,0)-shift(Byo,1,0))
  t2=(1./(2.*lamday))*(shift(Bxo,0,-1)-shift(Bxo,0,1))
  Jzo=(1./mu0)*(t1-t2)
  ss = size(Jzo)
  xs = ss[1]
  ys = ss[2]
  xw1 = 20
  xw2 = 30
  yw1 = 25
  yw2 = 30
  Jzo_sub = Jzo[xw1:xs-1-xw2, yw1:ys-1-yw2]
  noise_lv = 0.02   ;the noise level in unit of A m^(-2)
  index = where(Jzo_sub ge noise_lv or Jzo_sub le -1.0*noise_lv)  ;The index where |Jz| is greater than noise_lv A m^(-2)
  num = n_elements(index)
  Jz_tot_abs[i] = total(abs(Jzo_sub(index)))*lamdax*lamday
 ; index2 = where(Jzo_sub lt noise_lv and Jzo_sub gt -1.0*noise_lv)
 ; noise_lv2 = mean(abs(Jzo_sub[index2]))
 ; print,'noise_lv2:',noise_lv2
  temp_err = fltarr(50)
  for j= 0,49 do begin
    rand_err = randomn(seed,xs,ys)*noise_lv/3.0    ;We assume 3*sigma equals noise_lv, where sigma is the standard deviation of the random noise.
    rand_err_sub = rand_err[xw1:xs-1-xw2, yw1:ys-1-yw2]
    Jzo_sub = Jzo[xw1:xs-1-xw2, yw1:ys-1-yw2]
    Jzo_sub = Jzo_sub + rand_err_sub
    index = where(Jzo_sub ge noise_lv or Jzo_sub le -1.0*noise_lv)
    temp_err[j] = total(abs(Jzo_sub(index)))*lamdax*lamday
  endfor
  error_jz[i] = stddev(temp_err)  
  Jzo_sub = Jzo[xw1:xs-1-xw2, yw1:ys-1-yw2]
  ;current=-alog(abs(Jzo)/(abs(bzo)+cc))
  current=max(abs(Jzo))-abs(Jzo)          ;to reverse the color
  print,'max, min of currents:',max(Jzo),min(Jzo)
  mapcz = smapbz
  mapcz.data = current
  mapcz.id = 'Vertical Current'
  plot_map,mapcz,color=0,bcolor=0,charsize=2.0,dmax=0.08,dmin=0.0,xrange=xr1,yrange=yr1
  write_png,'./scz20100917__'+strtrim(string(i,format='(i4.4)'),2)+'.png',tvrd()
  mapcz.id = 'Vertical Current and Bxy'
  plot_map,mapcz,color=0,bcolor=0,charsize=2.0,dmax=0.08,dmin=0.0,xrange=xr1,yrange=yr1
  plot_vmap,/over,smapbx,smapby,limit=200,scale=0.0025,iskip=3,jskip=3,$
            v_color=255,axis_color=0,/Nolabels,v_thick=2.0,/Noaxis
  write_png,'./scz_bxy20100917__'+strtrim(string(i,format='(i4.4)'),2)+'.png',tvrd()
  
  bzo_sub = Bzo[xw1:xs-1-xw2, yw1:ys-1-yw2]
  index = where(bzo_sub lt 0.0)
  Jzo_sub_pos = Jzo_sub
  Jzo_sub_pos[index] = 0.0     ;The current Jz is set to 0.0 where Bz is negative
  index = where(Jzo_sub_pos ge noise_lv or Jzo_sub_pos le -1.0*noise_lv)    ;index for both postive Jz and negative Jz
  Jz_tot_pos[i] = total(abs(Jzo_sub_pos(index)))*lamdax*lamday
  index = where(bzo_sub gt 0.0)
  Jzo_sub_neg = Jzo_sub
  Jzo_sub_neg[index] = 0.0     ;The current Jz is set to 0.0 where Bz is positive
  index = where(Jzo_sub_neg ge noise_lv or Jzo_sub_neg le -1.0*noise_lv)    ;index for both postive Jz and negative Jz
  Jz_tot_neg[i] = total(abs(Jzo_sub_neg(index)))*lamdax*lamday  
endfor

col = 2
row = 2
!P.MULTI = [0, col, row]

lm=0.08 ;left margin
rm=0.05 ;right margin
bm=0.08 ;bottom margin
tm=0.04 ;top margin
cd=0.10 ;column distance
rd=0.08 ;row distance
cw=(1.0-lm-rm-(col-1)*cd)/col ;column width
rw=(1.0-bm-tm-(row-1)*rd)/row ;row width

wdef,1,1200,1000

max1 = max(smapbz.data)
smapbz.data[xw1,yw1:ys-yw2-1] = max1
smapbz.data[xs-xw2-1,yw1:ys-yw2-1] = max1
smapbz.data[xw1:xs-xw2-1,yw1] = max1
smapbz.data[xw1:xs-xw2-1,ys-yw2-1] = max1

i=0
j=0
x0p=lm+(cw+cd)*j
y0p=bm+(rw+rd)*(row-1-i)
x1p=lm+(cw+cd)*j+cw
y1p=bm+(rw+rd)*(row-1-i)+rw
smapbz.id='SDO/HMI Bz'
plot_map,smapbz,color=0,bcolor=0,charsize=1.5,dmax=800,dmin=-800,xrange=xr1,yrange=yr1,position = [x0p,y0p,x1p,y1p]

i=0
j=1
x0p=lm+(cw+cd)*j
y0p=bm+(rw+rd)*(row-1-i)
x1p=lm+(cw+cd)*j+cw
y1p=bm+(rw+rd)*(row-1-i)+rw
smapbz.id='SDO/HMI Bxyz'
plot_map,smapbz,color=0,bcolor=0,charsize=1.5,dmax=800,dmin=-800,xrange=xr1,yrange=yr1,position = [x0p,y0p,x1p,y1p]
plot_vmap,/over,smapbx,smapby,mapbz=smapbz,limit=200,scale=0.0025,iskip=3,jskip=3,$
          v_color=255,axis_color=0,/Nolabels,v_thick=2.0,/Noaxis

i=1
j=0
x0p=lm+(cw+cd)*j
y0p=bm+(rw+rd)*(row-1-i)
x1p=lm+(cw+cd)*j+cw
y1p=bm+(rw+rd)*(row-1-i)+rw
plot_image,Jzo,color=0,position = [x0p,y0p,x1p,y1p],charsize=2
plots,[xw1,xs-1-xw2,xs-1-xw2,xw1,xw1],[yw1,yw1,ys-1-yw2,ys-1-yw2,yw1],thick=4

i=1
j=1
x0p=lm+(cw+cd)*j
y0p=bm+(rw+rd)*(row-1-i)
x1p=lm+(cw+cd)*j+cw
y1p=bm+(rw+rd)*(row-1-i)+rw
tarray = anytim(['2010-09-17 00:00:00.000', '2010-09-17 01:00:00.000', '2010-09-17 02:00:00.000', $
                 '2010-09-17 03:00:00.000', '2010-09-17 04:00:00.000', '2010-09-17 05:00:00.000', $
                 '2010-09-17 06:00:00.000', '2010-09-17 08:00:00.000', '2010-09-17 09:00:00.000', $
                 '2010-09-17 09:36:00.000'])
utplot,tarray,jz_tot_abs*1.0e-12,'1979-01-01 00:00:00.000',timerange=['2010-09-17 00:00:00.000', '2010-09-17 09:36:00.000'], $
     xstyle=1,ystyle=1+2,color=0,ytitle='Integrated Electric Current (10!e12!n A)',position = [x0p,y0p,x1p,y1p],charsize=1.5,psym=-4,symsize=2
;errplot,tarray,(jz_tot_abs - error_jz)*1.0e-12,(jz_tot_abs + error_jz)*1.0e-12,color=0
ss1 = size(Jzo_sub)
;xs1 = ss1[1]
;ys1 = ss1[2]
;tvscl,congrid(Jzo_sub,8*xs1,8*ys1),100,100
write_png,'./currents_evolution.png',tvrd()



wdef,2,1200,500
col = 2
row = 1
!P.MULTI = [0, col, row]
utplot,tarray,jz_tot_pos*1.0e-12,'1979-01-01 00:00:00.000',timerange=['2010-09-17 00:00:00.000', '2010-09-17 09:36:00.000'], $
     xstyle=1,ystyle=1+2,color=0,ytitle='Integrated Electric Current (10!e12!n A)',psym = -1,charsize=1.5,symsize=2,title='Currents on the Positive Polarity'
utplot,tarray,jz_tot_neg*1.0e-12,'1979-01-01 00:00:00.000',timerange=['2010-09-17 00:00:00.000', '2010-09-17 09:36:00.000'], $
     xstyle=1,ystyle=1+2,color=0,ytitle='Integrated Electric Current (10!e12!n A)',psym = -2,charsize=1.5,symsize=2,title='Currents on the Negative Polarity'
write_png,'./currents_evolution_pos_neg.png',tvrd()

save,tarray,jz_tot_abs,jz_tot_pos,jz_tot_neg,error_jz,filename='currents_evolution.sav'
print,jz_tot_abs
print,error_jz
!p.multi=0
!p.background=0
end
