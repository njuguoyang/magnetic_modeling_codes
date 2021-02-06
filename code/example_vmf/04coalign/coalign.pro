;+
; NAME :
;   coalign
; PURPOSE :
;   Align map1 referred to map2 by a feature identification method.
; CATEGORY :
; CALLING SEQUENCE :
;   coalign,map1,map2
; INPUTS :
;   map1
;   map2
;   Map1 and map2 are two maps. For example, they should contain the 
;   contents as follows:
;** Structure <ae4578>, 13 tags, length=1580480, data length=1580474, refs=1:
;   DATA            DOUBLE    Array[483, 409]
;   XC              FLOAT           514.694
;   YC              FLOAT           135.609
;   DX              FLOAT          0.233491
;   DY              FLOAT          0.233491
;   TIME            STRING    '23-May-2008 17:12:29.450'
;   DUR             FLOAT           0.00000
;   ID              STRING    'THEMIS/MTR Fe 6302 Bz'
;   XUNITS          STRING    'arcsec'
;   YUNITS          STRING    'arcsec'
;   ROLL_ANGLE      FLOAT           0.00000
;   ROLL_CENTER     FLOAT     Array[2]
;   SOHO            INT              0
; OPTIONAL INPUT PARAMETERS :
; KEYWORD PARAMETERS :
; OUTPUTS :
;   map1 that aligned with map2
; MODIFICATION HISTORY :
;   2009.02 Guo Yang@LESIA, Observatoire de Paris
;-

pro coalign,map1,map2

tagname=tag_names(map1)
temp=where(tagname eq 'RTIME',ico)
if ico gt 0. then begin
  time_mtr=map1.rtime
endif else begin
  time_mtr=map1.time
endelse
dxmtr=map1.dx & dymtr=map1.dy
window,/free
plot_map,map1
window,/free
plot_map,map2,/log_scale
sub_map,map2,smap1,/noplot
sub_map,smap1,smap
plot_map,smap
rmap=drot_map(smap,time=time_mtr)
dxmdi=rmap.dx
dymdi=rmap.dy
setpts,p,map1.data,rmap.data
xmeanmdi=mean(p[0,0,*])
ymeanmdi=mean(p[1,0,*])
xmeanmtr=mean(p[0,1,*])
ymeanmtr=mean(p[1,1,*])
ssmdi=size(rmap.data)
xcmdi=ssmdi[1]/2.0
ycmdi=ssmdi[2]/2.0
xcmtr_ref2mdi=xmeanmtr-dxmdi/dxmtr*(xmeanmdi-xcmdi)
ycmtr_ref2mdi=ymeanmtr-dymdi/dymtr*(ymeanmdi-ycmdi)
ssmtr=size(map1.data)
xcmtr_index=ssmtr[1]/2.0
ycmtr_index=ssmtr[2]/2.0
xcmtr=(xcmtr_index-xcmtr_ref2mdi)*dxmtr+rmap.xc
ycmtr=(ycmtr_index-ycmtr_ref2mdi)*dymtr+rmap.yc
map1.xc=xcmtr
map1.yc=ycmtr
window,/free
plot_map,map1
dmap2=drot_map(map2,time=time_mtr)
plot_map,dmap2,/over

end
