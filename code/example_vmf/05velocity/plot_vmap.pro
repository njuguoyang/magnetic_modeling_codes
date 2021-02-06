pro plot_vmap,mapbx,mapby,mapbz=mapbz,_extra=_extra

if (n_elements(mapbx) ne 1) or (n_elements(mapby) ne 1) then begin
  print,'Error: Number of elements must be 1'
endif

if (mapbx.xc ne mapby.xc) or $
   (mapbx.yc ne mapby.yc) or $
   (mapbx.dx ne mapby.dx) or $
   (mapbx.dy ne mapby.dy) or $
   (total(size(mapbx.data) - size(mapby.data)) ne 0) then begin
  print,'Data Coordinate is not identical between Bx & By'
endif

sz=size(mapbx.data)
ix=sz(1)
jx=sz(2)

xcdt=mapbx.xc+mapbx.dx*findgen(ix)
xcdt=xcdt-(max(xcdt)-min(xcdt))/2
ycdt=mapbx.yc+mapbx.dy*findgen(jx)
ycdt=ycdt-(max(ycdt)-min(ycdt))/2

  fccnve ,/nofil,/nocon $
    ,xcfc=xcdt,ycfc=ycdt $
    ,mapbx.data,mapby.data,vz=mapbz.data,_extra=_extra


end
