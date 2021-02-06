pro fccnve, arg0,arg1,arg2,arg3,vz=vz $
     ,Xcord=xcord,Ycord=ycord $
     ,Xcfc=xcfc,Ycfc=ycfc,Xccn=xccn,Yccn=yccn,Xcve=xcve,Ycve=ycve $
     ,Xrange=xrange,Yrange=yrange $
     ,Nofilled=nofilled,Nocontour=nocontour,Novector=novector $
     ,Nolabels=nolabels,Noaxis=noaxis,Position=position $
     ,Levels1=levels1, Clr_index=clr_index, Out=out $
     ,Levels2=levels2 $
     ,C_colors=c_colors,C_linestyle=c_linestyle,C_thick=c_thick $
     ,Scale=scale,Limit=limit,Iskip=iskip,Jskip=jskip $
     ,Iplot=iplot,Jplot=jplot,Xplot=xplot,Yplot=yplot $
     ,V_color=v_color,V_thick=v_thick $
     ,Index_color=index_color,Index_size=index_size $
     ,Index_charsize=index_charsize $
     ,Sample=sample $
     ,Xindex=xindex, Yindex=yindex, Xarrow=xarrow $
     ,Bar=bar, B_position=b_position, B_max=b_max, B_min=b_min $
     ,Axis_color=axis_color, Xtitle=xtitle, Ytitle=ytitle $
     ,Title_charsize=title_charsize $
     ,Useold=useold,T3d=t3d,Zvalue=zvalue,noerase=noerase $
     ,over=over,No_arrow_head=No_arrow_head


;+
; NAME:
;   fccnve
; PURPOSE:
;   Overlay a filled-color map, a contour plot, and a vector plot.
;   This procedure accept non-uniform spacing.
;   (cf. 'tvcn' etc. only accept uniform spacing.)
; CATEGORY:
;   General graphics.
; CALLING SEQUENCE:
;   fccnve, a, b, vx, vy
; INPUTS:
;   a = 2 dimensional array to display (color map).
;   b = 2 dimensional array to display (contour).
;   vx = 2 dimensional array indicating the x-component of a vector filed
;   vy = 2 dimensional array indicating the y-component of a vector filed
; OPTIONAL INPUT PARAMETERS:
;
; KEYWORD PARAMETERS:
;   vz = 2 dimensional array indicating the z-component of a vector filed
;   xcord = 1 dimensional array indicating the x-coordinate
;   ycord = 1 dimensional array indicating the y-coordinate
;   xcfc =
;   ycfc =
;   xccn =
;   yccn =
;   xcve = 
;   ycve = 
;   xrange = 
;   yrange = 
;   nolabels = set to suppress displaying x & y labels
;   position = set the position in normalized cordinate [x0,y0,x1,y1]
;
;   levels1 = 1-d array for contour levels of color map 
;   clr_index = 1-d array for color map. 
;            Each element of this array cooresponds to one in levels1.
;       Note: [levels1(0) is a dummy]
;                        a < levels1(1) --> clr_index(0)
;          levels1(1) <= a < levels1(2) --> clr_index(1)
;          levels1(2) <= a < levels1(3) --> clr_index(2)
;          ....
;          levels1(mx-1) <= a           --> clr_index(mx-1)
;
;   out =  The level assumed outside the plotted region.
;
;   levels2 = 1-d array for levels of contour plot (for "b")
;   c_colors = color index for contour lines (for "b" and "c")
;   c_linestyle = line-style index for contour lines (for "b" and "c")
;   c_thick = thickness index for contour lines (for "b" and "c")
;
;   scale = scale size of the arrows in vector plot
;   limit = limitation for the magnitude of the arrows in vector plot,
;           if vector length is below this value no arrow appears.
;   iskip = arrows are plotted with this skip in x-direction
;   jskip = arrows are plotted with this skip in y-direction
;   iplot = 1-d array to detemine the x-position of arrows to be plotted
;   jplot = 1-d array to detemine the y-position of arrows to be plotted
;   xplot = 1-d array to detemine the x-position of arrows to be plotted
;   yplot = 1-d array to detemine the y-position of arrows to be plotted
;    NOTE: xplot,yplot is prior to iplot and jplot.
;          iplot,jplot is prior to iskip and jskip.
;   v_color = color index for vector plot
;   index_color = color index for index vector at the top-right
;   index_size = size of index vector at the top-right
;
;   bar = set to plot the color bar for color-map
;   b_position = set to give the color bar at the right-hand-side
;                Default is the left-hand-side
;   b_max = the maximum value for the color bar
;   b_min = the minimum value for the color bar
;
;   axis_color = the color index of axis
;   xtitle = title label string for x-axis
;   ytitle = title label string for y-axis
;   title_charsize = size of characters for title label
;   
; OUTPUTS:
;       No explicit outputs.
; COMMON BLOCKS:
;       none.
; SIDE EFFECTS:
;       The currently selected display is affected.
; RESTRICTIONS:
;       None that are obvious.
; PROCEDURE:
;
; MODIFICATION HISTORY:
;       T. Yokoyama (NAOJ) Oct. 8, 1996
;       TY 2002-12-30 fixed bug: calling argument (a,b,vx,vy)->(arg0,1,2,3)
;-

on_error,2                    ;Return to caller if an error occurs

c_mx = !d.table_size-1        ;Brightest color (white)
c_mn = 0b                     ;Darkest   color (black)

;;;;;;;;;;;;;;;; default parameters ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; Determine what to plot

fflg=0 & cflg=0 & vflg=0
if (keyword_set(nofilled) eq 0)  then fflg=1
if (keyword_set(nocontour) eq 0) then cflg=1
if (keyword_set(novector) eq 0)  then vflg=1

if (n_params() eq 0) then goto,error0
if (n_params() eq 1) and not keyword_set(nofilled) then begin
  fflg=1 & cflg=0 & vflg=0
endif
if (n_params() eq 1) and keyword_set(nofilled) then begin
  fflg=0 & cflg=1 & vflg=0
endif

case 1 of
  (fflg and (not cflg) and (not vflg)) : begin ; f
    a=arg0
  end
  ((not fflg) and cflg and (not vflg)) : begin ; c
    case 1 of
      (n_params() eq 1) : begin & b=arg0 & end
      (n_params() ge 4) : begin & b=arg1 & end
      else : goto,error0
    endcase
  end
  ((not fflg) and (not cflg) and vflg) : begin ; v
    case 1 of
      (n_params() eq 2) : begin & vy=arg1 & vx=arg0 & end
      (n_params() ge 4) : begin & vy=arg3 & vx=arg2 & end
      else : goto,error0
    endcase
  end
  (fflg and cflg and (not vflg))       : begin ; fc
    case 1 of
      (n_params() eq 2) : begin & a=arg0 & b=arg1 & end
      (n_params() ge 4) : begin & a=arg0 & b=arg1 & end
      else : goto,error0
    endcase
  end
  (fflg and (not cflg) and vflg)       : begin ; fv
    case 1 of
      (n_params() eq 3) : begin & a=arg0 & vx=arg1 & vy=arg2 & end
      (n_params() ge 4) : begin & a=arg0 & vx=arg2 & vy=arg3 & end
      else : goto,error0
    endcase
  end
  ((not fflg) and cflg and vflg)       : begin ; cv
    case 1 of
      (n_params() eq 3) : begin & b=arg0 & vx=arg1 & vy=arg2 & end
      (n_params() ge 4) : begin & b=arg1 & vx=arg2 & vy=arg3 & end
      else : goto,error0
    endcase
  end
  (fflg and cflg and vflg)             : begin ; fcv
    case 1 of
      (n_params() ge 4) : begin & a=arg0 & b=arg1 & vx=arg2 & vy=arg3 & end
      else : goto,error0
    endcase
  end
  else: goto,error0
endcase

;;; Determine coordinate

L100:
if (n_elements(xcord) ne 0) then begin
  if (n_elements(xcfc) eq 0) then xcfc=xcord
  if (n_elements(xccn) eq 0) then xccn=xcord
  if (n_elements(xcve) eq 0) then xcve=xcord
endif else begin
  if (n_elements(xcve) ne 0) then xcord=xcve
  if (n_elements(xccn) ne 0) then xcord=xccn
  if (n_elements(xcfc) ne 0) then xcord=xcfc
  if (n_elements(xcord) eq 0) then xcord=indgen((size(arg0))(1))
  goto,L100
endelse

L110:
if (n_elements(ycord) ne 0) then begin
  if (n_elements(ycfc) eq 0) then ycfc=ycord
  if (n_elements(yccn) eq 0) then yccn=ycord
  if (n_elements(ycve) eq 0) then ycve=ycord
endif else begin
  if (n_elements(ycve) ne 0) then ycord=ycve
  if (n_elements(yccn) ne 0) then ycord=yccn
  if (n_elements(ycfc) ne 0) then ycord=ycfc
  if (n_elements(ycord) eq 0) then ycord=indgen((size(arg0))(2))
  goto,L110
endelse

if keyword_set(over) then begin
  xrange=!x.crange
  yrange=!y.crange
  goto,L1000
endif

if (n_elements(xrange) eq 0) then begin
  xrange=fltarr(2)
  xrange(0)=min([xcfc(0),xccn(0),xcve(0)])
  xrange(1)=max([xcfc(n_elements(xcfc)-1),xccn(n_elements(xccn)-1) $
                ,xcve(n_elements(xcve)-1)])
endif
if (n_elements(yrange) eq 0) then begin
  yrange=fltarr(2)
  yrange(0)=min([ycfc(0),yccn(0),ycve(0)])
  yrange(1)=max([ycfc(n_elements(ycfc)-1),yccn(n_elements(yccn)-1) $
                ,ycve(n_elements(ycve)-1)])
endif

;if n_elements(position) eq 0 then begin
;  contour,[[0,0],[1,1]],xrange,yrange,/nodata, xstyle=5, ystyle = 5, t3d=t3d $
;    ,noerase=noerase,zvalue=zvalue
;endif else begin
  contour,[[0,0],[1,1]],xrange,yrange,/nodata, xstyle=5, ystyle = 5, t3d=t3d $
    ,noerase=noerase,zvalue=zvalue $
    ,position=position
;endelse

L1000:

px = !x.window * !d.x_vsize	;Get size of window in device units
py = !y.window * !d.y_vsize
dxwindow=(!x.window(1)-!x.window(0))*!d.x_vsize
swx = px(1)-px(0)		;Size in x in device units
swy = py(1)-py(0)		;Size in Y

;;;;;;;;; Filled Contour Plot ;;;;;;;;;;;;;;;
if keyword_set(fflg) then begin
  if n_elements(levels1) eq 0 then begin
    mx=max(a)
    mn=min(a)
    levels1=mn+(mx-mn)/16.*findgen(16)
  endif
  mlvx=n_elements(levels1)
  if n_elements(clr_index) eq 0 then begin
    ; clr_index=(!d.n_colors<256)/16*lindgen(16)
      clr_index=(!d.table_size-1)/mlvx*lindgen(mlvx)
  endif
  if n_elements(out) eq 0 then out=min(a)
  if n_elements(bar) eq 0 then bar=0
  if n_elements(b_position) eq 0 then b_position=0
  if n_elements(b_min) eq 0 then b_min=min(levels1)
  if n_elements(b_max) eq 0 then b_max=max(levels1)

  sza = size(a)             ;Size of image
  szlv = size(levels1)
  if (szlv(1) ge 30) then begin
    print,'# of Levels for filled contour is less eq than 29'
    return
  endif
  cn_itx=szlv(1)/29+1
  amari=szlv(1)-29*(cn_itx-1)
  szcl = size(clr_index)
  cl_itx=szlv(1)/szcl(1)+1
  cl_am=szlv(1)-szcl(1)*(cl_itx-1)
; clr_index0=bytarr(szlv(1)+1)
  clr_index0=lonarr(szlv(1)+1)
  if (cl_itx ge 2) then $
    for cl_it=0,cl_itx-2 do $
    clr_index0(cl_it*szcl(1):(cl_it+1)*szcl(1)-1)=clr_index
  cl_it=cl_itx-1
  clr_index0(cl_it*szcl(1):cl_it*szcl(1)+cl_am)=clr_index(0:cl_am)


  if ((long(strmid(!version.release,0,1)) le 3) or keyword_set(useold)) $
    then begin

    i0=(min(where(xcfc ge xrange(0)))-1)>0
    i1=(max(where(xcfc le xrange(1)))+1)<(sza(1)-1)
    j0=(min(where(ycfc ge yrange(0)))-1)>0
    j1=(max(where(ycfc le yrange(1)))+1)<(sza(2)-1)

    aa=a(i0:i1,j0:j1)
    szaa=size(aa)
    xcfca=xcfc(i0:i1) & ycfca=ycfc(j0:j1)

    data2=replicate(out,szaa(1)+2,szaa(2)+2)
    xx2=[xcfca(0)-(xcfca(1)-xcfca(0)),xcfca $
        ,xcfca(szaa(1)-1)+(xcfca(szaa(1)-1)-xcfca(szaa(1)-2))]
    yy2=[ycfca(0)-(ycfca(1)-ycfca(0)),ycfca $
        ,ycfca(szaa(2)-1)+(ycfca(szaa(2)-1)-ycfca(szaa(2)-2))]
    data2(1,1)=aa

    spawn,'echo $$',pid
    tmpfile=pid(0)+'tmp_cpath.dat'
    contour,/noerase,data2,xx2,yy2,levels=levels1,path_filename=tmpfile,$
         xst=5,yst=5,$
         pos = [px(0),py(0), px(0)+swx,py(0)+swy],/dev ,$
         xrange=xrange,yrange=yrange,t3d=t3d,zvalue=zvalue
    base_color=clr_index(0)
    polyfill,[px(0),px(0)+swx,px(0)+swx,px(0)] $
       ,[py(0),py(0),py(0)+swy,py(0)+swy],/dev $
       ,color=base_color
    polycontour,tmpfile,color_index=[clr_index0(0),clr_index0],/delete_file

   ;    make a sharp frame
    frame_color=c_mn
    if (!d.flags and 512) ne 0 then frame_color=c_mn
    dx0=1.4*(xx2(1)-xx2(0)) & dx1=1.4*(xx2(szaa(1)+1)-xx2(szaa(1)))
    dy0=1.7*(yy2(1)-yy2(0)) & dy1=1.7*(yy2(szaa(2)+1)-yy2(szaa(2)))
    polyfill,[xrange(0),xrange(0)-dx0,xrange(0)-dx0,xrange(0)],$
       [yrange(0),yrange(0),yrange(1),yrange(1)],$
       color=frame_color
    polyfill,[xrange(1),xrange(1)+dx1,xrange(1)+dx1,xrange(1)],$
       [yrange(0),yrange(0),yrange(1),yrange(1)],$
       color=frame_color
    polyfill,[xrange(0)-dx0,xrange(1)+dx1,xrange(1)+dx1,xrange(0)-dx0],$
       [yrange(0)-dy0,yrange(0)-dy0,yrange(0),yrange(0)],$
       color=frame_color
    polyfill,[xrange(0)-dx0,xrange(1)+dx1,xrange(1)+dx1,xrange(0)-dx0],$
       [yrange(1)+dy1,yrange(1)+dy1,yrange(1),yrange(1)],$
       color=frame_color

  endif else begin ;version 4 or later
      data2=replicate(out,sza(1)+2,sza(2)+2)
      xx2=[xcfc(0)-(xcfc(1)-xcfc(0)),xcfc $
          ,xcfc(sza(1)-1)+(xcfc(sza(1)-1)-xcfc(sza(1)-2))]
      yy2=[ycfc(0)-(ycfc(1)-ycfc(0)),ycfc $
          ,ycfc(sza(2)-1)+(ycfc(sza(2)-1)-ycfc(sza(2)-2))]
      data2(1,1)=a

      if (cn_itx ge 2) then begin
        cn_it=0
        levelstmp=levels1(cn_it*29+0:cn_it*29+28)
        clr_indextmp=clr_index0(cn_it*29+0:cn_it*29+29)
      endif else begin  ; cn_itx lt 2
        levelstmp=levels1
        clr_indextmp=clr_index0
      endelse  ; cn_itx lt 2

;     if ((!d.name ne 'Z') or (not keyword_set(t3d))) then begin
        contour,/noerase,/cell_fill,a,xcfc,ycfc,levels=levelstmp, $
          xst=5,yst=5,$
          pos = [px(0),py(0), px(0)+swx,py(0)+swy],/dev,$
          xrange=xrange,yrange=yrange, $
          c_colors=clr_indextmp, t3d=t3d, zvalue=zvalue
;     endif else begin  ; t3d
;       screen=!d.name
;       set_plot,'x',/copy
;       window,/free,/pixmap,xsize=480,ysize=480
;;       device,pseudo_color=8
;        contour,/fill,data2,xx2,yy2,levels=levelstmp, $
;          xst=5,yst=5,$
;          pos = [0,0,1,1], $
;          xrange=xrange,yrange=yrange, $
;          c_colors=clr_indextmp
;        img=tvrd()
;        wdelete
;    
;        set_plot,screen
;        x0=dev2nml(px(0)) & x1=dev2nml(px(0)+swx)
;        y0=dev2nml(py(0),/ycord) & y1=dev2nml(py(0)+swy,/ycord)
;        polyfill,/normal $
;          ,[[x0,y0,zvalue],[x1,y0,zvalue],[x1,y1,zvalue],[x0,y1,zvalue]] $
;          ,image_coord=[[0,0],[480,0],[480,480],[0,480]],pattern=img,/t3d $
;          ,color=clr_indextmp
;        contour,/noerase,/nodata,data2,xx2,yy2,levels=levelstmp, $
;          xst=5,yst=5,$
;          pos = [px(0),py(0), px(0)+swx,py(0)+swy],/dev,$
;          xrange=xrange,yrange=yrange, $
;          c_colors=clr_indextmp, t3d=t3d, zvalue=zvalue
;      endelse

  endelse ; version 4 or later


;;;;; color bar 
if (bar eq 1) then  begin

  pxb=px
  if (b_position eq 0) then begin
  pxb(0)=px(0)
  pxb(1)=px(0)+0.05*dxwindow
  px(0) =px(0)+0.15*dxwindow
  endif else begin
  pxb(0)=px(0)+0.95*dxwindow
  pxb(1)=px(1)
  px(1) =px(0)+0.85*dxwindow
  endelse
  pyb=py-8.
  swxb = pxb(1)-pxb(0)		;Size in x in device units
  swyb = pyb(1)-pyb(0)		;Size in Y

ib=20 & jb=256
db=float(b_max-b_min)
yb=b_min+db/jb*findgen(jb)
xb=findgen(ib)
bardata=fltarr(ib,jb) & for i=0,ib-1 do bardata(i,*)=yb
bardata2=replicate(out,ib+2,jb+2)
xb2=[xb(0)-(xb(1)-xb(0)),xb,xb(ib-1)+(xb(ib-1)-xb(ib-2))]
yb2=[yb(0)-(yb(1)-yb(0)),yb,yb(jb-1)+(yb(jb-1)-yb(jb-2))]
bardata2(1,1)=bardata
if (long(strmid(!version.release,0,1)) le 3) then begin
spawn,'echo $$',pid
tmpfile=pid(0)+'tmp_cpath.dat'
contour,bardata2,xb2,yb2,levels=levels1,path_filename=tmpfile,$
         xst=5,yst=5,$
         pos = [pxb(0),pyb(0), pxb(0)+swxb,pyb(0)+swyb],/dev,$
         xrange=[xb(0),xb(ib-1)],$
         yrange=[yb(0),yb(jb-1)],/noerase
base_color=clr_index(0)
polyfill,[xb(0),xb(ib-1),xb(ib-1),xb(0)],$
   [yb(0),yb(0),yb(jb-1),yb(jb-1)],$
   color=base_color
polycontour,tmpfile,color_index=clr_index
spawn,'/bin/rm '+tmpfile
endif else begin
contour,/fill,bardata2,xb2,yb2,levels=levels1, $
         xst=5,yst=5,$
         pos = [pxb(0),pyb(0), pxb(0)+swxb,pyb(0)+swyb],/dev,$
         xrange=[xb(0),xb(ib-1)],$
         yrange=[yb(0),yb(jb-1)],/noerase, $
         c_colors=clr_index
endelse

frame_color=c_mn
if (!d.flags and 512) ne 0 then frame_color=c_mn
polyfill,[xb(0),xb2(0),xb2(0),xb(0)],$
   [yb(0),yb(0),yb(jb-1),yb(jb-1)],$
   color=frame_color
polyfill,[xb(ib-1),xb2(ib+1),xb2(ib+1),xb(ib-1)],$
   [yb(0),yb(0),yb(jb-1),yb(jb-1)],$
   color=frame_color
polyfill,[xb2(0),xb2(ib+1),xb2(ib+1),xb2(0)],$
   [yb2(0),yb2(0),yb(0),yb(0)],$
   color=frame_color
polyfill,[xb2(0),xb2(ib+1),xb2(ib+1),xb2(0)],$
   [yb2(jb+1),yb2(jb+1),yb(jb-1),yb(jb-1)],$
   color=frame_color

   blank_arr=[' ',' ',' ',' ',' ',' ',' ',' ']
   if n_elements(axis_color) eq 0 then begin
     if (!d.flags and 512) ne 0 then axis_color=c_mn else axis_color=c_mx
   endif
   axis,yaxis=0                 ,color=axis_color,yst=1
   axis,yaxis=1,ytickn=blank_arr,color=axis_color,yst=1
   plots,[xb(0),xb(ib-1)],[yb(0),yb(0)],color=axis_color
   plots,[xb(0),xb(ib-1)],[yb(jb-1),yb(jb-1)],color=axis_color

endif

endif ; keyword(filled)

;;;;;;;;;;;;;;;; Vector plot ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if keyword_set(vflg) then begin
  if n_elements(v_color) eq 0 then begin
    v_color=c_mx
    if (!d.flags and 512) ne 0 then v_color=c_mn
  endif
  if n_elements(v_thick) eq 0 then v_thick=1
  if (n_elements(index_color) eq 0) then begin
    index_color=c_mx
    if (!d.flags and 512) ne 0 then index_color=c_mn
  endif
  if n_elements(index_size) eq 0 then index_size=0.
  if n_elements(index_charsize) eq 0 then index_charsize=!p.charsize
  if (n_elements(xarrow) eq 0) then xarrow=0.75
  if (n_elements(xindex) eq 0) then xindex=0.85
  if (n_elements(yindex) eq 0) then yindex=1.01
  if (keyword_set(No_arrow_head)) then begin  
    vfield,vx,vy,xcve,ycve,vz=vz,/noerase,xst=5,yst=5,$
        pos = [px(0),py(0), px(0)+swx,py(0)+swy],/dev,$
        scale=scale,limit=limit,iskip=iskip,jskip=jskip, $
        iplot=iplot,jplot=jplot,xplot=xplot,yplot=yplot, $
        color=v_color,index_color=index_color,index_size=index_size, $
        xindex=xindex,yindex=yindex,xarrow=xarrow,sample=sample, $
        xrange=xrange,yrange=yrange,index_charsize=index_charsize, $
        t3d=t3d,thick=v_thick,zvalue=zvalue,/No_arrow_head
  endif else begin
    vfield,vx,vy,xcve,ycve,vz=vz,/noerase,xst=5,yst=5,$
        pos = [px(0),py(0), px(0)+swx,py(0)+swy],/dev,$
        scale=scale,limit=limit,iskip=iskip,jskip=jskip, $
        iplot=iplot,jplot=jplot,xplot=xplot,yplot=yplot, $
        color=v_color,index_color=index_color,index_size=index_size, $
        xindex=xindex,yindex=yindex,xarrow=xarrow,sample=sample, $
        xrange=xrange,yrange=yrange,index_charsize=index_charsize, $
        t3d=t3d,thick=v_thick,zvalue=zvalue
  endelse
endif
;;;;;;;;; Contour Plot ;;;;;;;;;;;;;;;
if keyword_set(cflg) then begin
  if n_elements(levels2) eq 0 then begin
    mx=max(b)
    mn=min(b)
    levels2=mn+(mx-mn)/16.*findgen(16)
  endif
  if n_elements(c_colors) eq 0 then begin
    c_colors = c_mx
    if (!d.flags and 512) ne 0 then c_colors=c_mn
  endif
  if n_elements(c_thick) eq 0 then c_thick = 1.0
  if n_elements(c_linestyle) eq 0 then c_linestyle = 0

szlv = size(levels2)
cn_itx=szlv(1)/30+1
amari=szlv(1)-30*(cn_itx-1)

if (cn_itx ge 2) then begin
  for cn_it=0,cn_itx-2 do begin
    levelstmp=levels2(cn_it*30+0:cn_it*30+29)
    contour,b,xccn,yccn,/noerase,xstyle=5,ystyle=5,$	
        pos = [px(0),py(0), px(0)+swx,py(0)+swy],/dev,$
	   levels=levelstmp,$
        c_colors=c_colors,c_linestyle=c_linestyle,c_thick=c_thick, $
        xrange=xrange,yrange=yrange,t3d=t3d, zvalue=zvalue
  endfor
  if (amari ge 1) then begin
    cn_it=cn_itx-1
    levelstmp=levels2(cn_it*30+0:cn_it*30+amari-1)
    contour,b,xccn,yccn,/noerase,xstyle=5,ystyle=5,$	
        pos = [px(0),py(0), px(0)+swx,py(0)+swy],/dev,$
	   levels=levelstmp,$
        c_colors=c_colors,c_linestyle=c_linestyle,c_thick=c_thick, $
        xrange=xrange,yrange=yrange,t3d=t3d, zvalue=zvalue
  endif
endif else begin
    contour,b,xccn,yccn,/noerase,xstyle=5,ystyle=5,$	
        pos = [px(0),py(0), px(0)+swx,py(0)+swy],/dev,$
	   levels=levels2,$
        c_colors=c_colors,c_linestyle=c_linestyle,c_thick=c_thick, $
        xrange=xrange,yrange=yrange,t3d=t3d, zvalue=zvalue
endelse

endif
;;;;;;;;; Axises ;;;;;;;;;;;;;;;
; last for overwrite on images

if n_elements(axis_color) eq 0 then begin
  if (!d.flags and 512) ne 0 then axis_color=c_mn else axis_color=c_mx
endif
if n_elements(xtitle) eq 0 then xtitle=' '
if n_elements(ytitle) eq 0 then ytitle=' '
;if n_elements(title_charsize) eq 0 then title_charsize=1.

if keyword_set(noaxis) then return
blank_arr=[' ',' ',' ',' ',' ',' ',' ',' ']
if keyword_set(nolabels) then begin
   axis,xaxis=0,xtickn=blank_arr,color=axis_color,xst=1,t3d=t3d, zvalue=zvalue
   axis,xaxis=1,xtickn=blank_arr,color=axis_color,xst=1,t3d=t3d, zvalue=zvalue
   axis,yaxis=0,ytickn=blank_arr,color=axis_color,yst=1,t3d=t3d, zvalue=zvalue
   axis,yaxis=1,ytickn=blank_arr,color=axis_color,yst=1,t3d=t3d, zvalue=zvalue
endif else begin
   axis,xaxis=1,xtickn=blank_arr,color=axis_color,xst=1,t3d=t3d, zvalue=zvalue
   axis,yaxis=1,ytickn=blank_arr,color=axis_color,yst=1,t3d=t3d, zvalue=zvalue
 if (n_elements(title_charsize)) then begin
   axis,xaxis=0                 ,color=axis_color,xst=1 $
     ,xtitle=xtitle,charsize=title_charsize,t3d=t3d, zvalue=zvalue
   axis,yaxis=0                 ,color=axis_color,yst=1 $
     ,ytitle=ytitle,charsize=title_charsize,t3d=t3d, zvalue=zvalue
 endif else begin
   axis,xaxis=0                 ,color=axis_color,xst=1 $
     ,xtitle=xtitle,t3d=t3d, zvalue=zvalue
   axis,yaxis=0                 ,color=axis_color,yst=1 $
     ,ytitle=ytitle,t3d=t3d, zvalue=zvalue
 endelse
endelse

return

error0 :
   print,' number of data is inconsistent.'

end
