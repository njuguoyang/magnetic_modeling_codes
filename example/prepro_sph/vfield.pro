PRO VFIELD,U,V,X,Y,vz=vz, $
       Noerase = noerase, Position = position, $
       Xrange=xrange,Yrange=yrange, $
       Xstyle=xstyle,Ystyle=ystyle, Device = device, $
       Color = color, Title= title, $
       Limit=limit, Scale=scale, $
       Iskip = iskip, Jskip = jskip, $
       Iplot = iplot, Jplot = jplot, $
       Xplot = xplot, Yplot = yplot, $
       Index_color=index_color, Index_size=index_size, $
       Index_charsize=index_charsize, $
       Xindex=xindex, Yindex=yindex, Xarrow=xarrow, $
       F_color=f_color,V_color=v_color, $
       Xtitle=xtitle,Ytitle=ytitle, $
       I0skip=i0skip, J0skip=j0skip,Sample=sample, $
       Xticks=xticks,Yticks=yticks,Xtickv=xtickv,Ytickv=ytickv,$
       Xtickname=xtickname,Ytickname=ytickname,Xminor=xminor,Yminor=yminor, $
       T3d=t3d,thick=thick,Zvalue=zvalue,No_arrow_head=No_arrow_head

;+
; NAME:
;   vfield
; PURPOSE:
;   make a vector plot
; CATEGORY:
;   General graphics.
; CALLING SEQUENCE:
;   vfield, vx, vy, x, y
; INPUTS:
;   vx = 2 dimensional array indicating the x-component of a vector filed
;   vy = 2 dimensional array indicating the y-component of a vector filed
;   x = 1 dimensional array indicating the x-coordinate
;   y = 1 dimensional array indicating the y-coordinate
; OPTIONAL INPUT PARAMETERS:
;
; KEYWORD PARAMETERS:
;   noerase = 
;   position = set the position in normalized cordinate [x0,y0,x1,y1]
;   xrange = desired data range of the x-axis
;   yrange = desired data range of the y-axis
;   xstyle = 
;   ystyle = 
;   device =
;   color = 
;   title = 
;   scale = scale size of the arrows in vector plot
;   limit = limitation for the magnitude of the arrows in vector plot,
;           if vector length is below this value no arrow appears.
;   iskip = arrows are plotted with this skip in x-direction
;   jskip = arrows are plotted with this skip in y-direction
;   i0skip = set skip of starting point of arrow in x-direction
;   j0skip = set skip of starting point of arrow in y-direction
;   iplot = 1-d array to detemine the x-position of arrows to be plotted
;   jplot = 1-d array to detemine the y-position of arrows to be plotted
;    NOTE: iplot,jplot is prior to iskip and jskip.
;   color = color index for vector plot
;   index_color = color index for index vector at the top-right
;   index_size = size of index vector at the top-right
;   xindex = x-position of index in normalized unit
;   yindex = y-position of index in normalized unit
;   xarrow = x-position of index vector in normalized unit
;   xtitle = xtitle
;   ytitle = ytitle
;   v_color = The color index used to draw arrows
;   f_color = The color index used to draw frame
;   xticks = Same as the keyword of 'contour'
;   yticks = Same as the keyword of 'contour'
;   xtickv = Same as the keyword of 'contour'
;   ytickv = Same as the keyword of 'contour'
;   xtickname = Same as the keyword of 'contour'
;   ytickname = Same as the keyword of 'contour'
;   xminor = Same as the keyword of 'contour'
;   yminor = Same as the keyword of 'contour'
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
;    ver 0.0  T. yokoyama (Ntl Astron. Obs. Japan) 1995
;    Yngwie J. Matsumoto 1993.8.*
;    Yngwie J. Matsumoto 1993.10.8
;    matsu@a.phys.nagoya-u.ac.jp Jun 3 1995
;    matsu@a.phys.nagoya-u.ac.jp  Sat Feb 24 22:53:11 JST 1996
;    T. yokoyama (Ntl Astron. Obs. Japan) Oct. 10, 1996
;        stop supporting: keywords; margx,margy
;        new support: keywords; xrange,yrange,iplot,jplot,index_*
;-


;;;;;;;;; Error checking ;;;;;;;;;;;;;;;
        on_error,2                      ;Return to caller if an error occurs
        s = size(u)
        t = size(v)
        if s(0) ne 2 then begin
baduv:   message, 'U and V parameters must be 2D and same size.'
               endif
        if total(abs(s(0:2)-t(0:2))) ne 0 then goto,baduv
;
        if n_params(0) lt 3 then x = findgen(s(1)) else $
                if n_elements(x) ne s(1) then begin
badxy:                  message, 'X and Y arrays have incorrect size.'
                        endif
        if n_params(1) lt 4 then y = findgen(s(2)) else $
;               if n_elements(y) ne s(2) then goto,badxy
                if n_elements(y) ne s(2) then begin
print,n_elements(y),s(2)
goto,badxy
endif
;

;;;;;;;;; default parameters ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if (n_elements(x) eq 0) then x=lindgen(s(1))
if (n_elements(y) eq 0) then y=lindgen(s(2))
if (n_elements(xrange) eq 0) then xrange=[min(x),max(x)]
if (n_elements(yrange) eq 0) then yrange=[min(y),max(y)]

if n_elements(noerase) eq 0 then noerase = 0
if n_elements(xstyle) eq 0 then xstyle=1
if n_elements(ystyle) eq 0 then ystyle=1
if n_elements(device) eq 0 then device = 0
if n_elements(title) le 0 then title = ''

if n_elements(limit) eq 0 then limit=0.
if n_elements(scale) eq 0 then scale=0.5
if n_elements(iskip) eq 0 then iskip=1
if n_elements(jskip) eq 0 then jskip=1
if n_elements(i0skip) eq 0 then i0skip=0
if n_elements(j0skip) eq 0 then j0skip=0
if n_elements(iplot) eq 0 then iplot=i0skip+iskip*lindgen(s(1)/iskip)
if n_elements(jplot) eq 0 then jplot=j0skip+jskip*lindgen(s(2)/jskip)
if (n_elements(xplot) ne 0) then begin
  i0=min(where(xplot gt xrange(0)))
  i1=max(where(xplot lt xrange(1)))
  iplot=lindgen(i1-i0+1)
  for i=0,i1-i0 do begin
    whr=where(x ge xplot(i+i0),count)
    if (count ne 0) then  begin
      iplot(i)=min(whr)
    endif else begin
      goto,out
    endelse
  endfor
out: iplot=iplot(0:i-1)
endif
if (n_elements(yplot) ne 0) then begin
  j0=min(where(yplot gt yrange(0)))
  j1=max(where(yplot lt yrange(1)))
  jplot=lindgen(j1-j0+1)
  for j=0,j1-j0 do begin
    whr=where(y ge yplot(j+j0),count)
    if (count ne 0) then  begin
      jplot(j)=min(whr)
    endif else begin
      goto,out1
    endelse
  endfor
out1: jplot=jplot(0:j-1)
endif
if n_elements(index_size) eq 0 then index_size=5.
if n_elements(index_charsize) eq 0 then index_charsize=!p.charsize
if index_charsize eq 0 then index_charsize=1
if (index_size eq 0.) then index_size=5.
if (n_elements(xindex) eq 0) then xindex=0.85
if (n_elements(yindex) eq 0) then yindex=1.03
if (n_elements(xarrow) eq 0) then xarrow=0.75 

if n_elements(xtitle) le 0 then xtitle=''
if n_elements(ytitle) le 0 then ytitle=''
if n_elements(xticks) eq 0 then xticks=!x.ticks
if n_elements(yticks) eq 0 then yticks=!y.ticks
if n_elements(xtickv) eq 0 then xtickv=!x.tickv
if n_elements(ytickv) eq 0 then ytickv=!y.tickv
if n_elements(xminor) eq 0 then xminor=!x.minor
if n_elements(yminor) eq 0 then yminor=!y.minor
if n_elements(xtickname) eq 0 then xtickname=''
if n_elements(ytickname) eq 0 then ytickname=''
if (n_elements(color) eq 0) then begin
  if (!d.flags and 512) then color=0 else color=!d.table_size-1
endif
if (!d.flags and 512) then c_colors=c_mn
if (n_elements(f_color) eq 0) then f_color = color else f_color = !p.color
if (n_elements(v_color) eq 0) then v_color = color else v_color = !p.color
if (n_elements(index_color) eq 0) then index_color = color
if (n_elements(t3d) eq 0) then t3d=0
if (n_elements(thick) eq 0) then thick=1
if (n_elements(zvalue) eq 0) then zvalue=0

;;;;;;;;; Parameters ;;;;;;;;;;;;;;;

        sz=size(iplot) & iarx=sz(1)
        sz=size(jplot) & jarx=sz(1)

        x0 = xrange(0) & x1=xrange(1)
        y0 = yrange(0) & y1=yrange(1)


;;;;;;;;; Frame ;;;;;;;;;;;;;;;
        if n_elements(position) eq 0 then begin
            plot,[x0,x1],[y1,y0],/nodata,$
            xstyle=xstyle,ystyle=ystyle, $
              title=title,xtitle=xtitle,ytitle=ytitle,$
              noerase=noerase, color=f_color, $
              xticks=xticks,yticks=yticks,xtickv=xtickv,ytickv=ytickv,$
              xtickname=xtickname,ytickname=ytickname,$
              xminor=xminor,yminor=yminor,zvalue=zvalue
        endif else begin
            plot,[x0,x1],[y1,y0],/nodata,$
              xstyle=xstyle,ystyle=ystyle, $
              title=title,xtitle=xtitle,ytitle=ytitle,$
              noerase=noerase, color=f_color, position=position,device=device,$
              xticks=xticks,yticks=yticks,xtickv=xtickv,ytickv=ytickv,$
              xtickname=xtickname,ytickname=ytickname,$
              xminor=xminor,yminor=yminor,zvalue=zvalue
        endelse

        px = !x.window * !d.x_vsize
        py = !y.window * !d.y_vsize
        swx = px(1)-px(0)
        swy = py(1)-py(0)

;;;;;;;;; Arrows ;;;;;;;;;;;;;;;
        r = .3                          ;len of arrow head
        angle = 22.5 * !dtor            ;Angle of arrowhead
        st = r * sin(angle)             ;sin 22.5 degs * length of head
        ct = r * cos(angle)

        for iar=0,iarx-1 do begin
           i=iplot(iar)
           xp=x(i)
        for jar=0,jarx-1 do begin
           j=jplot(jar)
           yp=y(j)
           if (((xp ge x0) and (xp le x1)) and ((yp ge y0) and (yp le y1))) $
            then begin

           vecx=u(i,j)
           vecy=v(i,j)
           vec=sqrt(vecx*vecx+vecy*vecy)

           if vec ge limit then begin
              if (keyword_set(vz)) then begin
                if (vz[i,j] ge 0.0) then v_color=61 else v_color=255
              endif else begin
                v_color=0
              endelse
              vxp=vecx*scale
              vyp=vecy*scale

              xp1=xp+vxp
              yp1=yp+vyp
              if not (keyword_set(no_arrow_head)) then begin
                plots,[xp,xp1,xp1-(ct*vxp+st*vyp),xp1,xp1-(ct*vxp-st*vyp)],$  ;Plot sticks WITH arrow head
                      [yp,yp1,yp1-(ct*vyp-st*vxp),yp1,yp1-(ct*vyp+st*vxp)],$
                      color=v_color,t3d=t3d,thick=thick
              endif else begin
                plots,[xp,xp1],[yp,yp1],color=v_color,t3d=t3d,thick=thick    ;Plot sticks WITHOUT arrow head
              endelse
           endif
           endif
        endfor
        endfor

;;;;;;;;; Index arrow ;;;;;;;;;;;;;;;
if (keyword_set(sample) ne 0) then begin
              loadct,13
              vxp=index_size*scale
              vyp=0.*scale
;             xp=x0 $
;               +((px(0)+xarrow*swx)/!d.x_size-!x.window(0)) $
;               /(!x.window(1)-!x.window(0))*(x1-x0)
;              xp=x1-6*float(!d.x_ch_size)*index_charsize/!d.x_size $
;               *(x1-x0)/(!x.window(1)-!x.window(0))-vxp-71         ;!!!!!!!!!!!!!
;              yp=y0 $
;                +((py(0)+yindex*swy+0.5*!d.y_ch_size)/!d.y_size-!y.window(0)) $
;                /(!y.window(1)-!y.window(0))*(y1-y0)-19           ;!!!!!!!!!!!!!
              x_sample_position = 0.15          ; x- normalized position from the upper right corner
              y_sample_position = 0.05          ; y- normalized position from the upper right corner
              xp=x1 - vxp - x_sample_position*(x1-x0)/(!x.window(1)-!x.window(0))
              yp=y1 - (y_sample_position)*(y1-y0)/(!y.window(1)-!y.window(0)) + 0.5*!d.y_ch_size/!d.y_size*(y1-y0)
              print,'Sample arrow position in data coorinates:',xp,yp
              xp1=xp+vxp
              yp1=yp+vyp
              plots,[xp,xp1,xp1-(ct*vxp-st*vyp),xp1,xp1-(ct*vxp+st*vyp)],$
                    [yp,yp1,yp1-(ct*vyp-st*vxp),yp1,yp1-(ct*vyp+st*vxp)],$
                    color=index_color,thick=thick

;             xp=px(0)+xindex*swx
;             yp=py(0)+yindex*swy
;               xyouts,xp,yp,/device,' = '+ $
;               strmid(strcompress(string(index_size),/remove_all),0,3),$
;                   color=index_color,charsize=index_charsize
              xp=!x.window(1)-x_sample_position                               ;!!!!!!!!!!!!!
;             yp=!y.window(1)+0.2*float(!d.y_ch_size)/!d.y_size
              yp=!y.window(1)-y_sample_position                               ;!!!!!!!!!!!!!
                xyouts,alignment=0,/normal,xp,yp,' = '+ $
                strmid(strcompress(string(index_size),/remove_all),0,4)+' Gauss',$
                    color=index_color,charsize=index_charsize,charthick=3
endif

end
