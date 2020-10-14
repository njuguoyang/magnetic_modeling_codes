;+
; NAME :
;   footpoints
; PURPOSE :
;   Select footpoints as initial values to integrate the field lines
; CATEGORY :
;
; CALLING SEQUENCE :
;   Called by fieldline3d.pro, which is used for drawing 3 dimensional field lines
; INPUTS :
;
; OUTPUTS :
;
; COMMON BLOCKS :
;
; MODIFICATION HISTORY :
;   2007.09 Guo Yang
;-
pro footpoints_event,event
widget_control,event.top,get_uvalue=infoptr
info1=*infoptr
widget_control,event.id,get_uvalue=widget


case widget of
  'draw':begin
  widget_control,event.id,draw_motion_events=0
  if(event.press eq 1)then begin
    x=event.x
    y=event.y
    xyouts,x,y-2,'+',/device,alignment=0.5,color=0,charsize=1.4,charthick=1.2
    info1.number=info1.number+1L
    label_text=string(x,y,info1.number,$
    format='("X:",i3,1x,"Y:",i3,1x," Number of Footpoints:",i3)')
    widget_control,info1.label,set_value=label_text
    info1.footx[info1.number-1]=x/info1.scale0
    info1.footy[info1.number-1]=y/info1.scale0
    *infoptr=info1
    print,'FP No. and scaled XY Position',$
    (*infoptr).number,(*infoptr).footx[(*infoptr).number-1],(*infoptr).footy[(*infoptr).number-1]
  endif
  if(event.release eq 1)then begin
    widget_control,event.id,draw_motion_events=0
  endif
  end
  else:print,'unrecognized event:',widget
endcase
end

pro footpoints_OK,event
print,'detroying widgets'
widget_control,event.top,/destroy
end

pro control_cleanup,id
print,'cleaning up'
widget_control,id,get_uvalue=infoptr
end

pro footpoints,bz,infoptr,level_plus=level_plus,level_minus=level_minus,levels=levels

device,get_screen_size=ss
xoffset=ss[0]/8
yoffset=ss[1]/8
bzsize=size(bz)
scale0=min([ss[0]/bzsize[1],ss[1]/bzsize[2]])
print,'Initial Scale:',scale0
if (float(ss[0])/bzsize[1] le float(ss[1])/bzsize[2]) then begin
  while (scale0*bzsize[1]+16 ge ss[0]*0.88) do scale0=scale0-1
endif else begin
  while (scale0*bzsize[2]+88 ge ss[1]*0.88) do scale0=scale0-1 
endelse
print,'Final Scale:',scale0
xsize=scale0*bzsize[1]
ysize=scale0*bzsize[2]

tlb=widget_base(column=1,title='contour of Bz',xoffset=xoffset,yoffset=yoffset,tlb_frame_attr=1)
draw=widget_draw(tlb,xsize=xsize,ysize=ysize,uvalue='draw',/button_events)
label=widget_label(tlb,value='Click to select footpoints. Press OK button to end.',$
  /align_center,/dynamic_resize)
buttbase=widget_base(tlb,row=1,/align_center)
label=widget_label(buttbase,value='X:Y: Number of Footpoints:',/align_left,/dynamic_resize)
button=widget_button(buttbase,value='OK',uvalue='OK',event_pro='footpoints_OK')

widget_control,tlb,/realize
widget_control,draw,get_value=winid
wset,winid

device,decomposed=0
loadct,0
contour,bz,levels=level_plus,background=255,color=0,position=[0.,0.,1.0,1.0],min_value=-0.01,xstyle=1+4+16,ystyle=1+4+16
contour,bz,levels=level_minus,max_value=0.01,position=[0.,0.,1.0,1.0],xstyle=1+4+16,ystyle=1+4+16,/noerase,c_linestyle=4,color=0
loadct,1
if (not keyword_set(levels)) then begin
  contour,bz,levels=[0,0.01],position=[0,0,1,1],xstyle=1+4+16,ystyle=1+4+16,/noerase,c_linestyle=1,color=0
endif

;info={winid:0L,label:0L,number:0L,scale0:0.0,footx:fltarr(100),footy:fltarr(100)}
;infoptr=ptr_new(info)
(*infoptr).winid=winid & (*infoptr).label=label & (*infoptr).number=0L
(*infoptr).scale0=scale0
(*infoptr).footx=0.0 & (*infoptr).footy=0.0
widget_control,tlb,set_uvalue=infoptr
loadct,3
xmanager,'footpoints',tlb,cleanup='control_cleanup'
print,'done creating widgets'

end
