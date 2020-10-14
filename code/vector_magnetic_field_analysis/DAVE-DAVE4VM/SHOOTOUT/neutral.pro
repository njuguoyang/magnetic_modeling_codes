pro neutral

    DIR='SCHUCK/'
    file=dir+'dave4vm.fts'
;
    SM=7
;
    DATA=mrdfits(file,1)
    truth=mrdfits(dir+'dave4vm.fts',2)
    dave4vm=mrdfits(dir+'dave4vm.fts',3)
;
    i1=100
    i2=200
    NX=i2-i1+1
    NY=i2-i1+1
    x=dindgen(NX)*data.dx*1.d-1
    y=dindgen(NY)*data.dy*1.d-1
    x=x-x[NX/2]
    y=y-y[NY/2]
    xrange=[x[0],x[nx-1]]
    yrange=[y[0],y[ny-1]]
    xtitle='[Mm]'
    ytitle=xtitle
;
    set_plot,'ps'
    !P.font=0
    th=4
    cth=2
    cs=1.2
    extra={ticklen:-0.02,thick:th,xthick:th,ythick:th,charsize:cs,charthick:cth}
    device,file=DIR+'Neutral.eps',/color,bits_per_Pixel=8,xsize=5,ysize=5,/inches
    loadct,0
     mx=max(data.BZ,min=mn)
   print,mn,mx
    plotimage,bytscl(data.BZ[i1:i2,i1:i2],min=-3000,max=3000),/preserve_aspect,imgxrange=xrange,imgyrange=yrange,xtitle=xtitle,ytitle=ytitle,_extra=extra
    mx=max(truth.poynting[i1:i2,i1:i2],min=mn)
    print,mx,mn
    levels=(-1+dindgen(40)*0.5)
    dum=fsc_color(/allcolors,colorstructure=ps)    
    contour,smooth(data.BZ[i1:i2,i1:i2],SM),X,Y,/overplot,levels=[0],/follow,color=ps.green,_extra=extra
;    contour,truth.poynting[i1:i2,i1:i2],/overplot,levels=levels,/follow
    device,file=DIR+'Poynting.eps',/color,bits_per_Pixel=8
    loadct,0
    plotimage,bytscl(data.BZ[i1:i2,i1:i2],min=-3000,max=3000),/preserve_aspect,imgxrange=xrange,imgyrange=yrange,xtitle=xtitle,ytitle=ytitle,_extra=extra
    dum=fsc_color(/allcolors,colorstructure=ps)    
    contour,smooth(data.BZ[i1:i2,i1:i2],SM),X,Y,/overplot,levels=[0],/follow,color=ps.green,_extra=extra
    print,levels
    print
    linestyles=2*(levels lt 0.0)
    contour,smooth(truth.poynting[i1:i2,i1:i2],SM)/1.d5,X,Y,/overplot,levels=levels,/follow,color=ps.red,_extra=extra,c_linestyle=linestyles
    print,linestyles
    device,file=DIR+'Bz.eps',/color,bits_per_Pixel=8
    loadct,0
    plotimage,bytscl(data.BZ[i1:i2,i1:i2],min=-3000,max=3000),/preserve_aspect,imgxrange=xrange,imgyrange=yrange,xtitle=xtitle,ytitle=ytitle,_extra=extra
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    device,file=DIR+'Vector_lines.eps',/color,bits_per_Pixel=8
    loadct,0
    plotimage,bytscl(data.BZ[i1:i2,i1:i2],min=-3000,max=3000),/preserve_aspect,imgxrange=xrange,imgyrange=yrange,xtitle=xtitle,ytitle=ytitle,_extra=extra
    dum=fsc_color(/allcolors,colorstructure=ps)    
    contour,smooth(data.BZ[i1:i2,i1:i2],SM),X,Y,/overplot,levels=[0],/follow,color=ps.green,_extra=extra
    BX0=data.BX[i1:i2,i1:i2]
    BY0=data.BY[i1:i2,i1:i2]
    BX=BX0
    BY=BY0
    BX[*]=!values.d_nan
    BY[*]=!values.d_nan
;
    for i=0,NX-1L,3 do for j=0,NY-1L,3 do begin
        BX[i,j]=BX0[i,j]
        BY[i,j]=BY0[i,j]
    endfor
;
    velovect,BX,BY,X,Y,/overplot,length=4,color=ps.red
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    device,file=DIR+'Vector_lines_BW.eps',color=0,bits_per_Pixel=8
    loadct,0
    plotimage,bytscl(data.BZ[i1:i2,i1:i2],min=-3000,max=3000),/preserve_aspect,imgxrange=xrange,imgyrange=yrange,xtitle=xtitle,ytitle=ytitle,_extra=extra
;    dum=fsc_color(/allcolors,colorstructure=ps)    
    contour,smooth(data.BZ[i1:i2,i1:i2],SM),X,Y,/overplot,levels=[0],/follow,color=0,_extra=extra
    BX0=data.BX[i1:i2,i1:i2]
    BY0=data.BY[i1:i2,i1:i2]
    BX=BX0
    BY=BY0
    BX[*]=!values.d_nan
    BY[*]=!values.d_nan
;
    for i=0,NX-1L,3 do for j=0,NY-1L,3 do begin
        BX[i,j]=BX0[i,j]
        BY[i,j]=BY0[i,j]
    endfor
;
    BX1=BX
    BY1=BY
    BZ=data.BZ[i1:i2,i1:i2]
    POS=where(BZ gt 0,complement=NEG)
    BX1[NEG]=!values.d_nan
    BY1[NEG]=!values.d_nan
    velovect,BX1,BY1,X,Y,/overplot,length=4,color=0
    BX1=BX
    BY1=BY
    BX1[POS]=!values.d_nan
    BY1[POS]=!values.d_nan
    velovect,BX1,BY1,X,Y,/overplot,length=4,color=255
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    device,file=DIR+'Footpoint_lines.eps',/color,bits_per_Pixel=8
    loadct,0
    plotimage,bytscl(data.BZ[i1:i2,i1:i2],min=-3000,max=3000),/preserve_aspect,imgxrange=xrange,imgyrange=yrange,xtitle=xtitle,ytitle=ytitle,_extra=extra
    dum=fsc_color(/allcolors,colorstructure=ps)    
    contour,smooth(data.BZ[i1:i2,i1:i2],SM),X,Y,/overplot,levels=[0],/follow,color=ps.green,_extra=extra
    VX0=truth.ftvx[i1:i2,i1:i2]
    VY0=truth.ftvy[i1:i2,i1:i2]
    VX=VX0
    VY=VY0
    VX[*]=!values.d_nan
    VY[*]=!values.d_nan
;
    for i=0,NX-1L,3 do for j=0,NY-1L,3 do begin
        VX[i,j]=VX0[i,j]
        VY[i,j]=VY0[i,j]
    endfor
    mx_truth=max(vx^2+vy^2,/nan)
;
    velovect,VX,VY,X,Y,/overplot,length=4,color=ps.red
    VX0=dave4vm.ftvx[i1:i2,i1:i2]
    VY0=dave4vm.ftvy[i1:i2,i1:i2]
    VX=VX0
    VY=VY0
    VX[*]=!values.d_nan
    VY[*]=!values.d_nan
;
    for i=0,NX-1L,3 do for j=0,NY-1L,3 do begin
        VX[i,j]=VX0[i,j]
        VY[i,j]=VY0[i,j]
    endfor
    mx_dave4vm=max(vx^2+vy^2,/nan)
;
    velovect,VX,VY,X,Y,/overplot,length=4*mx_dave4vm/mx_truth,color=ps.blue
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    device,file=DIR+'Neutral_lines.eps',/color,bits_per_Pixel=8
    loadct,0
    plotimage,bytscl(data.BZ[i1:i2,i1:i2],min=-3000,max=3000),/preserve_aspect,imgxrange=xrange,imgyrange=yrange,xtitle=xtitle,ytitle=ytitle,_extra=extra
    dum=fsc_color(/allcolors,colorstructure=ps)    
    N=11L
    nevels=(dindgen(N)-N/2)*600
    c_linestyle=2*(nevels lt 0.0)
    contour,smooth(data.BZ[i1:i2,i1:i2],SM),X,Y,/overplot,levels=nevels,/follow,color=ps.green,_extra=extra,c_linestyle=c_linestyle
 ;   contour,smooth(truth.poynting[i1:i2,i1:i2],SM)/1.d5,X,Y,/overplot,levels=levels,/follow,color=ps.red,_extra=extra

    device,/close_file
    stop


end
