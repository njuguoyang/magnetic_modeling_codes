;+
; NAME :
;   fieldline3d
; PURPOSE :
;   Calculate and draw the field lines for a given 3-D vector field.
; CATEGORY :
;
; CALLING SEQUENCE :
;   fieldline3d,bx,by,bz,dx=dx,dy=dy,line_color=line_color,/nocontour,isn=isn,xindex,yindex
;
;   It calls the following procedures: footpoints.pro, corner.pro, xitp.pro, differential3.pro
;   Click to choose the footpoints for field line integration. The maximum
;   number of footpoints is 3000.
; INPUTS :
;   Bx[nx,ny,nz]
;   By[nx,ny,nz]
;   Bz[nx,ny,nz]
;   The dimensions should be equal for Bx, By, and Bz.
;
; OPTIONAL INPUT PARAMETERS (KEYWORD PARAMETERS):
;   dx,dy,dz  :  The grid space.
;   xyz0      :  Coordinates of the origin point (lower left and bottom corner). 
;                It must be an array with 3 elements. For example, 
;                xyz0=[100.0,200.0,0.0] means the origin is at [100.0,200.0,0.0].
;   line_color:  Set this keyword to choose the field line color. The color
;                is given by RGB 3 components. For example, line_color=[255,0,0]
;                sets lines red, line_color=[0,255,0] sets lines green, and so on.
;   line_thick:  Set this keyword to choose the field line thickness. The default
;                value is 1.6.
;   isn       :  Integration Step Number. The default is 3000, which means the
;                maximum integration step number is 3000 for each line. If this step
;                reaches, outputs give a warning "The maximum ISN (Integration Step
;                Number) reaches! Increase this number."
;   step_size :  Stepsize in the fraction of dz. It is 1/16 of dz by default.
;   xindex/yindex/zindex : X/Y/Z coordinates of the starting points of the field 
;                          lines to be plotted. If none is specified, the foot-points
;                          are selected by clicking the bottom boundary. If xindex 
;                          and yindex are specified and zindex in not, zindex is 0 
;                          for all the foot-points.                          
;   /sav_start:  Set this keyword to save the positions of the starting foot-points,
;                i.e. xindex/yindex/zindex in the file 'footpoints_position.sav'.
;   /sav_fl   :  Set this keyword to save the positions of the field lines in the file
;                'curve_infoxxxx.sav', where 'xxxx' is the field line No. At present,
;                this keyword only takes effect in cobination with the keyword 'both'.
;   /sav_pend :  Set this keyword to save the position of the positive ends of field
;                lines in the file 'positive_ends.sav', where positive is defined as 
;                the normal component of the magnetic field pointing inward of a volume. 
;                For example, on the bottom surface of a rectangular box, the positive 
;                end is where Bz > 0, but on the top surface it is where Bz < 0. The 
;                positive can be defined similarly on the lateral surfaces.
;   /sav_nend :  Set this keyword to save the position of the negative ends of field
;                lines in the file 'negative_ends.sav'.
;   linexyz   :  Set this keyword to save the position of the last magnetic field line       
;   /forward  :  Integration along the field lines forward.
;   /backwad  :  Integration along the field lines backward.
;   /both     :  Integration along the field lines in both directions. The priority of
;                these integration direction controlling keywords is /both > /forward 
;                > /backward. If none is specified, the integration direction is determined
;                by the Bz at the starting point. If Bz(0)<0, the direction is backward;
;                If Bz(0)>=0, the direction is upward.
;   levels    :  Specifies a vector containing the contour levels drawn by the CONTOUR 
;                procedure. A contour is drawn at each level in LEVELS. 
;   /percent  :  With this keyword, the positive and negative levels are regarded as 
;                fraction of maximum and minimum of the magnetic fields. 
;   /nocontour:  Set this keyword to hide the contour levels. It is useful
;                when you repeat drawing another group of field lines. Without
;                this keyword, the contour levels at the bottom will be drawn
;                again.
;   /noimage  :  Set this keyword to hide the background image.
;   /noplot   :  Set this keyword to prevent invoking iplot. At present, this keyword only 
;                takes effect in cobination with the keyword 'both'.
;
; OUTPUTS :
;   3-D magnetic field lines
; COMMON BLOCKS :
; SIDE EFFECTS :
; RESTRICTIONS :
; PROCEDURE :
; MODIFICATION HISTORY :
;   originally coded by M.T.Song by FORTRAN
;   2007.02 Guo Yang@Nanjing University, translated these codes from FORTRAN to IDL
;   2007.09 Guo Yang@Nanjing University, GUI and directly 3-D visualization without projection
;   2008.10 Guo Yang@Nanjing University & Observatoire de Paris, adopted the fourth-order Runge-Kutta method
;   2008.11 Li Hui  @PMO, added OPTIONAL INPUT PARAMETERS xindex, yindex
;   2009.11 Guo Yang@Nanjing University, added OPTIONAL INPUT PARAMETERS zindex, enabled to integrate the field
;           lines forward or backward along the field line, or along both directions.
;
;   This software is provided without any warranty. Permission to use,
;   copy, modify. Distributing modified or unmodified copies is granted,
;   provided this disclaimer and information are included unchanged.
;-

pro fieldline3d,Bx,By,Bz,dx=dx,dy=dy,dz=dz,xyz0=xyz0, $
                line_color=line_color,line_thick=line_thick, $
                isn=isn,step_size=step_size, $
                xindex=xindex,yindex=yindex,zindex=zindex, $
                sav_start=sav_start,sav_fl=sav_fl,sav_pend=sav_pend,sav_nend=sav_nend, $
                linexyz=linexyz,silent=silent, $
                forward=forward,backward=backward,both=both, $
                levels=levels,percent=percent, $
                nocontour=nocontour,noimage=noimage,noplot=noplot

common bvertex,xv,Bx0,By0,Bz0,nx,ny,nz
common prm1,orig,DELX1,DELY1,DELZ1,sig
common prm2,Binterp

If(n_params() Lt 3) Then Begin
  print,'ERROR!'
  print,'CALLING SEQUENCE :'
  print,'======================='
  print,'fieldline3d, Bx, By, Bz'
  print,'======================='
  Return
Endif
if (not keyword_set(silent)) then print,'+++ BEGIN +++'

Bx0=double(Bx) & By0=double(By) & Bz0=double(Bz)
;Bx=0.0 & By=0.0 & Bz=0.0
ss=size(Bx0)
nx=ss[1] & ny=ss[2] & nz=ss[3]
xv=dblarr(3,8) & DELX1=0.0 & DELY1=0.0 & DELZ1=0.0

if(keyword_set(isn)) then isn=isn else isn=3000    ;Integration Step Number
if(keyword_set(line_color)) then line_color=line_color else line_color=[255,0,0]
if(keyword_set(line_thick)) then line_thick=line_thick else line_thick=1.6
if(keyword_set(dx)) then DELX1=double(dx) else DELX1=1.0d
if(keyword_set(dy)) then DELY1=double(dy) else DELY1=1.0d
if(keyword_set(dz)) then DELZ1=double(dz) else DELZ1=1.0d
if(n_elements(xyz0) eq 3) then xyz0=double(xyz0) else xyz0=[0.d,0.d,0.d] 
if(keyword_set(step_size)) then h=double(step_size) else h=DELZ1/16.0    ;Integration Step Size
orig=xyz0
s=0.0    ;independent variable of the integration equation of field lines

;contour levels
BZMIN=MIN(Bz0[*,*,0])
BZMAX=MAX(Bz0[*,*,0])
if (keyword_set(levels) and not keyword_set(percent)) then begin
  level_minus=levels[where(levels le 0.0)]
  level_plus=levels[where(levels gt 0.0)]
endif 
if (keyword_set(levels) and keyword_set(percent)) then begin
  level_minus=-1.0*levels[where(levels le 0.0)]*BZMIN
  level_plus=levels[where(levels gt 0.0)]*BZMAX
endif
if (not keyword_set(levels)) then begin
CLP1=dblarr(33)  
CL2=BZMAX
CL11=-BZMIN
DEL=CL2/14.
DEL1=CL11/15.
CLP1(0)=-CL11+DEL1/2.
FOR I=1,14 DO BEGIN
  CLP1(I)=CLP1(0)+DEL1*I
ENDFOR
CLP1(15)=CLP1(14)/2.
CLP1(16)=CLP1(15)/2.
CLP1(17)=0.
CLP1(18)=-CLP1(16)
FOR I=19,32 DO BEGIN
  CLP1(I)=CLP1(18)+DEL*(I-18)
ENDFOR
level_minus=CLP1[[0,1,2,3,4,5,6,7,8,9]]
level_plus=CLP1[[24,25,26,27,28,29,30,31,32]]
endif

;select footpoints
if (not keyword_set(silent)) then print,'+++ FOOT POINT BEGINS +++'
if ((n_elements(xindex) gt 0) and (n_elements(xindex) eq n_elements(yindex))) then begin
  if (n_elements(zindex) ne n_elements(xindex)) then begin
    xcomm=xindex    ycomm=yindex
    zcomm=xindex*0.0+xyz0[2]  endif else begin
    xcomm=xindex    ycomm=yindex
    zcomm=zindex 
  endelse
  K1=n_elements(xindex)
endif else begin
  info={winid:0L,label:0L,number:0L,scale0:0.0,footx:fltarr(3000),footy:fltarr(3000)}
  infoptr=ptr_new(info)
  footpoints,Bz0[*,*,0],infoptr,level_plus=level_plus,level_minus=level_minus,levels=levels
  K1=(*infoptr).number
  xcomm=(*infoptr).footx[0:K1-1]
  ycomm=(*infoptr).footy[0:K1-1]
  ptr_free,infoptr
  xcomm=xcomm*DELX1+xyz0[0]
  ycomm=ycomm*DELY1+xyz0[1]
  zcomm=xcomm*0.0+xyz0[2]
endelse
if (not keyword_set(silent)) then print,'The Total Number of Foot Points:',K1
if (keyword_set(sav_pend)) then begin
  xpend=fltarr(K1)
  ypend=fltarr(K1)
  zpend=fltarr(K1)
endif
if (keyword_set(sav_nend)) then begin
  xnend=fltarr(K1)
  ynend=fltarr(K1)
  znend=fltarr(K1)
endif
if (not keyword_set(silent)) then print,'+++ FOOT POINT ENDS +++'
device,decomposed=1

;DRAWING THE CONTOUR AT THE BOTTOM
if (not keyword_set(silent)) then print,'+++ INTEGRATION BEGINS +++'
IF (keyword_set(both)) THEN BEGIN
  if (not keyword_set(silent)) then print,'Integration along the field lines in both directions.'
  FOR K=0,K1-1 DO BEGIN
  ;===================
  ;FORWARD INTEGRATION
  ;===================
    ln=0L
    xd=XCOMM[K]
    yd=YCOMM[K]
    zd=ZCOMM[K]
    tangentx=[0.0]
    tangenty=[0.0]
    tangentz=[0.0]
    bmag=[0.0]
    linex=xd
    liney=yd
    linez=zd
    xm=(xd-xyz0[0])/DELX1 		;convert to pixel position
    ym=(yd-xyz0[1])/DELY1 		;convert to pixel position
    zm=(zd-xyz0[2])/DELZ1 		;convert to pixel position
    ln=ln+1
    if (not keyword_set(silent)) then begin
      print,'==========    FOOTPOINT NO. ',string(k+1,format='(i4)'),'    =========='
      print,'FORWARD INTEGRATION'
      print,'                       X        Y        Z        Integration Steps'
      print,'Initial Position:',string(xd,format='(f9.2)'),string(yd,format='(f9.2)'),string(zd,format='(f9.2)'),string(ln-1,format='(i9)')
    endif
    s=0.0    ;independent variable of the integration equation of field lines   
    sig=+1.                    ;Forward
    while (xm GE 0.0 AND xm LE nx-1 AND $
           ym GE 0.0 AND ym LE ny-1 AND $
           zm GE 0.0 AND zm LE nz-1 AND $
           ln LE isn) do begin
      coords=[xd,yd,zd]
      dydx=differential3(s,coords)
      tangentx=[tangentx,dydx[0]]
      tangenty=[tangenty,dydx[1]]
      tangentz=[tangentz,dydx[2]]
      bmag=[bmag,Binterp]
      result=rk4(coords,dydx,s,h,'differential3')
      xd=result[0]
      yd=result[1]
      zd=result[2]
      linex=[linex,xd]
      liney=[liney,yd]
      linez=[linez,zd]
      xm=(xd-xyz0[0])/DELX1     ;convert to pixel position
      ym=(yd-xyz0[1])/DELY1     ;convert to pixel position
      zm=(zd-xyz0[2])/DELZ1     ;convert to pixel position
      ln=ln+1
      if (xd eq 0.0 and yd eq 0.0 and zd eq 0.0) then break
    endwhile
    if (ln-1 eq isn) then begin
      print,'++++++++++ Warning: The maximum ISN (Integration Step Number) reaches! Increase this number. ++++++++++'
    endif
    if (not keyword_set(silent)) then print,'Final   Position:',string(linex[ln-2],format='(f9.2)'),string(liney[ln-2],format='(f9.2)'),string(linez[ln-1],format='(f9.2)'),string(ln-1,format='(i9)')
    if (ln gt 2) then begin
      tangentx=[tangentx[2:ln-1],0.0]
      tangenty=[tangenty[2:ln-1],0.0]
      tangentz=[tangentz[2:ln-1],0.0]
      bmag=[bmag[2:ln-1],0.0]
    endif else begin
      tangentx=[0.0]
      tangenty=[0.0]
      tangentz=[0.0]
      bmag=[0.0]
    endelse 
    ;====================
    ;BACKWARD INTEGRATION
    ;==================== 
    xd=XCOMM[K]
    yd=YCOMM[K]
    zd=ZCOMM[K]
    xm=(xd-xyz0[0])/DELX1     ;convert to pixel position
    ym=(yd-xyz0[1])/DELY1     ;convert to pixel position
    zm=(zd-xyz0[2])/DELZ1     ;convert to pixel position
    if (not keyword_set(silent)) then begin
      print,'BACKWARD INTEGRATION'
      print,'                       X        Y        Z        Integration Steps'
      print,'Initial Position:',string(xd,format='(f9.2)'),string(yd,format='(f9.2)'),string(zd,format='(f9.2)'),string(ln-1,format='(i9)')
    endif
    sig=-1.                    ;backward
    while (xm GE 0.0 AND xm LE nx-1 AND $
           ym GE 0.0 AND ym LE ny-1 AND $
           zm GE 0.0 AND zm LE nz-1 AND $
           ln LE isn) do begin
      coords=[xd,yd,zd]
      dydx=differential3(s,coords)
      tangentx=[dydx[0]*sig,tangentx]    ;in this case the position is stored in a reversed order, so the tangent direction has also to be reversed
      tangenty=[dydx[1]*sig,tangenty]
      tangentz=[dydx[2]*sig,tangentz]
      bmag=[Binterp,bmag]
      result=rk4(coords,dydx,s,h,'differential3')
      xd=result[0]
      yd=result[1]
      zd=result[2]
      linex=[xd,linex]
      liney=[yd,liney]
      linez=[zd,linez]
      xm=(xd-xyz0[0])/DELX1     ;convert to pixel position
      ym=(yd-xyz0[1])/DELY1     ;convert to pixel position
      zm=(zd-xyz0[2])/DELZ1     ;convert to pixel position
      ln=ln+1
      if (xd eq 0.0 and yd eq 0.0 and zd eq 0.0) then break
    endwhile  
    if (ln-1 eq isn) then begin
      print,'++++++++++ Warning: The maximum ISN (Integration Step Number) reaches! Increase this number. ++++++++++'
    endif
    if (not keyword_set(silent)) then print,'Final   Position:',string(linex[1],format='(f9.2)'),string(liney[1],format='(f9.2)'),string(linez[1],format='(f9.2)'),string(ln-1,format='(i9)')
    if not (keyword_set(noplot)) then iplot,linex[0:ln-1],liney[0:ln-1],linez[0:ln-1],color=line_color,xtickfont_size=8,$
          ytickfont_size=8,ztickfont_size=8,thick=line_thick,/overplot
    linex=linex[1:ln-2]         ; The first and last point is derived by the last integration step, which is omitted here.
    liney=liney[1:ln-2]
    linez=linez[1:ln-2]
    tangentx=tangentx[0:ln-3]   ; The subscript is chosen in such a way to match the coordinates of the magnetic field line.
    tangenty=tangenty[0:ln-3]
    tangentz=tangentz[0:ln-3] 
    bmag=bmag[0:ln-3] 
    if (keyword_set(sav_fl)) then begin
      filename='curve_info'+string(k,format='(i4.4)')+'.sav'
      save,tangentx,tangenty,tangentz,linex,liney,linez,bmag,filename=filename
    endif
    if (keyword_set(sav_pend)) then begin
      xpend[K]=linex[0]
      ypend[K]=liney[0]
      zpend[K]=linez[0]
      save,xpend,ypend,zpend,filename='positive_ends.sav'
    endif
    if (keyword_set(sav_nend)) then begin
      xnend[K]=linex[ln-3]
      ynend[K]=liney[ln-3]
      znend[K]=linez[ln-3]
      save,xnend,ynend,znend,filename='negative_ends.sav'
    endif
    if (n_elements(linexyz) gt 0) then begin
      linexyz=fltarr(ln-2,3)
      linexyz[*,0]=linex
      linexyz[*,1]=liney
      linexyz[*,2]=linez
    endif
  ENDFOR
ENDIF ELSE BEGIN
  IF (keyword_set(forward) OR keyword_set(backward)) THEN BEGIN
    IF(keyword_set(forward)) THEN BEGIN
      sig=+1.
      print,'Integration along the field lines forward.'
    ENDIF ELSE BEGIN
      sig=-1.
      print,'Integration along the field lines backward.'
    ENDELSE
    FOR K=0,K1-1 DO BEGIN
      ln=0L
      xd=XCOMM[K]
      yd=YCOMM[K]
      zd=ZCOMM[K]
      linex=xd
      liney=yd
      linez=zd
      xm=(xd-xyz0[0])/DELX1     ;convert to pixel position
      ym=(yd-xyz0[1])/DELY1     ;convert to pixel position
      zm=(zd-xyz0[2])/DELZ1     ;convert to pixel position
      ln=ln+1
      print,'==========    FOOTPOINT NO. ',string(k+1,format='(i4)'),'    =========='
      if (sig eq 1.0) then print,'FORWARD INTEGRATION' else print,'BACKWARD INTEGRATION'
      print,'                       X        Y        Z        Integration Steps'
      print,'Initial Position:',string(xd,format='(f9.2)'),string(yd,format='(f9.2)'),string(zd,format='(f9.2)'),string(ln-1,format='(i9)')
      s=0.0    ;independent variable of the integration equation of field lines
      while (xm GE 0.0 AND xm LE nx-1 AND $
             ym GE 0.0 AND ym LE ny-1 AND $
             zm GE 0.0 AND zm LE nz-1 AND $
             ln LE isn) do begin
        coords=[xd,yd,zd]
        dydx=differential3(s,coords)
        tangentx=dydx[0]
        tangenty=dydx[1]
        tangentz=dydx[2]
        result=rk4(coords,dydx,s,h,'differential3')
        xd=result[0]
        yd=result[1]
        zd=result[2]
        linex=[linex,xd]
        liney=[liney,yd]
        linez=[linez,zd]
        xm=(xd-xyz0[0])/DELX1     ;convert to pixel position
        ym=(yd-xyz0[1])/DELY1     ;convert to pixel position
        zm=(zd-xyz0[2])/DELZ1     ;convert to pixel position
        ln=ln+1
        if (xd eq 0.0 and yd eq 0.0 and zd eq 0.0) then break
      endwhile
      if (ln-1 eq isn) then begin
        print,'++++++++++ Warning: The maximum ISN (Integration Step Number) reaches! Increase this number. ++++++++++'
      endif
      print,'Final   Position:',string(xd,format='(f9.2)'),string(yd,format='(f9.2)'),string(zd,format='(f9.2)'),string(ln-1,format='(i9)')
      iplot,linex[0:ln-1],liney[0:ln-1],linez[0:ln-1],color=line_color,xtickfont_size=8,$
            ytickfont_size=8,ztickfont_size=8,thick=line_thick,/overplot
    ENDFOR   
  ENDIF ELSE BEGIN
    print,'Integration along the field lines forward or backward according to Bz at the starting point.'
    FOR K=0,K1-1 DO BEGIN
      ln=0L
      xd=XCOMM[K]
      yd=YCOMM[K]
      zd=ZCOMM[K]
      linex=xd
      liney=yd
      linez=zd
      xm=(xd-xyz0[0])/DELX1     ;convert to pixel position
      ym=(yd-xyz0[1])/DELY1     ;convert to pixel position
      zm=(zd-xyz0[2])/DELZ1     ;convert to pixel position
      ln=ln+1
      IF(xv[2,0] LT 0.) THEN BEGIN
        sig=-1.
      ENDIF ELSE BEGIN
        sig=+1.
      ENDELSE
      print,'==========    FOOTPOINT NO. ',string(k+1,format='(i4)'),'    =========='
      if (sig eq 1.0) then print,'FORWARD INTEGRATION' else print,'BACKWARD INTEGRATION'
      print,'                       X        Y        Z        Integration Steps'
      print,'Initial Position:',string(xd,format='(f9.2)'),string(yd,format='(f9.2)'),string(zd,format='(f9.2)'),string(ln-1,format='(i9)')
      s=0.0    ;independent variable of the integration equation of field lines
      while (xm GE 0.0 AND xm LE nx-1 AND $
             ym GE 0.0 AND ym LE ny-1 AND $
             zm GE 0.0 AND zm LE nz-1 AND $
             ln LE isn) do begin
        coords=[xd,yd,zd]
        dydx=differential3(s,coords)
        tangentx=dydx[0]
        tangenty=dydx[1]
        tangentz=dydx[2]
        result=rk4(coords,dydx,s,h,'differential3')
        xd=result[0]
        yd=result[1]
        zd=result[2]
        linex=[linex,xd]
        liney=[liney,yd]
        linez=[linez,zd]
        xm=(xd-xyz0[0])/DELX1     ;convert to pixel position
        ym=(yd-xyz0[1])/DELY1     ;convert to pixel position
        zm=(zd-xyz0[2])/DELZ1     ;convert to pixel position
        ln=ln+1
        if (xd eq 0.0 and yd eq 0.0 and zd eq 0.0) then break
      endwhile
      if (ln-1 eq isn) then begin
        print,'++++++++++ Warning: The maximum ISN (Integration Step Number) reaches! Increase this number. ++++++++++'
      endif
      print,'Final   Position:',string(xd,format='(f9.2)'),string(yd,format='(f9.2)'),string(zd,format='(f9.2)'),string(ln-1,format='(i9)')
      iplot,linex[0:ln-1],liney[0:ln-1],linez[0:ln-1],color=line_color,xtickfont_size=8,$
            ytickfont_size=8,ztickfont_size=8,thick=line_thick,/overplot
    ENDFOR
  ENDELSE
ENDELSE

if not (keyword_set(noimage) and keyword_set(nocontour)) then begin
  iplot,[xyz0[0],xyz0[0]+nx*DELX1],[xyz0[1]+ny*DELY1,xyz0[1]+ny*DELY1],[xyz0[2],xyz0[2]],/overplot
  iplot,[xyz0[0]+nx*DELX1,xyz0[0]+nx*DELX1],[xyz0[1],xyz0[1]+ny*DELY1],[xyz0[2],xyz0[2]],/overplot
endif
x_cont=findgen(nx)*DELX1+xyz0[0]
y_cont=findgen(ny)*DELY1+xyz0[1]

if not (keyword_set(noimage)) then iimage,bz0[*,*,1],image_location=[xyz0[0],xyz0[1]], $
   image_dimensions=[nx*DELX1,ny*DELY1],overplot=1,zvalue=xyz0[2]
if not (keyword_set(nocontour)) then begin
  if (keyword_set(levels)) then begin
    icontour,Bz0[*,*,0],x_cont,y_cont,c_value=level_plus,min_value=-0.01,overplot=1,color=[50,50,50],zvalue=xyz0[2]
    icontour,Bz0[*,*,0],x_cont,y_cont,c_value=level_minus,max_value=0.01,overplot=1,color=[100,100,100],zvalue=xyz0[2]
  endif else begin
    icontour,Bz0[*,*,0],x_cont,y_cont,c_value=level_plus,min_value=-0.01,overplot=1,color=[50,50,50],zvalue=xyz0[2]
    icontour,Bz0[*,*,0],x_cont,y_cont,c_value=level_minus,max_value=0.01,overplot=1,color=[100,100,100],zvalue=xyz0[2]
    icontour,Bz0[*,*,0],x_cont,y_cont,c_value=[0,0.01],C_LINESTYLE=[3],color=[0,0,0],overplot=1,zvalue=xyz0[2]
  endelse
endif
if (keyword_set(sav_start)) then begin
  xindex=xcomm
  yindex=ycomm
  zindex=zcomm
  save,xindex,yindex,zindex,filename='footpoints_position.sav'
endif  

if (not keyword_set(silent)) then print,'+++ DRAWING FIELD LINES ENDS +++'
END
