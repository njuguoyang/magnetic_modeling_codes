pro derives2,dx,dy,psi,psi_xx,psi_yy,psi_xy,psi_x,psi_y,backward=backward
;
; assumes constant spacing  
;
  type=size(psi,/type)
  NAN=!values.F_nan
  if (type eq 5) then Nan=!values.D_nan
;
  sz=size(psi)
  nx=sz[1]
  ny=sz[2]
;
; compute second order derives of psi at half grid points
;  psi  :     1     2     3     4     5
;  psi_x:        1     2     3     4
; psi_xx:           1     2     3
;
;   deal with interior points
  if keyword_set(backward) then begin
     psi_x=(psi-shift(psi,1,0))/dx
     psi_y=(psi-shift(psi,0,1))/dy
     psi_x[0,*]=nan
     psi_y[*,0]=nan
  endif else begin
    psi_x=(shift(psi,-1,0)-psi)/dx
    psi_y=(shift(psi,0,-1)-psi)/dy
    psi_x[nx-1,*]=nan
    psi_y[*,ny-1]=nan
  endelse
;  stop
;
  psi_xy=((psi-shift(psi,-1,0)-shift(psi,0,-1)+shift(psi,-1,-1))/(dx*dy))
  psi_xy[nx-1,*]=nan
  psi_xy[*,ny-1]=nan
;
  psi_xx=(shift(psi,-1,0)-2*psi+shift(psi,1,0))/(dx^2)
  psi_yy=(shift(psi,0,-1)-2*psi+shift(psi,0,1))/(dy^2)
;
;   deal with boundary points pg. 914 Abramowitz
;
  psi_xx[0,*]=NAN
  psi_xx[nx-1,*]=NAN
;
  psi_yy[*,0]=NAN
  psi_yy[*,ny-1]=NAN
;
;  stop
;
;
return

end
