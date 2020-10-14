PRO DIV_FREE_SOLUTION_RING,mag,scale,a,mirror=mirror

; PURPOSE: Double precision version of the routine suggested by
; B. J. LaBonte. It calculates the solution a satisfying the gauge
; conditions of Chae (2001) by means of the fast Fourier transform
; FEATURES: Pads the original image with zeros (if /mirror is not
; set), or it pads it with a mirror image of the original if /mirror
; is set
; PROGRAMMERS: B. J. LaBonte & M. K. Georgoulis (JHU/APL, 10/13/05)

res=size(mag) & id1=res(1) & id2=res(2)
mag_ref=mag
if keyword_set(mirror) then EXPAND_IMAGE,mag_ref,mag,id1,id2,stx,sty,/mirror else $
                            EXPAND_IMAGE,mag_ref,mag,id1,id2,stx,sty

; Get the basic transform
nx = N_ELEMENTS(mag(*,0))
ny = N_ELEMENTS(mag(0,*))
ftm = FFT( mag, -1, /double)

; Multiply by i
ftr = CONJ( COMPLEX( IMAGINARY(ftm), FLOAT(ftm) ) )

; Generate wavenumber images
kx = REBIN( DINDGEN(nx), nx, ny ) - DOUBLE(nx-2)/2.
kx = SHIFT( kx, -(nx-2)/2, 0 )
kx = -kx * 2. * !DPI / DOUBLE(nx)
ky = REBIN( TRANSPOSE( DINDGEN(ny) ), nx, ny ) - DOUBLE(ny-2)/2.
ky = SHIFT( ky, 0, -(ny-2)/2 )
ky = -ky * 2. * !DPI / DOUBLE(ny)
k2 = kx^2 + ky^2
k2(0,0) = 1.

; Get vector potential from inverse transform
apx = DOUBLE( FFT( (ky/k2) * ftr, 1,/double ) )
apy = DOUBLE( FFT( (-kx/k2) * ftr, 1,/double ) )

a = [ [[apx]], [[apy]] ] * scale
a = [ [[apx(stx:stx+id1-1,sty:sty+id2-1)]], $
      [[apy(stx:stx+id1-1,sty:sty+id2-1)]] ] * scale
mag=mag_ref

RETURN
END
