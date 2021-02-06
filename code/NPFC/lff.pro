;+
;NAME:
;     LFF
;PURPOSE:
;     Computes a linear (constant alpha) force free magnetic field from
;     either the observed vertical or line-of-sight field.
;CATEGORY:
;CALLING SEQUENCE:
;     bff = lff(bz,b0=b0,pangle=pangle,lat=lat,cmd=cmd)
;INPUTS:
;     bz = line-of-sight or vertical magnetic field (see /normal below).
;          If bz is a structure it is assumed to be a standard magnetic
;          field structure from bfits.pro.  In that case, the LOS or
;          vertical field and the pointing information (b0, p, lat, cmd, 
;          pixel) is taken from the structure.
;OPTIONAL INPUT PARAMETERS:
;KEYWORD PARAMETERS
;     /normal = Nonzero to specify that Bz is the normal component (i.e.
;               heliographic B_z rather than image-plane (LOS) B_z).
;     /image = Return the field in the image coordinates rather than the 
;              heliographic coordinates
;     b0 = solar b0 angle in radians
;     pangle  = solar p angle in radians
;     lat = latitude in radians, can be an array if the array is the same 
;           size as the bz array.
;     cmd = cmd in radians, can be an array if the array is the same 
;           size as the bz array.
;     z = Array of heights (Mm) at which to evaluate force-free field. 
;         Default = [0.0].
;     alpha = Force free alpha in units of 1/(x pixels).  (default = 0.0)
;     pixel = Size of pixel in arcsec.  Two element array specifying
;             [dx,dy]; a scalar number applies to both. Default is 5.7".
;     guard = Add a guard ring around the input data.  The value gives the
;             multiplicative increase in size.  Default = 1 (no guard ring).
;             Guard must be integer.
;     xh = heliographic x-coordinates of pixels in Mm (Output).
;     yh = heliographic y-coordinates of pixels in Mm (Output).
;          (zh=0 on the photosphere)
;     dA = area of pixel projected on photosphere (Mm^2). This is the area of
;          the image plane pixel at the Sun divided by mu=cos(distance from 
;          center). Should also be d(xh)*d(yh)*(parallelogram factor) =
;          d(xh)*d(yh)*sin(corner angle of grid). Output.
;     /help = Help.
;     /quiet = work quietly.
;
;   Horizontal field for the large scale Fourier terms:  if not set, the
;   arbitrary constant is set to 0.0.
;   These assume that the 180 degree ambiguity is resolved.  If it is not,
;   they should not be specified.  If all are set, bx and by are used.
;   btrans and azim are inconsistent with /normal.
;
;     BX      = The heliographic x component of the field.  Used to compute
;               the large scale fourier terms when |alpha| is large.
;     BY      = The heliographic y component of the field.  Used to compute
;               the large scale fourier terms when |alpha| is large.
;     BTRANS  = The magnitude of the transverse magnetic field.  Used to 
;               compute the large scale fourier terms when |alpha| is large.
;     AZIM    = The azimuth of the transverse magnetic field in DEGREES
;               CCW from terrestrial north.  Used to compute the large scale 
;               fourier terms when |alpha| is large.
;
;OUTPUTS:
;     bff = x,y,z force-free field (nx,ny,nz,3)
;COMMON BLOCKS:
;SIDE EFFECTS:
;RESTRICTIONS:
;
;     Image-plane x, y, and z are in the terrestrial West and North
;     directions and towards the observer.
;
;     Heliographic x, y, z are in the directions of solar West
;     (increasing longitude), solar North (increasing latitude),
;     and vertically upwards from the photosphere.
;
;PROCEDURE:
;     The calculation is done in the image plane coordinates.
;EXAMPLE:
;     f = ff(mag.b_long,b0=mag.point.b0,p=mag.point.p,lat=mag.point.lat, $
;            cmd=mag.point.cmd,pixel=mag.point.pix_size,btrans=mag.b_trans, $
;            azim=mag.b_azim,alpha=get_alpha(mag))
;MODIFICATION HISTORY:
; 900xxx Yuhong Fan: Written.
; 9109xx A. McClymont: 
;    Added comments and refs, cleaned up code, checked transformations.
;    Added keywords and Biz, Bh output.  Note that output Bz may be
;    different from input Bz due to fudging for zero net flux.
;    For consistency the output value of B[ih]z should be used in
;    conjunction with the output B[ih][xy] values in calculations.
; 911014 ANM    Remove guard stuff, assume already done.
; 940428 T. Metcalf
;    Converted Sandy and Yuhong's potential92.pro to compute the linear force 
;    free field instead of the potential field.  Uses the equations of 
;    Alissandrakis, A&A, 100, 197, 1981 and Gary, ApJS, 69, 323, 1989.
; 960216 T. Metcalf
;    Put guard stuff back in as an option
;- 

; SHOULD USE REFLECTION TO MAKE FIELD LINES ENCLOSED RATHER THAN PERIODIC.
; IF Z!=0 AND IMAGE PLANE, IS THE PLANE COMPUTED ALONG THE LOS I.E. NOT
; VERTICALLY ABOVE THE PHOTOSPHERE?
; WHY DOES LAST Z TAKE SO MUCH LONGER THAN FIRST?

function lff,b0=b0,pangle=p,cmd=lcin,lat=bcin,BzIn, xh=xh,yh=yh,dA=dA, $
		HELP=help,NORMAL=normal,PIXEL=pixel,Z=z, $
		STOP=stop,QUIET=quiet,IMAGE=image,TICK=tick,alpha=inalpha, $
                bx=bxin,by=byin,btrans=btrans,azim=azim,guard=guard

;******************************************************************************

IF KEYWORD_SET(HELP) THEN BEGIN
   PRINT,'FUNCTION lff,b0=b0,pangle=pangle,lat=lat,cmd=cmd,Bz,xh=xh,yh=yh, $'
   PRINT,'dA=dA,HELP=help,NORMAL=normal,PIXEL=pixel,Z=z,alpha=alpha, $'
   PRINT,'image=image'
   RETURN,-1
ENDIF
ss = SIZE(BzIn)

if n_elements(lcin) GT 0 then lc = lcin
if n_elements(bcin) GT 0 then bc = bcin

if ss(n_elements(ss)-2) EQ 8 then begin  ; Is it a structure?
   if keyword_set(normal) then bz = BzIn.bz else bz = BzIn.b_long
   if n_elements(b0) LE 0 then b0 = BzIn.point.b0
   if n_elements(p) LE 0 then p = BzIn.point.p
   if tag_index(BzIn,'LATITUDE') GE 0 then begin
      if (n_elements(bc) LE 0) and (n_elements(bz) EQ n_elements(BzIn.latitude)) then $
         bc = BzIn.latitude*!dtor
   endif
   if tag_index(BzIn,'CMD') GE 0 then begin
      if (n_elements(lc) LE 0) AND (n_elements(bz) EQ n_elements(BzIn.cmd)) then $
         lc = BzIn.cmd*!dtor
   endif
   if n_elements(bc) LE 0 then bc = BzIn.point.lat
   if n_elements(lc) LE 0 then lc = BzIn.point.cmd
   if n_elements(pixel) LE 0 then pixel = BzIn.point.pix_size
endif else bz=bzin

ss = SIZE(Bz)
Nxin = ss(1) & Nyin = ss(2)
Nx1in=Nxin-1 & Ny1in=Nyin-1 ; TEMP

if n_elements(guard) LE 0 then guard = 1
use_guard_ring = (n_elements(guard) EQ 1L) AND (long(guard) GT 1L)

if use_guard_ring then begin
   if NOT keyword_set(quiet) then message,/info,'Guard factor is '+string(guard)
   bzbig = make_array(long(nxin*long(guard)), $
                      long(nyin*long(guard)), $
                      type=ss(n_elements(ss)-2),value=0)
   bzbig(0:nx1in,0:ny1in) = bz
   bz = temporary(bzbig)
   if n_elements(bc) EQ nxin*nyin then begin
      bcbig = make_array(long(nxin*long(guard)), $
                         long(nyin*long(guard)), $
                         /float,value=mean(bc))
      bcbig(0:nx1in,0:ny1in) = bc
      bc = temporary(bcbig)
   endif
   if n_elements(lc) EQ nxin*nyin then begin
      lcbig = make_array(long(nxin*long(guard)), $
                         long(nyin*long(guard)), $
                         /float,value=mean(lc))
      lcbig(0:nx1in,0:ny1in) = lc
      lc = temporary(lcbig)
   endif
endif

ss = SIZE(Bz)
Nx = ss(1) & Ny = ss(2)
Nx1=Nx & Ny1=Ny ; TEMP

; These are the transformation matrix elements given by eq. (1)
; of Gary and Hagyard, Solar Phys. 126, 21 (1990).
; The components of the magnetic field in heliographic coordinates (H)
; are given in terms of the image-plane components (I) by BH = A BI, where
; the B's are column vectors and A = [axx axy axz]
;				     [ayx ayy ayz]
;				     [azz azy azz]
; Note that A is orthogonal, so Inverse(A) = Transpose(A).
; Also note coordinates transform same as fields (see matrix c below).
; Also note last column of A gives direction cosines of LOS in helio-coords.

axx = - SIN(b0)*SIN(p)*SIN(lc) + COS(p)*COS(lc)
axy = + SIN(b0)*COS(p)*SIN(lc) + SIN(p)*COS(lc)
axz = - COS(b0)*SIN(lc)
ayx = - SIN(bc)*( SIN(b0)*SIN(p)*COS(lc) + COS(p)*SIN(lc) ) $
      - COS(bc)*COS(b0)*SIN(p)
ayy = + SIN(bc)*( SIN(b0)*COS(p)*COS(lc) - SIN(p)*SIN(lc) ) $
      + COS(bc)*COS(b0)*COS(p)
ayz = - COS(b0)*SIN(bc)*COS(lc) + SIN(b0)*COS(bc)
azx = + COS(bc)*( SIN(b0)*SIN(p)*COS(lc) + COS(p)*SIN(lc) ) $
      - SIN(bc)*COS(b0)*SIN(p)
azy = - COS(bc)*( SIN(b0)*COS(p)*COS(lc) - SIN(p)*SIN(lc) ) $
      + SIN(bc)*COS(b0)*COS(p)
azz = + COS(bc)*COS(b0)*COS(lc) + SIN(bc)*SIN(b0)

; These are the transformation matrix elements given by eq. (2)
; of Gary and Hagyard, Solar Phys. 126, 21 (1990).
; The image-plane coordinates xI are related to the heliographic coordinates
; xH by xI = C xH, where C = [cxx cxy]
;			     [cyx cyy]
; THIS ASSUMES zH=0 I.E. ON PHOTOSPHERE, so 2x2 matrix instead of 3x3.
; Note that zI = czx*xH + czy*yH, which they do not define, but
; it exists (and is nonzero).

cxx = axx & cxy = ayx & cyx = axy & cyy = ayy

;;IF NOT KEYWORD_SET(NORMAL) AND NOT KEYWORD_SET(quiet) THEN $
;;  PRINT, 'For this LOS data the net flux should be zero only if the field', $
;;	 'is really force-free with constant alpha.'

fa = FFT(bz,-1) ; Fourier representation of B_l (or B_n)
if NOT keyword_set(quiet) then begin
   print
   if keyword_set(normal) then type = 'normal' else type='LOS'
   message,/info,strcompress('Net '+type+' flux = '+string(fa(0,0)))
   print
endif
fa(0,0) = 0.0	; make sure net flux is zero

; image plane spatial frequencies of Fourier components (k_xi,k_yi)
kxi = 2*!PI/Nx1*SHIFT(INDGEN(Nx1)-(Nx1-1)/2,-(Nx1-1)/2)#REPLICATE(1,Ny1)
kyi = 2*!PI/Ny1*REPLICATE(1,Nx1)#SHIFT(INDGEN(Ny1)-(Ny1-1)/2,-(Ny1-1)/2)
; Originally in units of 2*!PI*(cycles per pixel) i.e. radians per pixel,
; so assumed square pixels. Convert to e.g. (Mm)^-1, so we know units for z.
CASE N_ELEMENTS(PIXEL) OF
0: BEGIN & dxi = 4*SQRT(2) & dyi = 4*SQRT(2) & END
1: BEGIN & dxi = pixel(0) & dyi = pixel(0) & END
2: BEGIN & dxi = pixel(0) & dyi = pixel(1) & END
ENDCASE
dxi = abs((149e3/206265.0)*dxi)
dyi = abs((149e3/206265.0)*dyi) ; pixel size in Mm.
kxi = kxi/dxi & kyi = kyi/dyi ; radians per Mm

if n_elements(inalpha) NE 1 then alpha = 0. $
else alpha=float(inalpha)/abs(dxi)  ; assume inalpha is in inverse pixels
alpha2 = complex(alpha^2,0.)

; heliocentric spatial frequencies (k_x,k_y,k_z)
kx = cxx*kxi + cyx*kyi & ky = cxy*kxi + cyy*kyi & kz2 = (kx^2+ky^2) 

k = sqrt(complex(kz2,0.0)-alpha2)

; spatial frequency in the data vector (B_l or B_n) direction
IF KEYWORD_SET(NORMAL)	THEN BEGIN
   kl2 = COMPLEX(0.0, kz2)
   kl2b = kl2 
ENDIF ELSE BEGIN
   kl2 = COMPLEX(axz*(k*kx-alpha*ky) + ayz*(k*ky+alpha*kx), azz*kz2)
   kl2b = COMPLEX(axz*(-k*kx-alpha*ky) + ayz*(-k*ky+alpha*kx), azz*kz2)
ENDELSE
kl2(0,0) = 1.0  ; (0,0) is singular, dont divide by zero
kl2b(0,0) = 1.0 ; (0,0) is singular, dont divide by zero

; Deal with large scale Fourier components with infinite energy.
; Requires an ambiguity resolved magnetogram to do this right.  If not, 
; just set C to 0.  See Gary, ApJS, 69, 323, 1989, although his equation
; for C is off a bit.  For imaginary k, there are two solutions to
; C for +\- k, however, the two solutions are equal when C is real.

bad = where(imaginary(k) NE 0.0,nbad)

; kz2(0,0) is always 0, but should not be considered bad.

if nbad GT 1 then begin
   bad = bad(1:nbad-1)
   nbad = nbad - 1
endif

if keyword_set(bxin) AND keyword_set(byin) then begin
   bx = bxin
   by = byin
endif else if keyword_set(btrans) AND keyword_set(azim) and $
          NOT keyword_set(normal) then begin
   ; Assume azimuth is CCW from terrestrial north
   bxi = -btrans*sin(azim*!dtor)
   byi = +btrans*cos(azim*!dtor)
   bx = axx*bxi + axy*byi + axz*bz(0:nx1-1,0:ny1-1) ;transform to heliographic field
   by = ayx*bxi + ayy*byi + ayz*bz(0:nx1-1,0:ny1-1)
   bxi=0   ; save some memory
   byi=0   ; save some memory
endif
if keyword_set(bx) and keyword_set(by) and $
   (nbad GT 0 or NOT keyword_set(quiet)) then begin
   if use_guard_ring then begin
      bxbig = make_array(long(nxin*long(guard)), $
                         long(nyin*long(guard)), $
                         type=ss(n_elements(ss)-2),value=0)
      bxbig(0:nx1in,0:ny1in) = bx
      bx = temporary(bxbig)
      bybig = make_array(long(nxin*long(guard)), $
                         long(nyin*long(guard)), $
                         type=ss(n_elements(ss)-2),value=0)
      bybig(0:nx1in,0:ny1in) = by
      by = temporary(bybig)
   endif
   fbx = FFT(bx,-1)
   fby = FFT(by,-1)
   junk = check_math(1,1) ;Turn math errors off avoiding over and underflows
   if nbad GT 0 then begin
      Cff = float((-kl2/kz2) * (fbx*kx + fby*ky)/(abs(k)*fa))
      Cff(0,0) = 0.0
   endif
   if NOT keyword_set(quiet) then begin
      calpha = float((kl2/kz2) * (fby*kx - fbx*ky)/fa)
      calpha(0,0) = 0.0
   endif
   junk = check_math(0,0) ; Turn math errors back on
   if nbad GT 0 then C = median(Cff(bad)) else C = 0.0
   ; If the linear force-free field is a good approximation, C and alpha
   ; should be constant.
   if NOT keyword_set(quiet) then begin
      message,/info,'If the linear force free approximation is good, alpha'
      if nbad GT 0 then message,/info,'and C should be constant' else $
         message,/info,'should be constant'
      message,/info,strcompress('alpha = '+ $
                    string(median(calpha*abs(dxi)))+' +\- '+ $
                    string(stdev(calpha*abs(dxi))))
      if nbad GT 0 then message,/info,'    C ='+ $
                    strcompress(' '+string(C)+' +\- '+ $
                    string(stdev(Cff))+', # inf. energy terms = '+ $
                    string(nbad)+ ' of '+string(n_elements(k)))
      print
   endif
   Cff = 0.0     ; save some memory
   calpha = 0.0  ; save some memory
endif else begin
   C = 0.0   ; default C is zero (presumably ambiguity not resolved)
endelse
if nbad GT 0 then begin
   a1 = complex(1.0,+C)/2.0   ; assumes C is real
   a2 = complex(1.0,-C)/2.0
endif

; Compute heliographic field

iphihat0 = fa/kl2
if nbad GT 0 then iphihat0B = fa(bad)/kl2b(bad)

IF N_ELEMENTS(z) EQ 0 THEN z = 0.0 ; photospheric height if not given
Nz = N_ELEMENTS(z) & B = FLTARR(Nx,Ny,Nz,3,/NOZERO)
FOR iz=0,Nz-1 DO BEGIN
  IF KEYWORD_SET(tick) THEN PRINT,STRING("15B),Nz-iz,FORM='(a,i,$)'
  iphihat = iphihat0*EXP(-k*z(iz)) ; get field at height z
  fbx = (k*kx-alpha*ky)*iphihat
  fby = (k*ky+alpha*kx)*iphihat
  fbz = complex(0.,1.)*(kz2)*iphihat
  if nbad GT 0 then begin
    ekz = EXP(+k(bad)*z(iz)) 
    fbx(bad) = a1*fbx(bad)+a2*(-k(bad)*kx(bad)-alpha*ky(bad))*iphihat0B*ekz
    fby(bad) = a1*fby(bad)+a2*(-k(bad)*ky(bad)+alpha*kx(bad))*iphihat0B*ekz
    fbz(bad) = a1*fbz(bad)+a2*complex(0.,1.)*(kz2(bad))*iphihat0B*ekz
  endif  
  B(*,*,iz,0) = FLOAT(FFT(fbx,1))
  B(*,*,iz,1) = FLOAT(FFT(fby,1))
  B(*,*,iz,2) = FLOAT(FFT(fbz,1))
ENDFOR
IF KEYWORD_SET(tick) THEN PRINT

xi = dxi*FINDGEN(Nx1)#REPLICATE(1,Ny1) & yi = dyi*REPLICATE(1,Nx1)#FINDGEN(Ny1)
zi = - ( azx*xi + azy*yi )/azz ; "height" of photosphere along LOS
xh = axx*xi + axy*yi + axz*zi & yh = ayx*xi + ayy*yi + ayz*zi
cosa = SIN(b0)*SIN(bc) + COS(b0)*COS(bc)*COS(lc) ; ang. dist. from sun centre
dA = dxi*dyi/cosa ; area of photosphere Mm^2 per pixel

;;; compute potential field energy
;;iz = WHERE(z EQ 0.0) & iz = iz(0)
;;IF NOT KEYWORD_SET(quiet) AND iz GE 0 THEN BEGIN
;;  Bhx = B(*,*,iz,0) & Bhy = B(*,*,iz,1) & Bhz = B(*,*,iz,2)
;;  ; "Low" method (NASA CP2374)
;;  ; forces on photosphere -- check validity of solution
;;  F0 = 0.5*TOTAL(Bhx^2+Bhy^2+Bhz^2)
;;  Fx=TOTAL(Bhx*Bhz) & Fy=TOTAL(Bhy*Bhz) & Fz=0.5*TOTAL(Bhz^2-Bhx^2-Bhy^2)
;;  PRINT,'Total forces of potential field on photosphere cf. average B^2:'
;;  PRINT,'	(to check validity of solution)'
;;  PRINT,'	[Fx,Fy,Fz]/F0=',Fx/F0,Fy/F0,Fz/F0
;;
;;  eLow = 1e-6/(4*!PI)*dA*TOTAL((xh*Bhx+yh*Bhy)*Bhz) ; in units of 1e30 ergs
;;  PRINT,'Potential field energy (Low method) =',eLow,' *1e30 erg'
;;
;;  ; "Semel" method
;;  phi = FLOAT(FFT(-COMPLEX(0,1)*iphihat0,1))
;;  eSemel = -1e-6/(8*!PI)*dA*TOTAL(phi*Bhz) ; hopefully in same units
;;  PRINT,'Potential field energy (Semel method) =',eSemel,' *1e30 erg'
;;  PRINT,'If these energies disagree significantly, it may mean that it is'
;;  PRINT,'necessary to put a wider guard ring around the observed region.'
;;ENDIF

IF KEYWORD_SET(stop) THEN STOP,'Stopped in lff.pro, type .c to continue'

IF KEYWORD_SET(image) THEN BEGIN ; image-plane components
   ; assume only wanted if Nz=1
   IF Nz GT 1 THEN MESSAGE,'/image inconsistent with nz GT 1'
   Bix = axx*B(*,*,0,0) + ayx*B(*,*,0,1) + azx*B(*,*,0,2)
   Biy = axy*B(*,*,0,0) + ayy*B(*,*,0,1) + azy*B(*,*,0,2)
   Biz = axz*B(*,*,0,0) + ayz*B(*,*,0,1) + azz*B(*,*,0,2)
   B(*,*,0,0) = Bix & B(*,*,0,1) = Biy & B(*,*,0,2) = Biz
ENDIF

; make B a 2+1 D array if z is scalar
IF NOT IsArray(z) THEN B = REFORM(B,Nx,Ny,3)

if use_guard_ring then begin
   B  = temporary(B(0:nx1in,0:ny1in,*,*))
   xh = temporary(xh(0:nx1in,0:ny1in))
   yh = temporary(yh(0:nx1in,0:ny1in))
;   dA = temporary(dA(0:nx1in,0:ny1in))
endif

return,B

END




