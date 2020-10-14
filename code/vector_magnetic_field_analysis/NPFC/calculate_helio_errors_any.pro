PRO CALCULATE_HELIO_ERRORS_ANY,Bl,Btr,phil,Bz,Bh,phi,B0,L0,P,Bc,Lc,sig_Bl,$
sig_Btr,sig_phil,sig_Bh,sig_Bz,sig_B

; PURPOSE: Using the uncertainties in the longitudinal and the
; transverse magnetic field components, this routine calculates the
; uncertainties in the vertical and the horizontal (heliographic)
; magnetic field components, sig_Bz and sig_Bh, respectively. This is
; done by error propagation in the transformation matrix introduced by
; Gary & Hagyard (1990). If desired, the routine also calculates the
; uncertainties in the individual horizontal magnetic field
; components, sig_Bx and sig_By, respectively.
; PROGRAMMER: Manolis K. Georgoulis (JHU/APL, 10/13/05)

r1=where(abs(Bl) gt 5.d3 or Btr gt 5.d3 or phil gt 2.*!dpi,ico)
if ico gt 0. then begin & Bl(r1)=0. & Btr(r1)=0. & phil(r1)=0. &endif

a11=-sin(B0)*sin(P)*sin(Lc-L0) + cos(P)*cos(Lc-L0)
a12=sin(B0)*cos(P)*sin(Lc-L0) + sin(P)*cos(Lc-L0)
a13=-cos(B0)*sin(Lc-L0)
a21=-sin(Bc)*(sin(B0)*sin(P)*cos(Lc-L0) + cos(P)*sin(Lc-L0)) $
             -cos(Bc)*cos(B0)*sin(P)
a22=sin(Bc)*(sin(B0)*cos(P)*cos(Lc-L0) - sin(P)*sin(Lc-L0)) $
                  +cos(Bc)*cos(B0)*cos(P)
a23=-cos(B0)*sin(Bc)*cos(Lc-L0) + sin(B0)*cos(Bc)
a31=cos(Bc)*(sin(B0)*sin(P)*cos(Lc-L0) + cos(P)*sin(Lc-L0)) $
            -sin(Bc)*cos(B0)*sin(P)
a32=-cos(Bc)*(sin(B0)*cos(P)*cos(Lc-L0) - sin(P)*sin(Lc-L0)) $
                   +sin(Bc)*cos(B0)*cos(P)
a33=cos(Bc)*cos(B0)*cos(Lc-L0) + sin(Bc)*sin(B0)
B=sqrt(Bl^2.+Btr^2)
sin_theta=Btr/B & cos_theta=Bl/B
r1=where(B eq 0.,ico) & if ico gt 0. then begin & sin_theta(r1)=0. & cos_theta(r1)=0. &endif
restructure_angle,sin_theta,cos_theta,thl

sig_B=(abs(Bl)/B)*sig_Bl + (Btr/B)*sig_Btr
r1=where(B eq 0.,ico) & if ico gt 0. then sig_B(r1)=0.
sig_Bksi=abs(cos(phil))*sig_Btr + Btr*abs(sin(phil))*sig_phil
sig_Beta=abs(sin(phil))*sig_Btr + Btr*abs(cos(phil))*sig_phil
sig_Bx=abs(a11)*sig_Bksi + abs(a12)*sig_Beta + abs(a13)*sig_Bl
sig_By=abs(a21)*sig_Bksi + abs(a22)*sig_Beta + abs(a23)*sig_Bl
sig_Bz=abs(a31)*sig_Bksi + abs(a32)*sig_Beta + abs(a33)*sig_Bl
r1=where(B eq 0.,ico) 
if ico gt 0. then begin & sig_Bx(r1)=0. & sig_By(r1)=0. & sig_Bz(r1)=0. &endif

Bksi=Btr*cos(phil) & Beta=Btr*sin(phil)
dBh_Bksi=(1./Bh)*((a11^2.+a21^2.)*Bksi + (a11*a13+a21*a23)*Bl + $
                  (a11*a12+a21*a22)*Beta)
dBh_Beta=(1./Bh)*((a12^2.+a22^2.)*Beta + (a11*a12+a21*a22)*Bksi + $
                  (a12*a13+a22*a23)*Bl)
dBh_Bl=(1./Bh)*((a13^2.+a23^2.)*Bl + (a11*a13+a21*a23)*Bksi + $
                (a12*a13+a22*a23)*Beta)
sig_Bh=abs(dBh_Bksi)*sig_Bksi + abs(dBh_Beta)*sig_Beta + abs(dBh_Bl)*sig_Bl
r1=where(sig_Bh*0. ne 0.,ico) & if ico gt 0. then sig_Bh(r1)=0.

RETURN
END





