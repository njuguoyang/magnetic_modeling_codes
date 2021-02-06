PRO CALCULATE_ERRORS,sig_Bz,sig_Bh,sig_B,Bz,Bh,err_dBz,err_Fz,err_Jmin,lamda

; PURPOSE: Here we calculate the relative errors in the vertical
; magnetic field gradient (err_dBz), the vertical Lorentz force
; (err_Fz) and the cross-field electric current density (err_Jhp)
; (1) The calculation is done with the heliographic field components
;     on the image plane
; (2) The Gaussian unit system is used
; (3) lamda is the pixel size in cm
; PROGRAMMER: Manolis K. Georgoulis (JHU/APL, 10/13/05)

B=sqrt(Bz^2.+Bh^2.)
err_B=sig_B/B 
r1=where(B eq 0.,ico) & if ico gt 0. then err_B(r1)=0.
err_dBz=sig_Bz/abs(Bz) + sig_Bh/Bh + err_B
r1=where(err_dBz*0. ne 0.,ico) & if ico gt 0. then err_dBz(r1)=0.

dB_z=2.*abs(Bz)*B/(lamda*Bh)
sig_dBz=(2.*B/(lamda*Bh))*sig_Bz + 2.*abs(Bz)/(lamda*Bh)*sig_B + $
        2.*abs(Bz)*B/(lamda*Bh^2.)*sig_Bh
Fz=(1./(!dpi*lamda))*Bh*abs(Bz) + (B/(4.*!dpi))*abs(dB_z)
sig_Fz=(1./(!dpi*lamda))*abs(Bz)*sig_Bh + (1./(!dpi*lamda))*Bh*sig_Bz + $
       (1./(4.*!dpi))*abs(dB_z)*sig_B + (1./(4.*!dpi))*B*sig_dBz
err_Fz=sig_Fz/Fz
r1=where(err_Fz*0. ne 0.,ico) & if ico gt 0. then err_Fz(r1)=0.

;err_Fz=err_B + err_dBz
err_Jmin=err_Fz+sig_Bh/Bh
r1=where(Bh eq 0.,ico) & if ico gt 1. then err_Jmin(r1)=0.

RETURN
END
