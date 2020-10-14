PRO LOS_IMAGE_SOLUTION,phi_res_im,phi1_im,phi2_im,phi_ref,phi_res_los,$
idim1,idim2

; PURPOSE: From the solution phi_res_im of the heliographic azimuth on
; the image plane, find the solution phi_res_los of the LOS azimuth,
; where (phi1_im,phi2_im) are the two possible solutions of the
; heliographic azimuth and phi_ref is the unresolved LOS azimuth
; PROGRAMMER: Manolis K. Georgoulis (JHU/APL, 10/13/05)

phi_res_los=dblarr(idim1,idim2)
r1=where(phi_res_im eq phi1_im,ic) 
if ic gt 0. then phi_res_los(r1)=phi_ref(r1)
r1=where(phi_res_im eq phi2_im,ic) 
if ic gt 0. then phi_res_los(r1)=phi_ref(r1)+!dpi
r1=where(phi_res_los gt 2.*!dpi,ic) 
if ic gt 0. then phi_res_los(r1)=phi_res_los(r1)-2.*!dpi

RETURN
END
