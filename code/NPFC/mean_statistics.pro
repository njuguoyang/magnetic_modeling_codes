PRO MEAN_STATISTICS,im,ibox

; PURPOSE: Find the average of the closest Cartesian neighborhood of a
; given image im where the neighborhood is defined by a linear
; distance imbox from each point of im
; PROGRAMMER: Manolis K. Georgoulis (JHU/APL, 10/13/05) 

res=size(im) & id1=res(1) & id2=res(2)
im_ref=im
for i=-ibox,ibox do for j=-ibox,ibox do im=im+shift(im_ref,i,j)
im=im/(2.*ibox+1.)^2.

RETURN
END



