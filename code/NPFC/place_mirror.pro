PRO PLACE_MIRROR,im,x1,x2,y1,y2,mr

; PURPOSE: Place an image mr in specified locations of an image
; im. The edge locations in im where mr is to be placed are
; (x1,y1) and (x2,y2)
; PROGRAMMER: Manolis K. Georgoulis (JHU/APL, 10/13/05)

nxa=intarr(2) & nya=intarr(2)
res=size(im(x1:x2,y1:y2)) & nxa(0)=res(1) & nya(0)=res(2)
res=size(mr) & nxa(1)=res(1) & nya(1)=res(2)

nx=min(nxa) & ny=min(nya)
im(x1:x1+nx-1,y1:y1+ny-1)=mr(0:nx-1,0:ny-1)

RETURN
END
