;+
; NAME :
;   CORNER
; PURPOSE :
;   Get the values of 8 vertices in the 3 dimensional data cube
; CATEGORY :
;
; CALLING SEQUENCE :
;   Called by fieldline3d.pro, which is used for drawing 3 dimensional field lines
; INPUTS :
;
; OUTPUTS :
;
; COMMON BLOCKS :
; MODIFICATION HISTORY :
;   originally coded by M.T.Song by FORTRAN
;   2007.02 Guo Yang translated these codes from FORTRAN to IDL
;-
PRO CORNER,n1,n2,n3,Bx,By,Bz,xv
n1p=n1+1
n2p=n2+1
n3p=n3+1
xv[0,0]=Bx[n1,n2,n3]
xv[0,1]=Bx[n1p,n2,n3]
xv[1,0]=By[n1,n2,n3]
xv[1,1]=By[n1p,n2,n3]
xv[2,0]=Bz[n1,n2,n3]
xv[2,1]=Bz[n1p,n2,n3]
xv[0,2]=Bx[n1,n2p,n3]
xv[0,3]=Bx[n1p,n2p,n3]
xv[1,2]=By[n1,n2p,n3]
xv[1,3]=By[n1p,n2p,n3]
xv[2,2]=Bz[n1,n2p,n3]
xv[2,3]=Bz[n1p,n2p,n3]
xv[0,4]=Bx[n1,n2,n3p]
xv[0,5]=Bx[n1p,n2,n3p]
xv[1,4]=By[n1,n2,n3p]
xv[1,5]=By[n1p,n2,n3p]
xv[2,4]=Bz[n1,n2,n3p]
xv[2,5]=Bz[n1p,n2,n3p]
xv[0,6]=Bx[n1,n2p,n3p]
xv[0,7]=Bx[n1p,n2p,n3p]
xv[1,6]=By[n1,n2p,n3p]
xv[1,7]=By[n1p,n2p,n3p]
xv[2,6]=Bz[n1,n2p,n3p]
xv[2,7]=Bz[n1p,n2p,n3p]
RETURN
END
