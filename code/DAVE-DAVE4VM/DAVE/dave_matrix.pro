function dave_matrix,Bz,Bzx,Bzy,Bzt,psf,psfx,psfy,psfxx,psfyy,psfxy,nu,double=double,float=float
;
   sz=size(Bz)
   A=make_array(7,7,sz[1],sz[2],double=double,float=float)
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;    construct matrix elements for LKA algorithm
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;    DAVE elements (depend only on Bz,Bzx,Bzy,Bzt)
;
   GG=Bz*Bz
   G=reform(convol(GG,psf,/center),1,sz[1],sz[2])                       ; 1
;
   GGx=Bz*Bzx
   Gx=reform(convol(GGx,psf,/center),1,sz[1],sz[2])                     ; 2
;
   xGx=reform(convol(GGx,psfx,/center),1,sz[1],sz[2])                   ; 3
;
   yGx=reform(convol(GGx,psfy,/center),1,sz[1],sz[2])                   ; 4
;
   GGy=Bz*Bzy
   Gy=reform(convol(GGy,psf,/center),1,sz[1],sz[2])                     ; 5
;  
   xGy=reform(convol(GGy,psfx,/center),1,sz[1],sz[2])                   ; 6
;
   yGy=reform(convol(GGy,psfy,/center),1,sz[1],sz[2])                   ; 7

   GGt=Bzt*Bz
   Ht=reform(convol(GGt,psf,/center),1,sz[1],sz[2])                     ; 8
;
   GGxx=Bzx*Bzx
   Gxx=reform(convol(GGxx,psf,/center),1,sz[1],sz[2])                   ; 9
;
   GGyy=Bzy*Bzy
   Gyy=reform(convol(GGyy,psf,/center),1,sz[1],sz[2])                   ; 10
;
   GGxy=Bzx*Bzy
   Gxy=reform(convol(GGxy,psf,/center),1,sz[1],sz[2])                   ; 11
;
   GGtx=Bzt*Bzx
   Gtx=reform(convol(GGtx,psf,/center),1,sz[1],sz[2])                   ; 12
;
   GGty=Bzt*Bzy
   Gty=reform(convol(GGty,psf,/center),1,sz[1],sz[2])                   ; 13
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
   xGxx=reform(convol(GGxx,psfx,/center),1,sz[1],sz[2])                 ; 14
;
   xGyy=reform(convol(GGyy,psfx,/center),1,sz[1],sz[2])                 ; 15
;
   xGxy=reform(convol(GGxy,psfx,/center),1,sz[1],sz[2])                 ; 16
;
   xGtx=reform(convol(GGtx,psfx,/center),1,sz[1],sz[2])                 ; 17
;
   xGty=reform(convol(GGty,psfx,/center),1,sz[1],sz[2])                 ; 18
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
   yGxx=reform(convol(GGxx,psfy,/center),1,sz[1],sz[2])                 ; 19
;
   yGyy=reform(convol(GGyy,psfy,/center),1,sz[1],sz[2])                 ; 20
;
   yGxy=reform(convol(GGxy,psfy,/center),1,sz[1],sz[2])                 ; 21
;
   yGtx=reform(convol(GGtx,psfy,/center),1,sz[1],sz[2])                 ; 22
;
   yGty=reform(convol(GGty,psfy,/center),1,sz[1],sz[2])                 ; 23
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
   xxGxx=reform(convol(GGxx,psfxx,/center),1,sz[1],sz[2])               ; 24
;
   xxGxy=reform(convol(GGxy,psfxx,/center),1,sz[1],sz[2])               ; 25
;
   xxGyy=reform(convol(GGyy,psfxx,/center),1,sz[1],sz[2])               ; 26
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
   xyGxx=reform(convol(GGxx,psfxy,/center),1,sz[1],sz[2])               ; 27
;
   xyGyy=reform(convol(GGyy,psfxy,/center),1,sz[1],sz[2])               ; 28
;
   xyGxy=reform(convol(GGxy,psfxy,/center),1,sz[1],sz[2])               ; 29
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
   yyGxx=reform(convol(GGxx,psfyy,/center),1,sz[1],sz[2])               ; 30
;
   yyGxy=reform(convol(GGxy,psfyy,/center),1,sz[1],sz[2])               ; 31
;
   yyGyy=reform(convol(GGyy,psfyy,/center),1,sz[1],sz[2])               ; 32
;
   GGtt=Bzt*Bzt
   Gtt=reform(convol(GGtt,psf,/center),1,sz[1],sz[2])                   ; 33
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   
;   
   A[*]=reform([Gxx,Gxy,nu*Gx+xGxx,nu*Gx+yGxy,yGxx,xGxy,Gtx,$
          Gxy,Gyy,nu*Gy+xGxy,nu*Gy+yGyy,yGxy,xGyy,Gty,$
          nu*Gx+xGxx,nu*Gy+xGxy,nu*(G*nu+2*xGx)+xxGxx,nu*(G*nu+xGx+yGy)+xyGxy,nu*yGx+xyGxx,nu*xGy+xxGxy,nu*Ht+xGtx,$   
          nu*Gx+yGxy,nu*Gy+yGyy,nu*(G*nu+xGx+yGy)+xyGxy,nu*(G*nu+2*yGy)+yyGyy,nu*yGx+yyGxy,nu*xGy+xyGyy,nu*Ht+yGty,$
          yGxx,yGxy,nu*yGx+xyGxx,nu*yGx+yyGxy,yyGxx,xyGxy,yGtx,$
          xGxy,xGyy,nu*xGy+xxGxy,nu*xGy+xyGyy,xyGxy,xxGyy,xGty,$
          Gtx,Gty,nu*Ht+xGtx,nu*Ht+yGty,yGtx,xGty,Gtt],7,7,sz[1],sz[2])
;
   return,A

end
