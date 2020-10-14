function dave4vm_matrix,Bx,Bxx,Bxy,By,Byx,Byy,Bz,Bzx,Bzy,Bzt,psf,psfx,psfy,psfxx,psfyy,psfxy,double=double,float=float
;
   sz=size(Bz)
   A=make_array(10,10,sz[1],sz[2],double=double,float=float)
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;    construct matrix elements for LKA algorithm
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;    DAVE elements (depend only on Bz,Bzx,Bzy,Bzt)
;
   G=reform(convol(Bz*Bz,psf,/center),1,sz[1],sz[2])                    ; 1
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
;  extra vector magnetogram terms
   BxBx=reform(convol(Bx*Bx,psf,/center),1,sz[1],sz[2]) 
   ByBy=reform(convol(By*By,psf,/center),1,sz[1],sz[2]) 
   BxBy=reform(convol(Bx*By,psf,/center),1,sz[1],sz[2]) 
   BzBx=reform(convol(Bz*Bx,psf,/center),1,sz[1],sz[2]) 
   BzBy=reform(convol(Bz*By,psf,/center),1,sz[1],sz[2])
; 
   BxBxx=reform(convol(Bx*Bxx,psf,/center),1,sz[1],sz[2]) 
   BxByy=reform(convol(Bx*Byy,psf,/center),1,sz[1],sz[2]) 
   BxxBxx=reform(convol(Bxx*Bxx,psf,/center),1,sz[1],sz[2]) 
   ByyByy=reform(convol(Byy*Byy,psf,/center),1,sz[1],sz[2]) 
   BxxByy=reform(convol(Bxx*Byy,psf,/center),1,sz[1],sz[2]) 
   ByBxx=reform(convol(By*Bxx,psf,/center),1,sz[1],sz[2])  
   ByByy=reform(convol(By*Byy,psf,/center),1,sz[1],sz[2]) 
;
   BzBxx=reform(convol(Bz*Bxx,psf,/center),1,sz[1],sz[2]) 
   BzByy=reform(convol(Bz*Byy,psf,/center),1,sz[1],sz[2])
;
   BztBxx=reform(convol(Bzt*Bxx,psf,/center),1,sz[1],sz[2]) 
   BztByy=reform(convol(Bzt*Byy,psf,/center),1,sz[1],sz[2])
; 
   BzxBx=reform(convol(Bzx*Bx,psf,/center),1,sz[1],sz[2]) 
   BzxBy=reform(convol(Bzx*By,psf,/center),1,sz[1],sz[2]) 
   BzxBxx=reform(convol(Bzx*Bxx,psf,/center),1,sz[1],sz[2]) 
   BzxByy=reform(convol(Bzx*Byy,psf,/center),1,sz[1],sz[2]) 
;
   BzyBx=reform(convol(Bzy*Bx,psf,/center),1,sz[1],sz[2]) 
   BzyBy=reform(convol(Bzy*By,psf,/center),1,sz[1],sz[2]) 
   BzyBxx=reform(convol(Bzy*Bxx,psf,/center),1,sz[1],sz[2]) 
   BzyByy=reform(convol(Bzy*Byy,psf,/center),1,sz[1],sz[2]) 
;
   BztBx=reform(convol(Bzt*Bx,psf,/center),1,sz[1],sz[2]) 
   BztBy=reform(convol(Bzt*By,psf,/center),1,sz[1],sz[2]) 
;
   xBzxBx=reform(convol(Bzx*Bx,psfx,/center),1,sz[1],sz[2]) 
   xBzxBy=reform(convol(Bzx*By,psfx,/center),1,sz[1],sz[2]) 
   xBzyBx=reform(convol(Bzy*Bx,psfx,/center),1,sz[1],sz[2]) 
   xBzyBy=reform(convol(Bzy*By,psfx,/center),1,sz[1],sz[2]) 
;
   yBzyBx=reform(convol(Bzy*Bx,psfy,/center),1,sz[1],sz[2]) 
   yBzyBy=reform(convol(Bzy*By,psfy,/center),1,sz[1],sz[2]) 
;
   yBzxBx=reform(convol(Bzx*Bx,psfy,/center),1,sz[1],sz[2]) 
   yBzxBy=reform(convol(Bzx*By,psfy,/center),1,sz[1],sz[2]) 
;
   yBxBxx=reform(convol(Bx*Bxx,psfy,/center),1,sz[1],sz[2]) 
   yBxByy=reform(convol(Bx*Byy,psfy,/center),1,sz[1],sz[2]) 
;
   yByBxx=reform(convol(By*Bxx,psfy,/center),1,sz[1],sz[2]) 
   yByByy=reform(convol(By*Byy,psfy,/center),1,sz[1],sz[2]) 
;
   xByBxx=reform(convol(By*Bxx,psfx,/center),1,sz[1],sz[2]) 
   xByByy=reform(convol(By*Byy,psfx,/center),1,sz[1],sz[2]) 
;
   xBzxBxx=reform(convol(Bzx*Bxx,psfx,/center),1,sz[1],sz[2]) 
   xBzxByy=reform(convol(Bzx*Byy,psfx,/center),1,sz[1],sz[2]) 
;
   yBzxBxx=reform(convol(Bzx*Bxx,psfy,/center),1,sz[1],sz[2]) 
   yBzxByy=reform(convol(Bzx*Byy,psfy,/center),1,sz[1],sz[2]) 
;
   xBxxBxx=reform(convol(Bxx*Bxx,psfx,/center),1,sz[1],sz[2]) 
   xBxxByy=reform(convol(Bxx*Byy,psfx,/center),1,sz[1],sz[2]) 
   xByyByy=reform(convol(Byy*Byy,psfx,/center),1,sz[1],sz[2]) 
;
   yBxxBxx=reform(convol(Bxx*Bxx,psfy,/center),1,sz[1],sz[2]) 
   yBxxByy=reform(convol(Bxx*Byy,psfy,/center),1,sz[1],sz[2]) 
   yByyByy=reform(convol(Byy*Byy,psfy,/center),1,sz[1],sz[2]) 
;
   xBxBxx=reform(convol(Bx*Bxx,psfx,/center),1,sz[1],sz[2]) 
   xBxByy=reform(convol(Bx*Byy,psfx,/center),1,sz[1],sz[2]) 
;
   xBzBxx=reform(convol(Bz*Bxx,psfx,/center),1,sz[1],sz[2]) 
   xBzByy=reform(convol(Bz*Byy,psfx,/center),1,sz[1],sz[2]) 
;
   xBztBxx=reform(convol(Bzt*Bxx,psfx,/center),1,sz[1],sz[2]) 
   xBztByy=reform(convol(Bzt*Byy,psfx,/center),1,sz[1],sz[2]) 
;
   yBztBxx=reform(convol(Bzt*Bxx,psfy,/center),1,sz[1],sz[2]) 
   yBztByy=reform(convol(Bzt*Byy,psfy,/center),1,sz[1],sz[2]) 
;
   xyBxxBxx=reform(convol(Bxx*Bxx,psfxy,/center),1,sz[1],sz[2]) 
   xyBxxByy=reform(convol(Bxx*Byy,psfxy,/center),1,sz[1],sz[2]) 
   xyByyByy=reform(convol(Byy*Byy,psfxy,/center),1,sz[1],sz[2]) 
;
   xyBzxBxx=reform(convol(Bzx*Bxx,psfxy,/center),1,sz[1],sz[2]) 
   xyBzxByy=reform(convol(Bzx*Byy,psfxy,/center),1,sz[1],sz[2]) 
   xyBzyBxx=reform(convol(Bzy*Bxx,psfxy,/center),1,sz[1],sz[2]) 
   xyBzyByy=reform(convol(Bzy*Byy,psfxy,/center),1,sz[1],sz[2]) 
;
   yBzBxx=reform(convol(Bz*Bxx,psfy,/center),1,sz[1],sz[2]) 
   yBzByy=reform(convol(Bz*Byy,psfy,/center),1,sz[1],sz[2]) 
;
   xBzyBxx=reform(convol(Bzy*Bxx,psfx,/center),1,sz[1],sz[2]) 
   xBzyByy=reform(convol(Bzy*Byy,psfx,/center),1,sz[1],sz[2]) 
   yBzyBxx=reform(convol(Bzy*Bxx,psfy,/center),1,sz[1],sz[2]) 
   yBzyByy=reform(convol(Bzy*Byy,psfy,/center),1,sz[1],sz[2]) 
;
   xxBxxBxx=reform(convol(Bxx*Bxx,psfxx,/center),1,sz[1],sz[2]) 
   xxBxxByy=reform(convol(Bxx*Byy,psfxx,/center),1,sz[1],sz[2]) 
   xxByyByy=reform(convol(Byy*Byy,psfxx,/center),1,sz[1],sz[2]) 
;
   xxBzxBxx=reform(convol(Bzx*Bxx,psfxx,/center),1,sz[1],sz[2]) 
   xxBzyBxx=reform(convol(Bzy*Bxx,psfxx,/center),1,sz[1],sz[2]) 
   xxBzxByy=reform(convol(Bzx*Byy,psfxx,/center),1,sz[1],sz[2]) 
   xxBzyByy=reform(convol(Bzy*Byy,psfxx,/center),1,sz[1],sz[2]) 
;
   yyBxxBxx=reform(convol(Bxx*Bxx,psfyy,/center),1,sz[1],sz[2]) 
   yyBxxByy=reform(convol(Bxx*Byy,psfyy,/center),1,sz[1],sz[2]) 
   yyByyByy=reform(convol(Byy*Byy,psfyy,/center),1,sz[1],sz[2]) 
;
   yyBzyBxx=reform(convol(Bzy*Bxx,psfyy,/center),1,sz[1],sz[2]) 
   yyBzyByy=reform(convol(Bzy*Byy,psfyy,/center),1,sz[1],sz[2]) 
;
   yyBzxBxx=reform(convol(Bzx*Bxx,psfyy,/center),1,sz[1],sz[2]) 
   yyBzxByy=reform(convol(Bzx*Byy,psfyy,/center),1,sz[1],sz[2]) 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   @vm_matrix
;   
   return,reform(A,10,10,sz[1],sz[2])

end
