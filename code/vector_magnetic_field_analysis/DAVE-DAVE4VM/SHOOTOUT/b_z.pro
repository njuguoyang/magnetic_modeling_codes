function B_Z,X,Z
   
     common loop,Z0,A
;
     R=abs(X) ; cylindrical "R"
     KC=Sqrt(4*A*R/((R+A)^2+(Z-Z0)^2))
     B_Z=1.d0/sqrt((R^2+A^2)+(Z-Z0)^2)*$
         (ELLIPTIC_K(KC)-(R^2-A^2+(Z-Z0)^2)/((R-A)^2+(Z-Z0)^2)*ELLIPTIC_E(KC))
return,B_Z

end

