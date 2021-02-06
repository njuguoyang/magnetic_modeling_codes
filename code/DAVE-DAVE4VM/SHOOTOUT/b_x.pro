function B_X,X,Z
   
     common loop,Z0,A
;
     R=abs(X) ; cylindrical "R"
     KC=Sqrt(4*A*R/((R+A)^2+(Z-Z0)^2))

     B_X=(z-z0)/X/sqrt((R^2+A^2)+(Z-Z0)^2)*$
         (-ELLIPTIC_K(KC)+(R^2+A^2+(Z-Z0)^2)/((R-A)^2+(Z-Z0)^2)*ELLIPTIC_E(KC))

     index=where(X eq 0,NZ)
     if (NZ gt 0) then B_X[index]=0.d0

return,B_X

end
