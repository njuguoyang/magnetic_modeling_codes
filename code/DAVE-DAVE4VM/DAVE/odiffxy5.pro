pro odiffxy5,image,dx,dy
;
;  gradient
   c1=0.12019d0
   c2=0.74038d0
   dx=(-c1*shift(image,-2,0)+c2*shift(image,-1,0)-c2*shift(image,1,0)+c1*shift(image,2,0))
   dy=(-c1*shift(image,0,-2)+c2*shift(image,0,-1)-c2*shift(image,0,1)+c1*shift(image,0,2))
  
;
return 

end
