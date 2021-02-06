pro diffxy5,image,dx,dy
;
;  gradient
   
   dx=(-shift(image,-2,0)+8*shift(image,-1,0)-8*shift(image,1,0)+shift(image,2,0))/12
   dy=(-shift(image,0,-2)+8*shift(image,0,-1)-8*shift(image,0,1)+shift(image,0,2))/12
;
return 

end
