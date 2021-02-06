pro diffxy,image,dx,dy
;
   sz=size(image)
;
;  gradient dx
    dx = (shift(image,-1,0) - shift(image,1,0))/2.
    dx[0,*] = (-3.0*image[0,*] + 4.0*image[1,*] - image[2,*])/2.
    dx[sz[1]-1,*] = (3.*image[sz[1]-1,*] - 4.*image[sz[1]-2,*] + image[sz[1]-3,*])/2.
;
;  gradient dy
    dy = (shift(image,0,-1) - shift(image,0,1))/2.
    dy[*,0] = (-3.0*image[*,0] + 4.0*image[*,1] - image[*,2])/2.
    dy[*,sz[2]-1] = (3.*image[*,sz[2]-1] - 4.*image[*,sz[2]-2] + image[*,sz[2]-3])/2.
;    stop
;
return 

end

