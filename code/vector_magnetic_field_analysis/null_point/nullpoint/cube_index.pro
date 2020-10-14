pro cube_index,bx,by,bz

xindex = dblarr(100)
yindex = dblarr(100)
zindex = dblarr(100)
; set the minmum value 
print,'**** calculating the boxes begin !!!'
;looking for the box where the null locate with the poincare principle.
cube_index = box(bx,by,bz)
print,'**** the boxes are calculated !!!'
save,cube_index,filename='cube_index.sav'
end
