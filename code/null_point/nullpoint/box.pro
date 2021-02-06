pro sub_field,bx,by,bz,obx,oby,obz,i,j,k

obx = dblarr(2,2,2)
oby = dblarr(2,2,2)
obz = dblarr(2,2,2)
;-----------------------------------------
obx[0,0,0] = bx[i,j,k]
obx[0,0,1] = bx[i,j,k+1]
obx[0,1,0] = bx[i,j+1,k]
obx[0,1,1] = bx[i,j+1,k+1]
obx[1,0,0] = bx[i+1,j,k]
obx[1,0,1] = bx[i+1,j,k+1]
obx[1,1,0] = bx[i+1,j+1,k]
obx[1,1,1] = bx[i+1,j+1,k+1]
;-----------------------------------------
oby[0,0,0] = by[i,j,k]
oby[0,0,1] = by[i,j,k+1]
oby[0,1,0] = by[i,j+1,k]
oby[0,1,1] = by[i,j+1,k+1]
oby[1,0,0] = by[i+1,j,k]
oby[1,0,1] = by[i+1,j,k+1]
oby[1,1,0] = by[i+1,j+1,k]
oby[1,1,1] = by[i+1,j+1,k+1]
;-----------------------------------------
obz[0,0,0] = bz[i,j,k]
obz[0,0,1] = bz[i,j,k+1]
obz[0,1,0] = bz[i,j+1,k]
obz[0,1,1] = bz[i,j+1,k+1]
obz[1,0,0] = bz[i+1,j,k]
obz[1,0,1] = bz[i+1,j,k+1]
obz[1,1,0] = bz[i+1,j+1,k]
obz[1,1,1] = bz[i+1,j+1,k+1]
;-----------------------------------------
end


;return boxes, which may contain a 3D null, to the main program null.pro
;the 3D magnetic field is needed for the calculation.
function box,bx,by,bz

ss = size(bz)
box_x = 0
box_y = 0
box_z = 0
n = 0

for i=0,ss[1]-2 do begin
	for j=0,ss[2]-2 do begin
		for k=0,ss[3]-2 do begin

			sub_field,bx,by,bz,obx,oby,obz,i,j,k
			flag = poincare(obx,oby,obz,3)
			if (flag ne 0) then begin
				box_x = [box_x,i]
				box_y = [box_y,j]
				box_z = [box_z,k]
				n = n+1
			endif
		endfor
	endfor
	print,'scaning process '+string(100.*i/ss[1],format='(f5.1)')+'% has been done.'
endfor

box = intarr(n+1,3)

box[*,0] = box_x
box[*,1] = box_y
box[*,2] = box_z

return,box
end
