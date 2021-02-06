;for 2D
;return the angle of the horzental field.
function angle,rxy1,rxy2

angle0 = acos(transpose(rxy1)#rxy2); dot product of two direction vector
direction0 = sgn(rxy2[0]*rxy1[1] - rxy1[0]*rxy2[1]);cross product of r2 \times r1
angle0 = angle0*direction
return,angle0
end

function cross,rxy1,rxy2
;calculate the cross product for two vectors 1 and 2
rxy3 = dblarr(3)
rxy3[0] = rxy1[1]*rxy2[2] - rxy1[2]*rxy2[1]
rxy3[1] = rxy1[2]*rxy2[0] - rxy1[0]*rxy2[2]
rxy3[2] = rxy1[0]*rxy2[1] - rxy1[1]*rxy2[0]
return,rxy3 
end

;for 3D
function area,rxy1,rxy2,rxy3
;calculate the angle between plane 1-2 and 1-3
r1=(transpose(rxy1)#rxy2)*rxy1
r2=(transpose(rxy1)#rxy3)*rxy1
h1 = rxy2 - r1
h2 = rxy3 - r2
area0 = acos((transpose(h1)#h2)/(sqrt(transpose(h1)#h1) * sqrt(transpose(h2)#h2)))
return,area0
end

function orientation,rxy1,rxy2,rxy3
;calculate the orientation for the new arraged sequence of the vector.
product = TRANSPOSE(cross(rxy2-rxy1,rxy3-rxy2))#rxy1
orientation = sgn(product)
;help,product,/str
return,orientation
end


pro sequence_v,bx,by,bz,v1,v2,v3
v1 = dblarr(3,12)
v2 = dblarr(3,12)
v3 = dblarr(3,12)

AA = [bx[0,0,0],by[0,0,0],bz[0,0,0]]
BB = [bx[1,0,0],by[1,0,0],bz[1,0,0]]
CC = [bx[1,1,0],by[1,1,0],bz[1,1,0]]
DD = [bx[0,1,0],by[0,1,0],bz[0,1,0]]
EE = [bx[0,0,1],by[0,0,1],bz[0,0,1]]
FF = [bx[1,0,1],by[1,0,1],bz[1,0,1]]
GG = [bx[1,1,1],by[1,1,1],bz[1,1,1]]
HH = [bx[0,1,1],by[0,1,1],bz[0,1,1]]
;-----------------------------------
v1[*,0] = AA
v2[*,0] = DD
v3[*,0] = BB

v1[*,1] = DD
v2[*,1] = CC
v3[*,1] = BB
;-----------------------------------
v1[*,2] = AA
v2[*,2] = BB
v3[*,2] = EE

v1[*,3] = BB
v2[*,3] = FF
v3[*,3] = EE
;-----------------------------------
v1[*,4] = AA
v2[*,4] = EE
v3[*,4] = DD

v1[*,5] = DD
v2[*,5] = EE
v3[*,5] = HH
;-----------------------------------
v1[*,6] = CC
v2[*,6] = GG
v3[*,6] = FF

v1[*,7] = BB
v2[*,7] = CC
v3[*,7] = FF
;-----------------------------------
v1[*,8] = GG
v2[*,8] = HH
v3[*,8] = FF

v1[*,9] = EE
v2[*,9] = FF
v3[*,9] = HH
;-----------------------------------
v1[*,10] = GG
v2[*,10] = CC
v3[*,10] = HH

v1[*,11] = CC
v2[*,11] = DD
v3[*,11] = HH
;-----------------------------------
end

;Return the poincare index for the identification of singularity in 2D and 3D vector field.
;Details of the theory and the algorithm can be found in the thesis of Dr. Zhao, Hui.
function poincare,bx,by,bz,dim
index = 0
bn = sqrt(bx^2 + by^2 + bz^2)
bx = bx/bn
by = by/bn
bz = bz/bn
msp = machar()
eps = msp.eps
if (dim ne 2 and dim ne 3) then begin
	print,'Dimension of data should be equal to 2 or 3 !'
	index = -999
endif

if (dim eq 2) then begin
	index = (angle([bx[0,0],by[0,0]],[bx[0,1],by[0,1]])) + (angle([bx[0,1],by[0,1]],[bx[1,1],by[1,1]])) + (angle([bx[1,1],by[1,1]],[bx[1,0],by[1,0]])) + (angle([bx[1,0],by[1,0]],[bx[0,0],by[0,0]]))
	index = index/(2*!dpi)
	if(abs(index) lt eps) then begin
		index = 0
	endif
endif

if(dim eq 3) then begin
area3 = 0.
sequence_v,bx,by,bz,v1,v2,v3
orientation0 = 1
	for i=0,11 do begin
		area0 = area(reform(v1[*,i]),reform(v2[*,i]),reform(v3[*,i]))
		area1 = area(reform(v2[*,i]),reform(v3[*,i]),reform(v1[*,i]))
		area2 = area(reform(v3[*,i]),reform(v1[*,i]),reform(v2[*,i]))
		orientation0 = orientation(reform(v1[*,i]),reform(v2[*,i]),reform(v3[*,i]))
		area3 = area3 + (area0 + area1 + area2 - !dpi)*orientation0
	endfor
area3 = area3/(4.*!dpi)
index = area3
	if(abs(index) lt eps) then begin
	index = 0
	endif
endif

return,index
end
