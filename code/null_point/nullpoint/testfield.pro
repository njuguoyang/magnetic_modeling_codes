pro testfield
;give a test field for testing the null scaning program
;we use the magnetic charges model, 4 charges with their magnitude of unit: 1(0,0,0),-1(0,1,0),-1(cos(7pi/6),sin(7pi/6),0),-1(cos(11pi/6),sin(11pi/6),0)
;domain dimension is (2*2*2)
;per pixel is 0.02
pixel = 0.02
bx = dblarr(200,200,200)
by = dblarr(200,200,200)
bz = dblarr(200,200,200)
origin=[99.5,99.5,0]
c1 = [0.,0.,-0.05];positive
c2 = [0.,1.,-0.05];negative
c3 = [cos(7*!dpi/6),sin(7*!dpi/6),-0.05];negative
c4 = [cos(11*!dpi/6),sin(11*!dpi/6),-0.05];negative
for i=0,199 do begin
	for j=0,199 do begin
		for k=0,199 do begin
			r = (1.*[i,j,k]-origin)*pixel
			r1 = r - c1
			r2 = r - c2
			r3 = r - c3
			r4 = r - c4
			b1 = 1.*r1/(sqrt(r1[0]^2 + r1[1]^2 + r1[2]^2)^3)
			b2 = -1.*r2/(sqrt(r2[0]^2 + r2[1]^2 + r2[2]^2)^3)
			b3 = -1.*r3/(sqrt(r3[0]^2 + r3[1]^2 + r3[2]^2)^3)
			b4 = -1.*r4/(sqrt(r4[0]^2 + r4[1]^2 + r4[2]^2)^3)
			bx[i,j,k] = b1[0] + b2[0] + b3[0] + b4[0]
			by[i,j,k] = b1[1] + b2[1] + b3[1] + b4[1]
			bz[i,j,k] = b1[2] + b2[2] + b3[2] + b4[2]
		endfor
	endfor
endfor

save,bx,by,bz,filename='testfield.sav'

print,'writing the vector field !'

ss=size(bz)

openw,1,'testfield.vtk'

printf,1,'# vtk DataFile Version 2.0'
printf,1,'Volume example'
printf,1,'ASCII'
printf,1,'DATASET STRUCTURED_POINTS'
printf,1,strjoin([strcompress('DIMENSIONS',/remove_all),strcompress(string(ss[1]),/remove_all),strcompress(string(ss[2]),/remove_all),strcompress(string(ss[3]),/remove_all)],' ')
printf,1,strjoin([strcompress('ASPECT_RATIO',/remove_all),strcompress(string(1.0),/remove_all),strcompress(string(1.0),/remove_all),strcompress(string(1.0),/remove_all)],' ')
printf,1,strjoin([strcompress('ORIGIN',/remove_all),strcompress(string(0),/remove_all),strcompress(string(0),/remove_all),strcompress(string(0),/remove_all)],' ')
printf,1,strjoin([strcompress('POINT_DATA',/remove_all),strcompress(string(ulong(ss[1]*ss[2]*ss[3]*1.0)),/remove_all)],' ')
printf,1,'VECTORS testfield FLOAT'

for k=0,ss[3]-1 do begin
    for j=0,ss[2]-1 do begin
        for i=0,ss[1]-1 do begin
        printf,1,bx[i,j,k],by[i,j,k],bz[i,j,k],format='(3E18.8)'
        endfor
   endfor
endfor
close,1
print,'io of data is completed !!!'

end
