;+
; NAME:
;   find_null
;
; PURPOSE: calculate the position, position error, null magnitude, 
;		   eigenvectors, eigenvalues, and the matrix of the null-point
;     	   by the Newton-Raphson scheme
;
; CALLING SEQUENCE:
;	find_null,bx,by,bz,'./null_para.sav','./cube_index.sav'
;
; INPUTS:
;	cubeindex : the box where the null locate with the poincare principle
;	bx,by,bz  : the magnetic field to look for a null
;
; PARAMETERS:
;         eps: threshold
;         n0 : the number of the step
;
; OUTPUTS:
;	save_name : contain the position, position error, null magnitude, 
;		        eigenvectors, eigenvalues, and the matrix of the null-point
;         
; MODIFICATION HISTORY:
; written by K. Yang 2014 Dec 3 NJU
; Modified by Z. Zhong (NJU, 2017 June 1) 
;-

pro find_null,bx,by,bz,save_name,cubeindex

;msp = machar()
;eps = msp.eps
eps = 1e-5
n0=99

xindex = dblarr(n0+1)
yindex = dblarr(n0+1)
zindex = dblarr(n0+1)
note1 = dblarr(n0+1)

restore,cubeindex
ss = size(cube_index)
null_num=0
if (ss[1] eq 1) then begin
	
	print,'*** There is no null inside the computational domain!!!'

endif else begin

	print,'There are '+strcompress(string(ss[1]-1,format='(I)'),/remove_all)+' null-points inside the computational domain !'
	null_posi = dblarr(ss[1]-1,3)
	null_posi_error = dblarr(ss[1]-1,3)
	null_matr = dblarr(ss[1]-1,3,3)
	posi_null = dblarr(3)
	note = dblarr(ss[1]-1)
	matr_null = dblarr(3,3)
	vector = complexarr(3,3)
	    
	for i=0,ss[1]-2 do begin
;		if ((i mod 5) eq 0) then begin
		print,'***** Calculating the '+strcompress(string(i),/remove_all)+' possible null !!!'
		print,'***** Newton iteration begin !'
;		endif
		
		for j=0,n0 do begin
			ran = randomU(seed,3)
      		posi_null = null_newton_raphson(cube_index[i+1,*]+ran,bx,by,bz,eps,matr_null,note0)
			note1[j] = note0
			xindex[j] = posi_null[0]
			yindex[j] = posi_null[1]
			zindex[j] = posi_null[2]
		endfor
		if (total(note1) eq 0 ) then begin
;		  print,'this is'+string(i)+'null point'
		  print, 'no null'
            x_mean = median(xindex)
            y_mean = median(yindex)
            z_mean = median(zindex)
            x_error = stddev(xindex)
            y_error = stddev(yindex)
            z_error = stddev(zindex)		  
		endif else begin 
            null_pos = where(note1 eq 1)
            print,n_elements(null_pos)
            print,xindex[null_pos],yindex[null_pos],zindex[null_pos]      
            x_mean = median(xindex[null_pos])
            y_mean = median(yindex[null_pos])
            z_mean = median(zindex[null_pos])
            x_error = stddev(xindex[null_pos])
            y_error = stddev(yindex[null_pos])
            z_error = stddev(zindex[null_pos])
    	endelse   
     	    ;mean position and error of the coordinate for the null
     	    null_posi[i,*] = [x_mean,y_mean,z_mean]
     	    null_posi_error[i,*] = [x_error,y_error,z_error]
     	    ;matrix value at the null
     	    matr = matrix_interp(bx,by,bz,x_mean,y_mean,z_mean)
     	    null_matr[i,*,*] = matr_null
     	    ;value of |B| at the null
     	    fieldb = field_interp(bx,by,bz,x_mean,y_mean,z_mean)
      	    note[i] = sqrt(fieldb[0]^2 + fieldb[1]^2 + fieldb[2]^2) 

	endfor

	null_num = where(note le eps)
	null_size = size(null_num)

	if (null_size[0] eq 0) then begin
		print,'*** There is no null inside the computational domain!!!'
	endif else begin 
	
    null_eigen_value = complexarr(null_size[1],3)
    null_eigen_vector = complexarr(null_size[1],3,3)
	null_position = dblarr(null_size[1],3)
	null_position_error = dblarr(null_size[1],3)
	null_matrix = dblarr(null_size[1],3,3)
	null_magnitude = dblarr(null_size[1])

	print,'**** begin to calculate the parameters of those null-points'
	for i=0,null_size[1]-1 do begin
		print,'**** No.'+strcompress(string(i),/remove_all)+' null point.'
		matrix_null=reform(null_matr(null_num[i],*,*))

		value = HQR(ELMHES(matrix_null),/double)

		if (value[0] eq 0 or value[1] eq 0 or value[2] eq 0) then begin
			vector = [[0.,0.,0.],[0,0,0],[0,0,0]]
		endif else begin
			vector = eigenvec(matrix_null,value)
		endelse
		null_eigen_value[i,*] = value
		null_eigen_vector[i,*,*] = vector
		null_position[i,*] = null_posi[null_num[i],*]
		null_position_error[i,*] = null_posi_error[null_num[i],*]
		null_matrix[i,*,*] = null_matr[null_num[i],*,*]
		null_magnitude[i] = note[null_num[i]]

	endfor
		save,null_position,null_position_error,null_magnitude,null_matrix,null_eigen_value,null_eigen_vector,filename=save_name
	endelse
endelse
end
