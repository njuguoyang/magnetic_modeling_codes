;+
;PURPOSE: To align the 9 selected magnetic field maps consecutively by feature identification
;         method, i.e., identifying common features of two images to estimate the offsets
;         between them.
;-


DIR='GUO2/'
file=find_file('../02ambiguity/04ambiguity_removed/mapbxyz*')
index=[11, 25, 27, 44, 63, 76, 92, 106, 118]
file=file[index]
;
nn=n_elements(file)
;nn=1    ;for test
for i=0,nn-2 do begin
  if (i eq 0) then begin
    restore,file[i]
    save,smapbz,smapbx,smapby,filename=DIR+'mapbxyz'+strtrim(string(i,format='(I03)'),2)+'.sav'
  endif else begin
    restore,DIR+'mapbxyz'+strtrim(string(i,format='(I03)'),2)+'.sav'
  endelse
  data_ref=smapbz.data
  dxmtr=smapbz.dx
  dymtr=smapbz.dy
  restore,file[i+1]
  dxmdi=smapbz.dx
  dymdi=smapbz.dy
  setpts,p,smapbz.data,data_ref
  xmeanmdi=mean(p[0,0,*])
  ymeanmdi=mean(p[1,0,*])
  xmeanmtr=mean(p[0,1,*])
  ymeanmtr=mean(p[1,1,*])
  ssmdi=size(data_ref)
  xcmdi=ssmdi[1]/2.0
  ycmdi=ssmdi[2]/2.0
  xcmtr_ref2mdi=xmeanmtr-dxmdi/dxmtr*(xmeanmdi-xcmdi)
  ycmtr_ref2mdi=ymeanmtr-dymdi/dymtr*(ymeanmdi-ycmdi)
  ssmtr=size(smapbz.data)
  xcmtr_index=ssmtr[1]/2.0
  ycmtr_index=ssmtr[2]/2.0
  xcshift=xcmtr_index-xcmtr_ref2mdi
  ycshift=ycmtr_index-ycmtr_ref2mdi
  xcshift=round(xcshift)
  ycshift=round(ycshift)
  print,'x- and y-shifts (pixel):',xcshift,ycshift
  data=shift(smapbz.data,xcshift,ycshift)
  smapbz.data=data
  data=shift(smapbx.data,xcshift,ycshift)
  smapbx.data=data
  data=shift(smapby.data,xcshift,ycshift)
  smapby.data=data
  save,smapbz,smapbx,smapby,filename=DIR+'mapbxyz'+strtrim(string(i+1,format='(I03)'),2)+'.sav'
endfor

end



