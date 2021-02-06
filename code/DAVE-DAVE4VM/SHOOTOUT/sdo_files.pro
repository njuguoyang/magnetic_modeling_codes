pro sdo_files

   DIR='SDO/'
   T_OBS00 = '2004-09-01T00:00:00.000
   T00=anytim(T_OBS00)

;
   verbose=1L
   base='hmi_lev1_hr_'
   file='./old_shootout/data_for_analysis.sav'
   restore,file,verbose=verbose,_extra=extra
   sz=size(BZ_START)
   TAU0=(t_stop[0]+t_start[0])/2
   writefits,'dum.fts',reform(BZ_START[*,*,0],sz[1],sz[2]),hd
   hd0=hd
   hd1=hd
   sxaddpar,hd0,'T_OBS',T_OBS00
   sxaddpar,hd1,'T_OBS',T_OBS00
;

   for i=0,sz[3]-1L do begin
       tags=['BX','BY','BZ']
;
       T0=t_start[i]-TAU0+T00
       T_OBS0=anytim(T0,/ccsds)
       TIME0=STRJOIN(STRSPLIT(strmid(T_OBS0,0,19), '._:-T',/EXTRACT))
       sxaddpar,hd0,'T_OBS',T_OBS0
;
       filename=dir+base+tags[0]+'_'+TIME0+'.fts'
;       print,filename & stop
       writefits,filename,reform(bx_start[*,*,i],sz[1],sz[2]),hd0
       filename=dir+base+tags[1]+'_'+TIME0+'.fts'
       writefits,filename,reform(by_start[*,*,i],sz[1],sz[2]),hd0
       filename=dir+base+tags[2]+'_'+TIME0+'.fts'
       writefits,filename,reform(bz_start[*,*,i],sz[1],sz[2]),hd0
;
       T1=t_stop[i]-TAU0+T00
       T_OBS1=anytim(T1,/ccsds)
       TIME1=STRJOIN(STRSPLIT(strmid(T_OBS1,0,19), '._:-T',/EXTRACT))
       sxaddpar,hd1,'T_OBS',T_OBS1

       DELTA_T=anytim(T_OBS1)-anytim(T_OBS0)
       print,i,' ',T_OBS1,'  ',T_OBS0,delta_t-(t_stop[i]-t_start[i])

       filename=dir+base+tags[0]+'_'+TIME1+'.fts'
       writefits,filename,reform(bx_stop[*,*,i],sz[1],sz[2]),hd1
       filename=dir+base+tags[1]+'_'+TIME1+'.fts'
       writefits,filename,reform(by_stop[*,*,i],sz[1],sz[2]),hd1
       filename=dir+base+tags[2]+'_'+TIME1+'.fts'
       writefits,filename,reform(bz_stop[*,*,i],sz[1],sz[2]),hd1
    endfor
;
   stop
;
end

pro check_files

     DIR='SDO/'
     initial=['hmi_lev1_hr_BX_20040831235754.fts','hmi_lev1_hr_BY_20040831235754.fts','hmi_lev1_hr_BZ_20040831235754.fts']
     final=['hmi_lev1_hr_BX_20040901000205.fts','hmi_lev1_hr_BY_20040901000205.fts','hmi_lev1_hr_BZ_20040901000205.fts']

     file='./old_shootout/data_for_analysis.sav'
     restore,file,verbose=verbose,_extra=extra
;
     BXI=mrdfits(DIR+initial[0])
     BXF=mrdfits(DIR+final[0])
     BYI=mrdfits(DIR+initial[1])
     BYF=mrdfits(DIR+final[1])
     BZI=mrdfits(DIR+initial[2])
     BZF=mrdfits(DIR+final[2])
     print,total((BXI-BX_START[*,*,0])^2,/double)
     print,total((BYI-BY_START[*,*,0])^2,/double)
     print,total((BZI-BZ_START[*,*,0])^2,/double)
     print,total((BXF-BX_STOP[*,*,0])^2,/double)
     print,total((BYF-BY_STOP[*,*,0])^2,/double)
     print,total((BZF-BZ_STOP[*,*,0])^2,/double)
     stop

end
