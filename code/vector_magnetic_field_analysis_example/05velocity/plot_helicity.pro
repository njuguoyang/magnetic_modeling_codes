;pro plot_helicity

DIR='./GUO/'
helicity=fltarr(118)
time=fltarr(118)
for i=0,117 do begin
  if (i eq 0) then begin
    time[i]=mag.dt/2.0
  endif else begin
    time[i]=mag.dt/2.0+t_half+time[i-1]
  endelse
  t_half=mag.dt/2.0
  restore,DIR+'dave4vm'+strtrim(string(i,format='(I03)'),2)+'.sav'
  sz=size(vel.u0)
  x1=vel.window_size[0]
  x2=sz[1]-1-vel.window_size[0]
  y1=vel.window_size[1]
  y2=sz[2]-1-vel.window_size[1]
  helicity[i]=total(vel.helicity[x1:x2,y1:y2])
endfor
file=find_file('../02ambiguity/alignment_02/mapbxyz001*')
restore,file
utplot,time,helicity,smapbz.time,background=255,color=0
end
