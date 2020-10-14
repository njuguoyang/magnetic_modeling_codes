function compute_los_velocity,mag,truth,angle
;
;    mag - magnetic data structure from shootout
;  angle - angle between z and los 
;
   sz=size(mag.bz)
;   
   n=[sin(angle),0.d0,cos(angle*!dpi/180.d0)]            ; normal vector
;
;     VX, VY and VZ should be time-centered
   VLOS=n[0]*truth.VX+n[1]*truth.VY+n[2]*truth.VZ
   mag=add_tag(mag,VLOS,'VLOS')
;
   NX=replicate(N[0],sz[1],sz[2])
   NY=replicate(N[1],sz[1],sz[2])
   NZ=replicate(N[2],sz[1],sz[2])
;
   mag=add_tag(mag,NX,'NX')
   mag=add_tag(mag,NY,'NY')
   mag=add_tag(mag,NZ,'NZ')
;
return,mag
;
end
