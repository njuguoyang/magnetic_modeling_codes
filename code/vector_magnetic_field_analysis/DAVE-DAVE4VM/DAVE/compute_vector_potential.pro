pro compute_vector_potential,mudpack,source,phi_structure,backward=backward,stop=stop
;
;   mudpack.loud=0
   sz=size(source)
   cxtxt=replicate(1.d0,sz[1],sz[2])
   cytyt=replicate(1.d0,sz[1],sz[2])
   cxtyt=replicate(0.d0,sz[1],sz[2])
   cxt=replicate(0.d0,sz[1],sz[2])
   cyt=replicate(0.d0,sz[1],sz[2])
   ce=replicate(0.d0,sz[1],sz[2])
;
   data=reform([[source],[cxtxt],[cytyt],[cxtyt],[cxt],[cyt],[ce]],sz[1],sz[2],7)
;
   mudpack_engine,mudpack,data,phi
;
   derives2,mudpack.dx,mudpack.dy,phi,phi_xx,phi_yy,phi_xy,phi_x,phi_y,$
            backward=backward
;   Ax=-phi_y   
;   Ay=phi_x   
   res=phi_xx+phi_yy-source
   if keyword_set(stop) then stop
   phi_structure={scalar:phi,x:phi_x,y:phi_y,xx:phi_xx,yy:phi_yy,xy:phi_xy,res:res}
;
   return
   
end
