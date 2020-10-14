function compute_fluxes,mag,vel
;
   common constants,speed_of_light
   if (N_ELEMENTS(speed_of_light) eq 0) then speed_of_light=2.998d10 ; cm/s
;
   dum=vel.U0
   dum[*]=0.d0
;
   vel=add_tag(vel,dum,'VB')
   vel=add_tag(vel,dum,'VPX')
   vel=add_tag(vel,dum,'VPY')
   vel=add_tag(vel,dum,'VPZ')
   vel=add_tag(vel,dum,'EX')
   vel=add_tag(vel,dum,'EY')
   vel=add_tag(vel,dum,'EZ')
   vel=add_tag(vel,dum,'FTVX')
   vel=add_tag(vel,dum,'FTVXX')
   vel=add_tag(vel,dum,'FTVY')
   vel=add_tag(vel,dum,'FTVYY')
   vel=add_tag(vel,dum,'POYNTING')
;   vel=add_tag(vel,dum,'Helicity')
;
;     flux transport vectors
   vel.FTVX=vel.U0*mag.Bz-vel.W0*mag.BX
   vel.FTVY=vel.V0*mag.Bz-vel.W0*mag.BY
;     derivatives
   odiffxy5,vel.FTVX/mag.dx,bzuxx,dum
   odiffxy5,vel.FTVY/mag.dy,dum,bzuyy
   vel.FTVXX=bzuxx
   vel.FTVYY=bzuyy
;
;     Poynting flux
   vel.Poynting=-(vel.FTVX*mag.Bx+vel.FTVY*mag.BY)/(4*!dpi)
;
;     Helicity flux (Demoulin and Berger)
;   vel.Helicity=-2*(vel.FTVX*mag.Ax+vel.FTVY*mag.AY)          ; [G^2 km^2/s]

;
;      parallel to the magnetic field plasma velocity
   vel.VB=mag.qx*vel.u0+mag.qy*vel.v0+mag.qz*vel.w0
;      perpendicular to the magnetic field plasma velocties
   vel.vpx=vel.u0-vel.vb*mag.qx
   vel.vpy=vel.v0-vel.vb*mag.qy
   vel.vpz=vel.w0-vel.vb*mag.qz
;
   vel.Ex=-(vel.V0*mag.bz-vel.W0*mag.by)/speed_of_light
   vel.Ey=-(vel.W0*mag.bx-vel.U0*mag.bz)/speed_of_light
   vel.Ez=-(vel.U0*mag.by-vel.V0*mag.bx)/speed_of_light
;

;
return,vel
;
end
