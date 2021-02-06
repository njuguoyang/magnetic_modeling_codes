function neutral_line,mag,vx,vy,vz
;
     FUX=vx*mag.bz
     FUY=vy*mag.bz
;
     FVX=-vz*mag.bx
     FVY=-vz*mag.by

     FU_PERP=FUX*MAG.BZX+FUY*MAG.BZY
     FU_PAR=-FUX*MAG.BZY+FUY*MAG.BZX
     FV_PERP=FVX*MAG.BZX+FVY*MAG.BZY
     FV_PAR=-FVX*MAG.BZY+FVY*MAG.BZX
     return,{FUX:FUX,FUY:FUY,FVX:FVX,FVY:FVY,$
                  FU_PERP:FU_PERP,FU_PAR:FU_PAR,$
                  FV_PERP:FV_PERP,FV_PAR:FV_PAR}

end

