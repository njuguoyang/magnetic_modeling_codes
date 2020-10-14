pro parallel
;
;    make figure to demonstrate how DAVE4VM estimates parallel
;    velocities is response to Pascal Demoulin's second salvo of
;    questions ;-)
;
     common loop,Z0,A
;
     N=11L
     I0=1.d0
     A=1.d0

     XM=2.5d0
     ZM=2.d0/1.d0
     Z0=-1.d0

     X=rebin((findgen(N)-N/2)*XM,N,N)/double(N/2)
     Z=transpose(rebin((findgen(N))*ZM,N,N))/double(N)

     BZ=B_Z(X,Z)
     BX=B_X(X,Z)
     vx=replicate(-1.d0,N,N)
     vz=replicate(1.d0,N,N)

     vdb=(vx*Bx+vz*Bz)/(Bx^2+Bz^2)
     vpx=vx-vdb*Bx
     vpz=vz-vdb*Bz
     set_plot,'ps'
     !P.font=0
     dummy=FSC_Color(/allcolors,colorstructure=ps)
     th=4
     cth=2
     cs=1.25
     dz=z[0,1]-z[0,0]
     ytitle='h/R'+sunsymbol()
     extra={xthick:th,ythick:th,charthick:cth,charsize:cs}

     device,xsize=6,ysize=3,filename='SCHUCK/f6.eps',/inches,/color,bits_per_pixel=8
     velovect,BX,BZ,X[*,0],Z[0,*],length=1,color=ps.black,$
              thick=th,_extra=extra,xtitle='X',ytitle=ytitle,ticklen=-0.02
     polyfill,[-1.4,-1.4,1.4,1.4],[-dz,dz,dz,-dz],color=ps.gray
     velovect,BX,BZ,X[*,0],Z[0,*],length=1,color=ps.black,$
              thick=th,_extra=extra,xtitle='X',ytitle=ytitle,ticklen=-0.02,/noerase
     mx=max(sqrt(vx^2+vz^2))
     mxp=max(sqrt(vpx^2+vpz^2))
     velovect,VX,VZ,X[*,0],Z[0,*],length=1,/overplot,color=ps.red,thick=th
     velovect,VPX,VPZ,X[*,0],Z[0,*],length=mxp/mx,/overplot,color=ps.blue,thick=th
     device,/close_file
;
;       Set up arrays to demonstrate equation~21 in manuscript
;
     i1=(where(X[*,0] eq -1))[0]          
     i2=(where(X[*,0] eq 1))[0]
     N=i2-i1+1L
     ind=i1+lindgen(N)

     A=dblarr(2,2*N)
     VP=dblarr(2*N)
     B2=Bx^2+Bz^2
     for j=0L,N-1 do begin
         i=2*J
         A[*,i]=[Bz[ind[j],0]^2/B2[ind[j],0],-Bx[ind[j],0]*Bz[ind[j],0]/B2[ind[j],0]]
         A[*,i+1]=[-Bx[ind[j],0]*Bz[ind[j],0]/B2[ind[j],0],Bx[ind[j],0]^2/B2[ind[j],0]]
         VP[i:i+1]=[VPX[ind[j],0],VPZ[ind[j],0]]
     endfor
     vt=LA_least_squares(A,VP,/double,method=3,rank=rank,rcondition=rcondition,residual=residual,status=status)

     print,'Total velocity inferred from data in Figure 6'
     print,vt,format='("total velocity=(",F8.2,",",F8.2,") which should be (-1,1)")


end

