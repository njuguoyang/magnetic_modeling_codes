pro mudpack_engine,mudpack,coefs,psi

   sz=size(work)
   spawn,'mktemp /tmp/map.sdf.XXXXXXXXX',infile
   spawn,'mktemp /tmp/phi.sdf.XXXXXXXXX',outfile
;
   print,'Input file name:',infile
   print,'Output file name:',outfile
;
   mudpack.loud=1
;
   iparms=long([mudpack.nx,mudpack.ny,mudpack.ixp,mudpack.jyq,mudpack.iex,$
           mudpack.jey,mudpack.iguess,mudpack.maxcy,mudpack.method,$
           mudpack.loud,mudpack.boundary])
;
   dparms=double([mudpack.xr,mudpack.yr,0])
;
   spawn,'\rm -f '+infile
   spawn,'\rm -f '+outfile
;
   sdf_write,infile,'iparms',iparms
   sdf_write,infile,'dparms',dparms
   sdf_write,infile,'coefs',coefs
;
   print,'spawning mudpack job (mudpack engine)'
;   cmd='$HOME/DAVE-DAVE4VM/MUD/elliptic2d '+infile+' '+outfile         ;comment by YG 20101011
   cmd='$HOME/idl/code/DAVE-DAVE4VM/MUD/elliptic2d '+infile+' '+outfile
   openw,1,'mybat.sh'
   printf,1,'#!/bin/csh'
   printf,1,'unlimit stacksize'
   printf,1,cmd
   printf,1,'exit 0'
   close,1
   spawn,'./mybat.sh'
   print,'spawning command:',cmd
;
   sdf_read,outfile,0,label0,iparms
   sdf_read,outfile,1,label1,ierror
   sdf_read,outfile,2,label2,dparms
   sdf_read,outfile,3,label3,psi
;
   spawn,'\rm -f '+infile
   spawn,'\rm -f '+outfile

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;      old way
;   openw,fn,infile,/f77_unformatted,/get_lun
;      writeu,fn,[mudpack.nx,mudpack.ny,mudpack.ixp,mudpack.jyq,mudpack.iex,mudpack.jey,mudpack.iguess,mudpack.maxcy,mudpack.method,mudpack.loud,mudpack.boundary]
;      writeu,fn,double([mudpack.xr,mudpack.yr,0])
;      writeu,fn,double(work)
;   close,fn
;   free_lun,fn
;
;   print,'spawning mudpack job (mudpack engine)'
;   print,'elliptic2d '+infile+' '+outfile
;   spawn,'$HOME/MUD/elliptic2d '+infile+' '+outfile
;
;   psi=dblarr(sz[1],sz[2])
;   iparms=lonarr(17+1)
;   fparms=dblarr(5)
;   openr,fn,outfile,/f77_unformatted,/get_lun
;      readu,fn,iparms
;      readu,fn,fparms
;      readu,fn,psi
;   close,fn
;   free_lun,fn
