PRO GET_EXTRAP_COMPONENTS,Bc,Bex

; PURPOSE: Get the extrapolated components Bc from the output Bex of
; the LFF routine
; PROGRAMMER: Manolis K. Georgoulis (JHU/APL,10/13/05)

if n_elements(size(Bex)) eq 6. then begin 
  Bc(0,*,*)=Bex(*,*,0)
  Bc(1,*,*)=Bex(*,*,1)
  Bc(2,*,*)=Bex(*,*,2)
endif else begin 
  Bc(0,*,*)=Bex(*,*,0,0)
  Bc(1,*,*)=Bex(*,*,0,1)
  Bc(2,*,*)=Bex(*,*,0,2)
endelse 

RETURN
END
