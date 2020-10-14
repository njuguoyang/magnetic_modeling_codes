pro make_all,no_elliptic2d=no_elliptic2d,clean=clean,skip_optimize=skip_optimize,pazo=pazo
;
;    CONTROL program for making figures for P.W. Schuck, "Tracking Vector
;    Magnetograms with the Magnetic Induction Equation", ApJ, 2008.
;
;    usage:
;      IDL>make_all/no_elliptic_2d
;
;    KEYWORDS        VALUE
;      NO_ELLIPTIC:
;                      0 - perform elliptic solves using mudpack
;                          This assumes that either the compiled
;                          version of MUDPACK distributed with the
;                          release works for your system OR you have 
;                          ABSOFT fortran and are able to compile it 
;                          from the source.
;                      1 - (RECOMMENDED) read in saved vector
;                          potentials instead of computing them with MUDPACK
;
;
;      CLEAN           0 - Don't erase files
;                      1 - Erase all plots, tables, and comparison
;                          files and generate them all from scratch
;
;      SKIP_OPTIMIZE   0 - Don't skip
;
;
;      PAZO            0 - (RECOMMENDED) PAZO math fonts are not
;                          present
;                      1 - PAZO math fonts present
;
;                          If you want the labels as produced for the
;                          paper you'll need to install the
;                          fonts on your machine AND make them
;                          accessible by IDL
;                           
;        1) get fonts
;           http://ftp.math.purdue.edu/mirrors/ctan.org/fonts/mathpazo/
;        2) edit resources/fonts/ps/font.map under your idl
;           distribution
;          
;           I have the following entries:
;      PazoMath                        PazoMath                        1
;      PazoMath-Italic                 PazoMathI                       1
;      PazoMath-Bold                   PazoMathBo                      1
;      PazoMath-BoldItalic             PazoMathBoI                     1
;  
;      AND the follwing AFM files
;
;      resource/fonts/ps/PazoMath.afm
;      resource/fonts/ps/PazoMathBlackBoardBo.afm
;      resource/fonts/ps/PazoMathBo.afm
;      resource/fonts/ps/PazoMathBoI.afm
;      resource/fonts/ps/PazoMathI.afm
;
;      generated with the idl routine PSAFM with the postscript AFM
;      files fplm*.afm at http://ftp.math.purdue.edu/mirrors/ctan.org/fonts/mathpazo/afm/
;
;
;
  if (N_ELEMENTS(clean) eq 0) then clean=0
  if (clean eq 1) then begin
     print,'CLEANING directories...'
     spawn,'\rm -f SCHUCK/*.eps'
     spawn,'\rm -f WELSCH/*.eps'
     spawn,'\rm -f SCHUCK/*.txt'
     spawn,'\rm -f WELSCH/*.txt'
     spawn,'\rm -f SCHUCK/*.tex'
     spawn,'\rm -f WELSCH/*.tex'
     print,'DONE'
     return
  endif
;
  if (N_ELEMENTS(skip_optimize) eq 0) then skip_optimize=0
;
;     accumulate results & make figures
;         OLD mask |Bz|>370~G
  print & print & print
  print,'compare DAVE and DAVE4VM on the original mask |Bz|>370~G'
  dave_results,/welsch,no_elliptic2d=no_elliptic2d
  dave4vm_results,/welsch,no_elliptic2d=no_elliptic2d
;
;         NEW mask |B|=sqrt(Bx^2+By^2+Bz^2)>370~G
  print & print & print
  print,'compare DAVE and DAVE4VM on the new mask |B|>370~G'
  dave_results,welsch=0,no_elliptic2d=no_elliptic2d
;
;     first run dave4vm with no horizontal magnetic fields
  dave4vm_results,welsch=0,/dave,no_elliptic2d=no_elliptic2d 
;
;     save results
  spawn,'cp SCHUCK/dave4vm.txt SCHUCK/dave4vm_Bx=By=0.txt'
;
;     now run it again with horizontal fields
  dave4vm_results,welsch=0,no_elliptic2d=no_elliptic2d
;
; verify code against original results
;     should have no differences (or minor ones cause by processor differences)
  spawn,'diff SCHUCK/dave.txt SCHUCK/dave.cmp  >SCHUCK/dave_err.txt'
;     should have no differences (or minor ones cause by processor differences)
  spawn,'diff SCHUCK/dave4vm.txt SCHUCK/dave4vm.cmp >SCHUCK/dave4vm_err.txt'
;     should have no differences (or minor ones cause by processor differences)
  spawn,'diff WELSCH/dave.txt WELSCH/dave.cmp  >WELSCH/dave_err.txt'
;     should have no differences (or minor ones cause by processor differences)
  spawn,'diff WELSCH/dave4vm.txt WELSCH/dave4vm.cmp >WELSCH/dave4vm_err.txt'
;
;       should have differences in *** TOTAL *** plasma velocities
  spawn,'diff SCHUCK/dave4vm_Bx=By=0.txt SCHUCK/dave.txt> dave_just_plasma_velocities.txt'
;     
  print & print & print
  if (skip_optimize eq 0) then begin
     print,'Perform WINDOW optimization calculations for DAVE and DAVE4VM'
     dave_optimize,no_elliptic2d=no_elliptic2d
     dave4vm_optimize,no_elliptic2d=no_elliptic2d
 endif else print,'*** Skipping WINDOW optimization calculations for DAVE and DAVE4VM'
;
;     make calibration curves for window size
  dave_calibration
  dave4vm_calibration

;     make tables
  print & print & print
  print,'Make tables DAVE and DAVE4VM on both masks'
  DIR='SCHUCK/'
  make_table,'DAVE',DIR+'dave.fts','DAVE4VM',DIR+'dave4vm.fts',DIR+'weak.tex'
  make_vtable,'DAVE',DIR+'dave.fts','DAVE4VM',DIR+'dave4vm.fts',DIR+'vweak.tex'
  DIR='WELSCH/'
  make_table,'DAVE',DIR+'dave.fts','DAVE4VM',DIR+'dave4vm.fts',DIR+'strong.tex'
  make_vtable,'DAVE',DIR+'dave.fts','DAVE4VM',DIR+'dave4vm.fts',DIR+'vstrong.tex'
;
  print & print & print
  print,'Verify codes against original results'
;       should have no differences (except perhaps tiny ones cause by
;       processor differences)
  print,'List of differences in DAVE     |B|>370 G: (minor differences okay)'
  spawn,'cat SCHUCK/dave_err.txt'
;       should have no differences (except perhaps tiny ones cause by
;       processor differences)
  print,'List of differences in DAVE4VM  |B|>370 G: (minor differences okay)'
  spawn,'cat SCHUCK/dave4vm_err.txt'
;       should have no differences (except perhaps tiny ones cause by
;       processor differences)
  print & print
  print,'List of differences in DAVE    |Bz|>370 G: (minor differences okay)'
  spawn,'cat WELSCH/dave_err.txt'
;       should have no differences (except perhaps tiny ones cause by
;       processor differences)
  print,'List of differences in DAVE4VM |Bz|>370 G: (minor differences okay)'
  spawn,'cat WELSCH/dave4vm_err.txt'
;       should have differences in TOTAL velocity 
;       see:  dave_just_plasma_velocities.cmp
  print & print
  print,'DAVE velocities are just perpendicular plasma velocities:'
  print,'This should show differences in total plasma velocity between'
  print,'DAVE and DAVE4VM'
  spawn,'cat dave_just_plasma_velocities.txt'
;
  print & print
  print,"make figures 8, 13, & 14"
  redundent,no_elliptic2d=no_elliptic2d,pazo=pazo
;
  print & print
  print,"make figure 6"
  parallel
;
  print,'**************************************'
  print,'*** DAVE-DAVE4VM figures completed ***'
  print,'***                                ***'
  print,'*** FLAGS:                         ***'
  print,'***   NO_ELLIPTIC2D: ',string(NO_ELLIPTIC2D,format='(i1)'),'             ***'
  print,'***   SKIP_OPTIMIZE: ',string(SKIP_OPTIMIZE,format='(i1)'),'             ***'
  print,'***            PAZO: ',string(PAZO,format='(i1)'),'             ***'

  print,'**************************************'
;     
  

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
