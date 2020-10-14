FC=ifort
#MKLROOT may be set automatically, but if not, uncomment the following, and insert appropriate path
#MKLROOT = /opt/intel/Compiler/11.1/046/mkl

ambig_files=ambig.f boxit.f buffer.f sortrx.f nacute5.f nacute6.f pacute.f potential.f printerror.f rd_field.f rd_fits.f rd_param.f rd_sotsp.f smooth.f cuspl.f cuspl2d.f cuspleval.f cuspleval2d.f transform.f update_fits.f write_field.f setup.f minimise_energy.f get_ran_pix.f ran3.f CalcE.f CalcDE.f CalcDE_reconfig.f reconfig.f

ambig: 	$(ambig_files) mkl_dfti.mod dfti_example_status_print.o
	$(FC) -extend_source -w -I/usr/include/cfitsio -O3 -ipo -axSSE4.2,SSE4.1,SSSE3 -mcmodel medium -shared-intel -o $@ $(ambig_files)  dfti_example_status_print.o -mkl -lcfitsio

clean:
	rm -f ambig *.mod *.o

mkl_dfti.mod: $(MKLROOT)/include/mkl_dfti.f90
	$(FC) -O3 -ipo -axSSE4.2,SSE4.1,SSSE3 -mcmodel medium -shared-intel -c $^

%.o: %.f90
	$(FC) -O3 -ipo -axSSE4.2,SSE4.1,SSSE3 -mcmodel medium -shared-intel -c $<
