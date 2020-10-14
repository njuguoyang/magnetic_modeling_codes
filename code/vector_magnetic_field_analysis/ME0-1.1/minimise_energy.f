c***********************************************************************************************************************************
      subroutine minimise_energy
      implicit none
c===================================================================================================================================
c This subroutine attempts to find the configuration of azimuthal angles that corresponds to the minimum of the "energy"
c summed over the region of interest. Simulated annealing is used to perform the search with a "temperature" that can vary from
c pixel to pixel. There are 3 main input parameters:
c tfac0 scales the initial temperature
c neq controls the number of reconfigurations that are attempted at each temperature setting
c tfactr is the cooling rate
c
c 2009, Ashley Crouch, ash@cora.nwra.com
c===================================================================================================================================
      include "include_file.sizes.f"
      include "include_file.bounds.f"
      include "include_file.seed.f"
      include "include_file.verb.f"
      include "include_file.anneal.f"
      include "include_file.ran_pix.f"
c
c Parameters tol_conv and nconv_min: If the relative changes in the energy is less than tol_conv for nconv_min consecutative
c temperature increments then the search has converged.
c
      real tol_conv
      parameter (tol_conv=1e-5)
      integer nconv_min
      parameter (nconv_min=10)
c
c Parameter tstop_par when the temperature is less than tstop_par times the initial temperature stop the search
c
      real tstop_par
      parameter (tstop_par=1e-7)

      integer i,j,idone,itemp,nconv,nsucc,nxnynz,neq_init,ii,jj
      real E,E_prev,de,ran3,max_de,tstop,t
      real tvar(nxmax,nymax)
c
c Initialisations
c
      call CalcE(E)

      nxny=dnx*dny
      nxnynz=nxny
      neq=neq*nxnynz
      neq_init=100*nxnynz

      nxjump=dnx/jump
      if (nxjump.lt.1) nxjump=1
      nyjump=dny/jump
      if (nyjump.lt.1) nyjump=1
      neq_init=neq_init/(nxjump*nyjump)
      neq=neq/(nxjump*nyjump)
      ia_prev=-1
      ja_prev=-1
c
c Calculate the initial temperature by sampling a very large number of reconfigurations (all accepted)
c
      do i=nxa,nxb
         do j=nya,nyb
            tvar(i,j)=0.
         enddo
      enddo
      max_de=0.
      do itemp=1,neq_init
c
c Loop over the sequence of pixels with random start point and separation: jump
c
         call get_ran_pix
         do ii=1,nxg
            i=ivec(ii)
            do jj=1,nyg
               j=jvec(jj)
c
c Calculate the change in the energy due the reconfiguration at the selected pixel
c
               call CalcDE(i,j,de)
               if (abs(de).gt.max_de) max_de=abs(de)
               if (abs(de).gt.tvar(i,j)) tvar(i,j)=abs(de)
c
c Update the total energy (for the initial temperature accept all reconfigurations)
c
               E=E+de
c
c Implement the reconfiguration
c
               call reconfig(i,j)
            enddo
         enddo
      enddo
c
c Scale the initial temperature according to input parameter tfac0
c
      do i=nxa,nxb
         do j=nya,nyb
            tvar(i,j)=tfac0*tvar(i,j)
         enddo
      enddo
      t=tfac0*max_de
c
c More initialisations
c
      tstop=tstop_par*t
      call CalcE(E)
      E_prev=E
      nconv=0
      idone=0
      if (iverb.eq.2) write(6,'(i9,3x,e15.8,3x,e15.8,3x,i3)')0,t*tstop_par/tstop,E,nconv
      do while (idone.eq.0)
         nsucc=0
         do itemp=1,neq
c
c Loop over the sequence of pixels with random start point and separation: jump
c
            call get_ran_pix
            do ii=1,nxg
               i=ivec(ii)
               do jj=1,nyg
                  j=jvec(jj)
c
c Calculate the change in the energy due the reconfiguration at the selected pixel
c
                  call CalcDE(i,j,de)
c
c Decide whether or not to accept the reconfiguration
c
                  if ((de.lt.0.).or.(ran3(seed).lt.exp(-de/tvar(i,j)))) then !metropolis criteria (simulated annealing, spatially dependent)
                     nsucc=nsucc+1
c
c Update the energy
c
                     E=E+de
c
c Implement the reconfiguration
c
                     call reconfig(i,j)
                  endif
               enddo
            enddo
         enddo
c
c Reduce the temperature according to input parameter tfactr
c
         do i=nxa,nxb
            do j=nya,nyb
               tvar(i,j)=tfactr*tvar(i,j)
            enddo
         enddo
         t=t*tfactr
c
c Convergence test
c
         if (abs(E-E_prev)/(abs(E)+abs(E_prev)).lt.tol_conv) then
            nconv=nconv+1
         else
            nconv=0
         endif
         E_prev=E
         if (iverb.eq.2) write(6,'(i9,3x,e15.8,3x,e15.8,3x,i3)')nsucc,t*tstop_par/tstop,E,nconv
c
c Stopping conditions:
c 1. no reconfigurations were accepted during this temperature increment
c 2. temperature is "small"
c 3. the change in energy is "small" over several consecutive temperature steps (converged)
c
         if ((nsucc.eq.0).or.(t.lt.tstop).or.(nconv.ge.nconv_min)) idone=1
      enddo
      return
      end
c***********************************************************************************************************************************
