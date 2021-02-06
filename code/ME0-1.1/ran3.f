c***********************************************************************************************************************************
      function ran3(seed)
c===================================================================================================================================
c     "Minimal standard" pseudo-random number generator of Park and
c     Miller.  Returns a uniform random deviate r s.t. 0 < r < 1.0.
c     Set seed to any non-zero integer value to initialize a sequence,
c     then do not change seed between calls for successive deviates
c     in the sequence.
c
c     References:
c        Park, S. and Miller, K., "Random Number Generators: Good Ones
c           are Hard to Find", Comm. ACM 31, 1192-1201 (Oct. 1988)
c        Park, S. and Miller, K., in "Remarks on Choosing and Imple-
c           menting Random Number Generators", Comm. ACM 36 No. 7,
c           105-110 (July 1993)
c
c This is ran0 from the genetic algorithm code PIKAIA developed by Paul Charbonneau & Barry Knapp
c===================================================================================================================================
c *** Declaration section ***
c
      implicit none
c
c     Input/Output:
      integer seed
c
c     Output:
      real ran3
c
c     Constants:
      integer A,M,Q,R
      parameter (A=48271,M=2147483647,Q=44488,R=3399)
      real SCALE,EPS,RNMX
      parameter (SCALE=1./M,EPS=1.2e-7,RNMX=1.-EPS)
c
c     Local:
      integer j
c
c *** Executable section ***
c
      j = seed/Q
      seed = A*(seed-j*Q)-R*j
      if (seed .lt. 0) seed = seed+M
      ran3 = min(seed*SCALE,RNMX)
      return
      end
c***********************************************************************************************************************************
