       PROGRAM FUN

       INTEGER N,INCX
       DOUBLE PRECISION ALPHA,TAU,SAFMIN,VNORM,CUTLO,CUTHI,DDIST2,DNRM2
       DOUBLE PRECISION XVECTOR(10)

       N=10
       INCX=1
       ALPHA = 2.0
       TAU = 2.0
       SAFMIN = 0.2225073858507201D-307
       CUTLO = 4.44089D-16/N
       CUTHI = 1.30438D19*10.0

       XVECTOR(1)=1.0;
       XVECTOR(2)=1.0;
       XVECTOR(3)= CUTHI;
       XVECTOR(4)=1.0;
       XVECTOR(5)=1.0;
       XVECTOR(6)= CUTHI;
       XVECTOR(7)=1.0;
       XVECTOR(8)= CUTHI;
       XVECTOR(9)=1.0;

       CALL HOUSG(N, ALPHA, XVECTOR, INCX, TAU, SAFMIN)
C       VNORM = DDIST2( N-1,XVECTOR,INCX,XVECTOR,INCX,1)
       VNORM = DNRM2( N-1,XVECTOR,INCX)
       
       write(6,10) -2.0*tau + tau*tau*(1.0+vnorm*vnorm)
   10  FORMAT (G22.16)
       STOP
       END
C
C234567==1=========2=========3==========4=========5=========6=========7==
C
      double precision function dnrm2 ( n, dx, incx)
      integer          next
      double precision   dx(1), cutlo, cuthi, hitest, sum, xmax,zero,one
      data   zero, one /0.0d0, 1.0d0/
c
c     euclidean norm of the n-vector stored in dx() with storage
c     increment incx .
c     if    n .le. 0 return with result = 0.
c     if n .ge. 1 then incx must be .ge. 1
c
c           c.l.lawson, 1978 jan 08
c
c     four phase method     using two built-in constants that are
c     hopefully applicable to all machines.
c         cutlo = maximum of  dsqrt(u/eps)  over all known machines.
c         cuthi = minimum of  dsqrt(v)      over all known machines.
c     where
c         eps = smallest no. such that eps + 1. .gt. 1.
c         u   = smallest positive no.   (underflow limit)
c         v   = largest  no.            (overflow  limit)
c
c     brief outline of algorithm..
c
c     phase 1    scans zero components.
c     move to phase 2 when a component is nonzero and .le. cutlo
c     move to phase 3 when a component is .gt. cutlo
c     move to phase 4 when a component is .ge. cuthi/m
c     where m = n for x() real and m = 2*n for complex.
c
c     values for cutlo and cuthi..
c     from the environmental parameters listed in the imsl converter
c     document the limiting values are as follows..
c     cutlo, s.p.   u/eps = 2**(-102) for  honeywell.  close seconds are
c                   univac and dec at 2**(-103)
c                   thus cutlo = 2**(-51) = 4.44089e-16
c     cuthi, s.p.   v = 2**127 for univac, honeywell, and dec.
c                   thus cuthi = 2**(63.5) = 1.30438e19
c     cutlo, d.p.   u/eps = 2**(-67) for honeywell and dec.
c                   thus cutlo = 2**(-33.5) = 8.23181d-11
c     cuthi, d.p.   same as s.p.  cuthi = 1.30438d19
c     data cutlo, cuthi / 8.232d-11,  1.304d19 /
c     data cutlo, cuthi / 4.441e-16,  1.304e19 /
      data cutlo, cuthi / 8.232d-11,  1.304d19 /
c
      if(n .gt. 0) go to 10
         dnrm2  = zero
         go to 300
c
   10 assign 30 to next
      sum = zero
      nn = n * incx
c                                                 begin main loop
      i = 1
   20    go to next,(30, 50, 70, 110)
   30 if( dabs(dx(i)) .gt. cutlo) go to 85
      assign 50 to next
      xmax = zero
c
c                        phase 1.  sum is zero
c
   50 if( dx(i) .eq. zero) go to 200
      if( dabs(dx(i)) .gt. cutlo) go to 85
c
c                                prepare for phase 2.
      assign 70 to next
      go to 105
c
c                                prepare for phase 4.
c
  100 i = j
      assign 110 to next
      sum = (sum / dx(i)) / dx(i)
  105 xmax = dabs(dx(i))
      go to 115
c
c                   phase 2.  sum is small.
c                             scale to avoid destructive underflow.
c
   70 if( dabs(dx(i)) .gt. cutlo ) go to 75
c
c                     common code for phases 2 and 4.
c                     in phase 4 sum is large.  scale to avoid overflow.
c
  110 if( dabs(dx(i)) .le. xmax ) go to 115
         sum = one + sum * (xmax / dx(i))**2
         xmax = dabs(dx(i))
         go to 200
c
  115 sum = sum + (dx(i)/xmax)**2
      go to 200
c
c
c                  prepare for phase 3.
c
   75 sum = (sum * xmax) * xmax
c
c
c     for real or d.p. set hitest = cuthi/n
c     for complex      set hitest = cuthi/(2*n)
c
   85 hitest = cuthi/float( n )
c
c                   phase 3.  sum is mid-range.  no scaling.
c
      do 95 j =i,nn,incx
      if(dabs(dx(j)) .ge. hitest) go to 100
   95    sum = sum + dx(j)**2
      dnrm2 = dsqrt( sum )
      go to 300
c
  200 continue
      i = i + incx
      if ( i .le. nn ) go to 20
c
c              end of main loop.
c
c              compute square root and adjust for scaling.
c
      dnrm2 = xmax * dsqrt(sum)
  300 continue
      return
      end
C
C234567==1=========2=========3=========4=========5=========6=========7==
C
C                            MANSUP         
C                           Version 3               March 15, 1996
C
C   A package of subroutines for support of MANPAK (VERS. 7) :
C
C   1. ODE Solvers
C      -----------
C
C   SUBROUTINE DOPSTA  ... Dormand-Prince order 5 RK-step routine
C                          for autonomous ODEs using reverse
C                          communication for function calls 
C
C   SUBROUTINE DOPSTN  ... Dormand-Prince order 5 RK-step routine
C                          for non-autonomous ODEs using reverse
C                          communication for function calls
C
C   2. Linear algebra routines
C      -----------------------
C
C   SUBROUTINE BIDIA   ... Reduces a general M by N matrix A, M >= N 
C                          to an M by N upper bidiagonal matrix
C
C   DOUBLE PRECISION FUNCTION DDIST2 
C                      ... Computes either the Euclidean distance 
C                          between two vectors X and Y or the 
C                          Euclidean norm of one such vector X
C
C   SUBROUTINE HOUSG   ... Generates a Householder reflector
C
C   SUBROUTINE HOUSL   ... Multiplies a matrix A from the left by 
C                          a Householder reflector
C
C   SUBROUTINE HOUSR   ... Multiplies a matrix A from the right by 
C                          a Householder reflector
C
C   SUBROUTINE LQF     ... Computes the LQ factorization with column 
C                          pivoting of an M by N matrix, M <= N.
C
C   SUBROUTINE LQAS    ... For an M x N matrix A with rank A = M <= N
C                          and given M-vector y, computes the solution 
C                          x of A x = y orthogonal to ker A
C
C   SUBROUTINE LUF     ... Computes an LU factorization with row 
C                          pivoting of an N-by-N matrix A
C
C   SUBROUTINE LUS1    ... Solves a linear system A * Z = X for an
C                          LU factored matrix A and a given vector X 
C
C   SUBROUTINE LUSK    ... Solves a linear system A * X = B for an
C                          LU factored matrix A and a given 
C                          matrix B on the right side
C
C   SUBROUTINE QRF     ... Computes the QR factorization with 
C                          column pivoting of a matrix A
C
C   SUBROUTINE QORG    ... Generate a matrix Q with orthonormal 
C                          columns as the product of given Householder 
C                          reflectors
C
C   SUBROUTINE QRS     ... Computes the least squares solution for
C                          an  M x N matrix A, M >= N, rank A = N, 
C                          and a given M-vector Y 
C
C   SUBROUTINE ROTG    ... Generates a plane rotation
C
C   SUBROUTINE ROTML   ... Multiplies a given matrix by a plane
C                          rotation from the left
C
C   SUBROUTINE ROTMR   ... Multiplies a given matrix by a plane
C                          rotation from the right
C
C   SUBROUTINE SVD     ... Computes the singular value decomposition 
C                          of an M x N, M >= N, upper bidiagonal matrix
C
C   SUBROUTINE SVD2D   ... Computes the singular value decomposition 
C                          of a 2-by-2 triangular matrix
C
C   3. Other support routines
C      ----------------------
C
C   SUBROUTINE MSGPRT  ... for uniform printing of messages
C
C   SUBROUTINE LUNIT   ... returns KL output-unit numbers 
C
C   SUBROUTINE DMACH   ... returns the relative machine precision 
C                          EPMACH and the safe minimum SAFMIN
C
C234567==1=========2=========3=========4=========5=========6=========7==
C
      SUBROUTINE DOPSTA( TASK,N,Y,YP,H,HMIN,HMAX,NMAX,ATOL,RTOL,ITOL,
     &                   W0,W1,W2,W3,W4,W5,W6,NSTEP,NACCPT,NREJCT )
C
      CHARACTER*6 TASK
      INTEGER ITOL,N,NACCPT,NMAX,NREJCT,NSTEP 
      DOUBLE PRECISION Y(*),YP(*),H,HMIN,HMAX,RTOL(*),ATOL(*)
      DOUBLE PRECISION W0(*),W1(*),W2(*),W3(*),W4(*),W5(*),W6(*)
C 
C234567--x---------x---------x---------x---------x---------x---------x--
C
C  Routine for taking a step with the Dormand-Prince Runge Kutta
C  method of order 5 in solving approximately the AUTONOMOUS 
C  initial value problem
C
C        y' = F(y),  y(0) = y0,
C
C  where
C        y    vector of dimension N,
C        F    mapping from R^N into R^N
C
C  The routine uses reverse communication and returns for all
C  evaluations of the user function for F. This is done under the
C  control of the character variable TASK.
C
C  To take a step compute YP = F(Y) for the given point Y and
C  call the routine with TASK = 'STEP'.
C  The routine returns with TASK = 'EVAL' signifying a request for
C  the evaluation of F at Y. 
C  Return to the routine with the computed vector YP and without  
C  changing TASK or any other variable in the calling sequence.
C  Alternately, without computing YP, a forced step reduction may be
C  requested by returning to the routine with TASK = 'REDUCE' and
C  with a suggested smaller stepsize in H. As default, when H is 
C  unchanged the step size is reduced by a factor 1/2. This
C  feature is included to allow for the solution of systems 
C  defined in terms of local coordinates on a manifold.
C
C  Upon final return TASK may have either one of the following 
C  values:
C    TASK = 'MINSTP'  step fell below minimum steplength HMIN
C    TASK = 'STPCNT'  maximal number of steps NMAX exceeded
C    TASK = 'ERROR'   error condition.
C
C  Variables in the calling sequence:
C  ----------------------------------
C  TASK   C  IN   Task identifier
C                 TASK = 'STEP'   take a new RK step
C                 TASK = 'EVAL'   YP = F(Y) for given Y has been
C                                 computed 
C                 TASK = 'REDUCE' force a step reduction
C            OUT  TASK = 'EVAL'   Compute YP = F(Y) for given Y
C                 TASK = 'DONE'   RK step successfully completed
C                 TASK = 'MINSTP' Step fell below allowed minimum
C                 TASK = 'STPCNT' Maximal number of steps exceeded
C                 TASK = 'ERROR'  Other error condition 
C  N      I  IN   Dimension of Y
C  Y      D  IN   Array of dimension N, the current point
C            OUT  The next computed point
C  YP     D  IN   Array of dimension N, function value at Y 
C            OUT  Function value at the next computed point 
C  H      D  IN   Given stepsize
C         D  OUT  Predicted next stepsize
C  NMAX   I  IN   Maximal number of steps
C  RTOL   D  IN   Relative error tolerance. May be
C                 an array of dimension 1 or N
C  ATOL   D  IN   Absolute error tolerance. May be
C                 an array of dimension 1 or N
C  ITOL   I  IN   Switch for RTOL and ATOL
C                 ITOL = 0: Both RTOL and ATOL are scalars
C                           The code keeps the local error of Y(I),
C                           roughly, below RTOL(1)*ABS(Y(I)) + ATOL(1) 
C                 ITOL = 1: Both RTOL and ATOL are vectors
C                           The code keeps the local error of Y(I),
C                           roughly, below RTOL(I)*ABS(Y(I)) + ATOL(I) 
C  W0-W6  D  WK   Seven work arrays of dimension N
C  NSTEP  I  I/O  Step counter
C  NACCPT I  I/O  Counter of accepted steps
C  NREJCT I  I/O  Counter of rejected steps (after the first step)
C
C234567--x---------x---------x---------x---------x---------x---------x--
C
C.....Step parameters
C
      DOUBLE PRECISION BETA,FAC1,FAC2,SAFE,EXPO1,FACMIN,FACIN1,FACIN2 
      PARAMETER( BETA=0.04D0, FAC1=0.2D0, FAC2=10.D0, SAFE=0.9D0,
     &           EXPO1=0.2D0 - BETA*0.75D0, FACMIN=1.0D-4, 
     &           FACIN1 = 1.0D0/FAC1, FACIN2 = 1.0D0/FAC2 )
C
C.....Dormand-Prince RK-5 coefficients for the autonomous case
C
      DOUBLE PRECISION A21,A31,A32,A41,A42,A43,A51,A52,A53,A54
      DOUBLE PRECISION A61,A62,A63,A64,A65,A71,A73,A74,A75,A76
      DOUBLE PRECISION E1,E3,E4,E5,E6,E7
C
      PARAMETER( A21 = 0.2D0,A31 = 3.0D0/40.0D0,
     &        A32 = 9.0D0/40.0D0, A41 = 44.0D0/45.0D0,
     &        A42 = -56.0D0/15.0D0, A43 = 32.0D0/9.0D0,
     &        A51 = 19372.0D0/6561.0D0, A52 = -25360.0D0/2187.0D0,
     &        A53 = 64448.0D0/6561.0D0, A54 = -212.0D0/729.0D0,
     &        A61 = 9017.0D0/3168.0D0, A62 = -355.0D0/33.0D0,
     &        A63 = 46732.0D0/5247.0D0, A64 = 49.0D0/176.0D0,
     &        A65 = -5103.0D0/18656.0D0, A71 = 35.0D0/384.0D0,
     &        A73 = 500.0D0/1113.0D0, A74 = 125.0D0/192.0D0,
     &        A75 = -2187.0D0/6784.0D0,A76 = 11.0D0/84.0D0 )
      PARAMETER( E1  = 71.0D0/57600.0D0, E3  = -71.0D0/16695.0D0,
     &        E4  = 71.0D0/1920.0D0, E5  = -17253.0D0/339200.0D0,
     &        E6  = 22.0D0/525.0D0,E7  = -1.0D0/40.0D0 )
C
C.....Other parameters
C
      DOUBLE PRECISION ZER,HALF,ONE
      PARAMETER( ZER=0.0D0, HALF=0.5D0, ONE=1.0D0 )
C 
C.....Fortran functions called
C
      DOUBLE PRECISION ABS,MAX,MIN,SQRT
C
C.....Local variables
C
      INTEGER I
      DOUBLE PRECISION ERR,ERRFAC,FAC,HNEW,SK,TMP
C
C.....Variables to be saved between calls
C
      LOGICAL REJECT
      INTEGER ICALL
      DOUBLE PRECISION ATOLI,FACOLD,HAB,POSNEG,RTOLI
      SAVE REJECT,ICALL,ATOLI,FACOLD,HAB,POSNEG,RTOLI
C
      DATA FACOLD/1.0D-4/
C
C.......................Executable statements...........................
C
C.....Check TASK
C
      IF (TASK .EQ. 'EVAL') THEN
         GOTO 100
C
      ELSEIF (TASK .EQ. 'STEP') THEN
C
C........Check data
C
         IF (ITOL .EQ. 0) THEN
            ATOLI  = ATOL(1)
            RTOLI  = RTOL(1)
         ENDIF
         REJECT = .FALSE.
         POSNEG = ONE
         HAB    = H
         IF (H .LT. ZER) THEN
            POSNEG = -POSNEG
            HAB = -HAB
         ENDIF
C
C........Store initial Y and YP
C
         DO 10 I = 1, N
            W0(I) = Y(I)
            W1(I) = YP(I)
   10    CONTINUE
         GOTO 50
C
      ELSEIF (TASK .EQ. 'REDUCE') THEN
C
C........Forced step reduction
C
         IF (ABS(H) .GE. HAB) H = HALF*H
         REJECT = .TRUE.
         IF(NACCPT .GE. 1)NREJCT = NREJCT + 1
C
      ELSE
C
C........TASK has an improper value
C
         TASK = 'ERROR'
         RETURN
      ENDIF
C
C.....Return point for rejected steps
C
   50 CONTINUE
C
C.....Check for maximal number of steps
C
      NSTEP = NSTEP + 1
      IF (NSTEP .GT. NMAX) THEN
         TASK = 'STPCNT'
         RETURN
      ENDIF
C
C.....Check for minimal stepsize
C
      IF (ABS(H) .LT. HMIN) THEN
         TASK = 'MINSTP'
         RETURN
      ENDIF
      TASK = 'EVAL'
C
C.....Stage Evaluations
C
      DO 90 I = 1,N
         Y(I)  = W0(I) + H*A21*W1(I)
   90 CONTINUE
      ICALL = 1
      RETURN
C
  100 CONTINUE
      IF (ICALL .EQ. 1) THEN
         DO 110 I=1,N
            W2(I) = YP(I) 
            Y(I)  = W0(I) + H*(A31*W1(I) + A32*W2(I))
  110    CONTINUE
         ICALL = 2
         RETURN
C
      ELSEIF (ICALL .EQ. 2) THEN
         DO 120 I=1,N 
            W3(I) = YP(I) 
            Y(I)  = W0(I) + H*(A41*W1(I) + A42*W2(I) + A43*W3(I))
  120    CONTINUE
         ICALL = 3
         RETURN
C
      ELSEIF (ICALL .EQ. 3) THEN
         DO 130 I=1,N 
            W4(I) = YP(I) 
            Y(I) = W0(I) + H*(A51*W1(I) + A52*W2(I) + A53*W3(I)
     &                                              + A54*W4(I))
  130    CONTINUE
         ICALL = 4
         RETURN
C
      ELSEIF (ICALL .EQ. 4) THEN
         DO 140 I=1,N 
            W5(I) = YP(I) 
            Y(I) = W0(I) + H*(A61*W1(I) + A62*W2(I) + A63*W3(I)
     &                                  + A64*W4(I) + A65*W5(I))
  140    CONTINUE
         ICALL = 5
         RETURN
C
      ELSEIF (ICALL .EQ. 5) THEN
         DO 150 I=1,N 
            W6(I) = YP(I) 
            Y(I) = W0(I) + H*(A71*W1(I) + A73*W3(I) + A74*W4(I)
     &                                  + A75*W5(I) + A76*W6(I))
  150    CONTINUE
         ICALL = 6
         RETURN
C
      ELSEIF (ICALL .EQ. 6) THEN
         DO 160 I=1,N
            W4(I) = H*(E1*W1(I) + E3*W3(I) + E4*W4(I) + E5*W5(I)
     &                                     + E6*W6(I) + E7*YP(I))
  160    CONTINUE
      ENDIF
C
C.....Error estimation
C
      ERR = ZER
      IF (ITOL .EQ. 0) THEN
         DO 200 I=1,N 
            SK  = ATOLI + RTOLI*MAX(ABS(W0(I)),ABS(Y(I)))
            TMP = W4(I)/SK
            ERR = ERR + TMP*TMP
  200    CONTINUE
      ELSE
         DO 210 I=1,N 
            SK  = ATOL(I) + RTOL(I)*MAX(ABS(W0(I)),ABS(Y(I)))
            TMP = W4(I)/SK
            ERR = ERR + TMP*TMP
  210    CONTINUE
      ENDIF
      ERR = SQRT(ERR/DBLE(N))
C
C.....Estimation of the next step including Lund stabilization
C.....It is required that FAC1 <= HNEW/H <= FAC2
C
      ERRFAC = ERR**EXPO1
      FAC    = ERRFAC/(FACOLD**BETA)
      FAC    = MAX(FACIN2, MIN(FACIN1,FAC/SAFE))
      HNEW   = H/FAC
C
C.....Error test
C
      IF (ERR .LE. ONE) THEN
C
C........Accepted step
C
         FACOLD = MAX(ERR,FACMIN)
         NACCPT = NACCPT+1
C
C........Set the next stepsize
C
         IF (ABS(HNEW) .GT. ABS(HMAX) 
     &         .AND. ABS(HMAX) .GT. ZER)HNEW = POSNEG*HMAX  
         IF(REJECT)HNEW = POSNEG*MIN(ABS(HNEW),ABS(H))
         H = HNEW
C
C........Successful return
C
         TASK = 'DONE'
         RETURN
C
      ELSE
C
C........Rejected step. Try smaller H
C
         HNEW   = H/MIN(FACIN1,ERRFAC/SAFE)
         REJECT = .TRUE.
         IF (NACCPT .GE. 1) NREJCT = NREJCT+1   
         H = HNEW
         GOTO 50
C
      ENDIF
C
C.....End of DOPSTA
C
      END
C 
C234567==1=========2=========3=========4=========5=========6=========7==
C
      SUBROUTINE DOPSTN( TASK,N,T,Y,YP,H,HMIN,HMAX,NMAX,
     &                   ATOL,RTOL,ITOL,W0,W1,W2,W3,W4,W5,W6,
     &                   JPOL,UINT,NSTEP,NACCPT,NREJCT )
C
      CHARACTER*6 TASK
      INTEGER ITOL,JPOL,N,NACCPT,NMAX,NREJCT,NSTEP 
      DOUBLE PRECISION T,Y(*),YP(*),H,HMIN,HMAX,RTOL(*),ATOL(*)
      DOUBLE PRECISION W0(*),W1(*),W2(*),W3(*),W4(*),W5(*),W6(*)
      DOUBLE PRECISION UINT(5,*)
C 
C234567--x---------x---------x---------x---------x---------x---------x--
C
C  Subroutine for taking a step with the Dormand-Prince Runge 
C  Kutta method of order 5 in solving the NONAUTONOMOUS initial  
C  value problem
C
C        y' = f(t,y),  y(t0) = y0,
C
C  where
C        t    the time
C        y    vector of dimension N,
C        F    mapping from R^N into R^N
C
C  The routine uses reverse communication and returns for all
C  evaluations of the user function for F. This is done under the
C  control of the character variable TASK.
C
C  To take a step compute YP = F(T,Y) for the given point (T,Y) and
C  call the routine with TASK = 'STEP'.
C  The routine returns with TASK = 'EVAL' signifying a request for
C  the evaluation of F at (T,Y). 
C  Return to the routine with the computed vector YP and without  
C  changing TASK or any other variable in the calling sequence.
C  Alternately, without computing YP, a forced step reduction may be
C  requested by returning to the routine with TASK = 'REDUCE' and
C  with a suggested smaller stepsize in H. As default, when H is 
C  unchanged the step size is reduced by a factor 1/2. This
C  feature is included to allow for the solution of systems 
C  defined in terms of local coordinates on a manifold.
C
C  Upon final return TASK may have either one of the following 
C  values:
C    TASK = 'MINSTP'  step fell below minimum steplength HMIN
C    TASK = 'STPCNT'  maximal number of steps NMAX exceeded
C    TASK = 'ERROR'   error condition.
C
C  Variables in the calling sequence:
C  ----------------------------------
C  TASK   C  IN   Task identifier
C                 TASK = 'STEP'   take a new RK step
C                 TASK = 'EVAL'   YP = F(Y,T) for given Y has been
C                                 computed 
C                 TASK = 'REDUCE' force a step reduction
C            OUT  TASK = 'EVAL'   Compute YP = F(Y,T) for given Y
C                 TASK = 'DONE'   RK step successfully completed
C                 TASK = 'MINSTP' Step fell below allowed minimum
C                 TASK = 'STPCNT' Maximal number of steps exceeded
C                 TASK = 'ERROR'  Other error condition 
C  N      I  IN   Dimension of Y
C  T      D  IN   T-value at the current point
C            OUT  T-value at the next point.
C  Y      D  IN   Array of dimension N, the current point
C            OUT  The next computed point
C  YP     D  IN   Array of dimension N, function value at Y 
C            OUT  Function value at the next computed point 
C  H      D  IN   Given stepsize
C         D  OUT  Predicted next stepsize
C  HMIN   D  IN   Minimal step
C  HMAX   D  IN   Maximal stepsize, HMAX > 0
C  NMAX   I  IN   Maximal number of steps
C  RTOL   D  IN   Relative error tolerance. May be
C                 an array of dimension 1 or N
C  ATOL   D  IN   Absolute error tolerance. May be
C                 an array of dimension 1 or N
C  ITOL   I  IN   Switch for RTOL and ATOL
C                 ITOL = 0: Both RTOL and ATOL are scalars
C                           The code keeps the local error of Y(I),
C                           roughly, below RTOL(1)*ABS(Y(I)) + ATOL(1) 
C                 ITOL = 1: Both RTOL and ATOL are vectors
C                           The code keeps the local error of Y(I),
C                           roughly, below RTOL(I)*ABS(Y(I)) + ATOL(I) 
C  W0-W6  D  WK   Seven work arrays of dimension N
C  JPOL   I  IN   Interpolation indicator
C                 JPOL = 0  No interpolation is requested
C                 JPOL = 1  Interpolation for continuous output is
C                           requested
C  UINT   D  OUT  Interpolation array of dimension 5 x N. 
C                 When JPOL = 0 then UINT is not referenced
C  NSTEP  I  I/O  Step counter
C  NACCPT I  I/O  Counter of accepted steps
C  NREJCT I  I/O  Counter of rejected steps (after the first step)
C
C234567--x---------x---------x---------x---------x---------x---------x--
C
C.....Step parameters
C
      DOUBLE PRECISION BETA,FAC1,FAC2,SAFE,EXPO1,FACMIN,FACIN1,FACIN2 
      PARAMETER( BETA=0.04D0, FAC1=0.2D0, FAC2=10.D0, SAFE=0.9D0,
     &           EXPO1=0.2D0 - BETA*0.75D0, FACMIN=1.0D-4, 
     &           FACIN1 = 1.0D0/FAC1, FACIN2 = 1.0D0/FAC2 )
C
C.....Dormand-Prince RK-5 coefficients for the non-autonomous case
C.....with interpolation
C
      DOUBLE PRECISION A21,A31,A32,A41,A42,A43,A51,A52,A53,A54
      DOUBLE PRECISION A61,A62,A63,A64,A65,A71,A73,A74,A75,A76
      DOUBLE PRECISION C2,C3,C4,C5,E1,E3,E4,E5,E6,E7
      DOUBLE PRECISION D21,D23,D24,D25,D26,D31,D33,D34,D35,D36
      DOUBLE PRECISION D41,D43,D44,D45,D46
C
      PARAMETER( A21 = 0.2D0,A31 = 3.0D0/40.0D0,
     &        A32 = 9.0D0/40.0D0, A41 = 44.0D0/45.0D0,
     &        A42 = -56.0D0/15.0D0, A43 = 32.0D0/9.0D0,
     &        A51 = 19372.0D0/6561.0D0, A52 = -25360.0D0/2187.0D0,
     &        A53 = 64448.0D0/6561.0D0, A54 = -212.0D0/729.0D0,
     &        A61 = 9017.0D0/3168.0D0, A62 = -355.0D0/33.0D0,
     &        A63 = 46732.0D0/5247.0D0, A64 = 49.0D0/176.0D0,
     &        A65 = -5103.0D0/18656.0D0, A71 = 35.0D0/384.0D0,
     &        A73 = 500.0D0/1113.0D0, A74 = 125.0D0/192.0D0,
     &        A75 = -2187.0D0/6784.0D0,A76 = 11.0D0/84.0D0 )
      PARAMETER( C2 = 0.2D0,C3 = 0.3D0,C4 = 0.8D0,C5 = 8.0D0/9.0D0)
      PARAMETER( E1  = 71.0D0/57600.0D0, E3  = -71.0D0/16695.0D0,
     &        E4  = 71.0D0/1920.0D0, E5  = -17253.0D0/339200.0D0,
     &        E6  = 22.0D0/525.0D0,E7  = -1.0D0/40.0D0 )
      PARAMETER( D21=-1337.0D0/480.0D0,D23=4216.0D0/1113.0D0,
     &        D24=-135.0D0/80.0D0,D25=-2187.0D0/8480.0D0,
     &        D26=66.0D0/70.0D0,D31=1039.0D0/360.0D0,
     &        D33=-468200.0D0/83475.0D0,D34=9.0D0/2.0D0,
     &        D35=400950.0D0/318000.0D0,D36=-638.0D0/210.0D0,
     &        D41=-1163.0D0/1152.0D0,D43=37900.0D0/16695.0D0,
     &        D44=-415.0D0/192.0D0,D45=-674325.0D0/508800.0D0,
     &        D46=374.0D0/168.0D0 )
C
C.....Other parameters
C
      DOUBLE PRECISION ZER,HALF,ONE
      PARAMETER( ZER=0.0D0, HALF=0.5D0, ONE=1.0D0 )
C 
C.....Fortran functions called
C
      DOUBLE PRECISION ABS,MAX,MIN,SQRT
C
C.....Local variables
C
      INTEGER I
      DOUBLE PRECISION ERR,ERRFAC,FAC,HNEW,SK,TMP
C
C.....Variables to be saved between calls
C
      LOGICAL REJECT
      INTEGER ICALL
      DOUBLE PRECISION ATOLI,FACOLD,HAB,POSNEG,RTOLI,T0
      SAVE REJECT,ICALL,ATOLI,FACOLD,HAB,POSNEG,RTOLI,T0
C
      DATA FACOLD/1.0D-4/
C
C.......................Executable statements...........................
C
C.....Check TASK
C
      IF (TASK .EQ. 'EVAL') THEN
         GOTO 100
C
      ELSEIF (TASK .EQ. 'STEP') THEN
C
C........Check data
C
         IF (ITOL .EQ. 0) THEN
            ATOLI  = ATOL(1)
            RTOLI  = RTOL(1)
         ENDIF
         REJECT = .FALSE.
         POSNEG = ONE
         HAB    = H
         IF (H .LT. ZER) THEN
            POSNEG = -POSNEG
            HAB = -HAB
         ENDIF
C
C........Store initial T, Y and YP
C
         T0 = T
         DO 10 I = 1, N
            W0(I) = Y(I)
            W1(I) = YP(I)
   10    CONTINUE
         GOTO 50
C
      ELSEIF (TASK .EQ. 'REDUCE') THEN
C
C........Forced step reduction
C
         IF (ABS(H) .GE. HAB) H = HALF*H
         REJECT = .TRUE.
         IF(NACCPT .GE. 1)NREJCT = NREJCT + 1
C
      ELSE
C
C........TASK has an improper value
C
         TASK = 'ERROR'
         RETURN
      ENDIF
C
C.....Restart point for rejected steps
C
   50 CONTINUE
C
C.....Check for maximal number of steps
C
      NSTEP = NSTEP + 1
      IF (NSTEP .GT. NMAX) THEN
         TASK = 'STPCNT'
         RETURN
      ENDIF
C
C.....Check for minimal stepsize
C
      IF (ABS(H) .LT. HMIN) THEN
         TASK = 'MINSTP'
         RETURN
      ENDIF
      T = T0
      TASK = 'EVAL'
C
C.....Stage Evaluations
C
      DO 90 I = 1,N
         Y(I)  = W0(I) + H*A21*W1(I)
   90 CONTINUE
      T = T0 + C2*H 
      ICALL = 1
      RETURN
C
  100 CONTINUE
      IF (ICALL .EQ. 1) THEN
         DO 110 I=1,N
            W2(I) = YP(I) 
            Y(I)  = W0(I) + H*(A31*W1(I) + A32*W2(I))
  110    CONTINUE
         T = T0 + C3*H
         ICALL = 2
         RETURN
C
      ELSEIF (ICALL .EQ. 2) THEN
         DO 120 I=1,N 
            W3(I) = YP(I) 
            Y(I)  = W0(I) + H*(A41*W1(I) + A42*W2(I) + A43*W3(I))
  120    CONTINUE
         T = T0 + C4*H
         ICALL = 3
         RETURN
C
      ELSEIF (ICALL .EQ. 3) THEN
         DO 130 I=1,N 
            W4(I) = YP(I) 
            Y(I) = W0(I) + H*(A51*W1(I) + A52*W2(I) + A53*W3(I)
     &                                              + A54*W4(I))
  130    CONTINUE
         T = T0 + C5*H
         ICALL = 4
         RETURN
C
      ELSEIF (ICALL .EQ. 4) THEN
         DO 140 I=1,N 
            W5(I) = YP(I) 
            Y(I) = W0(I) + H*(A61*W1(I) + A62*W2(I) + A63*W3(I)
     &                                  + A64*W4(I) + A65*W5(I))
  140    CONTINUE
         T = T0 + H
         ICALL = 5
         RETURN
C
      ELSEIF (ICALL .EQ. 5) THEN
         DO 150 I=1,N 
            W6(I) = YP(I) 
            Y(I) = W0(I) + H*(A71*W1(I) + A73*W3(I) + A74*W4(I)
     &                                  + A75*W5(I) + A76*W6(I))
  150    CONTINUE
C
C........If continuous output is requested calculate UINT0 - UINT4
C
         IF (JPOL .EQ. 1) THEN 
            DO 160 I = 1,N
               UINT(1,I) = W0(I)
               UINT(2,I) = W1(I)
               UINT(3,I) = D21*W1(I) + D23*W3(I) + D24*W4(I)
     &                               + D25*W5(I) + D26*W6(I)
               UINT(4,I) = D31*W1(I) + D33*W3(I) + D34*W4(I)
     &                               + D35*W5(I) + D36*W6(I)
               UINT(5,I) = D41*W1(I) + D43*W3(I) + D44*W4(I)
     &                               + D45*W5(I) + D46*W6(I)
  160       CONTINUE
         ENDIF
         ICALL = 6
         RETURN
C
      ELSEIF (ICALL .EQ. 6) THEN
         DO 170 I=1,N
            W4(I) = H*(E1*W1(I) + E3*W3(I) + E4*W4(I) + E5*W5(I)
     &                                     + E6*W6(I) + E7*YP(I))
  170    CONTINUE
      ENDIF
C
C.....Error estimation
C
      ERR = ZER
      IF (ITOL .EQ. 0) THEN
         DO 200 I=1,N 
            SK  = ATOLI + RTOLI*MAX(ABS(W0(I)),ABS(Y(I)))
            TMP = W4(I)/SK
            ERR = ERR + TMP*TMP
  200    CONTINUE
      ELSE
         DO 210 I=1,N 
            SK  = ATOL(I) + RTOL(I)*MAX(ABS(W0(I)),ABS(Y(I)))
            TMP = W4(I)/SK
            ERR = ERR + TMP*TMP
  210    CONTINUE
      ENDIF
      ERR = SQRT(ERR/DBLE(N))
C
C.....Estimation of the next step including Lund stabilization
C.....It is required that FAC1 <= HNEW/H <= FAC2
C
      ERRFAC = ERR**EXPO1
      FAC    = ERRFAC/(FACOLD**BETA)
      FAC    = MAX(FACIN2, MIN(FACIN1,FAC/SAFE))
      HNEW   = H/FAC
C
C.....Error test
C
      IF (ERR .LE. ONE) THEN
C
C........Accepted step
C
         FACOLD = MAX(ERR,FACMIN)
         NACCPT = NACCPT+1
C
C........Set the next stepsize
C
         IF (ABS(HNEW) .GT. ABS(HMAX) 
     &         .AND. ABS(HMAX) .GT. ZER)HNEW = POSNEG*HMAX  
         IF(REJECT)HNEW = POSNEG*MIN(ABS(HNEW),ABS(H))
         H = HNEW
C
C........Successful return
C
         TASK = 'DONE'
         RETURN
C
      ELSE
C
C........Rejected step. Try smaller H
C
         HNEW   = H/MIN(FACIN1,ERRFAC/SAFE)
         REJECT = .TRUE.
         IF (NACCPT .GE. 1) NREJCT = NREJCT+1   
         H = HNEW
         GOTO 50
C
      ENDIF
C
C.....End of DOPSTN
C
      END
C 
C234567==1=========2=========3=========4=========5=========6=========7==
C
      SUBROUTINE BIDIA( M,N,A,LDA,D,E,Q,LDQ,WRK,SAFMIN,IER )
C
      INTEGER IER,LDA,LDQ,M,N
      DOUBLE PRECISION A(LDA,*),D(*),E(*),Q(LDQ,*),WRK(*),SAFMIN
C
C234567--x---------x---------x---------x---------x---------x---------x--
C
C  Reduces a general M by N matrix A with M >= N to an M by N,
C  upper bidiagonal-matrix B by a transformation: Q' * A * P = B 
C  with an M x M orthogonal matrix Q and an N x N orthogonal 
C  matrix P. The matrix P is returned in the first N rows of the 
C  array A and the matrix Q' in the array Q.
C
C  Variables in the calling sequence
C  ---------------------------------
C  M      I   IN   The number of rows in the matrix A.  M >= 0.
C  N      I   IN   The number of columns in the matrix A.  N >= 0.
C  A      D   IN   The given M by N matrix, M >= N
C             OUT  The first N rows contain the N x N orthogonal
C                  matrix Q.
C  LDA    I   IN   The leading dimension of A. LDA >= max(1,M)
C  D      D   OUT  Array of dimension N, the diagonal elements 
C                  of the bidiagonal matrix B:  D(i) = A(i,i).
C  E      D   OUT  Array of dimension N-1, the off-diagonal 
C                  elements of the bidiagonal matrix B:
C                  E(i) = A(i,i+1) for i = 1,2,...,n-1;
C  Q      D   OUT  The transpose of the left orthogonal matrix.
C  LDQ    D   IN   The leading dimension of Q, LDQ >= max(1,M)
C  WRK    D   OUT  Work array of dimension N
C  SAFMIN D   IN   Safe minimum such that 1.0/SAFMIN does not
C                  overflow
C  IER    I   OUT  Error indicator
C                  IER = 0: successful exit
C                  IER = 1: input data error for BIDIA
C                  IER = 2: input data error from HOUSL or HOUSR
C
C234567--x---------x---------x---------x---------x---------x---------x--
C
C.....External subroutines
C
      EXTERNAL HOUSG,HOUSL,HOUSR
C
C.....Parameters
C
      DOUBLE PRECISION ZER, ONE
      PARAMETER ( ZER = 0.0D+0, ONE = 1.0D+0 )
C
C.....Local variables
C
      INTEGER I,IP1,J,K,MI
      DOUBLE PRECISION TAUI
C
C.......................Executable statements...........................
C
C.....Test for input errors
C
      IF( (M.LT.0) .OR. (N.LT.0) .OR. (N.GT.M) 
     &             .OR. (LDA.LT.1) .OR. (LDA.LT.M) 
     &             .OR. (LDQ.LT.1) .OR. (LDQ.LT.M) ) THEN
         CALL MSGPRT('BIDIA',' Error in the input-dimensions')
         IER = 1
         RETURN
      ENDIF
C
C.....Quick return if possible
C
      IER = 0
      IF( (M.EQ.0) .OR. (N.EQ.0) ) RETURN
C
C.....Initialize Q
C
      DO 20 J = 1, M
         DO 10 I = 1, M
            Q(I,J) = ZER
   10    CONTINUE
         Q(J,J) = ONE
   20 CONTINUE
C
C.....Reduction to upper bidiagonal form
C
      DO 30 I = 1, N
         IP1 = I + 1
         MI  = M - I + 1
C
C........Generate elementary reflector HQ(i) to annihilate A(i+1:m,i)
C
         TAUI = ZER
         IF( I .LT. M )CALL HOUSG( MI,A(I,I),A(IP1,I),1,TAUI,SAFMIN )
C
         D(I) = A(I,I)
C
C........Apply HQ(i) to A(i:m,i+1:n) from the left
C
         CALL HOUSL( MI,N-I,A(I,I),1,TAUI,A(I,IP1),LDA,IER )
         IF( IER .NE. 0) THEN
            CALL MSGPRT('BIDIA',' Input-dimension error in HOUSL') 
            IER = 2
            RETURN
         ENDIF
C
C........Apply HQ(I) to Q(1:M,I:N) from the right
C
         CALL HOUSR( M,MI,A(I,I),1,TAUI,Q(1,I),LDQ,IER )
         IF( IER .NE. 0) THEN
            CALL MSGPRT('BIDIA',' Input-dimension error in HOUSR') 
            IER = 2
            RETURN
         ENDIF
C
         IF( I .LT. N ) THEN
C
C...........Generate elementary reflector HP(i) to annihilate
C...........A(i,i+2:n)
C
            K = I + 2
            CALL HOUSG(N-I,A(I,IP1),A(I,K),LDA,WRK(I),SAFMIN)
            E(I) = A(I,IP1)
C
C...........Apply HP(i) to A(i+1:m,i+1:n) from the right
C
            CALL HOUSR(M-I,N-I,A(I,IP1),LDA,WRK(I),A(IP1,IP1),LDA,IER)
            IF( IER .NE. 0) THEN
               CALL MSGPRT('BIDIA',' Input-dimension error in HOUSR') 
               IER = 2
               RETURN
            ENDIF
         ELSE
            WRK(I) = ZER
         ENDIF
   30 CONTINUE
C
C.....Initialize the n-th row of A
C
      DO 40 J = 1, N
         A(N,J) = ZER
   40 CONTINUE
      A(N,N) = ONE
C
C.....Multiply A*HP(N-1)*...*HP(1)
C
      DO 60 I = N-1,1,-1
         IP1 = I + 1
C
C........Apply HP(I) to A(IP1:N,IP1:N) from the right
C
         CALL HOUSR(N-I,N-I,A(I,IP1),LDA,WRK(I),A(IP1,IP1),LDA,IER)
         IF( IER .NE. 0) THEN
            CALL MSGPRT('BIDIA',' Input-dimension error in HOUSR') 
            IER = 2
            RETURN
         ENDIF
C
C........Initialize the i-th row
C
         DO 50 J = 1, N
            A(I,J) = ZER
   50    CONTINUE
         A(I,I) = ONE
   60 CONTINUE
C
      RETURN
C
C.....End of BIDIA
C
      END
C
C234567==1=========2=========3=========4=========5=========6=========7==
C
      DOUBLE PRECISION FUNCTION DDIST2( N,X,INCX,Y,INCY,KVEC )
C
      INTEGER N, INCX, INCY, KVEC
      DOUBLE PRECISION X(*), Y(*)
C
C234567--1---------2---------3---------4---------5---------6---------7--
C
C  Computes either the Euclidean distance between two N-dimensional 
C  vectors X and Y or the Euclidean norm of one such vector X. 
C
C  Call the routine 
C  either
C     with KVEC = 2 and two vectors X and Y of dimension N stored
C     with storage-increments INCX and INCY, respectively.
C  or
C     with KVEC = 1 and one vector X of dimension N stored
C     with storage-increment INCX. In this case the second array 
C     is not referenced and can be a dummy array or simply the
C     array X again.
C
C  If N .LE. 0 then zero is returned, if N .GE. 1 then the
C  storage increments cannot be zero.
C
C  The algorithm follows the four-phase method of C. L. Lawson in the
C  LAPACK routine DNRM2.F.  As in DNRM2.F two built-in constants are 
C  used that are hopefully applicable to all machines.
C     CUTLO = maximum of  DSQRT(u/eps)  over all known machines.
C     CUTHI = minimum of  DSQRT(v)      over all known machines.
C  where
C     eps = smallest number such that 1.0d0 + eps .gt. 1.0d0
C     u   = smallest positive number  (underflow limit)
C     v   = largest  number           (overflow  limit)
C
C  Values for CUTLO and CUTHI listed in DNRM2.F are as follows:
C
C     CUTLO, s.p.  u/eps = 2**(-102) for honeywell. Close 
C                  seconds are univac and dec at 2**(-103)
C                  thus CUTLO = 2**(-51) = 4.44089e-16
C     CUTHI, s.p.  v = 2**127 for univac, honeywell, and dec.
C                  thus CUTHI = 2**(63.5) = 1.30438e19
C     CUTLO, d.p.  u/eps = 2**(-67) for honeywell and dec.
C                  thus CUTLO = 2**(-33.5) = 8.23181d-11
C     CUTHI, d.p.  same as s.p.  CUTHI = 1.30438d19
C
C     data cutlo, cuthi / 8.232d-11,  1.304d19 /
C     data cutlo, cuthi / 4.441e-16,  1.304e19 /                 
C
C  In line with the four phases of DNRM2.F, the algorithm uses 
C  four states identified by LEVEL = 0,1,2,3, respectively, 
C  which correspond to the following cases:
C
C     LEVEL = 0  only zero terms have been found so far
C     LEVEL = 1  all nonzero terms encountered so far do not
C                exceed CUTLO in modulus
C     LEVEL = 2  there are some terms that are larger than
C                CUTLO in modulus but none exceeds 
C                HITEST = CUTLO/DBLE(N) 
C     LEVEL = 3  there are terms that exceed HITEST in modulus.
C
C  All state transitions can only increase the LEVEL.
C
C  Variables in the calling sequence
C  ---------------------------------
C     N     I    IN   Dimension of the vectors X and Y
C     X     D    IN   The first vector of dimension N
C     INCX  I    IN   Storage increment of X
C     Y     D    IN   The second vector of dimension N
C     INYY  D    IN   Storage increment of Y
C     KVEC  I    IN   Number of vectors
C                     KVEC = 1  Only one vector, namely X, is given,
C                               the Euclidean norm of X is computed  
C                               and the Y array is not referenced
C                               This is the default
C                     KVEC = 2  Two vectors X and Y  are given,
C                               the Euclidean distance between 
C                               X and Y is computed
C
C234567--1---------2---------3---------4---------5---------6---------7--
C
C.....FORTRAN functions called
C
      DOUBLE PRECISION ABS, SQRT, DBLE
C
C.....Parameters and data
C
      DOUBLE PRECISION ZER, ONE, CUTLO, CUTHI
      PARAMETER( ZER = 0.0D0, ONE = 1.0D0, CUTLO = 8.232D-11, 
     &          CUTHI = 1.304D19 )
C
C.....Local variables
C
      INTEGER NN, IX, IY, J, LEVEL, JOB
      DOUBLE PRECISION DIFF, DX, DY, HITEST, SUM, TMP, TRM, XMAX
C
C.......................Executable statements...........................
C
C.....Return DDIST2 = ZER if N is not positive or if INCX = 0
C
      DDIST2 = ZER
      IF (N .LE. 0 .OR. INCX .EQ. 0) RETURN
C
      JOB = 1
      IF ( KVEC .EQ. 2 ) THEN
C
C........There are two vectors -- return DDIST2 = ZER if INCY = 0
C
         JOB = 2
         IF (INCY .EQ. 0) RETURN
      ENDIF
C
C.....Initializations
C
      NN = N
      HITEST = CUTHI/DBLE(NN)
      IF (HITEST .LT. CUTLO) HITEST = CUTLO
C
      IX = 1
      IF (INCX .LT. 0) IX = (-NN+1)*INCX + 1
      IY = 1
      IF (INCY .LT. 0) IY = (-NN+1)*INCY + 1
C
      SUM = ZER
      XMAX = ZER
      LEVEL = 0
C
C.....Loop over all N terms
C
      DO 10 J = 1,NN
         DX = X(IX)
         IX = IX + INCX
         IF (JOB .EQ. 1) THEN
C
C...........One vector
C
            DIFF = DX
         ELSE
C
C...........Two vectors
C
            DY = Y(IY)
            IY = IY + INCY
            IF (SIGN(ONE,DX) .NE. SIGN(ONE,DY)) THEN
               DIFF = DX - DY
            ELSE
               IF (ABS(DX) .LT. ABS(DY)) THEN
                  TMP = DY
                  DY  = DX
                  DX  = TMP
               ENDIF
               IF (DX .EQ. ZER) THEN
                  DIFF = ZER
               ELSE
                  DIFF = DX*(ONE - DY/DX)
               ENDIF
            ENDIF
         ENDIF
C
C........Summation of the squares of the nonzero terms
C
         IF (DIFF .NE. ZER) THEN
            TRM = ABS(DIFF) 
C
            IF (TRM .LE. CUTLO) THEN
C
C..............Small nonzero terms -- transition to level 1
C
               IF (LEVEL .EQ. 0) THEN
                  LEVEL = 1
                  XMAX = TRM
                  TMP = TRM/XMAX
                  SUM = SUM + TMP*TMP
               ELSE 
                  IF (LEVEL .EQ. 1) THEN
                     IF( TRM .GT. XMAX ) THEN
                        TMP = XMAX/TRM
                        SUM = ONE + SUM*TMP*TMP
                        XMAX = TRM
                     ELSE
                        TMP = TRM/XMAX
                        SUM = SUM + TMP*TMP
                     ENDIF
                  ELSE           
                     SUM = SUM + TRM*TRM
                  ENDIF
               ENDIF
C
            ELSE
C
C..............Mid-sized terms -- transition to level 2
C
               IF (LEVEL .EQ. 0) THEN
                  LEVEL = 2
               ELSE IF (LEVEL .EQ. 1) THEN
                  LEVEL = 2
                  SUM = (SUM * XMAX) * XMAX
               ENDIF
C
               IF (TRM .LE. HITEST) THEN
                  SUM = SUM + TRM*TRM
               ELSE
C
C.................Large terms
C
                  IF (LEVEL .EQ. 2) THEN
C
C                    Transition to level = 3
C
                     LEVEL = 3
                     SUM = (SUM / TRM) / TRM
                     XMAX = TRM
                  ENDIF
C
                  IF( TRM .GT. XMAX ) THEN
                     TMP = XMAX/TRM
                     SUM = ONE + SUM*TMP*TMP
                     XMAX = TRM
                  ELSE
                     TMP = TRM/XMAX
                     SUM = SUM + TMP*TMP
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
C
   10 CONTINUE
C
      IF (XMAX .EQ. ZER .OR. LEVEL.EQ.2) THEN
         DDIST2 = SQRT(SUM)
      ELSE
         DDIST2 = XMAX *SQRT(SUM)
      ENDIF
C
      RETURN
C
C.....End if DDIST2
C
      END
C
C234567==1=========2=========3=========4=========5=========6=========7==
C
      SUBROUTINE HOUSG( N,ALPHA,X,INCX,TAU,SAFMIN )
C
      INTEGER N,INCX
      DOUBLE PRECISION ALPHA,X(*),TAU,SAFMIN
C
C234567--x---------x---------x---------x---------x---------x---------x--
C
C  Generates an n-dimensional Householder reflector
C
C     H = I - tau*( 1 ) * ( 1 v' ),    H' * H = I,   tau scalar
C                 ( v ) 
C
C  such that
C
C     H * ( alpha ) = ( beta ),    alpha, beta scalars
C         (   x   )   (   0  )     x  (n-1)-dimensional vector
C
C  Because of H'* H = I it follows that
C 
C     alpha^2 + x^T x = beta^2    ==> beta = sqrt(alpha^2 + x^T x)
C
C     H'( beta )  = ( alpha )     ==> tau = (beta - alpha)/beta
C       (  0   )    (  x    )           v = xscal*x,  
C                                   xscal = 1/(beta - alpha)
C       
C     If x = 0, then tau = 0 and H = I, otherwise  1 <= tau <= 2.
C
C  This is an edited version of the LAPACK routine DLARFG
C
C  Variables in the calling sequence:
C  ----------------------------------
C  N      I   IN   Dimension of H
C  ALPHA  D   IN   The scalar alpha
C             OUT  The scalar beta
C  X      D   IN   The given vector x of dimension n - 1
C             OUT  The vector v
C  INCX   I   IN   The increment between elements of X, INCX .NE. 0
C  TAU    D   OUT  The scalar tau
C  SAFMIN D   IN   Safe minimum such that 1.0/SAFMIN does not
C                  overflow
C
C234567--x---------x---------x---------x---------x---------x---------x--
C
C.....Functions called
C
      DOUBLE PRECISION ABS, DDIST2, SQRT, DNRM2
C
C.....Parameters
C
      DOUBLE PRECISION ONE, ZER
      PARAMETER( ONE = 1.0D+0, ZER = 0.0D+0 )
C
C.....Local variables
C
      INTEGER I,IX,J,KNT,KX,NM1
      DOUBLE PRECISION A1,A2,BETA,RSAFMN,TMP,XNORM,XSCAL
C
C.......................Executable statements...........................
C
      IF( (N.LE.1) .OR. (INCX.EQ.0) )THEN
         TAU = ZER
         RETURN
      ENDIF
C
      NM1 = N - 1
      XNORM = DNRM2( NM1,X,INCX)
C      XNORM = DDIST2( NM1,X,INCX,X,INCX,1)
C
      IF( XNORM .EQ. ZER ) THEN
C
C........H is the identity
C
         TAU = ZER
      ELSE
C
C........General case
C
         KX = 1
         IF(INCX .LT. 0)KX = 1 - (N-1)*INCX
         A1 = XNORM
         A2 = ABS(ALPHA)
         IF ( A1 .LT. A2 ) THEN
            A1 = A2
            A2 = XNORM
         ENDIF
         IF ( A2 .EQ. ZER ) THEN
            BETA = A1
         ELSE
            TMP = A1/A2
            BETA = A2*SQRT( ONE + TMP*TMP )
         ENDIF
         IF (ALPHA .GT. ZER) BETA = -BETA
         IF( ABS(BETA) .GE. SAFMIN ) THEN
            TAU = (BETA - ALPHA) / BETA
            XSCAL = ONE / (ALPHA - BETA)
            ALPHA = BETA
         ELSE
C
C.......... XNORM, BETA may be inaccurate; scale X and recompute
C
            RSAFMN = ONE / SAFMIN
            KNT = 0
   10       CONTINUE
               KNT = KNT + 1
               IF( INCX .EQ. 1) THEN
                  DO 20 I = 1,NM1
                     X(I) = RSAFMN*X(I)
   20             CONTINUE
               ELSE
                  IX = KX
                  DO 30 I = 1,NM1
                     X(IX) = RSAFMN*X(IX)
                     IX = IX + INCX
   30             CONTINUE
               ENDIF
               BETA = BETA*RSAFMN
               ALPHA = ALPHA*RSAFMN
            IF( ABS( BETA ) .LT. SAFMIN )GOTO 10
C
C...........New BETA satisfies SAFMIN <= BETA <= 1.0
C
             XNORM = DNRM2( NM1,X,INCX)
C            XNORM = DDIST2( NM1, X, INCX,X, INCX,1)
            A1 = XNORM
            A2 = ABS(ALPHA)
            IF ( A1 .LT. A2 ) THEN
               A1 = A2
               A2 = XNORM
            ENDIF
            IF ( A2 .EQ. ZER ) THEN
               BETA = A1
            ELSE
               TMP = A1/A2
               BETA = A2*SQRT( ONE + TMP*TMP)
            ENDIF
            IF (ALPHA .GT. ZER) BETA = -BETA
            TAU = (BETA - ALPHA) / BETA
            XSCAL = ONE / (ALPHA - BETA)
            ALPHA = BETA
            DO 40 J = 1, KNT
               ALPHA = ALPHA*SAFMIN
   40       CONTINUE
         ENDIF
         IF( INCX .EQ. 1) THEN
            DO 50 I = 1,NM1
               X(I) = XSCAL*X(I)
   50       CONTINUE
         ELSE
            IX = KX
            DO 60 I = 1,NM1
               X(IX) = XSCAL*X(IX)
               IX = IX + INCX
   60       CONTINUE
         ENDIF
      ENDIF
C
      RETURN
C
C.....End of HOUSG
C
      END
C
C234567==1=========2=========3=========4=========5=========6=========7==
C
      SUBROUTINE HOUSL( M,N,X,INCX,TAU,A,LDA,IER )
C
      INTEGER M,N,INCX,LDA,IER
      DOUBLE PRECISION X(*),TAU,A(LDA,*)
C
C234567--x---------x---------x---------x---------x---------x---------x--
C
C  Multiplies a given M x N matrix A from the (L)eft by an 
C  M-dimensional Householder reflector
C
C     H = I - tau*x*x',   H' * H = I,  x = ( 1 )    
C                                          ( v ) 
C
C  that is, overwrite A by the product H * A
C
C  Variables in the calling sequence:
C  ----------------------------------
C  M     I   IN   The number of rows of A
C  N     I   IN   The number of columns of A
C  X     D   IN   The given vector x. x(1) = 1.0 is enforced
C  INCX  I   IN   The increment between elements of X, INCX .NE. 0
C  TAU   D   OUT  The scalar tau
C  A     D   IN   The given matrix of dimension M x N
C            OUT  The computed matrix product H * A
C  LDA   I   IN   The leading dimension of the array A, LDA >= M
C  IER   I   OUT  Error indicator
C                 IER = 0  no error
C                 IER = 1  input-data error
C
C234567--x---------x---------x---------x---------x---------x---------x--
C
C.....Parameter
C
      DOUBLE PRECISION ZER
      PARAMETER( ZER = 0.0D+0 )
C
C.....Local variables
C
      INTEGER I,IX,J,KX
      DOUBLE PRECISION SUM,TMP
C
C.......................Executable statements...........................
C
C.....Test input data
C
      IF( (M.LT.0) .OR. (N.LT.0) .OR. (INCX.EQ.0) .OR. (LDA.LT.M) )THEN
         CALL MSGPRT('HOUSL',' Error in the input-dimensions')
         IER = 1 
         RETURN
      ENDIF
C
C.....Quick return without error message
C
      IER = 0
      IF( (M.EQ.0) .OR. (N.EQ.0) .OR. (TAU.EQ.ZER) )RETURN
C
      IF( INCX .EQ. 1 )THEN
C
C........Form w = A' * x and  A := A - tau * x * w'
C
         DO 30 J = 1, N
            SUM = A(1,J)
            IF( M .GT. 1)THEN
               DO 10 I = 2, M
                  SUM = SUM + A(I,J)*X(I)
   10          CONTINUE
            ENDIF
            IF( SUM .NE. ZER )THEN
               TMP = -TAU*SUM
               A(1,J) = A(1,J) + TMP
               IF( M .GT. 1)THEN
                 DO 20 I = 2, M
                    A(I,J) = A(I,J) + X(I)*TMP
   20            CONTINUE
               ENDIF
            ENDIF
   30    CONTINUE
      ELSE
         KX = 1
         IF( INCX .LT. 0 )KX = 1 - (M-1)*INCX
C
C........Form w = A' * x and  A := A - tau * x * w'
C
         DO 60 J = 1, N
            SUM = A(1,J)
            IF( M .GT. 1)THEN
               IX  = KX
               DO 40 I = 2, M
                 IX  = IX + INCX
                 SUM = SUM + A(I,J)*X(IX)
   40          CONTINUE
            ENDIF
            IF( SUM .NE. ZER )THEN
               TMP = -TAU*SUM
               A(1,J) = A(1,J) + TMP
               IF( M .GT. 1)THEN
                  IX = KX
                  DO 50 I = 2,M
                     IX     = IX + INCX
                     A(I,J) = A(I,J) + X(IX)*TMP
   50             CONTINUE
               ENDIF
            ENDIF
   60    CONTINUE
      ENDIF
C
      RETURN
C
C.....End of HOUSL
C
      END
C
C234567==1=========2=========3=========4=========5=========6=========7==
C
      SUBROUTINE HOUSR( M,N,X,INCX,TAU,A,LDA,IER)
C
      INTEGER M,N,INCX,LDA,IER
      DOUBLE PRECISION X(*),TAU,A(LDA,*)
C
C234567--x---------x---------x---------x---------x---------x---------x--
C
C  Multiplies a given M x N matrix A from the (R)ight by an 
C  N-dimensional Householder reflector
C
C     H = I - tau * x * x',   H' * H = I,  x = ( 1 )    
C                                              ( v ) 
C
C  that is, overwrite A by the product A * H
C
C  Variables in the calling sequence:
C  ----------------------------------
C  M     I   IN   The number of rows of A
C  N     I   IN   The number of columns of A
C  X     D   IN   The given vector x. x(1) = 1.0 is enforced
C  INCV  I   IN   The increment between elements of V, INCV .NE. 0
C  TAU   D   OUT  The scalar tau
C  A     D   IN   The given matrix of dimension M x N
C            OUT  The computed matrix product A * H
C  LDA   I   IN   The leading dimension of the array A, LDA >= M
C  IER   I   OUT  Error indicator
C                 IER = 0  no error
C                 IER = 1  input-data error
C
C234567--x---------x---------x---------x---------x---------x---------x--
C
C.....Parameters
C
      DOUBLE PRECISION ZER
      PARAMETER( ZER = 0.0D+0 )
C
C.....Local variables
C
      INTEGER I,J,JX,KX
      DOUBLE PRECISION SUM,TMP
C
C.......................Executable statements...........................
C
C.....Test input data
C
      IF( (M.LT.0) .OR. (N.LT.0) .OR. (INCX.EQ.0) .OR. (LDA.LT.M) )THEN
         CALL MSGPRT('HOUSR',' Error in the input-dimensions')
         IER = 1 
         RETURN
      ENDIF
C
C.....Quick return without error message
C
      IER = 0
      IF( (M.EQ.0) .OR. (N.EQ.0) .OR. (TAU.EQ.ZER) )RETURN
C
      IF( INCX .EQ. 1 )THEN
C
C........Form w = A * x and A := A - tau * w * x'
C
         DO 30 I = 1, M
            SUM = A(I,1)
            IF( N .GT. 1)THEN
               DO 10 J = 2, N
                  SUM = SUM + A(I,J)*X(J)
   10          CONTINUE
            ENDIF
            IF( SUM .NE. ZER )THEN
               TMP = -TAU*SUM
               A(I,1) = A(I,1) + TMP
               IF( N .GT. 1)THEN
                  DO 20 J = 2, N
                     A(I,J) = A(I,J) + TMP*X(J)
   20             CONTINUE
               ENDIF
            ENDIF
   30    CONTINUE
C
      ELSE
C
         KX = 1
         IF( INCX .LE. 0 )KX = 1 - (N-1)*INCX
C
C........Form w = A * x and A := A - tau * w * x'
C
         DO 60 I = 1, M
            SUM = A(I,1)
            IF( N .GT. 1)THEN
               JX  = KX
               DO 40 J = 2, N
                  JX  = JX  + INCX
                  SUM = SUM + A(I,J)*X(JX)
   40          CONTINUE
            ENDIF
            IF( SUM .NE. ZER )THEN
               TMP = -TAU*SUM
               A(I,1) = A(I,1) + TMP
               IF( N .GT. 1)THEN
                  JX  = KX
                  DO 50 J = 2, N
                     JX  = JX  + INCX
                     A(I,J) = A(I,J) + TMP*X(JX)
   50             CONTINUE
               ENDIF
            ENDIF
   60    CONTINUE
      ENDIF
C
      RETURN
C
C.....End of HOUSL
C
      END
C
C234567==1=========2=========3=========4=========5=========6=========7==
C
      SUBROUTINE LQF( M,N,A,LDA,IPIV,TAU,WRK,SAFMIN,IER )
C
      INTEGER M,N,IPIV(*),LDA,IER
      DOUBLE PRECISION A(LDA,*),TAU(*),WRK(*),SAFMIN
C
C234567--x---------x---------x---------x---------x---------x---------x--
C
C  Computes the LQ factorization  P*A = L*Q  with column pivoting 
C  on the columns of an M by N matrix A under the assumption that 
C  M <= N. The matrix Q is represented as a product of elementary 
C  reflectors Q = H(1) H(2) . . . H(M), H(i) = I - tau * v * v'. 
C  Here tau is a real scalar, and v is a real vector with v(1:i-1) = 0 
C  and v(i) = 1. Upon return, the components v(i+1:n) are stored in 
C  A(i,i+1:n) and tau in the array TAU(i). The matrix P is stored in 
C  the array IPIV as follows: If JPVT(j) = i then the jth row of P is 
C  the ith canonical unit vector.
C
C  Variables in the calling sequence
C  ---------------------------------
C  M      I   IN   The number of rows of the matrix A. 0 <= M <= N.
C  N      I   IN   The number of columns of the matrix A.  N >= 0.
C  A      D   IN   The given M by N matrix
C             OUT  The elements on and below the diagonal contain
C                  the M by M lower triangular matrix L; the elements
C                  above the diagonal, together with the array TAU, 
C                  represent the orthogonal matrix Q as a product 
C                  of M elementary reflectors
C  LDA    I   IN   The leading dimension of A. LDA >= max(1,M).
C  IPIV   I   OUT  If IPIV(i) = k, then the i-th row of P*A was the 
C                  k-th row of A.
C  TAU    I   OUT  Array of dimension M, The scalar factors of 
C                  the elementary reflectors
C  WRK    D   WK   Work array of  dimension M
C  SAFMIN D   IN   Safe minimum such that 1.0/SAFMIN does not
C                  overflow
C  IER    I   OUT  Error indicator
C                  IER = 0   successful exit
C                  IER = 1   data-error in LQF
C                  IER = 2   data error reported from HOUSR
C
C234567--x---------x---------x---------x---------x---------x---------x--
C
C.....External subroutines
C
      EXTERNAL HOUSG,HOUSR
C
C.....Functions used
C
      DOUBLE PRECISION SQRT,DDIST2
C
C.....Parameters
C
      DOUBLE PRECISION FACT,ONE,ZER
      PARAMETER ( FACT = 0.05D+0, ONE = 1.0D+0, ZER = 0.0D+0 )
C
C.....Local variables
C
      INTEGER I,IMAX,IP1,ITMP,J
      DOUBLE PRECISION CNRM,CNRMJ,TAUI,TMP1,TMP2,TMP3
C
C.......................Executable statements...........................
C
C.....Test for input errors
C
      IF( (N.LT.0) .OR. (N.LT.M) .OR. (M.LT.0) 
     &             .OR. (LDA.LE.0) .OR. (LDA.LT.M) ) THEN
         CALL MSGPRT('LQF',' Error in the input-dimensions') 
         IER = 1
         RETURN
      ENDIF
C
C.....Quick return if possible
C
      IER = 0
      IF( (M .EQ.0) .OR. (N.EQ.0) )RETURN
C
C.....Initialize column norms and pivot array
C
      DO 10 I = 1, M
         CNRM = DDIST2(N,A(I,1),LDA,A(I,1),LDA,1)     
         TAU(I)  = CNRM
         WRK(I)  = CNRM
         IPIV(I) = I
   10 CONTINUE
C
C.....Main loop for the factorization
C
      DO 50 I = 1, M
         IP1 = I + 1
C
C........Determine ith pivot row
C
         IMAX = I
         CNRM = TAU(I)
         IF(I .LT. M) THEN
            DO 20 J = IP1,M
               IF (TAU(J) .GT. CNRM) THEN
                  IMAX = J
                  CNRM = TAU(J)
               ENDIF
   20       CONTINUE
C
C...........Swap the rows if necessary
C
            IF( IMAX .NE. I ) THEN
               DO 30 J = 1, N
                  TMP1      = A(IMAX,J)
                  A(IMAX,J) = A(I,J)
                  A(I,J)    = TMP1
   30          CONTINUE
               TAU(IMAX) = TAU(I)
               WRK(IMAX) = WRK(I)
               ITMP       = IPIV(IMAX)
               IPIV(IMAX) = IPIV(I)
               IPIV(I)    = ITMP
            ENDIF
         ENDIF
         IF( CNRM .NE. ZER) THEN
C
C...........Generate elementary reflector H(i)
C
            TAUI = ZER
            IF( I .LT. N )
     &         CALL HOUSG(N-I+1,A(I,I),A(I,IP1),LDA,TAUI,SAFMIN)
            IF( I .LT. M ) THEN
C
C..............Apply H(i) to A(i:n,i+1:m) from the right
C
               CALL HOUSR( M-I,N-I+1,A(I,I),LDA,TAUI,A(IP1,I),LDA,IER )
               IF( IER .NE. 0) THEN
                  CALL MSGPRT('LQF',' Input-dimension error in HOUSR') 
                  IER = 2
                  RETURN
               ENDIF
C
C..............Update row norms
C
               DO 40 J = IP1, M
                  CNRMJ = TAU(J)
                  IF( CNRMJ .NE. ZER ) THEN
                     TMP1 = ABS(A(J,I))/CNRMJ
                     TMP2 = ONE - TMP1*TMP1
                     IF( TMP2 .LE. ZER) THEN
                        TMP2 = ZER
                        TMP3 = ZER
                     ELSE
                        TMP3 = SQRT(TMP2)
                     ENDIF
                     TMP1 = CNRMJ/WRK(J)
                     TMP2 = ONE + FACT*TMP2*TMP1*TMP1
                     IF( TMP2 .EQ. ONE ) THEN
                        TAU(J) = DDIST2(N-I,A(J,IP1),LDA,A(J,IP1),LDA,1)
                        WRK(J) = TAU(J)
                     ELSE
                        TAU(J) = TAU(J)*TMP3
                     ENDIF
                  ENDIF
   40          CONTINUE
            ENDIF
            TAU(I) = TAUI
         ENDIF
   50 CONTINUE
C
      RETURN
C
C.....End of LQF
C
      END
C
C234567==1=========2=========3=========4=========5=========6=========7==
C
      SUBROUTINE LQAS( M,N,A,LDA,IPIV,TAU,Y,X,IER )
C
      INTEGER LDA,N,M,IPIV(*),IER
      DOUBLE PRECISION A(LDA,*),TAU(*),Y(*),X(*)
C
C234567--x---------x---------x---------x---------x---------x---------x--
C
C   For an  M x N matrix A with M <= N and rank A = M and a given
C   M-vector y, compute the unique solution x of A x = y which
C   is orthogonal to ker A. The algorithm 
C
C     1.  Solve L Y := P * Y
C     2.  X := (Y)
C              (0)
C     2.  X := Q^T*X
C
C   is used under the assumption that the arrays A, TAU, and 
C   IPIV contain the LQ-factorization of the matrix A.
C
C   Variables in the calling sequence
C   ---------------------------------
C   M    I   IN   Number of rows of the matrix A, M <= N
C   N    I   IN   Number of columns of the matrix A
C   A    D   IN   Array of dimension LDA x n, the factored matrix
C   LDA  I   IN   Leading dimension of A, LDA >= M
C   IPIV I   IN   Array of dimension M, the permutation of
C                 the LQ-factorization as returned by LQF
C   TAU  D   IN   Array of dimension M containing the scalar 
C                 factors of the elementary reflectors, as 
C                 returned by LQF
C   Y    D   IN   Array of dimension M, the given vector
C   X    D   OUT  Array of dimension N, the computed vector
C   IER  I   OUT  Error indicator
C                 IER = 1  input data error
C                 IER < 0  zero pivot encountered
C                          IER = -J signifies that the
C                          J-th diagonal element of the
C                          triangular matrix R is zero.
C
C234567--x---------x---------x---------x---------x---------x---------x--
C
C.....Parameter
C
      DOUBLE PRECISION ZER
      PARAMETER( ZER=0.0D0 )
C
C.....Local variables
C
      INTEGER I,IM1,IP,J,JP1
      DOUBLE PRECISION SUM,T
C
C.......................Executable statements...........................
C
C.....Test for data-input errors
C
      IF( (M.LT.0) .OR. (N.LT.0) .OR. (M.GT.N) .OR. (LDA.LT.M) )THEN
         CALL MSGPRT('LQAS',' Error in the input-dimensions') 
         IER = 1
         RETURN
      ENDIF
C
C.....Quick return if possible
C
      IER = 0
      IF( (M.EQ.0) .OR. (N.EQ.0) )RETURN
C
C.....Take care of the case N = 1
C
      IF (N .EQ. 1) THEN
         IF(A(1,1) .EQ. ZER) THEN
            CALL MSGPRT('LQAS',' Zero pivot encountered') 
            IER = -1
         ELSE
            X(1) = Y(1)/A(1,1)
         ENDIF
         RETURN
      ENDIF
C
C.....Solve L * Z = P * Y for Z and extend to an M-vector
C
      DO 20 I = 1,N
         IF (I .LE. M) THEN
            IF(A(I,I) .EQ. ZER) THEN
               CALL MSGPRT('LQAS',' Zero pivot encountered') 
               IER = -I
               RETURN
            ENDIF
            IP = IPIV(I)
            SUM = Y(IP)
            IF(I .GT. 1) THEN
               IM1 = I-1
               DO 10 J = 1, IM1
                  SUM = SUM - A(I,J)*X(J)
   10          CONTINUE
            ENDIF
            X(I) = SUM/A(I,I)
         ELSE
            X(I) = ZER
         ENDIF
   20 CONTINUE
C
C.....Form Q^T*X
C
      DO 50 J = M, 1, -1
         JP1 = J + 1
         SUM = X(J)
         IF( J .LT. N ) THEN
            DO 30 I = JP1, N
               SUM = SUM + A(J,I)*X(I)
   30       CONTINUE
         ENDIF
         IF( SUM .NE. ZER ) THEN
            T = - TAU(J)*SUM
            X(J) = X(J) + T
            IF( J .LT. N ) THEN
               DO 40 I = JP1, N
                  X(I) = X(I) + T*A(J,I)
   40          CONTINUE            
            ENDIF
         ENDIF
   50 CONTINUE
C
      RETURN
C
C.....End of LQAS
C
      END
C
C234567==1=========2=========3=========4=========5=========6=========7==
C
      SUBROUTINE LUF( N,A,LDA,IPIV,IER )
C
      INTEGER IER,IPIV(*),LDA,N
      DOUBLE PRECISION A(LDA,*)
C
C234567--x---------x---------x---------x---------x---------x---------x--
C
C  Computes an LU factorization of a general N-by-N matrix A
C  using partial pivoting with row interchanges:
C
C                   A = P * L * U
C
C  where P is a permutation matrix, L is lower triangular with unit
C  diagonal elements, and U is upper triangular.
C
C  Variables in the calling sequence
C  ---------------------------------
C  N    I   IN   The dimension of the square matrix A.  N >= 0.
C  A    D   IN   Array of dimension (LDA,N), the given matrix A
C           OUT  The factors L and U from the factorization
C                A = P*L*U; the unit diagonal elements of L are 
C                not stored.
C  LDA  I   IN   The leading dimension of the array A. LDA >= max(1,N)
C  IPIV I   OUT  Array of dimension N
C                The pivot indices; for 1 <= I <= N, row i of
C                the matrix was interchanged with row IPIV(I).
C  IER  I   OUT  Error indicator
C                IER = 0:  successful exit
C                IER = 1: input-data error
C                IER < 0:  if |IER| = K, U(K,K) is exactly zero. The 
C                          factorization has been completed, but the 
C                          factor U is exactly singular, and division 
C                          by zero will occur if it is used to solve a 
C                          system of equations.
C
C234567--x---------x---------x---------x---------x---------x---------x--
C
C.....Functions called
C
      DOUBLE PRECISION ABS
C
C.....Parameters
C
      DOUBLE PRECISION ONE, ZER
      PARAMETER ( ONE = 1.0D+0, ZER = 0.0D+0 )
C
C.....Local variables
C
      INTEGER I,J,K,KP1,KPIV
      DOUBLE PRECISION APIV,T
C
C.......................Executable statements...........................
C
C.....Test the input data
C
      IF( (N.LT.0) .OR. (LDA.LT.N)  ) THEN
         CALL MSGPRT('LUF',' Error in the input-dimensions') 
         IER = 1
         RETURN
      ENDIF
C
C.....Quick return if possible
C
      IER = 0
      IF( N .EQ. 0 )RETURN
C
C.....Loop over the columns
C
      DO 100 K = 1, N
         KP1 = K + 1
C
C........Find pivot index
C
         KPIV = K
         APIV = ABS(A(K,K))
         IF( K .LT. N) THEN
            DO 20 I = KP1,N
               IF (ABS(A(I,K)) .GT. APIV) THEN
                  KPIV = I
                  APIV = ABS(A(I,K))
               ENDIF
   20       CONTINUE
         ENDIF
         IPIV(K) = KPIV
         IF( APIV .EQ. ZER ) THEN
            CALL MSGPRT('LUF',' Zero pivot encountered') 
            IER = -K
         ELSE
C
C...........Apply the interchange to columns 1:N.
C
            IF( KPIV .NE. K )THEN
               DO 30 I = 1, N
                  T         = A(KPIV,I)
                  A(KPIV,I) = A(K,I)
                  A(K,I)    = T
   30          CONTINUE
            ENDIF
C
C...........Compute the multipliers
C
            IF( K .LT. N )THEN
               T = ONE/A(K,K)
               DO 40 I = KP1,N
                  A(I,K) = A(I,K)*T
   40          CONTINUE
C
C..............Update trailing submatrix
C
               DO 60 J = KP1,N
                  T = -A(K,J)
                  DO 50 I = KP1,N
                     A(I,J) = A(I,J) + A(I,K)*T
   50             CONTINUE
   60          CONTINUE
            ENDIF
         ENDIF
  100 CONTINUE
C
      RETURN
C
C.....End of LUF
C
      END
C
C234567==1=========2=========3=========4=========5=========6=========7==
C
      SUBROUTINE LUS1( N,A,LDA,IPIV,X,IER )
C
      INTEGER IER,IPIV(*),LDA,N
      DOUBLE PRECISION A(LDA,*),X(*)
C
C234567--x---------x---------x---------x---------x---------x---------x--
C
C  Solves a system of linear equations A * Z = X for an N-by-N matrix 
C  A and an N-vector X and returns the result in X. It is assumed
C  that the arrays A and IPIV contain the LU factorization returned 
C  by LUF
C
C  Variables in the calling sequence
C  ---------------------------------
C  N    I   IN   The dimension of the square matrix A.  N >= 0.
C  A    D   IN   Array of dimension (LDA,N), the factors L and U 
C                from the factorization A = P*L*U as computed by DLUF
C  LDA  I   IN   The leading dimension of the array A. LDA >= max(1,N)
C  IPIV I   IN   Array of dimension N, the pivot indices from DLUF; 
C                for 1<=i<=N, row i of the matrix was interchanged 
C                with row IPIV(i).
C  X    D   IN   Array of dimension N, the right hand side vector 
C           OUT  The solution vector X.
C  IER  I   OUT  Error indicator
C                IER = 0:  successful exit
C                IER = 1: input-data error
C                IER < 0:  if |IER| = K, the K-th pivot is zero.
C
C234567--x---------x---------x---------x---------x---------x---------x--
C
C.....Parameters
C
      DOUBLE PRECISION ZER
      PARAMETER ( ZER = 0.0D0)
C
C.....Local variables
C
      INTEGER I,IP,J,NM1
      DOUBLE PRECISION TMP
C
C.......................Executable statements...........................
C
C.....Test the input data
C
      IF( (N.LT.0) .OR. (LDA.LT.N)  ) THEN
         CALL MSGPRT('LUS1',' Error in the input-dimensions') 
         IER = 1
         RETURN
      ENDIF
C
C.....Quick return if possible
C
      IER = 0
      IF( N .EQ. 0 )RETURN
C
C.....TAKE CARE OF THE CASE N = 1
C
      IF( N .EQ. 1) THEN
         IF( A(1,1) .EQ. ZER )THEN
            CALL MSGPRT('LUS1',' Zero pivot encountered') 
            IER = -1
         ELSE         
            X(1) = X(1)/A(1,1)
         ENDIF
         RETURN
      ENDIF
      NM1 = N - 1
C
C.....Apply row interchanges to the right side.
C
      DO 10 I = 1, N
         IP = IPIV(I)
         IF( IP .NE. I )THEN
            TMP  = X(IP)
            X(IP) = X(I)
            X(I) = TMP
         ENDIF
   10 CONTINUE
C
C.....Solve L*Z = X, overwriting B with Z.
C
      DO 30 J = 1, NM1
         TMP = X(J)
         IF( TMP .NE. ZER )THEN
            DO 20 I = J+1, N
               X(I) = X(I) - TMP*A(I,J)
   20       CONTINUE
         ENDIF
   30 CONTINUE
C
C.....Solve U*Z = X, overwriting X with Z.
C
      DO 50 J = N, 1, -1
         IF( X(J) .NE. ZER )THEN
            IF( A(J,J) .EQ. ZER )THEN
               CALL MSGPRT('LUS1',' Zero pivot encountered') 
               IER = -J
            ELSE         
               X(J) = X(J)/A(J,J)
               IF( J .GT. 1) THEN
                  DO 40 I = 1, J - 1
                     X(I) = X(I) - X(J)*A(I,J)
   40             CONTINUE
               ENDIF
            ENDIF
         ENDIF
   50 CONTINUE
C
      RETURN
C
C.....End of LUS1
C
      END
C
C234567==1=========2=========3=========4=========5=========6=========7==
C
      SUBROUTINE LUSK( N,A,LDA,IPIV,B,LDB,K,IER )
C
      INTEGER IER,IPIV(*),LDA,LDB,N,K
      DOUBLE PRECISION A(LDA,*),B(LDB,*)
C
C234567--x---------x---------x---------x---------x---------x---------x--
C
C  Solves a system of linear equations A * X = B for an N-by-N matrix 
C  A and an N x K matrix B on the right side. It is assumed that the
C  array A contains the LU factorization of the system's matrix 
C  returned by LUF
C
C  Variables in the calling sequence
C  ---------------------------------
C  N    I   IN   The dimension of the square matrix A.  N >= 0.
C  A    D   IN   Array of dimension (LDA,N), the factors L and U 
C                from the factorization A = P*L*U as computed by DLUF
C  LDA  I   IN   The leading dimension of the array A. LDA >= max(1,N)
C  IPIV I   IN   Array of dimension N, the pivot indices from DLUF; 
C                for 1<=i<=N, row i of the matrix was interchanged 
C                with row IPIV(i).
C  B    D   IN   Array of dimension (LDB,NRHS), the right hand side 
C                matrix B.
C           OUT  The solution matrix X.
C  LDB  I   IN   The leading dimension of the array B. LDB >= max(1,N)
C  K    I   IN   The number of right hand sides, i.e., the number 
C                of columns of the matrix B.  K >= 1.
C  IER  I   OUT  Error indicator
C                IER = 0:  successful exit
C                IER = 1: input-data error
C                IER < 0:  if |IER| = K, the K-th pivot is zero.
C
C234567--x---------x---------x---------x---------x---------x---------x--
C
C.....Parameters
C
      DOUBLE PRECISION ZER
      PARAMETER ( ZER = 0.0D0)
C
C.....Local variables
C
      INTEGER I,IP,J,L,NM1
      DOUBLE PRECISION TMP
C
C.......................Executable statements...........................
C
C.....Test the input data
C
      IF( (N.LT.0) .OR. (K.LT.0) 
     &             .OR. (LDA.LT.N) .OR. (LDB.LT.N) ) THEN
         CALL MSGPRT('LUSK',' Error in the input-dimensions') 
         IER = 1
         RETURN
      ENDIF
C
C.....Quick return if possible
C
      IER = 0
      IF( (N.EQ.0) .OR. (K.EQ.0) )RETURN
C
C.....TAKE CARE OF THE CASE N = 1
C
      IF( N .EQ. 1) THEN
         IF( A(1,1) .EQ. ZER )THEN
            CALL MSGPRT('LUSK',' Zero pivot encountered') 
            IER = -1
         ELSE
            DO 10 L = 1, K         
               B(1,L) = B(1,L)/A(1,1)
   10       CONTINUE
         ENDIF
         RETURN
      ENDIF
      NM1 = N - 1
C
C.....Apply row interchanges to the right hand sides.
C
      DO 30 I = 1, N
         IP = IPIV( I )
         IF( IP .NE. I )THEN
            DO 20 L = 1, K
               TMP  = B(IP,L)
               B(IP,L) = B(I,L)
               B(I,L) = TMP
   20       CONTINUE
         ENDIF
   30 CONTINUE
C
C.....Solve L*X = B, overwriting B with X.
C
      DO 60, L = 1, K
         DO 50 I = 1,NM1
            IF( B(I,L) .NE. ZER )THEN
               DO 40 J = I+1, N
                  B(J,L) = B(J,L) - B(I,L)*A(J,I)
   40          CONTINUE
            ENDIF
   50    CONTINUE
   60 CONTINUE
C
C.....Solve U*X = B, overwriting B with the solution.
C
      DO 90 L = 1, K
         DO 80 I = N, 1, -1
            IF( B(I,L) .NE. ZER )THEN
               IF( A(I,I) .EQ. ZER )THEN
                  CALL MSGPRT('LUSK',' Zero pivot encountered') 
                  IER = -I
                  RETURN
               ELSE
                  B(I,L) = B(I,L)/A(I,I)
               ENDIF
               IF( I .GT. 1) THEN
                  DO 70 J = 1, I-1
                     B(J,L) = B(J,L) - B(I,L)*A(J,I)
   70             CONTINUE
               ENDIF
            ENDIF
   80    CONTINUE
   90 CONTINUE
C
      RETURN
C
C.....End of LUSK
C
      END
C
C234567==1=========2=========3=========4=========5=========6=========7==
C
      SUBROUTINE QRF( M,N,A,LDA,IPIV,TAU,WRK,SAFMIN,IER )
C
      INTEGER M,N,IPIV(*),LDA,IER
      DOUBLE PRECISION A(LDA,*),TAU(*),WRK(*),SAFMIN
C
C234567--x---------x---------x---------x---------x---------x---------x--
C
C  Compute the QR factorization with column pivoting on all columns 
C  of an M by N matrix A
C                      A*P = Q*R
C
C  The matrix Q is represented as a product of elementary reflectors
C
C     Q = H(1) H(2) . . . H(k),   where k = min(m,n).
C
C  Each H(i) has the form
C
C     H(i) = I - tau * v * v'
C
C  where tau is a real scalar, and v is a real vector with v(1:i-1) = 0
C  and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i), 
C  and tau in TAU(i).
C
C  The matrix P is represented in IPIV as follows: If IPIV(j) = i
C  then the jth column of P is the ith canonical unit vector.
C
C  Variables in the calling sequence
C  ---------------------------------
C  M      I   IN   The number of rows of the matrix A.  M >= 0.
C  N      I   IN   The number of columns of the matrix A.  N >= 0.
C  A      D   IN   The given M by N matrix
C             OUT  The elements on and above the diagonal contain
C                  the min(m,n) by n upper trapezoidal matrix R 
C                  (R is upper triangular if m >= n); the elements 
C                   belowthe diagonal, together with the array TAU, 
C                   represent the orthogonal matrix Q as a product 
C                   of min(m,n) elementary reflectors
C  LDA    I   IN   The leading dimension of A. LDA >= max(1,M).
C  IPIV   I   OUT  If IPIV(i) = k, then the i-th column of A*P was the 
C                   k-th column of A.
C  TAU    I   OUT  Array of dimension min(M,N), The scalar factors of 
C                  the elementary reflectors
C  WRK    D   WK   Work array of  dimension N
C  SAFMIN D   IN   Safe minimum such that 1.0/SAFMIN does not
C                  overflow
C  IER    I   OUT  Error indicator
C                  IER = 0   successful exit
C                  IER = 1   data-error
C                  IER = 2   input data error in HOUSL
C                  IER < 0   zero pivot encountered
C
C234567--x---------x---------x---------x---------x---------x---------x--
C
C.....External subroutines
C
      EXTERNAL HOUSG,HOUSL
C
C.....Functions used
C
      DOUBLE PRECISION SQRT,DDIST2
C
C.....Parameters
C
      DOUBLE PRECISION FACT,ONE,ZER
      PARAMETER ( FACT = 0.05D+0, ONE = 1.0D+0, ZER = 0.0D+0 )
C
C.....Local variables
C
      INTEGER I,IMAX,IP1,ITMP,J,MN
      DOUBLE PRECISION CNRM,CNRMJ,TAUI,TMP1,TMP2,TMP3
C
C.......................Executable statements...........................
C
C.....Test for input errors
C
      IF( (M.LT.0) .OR. (N.LT.0) .OR. (LDA.LT.M) ) THEN
         CALL MSGPRT('QRF',' Error in the input-dimensions') 
         IER = 1
         RETURN
      ENDIF
C
C.....Quick return if possible
C
      IER = 0
      IF( (M.EQ.0) .OR. (N.EQ.0) )RETURN
C
C.....Initialize column norms and pivot array
C
      DO 10 I = 1, N
         CNRM = DDIST2(M,A(1,I),1,A(1,I),1,1)
         TAU(I)  = CNRM
         WRK(I)  = CNRM
         IPIV(I) = I
   10 CONTINUE
C
      MN = N
      IF( M .LT. N)MN = M
C
C.....Main loop for the factorization
C
      DO 50 I = 1, MN
         IP1 = I + 1
C
C........Determine pivot column
C
         CNRM = TAU(I)
         IMAX = I
         IF(I .LT. N) THEN
            DO 20 J = IP1,N
               IF (TAU(J) .GT. CNRM) THEN
                  IMAX = J
                  CNRM = TAU(J)
               ENDIF
   20       CONTINUE
            IF( CNRM .EQ. ZER) THEN
               CALL MSGPRT('QRF',' Zero pivot encountered') 
               IER = -I
            ENDIF
C
C...........Swap the columns if necessary
C
            IF( IMAX .NE. I ) THEN
               DO 30 J = 1, M
                  TMP1      = A(J,IMAX)
                  A(J,IMAX) = A(J,I)
                  A(J,I)    = TMP1
   30          CONTINUE
               ITMP       = IPIV(IMAX)
               IPIV(IMAX) = IPIV(I)
               IPIV(I)    = ITMP
               TAU(IMAX)  = TAU(I)
               WRK(IMAX)  = WRK(I)
            ENDIF
         ENDIF
         IF( CNRM .EQ. ZER) THEN
            CALL MSGPRT('QRF',' Zero pivot encountered') 
            IER = -I
         ELSE
C
C...........Generate elementary reflector H(i)
C
            TAUI = ZER
            IF( I .LT. M )
     &         CALL HOUSG(M-I+1,A(I,I),A(IP1,I),1,TAUI,SAFMIN)
            IF( I .LT. N ) THEN
C
C..............Apply H(i) to A(i:m,i+1:n) from the left
C
               CALL HOUSL( M-I+1,N-I,A(I,I),1,TAUI,A(I,IP1),LDA,IER )
               IF( IER .NE. 0) THEN
                  CALL MSGPRT('QRF',' Input-dimension error in HOUSL') 
                  IER = 2
                  RETURN
               ENDIF
C
C..............Update column norms
C
               DO 40 J = IP1, N
                  CNRMJ = TAU(J)
                  IF( CNRMJ .NE. ZER ) THEN
                     TMP1 = ABS(A(I,J))/CNRMJ
                     TMP2 = ONE - TMP1*TMP1
                     IF( TMP2 .LE. ZER) THEN
                        TMP2 = ZER
                        TMP3 = ZER
                     ELSE
                        TMP3 = SQRT(TMP2)
                     ENDIF
                     TMP1 = CNRMJ/WRK(J)
                     TMP2 = ONE + FACT*TMP2*TMP1*TMP1
                     IF( TMP2 .EQ. ONE ) THEN
                        TAU(J) = DDIST2(M-I,A(IP1,J),1,A(IP1,J),1,1)
                        WRK(J) = TAU(J)
                     ELSE
                        TAU(J) = TAU(J)*TMP3
                     ENDIF
                  ENDIF
   40          CONTINUE
            ENDIF
            TAU(I) = TAUI
         ENDIF
   50 CONTINUE
C
      RETURN
C
C.....End of QRF
C
      END
C
C234567==1=========2=========3=========4=========5=========6=========7==
C
      SUBROUTINE QORG( M,N,K,L,A,LDA,TAU,Q,LDQ,IER )
C
      INTEGER IER,K,L,LDA,LDQ,M,N
      DOUBLE PRECISION A(LDA,*),TAU(*),Q(LDQ,*)
C
C234567--x---------x---------x---------x---------x---------x---------x--
C
C  Generate an M by NQ real matrix Q with NQ = L - K + 1 orthonormal 
C  columns formed as the product of N Householder reflectors 
C  H(1),...,H(N) of order M >= N as returned in the columns of the 
C  array A by QRF. 
C
C  Variables in the calling sequence:
C  ----------------------------------
C  
C  M    I  IN   The number of rows of the matrix A. M >= 0.
C  N    I  IN   The number of columns of the matrix A. M >= N >= 0.
C  K    I  IN   The lower column index, M >= K >= 1
C               K is set to 1 if K <= 0 
C  L    I  IN   The higher column index, M >= L >= K >= 1
C               L is set to M if L >= M 
C  A    D  IN   Array of dimension (LDA,N)
C               For IA = 1 the i-th column of A is assumed to contain
C                   the vector defining the elementary reflector 
C                   H(i), for i = 1,2,...,N, as returned by QRF 
C               For IA = 2 the i-th row of A is assumed to contain
C                   the vector defining the elementary reflector 
C                   H(i), for i = 1,2,...,M, as returned by QRF 
C  LDA  I  IN   The leading dimension of the array A, LDA >= max(1,M).
C  TAU  D  IN   Array of dimension (K) where TAU(i) contains the 
C               scalar factor of the elementary reflector H(i), as 
C               returned by QRF.
C  Q    D  OUT  The orthogonal matrix.
C  LDQ  D  IN   The leading dimension of the array Q, LDQ >= max(1,M) 
C  IER  I  OUT  Error indicator
C               IER = 0  successful exit
C               IER = 1  input-data error
C               IER = 2  input-data error in HOUSL
C
C  NOTE: The arrays A and Q cannot be identified
C
C234567--x---------x---------x---------x---------x---------x---------x--
C
C.....External subroutines called
C
      EXTERNAL HOUSL
C
C.....Parameters
C
      DOUBLE PRECISION ONE, ZER
      PARAMETER ( ONE = 1.0D+0, ZER = 0.0D+0 )
C
C.....Local variables
C
      INTEGER I,J1,J2,KM1,NL,NQ
C
C.......................Executable statements...........................
C
C.....Test for data-input errors
C
      IF( (M.LT.0) .OR. (N.LT.0) .OR. (N.GT.M) .OR. (K.GT.L)
     &             .OR. (LDA.LT.1) .OR. (LDA.LT.M) 
     &             .OR. (LDQ.LT.1) .OR. (LDQ.LT.M) ) THEN
         CALL MSGPRT('QORG',' Error in the input-dimensions') 
         IER = 1
         RETURN
      ENDIF
C
C.....Quick return if possible
C
      IER = 0
      IF( (N.EQ.0) .OR. (K.EQ.0) )RETURN
C
C.....Set defaults
C
      IF( K .LT. 0) K = 1
      IF( L .GT. M) L = M
C
C.....Initialize Q
C
      KM1 = K - 1
      NQ  = L - KM1
      DO 20 J1 = 1,NQ
         DO 10 I = 1, M
            Q(I,J1) = ZER
   10    CONTINUE
         Q(J1+KM1,J1) = ONE
   20 CONTINUE
C
C.....Apply H(1)*...*H(NL)*Q
C
      NL = L
      IF( L .GT. N) NL = L
      J2 = NQ 
      DO 50 J1 = NL, 1, -1
C
C........Apply H(J1) to Q(J1:M,J1:NQ) from the left
C
         CALL HOUSL(M-J1+1,NQ-J2+1,A(J1,J1),1,TAU(J1),Q(J1,J2),LDQ,IER)
         IF( IER .NE. 0) THEN
            CALL MSGPRT('QORG',' Input-dimension error in HOUSL') 
            IER = 2
            RETURN
         ENDIF
         J2 = J2 - 1
         IF( J2 .LT. 1 )J2 = 1
   50 CONTINUE
C
      RETURN
C
C.....End of QORG
C
      END
C234567==1=========2=========3=========4=========5=========6=========7==
C
      SUBROUTINE QRS( M,N,A,LDA,TAU,Y,X,IER )
C
      INTEGER LDA,N,M,IER
      DOUBLE PRECISION A(LDA,*),TAU(*),X(*),Y(*)
C
C234567--x---------x---------x---------x---------x---------x---------x--
C
C   For an  M x N matrix A with M >= N and rank A = N, and a given
C   M-vector Y, compute the least squares solution
C
C      min { || A*X - Y || ; X in R^N }
C
C   by the algorithm  
C
C       1.  Y := Q^T Y
C       2.  IF( M >= N ) THEN Solve R*Z := (I, 0)Y, Set X := Z
C                        ELSE Solve R*Z := Y, Set X = (Z,0)
C
C   under the assumption that the arrays A and TAU contain the
C   QR-factorization of an M x N matrix A.
C
C   The routine returns Q^T Y in the array Y and the least squares
C   solution in the array X. The array X may be identified with Y 
C   if Q^T*Y is not needed.
C
C   Variables in the calling sequence
C   ---------------------------------
C   M    I   IN   Number of rows of the matrix A, M >= N
C   N    I   IN   Number of columns of the matrix A
C   A    D   IN   Array of dimension LDA x n, the factored matrix
C   LDA  I   IN   Leading dimension of A, LDA >= M
C   TAU  D   IN   Array of dimension N containing the scalar 
C                 factor of the elementary reflectors, as 
C                 returned by QRF
C   Y    D   IN   Array of dimension M, the given vector Y 
C            OUT  The vector Q^T Y
C   X    D   OUT  Array of dimension N, the computed solution
C   IER  I   OUT  Error indicator
C                 IER = 0   no error
C                 IER = 1  input data error
C                 IER < 0  zero pivot encountered
C                          IER = -J signifies that the
C                          J-th diagonal element of the
C                          triangular matrix R is zero.
C
C234567--x---------x---------x---------x---------x---------x---------x--
C
C.....Parameter
C
      DOUBLE PRECISION ZER
      PARAMETER( ZER=0.0D0 )
C
C.....Local variables
C
      INTEGER I,J,JM1,JP1
      DOUBLE PRECISION SUM,T
C
C.......................Executable statements...........................
C
C.....Test for data-input errors
C
      IF( (M.LT.0) .OR. (N.LT.0) .OR. (M.LT.N)   
     &             .OR. (LDA.LT.1).OR. (LDA.LT.M) ) THEN
         CALL MSGPRT('QRS',' Error in the input-dimensions') 
         IER = 1
         RETURN
      ENDIF
C
C.....Quick returns if possible
C
      IER = 0
      IF( (M.EQ.0) .OR. (N.EQ.0) )RETURN      
C
C.....Take care of the case M = 1
C
      IF (M .EQ. 1) THEN
         IF( A(1,1) .EQ. ZER ) THEN
            CALL MSGPRT('QRS',' Zero pivot encountered') 
            IER = -1
         ELSE
            X(1) = Y(1)/A(1,1)
         ENDIF
         RETURN
      ENDIF
C
C.....Compute Q^T*Y
C
      DO 30 J = 1, N
         JP1 = J + 1
         SUM = Y(J)
         IF( J .LT. M ) THEN
            DO 10 I = JP1, M
               SUM = SUM + A(I,J)*Y(I)
   10       CONTINUE
         ENDIF
         IF( SUM .NE. ZER ) THEN
            T = - TAU(J)*SUM
            Y(J) = Y(J) + T
            IF( J .LT. M ) THEN
               DO 20 I = JP1, M
                  Y(I) = Y(I) + T*A(I,J)
   20          CONTINUE            
            ENDIF
         ENDIF
   30 CONTINUE
C
      DO 40 J = 1, N
         X(J) = Y(J)
   40 CONTINUE
C
C.....Solve R*Z := (I,0)Y, and set X := Z
C.....or solve  R*Z := Y, and set X := (Z,0)
C
      DO 60 J = N, 1, -1
         IF (A(J,J) .EQ. ZER) THEN
            CALL MSGPRT('QRS',' Zero pivot encountered') 
            IER = -J
            RETURN
         ENDIF
         X(J) = X(J)/A(J,J)
         IF (J .GT. 1) THEN
            JM1 = J - 1
            T = -X(J)
            DO 50 I = 1, JM1
               X(I) = X(I) + T*A(I,J)
   50       CONTINUE
         ENDIF
   60 CONTINUE
C
      RETURN
C
C.....End of QRS
C
      END
C
C234567==x=========x=========x=========x=========x=========x=========x==
C
      SUBROUTINE ROTG( F,G,CS,SN,R )
C
      DOUBLE PRECISION  CS,F,G,R,SN
C
C234567--x---------x---------x---------x---------x---------x---------x--
C
C  Generate a plane rotation so that
C
C     [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.
C     [ -SN  CS  ]     [ G ]     [ 0 ]
C
C  If G = 0, then CS = 1 and SN = 0.
C  If F = 0 and (G .ne. 0), then CS = 0 and SN = 1
C  F and G will remain unchanged
C
C  Variables in the calling sequence:
C  ----------------------------------
C  F    D    IN    The first component of vector to be rotated
C  G    D    IN    The second component of vector to be rotated
C  CS   D    OUT   The cosine term of the rotation
C  SN   D    OUT   The sine term of the rotation
C  R    D    OUT   The nonzero component of the rotated vector 
C
C234567--x---------x---------x---------x---------x---------x---------x--
C
C.....Intrinsic functions called
C
      DOUBLE PRECISION ABS, SQRT
C
C.....Local variables
C
      DOUBLE PRECISION T,TT
C
C.....Parameters
C
      DOUBLE PRECISION  ONE, ZER
      PARAMETER ( ONE = 1.0D0, ZER = 0.0D0 )
C
C.......................Executable statements...........................
C
      IF( G .EQ. ZER ) THEN
         CS = ONE
         SN = ZER
         R  = F
      ELSEIF( F .EQ. ZER ) THEN
         CS = ZER
         SN = ONE
         R = G
      ELSE
         IF( ABS(F) .GT. ABS(G) ) THEN
            T = G / F
            TT = SQRT( ONE+T*T )
            CS = ONE / TT
            SN = T*CS
            R  = F*TT
         ELSE
            T = F / G
            TT = SQRT( ONE+T*T )
            SN = ONE / TT
            CS = T*SN
            R = G*TT
         ENDIF
      ENDIF
C
      RETURN
C
      END
C
C234567==1=========2=========3=========4=========5=========6=========7==
C
      SUBROUTINE ROTML( M, N, A, LDA, J1, J2, C, S, IER )
C
      INTEGER M,N,LDA,J1,J2,IER
      DOUBLE PRECISION A(LDA,*),C,S
C
C234567--x---------x---------x---------x---------x---------x---------x--
C
C  Computes the product P * A of an M x M plane rotation P times 
C  an M x N matrix A, where
C
C                           [ 1 ..  0 .. 0 .. 0 ]
C                           [ . ..  . .. . .. . ]
C                           [ 0 ..  c .. s .. 0 ]  <-- j1
C     P = P(j1, j2, c, s) = [ . ..  . .. . .. . ]
C                           [ 0 .. -s .. c .. 0 ]  <-- j2
C                           [ . ..  . .. . .. . ]
C                           [ 0 ..  0 .. 0 .. 1 ]
C
C  Variables in the calling sequence
C  ---------------------------------
C  M    I   IN   Number of rows of the matrix A
C  N    I   IN   Number of columns of the matrix A
C  A    D   IN   Array of dimension LDA x N containing the matrix A
C           OUT  The product P * A
C  LDA  I   IN   The leading dimension of A, LDA >= M
C  J1   I   IN   The index of the first variable defining the plane
C                of rotation, 1 <= J1 <= M
C  J2   I   IN   The index of the second variable defining the plane
C                of rotation, 1 <= J2 <= M
C  C    D   IN   The cosine value of the rotation
C  S    D   IN   The sine value of the rotation
C  IER  I   OUT  Error indicator
C                IER = 0  no error
C                IER = 1  input data error
C
C234567--x---------x---------x---------x---------x---------x---------x--
C
C.....Parameters
C
      DOUBLE PRECISION ONE, ZER
      PARAMETER ( ONE = 1.0D+0, ZER = 0.0D+0 )
C
C.....Local variables
C
      INTEGER I
      DOUBLE PRECISION TMP
C
C.......................Executable statements...........................
C
C.....Test the input data
C
      IF( (M.LT.0) .OR. (N.LT.0) .OR. (LDA.LT.M)
     &             .OR. (J1.LT.1) .OR. (J1.GT.M)
     &             .OR. (J2.LT.1) .OR. (J2.GT.M)
     &             .OR. (J1.EQ.J2) )THEN
         CALL MSGPRT('ROTML',' Zero pivot encountered') 
         IER = 1
         RETURN
      ENDIF
C
C.....Quick return if possible
C
      IER = 0
      IF( ( M.EQ.0 ) .OR. ( N.EQ.0 ) )RETURN
C
C.....Compute the product
C
      IF( ( C.NE.ONE ) .OR. ( S.NE.ZER ) ) THEN
         DO 10 I = 1, N
            TMP = A(J2,I)
            A(J2,I) = C*TMP - S*A(J1,I)
            A(J1,I) = S*TMP + C*A(J1,I)
   10    CONTINUE
      ENDIF
C
      RETURN
C
C.....End of ROTML
C
      END
C
C234567==1=========2=========3=========4=========5=========6=========7==
C
      SUBROUTINE ROTMR( M, N, A, LDA, J1, J2, C, S, IER )
C
      INTEGER M,N,LDA,J1,J2,IER
      DOUBLE PRECISION A(LDA,*),C,S
C
C234567--x---------x---------x---------x---------x---------x---------x--
C
C  Computes the product A * P' of an M x N matrix A times the
C  transpose P' of an M x M plane rotation P, where
C
C                           [ 1 ..  0 .. 0 .. 0 ]
C                           [ . ..  . .. . .. . ]
C                           [ 0 ..  c .. s .. 0 ]  <-- j1
C     P = P(j1, j2, c, s) = [ . ..  . .. . .. . ]
C                           [ 0 .. -s .. c .. 0 ]  <-- j2
C                           [ . ..  . .. . .. . ]
C                           [ 0 ..  0 .. 0 .. 1 ]
C
C  Variables in the calling sequence
C  ---------------------------------
C  M    I   IN   Number of rows of the matrix A
C  N    I   IN   Number of columns of the matrix A
C  A    D   IN   Array of dimension LDA x N containing the matrix A
C           OUT  The product P * A
C  LDA  I   IN   The leading dimension of A, LDA >= M
C  J1   I   IN   The index of the first variable defining the plane
C                of rotation, 1 <= J1 <= M
C  J2   I   IN   The index of the second variable defining the plane
C                of rotation, 1 <= J2 <= M
C  C    D   IN   The cosine value of the rotation
C  S    D   IN   The sine value of the rotation
C  IER  I   OUT  Error indicator
C                IER = 0  no error
C                IER = 1  input data error
C
C234567--x---------x---------x---------x---------x---------x---------x--
C
C.....Parameters
C
      DOUBLE PRECISION ONE, ZER
      PARAMETER ( ONE = 1.0D+0, ZER = 0.0D+0 )
C
C.....Local variables
C
      INTEGER I
      DOUBLE PRECISION TMP
C
C.......................Executable statements...........................
C
C.....Test the input data
C
      IF( (M.LT.0) .OR. (N.LT.0) .OR. (LDA.LT.M)
     &             .OR. (J1.LT.1) .OR. (J1.GT.M)
     &             .OR. (J2.LT.1) .OR. (J2.GT.M)
     &             .OR. (J1.EQ.J2) )THEN
         CALL MSGPRT('ROTMR',' Zero pivot encountered') 
         IER = 1
         RETURN
      ENDIF
C
C.....Quick return if possible
C
      IER = 0
      IF( ( M.EQ.0 ) .OR. ( N.EQ.0 ) )RETURN
C
C.....Compute the product
C
      IF( ( C.NE.ONE ) .OR. ( S.NE.ZER ) ) THEN
         DO 10 I = 1, N
            TMP = A(I,J2)
            A(I,J2) = C*TMP - S*A(I,J1)
            A(I,J1) = S*TMP + C*A(I,J1)
   10    CONTINUE
      ENDIF
C
      RETURN
C
C.....End of ROTML
C
      END
C
C234567==1=========2=========3=========4=========5=========6=========7==
C
      SUBROUTINE SVD(M,N,D,E,U,LDU,VT,LDV,EPMACH,SAFMIN,JOB,IER)
C
      INTEGER M,N,LDU,LDV,JOB,IER
      DOUBLE PRECISION D(*),E(*),U(LDU,*),VT(LDV,*)
      DOUBLE PRECISION EPMACH,SAFMIN
C
C234567--x---------x---------x---------x---------x---------x---------x--
C
C  Compute the singular value decomposition (SVD) of a real
C  M x N, M >= N, upper bidiagonal matrix
C
C             B = Q * S * P'
C
C  where P' is the transpose of P, S is an M x N diagonal matrix 
C  with non-negative diagonal elements (the singular values of B),
C  and the diagonal D and upper diagonal E of the matrix B are given.
C
C  The routine computes S, and optionally the matrices U * Q, and
C  P' * VT, for given real input matrices U, and VT.
C
C  The routine is based on the Lapack routine DBDSQR and the
C  algorithm discussed in "Computing  Small Singular Values of 
C  Bidiagonal Matrices With Guaranteed High Relative Accuracy," by 
C  J. Demmel and W. Kahan, LAPACK Working Note #3. 
C
C  Variables in the calling sequence
C  ---------------------------------
C  M      I   IN    Number of rows of B
C  N      I   IN    Number of columns of B, M >= N
C  D      D   IN    Array of dimension N containing the N diagonal 
C                   elements of the bidiagonal matrix B
C             OUT   The singular values of B in decreasing order
C  E      D   IN    Array of dimension N-1 off-diagonal elements of 
C                   the bidiagonal matrix B.
C             OUT   The matrix is destroyed
C  VT     D   IN    Array of dimension LDV by N, LDV >= N, containing
C                   a given N x N matrix VT
C             OUT   The product P' * VT
C  LDV    I   IN    The leading dimension of the array VT,
C                   LDV >= max(1,N)
C  U      D   IN    Array of dimension LDU x M, LDU >= M, containing
C                   a given M x M matrix U
C             OUT   The product U * Q
C  LDU    I   IN    The leading dimension of the array U,
C                   LDU >= max(1,M)
C  EPMACH D   IN    Machine epsilon, the smallest number such
C                   that 1.0D0 + EPMACH .NE. 1.0D0
C  SAFMIN D   IN    Safe minimum such that 1.0D0/SAFMIN does 
C                   not overflow    
C  JOB    I   IN    Job indicator
C                   JOB = 0  No singular vectors are computed
C                   JOB = 1  The singular vectors are computed
C  IER    I   OUT   Error indicator
C                   IER = 0:  successful exit
C                   IER = 1:  illegal data entry
C                   IER < 0:  the algorithm did not converge; 
C                             D and E contain the elements of a 
C                             bidiagonal matrix which is orthogonally
C                             similar to the input matrix B;  if 
C                             | IER | = i, then i elements of E have 
C                             not converged to zero.
C
C234567--x---------x---------x---------x---------x---------x---------x--
C
C.....Parameters
C
      INTEGER MAXITR
      PARAMETER( MAXITR = 6 )
      DOUBLE PRECISION  ZER,ONE,TEN,HNDRD,HNDRTH,MEIGTH
      PARAMETER ( ZER = 0.0D0, ONE = 1.0D0, TEN = 1.0D1,
     &         HNDRD = 1.0D2, HNDRTH = 1.0D-2, MEIGTH = -0.125D0 )
C
C.....Local variables
C
      LOGICAL SINGVC
      INTEGER I,IDIR,ISUB,ITER,J,JJ,LL,LLL, 
     &        MPT,MAXIT,NM1,NP1MI,OLDLL,OLDM
      DOUBLE PRECISION ABSE,ABSS,COSL,COSR,CS,DN,F,G,GAP,GMAX,
     &                 H,MU,OLDCS,OLDSN,R,SAFBD,SHIFT,SIGMN,
     &                 SIGMX,SINL,SINR,SLL,SMAX,SMIN,SMINL,
     &                 SMINLO,SMINOA,SN,TMP,THRESH,TOL,TOLMUL
C
      DATA OLDSN/0.0D0/
C
C.......................Executable statements...........................
C
C.....Test input data
C
      IER = 0
      IF( (M.LT.0) .OR. (N.LT.0) .OR. (M.LT.N) 
     &             .OR. (LDV.LT.N) .OR. (LDU.LT.M) ) THEN
         CALL MSGPRT('SVD',' Zero pivot encountered') 
         IER = 1
         RETURN
      ENDIF
C
C.....Quick return if possible
C
      IER = 0
      IF( ( M.EQ.0 ) .OR. ( N.EQ.0 ) )RETURN
      IF( N .EQ. 1 ) GOTO 200
C
C.....SINGVC is true if any singular vectors desired, false otherwise
C
      SINGVC = JOB .NE. 0
      DN = DBLE(N)
C
C.....Set tolerance
C
      TOLMUL = EPMACH**MEIGTH
      IF( TOLMUL .GT. HNDRD ) TOLMUL = HNDRD
      IF( TOLMUL .LT. TEN ) TOLMUL = TEN
      TOL = TOLMUL*EPMACH
C
C.....Compute approximate maximum, minimum singular values
C
      SMAX = ABS(D(N))
      NM1 = N - 1
      IF( N .GT. 1)THEN
         DO 20 I = 1, NM1
            IF( ABS(D(I)) .GT. SMAX) SMAX = ABS(D(I))
            IF( ABS(E(I)) .GT. SMAX) SMAX = ABS(E(I))
   20    CONTINUE
      ENDIF
      SMINL = ZER
      SMINOA = ABS(D(1))
      IF( SMINOA .NE. ZER )THEN
         MU = SMINOA
         IF( N .GT. 1)THEN
            DO 30 I = 2, N
               MU = ABS(D(I))*( MU / (MU + ABS(E(I-1) )))
               IF( MU .LT. SMINOA) SMINOA = MU
               IF( SMINOA .EQ. ZER )GO TO 40
   30       CONTINUE
         ENDIF
      ENDIF
   40 CONTINUE
      SMINOA = SMINOA / SQRT(DN)
C
C.....Prepare for main iteration loop for the singular values
C
      MAXIT = MAXITR*N*N
      SAFBD = DBLE(MAXIT)*SAFMIN
      ITER = 0
      OLDLL = -1
      OLDM = -1
      THRESH = TOL*SMINOA
      IF( THRESH .LT. SAFBD) THRESH = SAFBD
C
C.....MPT points to last entry of unconverged part of matrix
C
      MPT = N
C
C.....Begin main iteration loop
C
   50 CONTINUE
C
C.....Check for convergence
C
      IF( MPT .LE. 1 )GOTO 200
C
C.....Check for maximal iteration count
C
      IF( ITER .GT. MAXIT )THEN
C
C........Maximum number of iterations exceeded, 
C........Convergence failure
C
         IER = 0
         DO 60 I = 1, NM1
            IF( E(I) .NE. ZER )IER = IER-1
   60    CONTINUE
         CALL MSGPRT('SVD',' The algorithm failed to converge') 
         RETURN
      ENDIF
C
C.....Find diagonal block of matrix to work on
C
      SMAX = ABS(D(MPT))
      SMIN = SMAX
      DO 70 LLL = 1, MPT
         LL = MPT - LLL
         IF( LL .EQ. 0 ) GOTO 90
         ABSS = ABS(D(LL))
         IF( ABSS .LT. SMIN )SMIN = ABSS
         IF( ABSS .GT. SMAX )SMAX = ABSS
         ABSE = ABS(E(LL))
         IF( ABSE .LE. THRESH ) GOTO 80
         IF( ABSE .GT. SMAX )SMAX = ABSE
   70 CONTINUE
C
   80 CONTINUE
      E(LL) = ZER
C
C.....Matrix splits since E(LL) = 0
C
      IF( LL .EQ. MPT-1 ) THEN
C
C........Convergence of bottom singular value, return to top of loop
C
         MPT = MPT - 1
         GOTO 50
      ENDIF
C
   90 CONTINUE
      LL = LL + 1
C
C.....E(LL) through E(MPT-1) are nonzero, E(LL-1) is zero
C
      IF( LL .EQ. MPT-1 ) THEN
C
C........2 by 2 block, handle separately
C
         CALL SVD2D( D(MPT-1),E(MPT-1),D(MPT),SIGMN,SIGMX,
     &               SINR,COSR,SINL,COSL,EPMACH,1 )
         D(MPT-1) = SIGMX
         E(MPT-1) = ZER
         D(MPT)   = SIGMN
C
C........Compute singular vectors, if desired
C
         IF( SINGVC ) THEN
            CALL ROTML(N,N,VT,LDV,MPT-1,MPT,COSR,SINR,IER)
            CALL ROTMR(M,M,U ,LDU,MPT-1,MPT,COSL,SINL,IER)
            IF( IER .EQ. 1) RETURN
         ENDIF
         MPT = MPT - 2
         GOTO 50
      ENDIF
C
C.....If working on new submatrix, choose shift direction
C.....(from larger end diagonal entry towards smaller)
C
      IF( LL.GT.OLDM .OR. MPT.LT.OLDLL ) THEN
         IF( ABS(D(LL)) .GE. ABS(D(MPT)) ) THEN
C
C...........Chase bulge from top (big end) to bottom (small end)
C
            IDIR = 1
         ELSE
C
C...........Chase bulge from bottom (big end) to top (small end)
C
            IDIR = 2
         ENDIF
      ENDIF
C
C.....Apply convergence tests
C
      IF( IDIR .EQ. 1 ) THEN
C
C........Run convergence test in forward direction
C........First apply standard test to bottom of matrix
C
         IF( (ABS(E(MPT-1)) .LE. ABS(TOL)*ABS(D(MPT))) ) THEN
            E(MPT-1) = ZER
            GOTO 50
         ENDIF
C
C........Apply convergence criterion forward
C
         MU = ABS(D(LL))
         SMINL = MU
         DO 100 LLL = LL, MPT - 1
            IF( ABS(E(LLL) ) .LE. TOL*MU ) THEN
               E(LLL) = ZER
               GOTO 50
            ENDIF
            SMINLO = SMINL
            MU = ABS(D(LLL+1))*( MU / ( MU+ABS(E(LLL) )))
            IF( SMINL .GT. MU) SMINL = MU
  100    CONTINUE
C
C........If singular values only wanted, apply gap test to bottom
C........end of matrix
C
         IF( .NOT. SINGVC ) THEN
            GAP = SMINLO / SQRT(DBLE(MPT-LL)) - ABS(D(MPT))
            IF( GAP .GT. ZER ) THEN
               GMAX = GAP
               ABSS = ABS(D(MPT))
               IF( GMAX .LT. ABSS) GMAX = ABSS
               ABSE = ABS(E(MPT-1))
               IF( GMAX .LT. ABSE) GMAX = ABSE
               IF( ( ABSE/GMAX )**2 
     &               .LE. TOL*( GAP/GMAX )*( ABSS/GMAX ) ) THEN
                  E(MPT-1) = ZER
                  GOTO 50
               ENDIF
            ENDIF
         ENDIF
      ELSE
C
C........Run convergence test in backward direction
C........First apply standard test to top of matrix
C
         IF( ( ABS(E(LL)) .LE. ABS(TOL)*ABS(D(LL)) )) THEN
            E(LL) = ZER
            GOTO 50
         ENDIF
         MU = ABS(D(MPT))
         SMINL = MU
         DO 110 LLL = MPT-1, LL, -1
            IF( ABS(E(LLL)) .LE. TOL*MU ) THEN
               E(LLL) = ZER
               GOTO 50
            ENDIF
            SMINLO = SMINL
            MU = ABS(D(LLL ))*( MU / ( MU+ABS(E(LLL)) ) )
            IF( SMINL .GT. MU) SMINL = MU
  110    CONTINUE
C
C........If singular values only wanted, apply gap test to top
C........end of matrix
C
         IF( .NOT. SINGVC ) THEN
            GAP = SMINLO / SQRT( DBLE(MPT-LL) ) - ABS(D(LL))
            IF( GAP .GT. ZER ) THEN
               ABSS = ABS(D(LL))
               IF( GAP .LT. ABSS) GAP = ABSS
               ABSE = ABS(E(LL))
               IF( GAP .LT. ABSE) GAP = ABSE
               IF( (ABSE / GMAX )**2 .LE. TOL*(GAP/GMAX)*
     &                (ABSS/GMAX) ) THEN
                   E(LL) = ZER
                   GOTO 50
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      OLDLL = LL
      OLDM = MPT
C
C.....Compute shift.  First, test if shifting would ruin relative
C.....accuracy, and if so set the shift to zero.
C
      IF( DN*TOL*(SMINL/SMAX) .LE. MAX(EPMACH,HNDRTH*TOL) )THEN
C
C.........Use a zero shift to avoid loss of relative accuracy
C
         SHIFT = ZER
      ELSE
C
C........Compute the shift from 2-by-2 block at end of matrix
C
         IF( IDIR .EQ. 1 ) THEN
            SLL = ABS(D(LL))
            JJ = MPT-1
         ELSE
            SLL = ABS(D(MPT))
            JJ = LL
         ENDIF
         CALL SVD2D( D(JJ),E(JJ),D(JJ+1),SHIFT,R,
     &               TMP,TMP,TMP,TMP,EPMACH,0 )
C
C........Test if shift negligible, and if so set to zero
C
         IF( SLL .GT. ZER ) THEN
            IF( ( SHIFT/SLL )**2 .LT. EPMACH ) SHIFT = ZER
         ENDIF
      ENDIF
C
C.....Increment iteration count
C
      ITER = ITER + MPT - LL
C
C.....If SHIFT = 0, do simplified QR iteration
C
      IF( SHIFT .EQ. ZER ) THEN
         IF( IDIR .EQ. 1 ) THEN
C
C...........Chase bulge from top to bottom and update
C...........singular vectors if desired
C
            CS = ONE
            OLDCS = ONE
            DO 120 I = LL, MPT-1
               CALL ROTG( D(I)*CS, E(I), CS, SN, R )
               IF( I .GT. LL) E(I-1) = OLDSN*R
               CALL ROTG( OLDCS*R, D(I+1)*SN, OLDCS, OLDSN, D(I))
               IF( SINGVC ) THEN
                  CALL ROTML(N,N,VT,LDV,I,I+1,CS,SN,IER)
                  CALL ROTMR(M,M,U ,LDU,I,I+1,OLDCS,OLDSN,IER)
                  IF( IER .EQ. 1) RETURN
               ENDIF
  120       CONTINUE
            H = D(MPT)*CS
            D(MPT) = H*OLDCS
            E(MPT-1) = H*OLDSN
C
C...........Test convergence
C
            IF( ABS( E(MPT-1) ) .LE. THRESH ) E(MPT-1) = ZER
C
         ELSE
C
C...........Chase bulge from bottom to top and update
C...........singular vectors if desired
C
            CS = ONE
            OLDCS = ONE
            DO 130 I = MPT, LL+1, -1
               CALL ROTG( D(I)*CS, E(I-1), CS, SN, R )
               IF( I .LT. MPT )E(I) = OLDSN*R
               CALL ROTG( OLDCS*R, D(I-1)*SN, OLDCS, OLDSN, D(I))
               IF( SINGVC ) THEN
                  CALL ROTML(N,N,VT,LDV,I,I+1,OLDCS,-OLDSN,IER)
                  CALL ROTMR(M,M,U ,LDU,I,I+1,CS,-SN,IER)
                  IF( IER .EQ. 1) RETURN
               ENDIF
  130       CONTINUE
            H = D(LL)*CS
            D(LL) = H*OLDCS
            E(LL) = H*OLDSN
C
C...........Test convergence
C
            IF( ABS( E(LL) ) .LE. THRESH )E(LL) = ZER
         ENDIF
C
      ELSE
C
C........Use nonzero shift
C
         IF( IDIR .EQ. 1 ) THEN
C
C...........Chase bulge from top to bottom and update
C...........singular vectors if desired
C
            F = (ABS(D(LL))-SHIFT)*( SIGN(ONE,D(LL)) + SHIFT/D(LL) )
            G = E(LL)
C
            DO 140 I = LL, MPT-1
               CALL ROTG( F, G, COSR, SINR, R )
               IF( I .GT. LL ) E(I-1) = R
               F = COSR*D(I) + SINR*E(I)
               E(I) = COSR*E(I) - SINR*D(I)
               G = SINR*D(I+1)
               D(I+1) = COSR*D(I+1)
               CALL ROTG( F, G, COSL, SINL, R )
               D(I) = R
               F = COSL*E(I) + SINL*D(I+1)
               D(I+1) = COSL*D(I+1) - SINL*E(I)
               IF( I .LT. MPT-1 )THEN               
                  G = SINL*E(I+1)
                  E(I+1) = COSL*E(I+1)
               ELSE
                  E(MPT-1) = F
               ENDIF
               IF( SINGVC ) THEN
                  CALL ROTML(N,N,VT,LDV,I,I+1,COSR,SINR,IER)
                  CALL ROTMR(M,M,U ,LDU,I,I+1,COSL,SINL,IER)
                  IF( IER .EQ. 1) RETURN
               ENDIF
  140       CONTINUE
C
C...........Test convergence
C
            IF( ABS( E(MPT-1) ).LE.THRESH ) E(MPT-1) = ZER
C
         ELSE
C
C...........Chase bulge from bottom to top and
C...........update singular vectors if desired
C
            F = (ABS(D(MPT)) - SHIFT)*(SIGN(ONE,D(MPT)) + SHIFT/D(MPT))
            G = E(MPT-1)
C
            DO 150 I = MPT, LL+1, -1
               CALL ROTG( F, G, COSR, SINR, R )
               IF( I .LT. MPT ) E(I) = R
               F = COSR*D(I) + SINR*E(I-1)
               E(I-1) = COSR*E(I-1) - SINR*D(I)
               G = SINR*D(I-1)
               D(I-1) = COSR*D(I-1)
               CALL ROTG( F, G, COSL, SINL, R )
               D(I) = R
               F = COSL*E(I-1) + SINL*D(I-1)
               D(I-1) = COSL*D(I-1) - SINL*E(I-1)
               IF( I .GT. LL+1 ) THEN
                  G = SINL*E(I-2)
                  E(I-2) = COSL*E(I-2)
               ELSE
                  E(LL) = F
               ENDIF
               IF( SINGVC ) THEN
                  CALL ROTML(N,N,VT,LDV,I,I+1,COSL,-SINL,IER)
                  CALL ROTMR(M,M,U ,LDU,I,I+1,COSR,-SINR,IER)
                  IF( IER .EQ. 1) RETURN
               ENDIF
  150       CONTINUE
C
C...........Test convergence
C
            IF( ABS( E(LL) ) .LE. THRESH ) E(LL) = ZER
C
         ENDIF
      ENDIF
C
C.....QR iteration finished, go back and check convergence
C
      GOTO 50
C
C.....All singular values converged, so make them positive
C
  200 CONTINUE
C
      DO 220 I = 1, N
         IF( D(I) .LT. ZER ) THEN
            D(I) = -D(I)
C
C...........Change sign of singular vectors, if desired
C
            IF( SINGVC )THEN
               DO 210 JJ = 1,N
                  VT(I,JJ) = -VT(I,JJ)
  210          CONTINUE
            ENDIF
         ENDIF
  220 CONTINUE
C
C.....Sort the singular values into decreasing order (insertion sort on
C.....singular values, but only one transposition per singular vector)
C
      DO 260 I = 1, NM1
C
C........Scan for smallest D(I)
C
         ISUB = 1
         SMIN = D(1)
         NP1MI = N + 1 - I
         DO 230 J = 2, NP1MI
            IF( D(J) .LE. SMIN ) THEN
               ISUB = J
               SMIN = D(J)
            ENDIF
  230    CONTINUE
         IF( ISUB .NE. NP1MI ) THEN
C
C...........Swap singular values and vectors
C
            D(ISUB)  = D(NP1MI)
            D(NP1MI) = SMIN
            IF( SINGVC ) THEN
               DO 240 JJ = 1, N
                  TMP          = VT(ISUB,JJ)
                  VT(ISUB,JJ)  = VT(NP1MI,JJ)
                  VT(NP1MI,JJ) = TMP
  240          CONTINUE
               DO 250 JJ = 1, M
                  TMP         = U(JJ,ISUB)
                  U(JJ,ISUB)  = U(JJ,NP1MI)
                  U(JJ,NP1MI) = TMP
  250          CONTINUE
            ENDIF
         ENDIF
  260 CONTINUE
C
      RETURN
C
C.....End of SVD
C
      END
C
C234567==1=========2=========3=========4=========5=========6=========7==
C
      SUBROUTINE SVD2D( F,G,H,SSMIN,SSMAX,SNR,CSR,SNL,CSL,EPMACH,JOB )
C
      INTEGER JOB
      DOUBLE PRECISION CSL,CSR,F,G,H,SNL,SNR,SSMAX,SSMIN,EPMACH
C
C234567--x---------x---------x---------x---------x---------x---------x--
C
C  Computes the singular value decomposition of a 2-by-2 triangular 
C  matrix
C                [  F   G  ]
C                [  0   H  ].
C  On return, abs(SSMAX) is the larger singular value, abs(SSMIN) 
C  is the smaller singular value, and (CSL,SNL) and (CSR,SNR) are 
C  the left and right singular vectors for abs(SSMAX) in the 
C  decomposition
C
C     [ CSL  SNL ] [  F   G  ] [ CSR -SNR ]  =  [ SSMAX   0   ]
C     [-SNL  CSL ] [  0   H  ] [ SNR  CSR ]     [  0    SSMIN ].
C
C  Any input parameter may be aliased with any output parameter.
C
C  Barring over/underflow and assuming a guard digit in subtraction,
C  all output quantities are correct to within a few units in the 
C  last place (ulps).
C
C  In IEEE arithmetic, the code works correctly if one matrix entry
C  is infinite.
C
C  Overflow will not occur unless the largest singular value itself
C  overflows or is within a few ulps of overflow. (On machines with
C  partial overflow, like the Cray, overflow may occur if the largest
C  singular value is within a factor of 2 of overflow.)
C
C  Underflow is harmless if underflow is gradual. Otherwise, results
C  may correspond to a matrix modified by perturbations of size near
C  the underflow threshold.
C
C  Edited Version of LAPACK routine DLASV2.
C
C  Variables in the calling sequence
C  ---------------------------------
C  F      D   IN   The (1,1) entry of the 2-by-2 matrix.
C  G      D   IN   The (1,2) entry of the 2-by-2 matrix.
C  H      D   IN   The (1,1) entry of the 2-by-2 matrix.
C  SSMIN  D   OUT  abs(SSMIN) is the smaller singular value.
C  SSMAX  D   OUT  abs(SSMAX) is the larger singular value.
C  SNL    D   OUT  first component of a unit left singular vector
C                  for the singular value abs(SSMAX)
C  CSL    D   OUT  second component of a unit left singular vector
C                  for the singular value abs(SSMAX)
C  SNR    D   OUT  first component of a unit right singular vector
C                  for the singular value abs(SSMAX)
C  CSR    D   OUT  first component of a unit right singular vector
C                  for the singular value abs(SSMAX)
C  EPMACH D   IN   Machine epsilon, the smallest number such that
C                  1.0D0 + EPMACH .NE. 1.0D0
C  JOB    I   IN   Job indicator
C                  JOB = 0 compute singular values SSMIN,SSMAX only
C                  JOB = 1 compute singular values and vectors
C
C234567--x---------x---------x---------x---------x---------x---------x--
C
C.....Intrinsic functions called
C
      DOUBLE PRECISION ABS,SIGN,SQRT
C
C.....Parameters
C
      DOUBLE PRECISION ZER,HALF,ONE,TWO,FOUR
      PARAMETER  ( ZER = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, 
     &             TWO = 2.0D0, FOUR = 4.0D0 )
C
C.....Local variables
C
      LOGICAL GASMAL,SINGVC,SWAP
      INTEGER PMAX
      DOUBLE PRECISION A,CLT,CRT,D,FA,FT,GA,GT,HA,HT,L,M,
     &                 MM,R,S,SLT,SRT,T,TMP,TSIGN,TT
C
C.......................Executable statements...........................
C
      SINGVC = (JOB .EQ. 0)
      FT = F
      GT = G
      HT = H
      FA = ABS(FT)
      GA = ABS(GT)
      HA = ABS(HT)
C
C.....PMAX points to the maximum absolute entry of matrix
C..........PMAX = 1 if F largest in absolute values
C..........PMAX = 2 if G largest in absolute values
C..........PMAX = 3 if H largest in absolute values
C
      PMAX = 1
      SWAP = ( HA .GT. FA )
      IF( SWAP ) THEN
C
C........Ensure FA .ge. HA
C
         PMAX = 3
         TMP  = FT
         FT   = HT
         HT   = TMP
         TMP  = FA
         FA   = HA
         HA   = TMP
      ENDIF
C
      IF( GA .EQ. ZER ) THEN
C
C........Diagonal matrix
C
         SSMIN = HA
         SSMAX = FA
C
C........Return if no singular vectors are wanted
C
         IF( SINGVC ) RETURN
C
         CLT = ONE
         CRT = ONE
         SLT = ZER
         SRT = ZER
      ELSE
         GASMAL = .TRUE.
         IF( GA .GT. FA ) THEN
            PMAX = 2
            IF( ( FA/GA ) .LT. EPMACH ) THEN
C
C..............Case of very large GA
C
               GASMAL = .FALSE.
               SSMAX = GA
               IF( HA .GT. ONE ) THEN
                  SSMIN = FA / (GA/HA)
               ELSE
                  SSMIN = ( FA/GA )*HA
               ENDIF
C
C..............Return if no singular vectors are wanted
C
               IF( SINGVC ) RETURN
C
               CLT = ONE
               SLT = HT/GT
               SRT = ONE
               CRT = FT/GT
            ENDIF
         ENDIF
C
         IF( GASMAL ) THEN
C
C...........Normal case
C
            D = FA - HA
            IF( D .EQ. FA ) THEN
C
C..............Copes with infinite F or H
C
               L = ONE
            ELSE
               L = D/FA
            ENDIF
C
C...........Note that 0 .le. L .le. 1
C
            M = GT/FT
C
C...........Note that ABS(M) .LE. 1.0D0/EPMACH
C
            T = TWO - L
C
C...........Note that T .ge. 1
C
            MM = M*M
            TT = T*T
            S = SQRT(TT + MM)
C
C...........Note that 1.0D0 .LE. S .LE. 1.0D0 + 1.0D0/EPMACH
C
            IF( L .EQ. ZER ) THEN
               R = ABS(M)
            ELSE
               R = SQRT(L*L+MM)
            ENDIF
C
C...........Note that 0 .LE. R .LE. 1.0D0 + 1.0D0/EPMACH
C
            A = HALF*(S + R)
C
C...........Note that 1 .le. A .le. 1 + abs(M)
C
            SSMIN = HA/A
            SSMAX = FA*A
C
C...........Return if no singular vectors are wanted
C
            IF( SINGVC ) RETURN
C
            IF( MM .EQ. ZER ) THEN
C
C..............Note that M is very tiny
C
               IF( L .EQ. ZER ) THEN
                  T = SIGN(TWO,FT )*SIGN(ONE,GT )
               ELSE
                  T = GT/SIGN(D,FT ) + M/T
               ENDIF
            ELSE
               T = ( M/(S+T) + M/(R+L) )*(ONE + A)
            ENDIF
            L = SQRT(T*T+FOUR)
            CRT = TWO/L
            SRT = T/L
            CLT = ( CRT+SRT*M )/A
            SLT = ( HT/ FT)*SRT/A
         ENDIF
      ENDIF
      IF( SWAP ) THEN
         CSL = SRT
         SNL = CRT
         CSR = SLT
         SNR = CLT
      ELSE
         CSL = CLT
         SNL = SLT
         CSR = CRT
         SNR = SRT
      ENDIF
C
C.....Correct signs of SSMAX and SSMIN
C
      IF( PMAX .EQ. 1 )
     &          TSIGN = SIGN(ONE,CSR )*SIGN(ONE,CSL)*SIGN(ONE,F)
      IF( PMAX .EQ. 2 )
     &          TSIGN = SIGN(ONE,SNR )*SIGN(ONE,CSL )*SIGN(ONE,G)
      IF( PMAX.EQ.3 )
     &          TSIGN = SIGN(ONE,SNR)*SIGN(ONE,SNL)*SIGN(ONE,H)
      SSMAX = SIGN(SSMAX,TSIGN)
      SSMIN = SIGN(SSMIN,TSIGN*SIGN(ONE,F)*SIGN(ONE,H))
C
      RETURN
C
C.....End of SVD2D
C
      END
C
C234567==1=========2=========3=========4=========5=========6=========7==
C
      SUBROUTINE MSGPRT(PROGR,MESSG)
C
      CHARACTER*(*) PROGR,MESSG
C
C234567--x---------x---------x---------x---------x---------x---------x--
C
C     This routine processes and prints a diagnostic message which may 
C     or may not represent an error. It follows closely the design of 
C     the SLATEC error processing routine XERMSG. However, MSGPRT only 
C     prints the supplied message and then returns control to the 
C     calling routine.
C
C     At the first call of MSGPRT the routine calls the integer
C     function LUNIT to obtain an output unit number which will 
C     be used for all subsequent output from the routine. 
C     None of the arguments to MSGPRT are modified by the routine
C
C  Variables in the calling sequence
C  ---------------------------------
C     PROGR  C  IN  Character constant or variable with the name 
C                   of the routine that sends the message.
C     MESSG  C  IN  Character constant or variable with the text of 
C                   the message. In the example below, the message 
C                   is a character constant that contains a generic 
C                   message.
C
C                    CALL MESSG ('PROGR',
C                   &   'The matrix-order exceeds the row dimension')
C
C                   where 'PROGR' is the name of the calling program.
C                   It is possible to generate a message that contains 
C                   numeric values. Specific numeric values can be 
C                   converted into character strings using formatted 
C                   WRITE statements into character variables. This is 
C                   called standard Fortran internal file I/O and is 
C                   exemplified in the first three lines of the 
C                   following example which also shows the use of // 
C                   for the concatenation of substrings of characters
C
C                      CHARACTER*5 CHARN, CHARL
C                      WRITE (CHARN,10) N
C                      WRITE (CHARL,10) LDA
C                   10 FORMAT(I5)
C                      CALL MESSG ('PROGR', 'The order'//CHARN//
C                     &   ' of the matrix exceeds the row dimension'//
C                     &   CHARL)
C
C                   By judicious use of // it can be avoided that some
C                   character constant is continued to the next line.
C                   Morover, concatenation of several parts of a 
C                   message has the advantage, over the formation of 
C                   a long message as one large character variable, 
C                   that it avoids the need for knowing how long the 
C                   message will be in order to declare an adequate 
C                   length for that large a character variable.  
C                   Long messages will be broken into pieces of 72 
C                   characters (set by NWRAP) for printing on multiple 
C                   lines. The error message is scanned backwards to 
C                   ignore trailing blanks; moreover, the substring 
C                   '$$' is treated as a new line sentinel. Thus the 
C                   use of the separator '$$' can avoid the need for
C                   counting multiples of 72 characters; but, of
C                   course, '$$' must occur within 72 characters of
C                   the start of each line to have its intended effect.
C
C234567--x---------x---------x---------x---------x---------x---------x--
C
C.....Subroutine called
C
      EXTERNAL LUNIT
C
C.....Local variables
C
      INTEGER I,INDEX,IDELTA,LEN,LPIECE,LWRAP,LENMSG
      INTEGER MAX,MIN,N,NEXTC
      CHARACTER*148 CBUFF
C
C.....Variables to be saved between calls
C
      LOGICAL NOPRNT,FIRST
      INTEGER KL,LOUT(5)
      SAVE NOPRNT,FIRST,KL,LOUT
C
C.....Data
C
      INTEGER NWRAP
      CHARACTER*2 NEWLIN
      DATA FIRST/.TRUE./, NWRAP/72/, NEWLIN/'$$'/
C
C.......................Executable statements...........................
C
C.....Determine the output units
C
      IF (FIRST) THEN
         FIRST = .FALSE.
         CALL LUNIT(KL, LOUT)
         NOPRNT = ( KL .LT. 1 .OR. KL .GT. 5 )
      ENDIF
C
      IF(NOPRNT)RETURN
C
C.....Get the name of the subroutine reporting the message
C
      CBUFF(1:14) = ' Message from '
      I = MIN(LEN(PROGR),10)
      N = 14 + I      
      CBUFF(15:N) = PROGR(1:I)
C
C.....Print the name of the subroutine
C
      DO 10 I = 1, KL
         WRITE(LOUT(I),200) CBUFF(1:N)
   10 CONTINUE
C     
C.....LWRAP is the maximum number of characters to be taken at one
C.....time from MESSG to print on one line
C
      LWRAP = MAX(16, MIN(132, NWRAP))
C
C.....Set LENMSG to the length of MESSG, ignore any trailing blanks
C
      LENMSG = LEN(MESSG)
      N = LENMSG
      DO 20 I=1,N
         IF (MESSG(LENMSG:LENMSG) .NE. ' ') GOTO 30
         LENMSG = LENMSG - 1
   20 CONTINUE
C
C.....If the message is all blanks, then print a blank line
C
   30 CONTINUE
      IF (LENMSG .EQ. 0) THEN
         CBUFF(1:1) = ' '
         DO 40 I = 1, KL
            WRITE(LOUT(I),200) CBUFF(1:1)
   40    CONTINUE
         RETURN
      ENDIF
C
C.....Set NEXTC to the position in MESSG where the next substring
C.....starts. From this position the line is scanned for the new line
C.....sentinel. The intrinsic function INDEX returns zero if there is
C ....no occurrence or if the length of the second argument.
C
C.....Control loops back to label 50 until all pieces are printed.
C.....When NEXTC exceeds LENMSG, there is no more to print.
C
C.....The variable LPIECE gives the number of characters to be
C.....taken from MESSG starting at the position NEXTC. Several 
C.....cases of LPIECE have to be checked in the following order.
C
C.....LPIECE .EQ. 0  The new line sentinel does not occur in the
C                    remainder of the string. LPIECE is set to
C                    LWRAP or LENMSG+1-NEXTC whichever is less.
C
C.....LPIECE .EQ. 1  The new line sentinel starts at
C                    MESSG(NEXTC:NEXTC). LPIECE is effectively zero.
C                    Nothing is printed and NEXTC is incremented
C                    by 2. This takes care of the case when there 
C                    is a message with exactly 72 characters followed
C                    by a new line sentinel followed by more characters. 
C
C.....LPIECE .GT. LWRAP+1  Reduce LPIECE to LWRAP.
C
C.....Otherwise      This means that 2 .LE. LPIECE .LE. LWRAP+1
C                    Reset LPIECE = LPIECE-1.  
C                    Note that this properly handles the end case
C                    where LPIECE .EQ. LWRAP+1; that is, when the
C                    sentinel falls exactly at the end of a line.
C
      NEXTC = 0
C
C.....Loop point
C
  100 CONTINUE
C 
      LPIECE = INDEX(MESSG(NEXTC:LENMSG), NEWLIN)
      IF (LPIECE .EQ. 0) THEN
C
C........No new line sentinel found
C
         IDELTA = 0
         LPIECE = MIN(LWRAP, LENMSG+1-NEXTC)
         IF (LPIECE .LT. LENMSG+1-NEXTC) THEN
            DO 110 I=LPIECE+1,2,-1
               IF (MESSG(NEXTC+I-1:NEXTC+I-1) .EQ. ' ') THEN
                  LPIECE = I-1
                  IDELTA = 1
                  GOTO 120
               ENDIF
  110        CONTINUE
         ENDIF
  120    CBUFF(1:LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC = NEXTC + LPIECE + IDELTA
      ELSEIF (LPIECE .EQ. 1) THEN
C
C........a new line sentinel found at MESSG(NEXTC:NEXTC+1)
C........do not print a blank line
C
         NEXTC = NEXTC + 2
         GO TO 100
C
      ELSEIF (LPIECE .GT. LWRAP+1) THEN
C
C........lpiece has to be reduced to LWRAP
C
         IDELTA = 0
         LPIECE = LWRAP
         DO 140 I=LPIECE+1,2,-1
            IF (MESSG(NEXTC+I-1:NEXTC+I-1) .EQ. ' ') THEN
               LPIECE = I-1
               IDELTA = 1
               GOTO 150
            ENDIF
  140    CONTINUE
  150    CBUFF(1:LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC = NEXTC + LPIECE + IDELTA
C
      ELSE
C
C........now  2 .LE. LPIECE .LE. LWRAP+1
C........decrement LPIECE by one.
C
         LPIECE = LPIECE - 1
         CBUFF(1:LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC  = NEXTC + LPIECE + 2
      ENDIF
C
C.....Print
C
      DO 160 I = 1,KL
         WRITE(LOUT(I),200) CBUFF(1:LPIECE)
  160 CONTINUE
C
  200 FORMAT(A)
C
C.....Go back to the loop point if there is more to be printed
C
      IF (NEXTC .LE. LENMSG) GOTO 100
      RETURN
C
C.....End of MSGPRT
C
      END
C
C234567==1=========2=========3=========4=========5=========6=========7==
C
      SUBROUTINE LUNIT(KL, LOUT)
C
      INTEGER KL, LOUT(*)
C
C234567--x---------x---------x---------x---------x---------x---------x--
C
C   Generic function to return KL output-unit numbers for use by
C   by the message routine MSGPRT. For KL <= 0 or KL > 5 all
C   printout by MSGPRT will be suppressed.
C   As written the routine returns KL = 1 and LOUT(1) = 6. Modify 
C   the routine if more or different output units are to be used. 
C
C  Variables in the calling sequence
C  ---------------------------------
C   KL   I   OUT  Number of different output units to be used
C                 by MSGPRT. For KL <= 0 and KL > 5 all printout
C                 by MSGPRT is suppressed.
C   LOUT I   OUT  Array of dimension KL, 1 <= KL<= 5, for the
C                 KL output-units to be used by MSGPRT.
C
C234567--x---------x---------x---------x---------x---------x---------x--
C
      KL = 1
      LOUT(1) = 6
C
C.....End of LUNIT
C
      RETURN
      END
C
C234567==1=========2=========3=========4=========5=========6=========7==
C
      SUBROUTINE DMACH( EPMACH, SAFMIN )
C
      DOUBLE PRECISION EPMACH, SAFMIN
C
C234567--x---------x---------x---------x---------x---------x---------x--
C
C  To determine the relative machine precision EPMACH and the 
C  safe minimum SAFMIN.
C
C  Variables in the calling sequence
C  ---------------------------------
C  EPMACH D   OUT   The smallest floating point number such that
C                   1.0D0+EPS .NE. 1.0D0
C  SAFMIN D   OUT   SAFMIN = BETA**M where M is the largest in 
C                   magnitude negative integer M such that BETA**M 
C                   is positive and normalized
C
C234567--x---------x---------x---------x---------x---------x---------x--
C
C.....External function
C
      DOUBLE PRECISION ABS
C
C.....Parameters
C
      DOUBLE PRECISION ZER, ONE
      PARAMETER( ZER=0.0D0, ONE=1.0D0 )
C
C.....Local variables
C
      INTEGER I,IBETA,IT,N
      DOUBLE PRECISION A,B,BETA,BETAIN,C,D,T,Y,Z
C
C.....Variables to be saved
C
      LOGICAL DONE
      DOUBLE PRECISION EPMCH,SFMIN
      SAVE DONE,EPMCH,SFMIN
C
C.....Data statements
C
      DATA DONE/.FALSE./, EPMCH/0.0D0/, SFMIN/0.0D0/
C
C.......................Executable statements...........................
C
C.....If this is not the first call, return the previously
C.....computed values
C
      IF( DONE ) THEN
         EPMACH = EPMCH
         SAFMIN = SFMIN
         RETURN
      ENDIF
      DONE = .TRUE.
C
C.....Determine the base of the floating point representation
C
      A = ONE
   10 CONTINUE
         A = A + A
         B = A + ONE
         C = B - A
      IF (C-ONE .EQ. ZER) GOTO 10
      B = ONE
   20 CONTINUE
         B = B + B
         C = A + B
         IBETA = INT(C - A)
      IF (IBETA .EQ. 0) GOTO 20
      BETA = DBLE(IBETA)
C
C.....Determine the number of digits in the representation
C
      IT = 0
      B = ONE
   30 CONTINUE
         IT = IT + 1
         B = B * BETA
         C = B + ONE
         D = C - B
      IF (D-ONE .EQ. ZER) GOTO 30
C
C.....Determine the smallest positive floating-point number 
C.....such that 1.0D0 + EPS .NE. 1.0D0
C
      N = IT + 3
      BETAIN = ONE/BETA
      A = ONE
      DO 40 I = 1, N
         A = A * BETAIN
   40 CONTINUE
C
   50 CONTINUE
         B = ONE + A
         IF (B - ONE .NE. ZER) GOTO 60
         A = A * BETA
         GOTO 50
   60 CONTINUE 
      EPMCH = A
C
C.....Determine the largest in magnitude negative integer M 
C.....such that BETA**M is positive and normalized and set 
C.....SAFMIN = BETA**M.
C
      Z = BETAIN
      T = ONE + EPMCH
   70 CONTINUE
         Y = Z
         Z = Y * Y
         A = Z * ONE
         C = Z * T
         IF ((A+A .EQ. ZER) .OR. (ABS(Z) .GE. Y)) GOTO 80
         D = C * BETAIN
         IF (D*BETA .EQ. Z) GOTO 80
         GOTO 70
   80 CONTINUE
         SFMIN = Y
         Y = Y * BETAIN
         A = Y * ONE
         C = Y * T
         IF (((A+A) .EQ. ZER) .OR. (ABS(Y) .GE. SFMIN)) GOTO 90
         D = C * BETAIN
         IF ((D*BETA .NE. Y) .OR. (C .EQ. Y)) GOTO 80
         SFMIN = Y
   90 CONTINUE
C
      EPMACH = EPMCH
      SAFMIN = SFMIN
C
      RETURN
C
C.....End of DMACH
C
      END
C
C234567==1=========2=========3=========4=========5=========6=========7==
