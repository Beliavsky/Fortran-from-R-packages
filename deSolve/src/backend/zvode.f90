! Original authors: Peter N. Brown, Alan C. Hindmarsh,
!   Geore D. Byrne (see original author statement below)
!
!  Adapted for use in R package deSolve by the deSolve authors.
!

!DECK ZVODE
SUBROUTINE ZVODE (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK, &
& ISTATE, IOPT, ZWORK, LZW, RWORK, LRW, IWORK, LIW, &
& JAC, MF, RPAR, IPAR)
EXTERNAL F, JAC
COMPLEX(KIND=KIND(0.0d0)) Y, ZWORK
DOUBLE PRECISION T, TOUT, RTOL, ATOL, RWORK
INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LZW, LRW, IWORK, LIW, &
& MF, IPAR
DIMENSION Y(*), RTOL(*), ATOL(*), ZWORK(LZW), RWORK(LRW), &
& IWORK(LIW), RPAR(*), IPAR(*)
!-----------------------------------------------------------------------
! ZVODE: Variable-coefficient Ordinary Differential Equation solver,
! with fixed-leading-coefficient implementation.
! This version is in complex double precision.
!
! ZVODE solves the initial value problem for stiff or nonstiff
! systems of first order ODEs,
!     dy/dt = f(t,y) ,  or, in component form,
!     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(NEQ)) (i = 1,...,NEQ).
! Here the y vector is treated as complex.
! ZVODE is a package based on the EPISODE and EPISODEB packages, and
! on the ODEPACK user interface standard, with minor modifications.
!
! NOTE: When using ZVODE for a stiff system, it should only be used for
! the case in which the function f is analytic, that is, when each f(i)
! is an analytic function of each y(j).  Analyticity means that the
! partial derivative df(i)/dy(j) is a unique complex number, and this
! fact is critical in the way ZVODE solves the dense or banded linear
! systems that arise in the stiff case.  For a complex stiff ODE system
! in which f is not analytic, ZVODE is likely to have convergence
! failures, and for this problem one should instead use DVODE on the
! equivalent real system (in the real and imaginary parts of y).
!-----------------------------------------------------------------------
! Authors:
!               Peter N. Brown and Alan C. Hindmarsh
!               Center for Applied Scientific Computing
!               Lawrence Livermore National Laboratory
!               Livermore, CA 94551
! and
!               George D. Byrne (Prof. Emeritus)
!               Illinois Institute of Technology
!               Chicago, IL 60616
!-----------------------------------------------------------------------
! For references, see DVODE.
!-----------------------------------------------------------------------
! Summary of usage.
!
! Communication between the user and the ZVODE package, for normal
! situations, is summarized here.  This summary describes only a subset
! of the full set of options available.  See the full description for
! details, including optional communication, nonstandard options,
! and instructions for special situations.  See also the example
! problem (with program and output) following this summary.
!
! A. First provide a subroutine of the form:
!           SUBROUTINE F (NEQ, T, Y, YDOT, RPAR, IPAR)
!           DOUBLE COMPLEX Y(NEQ), YDOT(NEQ)
!           DOUBLE PRECISION T
! which supplies the vector function f by loading YDOT(i) with f(i).
!
! B. Next determine (or guess) whether or not the problem is stiff.
! Stiffness occurs when the Jacobian matrix df/dy has an eigenvalue
! whose real part is negative and large in magnitude, compared to the
! reciprocal of the t span of interest.  If the problem is nonstiff,
! use a method flag MF = 10.  If it is stiff, there are four standard
! choices for MF (21, 22, 24, 25), and ZVODE requires the Jacobian
! matrix in some form.  In these cases (MF .gt. 0), ZVODE will use a
! saved copy of the Jacobian matrix.  If this is undesirable because of
! storage limitations, set MF to the corresponding negative value
! (-21, -22, -24, -25).  (See full description of MF below.)
! The Jacobian matrix is regarded either as full (MF = 21 or 22),
! or banded (MF = 24 or 25).  In the banded case, ZVODE requires two
! half-bandwidth parameters ML and MU.  These are, respectively, the
! widths of the lower and upper parts of the band, excluding the main
! diagonal.  Thus the band consists of the locations (i,j) with
! i-ML .le. j .le. i+MU, and the full bandwidth is ML+MU+1.
!
! C. If the problem is stiff, you are encouraged to supply the Jacobian
! directly (MF = 21 or 24), but if this is not feasible, ZVODE will
! compute it internally by difference quotients (MF = 22 or 25).
! If you are supplying the Jacobian, provide a subroutine of the form:
!           SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD, RPAR, IPAR)
!           DOUBLE COMPLEX Y(NEQ), PD(NROWPD,NEQ)
!           DOUBLE PRECISION T
! which supplies df/dy by loading PD as follows:
!     For a full Jacobian (MF = 21), load PD(i,j) with df(i)/dy(j),
! the partial derivative of f(i) with respect to y(j).  (Ignore the
! ML and MU arguments in this case.)
!     For a banded Jacobian (MF = 24), load PD(i-j+MU+1,j) with
! df(i)/dy(j), i.e. load the diagonal lines of df/dy into the rows of
! PD from the top down.
!     In either case, only nonzero elements need be loaded.
!
! D. Write a main program which calls subroutine ZVODE once for
! each point at which answers are desired.  This should also provide
! for possible use of logical unit 6 for output of error messages
! by ZVODE.  On the first call to ZVODE, supply arguments as follows:
! F      = Name of subroutine for right-hand side vector f.
!          This name must be declared external in calling program.
! NEQ    = Number of first order ODEs.
! Y      = Double complex array of initial values, of length NEQ.
! T      = The initial value of the independent variable.
! TOUT   = First point where output is desired (.ne. T).
! ITOL   = 1 or 2 according as ATOL (below) is a scalar or array.
! RTOL   = Relative tolerance parameter (scalar).
! ATOL   = Absolute tolerance parameter (scalar or array).
!          The estimated local error in Y(i) will be controlled so as
!          to be roughly less (in magnitude) than
!             EWT(i) = RTOL*abs(Y(i)) + ATOL     if ITOL = 1, or
!             EWT(i) = RTOL*abs(Y(i)) + ATOL(i)  if ITOL = 2.
!          Thus the local error test passes if, in each component,
!          either the absolute error is less than ATOL (or ATOL(i)),
!          or the relative error is less than RTOL.
!          Use RTOL = 0.0 for pure absolute error control, and
!          use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error
!          control.  Caution: Actual (global) errors may exceed these
!          local tolerances, so choose them conservatively.
! ITASK  = 1 for normal computation of output values of Y at t = TOUT.
! ISTATE = Integer flag (input and output).  Set ISTATE = 1.
! IOPT   = 0 to indicate no optional input used.
! ZWORK  = Double precision complex work array of length at least:
!             15*NEQ                      for MF = 10,
!             8*NEQ + 2*NEQ**2            for MF = 21 or 22,
!             10*NEQ + (3*ML + 2*MU)*NEQ  for MF = 24 or 25.
! LZW    = Declared length of ZWORK (in user's DIMENSION statement).
! RWORK  = Real work array of length at least 20 + NEQ.
! LRW    = Declared length of RWORK (in user's DIMENSION statement).
! IWORK  = Integer work array of length at least:
!             30        for MF = 10,
!             30 + NEQ  for MF = 21, 22, 24, or 25.
!          If MF = 24 or 25, input in IWORK(1),IWORK(2) the lower
!          and upper half-bandwidths ML,MU.
! LIW    = Declared length of IWORK (in user's DIMENSION statement).
! JAC    = Name of subroutine for Jacobian matrix (MF = 21 or 24).
!          If used, this name must be declared external in calling
!          program.  If not used, pass a dummy name.
! MF     = Method flag.  Standard values are:
!          10 for nonstiff (Adams) method, no Jacobian used.
!          21 for stiff (BDF) method, user-supplied full Jacobian.
!          22 for stiff method, internally generated full Jacobian.
!          24 for stiff method, user-supplied banded Jacobian.
!          25 for stiff method, internally generated banded Jacobian.
! RPAR   = user-defined real or complex array passed to F and JAC.
! IPAR   = user-defined integer array passed to F and JAC.
! Note that the main program must declare arrays Y, ZWORK, RWORK, IWORK,
! and possibly ATOL, RPAR, and IPAR.  RPAR may be declared REAL, DOUBLE,
! COMPLEX, or DOUBLE COMPLEX, depending on the user's needs.
!
! E. The output from the first call (or any call) is:
!      Y = Array of computed values of y(t) vector.
!      T = Corresponding value of independent variable (normally TOUT).
! ISTATE = 2  if ZVODE was successful, negative otherwise.
!          -1 means excess work done on this call. (Perhaps wrong MF.)
!          -2 means excess accuracy requested. (Tolerances too small.)
!          -3 means illegal input detected. (See printed message.)
!          -4 means repeated error test failures. (Check all input.)
!          -5 means repeated convergence failures. (Perhaps bad
!             Jacobian supplied or wrong choice of MF or tolerances.)
!          -6 means error weight became zero during problem. (Solution
!             component i vanished, and ATOL or ATOL(i) = 0.)
!
! F. To continue the integration after a successful return, simply
! reset TOUT and call ZVODE again.  No other parameters need be reset.
!
!-----------------------------------------------------------------------
! EXAMPLE PROBLEM
!
! The program below uses ZVODE to solve the following system of 2 ODEs:
! dw/dt = -i*w*w*z, dz/dt = i*z; w(0) = 1/2.1, z(0) = 1; t = 0 to 2*pi.
! Solution: w = 1/(z + 1.1), z = exp(it).  As z traces the unit circle,
! w traces a circle of radius 10/2.1 with center at 11/2.1.
! For convenience, Main passes RPAR = (imaginary unit i) to FEX and JEX.
!
!     EXTERNAL FEX, JEX
!     DOUBLE COMPLEX Y(2), ZWORK(24), RPAR, WTRU, ERR
!     DOUBLE PRECISION ABERR, AEMAX, ATOL, RTOL, RWORK(22), T, TOUT
!     DIMENSION IWORK(32)
!     NEQ = 2
!     Y(1) = 1.0D0/2.1D0
!     Y(2) = 1.0D0
!     T = 0.0D0
!     DTOUT = 0.1570796326794896D0
!     TOUT = DTOUT
!     ITOL = 1
!     RTOL = 1.D-9
!     ATOL = 1.D-8
!     ITASK = 1
!     ISTATE = 1
!     IOPT = 0
!     LZW = 24
!     LRW = 22
!     LIW = 32
!     MF = 21
!     RPAR = DCMPLX(0.0D0,1.0D0)
!     AEMAX = 0.0D0
!     WRITE(6,10)
! 10  FORMAT('   t',11X,'w',26X,'z')
!     DO 40 IOUT = 1,40
!       CALL ZVODE(FEX,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,IOPT,
!    1             ZWORK,LZW,RWORK,LRW,IWORK,LIW,JEX,MF,RPAR,IPAR)
!       WTRU = 1.0D0/DCMPLX(COS(T) + 1.1D0, SIN(T))
!       ERR = Y(1) - WTRU
!       ABERR = ABS(DREAL(ERR)) + ABS(DIMAG(ERR))
!       AEMAX = MAX(AEMAX,ABERR)
!       WRITE(6,20) T, DREAL(Y(1)),DIMAG(Y(1)), DREAL(Y(2)),DIMAG(Y(2))
! 20    FORMAT(F9.5,2X,2F12.7,3X,2F12.7)
!       IF (ISTATE .LT. 0) THEN
!         WRITE(6,30) ISTATE
! 30      FORMAT(//'***** Error halt.  ISTATE =',I3)
!         STOP
!         ENDIF
! 40    TOUT = TOUT + DTOUT
!     WRITE(6,50) IWORK(11), IWORK(12), IWORK(13), IWORK(20),
!    1            IWORK(21), IWORK(22), IWORK(23), AEMAX
! 50  FORMAT(/' No. steps =',I4,'   No. f-s =',I5,
!    1        '   No. J-s =',I4,'   No. LU-s =',I4/
!    2        ' No. nonlinear iterations =',I4/
!    3        ' No. nonlinear convergence failures =',I4/
!    4        ' No. error test failures =',I4/
!    5        ' Max. abs. error in w =',D10.2)
!     STOP
!     END
!
!     SUBROUTINE FEX (NEQ, T, Y, YDOT, RPAR, IPAR)
!     DOUBLE COMPLEX Y(NEQ), YDOT(NEQ), RPAR
!     DOUBLE PRECISION T
!     YDOT(1) = -RPAR*Y(1)*Y(1)*Y(2)
!     YDOT(2) = RPAR*Y(2)
!     RETURN
!     END
!
!     SUBROUTINE JEX (NEQ, T, Y, ML, MU, PD, NRPD, RPAR, IPAR)
!     DOUBLE COMPLEX Y(NEQ), PD(NRPD,NEQ), RPAR
!     DOUBLE PRECISION T
!     PD(1,1) = -2.0D0*RPAR*Y(1)*Y(2)
!     PD(1,2) = -RPAR*Y(1)*Y(1)
!     PD(2,2) = RPAR
!     RETURN
!     END
!
! The output of this example program is as follows:
!
!   t           w                          z
!  0.15708     0.4763242  -0.0356919      0.9876884   0.1564345
!  0.31416     0.4767322  -0.0718256      0.9510565   0.3090170
!  0.47124     0.4774351  -0.1088651      0.8910065   0.4539906
!  0.62832     0.4784699  -0.1473206      0.8090170   0.5877853
!  0.78540     0.4798943  -0.1877789      0.7071067   0.7071069
!  0.94248     0.4817938  -0.2309414      0.5877852   0.8090171
!  1.09956     0.4842934  -0.2776778      0.4539904   0.8910066
!  1.25664     0.4875766  -0.3291039      0.3090169   0.9510566
!  1.41372     0.4919177  -0.3866987      0.1564343   0.9876884
!  1.57080     0.4977376  -0.4524889     -0.0000001   1.0000000
!  1.72788     0.5057044  -0.5293524     -0.1564346   0.9876883
!  1.88496     0.5169274  -0.6215400     -0.3090171   0.9510565
!  2.04204     0.5333540  -0.7356275     -0.4539906   0.8910065
!  2.19911     0.5586542  -0.8823669     -0.5877854   0.8090169
!  2.35619     0.6004188  -1.0806013     -0.7071069   0.7071067
!  2.51327     0.6764486  -1.3664281     -0.8090171   0.5877851
!  2.67035     0.8366909  -1.8175245     -0.8910066   0.4539904
!  2.82743     1.2657121  -2.6260146     -0.9510566   0.3090168
!  2.98451     3.0284506  -4.2182180     -0.9876884   0.1564343
!  3.14159    10.0000699   0.0000663     -1.0000000  -0.0000002
!  3.29867     3.0284170   4.2182053     -0.9876883  -0.1564346
!  3.45575     1.2657041   2.6260067     -0.9510565  -0.3090172
!  3.61283     0.8366878   1.8175205     -0.8910064  -0.4539907
!  3.76991     0.6764469   1.3664259     -0.8090169  -0.5877854
!  3.92699     0.6004178   1.0806000     -0.7071066  -0.7071069
!  4.08407     0.5586535   0.8823662     -0.5877851  -0.8090171
!  4.24115     0.5333535   0.7356271     -0.4539903  -0.8910066
!  4.39823     0.5169271   0.6215398     -0.3090168  -0.9510566
!  4.55531     0.5057041   0.5293523     -0.1564343  -0.9876884
!  4.71239     0.4977374   0.4524890      0.0000002  -1.0000000
!  4.86947     0.4919176   0.3866988      0.1564347  -0.9876883
!  5.02655     0.4875765   0.3291040      0.3090172  -0.9510564
!  5.18363     0.4842934   0.2776780      0.4539907  -0.8910064
!  5.34071     0.4817939   0.2309415      0.5877854  -0.8090169
!  5.49779     0.4798944   0.1877791      0.7071069  -0.7071066
!  5.65487     0.4784700   0.1473208      0.8090171  -0.5877850
!  5.81195     0.4774352   0.1088652      0.8910066  -0.4539903
!  5.96903     0.4767324   0.0718257      0.9510566  -0.3090168
!  6.12611     0.4763244   0.0356920      0.9876884  -0.1564342
!  6.28319     0.4761907   0.0000000      1.0000000   0.0000003
!
! No. steps = 542   No. f-s =  610   No. J-s =  10   No. LU-s =  47
! No. nonlinear iterations = 607
! No. nonlinear convergence failures =   0
! No. error test failures =  13
! Max. abs. error in w =  0.13E-03
!
!-----------------------------------------------------------------------
! Full description of user interface to ZVODE.
!
! The user interface to ZVODE consists of the following parts.
!
! i.   The call sequence to subroutine ZVODE, which is a driver
!      routine for the solver.  This includes descriptions of both
!      the call sequence arguments and of user-supplied routines.
!      Following these descriptions is
!        * a description of optional input available through the
!          call sequence,
!        * a description of optional output (in the work arrays), and
!        * instructions for interrupting and restarting a solution.
!
! ii.  Descriptions of other routines in the ZVODE package that may be
!      (optionally) called by the user.  These provide the ability to
!      alter error message handling, save and restore the internal
!      COMMON, and obtain specified derivatives of the solution y(t).
!
! iii. Descriptions of COMMON blocks to be declared in overlay
!      or similar environments.
!
! iv.  Description of two routines in the ZVODE package, either of
!      which the user may replace with his own version, if desired.
!      these relate to the measurement of errors.
!
!-----------------------------------------------------------------------
! Part i.  Call Sequence.
!
! The call sequence parameters used for input only are
!     F, NEQ, TOUT, ITOL, RTOL, ATOL, ITASK, IOPT, LRW, LIW, JAC, MF,
! and those used for both input and output are
!     Y, T, ISTATE.
! The work arrays ZWORK, RWORK, and IWORK are also used for conditional
! and optional input and optional output.  (The term output here refers
! to the return from subroutine ZVODE to the user's calling program.)
!
! The legality of input parameters will be thoroughly checked on the
! initial call for the problem, but not checked thereafter unless a
! change in input parameters is flagged by ISTATE = 3 in the input.
!
! The descriptions of the call arguments are as follows.
!
! F      = The name of the user-supplied subroutine defining the
!          ODE system.  The system must be put in the first-order
!          form dy/dt = f(t,y), where f is a vector-valued function
!          of the scalar t and the vector y.  Subroutine F is to
!          compute the function f.  It is to have the form
!               SUBROUTINE F (NEQ, T, Y, YDOT, RPAR, IPAR)
!               DOUBLE COMPLEX Y(NEQ), YDOT(NEQ)
!               DOUBLE PRECISION T
!          where NEQ, T, and Y are input, and the array YDOT = f(t,y)
!          is output.  Y and YDOT are double complex arrays of length
!          NEQ.  Subroutine F should not alter Y(1),...,Y(NEQ).
!          F must be declared EXTERNAL in the calling program.
!
!          Subroutine F may access user-defined real/complex and
!          integer work arrays RPAR and IPAR, which are to be
!          dimensioned in the calling program.
!
!          If quantities computed in the F routine are needed
!          externally to ZVODE, an extra call to F should be made
!          for this purpose, for consistent and accurate results.
!          If only the derivative dy/dt is needed, use ZVINDY instead.
!
! NEQ    = The size of the ODE system (number of first order
!          ordinary differential equations).  Used only for input.
!          NEQ may not be increased during the problem, but
!          can be decreased (with ISTATE = 3 in the input).
!
! Y      = A double precision complex array for the vector of dependent
!          variables, of length NEQ or more.  Used for both input and
!          output on the first call (ISTATE = 1), and only for output
!          on other calls.  On the first call, Y must contain the
!          vector of initial values.  In the output, Y contains the
!          computed solution evaluated at T.  If desired, the Y array
!          may be used for other purposes between calls to the solver.
!
!          This array is passed as the Y argument in all calls to
!          F and JAC.
!
! T      = The independent variable.  In the input, T is used only on
!          the first call, as the initial point of the integration.
!          In the output, after each call, T is the value at which a
!          computed solution Y is evaluated (usually the same as TOUT).
!          On an error return, T is the farthest point reached.
!
! TOUT   = The next value of t at which a computed solution is desired.
!          Used only for input.
!
!          When starting the problem (ISTATE = 1), TOUT may be equal
!          to T for one call, then should .ne. T for the next call.
!          For the initial T, an input value of TOUT .ne. T is used
!          in order to determine the direction of the integration
!          (i.e. the algebraic sign of the step sizes) and the rough
!          scale of the problem.  Integration in either direction
!          (forward or backward in t) is permitted.
!
!          If ITASK = 2 or 5 (one-step modes), TOUT is ignored after
!          the first call (i.e. the first call with TOUT .ne. T).
!          Otherwise, TOUT is required on every call.
!
!          If ITASK = 1, 3, or 4, the values of TOUT need not be
!          monotone, but a value of TOUT which backs up is limited
!          to the current internal t interval, whose endpoints are
!          TCUR - HU and TCUR.  (See optional output, below, for
!          TCUR and HU.)
!
! ITOL   = An indicator for the type of error control.  See
!          description below under ATOL.  Used only for input.
!
! RTOL   = A relative error tolerance parameter, either a scalar or
!          an array of length NEQ.  See description below under ATOL.
!          Input only.
!
! ATOL   = An absolute error tolerance parameter, either a scalar or
!          an array of length NEQ.  Input only.
!
!          The input parameters ITOL, RTOL, and ATOL determine
!          the error control performed by the solver.  The solver will
!          control the vector e = (e(i)) of estimated local errors
!          in Y, according to an inequality of the form
!                      rms-norm of ( e(i)/EWT(i) )   .le.   1,
!          where       EWT(i) = RTOL(i)*abs(Y(i)) + ATOL(i),
!          and the rms-norm (root-mean-square norm) here is
!          rms-norm(v) = sqrt(sum v(i)**2 / NEQ).  Here EWT = (EWT(i))
!          is a vector of weights which must always be positive, and
!          the values of RTOL and ATOL should all be non-negative.
!          The following table gives the types (scalar/array) of
!          RTOL and ATOL, and the corresponding form of EWT(i).
!
!             ITOL    RTOL       ATOL          EWT(i)
!              1     scalar     scalar     RTOL*ABS(Y(i)) + ATOL
!              2     scalar     array      RTOL*ABS(Y(i)) + ATOL(i)
!              3     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL
!              4     array      array      RTOL(i)*ABS(Y(i)) + ATOL(i)
!
!          When either of these parameters is a scalar, it need not
!          be dimensioned in the user's calling program.
!
!          If none of the above choices (with ITOL, RTOL, and ATOL
!          fixed throughout the problem) is suitable, more general
!          error controls can be obtained by substituting
!          user-supplied routines for the setting of EWT and/or for
!          the norm calculation.  See Part iv below.
!
!          If global errors are to be estimated by making a repeated
!          run on the same problem with smaller tolerances, then all
!          components of RTOL and ATOL (i.e. of EWT) should be scaled
!          down uniformly.
!
! ITASK  = An index specifying the task to be performed.
!          Input only.  ITASK has the following values and meanings.
!          1  means normal computation of output values of y(t) at
!             t = TOUT (by overshooting and interpolating).
!          2  means take one step only and return.
!          3  means stop at the first internal mesh point at or
!             beyond t = TOUT and return.
!          4  means normal computation of output values of y(t) at
!             t = TOUT but without overshooting t = TCRIT.
!             TCRIT must be input as RWORK(1).  TCRIT may be equal to
!             or beyond TOUT, but not behind it in the direction of
!             integration.  This option is useful if the problem
!             has a singularity at or beyond t = TCRIT.
!          5  means take one step, without passing TCRIT, and return.
!             TCRIT must be input as RWORK(1).
!
!          Note:  If ITASK = 4 or 5 and the solver reaches TCRIT
!          (within roundoff), it will return T = TCRIT (exactly) to
!          indicate this (unless ITASK = 4 and TOUT comes before TCRIT,
!          in which case answers at T = TOUT are returned first).
!
! ISTATE = an index used for input and output to specify the
!          the state of the calculation.
!
!          In the input, the values of ISTATE are as follows.
!          1  means this is the first call for the problem
!             (initializations will be done).  See note below.
!          2  means this is not the first call, and the calculation
!             is to continue normally, with no change in any input
!             parameters except possibly TOUT and ITASK.
!             (If ITOL, RTOL, and/or ATOL are changed between calls
!             with ISTATE = 2, the new values will be used but not
!             tested for legality.)
!          3  means this is not the first call, and the
!             calculation is to continue normally, but with
!             a change in input parameters other than
!             TOUT and ITASK.  Changes are allowed in
!             NEQ, ITOL, RTOL, ATOL, IOPT, LRW, LIW, MF, ML, MU,
!             and any of the optional input except H0.
!             (See IWORK description for ML and MU.)
!          Note:  A preliminary call with TOUT = T is not counted
!          as a first call here, as no initialization or checking of
!          input is done.  (Such a call is sometimes useful to include
!          the initial conditions in the output.)
!          Thus the first call for which TOUT .ne. T requires
!          ISTATE = 1 in the input.
!
!          In the output, ISTATE has the following values and meanings.
!           1  means nothing was done, as TOUT was equal to T with
!              ISTATE = 1 in the input.
!           2  means the integration was performed successfully.
!          -1  means an excessive amount of work (more than MXSTEP
!              steps) was done on this call, before completing the
!              requested task, but the integration was otherwise
!              successful as far as T.  (MXSTEP is an optional input
!              and is normally 500.)  To continue, the user may
!              simply reset ISTATE to a value .gt. 1 and call again.
!              (The excess work step counter will be reset to 0.)
!              In addition, the user may increase MXSTEP to avoid
!              this error return.  (See optional input below.)
!          -2  means too much accuracy was requested for the precision
!              of the machine being used.  This was detected before
!              completing the requested task, but the integration
!              was successful as far as T.  To continue, the tolerance
!              parameters must be reset, and ISTATE must be set
!              to 3.  The optional output TOLSF may be used for this
!              purpose.  (Note: If this condition is detected before
!              taking any steps, then an illegal input return
!              (ISTATE = -3) occurs instead.)
!          -3  means illegal input was detected, before taking any
!              integration steps.  See written message for details.
!              Note:  If the solver detects an infinite loop of calls
!              to the solver with illegal input, it will cause
!              the run to stop.
!          -4  means there were repeated error test failures on
!              one attempted step, before completing the requested
!              task, but the integration was successful as far as T.
!              The problem may have a singularity, or the input
!              may be inappropriate.
!          -5  means there were repeated convergence test failures on
!              one attempted step, before completing the requested
!              task, but the integration was successful as far as T.
!              This may be caused by an inaccurate Jacobian matrix,
!              if one is being used.
!          -6  means EWT(i) became zero for some i during the
!              integration.  Pure relative error control (ATOL(i)=0.0)
!              was requested on a variable which has now vanished.
!              The integration was successful as far as T.
!
!          Note:  Since the normal output value of ISTATE is 2,
!          it does not need to be reset for normal continuation.
!          Also, since a negative input value of ISTATE will be
!          regarded as illegal, a negative output value requires the
!          user to change it, and possibly other input, before
!          calling the solver again.
!
! IOPT   = An integer flag to specify whether or not any optional
!          input is being used on this call.  Input only.
!          The optional input is listed separately below.
!          IOPT = 0 means no optional input is being used.
!                   Default values will be used in all cases.
!          IOPT = 1 means optional input is being used.
!
! ZWORK  = A double precision complex working array.
!          The length of ZWORK must be at least
!             NYH*(MAXORD + 1) + 2*NEQ + LWM    where
!          NYH    = the initial value of NEQ,
!          MAXORD = 12 (if METH = 1) or 5 (if METH = 2) (unless a
!                   smaller value is given as an optional input),
!          LWM = length of work space for matrix-related data:
!          LWM = 0             if MITER = 0,
!          LWM = 2*NEQ**2      if MITER = 1 or 2, and MF.gt.0,
!          LWM = NEQ**2        if MITER = 1 or 2, and MF.lt.0,
!          LWM = NEQ           if MITER = 3,
!          LWM = (3*ML+2*MU+2)*NEQ     if MITER = 4 or 5, and MF.gt.0,
!          LWM = (2*ML+MU+1)*NEQ       if MITER = 4 or 5, and MF.lt.0.
!          (See the MF description for METH and MITER.)
!          Thus if MAXORD has its default value and NEQ is constant,
!          this length is:
!             15*NEQ                    for MF = 10,
!             15*NEQ + 2*NEQ**2         for MF = 11 or 12,
!             15*NEQ + NEQ**2           for MF = -11 or -12,
!             16*NEQ                    for MF = 13,
!             17*NEQ + (3*ML+2*MU)*NEQ  for MF = 14 or 15,
!             16*NEQ + (2*ML+MU)*NEQ    for MF = -14 or -15,
!              8*NEQ                    for MF = 20,
!              8*NEQ + 2*NEQ**2         for MF = 21 or 22,
!              8*NEQ + NEQ**2           for MF = -21 or -22,
!              9*NEQ                    for MF = 23,
!             10*NEQ + (3*ML+2*MU)*NEQ  for MF = 24 or 25.
!              9*NEQ + (2*ML+MU)*NEQ    for MF = -24 or -25.
!
! LZW    = The length of the array ZWORK, as declared by the user.
!          (This will be checked by the solver.)
!
! RWORK  = A real working array (double precision).
!          The length of RWORK must be at least 20 + NEQ.
!          The first 20 words of RWORK are reserved for conditional
!          and optional input and optional output.
!
!          The following word in RWORK is a conditional input:
!            RWORK(1) = TCRIT = critical value of t which the solver
!                       is not to overshoot.  Required if ITASK is
!                       4 or 5, and ignored otherwise.  (See ITASK.)
!
! LRW    = The length of the array RWORK, as declared by the user.
!          (This will be checked by the solver.)
!
! IWORK  = An integer work array.  The length of IWORK must be at least
!             30        if MITER = 0 or 3 (MF = 10, 13, 20, 23), or
!             30 + NEQ  otherwise (abs(MF) = 11,12,14,15,21,22,24,25).
!          The first 30 words of IWORK are reserved for conditional and
!          optional input and optional output.
!
!          The following 2 words in IWORK are conditional input:
!            IWORK(1) = ML     These are the lower and upper
!            IWORK(2) = MU     half-bandwidths, respectively, of the
!                       banded Jacobian, excluding the main diagonal.
!                       The band is defined by the matrix locations
!                       (i,j) with i-ML .le. j .le. i+MU.  ML and MU
!                       must satisfy  0 .le.  ML,MU  .le. NEQ-1.
!                       These are required if MITER is 4 or 5, and
!                       ignored otherwise.  ML and MU may in fact be
!                       the band parameters for a matrix to which
!                       df/dy is only approximately equal.
!
! LIW    = the length of the array IWORK, as declared by the user.
!          (This will be checked by the solver.)
!
! Note:  The work arrays must not be altered between calls to ZVODE
! for the same problem, except possibly for the conditional and
! optional input, and except for the last 2*NEQ words of ZWORK and
! the last NEQ words of RWORK.  The latter space is used for internal
! scratch space, and so is available for use by the user outside ZVODE
! between calls, if desired (but not for use by F or JAC).
!
! JAC    = The name of the user-supplied routine (MITER = 1 or 4) to
!          compute the Jacobian matrix, df/dy, as a function of
!          the scalar t and the vector y.  It is to have the form
!               SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD,
!                               RPAR, IPAR)
!               DOUBLE COMPLEX Y(NEQ), PD(NROWPD,NEQ)
!               DOUBLE PRECISION T
!          where NEQ, T, Y, ML, MU, and NROWPD are input and the array
!          PD is to be loaded with partial derivatives (elements of the
!          Jacobian matrix) in the output.  PD must be given a first
!          dimension of NROWPD.  T and Y have the same meaning as in
!          Subroutine F.
!               In the full matrix case (MITER = 1), ML and MU are
!          ignored, and the Jacobian is to be loaded into PD in
!          columnwise manner, with df(i)/dy(j) loaded into PD(i,j).
!               In the band matrix case (MITER = 4), the elements
!          within the band are to be loaded into PD in columnwise
!          manner, with diagonal lines of df/dy loaded into the rows
!          of PD. Thus df(i)/dy(j) is to be loaded into PD(i-j+MU+1,j).
!          ML and MU are the half-bandwidth parameters. (See IWORK).
!          The locations in PD in the two triangular areas which
!          correspond to nonexistent matrix elements can be ignored
!          or loaded arbitrarily, as they are overwritten by ZVODE.
!               JAC need not provide df/dy exactly.  A crude
!          approximation (possibly with a smaller bandwidth) will do.
!               In either case, PD is preset to zero by the solver,
!          so that only the nonzero elements need be loaded by JAC.
!          Each call to JAC is preceded by a call to F with the same
!          arguments NEQ, T, and Y.  Thus to gain some efficiency,
!          intermediate quantities shared by both calculations may be
!          saved in a user COMMON block by F and not recomputed by JAC,
!          if desired.  Also, JAC may alter the Y array, if desired.
!          JAC must be declared external in the calling program.
!               Subroutine JAC may access user-defined real/complex and
!          integer work arrays, RPAR and IPAR, whose dimensions are set
!          by the user in the calling program.
!
! MF     = The method flag.  Used only for input.  The legal values of
!          MF are 10, 11, 12, 13, 14, 15, 20, 21, 22, 23, 24, 25,
!          -11, -12, -14, -15, -21, -22, -24, -25.
!          MF is a signed two-digit integer, MF = JSV*(10*METH + MITER).
!          JSV = SIGN(MF) indicates the Jacobian-saving strategy:
!            JSV =  1 means a copy of the Jacobian is saved for reuse
!                     in the corrector iteration algorithm.
!            JSV = -1 means a copy of the Jacobian is not saved
!                     (valid only for MITER = 1, 2, 4, or 5).
!          METH indicates the basic linear multistep method:
!            METH = 1 means the implicit Adams method.
!            METH = 2 means the method based on backward
!                     differentiation formulas (BDF-s).
!          MITER indicates the corrector iteration method:
!            MITER = 0 means functional iteration (no Jacobian matrix
!                      is involved).
!            MITER = 1 means chord iteration with a user-supplied
!                      full (NEQ by NEQ) Jacobian.
!            MITER = 2 means chord iteration with an internally
!                      generated (difference quotient) full Jacobian
!                      (using NEQ extra calls to F per df/dy value).
!            MITER = 3 means chord iteration with an internally
!                      generated diagonal Jacobian approximation
!                      (using 1 extra call to F per df/dy evaluation).
!            MITER = 4 means chord iteration with a user-supplied
!                      banded Jacobian.
!            MITER = 5 means chord iteration with an internally
!                      generated banded Jacobian (using ML+MU+1 extra
!                      calls to F per df/dy evaluation).
!          If MITER = 1 or 4, the user must supply a subroutine JAC
!          (the name is arbitrary) as described above under JAC.
!          For other values of MITER, a dummy argument can be used.
!
! RPAR     User-specified array used to communicate real or complex
!          parameters to user-supplied subroutines.  If RPAR is an
!          array, it must be dimensioned in the user's calling program;
!          if it is unused or it is a scalar, then it need not be
!          dimensioned.  The type of RPAR may be REAL, DOUBLE, COMPLEX,
!          or DOUBLE COMPLEX, depending on the user program's needs.
!          RPAR is not type-declared within ZVODE, but simply passed
!          (by address) to the user's F and JAC routines.
!
! IPAR     User-specified array used to communicate integer parameter
!          to user-supplied subroutines.  If IPAR is an array, it must
!          be dimensioned in the user's calling program.
!-----------------------------------------------------------------------
! Optional Input.
!
! The following is a list of the optional input provided for in the
! call sequence.  (See also Part ii.)  For each such input variable,
! this table lists its name as used in this documentation, its
! location in the call sequence, its meaning, and the default value.
! The use of any of this input requires IOPT = 1, and in that
! case all of this input is examined.  A value of zero for any
! of these optional input variables will cause the default value to be
! used.  Thus to use a subset of the optional input, simply preload
! locations 5 to 10 in RWORK and IWORK to 0.0 and 0 respectively, and
! then set those of interest to nonzero values.
!
! NAME    LOCATION      MEANING AND DEFAULT VALUE
!
! H0      RWORK(5)  The step size to be attempted on the first step.
!                   The default value is determined by the solver.
!
! HMAX    RWORK(6)  The maximum absolute step size allowed.
!                   The default value is infinite.
!
! HMIN    RWORK(7)  The minimum absolute step size allowed.
!                   The default value is 0.  (This lower bound is not
!                   enforced on the final step before reaching TCRIT
!                   when ITASK = 4 or 5.)
!
! MAXORD  IWORK(5)  The maximum order to be allowed.  The default
!                   value is 12 if METH = 1, and 5 if METH = 2.
!                   If MAXORD exceeds the default value, it will
!                   be reduced to the default value.
!                   If MAXORD is changed during the problem, it may
!                   cause the current order to be reduced.
!
! MXSTEP  IWORK(6)  Maximum number of (internally defined) steps
!                   allowed during one call to the solver.
!                   The default value is 500.
!
! MXHNIL  IWORK(7)  Maximum number of messages printed (per problem)
!                   warning that T + H = T on a step (H = step size).
!                   This must be positive to result in a non-default
!                   value.  The default value is 10.
!
!-----------------------------------------------------------------------
! Optional Output.
!
! As optional additional output from ZVODE, the variables listed
! below are quantities related to the performance of ZVODE
! which are available to the user.  These are communicated by way of
! the work arrays, but also have internal mnemonic names as shown.
! Except where stated otherwise, all of this output is defined
! on any successful return from ZVODE, and on any return with
! ISTATE = -1, -2, -4, -5, or -6.  On an illegal input return
! (ISTATE = -3), they will be unchanged from their existing values
! (if any), except possibly for TOLSF, LENZW, LENRW, and LENIW.
! On any error return, output relevant to the error will be defined,
! as noted below.
!
! NAME    LOCATION      MEANING
!
! HU      RWORK(11) The step size in t last used (successfully).
!
! HCUR    RWORK(12) The step size to be attempted on the next step.
!
! TCUR    RWORK(13) The current value of the independent variable
!                   which the solver has actually reached, i.e. the
!                   current internal mesh point in t.  In the output,
!                   TCUR will always be at least as far from the
!                   initial value of t as the current argument T,
!                   but may be farther (if interpolation was done).
!
! TOLSF   RWORK(14) A tolerance scale factor, greater than 1.0,
!                   computed when a request for too much accuracy was
!                   detected (ISTATE = -3 if detected at the start of
!                   the problem, ISTATE = -2 otherwise).  If ITOL is
!                   left unaltered but RTOL and ATOL are uniformly
!                   scaled up by a factor of TOLSF for the next call,
!                   then the solver is deemed likely to succeed.
!                   (The user may also ignore TOLSF and alter the
!                   tolerance parameters in any other way appropriate.)
!
! NST     IWORK(11) The number of steps taken for the problem so far.
!
! NFE     IWORK(12) The number of f evaluations for the problem so far.
!
! NJE     IWORK(13) The number of Jacobian evaluations so far.
!
! NQU     IWORK(14) The method order last used (successfully).
!
! NQCUR   IWORK(15) The order to be attempted on the next step.
!
! IMXER   IWORK(16) The index of the component of largest magnitude in
!                   the weighted local error vector ( e(i)/EWT(i) ),
!                   on an error return with ISTATE = -4 or -5.
!
! LENZW   IWORK(17) The length of ZWORK actually required.
!                   This is defined on normal returns and on an illegal
!                   input return for insufficient storage.
!
! LENRW   IWORK(18) The length of RWORK actually required.
!                   This is defined on normal returns and on an illegal
!                   input return for insufficient storage.
!
! LENIW   IWORK(19) The length of IWORK actually required.
!                   This is defined on normal returns and on an illegal
!                   input return for insufficient storage.
!
! NLU     IWORK(20) The number of matrix LU decompositions so far.
!
! NNI     IWORK(21) The number of nonlinear (Newton) iterations so far.
!
! NCFN    IWORK(22) The number of convergence failures of the nonlinear
!                   solver so far.
!
! NETF    IWORK(23) The number of error test failures of the integrator
!                   so far.
!
! The following two arrays are segments of the ZWORK array which
! may also be of interest to the user as optional output.
! For each array, the table below gives its internal name,
! its base address in ZWORK, and its description.
!
! NAME    BASE ADDRESS      DESCRIPTION
!
! YH       1             The Nordsieck history array, of size NYH by
!                        (NQCUR + 1), where NYH is the initial value
!                        of NEQ.  For j = 0,1,...,NQCUR, column j+1
!                        of YH contains HCUR**j/factorial(j) times
!                        the j-th derivative of the interpolating
!                        polynomial currently representing the
!                        solution, evaluated at t = TCUR.
!
! ACOR     LENZW-NEQ+1   Array of size NEQ used for the accumulated
!                        corrections on each step, scaled in the output
!                        to represent the estimated local error in Y
!                        on the last step.  This is the vector e in
!                        the description of the error control.  It is
!                        defined only on a successful return from ZVODE.
!
!-----------------------------------------------------------------------
! Interrupting and Restarting
!
! If the integration of a given problem by ZVODE is to be
! interrrupted and then later continued, such as when restarting
! an interrupted run or alternating between two or more ODE problems,
! the user should save, following the return from the last ZVODE call
! prior to the interruption, the contents of the call sequence
! variables and internal COMMON blocks, and later restore these
! values before the next ZVODE call for that problem.  To save
! and restore the COMMON blocks, use subroutine ZVSRCO, as
! described below in part ii.
!
! In addition, if non-default values for either LUN or MFLAG are
! desired, an extra call to XSETUN and/or XSETF should be made just
! before continuing the integration.  See Part ii below for details.
!
!-----------------------------------------------------------------------
! Part ii.  Other Routines Callable.
!
! The following are optional calls which the user may make to
! gain additional capabilities in conjunction with ZVODE.
! (The routines XSETUN and XSETF are designed to conform to the
! SLATEC error handling package.)
!
!     FORM OF CALL                  FUNCTION
!  CALL XSETUN(LUN)           Set the logical unit number, LUN, for
!                             output of messages from ZVODE, if
!                             the default is not desired.
!                             The default value of LUN is 6.
!
!  CALL XSETF(MFLAG)          Set a flag to control the printing of
!                             messages by ZVODE.
!                             MFLAG = 0 means do not print. (Danger:
!                             This risks losing valuable information.)
!                             MFLAG = 1 means print (the default).
!
!                             Either of the above calls may be made at
!                             any time and will take effect immediately.
!
!  CALL ZVSRCO(RSAV,ISAV,JOB) Saves and restores the contents of
!                             the internal COMMON blocks used by
!                             ZVODE. (See Part iii below.)
!                             RSAV must be a real array of length 51
!                             or more, and ISAV must be an integer
!                             array of length 40 or more.
!                             JOB=1 means save COMMON into RSAV/ISAV.
!                             JOB=2 means restore COMMON from RSAV/ISAV.
!                                ZVSRCO is useful if one is
!                             interrupting a run and restarting
!                             later, or alternating between two or
!                             more problems solved with ZVODE.
!
!  CALL ZVINDY(,,,,,)         Provide derivatives of y, of various
!        (See below.)         orders, at a specified point T, if
!                             desired.  It may be called only after
!                             a successful return from ZVODE.
!
! The detailed instructions for using ZVINDY are as follows.
! The form of the call is:
!
!  CALL ZVINDY (T, K, ZWORK, NYH, DKY, IFLAG)
!
! The input parameters are:
!
! T         = Value of independent variable where answers are desired
!             (normally the same as the T last returned by ZVODE).
!             For valid results, T must lie between TCUR - HU and TCUR.
!             (See optional output for TCUR and HU.)
! K         = Integer order of the derivative desired.  K must satisfy
!             0 .le. K .le. NQCUR, where NQCUR is the current order
!             (see optional output).  The capability corresponding
!             to K = 0, i.e. computing y(T), is already provided
!             by ZVODE directly.  Since NQCUR .ge. 1, the first
!             derivative dy/dt is always available with ZVINDY.
! ZWORK     = The history array YH.
! NYH       = Column length of YH, equal to the initial value of NEQ.
!
! The output parameters are:
!
! DKY       = A double complex array of length NEQ containing the
!             computed value of the K-th derivative of y(t).
! IFLAG     = Integer flag, returned as 0 if K and T were legal,
!             -1 if K was illegal, and -2 if T was illegal.
!             On an error return, a message is also written.
!-----------------------------------------------------------------------
! Part iii.  COMMON Blocks.
! If ZVODE is to be used in an overlay situation, the user
! must declare, in the primary overlay, the variables in:
!   (1) the call sequence to ZVODE,
!   (2) the two internal COMMON blocks
!         /ZVOD01/  of length  83  (50 double precision words
!                         followed by 33 integer words),
!         /ZVOD02/  of length  9  (1 double precision word
!                         followed by 8 integer words),
!
! If ZVODE is used on a system in which the contents of internal
! COMMON blocks are not preserved between calls, the user should
! declare the above two COMMON blocks in his calling program to insure
! that their contents are preserved.
!
!-----------------------------------------------------------------------
! Part iv.  Optionally Replaceable Solver Routines.
!
! Below are descriptions of two routines in the ZVODE package which
! relate to the measurement of errors.  Either routine can be
! replaced by a user-supplied version, if desired.  However, since such
! a replacement may have a major impact on performance, it should be
! done only when absolutely necessary, and only with great caution.
! (Note: The means by which the package version of a routine is
! superseded by the user's version may be system-dependent.)
!
! (a) ZEWSET.
! The following subroutine is called just before each internal
! integration step, and sets the array of error weights, EWT, as
! described under ITOL/RTOL/ATOL above:
!     SUBROUTINE ZEWSET (NEQ, ITOL, RTOL, ATOL, YCUR, EWT)
! where NEQ, ITOL, RTOL, and ATOL are as in the ZVODE call sequence,
! YCUR contains the current (double complex) dependent variable vector,
! and EWT is the array of weights set by ZEWSET.
!
! If the user supplies this subroutine, it must return in EWT(i)
! (i = 1,...,NEQ) a positive quantity suitable for comparison with
! errors in Y(i).  The EWT array returned by ZEWSET is passed to the
! ZVNORM routine (See below.), and also used by ZVODE in the computation
! of the optional output IMXER, the diagonal Jacobian approximation,
! and the increments for difference quotient Jacobians.
!
! In the user-supplied version of ZEWSET, it may be desirable to use
! the current values of derivatives of y.  Derivatives up to order NQ
! are available from the history array YH, described above under
! Optional Output.  In ZEWSET, YH is identical to the YCUR array,
! extended to NQ + 1 columns with a column length of NYH and scale
! factors of h**j/factorial(j).  On the first call for the problem,
! given by NST = 0, NQ is 1 and H is temporarily set to 1.0.
! NYH is the initial value of NEQ.  The quantities NQ, H, and NST
! can be obtained by including in ZEWSET the statements:
!     DOUBLE PRECISION RVOD, H, HU
!     COMMON /ZVOD01/ RVOD(50), IVOD(33)
!     COMMON /ZVOD02/ HU, NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
!     NQ = IVOD(28)
!     H = RVOD(21)
! Thus, for example, the current value of dy/dt can be obtained as
! YCUR(NYH+i)/H  (i=1,...,NEQ)  (and the division by H is
! unnecessary when NST = 0).
!
! (b) ZVNORM.
! The following is a real function routine which computes the weighted
! root-mean-square norm of a vector v:
!     D = ZVNORM (N, V, W)
! where:
!   N = the length of the vector,
!   V = double complex array of length N containing the vector,
!   W = real array of length N containing weights,
!   D = sqrt( (1/N) * sum(abs(V(i))*W(i))**2 ).
! ZVNORM is called with N = NEQ and with W(i) = 1.0/EWT(i), where
! EWT is as set by subroutine ZEWSET.
!
! If the user supplies this function, it should return a non-negative
! value of ZVNORM suitable for use in the error control in ZVODE.
! None of the arguments should be altered by ZVNORM.
! For example, a user-supplied ZVNORM routine might:
!   -substitute a max-norm of (V(i)*W(i)) for the rms-norm, or
!   -ignore some components of V in the norm, with the effect of
!    suppressing the error control on those components of Y.
!-----------------------------------------------------------------------
! REVISION HISTORY (YYYYMMDD)
!  20060517  DATE WRITTEN, modified from DVODE of 20020430.
!  20061227  Added note on use for analytic f.
!-----------------------------------------------------------------------
! Other Routines in the ZVODE Package.
!
! In addition to Subroutine ZVODE, the ZVODE package includes the
! following subroutines and function routines:
!  ZVHIN     computes an approximate step size for the initial step.
!  ZVINDY    computes an interpolated value of the y vector at t = TOUT.
!  ZVSTEP    is the core integrator, which does one step of the
!            integration and the associated error control.
!  ZVSET     sets all method coefficients and test constants.
!  ZVNLSD    solves the underlying nonlinear system -- the corrector.
!  ZVJAC     computes and preprocesses the Jacobian matrix J = df/dy
!            and the Newton iteration matrix P = I - (h/l1)*J.
!  ZVSOL     manages solution of linear system in chord iteration.
!  ZVJUST    adjusts the history array on a change of order.
!  ZEWSET    sets the error weight vector EWT before each step.
!  ZVNORM    computes the weighted r.m.s. norm of a vector.
!  ZABSSQ    computes the squared absolute value of a double complex z.
!  ZVSRCO    is a user-callable routine to save and restore
!            the contents of the internal COMMON blocks.
!  ZACOPY    is a routine to copy one two-dimensional array to another.
!  ZGEFA and ZGESL   are routines from LINPACK for solving full
!            systems of linear algebraic equations.
!  ZGBFA and ZGBSL   are routines from LINPACK for solving banded
!            linear systems.
!  DZSCAL    scales a double complex array by a double prec. scalar.
!  DZAXPY    adds a D.P. scalar times one complex vector to another.
!  ZCOPY     is a basic linear algebra module from the BLAS.
!  DUMACH    sets the unit roundoff of the machine.
!  XERRWD, XSETUN, XSETF, IXSAV, and IUMACH handle the printing of all
!            error messages and warnings.  XERRWD is machine-dependent.
! Note: ZVNORM, ZABSSQ, DUMACH, IXSAV, and IUMACH are function routines.
! All the others are subroutines.
! The intrinsic functions called with double precision complex arguments
! are: ABS, DREAL, and DIMAG.  All of these are expected to return
! double precision real values.
!
!-----------------------------------------------------------------------
!
! Type declarations for labeled COMMON block ZVOD01 --------------------
!
DOUBLE PRECISION ACNRM, CCMXJ, CONP, CRATE, DRC, EL, &
& ETA, ETAMAX, H, HMIN, HMXI, HNEW, HRL1, HSCAL, PRL1, &
& RC, RL1, SRUR, TAU, TQ, TN, UROUND
INTEGER ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH, &
& L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, &
& LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP, &
& N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ, &
& NSLP, NYH
!
! Type declarations for labeled COMMON block ZVOD02 --------------------
!
DOUBLE PRECISION HU
INTEGER NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
!
! Type declarations for local variables --------------------------------
!
EXTERNAL ZVNLSD
LOGICAL IHIT
DOUBLE PRECISION ATOLI, BIG, EWTI, FOUR, H0, HMAX, HMX, HUN, ONE, &
& PT2, RH, RTOLI, SIZE, TCRIT, TNEXT, TOLSF, TP, TWO, ZERO
INTEGER I, IER, IFLAG, IMXER, JCO, KGO, LENIW, LENJ, LENP, LENZW, &
& LENRW, LENWM, LF0, MBAND, MFA, ML, MORD, MU, MXHNL0, MXSTP0, &
& NITER, NSLAST
CHARACTER(LEN=80) MSG
!
! Type declaration for function subroutines called ---------------------
!
DOUBLE PRECISION DUMACH, ZVNORM
!
DIMENSION MORD(2)
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to ZVODE.
!-----------------------------------------------------------------------
SAVE MORD, MXHNL0, MXSTP0
SAVE ZERO, ONE, TWO, FOUR, PT2, HUN
!-----------------------------------------------------------------------
! The following internal COMMON blocks contain variables which are
! communicated between subroutines in the ZVODE package, or which are
! to be saved between calls to ZVODE.
! In each block, real variables precede integers.
! The block /ZVOD01/ appears in subroutines ZVODE, ZVINDY, ZVSTEP,
! ZVSET, ZVNLSD, ZVJAC, ZVSOL, ZVJUST and ZVSRCO.
! The block /ZVOD02/ appears in subroutines ZVODE, ZVINDY, ZVSTEP,
! ZVNLSD, ZVJAC, and ZVSRCO.
!
! The variables stored in the internal COMMON blocks are as follows:
!
! ACNRM  = Weighted r.m.s. norm of accumulated correction vectors.
! CCMXJ  = Threshhold on DRC for updating the Jacobian. (See DRC.)
! CONP   = The saved value of TQ(5).
! CRATE  = Estimated corrector convergence rate constant.
! DRC    = Relative change in H*RL1 since last ZVJAC call.
! EL     = Real array of integration coefficients.  See ZVSET.
! ETA    = Saved tentative ratio of new to old H.
! ETAMAX = Saved maximum value of ETA to be allowed.
! H      = The step size.
! HMIN   = The minimum absolute value of the step size H to be used.
! HMXI   = Inverse of the maximum absolute value of H to be used.
!          HMXI = 0.0 is allowed and corresponds to an infinite HMAX.
! HNEW   = The step size to be attempted on the next step.
! HRL1   = Saved value of H*RL1.
! HSCAL  = Stepsize in scaling of YH array.
! PRL1   = The saved value of RL1.
! RC     = Ratio of current H*RL1 to value on last ZVJAC call.
! RL1    = The reciprocal of the coefficient EL(1).
! SRUR   = Sqrt(UROUND), used in difference quotient algorithms.
! TAU    = Real vector of past NQ step sizes, length 13.
! TQ     = A real vector of length 5 in which ZVSET stores constants
!          used for the convergence test, the error test, and the
!          selection of H at a new order.
! TN     = The independent variable, updated on each step taken.
! UROUND = The machine unit roundoff.  The smallest positive real number
!          such that  1.0 + UROUND .ne. 1.0
! ICF    = Integer flag for convergence failure in ZVNLSD:
!            0 means no failures.
!            1 means convergence failure with out of date Jacobian
!                   (recoverable error).
!            2 means convergence failure with current Jacobian or
!                   singular matrix (unrecoverable error).
! INIT   = Saved integer flag indicating whether initialization of the
!          problem has been done (INIT = 1) or not.
! IPUP   = Saved flag to signal updating of Newton matrix.
! JCUR   = Output flag from ZVJAC showing Jacobian status:
!            JCUR = 0 means J is not current.
!            JCUR = 1 means J is current.
! JSTART = Integer flag used as input to ZVSTEP:
!            0  means perform the first step.
!            1  means take a new step continuing from the last.
!            -1 means take the next step with a new value of MAXORD,
!                  HMIN, HMXI, N, METH, MITER, and/or matrix parameters.
!          On return, ZVSTEP sets JSTART = 1.
! JSV    = Integer flag for Jacobian saving, = sign(MF).
! KFLAG  = A completion code from ZVSTEP with the following meanings:
!               0      the step was succesful.
!              -1      the requested error could not be achieved.
!              -2      corrector convergence could not be achieved.
!              -3, -4  fatal error in VNLS (can not occur here).
! KUTH   = Input flag to ZVSTEP showing whether H was reduced by the
!          driver.  KUTH = 1 if H was reduced, = 0 otherwise.
! L      = Integer variable, NQ + 1, current order plus one.
! LMAX   = MAXORD + 1 (used for dimensioning).
! LOCJS  = A pointer to the saved Jacobian, whose storage starts at
!          WM(LOCJS), if JSV = 1.
! LYH, LEWT, LACOR, LSAVF, LWM, LIWM = Saved integer pointers
!          to segments of ZWORK, RWORK, and IWORK.
! MAXORD = The maximum order of integration method to be allowed.
! METH/MITER = The method flags.  See MF.
! MSBJ   = The maximum number of steps between J evaluations, = 50.
! MXHNIL = Saved value of optional input MXHNIL.
! MXSTEP = Saved value of optional input MXSTEP.
! N      = The number of first-order ODEs, = NEQ.
! NEWH   = Saved integer to flag change of H.
! NEWQ   = The method order to be used on the next step.
! NHNIL  = Saved counter for occurrences of T + H = T.
! NQ     = Integer variable, the current integration method order.
! NQNYH  = Saved value of NQ*NYH.
! NQWAIT = A counter controlling the frequency of order changes.
!          An order change is about to be considered if NQWAIT = 1.
! NSLJ   = The number of steps taken as of the last Jacobian update.
! NSLP   = Saved value of NST as of last Newton matrix update.
! NYH    = Saved value of the initial value of NEQ.
! HU     = The step size in t last used.
! NCFN   = Number of nonlinear convergence failures so far.
! NETF   = The number of error test failures of the integrator so far.
! NFE    = The number of f evaluations for the problem so far.
! NJE    = The number of Jacobian evaluations so far.
! NLU    = The number of matrix LU decompositions so far.
! NNI    = Number of nonlinear iterations so far.
! NQU    = The method order last used.
! NST    = The number of steps taken for the problem so far.
!-----------------------------------------------------------------------
COMMON /ZVOD01/ ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13), ETA, &
& ETAMAX, H, HMIN, HMXI, HNEW, HRL1, HSCAL, PRL1, &
& RC, RL1, SRUR, TAU(13), TQ(5), TN, UROUND, &
& ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH, &
& L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, &
& LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP, &
& N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ, &
& NSLP, NYH
COMMON /ZVOD02/ HU, NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
!
DATA  MORD(1) /12/, MORD(2) /5/, MXSTP0 /500/, MXHNL0 /10/
DATA ZERO /0.0D0/, ONE /1.0D0/, TWO /2.0D0/, FOUR /4.0D0/, &
& PT2 /0.2D0/, HUN /100.0D0/
!-----------------------------------------------------------------------
! Block A.
! This code block is executed on every call.
! It tests ISTATE and ITASK for legality and branches appropriately.
! If ISTATE .gt. 1 but the flag INIT shows that initialization has
! not yet been done, an error return occurs.
! If ISTATE = 1 and TOUT = T, return immediately.
!-----------------------------------------------------------------------
! KARLINE: INITIALISED IHIT TO AVOID COMPILER WARNINGS - SHOULD HAVE NO EFFEXT
IHIT = .TRUE.
IF (ISTATE .LT. 1 .OR. ISTATE .GT. 3) GO TO 601
IF (ITASK .LT. 1 .OR. ITASK .GT. 5) GO TO 602
IF (ISTATE .EQ. 1) GO TO 10
IF (INIT .NE. 1) GO TO 603
IF (ISTATE .EQ. 2) GO TO 200
GO TO 20
10 INIT = 0
IF (TOUT .EQ. T) RETURN
!-----------------------------------------------------------------------
! Block B.
! The next code block is executed for the initial call (ISTATE = 1),
! or for a continuation call with parameter changes (ISTATE = 3).
! It contains checking of all input and various initializations.
!
! First check legality of the non-optional input NEQ, ITOL, IOPT,
! MF, ML, and MU.
!-----------------------------------------------------------------------
20 IF (NEQ .LE. 0) GO TO 604
IF (ISTATE .EQ. 1) GO TO 25
IF (NEQ .GT. N) GO TO 605
25 N = NEQ
IF (ITOL .LT. 1 .OR. ITOL .GT. 4) GO TO 606
IF (IOPT .LT. 0 .OR. IOPT .GT. 1) GO TO 607
JSV = SIGN(1,MF)
MFA = ABS(MF)
METH = MFA/10
MITER = MFA - 10*METH
IF (METH .LT. 1 .OR. METH .GT. 2) GO TO 608
IF (MITER .LT. 0 .OR. MITER .GT. 5) GO TO 608
IF (MITER .LE. 3) GO TO 30
ML = IWORK(1)
MU = IWORK(2)
IF (ML .LT. 0 .OR. ML .GE. N) GO TO 609
IF (MU .LT. 0 .OR. MU .GE. N) GO TO 610
30 CONTINUE
! Next process and check the optional input. ---------------------------
IF (IOPT .EQ. 1) GO TO 40
MAXORD = MORD(METH)
MXSTEP = MXSTP0
MXHNIL = MXHNL0
IF (ISTATE .EQ. 1) H0 = ZERO
HMXI = ZERO
HMIN = ZERO
GO TO 60
40 MAXORD = IWORK(5)
IF (MAXORD .LT. 0) GO TO 611
IF (MAXORD .EQ. 0) MAXORD = 100
MAXORD = MIN(MAXORD,MORD(METH))
MXSTEP = IWORK(6)
IF (MXSTEP .LT. 0) GO TO 612
IF (MXSTEP .EQ. 0) MXSTEP = MXSTP0
MXHNIL = IWORK(7)
IF (MXHNIL .LT. 0) GO TO 613
IF (MXHNIL .EQ. 0) MXHNIL = MXHNL0
IF (ISTATE .NE. 1) GO TO 50
H0 = RWORK(5)
IF ((TOUT - T)*H0 .LT. ZERO) GO TO 614
50 HMAX = RWORK(6)
IF (HMAX .LT. ZERO) GO TO 615
HMXI = ZERO
IF (HMAX .GT. ZERO) HMXI = ONE/HMAX
HMIN = RWORK(7)
IF (HMIN .LT. ZERO) GO TO 616
!-----------------------------------------------------------------------
! Set work array pointers and check lengths LZW, LRW, and LIW.
! Pointers to segments of ZWORK, RWORK, and IWORK are named by prefixing
! L to the name of the segment.  E.g., segment YH starts at ZWORK(LYH).
! Segments of ZWORK (in order) are denoted  YH, WM, SAVF, ACOR.
! Besides optional inputs/outputs, RWORK has only the segment EWT.
! Within WM, LOCJS is the location of the saved Jacobian (JSV .gt. 0).
!-----------------------------------------------------------------------
60 LYH = 1
IF (ISTATE .EQ. 1) NYH = N
LWM = LYH + (MAXORD + 1)*NYH
JCO = MAX(0,JSV)
IF (MITER .EQ. 0) LENWM = 0
IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
  LENWM = (1 + JCO)*N*N
  LOCJS = N*N + 1
ENDIF
IF (MITER .EQ. 3) LENWM = N
IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
  MBAND = ML + MU + 1
  LENP = (MBAND + ML)*N
  LENJ = MBAND*N
  LENWM = LENP + JCO*LENJ
  LOCJS = LENP + 1
  ENDIF
LSAVF = LWM + LENWM
LACOR = LSAVF + N
LENZW = LACOR + N - 1
IWORK(17) = LENZW
LEWT = 21
LENRW = 20 + N
IWORK(18) = LENRW
LIWM = 1
LENIW = 30 + N
IF (MITER .EQ. 0 .OR. MITER .EQ. 3) LENIW = 30
IWORK(19) = LENIW
IF (LENZW .GT. LZW) GO TO 628
IF (LENRW .GT. LRW) GO TO 617
IF (LENIW .GT. LIW) GO TO 618
! Check RTOL and ATOL for legality. ------------------------------------
RTOLI = RTOL(1)
ATOLI = ATOL(1)
DO 70 I = 1,N
  IF (ITOL .GE. 3) RTOLI = RTOL(I)
  IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
  IF (RTOLI .LT. ZERO) GO TO 619
  IF (ATOLI .LT. ZERO) GO TO 620
70 CONTINUE
IF (ISTATE .EQ. 1) GO TO 100
! If ISTATE = 3, set flag to signal parameter changes to ZVSTEP. -------
JSTART = -1
IF (NQ .LE. MAXORD) GO TO 200
! MAXORD was reduced below NQ.  Copy YH(*,MAXORD+2) into SAVF. ---------
CALL ZCOPY (N, ZWORK(LWM), 1, ZWORK(LSAVF), 1)
GO TO 200
!-----------------------------------------------------------------------
! Block C.
! The next block is for the initial call only (ISTATE = 1).
! It contains all remaining initializations, the initial call to F,
! and the calculation of the initial step size.
! The error weights in EWT are inverted after being loaded.
!-----------------------------------------------------------------------
100 UROUND = DUMACH()
TN = T
IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 110
TCRIT = RWORK(1)
IF ((TCRIT - TOUT)*(TOUT - T) .LT. ZERO) GO TO 625
IF (H0 .NE. ZERO .AND. (T + H0 - TCRIT)*H0 .GT. ZERO) &
& H0 = TCRIT - T
110 JSTART = 0
IF (MITER .GT. 0) SRUR = SQRT(UROUND)
CCMXJ = PT2
MSBJ = 50
NHNIL = 0
NST = 0
NJE = 0
NNI = 0
NCFN = 0
NETF = 0
NLU = 0
NSLJ = 0
NSLAST = 0
HU = ZERO
NQU = 0
! Initial call to F.  (LF0 points to YH(*,2).) -------------------------
LF0 = LYH + NYH
CALL F (N, T, Y, ZWORK(LF0), RPAR, IPAR)
NFE = 1
! Load the initial value vector in YH. ---------------------------------
CALL ZCOPY (N, Y, 1, ZWORK(LYH), 1)
! Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
NQ = 1
H = ONE
CALL ZEWSET (N, ITOL, RTOL, ATOL, ZWORK(LYH), RWORK(LEWT))
DO 120 I = 1,N
  IF (RWORK(I+LEWT-1) .LE. ZERO) GO TO 621
  RWORK(I+LEWT-1) = ONE/RWORK(I+LEWT-1)
120 CONTINUE
IF (H0 .NE. ZERO) GO TO 180
! Call ZVHIN to set initial step size H0 to be attempted. --------------
CALL ZVHIN (N, T, ZWORK(LYH), ZWORK(LF0), F, RPAR, IPAR, TOUT, &
& UROUND, RWORK(LEWT), ITOL, ATOL, Y, ZWORK(LACOR), H0, &
& NITER, IER)
NFE = NFE + NITER
IF (IER .NE. 0) GO TO 622
! Adjust H0 if necessary to meet HMAX bound. ---------------------------
180 RH = ABS(H0)*HMXI
IF (RH .GT. ONE) H0 = H0/RH
! Load H with H0 and scale YH(*,2) by H0. ------------------------------
H = H0
CALL DZSCAL (N, H0, ZWORK(LF0), 1)
GO TO 270
!-----------------------------------------------------------------------
! Block D.
! The next code block is for continuation calls only (ISTATE = 2 or 3)
! and is to check stop conditions before taking a step.
!-----------------------------------------------------------------------
200 NSLAST = NST
KUTH = 0
IF (ITASK .EQ. 1) THEN
  GOTO 210
ELSE IF (ITASK .EQ. 2) THEN
  GOTO 250
ELSE IF (ITASK .EQ. 3) THEN
  GOTO 220
ELSE IF (ITASK .EQ. 4) THEN
  GOTO 230
ELSE IF (ITASK .EQ. 5) THEN
  GOTO 240
ENDIF

!      GO TO (210, 250, 220, 230, 240), ITASK
210 IF ((TN - TOUT)*H .LT. ZERO) GO TO 250
CALL ZVINDY (TOUT, 0, ZWORK(LYH), NYH, Y, IFLAG)
IF (IFLAG .NE. 0) GO TO 627
T = TOUT
GO TO 420
220 TP = TN - HU*(ONE + HUN*UROUND)
IF ((TP - TOUT)*H .GT. ZERO) GO TO 623
IF ((TN - TOUT)*H .LT. ZERO) GO TO 250
GO TO 400
230 TCRIT = RWORK(1)
IF ((TN - TCRIT)*H .GT. ZERO) GO TO 624
IF ((TCRIT - TOUT)*H .LT. ZERO) GO TO 625
IF ((TN - TOUT)*H .LT. ZERO) GO TO 245
CALL ZVINDY (TOUT, 0, ZWORK(LYH), NYH, Y, IFLAG)
IF (IFLAG .NE. 0) GO TO 627
T = TOUT
GO TO 420
240 TCRIT = RWORK(1)
IF ((TN - TCRIT)*H .GT. ZERO) GO TO 624
245 HMX = ABS(TN) + ABS(H)
IHIT = ABS(TN - TCRIT) .LE. HUN*UROUND*HMX
IF (IHIT) GO TO 400
TNEXT = TN + HNEW*(ONE + FOUR*UROUND)
IF ((TNEXT - TCRIT)*H .LE. ZERO) GO TO 250
H = (TCRIT - TN)*(ONE - FOUR*UROUND)
KUTH = 1
!-----------------------------------------------------------------------
! Block E.
! The next block is normally executed for all calls and contains
! the call to the one-step core integrator ZVSTEP.
!
! This is a looping point for the integration steps.
!
! First check for too many steps being taken, update EWT (if not at
! start of problem), check for too much accuracy being requested, and
! check for H below the roundoff level in T.
!-----------------------------------------------------------------------
250 CONTINUE
IF ((NST-NSLAST) .GE. MXSTEP) GO TO 500
CALL ZEWSET (N, ITOL, RTOL, ATOL, ZWORK(LYH), RWORK(LEWT))
DO 260 I = 1,N
  IF (RWORK(I+LEWT-1) .LE. ZERO) GO TO 510
  RWORK(I+LEWT-1) = ONE/RWORK(I+LEWT-1)
260 CONTINUE
270 TOLSF = UROUND*ZVNORM (N, ZWORK(LYH), RWORK(LEWT))
IF (TOLSF .LE. ONE) GO TO 280
TOLSF = TOLSF*TWO
IF (NST .EQ. 0) GO TO 626
GO TO 520
280 IF ((TN + H) .NE. TN) GO TO 290
NHNIL = NHNIL + 1
IF (NHNIL .GT. MXHNIL) GO TO 290
MSG = 'ZVODE--  Warning: internal T (=R1) and H (=R2) are'
CALL XERRWD (MSG, 50, 101, 1, 0, 0, 0, 0, ZERO, ZERO)
MSG='      such that in the machine, T + H = T on the next step  '
CALL XERRWD (MSG, 60, 101, 1, 0, 0, 0, 0, ZERO, ZERO)
MSG = '      (H = step size). solver will continue anyway'
CALL XERRWD (MSG, 50, 101, 1, 0, 0, 0, 2, TN, H)
IF (NHNIL .LT. MXHNIL) GO TO 290
MSG = 'ZVODE--  Above warning has been issued I1 times.  '
CALL XERRWD (MSG, 50, 102, 1, 0, 0, 0, 0, ZERO, ZERO)
MSG = '      it will not be issued again for this problem'
CALL XERRWD (MSG, 50, 102, 1, 1, MXHNIL, 0, 0, ZERO, ZERO)
290 CONTINUE
!-----------------------------------------------------------------------
! CALL ZVSTEP (Y, YH, NYH, YH, EWT, SAVF, VSAV, ACOR,
!              WM, IWM, F, JAC, F, ZVNLSD, RPAR, IPAR)
!-----------------------------------------------------------------------
CALL ZVSTEP (Y, ZWORK(LYH), NYH, ZWORK(LYH), RWORK(LEWT), &
& ZWORK(LSAVF), Y, ZWORK(LACOR), ZWORK(LWM), IWORK(LIWM), &
& F, JAC, F, ZVNLSD, RPAR, IPAR)
KGO = 1 - KFLAG
! Branch on KFLAG.  Note: In this version, KFLAG can not be set to -3.
!  KFLAG .eq. 0,   -1,  -2
IF (KGO .EQ. 1) THEN
  GOTO 300
ELSE IF (KGO .EQ. 2) THEN
  GOTO 530
ELSE IF (KGO .EQ. 3) THEN
  GOTO 540
ENDIF
!      GO TO (300, 530, 540), KGO
!-----------------------------------------------------------------------
! Block F.
! The following block handles the case of a successful return from the
! core integrator (KFLAG = 0).  Test for stop conditions.
!-----------------------------------------------------------------------
300 INIT = 1
KUTH = 0
IF (ITASK .EQ. 1) THEN
  GOTO 310
ELSE IF (ITASK .EQ. 2) THEN
  GOTO 400
ELSE IF (ITASK .EQ. 3) THEN
  GOTO 330
ELSE IF (ITASK .EQ. 4) THEN
  GOTO 340
ELSE IF (ITASK .EQ. 5) THEN
  GOTO 350
ENDIF

!      GO TO (310, 400, 330, 340, 350), ITASK
! ITASK = 1.  If TOUT has been reached, interpolate. -------------------
310 IF ((TN - TOUT)*H .LT. ZERO) GO TO 250
CALL ZVINDY (TOUT, 0, ZWORK(LYH), NYH, Y, IFLAG)
T = TOUT
GO TO 420
! ITASK = 3.  Jump to exit if TOUT was reached. ------------------------
330 IF ((TN - TOUT)*H .GE. ZERO) GO TO 400
GO TO 250
! ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary.
340 IF ((TN - TOUT)*H .LT. ZERO) GO TO 345
CALL ZVINDY (TOUT, 0, ZWORK(LYH), NYH, Y, IFLAG)
T = TOUT
GO TO 420
345 HMX = ABS(TN) + ABS(H)
IHIT = ABS(TN - TCRIT) .LE. HUN*UROUND*HMX
IF (IHIT) GO TO 400
TNEXT = TN + HNEW*(ONE + FOUR*UROUND)
IF ((TNEXT - TCRIT)*H .LE. ZERO) GO TO 250
H = (TCRIT - TN)*(ONE - FOUR*UROUND)
KUTH = 1
GO TO 250
! ITASK = 5.  See if TCRIT was reached and jump to exit. ---------------
350 HMX = ABS(TN) + ABS(H)
IHIT = ABS(TN - TCRIT) .LE. HUN*UROUND*HMX
!-----------------------------------------------------------------------
! Block G.
! The following block handles all successful returns from ZVODE.
! If ITASK .ne. 1, Y is loaded from YH and T is set accordingly.
! ISTATE is set to 2, and the optional output is loaded into the work
! arrays before returning.
!-----------------------------------------------------------------------
400 CONTINUE
CALL ZCOPY (N, ZWORK(LYH), 1, Y, 1)
T = TN
IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 420
IF (IHIT) T = TCRIT
420 ISTATE = 2
RWORK(11) = HU
RWORK(12) = HNEW
RWORK(13) = TN
IWORK(11) = NST
IWORK(12) = NFE
IWORK(13) = NJE
IWORK(14) = NQU
IWORK(15) = NEWQ
IWORK(20) = NLU
IWORK(21) = NNI
IWORK(22) = NCFN
IWORK(23) = NETF
RETURN
!-----------------------------------------------------------------------
! Block H.
! The following block handles all unsuccessful returns other than
! those for illegal input.  First the error message routine is called.
! if there was an error test or convergence test failure, IMXER is set.
! Then Y is loaded from YH, and T is set to TN.
! The optional output is loaded into the work arrays before returning.
!-----------------------------------------------------------------------
! The maximum number of steps was taken before reaching TOUT. ----------
500 MSG = 'ZVODE--  At current T (=R1), MXSTEP (=I1) steps   '
CALL XERRWD (MSG, 50, 201, 1, 0, 0, 0, 0, ZERO, ZERO)
MSG = '      taken on this call before reaching TOUT     '
CALL XERRWD (MSG, 50, 201, 1, 1, MXSTEP, 0, 1, TN, ZERO)
ISTATE = -1
GO TO 580
! EWT(i) .le. 0.0 for some i (not at start of problem). ----------------
510 EWTI = RWORK(LEWT+I-1)
MSG = 'ZVODE--  At T (=R1), EWT(I1) has become R2 .LE. 0.'
CALL XERRWD (MSG, 50, 202, 1, 1, I, 0, 2, TN, EWTI)
ISTATE = -6
GO TO 580
! Too much accuracy requested for machine precision. -------------------
520 MSG = 'ZVODE--  At T (=R1), too much accuracy requested  '
CALL XERRWD (MSG, 50, 203, 1, 0, 0, 0, 0, ZERO, ZERO)
MSG = '      for precision of machine:   see TOLSF (=R2) '
CALL XERRWD (MSG, 50, 203, 1, 0, 0, 0, 2, TN, TOLSF)
RWORK(14) = TOLSF
ISTATE = -2
GO TO 580
! KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
530 MSG = 'ZVODE--  At T(=R1) and step size H(=R2), the error'
CALL XERRWD (MSG, 50, 204, 1, 0, 0, 0, 0, ZERO, ZERO)
MSG = '      test failed repeatedly or with abs(H) = HMIN'
CALL XERRWD (MSG, 50, 204, 1, 0, 0, 0, 2, TN, H)
ISTATE = -4
GO TO 560
! KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ----
540 MSG = 'ZVODE--  At T (=R1) and step size H (=R2), the    '
CALL XERRWD (MSG, 50, 205, 1, 0, 0, 0, 0, ZERO, ZERO)
MSG = '      corrector convergence failed repeatedly     '
CALL XERRWD (MSG, 50, 205, 1, 0, 0, 0, 0, ZERO, ZERO)
MSG = '      or with abs(H) = HMIN   '
CALL XERRWD (MSG, 30, 205, 1, 0, 0, 0, 2, TN, H)
ISTATE = -5
! Compute IMXER if relevant. -------------------------------------------
560 BIG = ZERO
IMXER = 1
DO 570 I = 1,N
  SIZE = ABS(ZWORK(I+LACOR-1))*RWORK(I+LEWT-1)
  IF (BIG .GE. SIZE) GO TO 570
  BIG = SIZE
  IMXER = I
570 CONTINUE
IWORK(16) = IMXER
! Set Y vector, T, and optional output. --------------------------------
580 CONTINUE
CALL ZCOPY (N, ZWORK(LYH), 1, Y, 1)
T = TN
RWORK(11) = HU
RWORK(12) = H
RWORK(13) = TN
IWORK(11) = NST
IWORK(12) = NFE
IWORK(13) = NJE
IWORK(14) = NQU
IWORK(15) = NQ
IWORK(20) = NLU
IWORK(21) = NNI
IWORK(22) = NCFN
IWORK(23) = NETF
RETURN
!-----------------------------------------------------------------------
! Block I.
! The following block handles all error returns due to illegal input
! (ISTATE = -3), as detected before calling the core integrator.
! First the error message routine is called.   If the illegal input
! is a negative ISTATE, the run is aborted (apparent infinite loop).
!-----------------------------------------------------------------------
601 MSG = 'ZVODE--  ISTATE (=I1) illegal '
CALL XERRWD (MSG, 30, 1, 1, 1, ISTATE, 0, 0, ZERO, ZERO)
IF (ISTATE .LT. 0) GO TO 800
GO TO 700
602 MSG = 'ZVODE--  ITASK (=I1) illegal  '
CALL XERRWD (MSG, 30, 2, 1, 1, ITASK, 0, 0, ZERO, ZERO)
GO TO 700
603 MSG='ZVODE--  ISTATE (=I1) .GT. 1 but ZVODE not initialized      '
CALL XERRWD (MSG, 60, 3, 1, 1, ISTATE, 0, 0, ZERO, ZERO)
GO TO 700
604 MSG = 'ZVODE--  NEQ (=I1) .LT. 1     '
CALL XERRWD (MSG, 30, 4, 1, 1, NEQ, 0, 0, ZERO, ZERO)
GO TO 700
605 MSG = 'ZVODE--  ISTATE = 3 and NEQ increased (I1 to I2)  '
CALL XERRWD (MSG, 50, 5, 1, 2, N, NEQ, 0, ZERO, ZERO)
GO TO 700
606 MSG = 'ZVODE--  ITOL (=I1) illegal   '
CALL XERRWD (MSG, 30, 6, 1, 1, ITOL, 0, 0, ZERO, ZERO)
GO TO 700
607 MSG = 'ZVODE--  IOPT (=I1) illegal   '
CALL XERRWD (MSG, 30, 7, 1, 1, IOPT, 0, 0, ZERO, ZERO)
GO TO 700
608 MSG = 'ZVODE--  MF (=I1) illegal     '
CALL XERRWD (MSG, 30, 8, 1, 1, MF, 0, 0, ZERO, ZERO)
GO TO 700
609 MSG = 'ZVODE--  ML (=I1) illegal:  .LT.0 or .GE.NEQ (=I2)'
CALL XERRWD (MSG, 50, 9, 1, 2, ML, NEQ, 0, ZERO, ZERO)
GO TO 700
610 MSG = 'ZVODE--  MU (=I1) illegal:  .LT.0 or .GE.NEQ (=I2)'
CALL XERRWD (MSG, 50, 10, 1, 2, MU, NEQ, 0, ZERO, ZERO)
GO TO 700
611 MSG = 'ZVODE--  MAXORD (=I1) .LT. 0  '
CALL XERRWD (MSG, 30, 11, 1, 1, MAXORD, 0, 0, ZERO, ZERO)
GO TO 700
612 MSG = 'ZVODE--  MXSTEP (=I1) .LT. 0  '
CALL XERRWD (MSG, 30, 12, 1, 1, MXSTEP, 0, 0, ZERO, ZERO)
GO TO 700
613 MSG = 'ZVODE--  MXHNIL (=I1) .LT. 0  '
CALL XERRWD (MSG, 30, 13, 1, 1, MXHNIL, 0, 0, ZERO, ZERO)
GO TO 700
614 MSG = 'ZVODE--  TOUT (=R1) behind T (=R2)      '
CALL XERRWD (MSG, 40, 14, 1, 0, 0, 0, 2, TOUT, T)
MSG = '      integration direction is given by H0 (=R1)  '
CALL XERRWD (MSG, 50, 14, 1, 0, 0, 0, 1, H0, ZERO)
GO TO 700
615 MSG = 'ZVODE--  HMAX (=R1) .LT. 0.0  '
CALL XERRWD (MSG, 30, 15, 1, 0, 0, 0, 1, HMAX, ZERO)
GO TO 700
616 MSG = 'ZVODE--  HMIN (=R1) .LT. 0.0  '
CALL XERRWD (MSG, 30, 16, 1, 0, 0, 0, 1, HMIN, ZERO)
GO TO 700
617 CONTINUE
MSG='ZVODE--  RWORK length needed, LENRW (=I1), exceeds LRW (=I2)'
CALL XERRWD (MSG, 60, 17, 1, 2, LENRW, LRW, 0, ZERO, ZERO)
GO TO 700
618 CONTINUE
MSG='ZVODE--  IWORK length needed, LENIW (=I1), exceeds LIW (=I2)'
CALL XERRWD (MSG, 60, 18, 1, 2, LENIW, LIW, 0, ZERO, ZERO)
GO TO 700
619 MSG = 'ZVODE--  RTOL(I1) is R1 .LT. 0.0        '
CALL XERRWD (MSG, 40, 19, 1, 1, I, 0, 1, RTOLI, ZERO)
GO TO 700
620 MSG = 'ZVODE--  ATOL(I1) is R1 .LT. 0.0        '
CALL XERRWD (MSG, 40, 20, 1, 1, I, 0, 1, ATOLI, ZERO)
GO TO 700
621 EWTI = RWORK(LEWT+I-1)
MSG = 'ZVODE--  EWT(I1) is R1 .LE. 0.0         '
CALL XERRWD (MSG, 40, 21, 1, 1, I, 0, 1, EWTI, ZERO)
GO TO 700
622 CONTINUE
MSG='ZVODE--  TOUT (=R1) too close to T(=R2) to start integration'
CALL XERRWD (MSG, 60, 22, 1, 0, 0, 0, 2, TOUT, T)
GO TO 700
623 CONTINUE
MSG='ZVODE--  ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2)  '
CALL XERRWD (MSG, 60, 23, 1, 1, ITASK, 0, 2, TOUT, TP)
GO TO 700
624 CONTINUE
MSG='ZVODE--  ITASK = 4 or 5 and TCRIT (=R1) behind TCUR (=R2)   '
CALL XERRWD (MSG, 60, 24, 1, 0, 0, 0, 2, TCRIT, TN)
GO TO 700
625 CONTINUE
MSG='ZVODE--  ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)   '
CALL XERRWD (MSG, 60, 25, 1, 0, 0, 0, 2, TCRIT, TOUT)
GO TO 700
626 MSG = 'ZVODE--  At start of problem, too much accuracy   '
CALL XERRWD (MSG, 50, 26, 1, 0, 0, 0, 0, ZERO, ZERO)
MSG='      requested for precision of machine:   see TOLSF (=R1) '
CALL XERRWD (MSG, 60, 26, 1, 0, 0, 0, 1, TOLSF, ZERO)
RWORK(14) = TOLSF
GO TO 700
627 MSG='ZVODE--  Trouble from ZVINDY.  ITASK = I1, TOUT = R1.       '
CALL XERRWD (MSG, 60, 27, 1, 1, ITASK, 0, 1, TOUT, ZERO)
GO TO 700
628 CONTINUE
MSG='ZVODE--  ZWORK length needed, LENZW (=I1), exceeds LZW (=I2)'
CALL XERRWD (MSG, 60, 17, 1, 2, LENZW, LZW, 0, ZERO, ZERO)
!
700 CONTINUE
ISTATE = -3
RETURN
!
800 MSG = 'ZVODE--  Run aborted:  apparent infinite loop     '
CALL XERRWD (MSG, 50, 303, 2, 0, 0, 0, 0, ZERO, ZERO)
RETURN
!----------------------- End of Subroutine ZVODE -----------------------
END
!DECK ZVHIN
SUBROUTINE ZVHIN (N, T0, Y0, YDOT, F, RPAR, IPAR, TOUT, UROUND, &
& EWT, ITOL, ATOL, Y, TEMP, H0, NITER, IER)
EXTERNAL F
COMPLEX(KIND=KIND(0.0d0)) Y0, YDOT, Y, TEMP
DOUBLE PRECISION T0, TOUT, UROUND, EWT, ATOL, H0
INTEGER N, IPAR, ITOL, NITER, IER
DIMENSION Y0(*), YDOT(*), EWT(*), ATOL(*), Y(*), &
& TEMP(*), RPAR(*), IPAR(*)
!-----------------------------------------------------------------------
! Call sequence input -- N, T0, Y0, YDOT, F, RPAR, IPAR, TOUT, UROUND,
!                        EWT, ITOL, ATOL, Y, TEMP
! Call sequence output -- H0, NITER, IER
! COMMON block variables accessed -- None
!
! Subroutines called by ZVHIN:  F
! Function routines called by ZVHIN: ZVNORM
!-----------------------------------------------------------------------
! This routine computes the step size, H0, to be attempted on the
! first step, when the user has not supplied a value for this.
!
! First we check that TOUT - T0 differs significantly from zero.  Then
! an iteration is done to approximate the initial second derivative
! and this is used to define h from w.r.m.s.norm(h**2 * yddot / 2) = 1.
! A bias factor of 1/2 is applied to the resulting h.
! The sign of H0 is inferred from the initial values of TOUT and T0.
!
! Communication with ZVHIN is done with the following variables:
!
! N      = Size of ODE system, input.
! T0     = Initial value of independent variable, input.
! Y0     = Vector of initial conditions, input.
! YDOT   = Vector of initial first derivatives, input.
! F      = Name of subroutine for right-hand side f(t,y), input.
! RPAR, IPAR = User's real/complex and integer work arrays.
! TOUT   = First output value of independent variable
! UROUND = Machine unit roundoff
! EWT, ITOL, ATOL = Error weights and tolerance parameters
!                   as described in the driver routine, input.
! Y, TEMP = Work arrays of length N.
! H0     = Step size to be attempted, output.
! NITER  = Number of iterations (and of f evaluations) to compute H0,
!          output.
! IER    = The error flag, returned with the value
!          IER = 0  if no trouble occurred, or
!          IER = -1 if TOUT and T0 are considered too close to proceed.
!-----------------------------------------------------------------------
!
! Type declarations for local variables --------------------------------
!
DOUBLE PRECISION AFI, ATOLI, DELYI, H, HALF, HG, HLB, HNEW, HRAT, &
& HUB, HUN, PT1, T1, TDIST, TROUND, TWO, YDDNRM
INTEGER I, ITER
!
! Type declaration for function subroutines called ---------------------
!
DOUBLE PRECISION ZVNORM
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
SAVE HALF, HUN, PT1, TWO
DATA HALF /0.5D0/, HUN /100.0D0/, PT1 /0.1D0/, TWO /2.0D0/
!
NITER = 0
TDIST = ABS(TOUT - T0)
TROUND = UROUND*MAX(ABS(T0),ABS(TOUT))
IF (TDIST .LT. TWO*TROUND) GO TO 100
!
! Set a lower bound on h based on the roundoff level in T0 and TOUT. ---
HLB = HUN*TROUND
! Set an upper bound on h based on TOUT-T0 and the initial Y and YDOT. -
HUB = PT1*TDIST
ATOLI = ATOL(1)
DO 10 I = 1, N
  IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
  DELYI = PT1*ABS(Y0(I)) + ATOLI
  AFI = ABS(YDOT(I))
  IF (AFI*HUB .GT. DELYI) HUB = DELYI/AFI
10 CONTINUE
!
! Set initial guess for h as geometric mean of upper and lower bounds. -
ITER = 0
HG = SQRT(HLB*HUB)
! If the bounds have crossed, exit with the mean value. ----------------
IF (HUB .LT. HLB) THEN
  H0 = HG
  GO TO 90
ENDIF
!
! Looping point for iteration. -----------------------------------------
50 CONTINUE
! Estimate the second derivative as a difference quotient in f. --------
H = SIGN (HG, TOUT - T0)
T1 = T0 + H
DO 60 I = 1, N
  Y(I) = Y0(I) + H*YDOT(I)
60 CONTINUE
CALL F (N, T1, Y, TEMP, RPAR, IPAR)
DO 70 I = 1, N
  TEMP(I) = (TEMP(I) - YDOT(I))/H
70 CONTINUE
YDDNRM = ZVNORM (N, TEMP, EWT)
! Get the corresponding new value of h. --------------------------------
IF (YDDNRM*HUB*HUB .GT. TWO) THEN
  HNEW = SQRT(TWO/YDDNRM)
ELSE
  HNEW = SQRT(HG*HUB)
ENDIF
ITER = ITER + 1
!-----------------------------------------------------------------------
! Test the stopping conditions.
! Stop if the new and previous h values differ by a factor of .lt. 2.
! Stop if four iterations have been done.  Also, stop with previous h
! if HNEW/HG .gt. 2 after first iteration, as this probably means that
! the second derivative value is bad because of cancellation error.
!-----------------------------------------------------------------------
IF (ITER .GE. 4) GO TO 80
HRAT = HNEW/HG
IF ( (HRAT .GT. HALF) .AND. (HRAT .LT. TWO) ) GO TO 80
IF ( (ITER .GE. 2) .AND. (HNEW .GT. TWO*HG) ) THEN
  HNEW = HG
  GO TO 80
ENDIF
HG = HNEW
GO TO 50
!
! Iteration done.  Apply bounds, bias factor, and sign.  Then exit. ----
80 H0 = HNEW*HALF
IF (H0 .LT. HLB) H0 = HLB
IF (H0 .GT. HUB) H0 = HUB
90 H0 = SIGN(H0, TOUT - T0)
NITER = ITER
IER = 0
RETURN
! Error return for TOUT - T0 too small. --------------------------------
100 IER = -1
RETURN
!----------------------- End of Subroutine ZVHIN -----------------------
END
!DECK ZVINDY
SUBROUTINE ZVINDY (T, K, YH, LDYH, DKY, IFLAG)
COMPLEX(KIND=KIND(0.0d0)) YH, DKY
DOUBLE PRECISION T
INTEGER K, LDYH, IFLAG
DIMENSION YH(LDYH,*), DKY(*)
!-----------------------------------------------------------------------
! Call sequence input -- T, K, YH, LDYH
! Call sequence output -- DKY, IFLAG
! COMMON block variables accessed:
!     /ZVOD01/ --  H, TN, UROUND, L, N, NQ
!     /ZVOD02/ --  HU
!
! Subroutines called by ZVINDY: DZSCAL, XERRWD
! Function routines called by ZVINDY: None
!-----------------------------------------------------------------------
! ZVINDY computes interpolated values of the K-th derivative of the
! dependent variable vector y, and stores it in DKY.  This routine
! is called within the package with K = 0 and T = TOUT, but may
! also be called by the user for any K up to the current order.
! (See detailed instructions in the usage documentation.)
!-----------------------------------------------------------------------
! The computed values in DKY are gotten by interpolation using the
! Nordsieck history array YH.  This array corresponds uniquely to a
! vector-valued polynomial of degree NQCUR or less, and DKY is set
! to the K-th derivative of this polynomial at T.
! The formula for DKY is:
!              q
!  DKY(i)  =  sum  c(j,K) * (T - TN)**(j-K) * H**(-j) * YH(i,j+1)
!             j=K
! where  c(j,K) = j*(j-1)*...*(j-K+1), q = NQCUR, TN = TCUR, H = HCUR.
! The quantities  NQ = NQCUR, L = NQ+1, N, TN, and H are
! communicated by COMMON.  The above sum is done in reverse order.
! IFLAG is returned negative if either K or T is out of bounds.
!
! Discussion above and comments in driver explain all variables.
!-----------------------------------------------------------------------
!
! Type declarations for labeled COMMON block ZVOD01 --------------------
!
DOUBLE PRECISION ACNRM, CCMXJ, CONP, CRATE, DRC, EL, &
& ETA, ETAMAX, H, HMIN, HMXI, HNEW, HRL1, HSCAL, PRL1, &
& RC, RL1, SRUR, TAU, TQ, TN, UROUND
INTEGER ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH, &
& L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, &
& LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP, &
& N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ, &
& NSLP, NYH
!
! Type declarations for labeled COMMON block ZVOD02 --------------------
!
DOUBLE PRECISION HU
INTEGER NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
!
! Type declarations for local variables --------------------------------
!
DOUBLE PRECISION C, HUN, R, S, TFUZZ, TN1, TP, ZERO
INTEGER I, IC, J, JB, JB2, JJ, JJ1, JP1
CHARACTER(LEN=80) MSG
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
SAVE HUN, ZERO
!
COMMON /ZVOD01/ ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13), ETA, &
& ETAMAX, H, HMIN, HMXI, HNEW, HRL1, HSCAL, PRL1, &
& RC, RL1, SRUR, TAU(13), TQ(5), TN, UROUND, &
& ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH, &
& L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, &
& LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP, &
& N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ, &
& NSLP, NYH
COMMON /ZVOD02/ HU, NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
!
DATA HUN /100.0D0/, ZERO /0.0D0/
!
IFLAG = 0
IF (K .LT. 0 .OR. K .GT. NQ) GO TO 80
TFUZZ = HUN*UROUND*SIGN(ABS(TN) + ABS(HU), HU)
TP = TN - HU - TFUZZ
TN1 = TN + TFUZZ
IF ((T-TP)*(T-TN1) .GT. ZERO) GO TO 90
!
S = (T - TN)/H
IC = 1
IF (K .EQ. 0) GO TO 15
JJ1 = L - K
DO 10 JJ = JJ1, NQ
  IC = IC*JJ
10 CONTINUE
15 C = REAL(IC)
DO 20 I = 1, N
  DKY(I) = C*YH(I,L)
20 CONTINUE
IF (K .EQ. NQ) GO TO 55
JB2 = NQ - K
DO 50 JB = 1, JB2
  J = NQ - JB
  JP1 = J + 1
  IC = 1
  IF (K .EQ. 0) GO TO 35
  JJ1 = JP1 - K
  DO 30 JJ = JJ1, J
    IC = IC*JJ
30 CONTINUE
35 C = REAL(IC)
  DO 40 I = 1, N
    DKY(I) = C*YH(I,JP1) + S*DKY(I)
40 CONTINUE
50 CONTINUE
IF (K .EQ. 0) RETURN
55 R = H**(-K)
CALL DZSCAL (N, R, DKY, 1)
RETURN
!
80 MSG = 'ZVINDY-- K (=I1) illegal      '
CALL XERRWD (MSG, 30, 51, 1, 1, K, 0, 0, ZERO, ZERO)
IFLAG = -1
RETURN
90 MSG = 'ZVINDY-- T (=R1) illegal      '
CALL XERRWD (MSG, 30, 52, 1, 0, 0, 0, 1, T, ZERO)
MSG='      T not in interval TCUR - HU (= R1) to TCUR (=R2)      '
CALL XERRWD (MSG, 60, 52, 1, 0, 0, 0, 2, TP, TN)
IFLAG = -2
RETURN
!----------------------- End of Subroutine ZVINDY ----------------------
END
!DECK ZVSTEP
SUBROUTINE ZVSTEP (Y, YH, LDYH, YH1, EWT, SAVF, VSAV, ACOR, &
& WM, IWM, F, JAC, PSOL, VNLS, RPAR, IPAR)
EXTERNAL F, JAC, PSOL, VNLS
COMPLEX(KIND=KIND(0.0d0)) Y, YH, YH1, SAVF, VSAV, ACOR, WM
DOUBLE PRECISION EWT
INTEGER LDYH, IWM, IPAR
DIMENSION Y(*), YH(LDYH,*), YH1(*), EWT(*), SAVF(*), VSAV(*), &
& ACOR(*), WM(*), IWM(*), RPAR(*), IPAR(*)
!-----------------------------------------------------------------------
! Call sequence input -- Y, YH, LDYH, YH1, EWT, SAVF, VSAV,
!                        ACOR, WM, IWM, F, JAC, PSOL, VNLS, RPAR, IPAR
! Call sequence output -- YH, ACOR, WM, IWM
! COMMON block variables accessed:
!     /ZVOD01/  ACNRM, EL(13), H, HMIN, HMXI, HNEW, HSCAL, RC, TAU(13),
!               TQ(5), TN, JCUR, JSTART, KFLAG, KUTH,
!               L, LMAX, MAXORD, N, NEWQ, NQ, NQWAIT
!     /ZVOD02/  HU, NCFN, NETF, NFE, NQU, NST
!
! Subroutines called by ZVSTEP: F, DZAXPY, ZCOPY, DZSCAL,
!                               ZVJUST, VNLS, ZVSET
! Function routines called by ZVSTEP: ZVNORM
!-----------------------------------------------------------------------
! ZVSTEP performs one step of the integration of an initial value
! problem for a system of ordinary differential equations.
! ZVSTEP calls subroutine VNLS for the solution of the nonlinear system
! arising in the time step.  Thus it is independent of the problem
! Jacobian structure and the type of nonlinear system solution method.
! ZVSTEP returns a completion flag KFLAG (in COMMON).
! A return with KFLAG = -1 or -2 means either ABS(H) = HMIN or 10
! consecutive failures occurred.  On a return with KFLAG negative,
! the values of TN and the YH array are as of the beginning of the last
! step, and H is the last step size attempted.
!
! Communication with ZVSTEP is done with the following variables:
!
! Y      = An array of length N used for the dependent variable vector.
! YH     = An LDYH by LMAX array containing the dependent variables
!          and their approximate scaled derivatives, where
!          LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate
!          j-th derivative of y(i), scaled by H**j/factorial(j)
!          (j = 0,1,...,NQ).  On entry for the first step, the first
!          two columns of YH must be set from the initial values.
! LDYH   = A constant integer .ge. N, the first dimension of YH.
!          N is the number of ODEs in the system.
! YH1    = A one-dimensional array occupying the same space as YH.
! EWT    = An array of length N containing multiplicative weights
!          for local error measurements.  Local errors in y(i) are
!          compared to 1.0/EWT(i) in various error tests.
! SAVF   = An array of working storage, of length N.
!          also used for input of YH(*,MAXORD+2) when JSTART = -1
!          and MAXORD .lt. the current order NQ.
! VSAV   = A work array of length N passed to subroutine VNLS.
! ACOR   = A work array of length N, used for the accumulated
!          corrections.  On a successful return, ACOR(i) contains
!          the estimated one-step local error in y(i).
! WM,IWM = Complex and integer work arrays associated with matrix
!          operations in VNLS.
! F      = Dummy name for the user-supplied subroutine for f.
! JAC    = Dummy name for the user-supplied Jacobian subroutine.
! PSOL   = Dummy name for the subroutine passed to VNLS, for
!          possible use there.
! VNLS   = Dummy name for the nonlinear system solving subroutine,
!          whose real name is dependent on the method used.
! RPAR, IPAR = User's real/complex and integer work arrays.
!-----------------------------------------------------------------------
!
! Type declarations for labeled COMMON block ZVOD01 --------------------
!
DOUBLE PRECISION ACNRM, CCMXJ, CONP, CRATE, DRC, EL, &
& ETA, ETAMAX, H, HMIN, HMXI, HNEW, HRL1, HSCAL, PRL1, &
& RC, RL1, SRUR, TAU, TQ, TN, UROUND
INTEGER ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH, &
& L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, &
& LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP, &
& N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ, &
& NSLP, NYH
!
! Type declarations for labeled COMMON block ZVOD02 --------------------
!
DOUBLE PRECISION HU
INTEGER NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
!
! Type declarations for local variables --------------------------------
!
DOUBLE PRECISION ADDON, BIAS1,BIAS2,BIAS3, CNQUOT, DDN, DSM, DUP, &
& ETACF, ETAMIN, ETAMX1, ETAMX2, ETAMX3, ETAMXF, &
& ETAQ, ETAQM1, ETAQP1, FLOTL, ONE, ONEPSM, &
& R, THRESH, TOLD, ZERO
INTEGER I, I1, I2, IBACK, J, JB, KFC, KFH, MXNCF, NCF, NFLAG
!
! Type declaration for function subroutines called ---------------------
!
DOUBLE PRECISION ZVNORM
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
SAVE ADDON, BIAS1, BIAS2, BIAS3, &
& ETACF, ETAMIN, ETAMX1, ETAMX2, ETAMX3, ETAMXF, ETAQ, ETAQM1, &
& KFC, KFH, MXNCF, ONEPSM, THRESH, ONE, ZERO
!-----------------------------------------------------------------------
COMMON /ZVOD01/ ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13), ETA, &
& ETAMAX, H, HMIN, HMXI, HNEW, HRL1, HSCAL, PRL1, &
& RC, RL1, SRUR, TAU(13), TQ(5), TN, UROUND, &
& ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH, &
& L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, &
& LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP, &
& N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ, &
& NSLP, NYH
COMMON /ZVOD02/ HU, NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
!
DATA KFC/-3/, KFH/-7/, MXNCF/10/
DATA ADDON  /1.0D-6/,    BIAS1  /6.0D0/,     BIAS2  /6.0D0/, &
& BIAS3  /10.0D0/,    ETACF  /0.25D0/,    ETAMIN /0.1D0/, &
& ETAMXF /0.2D0/,     ETAMX1 /1.0D4/,     ETAMX2 /10.0D0/, &
& ETAMX3 /10.0D0/,    ONEPSM /1.00001D0/, THRESH /1.5D0/
DATA ONE/1.0D0/, ZERO/0.0D0/
!
KFLAG = 0
TOLD = TN
NCF = 0
JCUR = 0
NFLAG = 0
IF (JSTART .GT. 0) GO TO 20
IF (JSTART .EQ. -1) GO TO 100
!-----------------------------------------------------------------------
! On the first call, the order is set to 1, and other variables are
! initialized.  ETAMAX is the maximum ratio by which H can be increased
! in a single step.  It is normally 10, but is larger during the
! first step to compensate for the small initial H.  If a failure
! occurs (in corrector convergence or error test), ETAMAX is set to 1
! for the next increase.
!-----------------------------------------------------------------------
LMAX = MAXORD + 1
NQ = 1
L = 2
NQNYH = NQ*LDYH
TAU(1) = H
PRL1 = ONE
RC = ZERO
ETAMAX = ETAMX1
NQWAIT = 2
HSCAL = H
GO TO 200
!-----------------------------------------------------------------------
! Take preliminary actions on a normal continuation step (JSTART.GT.0).
! If the driver changed H, then ETA must be reset and NEWH set to 1.
! If a change of order was dictated on the previous step, then
! it is done here and appropriate adjustments in the history are made.
! On an order decrease, the history array is adjusted by ZVJUST.
! On an order increase, the history array is augmented by a column.
! On a change of step size H, the history array YH is rescaled.
!-----------------------------------------------------------------------
20 CONTINUE
IF (KUTH .EQ. 1) THEN
  ETA = MIN(ETA,H/HSCAL)
  NEWH = 1
  ENDIF
50 IF (NEWH .EQ. 0) GO TO 200
IF (NEWQ .EQ. NQ) GO TO 150
IF (NEWQ .LT. NQ) THEN
  CALL ZVJUST (YH, LDYH, -1)
  NQ = NEWQ
  L = NQ + 1
  NQWAIT = L
  GO TO 150
  ENDIF
IF (NEWQ .GT. NQ) THEN
  CALL ZVJUST (YH, LDYH, 1)
  NQ = NEWQ
  L = NQ + 1
  NQWAIT = L
  GO TO 150
ENDIF
!-----------------------------------------------------------------------
! The following block handles preliminaries needed when JSTART = -1.
! If N was reduced, zero out part of YH to avoid undefined references.
! If MAXORD was reduced to a value less than the tentative order NEWQ,
! then NQ is set to MAXORD, and a new H ratio ETA is chosen.
! Otherwise, we take the same preliminary actions as for JSTART .gt. 0.
! In any case, NQWAIT is reset to L = NQ + 1 to prevent further
! changes in order for that many steps.
! The new H ratio ETA is limited by the input H if KUTH = 1,
! by HMIN if KUTH = 0, and by HMXI in any case.
! Finally, the history array YH is rescaled.
!-----------------------------------------------------------------------
100 CONTINUE
LMAX = MAXORD + 1
IF (N .EQ. LDYH) GO TO 120
I1 = 1 + (NEWQ + 1)*LDYH
I2 = (MAXORD + 1)*LDYH
IF (I1 .GT. I2) GO TO 120
DO 110 I = I1, I2
  YH1(I) = ZERO
110 CONTINUE
120 IF (NEWQ .LE. MAXORD) GO TO 140
FLOTL = REAL(LMAX)
IF (MAXORD .LT. NQ-1) THEN
  DDN = ZVNORM (N, SAVF, EWT)/TQ(1)
  ETA = ONE/((BIAS1*DDN)**(ONE/FLOTL) + ADDON)
  ENDIF
IF (MAXORD .EQ. NQ .AND. NEWQ .EQ. NQ+1) ETA = ETAQ
IF (MAXORD .EQ. NQ-1 .AND. NEWQ .EQ. NQ+1) THEN
  ETA = ETAQM1
  CALL ZVJUST (YH, LDYH, -1)
  ENDIF
IF (MAXORD .EQ. NQ-1 .AND. NEWQ .EQ. NQ) THEN
  DDN = ZVNORM (N, SAVF, EWT)/TQ(1)
  ETA = ONE/((BIAS1*DDN)**(ONE/FLOTL) + ADDON)
  CALL ZVJUST (YH, LDYH, -1)
  ENDIF
ETA = MIN(ETA,ONE)
NQ = MAXORD
L = LMAX
140 IF (KUTH .EQ. 1) ETA = MIN(ETA,ABS(H/HSCAL))
IF (KUTH .EQ. 0) ETA = MAX(ETA,HMIN/ABS(HSCAL))
ETA = ETA/MAX(ONE,ABS(HSCAL)*HMXI*ETA)
NEWH = 1
NQWAIT = L
IF (NEWQ .LE. MAXORD) GO TO 50
! Rescale the history array for a change in H by a factor of ETA. ------
150 R = ONE
DO 180 J = 2, L
  R = R*ETA
  CALL DZSCAL (N, R, YH(1,J), 1 )
180 CONTINUE
H = HSCAL*ETA
HSCAL = H
RC = RC*ETA
NQNYH = NQ*LDYH
!-----------------------------------------------------------------------
! This section computes the predicted values by effectively
! multiplying the YH array by the Pascal triangle matrix.
! ZVSET is called to calculate all integration coefficients.
! RC is the ratio of new to old values of the coefficient H/EL(2)=h/l1.
!-----------------------------------------------------------------------
200 TN = TN + H
I1 = NQNYH + 1
DO 220 JB = 1, NQ
  I1 = I1 - LDYH
  DO 210 I = I1, NQNYH
    YH1(I) = YH1(I) + YH1(I+LDYH)
210 CONTINUE
220 CONTINUE
CALL ZVSET
RL1 = ONE/EL(2)
RC = RC*(RL1/PRL1)
PRL1 = RL1
!
! Call the nonlinear system solver. ------------------------------------
!
CALL VNLS (Y, YH, LDYH, VSAV, SAVF, EWT, ACOR, IWM, WM, &
& F, JAC, PSOL, NFLAG, RPAR, IPAR)
!
IF (NFLAG .EQ. 0) GO TO 450
!-----------------------------------------------------------------------
! The VNLS routine failed to achieve convergence (NFLAG .NE. 0).
! The YH array is retracted to its values before prediction.
! The step size H is reduced and the step is retried, if possible.
! Otherwise, an error exit is taken.
!-----------------------------------------------------------------------
  NCF = NCF + 1
  NCFN = NCFN + 1
  ETAMAX = ONE
  TN = TOLD
  I1 = NQNYH + 1
  DO 430 JB = 1, NQ
    I1 = I1 - LDYH
    DO 420 I = I1, NQNYH
      YH1(I) = YH1(I) - YH1(I+LDYH)
420 CONTINUE
430 CONTINUE
  IF (NFLAG .LT. -1) GO TO 680
  IF (ABS(H) .LE. HMIN*ONEPSM) GO TO 670
  IF (NCF .EQ. MXNCF) GO TO 670
  ETA = ETACF
  ETA = MAX(ETA,HMIN/ABS(H))
  NFLAG = -1
  GO TO 150
!-----------------------------------------------------------------------
! The corrector has converged (NFLAG = 0).  The local error test is
! made and control passes to statement 500 if it fails.
!-----------------------------------------------------------------------
450 CONTINUE
DSM = ACNRM/TQ(2)
IF (DSM .GT. ONE) GO TO 500
!-----------------------------------------------------------------------
! After a successful step, update the YH and TAU arrays and decrement
! NQWAIT.  If NQWAIT is then 1 and NQ .lt. MAXORD, then ACOR is saved
! for use in a possible order increase on the next step.
! If ETAMAX = 1 (a failure occurred this step), keep NQWAIT .ge. 2.
!-----------------------------------------------------------------------
KFLAG = 0
NST = NST + 1
HU = H
NQU = NQ
DO 470 IBACK = 1, NQ
  I = L - IBACK
  TAU(I+1) = TAU(I)
470 CONTINUE
TAU(1) = H
DO 480 J = 1, L
  CALL DZAXPY (N, EL(J), ACOR, 1, YH(1,J), 1 )
480 CONTINUE
NQWAIT = NQWAIT - 1
IF ((L .EQ. LMAX) .OR. (NQWAIT .NE. 1)) GO TO 490
CALL ZCOPY (N, ACOR, 1, YH(1,LMAX), 1 )
CONP = TQ(5)
490 IF (ETAMAX .NE. ONE) GO TO 560
IF (NQWAIT .LT. 2) NQWAIT = 2
NEWQ = NQ
NEWH = 0
ETA = ONE
HNEW = H
GO TO 690
!-----------------------------------------------------------------------
! The error test failed.  KFLAG keeps track of multiple failures.
! Restore TN and the YH array to their previous values, and prepare
! to try the step again.  Compute the optimum step size for the
! same order.  After repeated failures, H is forced to decrease
! more rapidly.
!-----------------------------------------------------------------------
500 KFLAG = KFLAG - 1
NETF = NETF + 1
NFLAG = -2
TN = TOLD
I1 = NQNYH + 1
DO 520 JB = 1, NQ
  I1 = I1 - LDYH
  DO 510 I = I1, NQNYH
    YH1(I) = YH1(I) - YH1(I+LDYH)
510 CONTINUE
520 CONTINUE
IF (ABS(H) .LE. HMIN*ONEPSM) GO TO 660
ETAMAX = ONE
IF (KFLAG .LE. KFC) GO TO 530
! Compute ratio of new H to current H at the current order. ------------
FLOTL = REAL(L)
ETA = ONE/((BIAS2*DSM)**(ONE/FLOTL) + ADDON)
ETA = MAX(ETA,HMIN/ABS(H),ETAMIN)
IF ((KFLAG .LE. -2) .AND. (ETA .GT. ETAMXF)) ETA = ETAMXF
GO TO 150
!-----------------------------------------------------------------------
! Control reaches this section if 3 or more consecutive failures
! have occurred.  It is assumed that the elements of the YH array
! have accumulated errors of the wrong order.  The order is reduced
! by one, if possible.  Then H is reduced by a factor of 0.1 and
! the step is retried.  After a total of 7 consecutive failures,
! an exit is taken with KFLAG = -1.
!-----------------------------------------------------------------------
530 IF (KFLAG .EQ. KFH) GO TO 660
IF (NQ .EQ. 1) GO TO 540
ETA = MAX(ETAMIN,HMIN/ABS(H))
CALL ZVJUST (YH, LDYH, -1)
L = NQ
NQ = NQ - 1
NQWAIT = L
GO TO 150
540 ETA = MAX(ETAMIN,HMIN/ABS(H))
H = H*ETA
HSCAL = H
TAU(1) = H
CALL F (N, TN, Y, SAVF, RPAR, IPAR)
NFE = NFE + 1
DO 550 I = 1, N
  YH(I,2) = H*SAVF(I)
550 CONTINUE
NQWAIT = 10
GO TO 200
!-----------------------------------------------------------------------
! If NQWAIT = 0, an increase or decrease in order by one is considered.
! Factors ETAQ, ETAQM1, ETAQP1 are computed by which H could
! be multiplied at order q, q-1, or q+1, respectively.
! The largest of these is determined, and the new order and
! step size set accordingly.
! A change of H or NQ is made only if H increases by at least a
! factor of THRESH.  If an order change is considered and rejected,
! then NQWAIT is set to 2 (reconsider it after 2 steps).
!-----------------------------------------------------------------------
! Compute ratio of new H to current H at the current order. ------------
560 FLOTL = REAL(L)
ETAQ = ONE/((BIAS2*DSM)**(ONE/FLOTL) + ADDON)
IF (NQWAIT .NE. 0) GO TO 600
NQWAIT = 2
ETAQM1 = ZERO
IF (NQ .EQ. 1) GO TO 570
! Compute ratio of new H to current H at the current order less one. ---
DDN = ZVNORM (N, YH(1,L), EWT)/TQ(1)
ETAQM1 = ONE/((BIAS1*DDN)**(ONE/(FLOTL - ONE)) + ADDON)
570 ETAQP1 = ZERO
IF (L .EQ. LMAX) GO TO 580
! Compute ratio of new H to current H at current order plus one. -------
CNQUOT = (TQ(5)/CONP)*(H/TAU(2))**L
DO 575 I = 1, N
  SAVF(I) = ACOR(I) - CNQUOT*YH(I,LMAX)
575 CONTINUE
DUP = ZVNORM (N, SAVF, EWT)/TQ(3)
ETAQP1 = ONE/((BIAS3*DUP)**(ONE/(FLOTL + ONE)) + ADDON)
580 IF (ETAQ .GE. ETAQP1) GO TO 590
IF (ETAQP1 .GT. ETAQM1) GO TO 620
GO TO 610
590 IF (ETAQ .LT. ETAQM1) GO TO 610
600 ETA = ETAQ
NEWQ = NQ
GO TO 630
610 ETA = ETAQM1
NEWQ = NQ - 1
GO TO 630
620 ETA = ETAQP1
NEWQ = NQ + 1
CALL ZCOPY (N, ACOR, 1, YH(1,LMAX), 1)
! Test tentative new H against THRESH, ETAMAX, and HMXI, then exit. ----
630 IF (ETA .LT. THRESH .OR. ETAMAX .EQ. ONE) GO TO 640
ETA = MIN(ETA,ETAMAX)
ETA = ETA/MAX(ONE,ABS(H)*HMXI*ETA)
NEWH = 1
HNEW = H*ETA
GO TO 690
640 NEWQ = NQ
NEWH = 0
ETA = ONE
HNEW = H
GO TO 690
!-----------------------------------------------------------------------
! All returns are made through this section.
! On a successful return, ETAMAX is reset and ACOR is scaled.
!-----------------------------------------------------------------------
660 KFLAG = -1
GO TO 720
670 KFLAG = -2
GO TO 720
680 IF (NFLAG .EQ. -2) KFLAG = -3
IF (NFLAG .EQ. -3) KFLAG = -4
GO TO 720
690 ETAMAX = ETAMX3
IF (NST .LE. 10) ETAMAX = ETAMX2
R = ONE/TQ(2)
CALL DZSCAL (N, R, ACOR, 1)
720 JSTART = 1
RETURN
!----------------------- End of Subroutine ZVSTEP ----------------------
END
!DECK ZVSET
SUBROUTINE ZVSET
!-----------------------------------------------------------------------
! Call sequence communication: None
! COMMON block variables accessed:
!     /ZVOD01/ -- EL(13), H, TAU(13), TQ(5), L(= NQ + 1),
!                 METH, NQ, NQWAIT
!
! Subroutines called by ZVSET: None
! Function routines called by ZVSET: None
!-----------------------------------------------------------------------
! ZVSET is called by ZVSTEP and sets coefficients for use there.
!
! For each order NQ, the coefficients in EL are calculated by use of
!  the generating polynomial lambda(x), with coefficients EL(i).
!      lambda(x) = EL(1) + EL(2)*x + ... + EL(NQ+1)*(x**NQ).
! For the backward differentiation formulas,
!                                     NQ-1
!      lambda(x) = (1 + x/xi*(NQ)) * product (1 + x/xi(i) ) .
!                                     i = 1
! For the Adams formulas,
!                              NQ-1
!      (d/dx) lambda(x) = c * product (1 + x/xi(i) ) ,
!                              i = 1
!      lambda(-1) = 0,    lambda(0) = 1,
! where c is a normalization constant.
! In both cases, xi(i) is defined by
!      H*xi(i) = t sub n  -  t sub (n-i)
!              = H + TAU(1) + TAU(2) + ... TAU(i-1).
!
!
! In addition to variables described previously, communication
! with ZVSET uses the following:
!   TAU    = A vector of length 13 containing the past NQ values
!            of H.
!   EL     = A vector of length 13 in which vset stores the
!            coefficients for the corrector formula.
!   TQ     = A vector of length 5 in which vset stores constants
!            used for the convergence test, the error test, and the
!            selection of H at a new order.
!   METH   = The basic method indicator.
!   NQ     = The current order.
!   L      = NQ + 1, the length of the vector stored in EL, and
!            the number of columns of the YH array being used.
!   NQWAIT = A counter controlling the frequency of order changes.
!            An order change is about to be considered if NQWAIT = 1.
!-----------------------------------------------------------------------
!
! Type declarations for labeled COMMON block ZVOD01 --------------------
!
DOUBLE PRECISION ACNRM, CCMXJ, CONP, CRATE, DRC, EL, &
& ETA, ETAMAX, H, HMIN, HMXI, HNEW, HRL1, HSCAL, PRL1, &
& RC, RL1, SRUR, TAU, TQ, TN, UROUND
INTEGER ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH, &
& L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, &
& LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP, &
& N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ, &
& NSLP, NYH
!
! Type declarations for local variables --------------------------------
!
DOUBLE PRECISION AHATN0, ALPH0, CNQM1, CORTES, CSUM, ELP, EM, &
& EM0, FLOTI, FLOTL, FLOTNQ, HSUM, ONE, RXI, RXIS, S, SIX, &
& T1, T2, T3, T4, T5, T6, TWO, XI, ZERO
INTEGER I, IBACK, J, JP1, NQM1, NQM2
!
DIMENSION EM(13)
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
SAVE CORTES, ONE, SIX, TWO, ZERO
!
COMMON /ZVOD01/ ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13), ETA, &
& ETAMAX, H, HMIN, HMXI, HNEW, HRL1, HSCAL, PRL1, &
& RC, RL1, SRUR, TAU(13), TQ(5), TN, UROUND, &
& ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH, &
& L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, &
& LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP, &
& N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ, &
& NSLP, NYH
!
DATA CORTES /0.1D0/
DATA ONE  /1.0D0/, SIX /6.0D0/, TWO /2.0D0/, ZERO /0.0D0/
!
FLOTL = REAL(L)
NQM1 = NQ - 1
NQM2 = NQ - 2
IF (METH .EQ. 1) THEN
  GOTO 100
ELSE IF (METH .EQ. 2) THEN
  GOTO 200
ENDIF
!      GO TO (100, 200), METH
!
! Set coefficients for Adams methods. ----------------------------------
100 IF (NQ .NE. 1) GO TO 110
EL(1) = ONE
EL(2) = ONE
TQ(1) = ONE
TQ(2) = TWO
TQ(3) = SIX*TQ(2)
TQ(5) = ONE
GO TO 300
110 HSUM = H
EM(1) = ONE
FLOTNQ = FLOTL - ONE
DO 115 I = 2, L
  EM(I) = ZERO
115 CONTINUE
DO 150 J = 1, NQM1
  IF ((J .NE. NQM1) .OR. (NQWAIT .NE. 1)) GO TO 130
  S = ONE
  CSUM = ZERO
  DO 120 I = 1, NQM1
    CSUM = CSUM + S*EM(I)/REAL(I+1)
    S = -S
120 CONTINUE
  TQ(1) = EM(NQM1)/(FLOTNQ*CSUM)
130 RXI = H/HSUM
  DO 140 IBACK = 1, J
    I = (J + 2) - IBACK
    EM(I) = EM(I) + EM(I-1)*RXI
140 CONTINUE
  HSUM = HSUM + TAU(J)
150 CONTINUE
! Compute integral from -1 to 0 of polynomial and of x times it. -------
S = ONE
EM0 = ZERO
CSUM = ZERO
DO 160 I = 1, NQ
  FLOTI = REAL(I)
  EM0 = EM0 + S*EM(I)/FLOTI
  CSUM = CSUM + S*EM(I)/(FLOTI+ONE)
  S = -S
160 CONTINUE
! In EL, form coefficients of normalized integrated polynomial. --------
S = ONE/EM0
EL(1) = ONE
DO 170 I = 1, NQ
  EL(I+1) = S*EM(I)/REAL(I)
170 CONTINUE
XI = HSUM/H
TQ(2) = XI*EM0/CSUM
TQ(5) = XI/EL(L)
IF (NQWAIT .NE. 1) GO TO 300
! For higher order control constant, multiply polynomial by 1+x/xi(q). -
RXI = ONE/XI
DO 180 IBACK = 1, NQ
  I = (L + 1) - IBACK
  EM(I) = EM(I) + EM(I-1)*RXI
180 CONTINUE
! Compute integral of polynomial. --------------------------------------
S = ONE
CSUM = ZERO
DO 190 I = 1, L
  CSUM = CSUM + S*EM(I)/REAL(I+1)
  S = -S
190 CONTINUE
TQ(3) = FLOTL*EM0/CSUM
GO TO 300
!
! Set coefficients for BDF methods. ------------------------------------
200 DO 210 I = 3, L
  EL(I) = ZERO
210 CONTINUE
EL(1) = ONE
EL(2) = ONE
ALPH0 = -ONE
AHATN0 = -ONE
HSUM = H
RXI = ONE
RXIS = ONE
IF (NQ .EQ. 1) GO TO 240
DO 230 J = 1, NQM2
! In EL, construct coefficients of (1+x/xi(1))*...*(1+x/xi(j+1)). ------
  HSUM = HSUM + TAU(J)
  RXI = H/HSUM
  JP1 = J + 1
  ALPH0 = ALPH0 - ONE/REAL(JP1)
  DO 220 IBACK = 1, JP1
    I = (J + 3) - IBACK
    EL(I) = EL(I) + EL(I-1)*RXI
220 CONTINUE
230 CONTINUE
ALPH0 = ALPH0 - ONE/REAL(NQ)
RXIS = -EL(2) - ALPH0
HSUM = HSUM + TAU(NQM1)
RXI = H/HSUM
AHATN0 = -EL(2) - RXI
DO 235 IBACK = 1, NQ
  I = (NQ + 2) - IBACK
  EL(I) = EL(I) + EL(I-1)*RXIS
235 CONTINUE
240 T1 = ONE - AHATN0 + ALPH0
T2 = ONE + REAL(NQ)*T1
TQ(2) = ABS(ALPH0*T2/T1)
TQ(5) = ABS(T2/(EL(L)*RXI/RXIS))
IF (NQWAIT .NE. 1) GO TO 300
CNQM1 = RXIS/EL(L)
T3 = ALPH0 + ONE/REAL(NQ)
T4 = AHATN0 + RXI
ELP = T3/(ONE - T4 + T3)
TQ(1) = ABS(ELP/CNQM1)
HSUM = HSUM + TAU(NQ)
RXI = H/HSUM
T5 = ALPH0 - ONE/REAL(NQ+1)
T6 = AHATN0 - RXI
ELP = T2/(ONE - T6 + T5)
TQ(3) = ABS(ELP*RXI*(FLOTL + ONE)*T5)
300 TQ(4) = CORTES*TQ(2)
RETURN
!----------------------- End of Subroutine ZVSET -----------------------
END
!DECK ZVJUST
SUBROUTINE ZVJUST (YH, LDYH, IORD)
COMPLEX(KIND=KIND(0.0d0)) YH
INTEGER LDYH, IORD
DIMENSION YH(LDYH,*)
!-----------------------------------------------------------------------
! Call sequence input -- YH, LDYH, IORD
! Call sequence output -- YH
! COMMON block input -- NQ, METH, LMAX, HSCAL, TAU(13), N
! COMMON block variables accessed:
!     /ZVOD01/ -- HSCAL, TAU(13), LMAX, METH, N, NQ,
!
! Subroutines called by ZVJUST: DZAXPY
! Function routines called by ZVJUST: None
!-----------------------------------------------------------------------
! This subroutine adjusts the YH array on reduction of order,
! and also when the order is increased for the stiff option (METH = 2).
! Communication with ZVJUST uses the following:
! IORD  = An integer flag used when METH = 2 to indicate an order
!         increase (IORD = +1) or an order decrease (IORD = -1).
! HSCAL = Step size H used in scaling of Nordsieck array YH.
!         (If IORD = +1, ZVJUST assumes that HSCAL = TAU(1).)
! See References 1 and 2 for details.
!-----------------------------------------------------------------------
!
! Type declarations for labeled COMMON block ZVOD01 --------------------
!
DOUBLE PRECISION ACNRM, CCMXJ, CONP, CRATE, DRC, EL, &
& ETA, ETAMAX, H, HMIN, HMXI, HNEW, HRL1, HSCAL, PRL1, &
& RC, RL1, SRUR, TAU, TQ, TN, UROUND
INTEGER ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH, &
& L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, &
& LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP, &
& N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ, &
& NSLP, NYH
!
! Type declarations for local variables --------------------------------
!
DOUBLE PRECISION ALPH0, ALPH1, HSUM, ONE, PROD, T1, XI,XIOLD, ZERO
INTEGER I, IBACK, J, JP1, LP1, NQM1, NQM2, NQP1
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
SAVE ONE, ZERO
!
COMMON /ZVOD01/ ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13), ETA, &
& ETAMAX, H, HMIN, HMXI, HNEW, HRL1, HSCAL, PRL1, &
& RC, RL1, SRUR, TAU(13), TQ(5), TN, UROUND, &
& ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH, &
& L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, &
& LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP, &
& N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ, &
& NSLP, NYH
!
DATA ONE /1.0D0/, ZERO /0.0D0/
!
IF ((NQ .EQ. 2) .AND. (IORD .NE. 1)) RETURN
NQM1 = NQ - 1
NQM2 = NQ - 2
IF (METH .EQ. 1) THEN
  GOTO 100
ELSE IF (METH .EQ. 2) THEN
  GOTO 200
ENDIF
!      GO TO (100, 200), METH
!-----------------------------------------------------------------------
! Nonstiff option...
! Check to see if the order is being increased or decreased.
!-----------------------------------------------------------------------
100 CONTINUE
IF (IORD .EQ. 1) GO TO 180
! Order decrease. ------------------------------------------------------
DO 110 J = 1, LMAX
  EL(J) = ZERO
110 CONTINUE
EL(2) = ONE
HSUM = ZERO
DO 130 J = 1, NQM2
! Construct coefficients of x*(x+xi(1))*...*(x+xi(j)). -----------------
  HSUM = HSUM + TAU(J)
  XI = HSUM/HSCAL
  JP1 = J + 1
  DO 120 IBACK = 1, JP1
    I = (J + 3) - IBACK
    EL(I) = EL(I)*XI + EL(I-1)
120 CONTINUE
130 CONTINUE
! Construct coefficients of integrated polynomial. ---------------------
DO 140 J = 2, NQM1
  EL(J+1) = REAL(NQ)*EL(J)/REAL(J)
140 CONTINUE
! Subtract correction terms from YH array. -----------------------------
DO 170 J = 3, NQ
  DO 160 I = 1, N
   YH(I,J) = YH(I,J) - YH(I,L)*EL(J)
160 CONTINUE
170 CONTINUE
RETURN
! Order increase. ------------------------------------------------------
! Zero out next column in YH array. ------------------------------------
180 CONTINUE
LP1 = L + 1
DO 190 I = 1, N
  YH(I,LP1) = ZERO
190 CONTINUE
RETURN
!-----------------------------------------------------------------------
! Stiff option...
! Check to see if the order is being increased or decreased.
!-----------------------------------------------------------------------
200 CONTINUE
IF (IORD .EQ. 1) GO TO 300
! Order decrease. ------------------------------------------------------
DO 210 J = 1, LMAX
  EL(J) = ZERO
210 CONTINUE
EL(3) = ONE
HSUM = ZERO
DO 230 J = 1,NQM2
! Construct coefficients of x*x*(x+xi(1))*...*(x+xi(j)). ---------------
  HSUM = HSUM + TAU(J)
  XI = HSUM/HSCAL
  JP1 = J + 1
  DO 220 IBACK = 1, JP1
    I = (J + 4) - IBACK
    EL(I) = EL(I)*XI + EL(I-1)
220 CONTINUE
230 CONTINUE
! Subtract correction terms from YH array. -----------------------------
DO 250 J = 3,NQ
  DO 240 I = 1, N
    YH(I,J) = YH(I,J) - YH(I,L)*EL(J)
240 CONTINUE
250 CONTINUE
RETURN
! Order increase. ------------------------------------------------------
300 DO 310 J = 1, LMAX
  EL(J) = ZERO
310 CONTINUE
EL(3) = ONE
ALPH0 = -ONE
ALPH1 = ONE
PROD = ONE
XIOLD = ONE
HSUM = HSCAL
IF (NQ .EQ. 1) GO TO 340
DO 330 J = 1, NQM1
! Construct coefficients of x*x*(x+xi(1))*...*(x+xi(j)). ---------------
  JP1 = J + 1
  HSUM = HSUM + TAU(JP1)
  XI = HSUM/HSCAL
  PROD = PROD*XI
  ALPH0 = ALPH0 - ONE/REAL(JP1)
  ALPH1 = ALPH1 + ONE/XI
  DO 320 IBACK = 1, JP1
    I = (J + 4) - IBACK
    EL(I) = EL(I)*XIOLD + EL(I-1)
320 CONTINUE
  XIOLD = XI
330 CONTINUE
340 CONTINUE
T1 = (-ALPH0 - ALPH1)/PROD
! Load column L + 1 in YH array. ---------------------------------------
LP1 = L + 1
DO 350 I = 1, N
  YH(I,LP1) = T1*YH(I,LMAX)
350 CONTINUE
! Add correction terms to YH array. ------------------------------------
NQP1 = NQ + 1
DO 370 J = 3, NQP1
  CALL DZAXPY (N, EL(J), YH(1,LP1), 1, YH(1,J), 1 )
370 CONTINUE
RETURN
!----------------------- End of Subroutine ZVJUST ----------------------
END
!DECK ZVNLSD
SUBROUTINE ZVNLSD (Y, YH, LDYH, VSAV, SAVF, EWT, ACOR, IWM, WM, &
& F, JAC, PDUM, NFLAG, RPAR, IPAR)
EXTERNAL F, JAC, PDUM
COMPLEX(KIND=KIND(0.0d0)) Y, YH, VSAV, SAVF, ACOR, WM
DOUBLE PRECISION EWT
INTEGER LDYH, IWM, NFLAG, IPAR
DIMENSION Y(*), YH(LDYH,*), VSAV(*), SAVF(*), EWT(*), ACOR(*), &
& IWM(*), WM(*), RPAR(*), IPAR(*)
!-----------------------------------------------------------------------
! Call sequence input -- Y, YH, LDYH, SAVF, EWT, ACOR, IWM, WM,
!                        F, JAC, NFLAG, RPAR, IPAR
! Call sequence output -- YH, ACOR, WM, IWM, NFLAG
! COMMON block variables accessed:
!     /ZVOD01/ ACNRM, CRATE, DRC, H, RC, RL1, TQ(5), TN, ICF,
!                JCUR, METH, MITER, N, NSLP
!     /ZVOD02/ HU, NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
!
! Subroutines called by ZVNLSD: F, DZAXPY, ZCOPY, DZSCAL, ZVJAC, ZVSOL
! Function routines called by ZVNLSD: ZVNORM
!-----------------------------------------------------------------------
! Subroutine ZVNLSD is a nonlinear system solver, which uses functional
! iteration or a chord (modified Newton) method.  For the chord method
! direct linear algebraic system solvers are used.  Subroutine ZVNLSD
! then handles the corrector phase of this integration package.
!
! Communication with ZVNLSD is done with the following variables. (For
! more details, please see the comments in the driver subroutine.)
!
! Y          = The dependent variable, a vector of length N, input.
! YH         = The Nordsieck (Taylor) array, LDYH by LMAX, input
!              and output.  On input, it contains predicted values.
! LDYH       = A constant .ge. N, the first dimension of YH, input.
! VSAV       = Unused work array.
! SAVF       = A work array of length N.
! EWT        = An error weight vector of length N, input.
! ACOR       = A work array of length N, used for the accumulated
!              corrections to the predicted y vector.
! WM,IWM     = Complex and integer work arrays associated with matrix
!              operations in chord iteration (MITER .ne. 0).
! F          = Dummy name for user-supplied routine for f.
! JAC        = Dummy name for user-supplied Jacobian routine.
! PDUM       = Unused dummy subroutine name.  Included for uniformity
!              over collection of integrators.
! NFLAG      = Input/output flag, with values and meanings as follows:
!              INPUT
!                  0 first call for this time step.
!                 -1 convergence failure in previous call to ZVNLSD.
!                 -2 error test failure in ZVSTEP.
!              OUTPUT
!                  0 successful completion of nonlinear solver.
!                 -1 convergence failure or singular matrix.
!                 -2 unrecoverable error in matrix preprocessing
!                    (cannot occur here).
!                 -3 unrecoverable error in solution (cannot occur
!                    here).
! RPAR, IPAR = User's real/complex and integer work arrays.
!
! IPUP       = Own variable flag with values and meanings as follows:
!              0,            do not update the Newton matrix.
!              MITER .ne. 0, update Newton matrix, because it is the
!                            initial step, order was changed, the error
!                            test failed, or an update is indicated by
!                            the scalar RC or step counter NST.
!
! For more details, see comments in driver subroutine.
!-----------------------------------------------------------------------
! Type declarations for labeled COMMON block ZVOD01 --------------------
!
DOUBLE PRECISION ACNRM, CCMXJ, CONP, CRATE, DRC, EL, &
& ETA, ETAMAX, H, HMIN, HMXI, HNEW, HRL1, HSCAL, PRL1, &
& RC, RL1, SRUR, TAU, TQ, TN, UROUND
INTEGER ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH, &
& L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, &
& LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP, &
& N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ, &
& NSLP, NYH
!
! Type declarations for labeled COMMON block ZVOD02 --------------------
!
DOUBLE PRECISION HU
INTEGER NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
!
! Type declarations for local variables --------------------------------
!
DOUBLE PRECISION CCMAX, CRDOWN, CSCALE, DCON, DEL, DELP, ONE, &
& RDIV, TWO, ZERO
INTEGER I, IERPJ, IERSL, M, MAXCOR, MSBP
!
! Type declaration for function subroutines called ---------------------
!
DOUBLE PRECISION ZVNORM
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
SAVE CCMAX, CRDOWN, MAXCOR, MSBP, RDIV, ONE, TWO, ZERO
!
COMMON /ZVOD01/ ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13), ETA, &
& ETAMAX, H, HMIN, HMXI, HNEW, HRL1, HSCAL, PRL1, &
& RC, RL1, SRUR, TAU(13), TQ(5), TN, UROUND, &
& ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH, &
& L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, &
& LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP, &
& N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ, &
& NSLP, NYH
COMMON /ZVOD02/ HU, NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
!
DATA CCMAX /0.3D0/, CRDOWN /0.3D0/, MAXCOR /3/, MSBP /20/, &
& RDIV  /2.0D0/
DATA ONE /1.0D0/, TWO /2.0D0/, ZERO /0.0D0/
!-----------------------------------------------------------------------
! On the first step, on a change of method order, or after a
! nonlinear convergence failure with NFLAG = -2, set IPUP = MITER
! to force a Jacobian update when MITER .ne. 0.
!-----------------------------------------------------------------------
IF (JSTART .EQ. 0) NSLP = 0
IF (NFLAG .EQ. 0) ICF = 0
IF (NFLAG .EQ. -2) IPUP = MITER
IF ( (JSTART .EQ. 0) .OR. (JSTART .EQ. -1) ) IPUP = MITER
! If this is functional iteration, set CRATE .eq. 1 and drop to 220
IF (MITER .EQ. 0) THEN
  CRATE = ONE
  GO TO 220
ENDIF
!-----------------------------------------------------------------------
! RC is the ratio of new to old values of the coefficient H/EL(2)=h/l1.
! When RC differs from 1 by more than CCMAX, IPUP is set to MITER
! to force ZVJAC to be called, if a Jacobian is involved.
! In any case, ZVJAC is called at least every MSBP steps.
!-----------------------------------------------------------------------
DRC = ABS(RC-ONE)
IF (DRC .GT. CCMAX .OR. NST .GE. NSLP+MSBP) IPUP = MITER
!-----------------------------------------------------------------------
! Up to MAXCOR corrector iterations are taken.  A convergence test is
! made on the r.m.s. norm of each correction, weighted by the error
! weight vector EWT.  The sum of the corrections is accumulated in the
! vector ACOR(i).  The YH array is not altered in the corrector loop.
!-----------------------------------------------------------------------
220 M = 0
DELP = ZERO
CALL ZCOPY (N, YH(1,1), 1, Y, 1 )
CALL F (N, TN, Y, SAVF, RPAR, IPAR)
NFE = NFE + 1
IF (IPUP .LE. 0) GO TO 250
!-----------------------------------------------------------------------
! If indicated, the matrix P = I - h*rl1*J is reevaluated and
! preprocessed before starting the corrector iteration.  IPUP is set
! to 0 as an indicator that this has been done.
!-----------------------------------------------------------------------
CALL ZVJAC (Y, YH, LDYH, EWT, ACOR, SAVF, WM, IWM, F, JAC, IERPJ, &
& RPAR, IPAR)
IPUP = 0
RC = ONE
DRC = ZERO
CRATE = ONE
NSLP = NST
! If matrix is singular, take error return to force cut in step size. --
IF (IERPJ .NE. 0) GO TO 430
250 DO 260 I = 1,N
  ACOR(I) = ZERO
260 CONTINUE
! This is a looping point for the corrector iteration. -----------------
270 IF (MITER .NE. 0) GO TO 350
!-----------------------------------------------------------------------
! In the case of functional iteration, update Y directly from
! the result of the last function evaluation.
!-----------------------------------------------------------------------
DO 280 I = 1,N
  SAVF(I) = RL1*(H*SAVF(I) - YH(I,2))
280 CONTINUE
DO 290 I = 1,N
  Y(I) = SAVF(I) - ACOR(I)
290 CONTINUE
DEL = ZVNORM (N, Y, EWT)
DO 300 I = 1,N
  Y(I) = YH(I,1) + SAVF(I)
300 CONTINUE
CALL ZCOPY (N, SAVF, 1, ACOR, 1)
GO TO 400
!-----------------------------------------------------------------------
! In the case of the chord method, compute the corrector error,
! and solve the linear system with that as right-hand side and
! P as coefficient matrix.  The correction is scaled by the factor
! 2/(1+RC) to account for changes in h*rl1 since the last ZVJAC call.
!-----------------------------------------------------------------------
350 DO 360 I = 1,N
  Y(I) = (RL1*H)*SAVF(I) - (RL1*YH(I,2) + ACOR(I))
360 CONTINUE
CALL ZVSOL (WM, IWM, Y, IERSL)
NNI = NNI + 1
IF (IERSL .GT. 0) GO TO 410
IF (METH .EQ. 2 .AND. RC .NE. ONE) THEN
  CSCALE = TWO/(ONE + RC)
  CALL DZSCAL (N, CSCALE, Y, 1)
ENDIF
DEL = ZVNORM (N, Y, EWT)
CALL DZAXPY (N, ONE, Y, 1, ACOR, 1)
DO 380 I = 1,N
  Y(I) = YH(I,1) + ACOR(I)
380 CONTINUE
!-----------------------------------------------------------------------
! Test for convergence.  If M .gt. 0, an estimate of the convergence
! rate constant is stored in CRATE, and this is used in the test.
!-----------------------------------------------------------------------
400 IF (M .NE. 0) CRATE = MAX(CRDOWN*CRATE,DEL/DELP)
DCON = DEL*MIN(ONE,CRATE)/TQ(4)
IF (DCON .LE. ONE) GO TO 450
M = M + 1
IF (M .EQ. MAXCOR) GO TO 410
IF (M .GE. 2 .AND. DEL .GT. RDIV*DELP) GO TO 410
DELP = DEL
CALL F (N, TN, Y, SAVF, RPAR, IPAR)
NFE = NFE + 1
GO TO 270
!
410 IF (MITER .EQ. 0 .OR. JCUR .EQ. 1) GO TO 430
ICF = 1
IPUP = MITER
GO TO 220
!
430 CONTINUE
NFLAG = -1
ICF = 2
IPUP = MITER
RETURN
!
! Return for successful step. ------------------------------------------
450 NFLAG = 0
JCUR = 0
ICF = 0
IF (M .EQ. 0) ACNRM = DEL
IF (M .GT. 0) ACNRM = ZVNORM (N, ACOR, EWT)
RETURN
!----------------------- End of Subroutine ZVNLSD ----------------------
END
!DECK ZVJAC
SUBROUTINE ZVJAC (Y, YH, LDYH, EWT, FTEM, SAVF, WM, IWM, F, JAC, &
& IERPJ, RPAR, IPAR)
EXTERNAL F, JAC
COMPLEX(KIND=KIND(0.0d0)) Y, YH, FTEM, SAVF, WM
DOUBLE PRECISION EWT
INTEGER LDYH, IWM, IERPJ, IPAR
DIMENSION Y(*), YH(LDYH,*), EWT(*), FTEM(*), SAVF(*), &
& WM(*), IWM(*), RPAR(*), IPAR(*)
!-----------------------------------------------------------------------
! Call sequence input -- Y, YH, LDYH, EWT, FTEM, SAVF, WM, IWM,
!                        F, JAC, RPAR, IPAR
! Call sequence output -- WM, IWM, IERPJ
! COMMON block variables accessed:
!     /ZVOD01/  CCMXJ, DRC, H, HRL1, RL1, SRUR, TN, UROUND, ICF, JCUR,
!               LOCJS, MITER, MSBJ, N, NSLJ
!     /ZVOD02/  NFE, NST, NJE, NLU
!
! Subroutines called by ZVJAC: F, JAC, ZACOPY, ZCOPY, ZGBFA, ZGEFA,
!                              DZSCAL
! Function routines called by ZVJAC: ZVNORM
!-----------------------------------------------------------------------
! ZVJAC is called by ZVNLSD to compute and process the matrix
! P = I - h*rl1*J , where J is an approximation to the Jacobian.
! Here J is computed by the user-supplied routine JAC if
! MITER = 1 or 4, or by finite differencing if MITER = 2, 3, or 5.
! If MITER = 3, a diagonal approximation to J is used.
! If JSV = -1, J is computed from scratch in all cases.
! If JSV = 1 and MITER = 1, 2, 4, or 5, and if the saved value of J is
! considered acceptable, then P is constructed from the saved J.
! J is stored in wm and replaced by P.  If MITER .ne. 3, P is then
! subjected to LU decomposition in preparation for later solution
! of linear systems with P as coefficient matrix. This is done
! by ZGEFA if MITER = 1 or 2, and by ZGBFA if MITER = 4 or 5.
!
! Communication with ZVJAC is done with the following variables.  (For
! more details, please see the comments in the driver subroutine.)
! Y          = Vector containing predicted values on entry.
! YH         = The Nordsieck array, an LDYH by LMAX array, input.
! LDYH       = A constant .ge. N, the first dimension of YH, input.
! EWT        = An error weight vector of length N.
! SAVF       = Array containing f evaluated at predicted y, input.
! WM         = Complex work space for matrices.  In the output, it
!              contains the inverse diagonal matrix if MITER = 3 and
!              the LU decomposition of P if MITER is 1, 2 , 4, or 5.
!              Storage of the saved Jacobian starts at WM(LOCJS).
! IWM        = Integer work space containing pivot information,
!              starting at IWM(31), if MITER is 1, 2, 4, or 5.
!              IWM also contains band parameters ML = IWM(1) and
!              MU = IWM(2) if MITER is 4 or 5.
! F          = Dummy name for the user-supplied subroutine for f.
! JAC        = Dummy name for the user-supplied Jacobian subroutine.
! RPAR, IPAR = User's real/complex and integer work arrays.
! RL1        = 1/EL(2) (input).
! IERPJ      = Output error flag,  = 0 if no trouble, 1 if the P
!              matrix is found to be singular.
! JCUR       = Output flag to indicate whether the Jacobian matrix
!              (or approximation) is now current.
!              JCUR = 0 means J is not current.
!              JCUR = 1 means J is current.
!-----------------------------------------------------------------------
!
! Type declarations for labeled COMMON block ZVOD01 --------------------
!
DOUBLE PRECISION ACNRM, CCMXJ, CONP, CRATE, DRC, EL, &
& ETA, ETAMAX, H, HMIN, HMXI, HNEW, HRL1, HSCAL, PRL1, &
& RC, RL1, SRUR, TAU, TQ, TN, UROUND
INTEGER ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH, &
& L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, &
& LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP, &
& N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ, &
& NSLP, NYH
!
! Type declarations for labeled COMMON block ZVOD02 --------------------
!
DOUBLE PRECISION HU
INTEGER NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
!
! Type declarations for local variables --------------------------------
!
COMPLEX(KIND=KIND(0.0d0)) DI, R1, YI, YJ, YJJ
DOUBLE PRECISION CON, FAC, ONE, PT1, R, R0, THOU, ZERO
INTEGER I, I1, I2, IER, II, J, J1, JJ, JOK, LENP, MBA, MBAND, &
& MEB1, MEBAND, ML, ML1, MU, NP1
!
! Type declaration for function subroutines called ---------------------
!
DOUBLE PRECISION ZVNORM
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this subroutine.
!-----------------------------------------------------------------------
SAVE ONE, PT1, THOU, ZERO
!-----------------------------------------------------------------------
COMMON /ZVOD01/ ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13), ETA, &
& ETAMAX, H, HMIN, HMXI, HNEW, HRL1, HSCAL, PRL1, &
& RC, RL1, SRUR, TAU(13), TQ(5), TN, UROUND, &
& ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH, &
& L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, &
& LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP, &
& N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ, &
& NSLP, NYH
COMMON /ZVOD02/ HU, NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
!
DATA ONE /1.0D0/, THOU /1000.0D0/, ZERO /0.0D0/, PT1 /0.1D0/
!
IERPJ = 0
HRL1 = H*RL1
! See whether J should be evaluated (JOK = -1) or not (JOK = 1). -------
JOK = JSV
IF (JSV .EQ. 1) THEN
  IF (NST .EQ. 0 .OR. NST .GT. NSLJ+MSBJ) JOK = -1
  IF (ICF .EQ. 1 .AND. DRC .LT. CCMXJ) JOK = -1
  IF (ICF .EQ. 2) JOK = -1
ENDIF
! End of setting JOK. --------------------------------------------------
!
IF (JOK .EQ. -1 .AND. MITER .EQ. 1) THEN
! If JOK = -1 and MITER = 1, call JAC to evaluate Jacobian. ------------
NJE = NJE + 1
NSLJ = NST
JCUR = 1
LENP = N*N
DO 110 I = 1,LENP
  WM(I) = ZERO
110 CONTINUE
CALL JAC (N, TN, Y, 0, 0, WM, N, RPAR, IPAR)
IF (JSV .EQ. 1) CALL ZCOPY (LENP, WM, 1, WM(LOCJS), 1)
ENDIF
!
IF (JOK .EQ. -1 .AND. MITER .EQ. 2) THEN
! If MITER = 2, make N calls to F to approximate the Jacobian. ---------
NJE = NJE + 1
NSLJ = NST
JCUR = 1
FAC = ZVNORM (N, SAVF, EWT)
R0 = THOU*ABS(H)*UROUND*REAL(N)*FAC
IF (R0 .EQ. ZERO) R0 = ONE
J1 = 0
DO 230 J = 1,N
  YJ = Y(J)
  R = MAX(SRUR*ABS(YJ),R0/EWT(J))
  Y(J) = Y(J) + R
  FAC = ONE/R
  CALL F (N, TN, Y, FTEM, RPAR, IPAR)
  DO 220 I = 1,N
    WM(I+J1) = (FTEM(I) - SAVF(I))*FAC
220 CONTINUE
  Y(J) = YJ
  J1 = J1 + N
230 CONTINUE
NFE = NFE + N
LENP = N*N
IF (JSV .EQ. 1) CALL ZCOPY (LENP, WM, 1, WM(LOCJS), 1)
ENDIF
!
IF (JOK .EQ. 1 .AND. (MITER .EQ. 1 .OR. MITER .EQ. 2)) THEN
JCUR = 0
LENP = N*N
CALL ZCOPY (LENP, WM(LOCJS), 1, WM, 1)
ENDIF
!
IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
! Multiply Jacobian by scalar, add identity, and do LU decomposition. --
CON = -HRL1
CALL DZSCAL (LENP, CON, WM, 1)
J = 1
NP1 = N + 1
DO 250 I = 1,N
  WM(J) = WM(J) + ONE
  J = J + NP1
250 CONTINUE
NLU = NLU + 1
CALL ZGEFA (WM, N, N, IWM(31), IER)
IF (IER .NE. 0) IERPJ = 1
RETURN
ENDIF
! End of code block for MITER = 1 or 2. --------------------------------
!
IF (MITER .EQ. 3) THEN
! If MITER = 3, construct a diagonal approximation to J and P. ---------
NJE = NJE + 1
JCUR = 1
R = RL1*PT1
DO 310 I = 1,N
  Y(I) = Y(I) + R*(H*SAVF(I) - YH(I,2))
310 CONTINUE
CALL F (N, TN, Y, WM, RPAR, IPAR)
NFE = NFE + 1
DO 320 I = 1,N
  R1 = H*SAVF(I) - YH(I,2)
  DI = PT1*R1 - H*(WM(I) - SAVF(I))
  WM(I) = ONE
  IF (ABS(R1) .LT. UROUND/EWT(I)) GO TO 320
  IF (ABS(DI) .EQ. ZERO) GO TO 330
  WM(I) = PT1*R1/DI
320 CONTINUE
RETURN
330 IERPJ = 1
RETURN
ENDIF
! End of code block for MITER = 3. -------------------------------------
!
! Set constants for MITER = 4 or 5. ------------------------------------
ML = IWM(1)
MU = IWM(2)
ML1 = ML + 1
MBAND = ML + MU + 1
MEBAND = MBAND + ML
LENP = MEBAND*N
!
IF (JOK .EQ. -1 .AND. MITER .EQ. 4) THEN
! If JOK = -1 and MITER = 4, call JAC to evaluate Jacobian. ------------
NJE = NJE + 1
NSLJ = NST
JCUR = 1
DO 410 I = 1,LENP
  WM(I) = ZERO
410 CONTINUE
CALL JAC (N, TN, Y, ML, MU, WM(ML1), MEBAND, RPAR, IPAR)
IF (JSV .EQ. 1) &
& CALL ZACOPY (MBAND, N, WM(ML1), MEBAND, WM(LOCJS), MBAND)
ENDIF
!
IF (JOK .EQ. -1 .AND. MITER .EQ. 5) THEN
! If MITER = 5, make ML+MU+1 calls to F to approximate the Jacobian. ---
NJE = NJE + 1
NSLJ = NST
JCUR = 1
MBA = MIN(MBAND,N)
MEB1 = MEBAND - 1
FAC = ZVNORM (N, SAVF, EWT)
R0 = THOU*ABS(H)*UROUND*REAL(N)*FAC
IF (R0 .EQ. ZERO) R0 = ONE
DO 560 J = 1,MBA
  DO 530 I = J,N,MBAND
    YI = Y(I)
    R = MAX(SRUR*ABS(YI),R0/EWT(I))
    Y(I) = Y(I) + R
530 CONTINUE
  CALL F (N, TN, Y, FTEM, RPAR, IPAR)
  DO 550 JJ = J,N,MBAND
    Y(JJ) = YH(JJ,1)
    YJJ = Y(JJ)
    R = MAX(SRUR*ABS(YJJ),R0/EWT(JJ))
    FAC = ONE/R
    I1 = MAX(JJ-MU,1)
    I2 = MIN(JJ+ML,N)
    II = JJ*MEB1 - ML
    DO 540 I = I1,I2
      WM(II+I) = (FTEM(I) - SAVF(I))*FAC
540 CONTINUE
550 CONTINUE
560 CONTINUE
NFE = NFE + MBA
IF (JSV .EQ. 1) &
& CALL ZACOPY (MBAND, N, WM(ML1), MEBAND, WM(LOCJS), MBAND)
ENDIF
!
IF (JOK .EQ. 1) THEN
JCUR = 0
CALL ZACOPY (MBAND, N, WM(LOCJS), MBAND, WM(ML1), MEBAND)
ENDIF
!
! Multiply Jacobian by scalar, add identity, and do LU decomposition.
CON = -HRL1
CALL DZSCAL (LENP, CON, WM, 1 )
II = MBAND
DO 580 I = 1,N
  WM(II) = WM(II) + ONE
  II = II + MEBAND
580 CONTINUE
NLU = NLU + 1
CALL ZGBFA (WM, MEBAND, N, ML, MU, IWM(31), IER)
IF (IER .NE. 0) IERPJ = 1
RETURN
! End of code block for MITER = 4 or 5. --------------------------------
!
!----------------------- End of Subroutine ZVJAC -----------------------
END
!DECK ZACOPY
SUBROUTINE ZACOPY (NROW, NCOL, A, NROWA, B, NROWB)
COMPLEX(KIND=KIND(0.0d0)) A, B
INTEGER NROW, NCOL, NROWA, NROWB
DIMENSION A(NROWA,NCOL), B(NROWB,NCOL)
!-----------------------------------------------------------------------
! Call sequence input -- NROW, NCOL, A, NROWA, NROWB
! Call sequence output -- B
! COMMON block variables accessed -- None
!
! Subroutines called by ZACOPY: ZCOPY
! Function routines called by ZACOPY: None
!-----------------------------------------------------------------------
! This routine copies one rectangular array, A, to another, B,
! where A and B may have different row dimensions, NROWA and NROWB.
! The data copied consists of NROW rows and NCOL columns.
!-----------------------------------------------------------------------
INTEGER IC
!
DO 20 IC = 1,NCOL
  CALL ZCOPY (NROW, A(1,IC), 1, B(1,IC), 1)
20 CONTINUE
!
RETURN
!----------------------- End of Subroutine ZACOPY ----------------------
END
!DECK ZVSOL
SUBROUTINE ZVSOL (WM, IWM, X, IERSL)
COMPLEX(KIND=KIND(0.0d0)) WM, X
INTEGER IWM, IERSL
DIMENSION WM(*), IWM(*), X(*)
!-----------------------------------------------------------------------
! Call sequence input -- WM, IWM, X
! Call sequence output -- X, IERSL
! COMMON block variables accessed:
!     /ZVOD01/ -- H, HRL1, RL1, MITER, N
!
! Subroutines called by ZVSOL: ZGESL, ZGBSL
! Function routines called by ZVSOL: None
!-----------------------------------------------------------------------
! This routine manages the solution of the linear system arising from
! a chord iteration.  It is called if MITER .ne. 0.
! If MITER is 1 or 2, it calls ZGESL to accomplish this.
! If MITER = 3 it updates the coefficient H*RL1 in the diagonal
! matrix, and then computes the solution.
! If MITER is 4 or 5, it calls ZGBSL.
! Communication with ZVSOL uses the following variables:
! WM    = Real work space containing the inverse diagonal matrix if
!         MITER = 3 and the LU decomposition of the matrix otherwise.
! IWM   = Integer work space containing pivot information, starting at
!         IWM(31), if MITER is 1, 2, 4, or 5.  IWM also contains band
!         parameters ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.
! X     = The right-hand side vector on input, and the solution vector
!         on output, of length N.
! IERSL = Output flag.  IERSL = 0 if no trouble occurred.
!         IERSL = 1 if a singular matrix arose with MITER = 3.
!-----------------------------------------------------------------------
!
! Type declarations for labeled COMMON block ZVOD01 --------------------
!
DOUBLE PRECISION ACNRM, CCMXJ, CONP, CRATE, DRC, EL, &
& ETA, ETAMAX, H, HMIN, HMXI, HNEW, HRL1, HSCAL, PRL1, &
& RC, RL1, SRUR, TAU, TQ, TN, UROUND
INTEGER ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH, &
& L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, &
& LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP, &
& N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ, &
& NSLP, NYH
!
! Type declarations for local variables --------------------------------
!
COMPLEX(KIND=KIND(0.0d0)) DI
DOUBLE PRECISION ONE, PHRL1, R, ZERO
INTEGER I, MEBAND, ML, MU
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
SAVE ONE, ZERO
!
COMMON /ZVOD01/ ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13), ETA, &
& ETAMAX, H, HMIN, HMXI, HNEW, HRL1, HSCAL, PRL1, &
& RC, RL1, SRUR, TAU(13), TQ(5), TN, UROUND, &
& ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH, &
& L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, &
& LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP, &
& N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ, &
& NSLP, NYH
!
DATA ONE /1.0D0/, ZERO /0.0D0/
!
IERSL = 0
IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
  GOTO 100
ELSE IF (MITER .EQ. 3) THEN
  GOTO 300
ELSE IF (MITER .LE. 5) THEN
  GOTO 400
ENDIF

!      GO TO (100, 100, 300, 400, 400), MITER
100 CALL ZGESL (WM, N, N, IWM(31), X, 0)
RETURN
!
300 PHRL1 = HRL1
HRL1 = H*RL1
IF (HRL1 .EQ. PHRL1) GO TO 330
R = HRL1/PHRL1
DO 320 I = 1,N
  DI = ONE - R*(ONE - ONE/WM(I))
  IF (ABS(DI) .EQ. ZERO) GO TO 390
  WM(I) = ONE/DI
320 CONTINUE
!
330 DO 340 I = 1,N
  X(I) = WM(I)*X(I)
340 CONTINUE
RETURN
390 IERSL = 1
RETURN
!
400 ML = IWM(1)
MU = IWM(2)
MEBAND = 2*ML + MU + 1
CALL ZGBSL (WM, MEBAND, N, ML, MU, IWM(31), X, 0)
RETURN
!----------------------- End of Subroutine ZVSOL -----------------------
END
!DECK ZVSRCO
SUBROUTINE ZVSRCO (RSAV, ISAV, JOB)
DOUBLE PRECISION RSAV
INTEGER ISAV, JOB
DIMENSION RSAV(*), ISAV(*)
!-----------------------------------------------------------------------
! Call sequence input -- RSAV, ISAV, JOB
! Call sequence output -- RSAV, ISAV
! COMMON block variables accessed -- All of /ZVOD01/ and /ZVOD02/
!
! Subroutines/functions called by ZVSRCO: None
!-----------------------------------------------------------------------
! This routine saves or restores (depending on JOB) the contents of the
! COMMON blocks ZVOD01 and ZVOD02, which are used internally by ZVODE.
!
! RSAV = real array of length 51 or more.
! ISAV = integer array of length 41 or more.
! JOB  = flag indicating to save or restore the COMMON blocks:
!        JOB  = 1 if COMMON is to be saved (written to RSAV/ISAV).
!        JOB  = 2 if COMMON is to be restored (read from RSAV/ISAV).
!        A call with JOB = 2 presumes a prior call with JOB = 1.
!-----------------------------------------------------------------------
DOUBLE PRECISION RVOD1, RVOD2
INTEGER IVOD1, IVOD2
INTEGER I, LENIV1, LENIV2, LENRV1, LENRV2
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
SAVE LENRV1, LENIV1, LENRV2, LENIV2
!
COMMON /ZVOD01/ RVOD1(50), IVOD1(33)
COMMON /ZVOD02/ RVOD2(1), IVOD2(8)
DATA LENRV1/50/, LENIV1/33/, LENRV2/1/, LENIV2/8/
!
IF (JOB .EQ. 2) GO TO 100
DO 10 I = 1,LENRV1
  RSAV(I) = RVOD1(I)
10 CONTINUE
DO 15 I = 1,LENRV2
  RSAV(LENRV1+I) = RVOD2(I)
15 CONTINUE
!
DO 20 I = 1,LENIV1
  ISAV(I) = IVOD1(I)
20 CONTINUE
DO 25 I = 1,LENIV2
  ISAV(LENIV1+I) = IVOD2(I)
25 CONTINUE
!
RETURN
!
100 CONTINUE
DO 110 I = 1,LENRV1
   RVOD1(I) = RSAV(I)
110 CONTINUE
DO 115 I = 1,LENRV2
   RVOD2(I) = RSAV(LENRV1+I)
115 CONTINUE
!
DO 120 I = 1,LENIV1
   IVOD1(I) = ISAV(I)
120 CONTINUE
DO 125 I = 1,LENIV2
   IVOD2(I) = ISAV(LENIV1+I)
125 CONTINUE
!
RETURN
!----------------------- End of Subroutine ZVSRCO ----------------------
END
!DECK ZEWSET
SUBROUTINE ZEWSET (N, ITOL, RTOL, ATOL, YCUR, EWT)
!***BEGIN PROLOGUE  ZEWSET
!***SUBSIDIARY
!***PURPOSE  Set error weight vector.
!***TYPE      DOUBLE PRECISION (SEWSET-S, DEWSET-D, ZEWSET-Z)
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!  This subroutine sets the error weight vector EWT according to
!      EWT(i) = RTOL(i)*ABS(YCUR(i)) + ATOL(i),  i = 1,...,N,
!  with the subscript on RTOL and/or ATOL possibly replaced by 1 above,
!  depending on the value of ITOL.
!
!***SEE ALSO  DLSODE
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   060502  DATE WRITTEN, modified from DEWSET of 930809.
!***END PROLOGUE  ZEWSET
COMPLEX(KIND=KIND(0.0d0)) YCUR
DOUBLE PRECISION RTOL, ATOL, EWT
INTEGER N, ITOL
INTEGER I
DIMENSION RTOL(*), ATOL(*), YCUR(N), EWT(N)
!
!***FIRST EXECUTABLE STATEMENT  ZEWSET
IF (ITOL .EQ. 1) THEN
  GOTO 10
ELSE IF (ITOL .EQ. 2) THEN
  GOTO 20
ELSE IF (ITOL .EQ. 3) THEN
  GOTO 30
ELSE IF (ITOL .EQ. 4) THEN
  GOTO 40
ENDIF
!      GO TO (10, 20, 30, 40), ITOL
10 CONTINUE
DO 15 I = 1,N
  EWT(I) = RTOL(1)*ABS(YCUR(I)) + ATOL(1)
15 CONTINUE
RETURN
20 CONTINUE
DO 25 I = 1,N
  EWT(I) = RTOL(1)*ABS(YCUR(I)) + ATOL(I)
25 CONTINUE
RETURN
30 CONTINUE
DO 35 I = 1,N
  EWT(I) = RTOL(I)*ABS(YCUR(I)) + ATOL(1)
35 CONTINUE
RETURN
40 CONTINUE
DO 45 I = 1,N
  EWT(I) = RTOL(I)*ABS(YCUR(I)) + ATOL(I)
45 CONTINUE
RETURN
!----------------------- END OF SUBROUTINE ZEWSET ----------------------
END
!DECK ZVNORM
DOUBLE PRECISION FUNCTION ZVNORM (N, V, W)
!***BEGIN PROLOGUE  ZVNORM
!***SUBSIDIARY
!***PURPOSE  Weighted root-mean-square vector norm.
!***TYPE      DOUBLE COMPLEX (SVNORM-S, DVNORM-D, ZVNORM-Z)
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!  This function routine computes the weighted root-mean-square norm
!  of the vector of length N contained in the double complex array V,
!  with weights contained in the array W of length N:
!    ZVNORM = SQRT( (1/N) * SUM( abs(V(i))**2 * W(i)**2 )
!  The squared absolute value abs(v)**2 is computed by ZABSSQ.
!
!***SEE ALSO  DLSODE
!***ROUTINES CALLED  ZABSSQ
!***REVISION HISTORY  (YYMMDD)
!   060502  DATE WRITTEN, modified from DVNORM of 930809.
!***END PROLOGUE  ZVNORM
COMPLEX(KIND=KIND(0.0d0)) V
DOUBLE PRECISION W,   SUM, ZABSSQ
INTEGER N,   I
DIMENSION V(N), W(N)
!
!***FIRST EXECUTABLE STATEMENT  ZVNORM
SUM = 0.0D0
DO 10 I = 1,N
  SUM = SUM + ZABSSQ(V(I)) * W(I)**2
10 CONTINUE
ZVNORM = SQRT(SUM/N)
RETURN
!----------------------- END OF FUNCTION ZVNORM ------------------------
END
!DECK ZABSSQ
DOUBLE PRECISION FUNCTION ZABSSQ(Z)
!***BEGIN PROLOGUE  ZABSSQ
!***SUBSIDIARY
!***PURPOSE  Squared absolute value of a double complex number.
!***TYPE      DOUBLE PRECISION (ZABSSQ-Z)
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!  This function routine computes the square of the absolute value of
!  a double precision complex number Z,
!    ZABSSQ = DREAL(Z)**2 * DIMAG(Z)**2
!***REVISION HISTORY  (YYMMDD)
!   060502  DATE WRITTEN.
!***END PROLOGUE  ZABSSQ
COMPLEX(KIND=KIND(0.0d0)) Z
ZABSSQ = REAL(Z)**2 + AIMAG(Z)**2
RETURN
!----------------------- END OF FUNCTION ZABSSQ ------------------------
END
!DECK DZSCAL
SUBROUTINE DZSCAL(N, DA, ZX, INCX)
!***BEGIN PROLOGUE  DZSCAL
!***SUBSIDIARY
!***PURPOSE  Scale a double complex vector by a double prec. constant.
!***TYPE      DOUBLE PRECISION (DZSCAL-Z)
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!  Scales a double complex vector by a double precision constant.
!  Minor modification of BLAS routine ZSCAL.
!***REVISION HISTORY  (YYMMDD)
!   060530  DATE WRITTEN.
!***END PROLOGUE  DZSCAL
COMPLEX(KIND=KIND(0.0d0)) ZX(*)
DOUBLE PRECISION DA
INTEGER I,INCX,IX,N
!
IF( N.LE.0 .OR. INCX.LE.0 )RETURN
IF(INCX.EQ.1)GO TO 20
! Code for increment not equal to 1
IX = 1
DO 10 I = 1,N
  ZX(IX) = DA*ZX(IX)
  IX = IX + INCX
10 CONTINUE
RETURN
! Code for increment equal to 1
20 DO 30 I = 1,N
  ZX(I) = DA*ZX(I)
30 CONTINUE
RETURN
END
!DECK DZAXPY
SUBROUTINE DZAXPY(N, DA, ZX, INCX, ZY, INCY)
!***BEGIN PROLOGUE  DZAXPY
!***PURPOSE  Real constant times a complex vector plus a complex vector.
!***TYPE      DOUBLE PRECISION (DZAXPY-Z)
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!  Add a D.P. real constant times a complex vector to a complex vector.
!  Minor modification of BLAS routine ZAXPY.
!***REVISION HISTORY  (YYMMDD)
!   060530  DATE WRITTEN.
!***END PROLOGUE  DZAXPY
COMPLEX(KIND=KIND(0.0d0)) ZX(*),ZY(*)
DOUBLE PRECISION DA
INTEGER I,INCX,INCY,IX,IY,N
IF(N.LE.0)RETURN
IF (ABS(DA) .EQ. 0.0D0) RETURN
IF (INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
! Code for unequal increments or equal increments not equal to 1
IX = 1
IY = 1
IF(INCX.LT.0)IX = (-N+1)*INCX + 1
IF(INCY.LT.0)IY = (-N+1)*INCY + 1
DO 10 I = 1,N
  ZY(IY) = ZY(IY) + DA*ZX(IX)
  IX = IX + INCX
  IY = IY + INCY
10 CONTINUE
RETURN
! Code for both increments equal to 1
20 DO 30 I = 1,N
  ZY(I) = ZY(I) + DA*ZX(I)
30 CONTINUE
RETURN
END


subroutine zgesl(a,lda,n,ipvt,b,job)
integer lda,n,ipvt(1),job
COMPLEX(KIND=KIND(0.0d0)) a(lda,*),b(*)
!
!     zgesl solves the COMPLEX(KIND=8) system
!     a * x = b  or  ctrans(a) * x = b
!     using the factors computed by zgeco or zgefa.
!
!     on entry
!
!        a       COMPLEX(KIND=8)(lda, n)
!                the output from zgeco or zgefa.
!
!        lda     integer
!                the leading dimension of the array  a .
!
!        n       integer
!                the order of the matrix  a .
!
!        ipvt    integer(n)
!                the pivot vector from zgeco or zgefa.
!
!        b       COMPLEX(KIND=8)(n)
!                the right hand side vector.
!
!        job     integer
!                = 0         to solve  a*x = b ,
!                = nonzero   to solve  ctrans(a)*x = b  where
!                            ctrans(a)  is the conjugate transpose.
!
!     on return
!
!        b       the solution vector  x .
!
!     error condition
!
!        a division by zero will occur if the input factor contains a
!        zero on the diagonal.  technically this indicates singularity
!        but it is often caused by improper arguments or improper
!        setting of lda .  it will not occur if the subroutines are
!        called correctly and if zgeco has set rcond .gt. 0.0
!        or zgefa has set info .eq. 0 .
!
!     to compute  inverse(a) * c  where  c  is a matrix
!     with  p  columns
!           call zgeco(a,lda,n,ipvt,rcond,z)
!           if (rcond is too small) go to ...
!           do 10 j = 1, p
!              call zgesl(a,lda,n,ipvt,c(1,j),0)
!        10 continue
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     blas zaxpy,zdotc
!     fortran dconjg
!
!     internal variables
!
COMPLEX(KIND=KIND(0.0d0)) zdotc,t
integer k,kb,l,nm1
! KS      double precision dreal,dimag
! KS      COMPLEX(KIND=8) zdumr,zdumi
! KS      dreal(zdumr) = zdumr
! KS      dimag(zdumi) = (0.0d0,-1.0d0)*zdumi
!
nm1 = n - 1
if (job .NE. 0) go to 50
!
!        job = 0 , solve  a * x = b
!        first solve  l*y = b
!
   if (nm1 .LT. 1) go to 30
   do 20 k = 1, nm1
      l = ipvt(k)
      t = b(l)
      if (l .EQ. k) go to 10
         b(l) = b(k)
         b(k) = t
10 continue
      call zaxpy(n-k,t,a(k+1,k),1,b(k+1),1)
20 continue
30 continue
!
!        now solve  u*x = y
!
   do 40 kb = 1, n
      k = n + 1 - kb
      b(k) = b(k)/a(k,k)
      t = -b(k)
      call zaxpy(k-1,t,a(1,k),1,b(1),1)
40 continue
go to 100
50 continue
!
!        job = nonzero, solve  ctrans(a) * x = b
!        first solve  ctrans(u)*y = b
!
   do 60 k = 1, n
      t = zdotc(k-1,a(1,k),1,b(1),1)
      b(k) = (b(k) - t)/conjg(a(k,k))
60 continue
!
!        now solve ctrans(l)*x = y
!
   if (nm1 .LT. 1) go to 90
   do 80 kb = 1, nm1
      k = n - kb
      b(k) = b(k) + zdotc(n-k,a(k+1,k),1,b(k+1),1)
      l = ipvt(k)
      if (l .EQ. k) go to 70
         t = b(l)
         b(l) = b(k)
         b(k) = t
70 continue
80 continue
90 continue
100 continue
return
end

subroutine zgbfa(abd,lda,n,ml,mu,ipvt,info)
integer lda,n,ml,mu,ipvt(*),info
COMPLEX(KIND=KIND(0.0d0)) abd(lda,*)
!
!     zgbfa factors a COMPLEX(KIND=8) band matrix by elimination.
!
!     zgbfa is usually called by zgbco, but it can be called
!     directly with a saving in time if  rcond  is not needed.
!
!     on entry
!
!        abd     COMPLEX(KIND=8)(lda, n)
!                contains the matrix in band storage.  the columns
!                of the matrix are stored in the columns of  abd  and
!                the diagonals of the matrix are stored in rows
!                ml+1 through 2*ml+mu+1 of  abd .
!                see the comments below for details.
!
!        lda     integer
!                the leading dimension of the array  abd .
!                lda must be .ge. 2*ml + mu + 1 .
!
!        n       integer
!                the order of the original matrix.
!
!        ml      integer
!                number of diagonals below the main diagonal.
!                0 .le. ml .lt. n .
!
!        mu      integer
!                number of diagonals above the main diagonal.
!                0 .le. mu .lt. n .
!                more efficient if  ml .le. mu .
!     on return
!
!        abd     an upper triangular matrix in band storage and
!                the multipliers which were used to obtain it.
!                the factorization can be written  a = l*u  where
!                l  is a product of permutation and unit lower
!                triangular matrices and  u  is upper triangular.
!
!        ipvt    integer(n)
!                an integer vector of pivot indices.
!
!        info    integer
!                = 0  normal value.
!                = k  if  u(k,k) .eq. 0.0 .  this is not an error
!                     condition for this subroutine, but it does
!                     indicate that zgbsl will divide by zero if
!                     called.  use  rcond  in zgbco for a reliable
!                     indication of singularity.
!
!     band storage
!
!           if  a  is a band matrix, the following program segment
!           will set up the input.
!
!                   ml = (band width below the diagonal)
!                   mu = (band width above the diagonal)
!                   m = ml + mu + 1
!                   do 20 j = 1, n
!                      i1 = max0(1, j-mu)
!                      i2 = min0(n, j+ml)
!                      do 10 i = i1, i2
!                         k = i - j + m
!                         abd(k,j) = a(i,j)
!                10    continue
!                20 continue
!
!           this uses rows  ml+1  through  2*ml+mu+1  of  abd .
!           in addition, the first  ml  rows in  abd  are used for
!           elements generated during the triangularization.
!           the total number of rows needed in  abd  is  2*ml+mu+1 .
!           the  ml+mu by ml+mu  upper left triangle and the
!           ml by ml  lower right triangle are not referenced.
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     blas zaxpy,zscal,izamax
!     fortran dabs,max0,min0
!
!     internal variables
!
COMPLEX(KIND=KIND(0.0d0)) t
integer i,izamax,i0,j,ju,jz,j0,j1,k,kp1,l,lm,m,mm,nm1
!
!KS      COMPLEX(KIND=8) zdum
double precision cabs1
!      double precision dreal,dimag
!      COMPLEX(KIND=8) zdumr,zdumi
!      dreal(zdumr) = zdumr
!      dimag(zdumi) = (0.0d0,-1.0d0)*zdumi
!KS      cabs1(zdum) = dabs(dreal(zdum)) + dabs(dimag(zdum))
!
m = ml + mu + 1
info = 0
!
!     zero initial fill-in columns
!
j0 = mu + 2
j1 = min0(n,m) - 1
if (j1 .LT. j0) go to 30
do 20 jz = j0, j1
   i0 = m + 1 - jz
   do 10 i = i0, ml
      abd(i,jz) = (0.0d0,0.0d0)
10 continue
20 continue
30 continue
jz = j1
ju = 0
!
!     gaussian elimination with partial pivoting
!
nm1 = n - 1
if (nm1 .LT. 1) go to 130
do 120 k = 1, nm1
   kp1 = k + 1
!
!        zero next fill-in column
!
   jz = jz + 1
   if (jz .GT. n) go to 50
   if (ml .LT. 1) go to 50
      do 40 i = 1, ml
         abd(i,jz) = (0.0d0,0.0d0)
40 continue
50 continue
!
!        find l = pivot index
!
   lm = min0(ml,n-k)
   l = izamax(lm+1,abd(m,k),1) + m - 1
   ipvt(k) = l + k - m
!
!        zero pivot implies this column already triangularized
!
   if (cabs1(abd(l,k)) .EQ. 0.0d0) go to 100
!
!           interchange if necessary
!
      if (l .EQ. m) go to 60
         t = abd(l,k)
         abd(l,k) = abd(m,k)
         abd(m,k) = t
60 continue
!
!           compute multipliers
!
      t = -(1.0d0,0.0d0)/abd(m,k)
      call zscal(lm,t,abd(m+1,k),1)
!
!           row elimination with column indexing
!
      ju = min0(max0(ju,mu+ipvt(k)),n)
      mm = m
      if (ju .LT. kp1) go to 90
      do 80 j = kp1, ju
         l = l - 1
         mm = mm - 1
         t = abd(l,j)
         if (l .EQ. mm) go to 70
            abd(l,j) = abd(mm,j)
            abd(mm,j) = t
70 continue
         call zaxpy(lm,t,abd(m+1,k),1,abd(mm+1,j),1)
80 continue
90 continue
   go to 110
100 continue
      info = k
110 continue
120 continue
130 continue
ipvt(n) = n
if (cabs1(abd(m,n)) .EQ. 0.0d0) info = n
return
end

subroutine zgbsl(abd,lda,n,ml,mu,ipvt,b,job)
integer lda,n,ml,mu,ipvt(1),job
COMPLEX(KIND=KIND(0.0d0)) abd(lda,*),b(*)
!
!     zgbsl solves the COMPLEX(KIND=8) band system
!     a * x = b  or  ctrans(a) * x = b
!     using the factors computed by zgbco or zgbfa.
!
!     on entry
!
!        abd     COMPLEX(KIND=8)(lda, n)
!                the output from zgbco or zgbfa.
!
!        lda     integer
!                the leading dimension of the array  abd .
!
!        n       integer
!                the order of the original matrix.
!
!        ml      integer
!                number of diagonals below the main diagonal.
!
!        mu      integer
!                number of diagonals above the main diagonal.
!
!        ipvt    integer(n)
!                the pivot vector from zgbco or zgbfa.
!
!        b       COMPLEX(KIND=8)(n)
!                the right hand side vector.
!
!        job     integer
!                = 0         to solve  a*x = b ,
!                = nonzero   to solve  ctrans(a)*x = b , where
!                            ctrans(a)  is the conjugate transpose.
!
!     on return
!
!        b       the solution vector  x .
!
!     error condition
!
!        a division by zero will occur if the input factor contains a
!        zero on the diagonal.  technically this indicates singularity
!        but it is often caused by improper arguments or improper
!        setting of lda .  it will not occur if the subroutines are
!        called correctly and if zgbco has set rcond .gt. 0.0
!        or zgbfa has set info .eq. 0 .
!
!     to compute  inverse(a) * c  where  c  is a matrix
!     with  p  columns
!           call zgbco(abd,lda,n,ml,mu,ipvt,rcond,z)
!           if (rcond is too small) go to ...
!           do 10 j = 1, p
!              call zgbsl(abd,lda,n,ml,mu,ipvt,c(1,j),0)
!        10 continue
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     blas zaxpy,zdotc
!     fortran dconjg,min0
!
!     internal variables
!
COMPLEX(KIND=KIND(0.0d0)) zdotc,t
integer k,kb,l,la,lb,lm,m,nm1
!      double precision dreal,dimag
!      COMPLEX(KIND=8) zdumr,zdumi
!      dreal(zdumr) = zdumr
!      dimag(zdumi) = (0.0d0,-1.0d0)*zdumi
!
m = mu + ml + 1
nm1 = n - 1
if (job .NE. 0) go to 50
!
!        job = 0 , solve  a * x = b
!        first solve l*y = b
!
   if (ml .EQ. 0) go to 30
   if (nm1 .LT. 1) go to 30
      do 20 k = 1, nm1
         lm = min0(ml,n-k)
         l = ipvt(k)
         t = b(l)
         if (l .EQ. k) go to 10
            b(l) = b(k)
            b(k) = t
10 continue
         call zaxpy(lm,t,abd(m+1,k),1,b(k+1),1)
20 continue
30 continue
!
!        now solve  u*x = y
!
   do 40 kb = 1, n
      k = n + 1 - kb
      b(k) = b(k)/abd(m,k)
      lm = min0(k,m) - 1
      la = m - lm
      lb = k - lm
      t = -b(k)
      call zaxpy(lm,t,abd(la,k),1,b(lb),1)
40 continue
go to 100
50 continue
!
!        job = nonzero, solve  ctrans(a) * x = b
!        first solve  ctrans(u)*y = b
!
   do 60 k = 1, n
      lm = min0(k,m) - 1
      la = m - lm
      lb = k - lm
      t = zdotc(lm,abd(la,k),1,b(lb),1)
      b(k) = (b(k) - t)/conjg(abd(m,k))
60 continue
!
!        now solve ctrans(l)*x = y
!
   if (ml .EQ. 0) go to 90
   if (nm1 .LT. 1) go to 90
      do 80 kb = 1, nm1
         k = n - kb
         lm = min0(ml,n-k)
         b(k) = b(k) + zdotc(lm,abd(m+1,k),1,b(k+1),1)
         l = ipvt(k)
         if (l .EQ. k) go to 70
            t = b(l)
            b(l) = b(k)
            b(k) = t
70 continue
80 continue
90 continue
100 continue
return
end
! KARLINE: created true functions out of these statement functions
! Thomas:  removed function definitions for dreal and dimag,
!          as 'real' and 'aimag' are now standards.
double precision function cabs1(zdum)
COMPLEX(KIND=KIND(0.0d0)) zdum
  cabs1 = dabs(real(zdum)) + dabs(aimag(zdum))
end function
! KARLINE: end new functions


subroutine zgefa(a,lda,n,ipvt,info)
integer lda,n,ipvt(*),info
COMPLEX(KIND=KIND(0.0d0)) a(lda,*)
!
!     zgefa factors a COMPLEX(KIND=8) matrix by gaussian elimination.
!
!     zgefa is usually called by zgeco, but it can be called
!     directly with a saving in time if  rcond  is not needed.
!     (time for zgeco) = (1 + 9/n)*(time for zgefa) .
!
!     on entry
!
!        a       COMPLEX(KIND=8)(lda, n)
!                the matrix to be factored.
!
!        lda     integer
!                the leading dimension of the array  a .
!
!        n       integer
!                the order of the matrix  a .
!
!     on return
!
!        a       an upper triangular matrix and the multipliers
!                which were used to obtain it.
!                the factorization can be written  a = l*u  where
!                l  is a product of permutation and unit lower
!                triangular matrices and  u  is upper triangular.
!
!        ipvt    integer(n)
!                an integer vector of pivot indices.
!
!        info    integer
!                = 0  normal value.
!                = k  if  u(k,k) .eq. 0.0 .  this is not an error
!                     condition for this subroutine, but it does
!                     indicate that zgesl or zgedi will divide by zero
!                     if called.  use  rcond  in zgeco for a reliable
!                     indication of singularity.
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     blas zaxpy,zscal,izamax
!     fortran dabs
!
!     internal variables
!
COMPLEX(KIND=KIND(0.0d0)) t
integer izamax,j,k,kp1,l,nm1
!
! KS    COMPLEX(KIND=8) zdum
double precision cabs1
!     double precision dreal,dimag
! KS    COMPLEX(KIND=8) zdumr,zdumi
! Karline: next three statement functions replaced with true functions above
!      dreal(zdumr) = zdumr
!      dimag(zdumi) = (0.0d0,-1.0d0)*zdumi
!      cabs1(zdum) = dabs(dreal(zdum)) + dabs(dimag(zdum))
!
!     gaussian elimination with partial pivoting
!
info = 0
nm1 = n - 1
if (nm1 .LT. 1) go to 70
do 60 k = 1, nm1
   kp1 = k + 1
!
!        find l = pivot index
!
   l = izamax(n-k+1,a(k,k),1) + k - 1
   ipvt(k) = l
!
!        zero pivot implies this column already triangularized
!
   if (cabs1(a(l,k)) .EQ. 0.0d0) go to 40
!
!           interchange if necessary
!
      if (l .EQ. k) go to 10
         t = a(l,k)
         a(l,k) = a(k,k)
         a(k,k) = t
10 continue
!
!           compute multipliers
!
      t = -(1.0d0,0.0d0)/a(k,k)
      call zscal(n-k,t,a(k+1,k),1)
!
!           row elimination with column indexing
!
      do 30 j = kp1, n
         t = a(l,j)
         if (l .EQ. k) go to 20
            a(l,j) = a(k,j)
            a(k,j) = t
20 continue
         call zaxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
30 continue
   go to 50
40 continue
      info = k
50 continue
60 continue
70 continue
ipvt(n) = n
if (cabs1(a(n,n)) .EQ. 0.0d0) info = n
return
end

