!  The code in this file is was taken from daspk.tgz from
!    https://www.netlib.org/ode/
!  Original authors: Linda R. Petzold, Peter N. Brown, 
!    Alan C. Hindmarsh, and Clement W. Ulrich
!  Adapted for use in R package deSolve by the deSolve authors.
!

! Work performed under the auspices of the U.S. Department of Energy
! by Lawrence Livermore National Laboratory under contract number 
! W-7405-Eng-48.
! karline: changed INFO, to also pass the index of the variables
!          error scaling ~ index of variables            
SUBROUTINE DDASPK (RES, NEQ, T, Y, YPRIME, TOUT, INFO, RTOL, ATOL, &
& IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC, PSOL)
!
!***BEGIN PROLOGUE  DDASPK
!***DATE WRITTEN   890101   (YYMMDD)
!***REVISION DATE  910624   
!***REVISION DATE  920929   (CJ in RES call, RES counter fix.)
!***REVISION DATE  921215   (Warnings on poor iteration performance)
!***REVISION DATE  921216   (NRMAX as optional input)
!***REVISION DATE  930315   (Name change: DDINI to DDINIT)
!***REVISION DATE  940822   (Replaced initial condition calculation)
!***REVISION DATE  941101   (Added linesearch in I.C. calculations)
!***REVISION DATE  941220   (Misc. corrections throughout)
!***REVISION DATE  950125   (Added DINVWT routine)
!***REVISION DATE  950714   (Misc. corrections throughout)
!***REVISION DATE  950802   (Default NRMAX = 5, based on tests.)
!***REVISION DATE  950808   (Optional error test added.)
!***REVISION DATE  950814   (Added I.C. constraints and INFO(14))
!***REVISION DATE  950828   (Various minor corrections.)
!***REVISION DATE  951006   (Corrected WT scaling in DFNRMK.)
!***REVISION DATE  960129   (Corrected RL bug in DLINSD, DLINSK.)
!***REVISION DATE  960301   (Added NONNEG to SAVE statement.)
!***CATEGORY NO.  I1A2
!***KEYWORDS  DIFFERENTIAL/ALGEBRAIC, BACKWARD DIFFERENTIATION FORMULAS,
!             IMPLICIT DIFFERENTIAL SYSTEMS, KRYLOV ITERATION
!***AUTHORS   Linda R. Petzold, Peter N. Brown, Alan C. Hindmarsh, and
!                  Clement W. Ulrich
!             Center for Computational Sciences & Engineering, L-316
!             Lawrence Livermore National Laboratory
!             P.O. Box 808,
!             Livermore, CA 94551
!***PURPOSE  This code solves a system of differential/algebraic 
!            equations of the form 
!               G(t,y,y') = 0 , 
!            using a combination of Backward Differentiation Formula 
!            (BDF) methods and a choice of two linear system solution 
!            methods: direct (dense or band) or Krylov (iterative).
!            This version is in double precision.
!-----------------------------------------------------------------------
!***DESCRIPTION
!
! *Usage:
!
!      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!      INTEGER NEQ, INFO(N), IDID, LRW, LIW, IWORK(LIW), IPAR(*)
!      DOUBLE PRECISION T, Y(*), YPRIME(*), TOUT, RTOL(*), ATOL(*),
!         RWORK(LRW), RPAR(*)
!      EXTERNAL  RES, JAC, PSOL
!
!      CALL DDASPK (RES, NEQ, T, Y, YPRIME, TOUT, INFO, RTOL, ATOL,
!     *   IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC, PSOL)
!
!  Quantities which may be altered by the code are:
!     T, Y(*), YPRIME(*), INFO(1), RTOL, ATOL, IDID, RWORK(*), IWORK(*)
!
!
! *Arguments:
!
!  RES:EXT          This is the name of a subroutine which you
!                   provide to define the residual function G(t,y,y')
!                   of the differential/algebraic system.
!
!  NEQ:IN           This is the number of equations in the system.
!
!  T:INOUT          This is the current value of the independent 
!                   variable.
!
!  Y(*):INOUT       This array contains the solution components at T.
!
!  YPRIME(*):INOUT  This array contains the derivatives of the solution
!                   components at T.
!
!  TOUT:IN          This is a point at which a solution is desired.
!
!  INFO(N):IN       This is an integer array used to communicate details
!                   of how the solution is to be carried out, such as
!                   tolerance type, matrix structure, step size and
!                   order limits, and choice of nonlinear system method.
!                   N must be at least 20.
!
!  RTOL,ATOL:INOUT  These quantities represent absolute and relative
!                   error tolerances (on local error) which you provide
!                   to indicate how accurately you wish the solution to
!                   be computed.  You may choose them to be both scalars
!                   or else both arrays of length NEQ.
!
!  IDID:OUT         This integer scalar is an indicator reporting what
!                   the code did.  You must monitor this variable to
!                   decide what action to take next.
!
!  RWORK:WORK       A real work array of length LRW which provides the
!                   code with needed storage space.
!
!  LRW:IN           The length of RWORK.
!
!  IWORK:WORK       An integer work array of length LIW which provides
!                   the code with needed storage space.
!
!  LIW:IN           The length of IWORK.
!
!  RPAR,IPAR:IN     These are real and integer parameter arrays which
!                   you can use for communication between your calling
!                   program and the RES, JAC, and PSOL subroutines.
!
!  JAC:EXT          This is the name of a subroutine which you may
!                   provide (optionally) for calculating Jacobian 
!                   (partial derivative) data involved in solving linear
!                   systems within DDASPK.
!
!  PSOL:EXT         This is the name of a subroutine which you must
!                   provide for solving linear systems if you selected
!                   a Krylov method.  The purpose of PSOL is to solve
!                   linear systems involving a left preconditioner P.
!
! *Overview
!
!  The DDASPK solver uses the backward differentiation formulas of
!  orders one through five to solve a system of the form G(t,y,y') = 0
!  for y = Y and y' = YPRIME.  Values for Y and YPRIME at the initial 
!  time must be given as input.  These values should be consistent, 
!  that is, if T, Y, YPRIME are the given initial values, they should 
!  satisfy G(T,Y,YPRIME) = 0.  However, if consistent values are not
!  known, in many cases you can have DDASPK solve for them -- see INFO(11).
!  (This and other options are described in more detail below.)
!
!  Normally, DDASPK solves the system from T to TOUT.  It is easy to
!  continue the solution to get results at additional TOUT.  This is
!  the interval mode of operation.  Intermediate results can also be
!  obtained easily by specifying INFO(3).
!
!  On each step taken by DDASPK, a sequence of nonlinear algebraic  
!  systems arises.  These are solved by one of two types of
!  methods:
!    * a Newton iteration with a direct method for the linear
!      systems involved (INFO(12) = 0), or
!    * a Newton iteration with a preconditioned Krylov iterative 
!      method for the linear systems involved (INFO(12) = 1).
!
!  The direct method choices are dense and band matrix solvers, 
!  with either a user-supplied or an internal difference quotient 
!  Jacobian matrix, as specified by INFO(5) and INFO(6).
!  In the band case, INFO(6) = 1, you must supply half-bandwidths
!  in IWORK(1) and IWORK(2).
!
!  The Krylov method is the Generalized Minimum Residual (GMRES) 
!  method, in either complete or incomplete form, and with 
!  scaling and preconditioning.  The method is implemented
!  in an algorithm called SPIGMR.  Certain options in the Krylov 
!  method case are specified by INFO(13) and INFO(15).
!
!  If the Krylov method is chosen, you may supply a pair of routines,
!  JAC and PSOL, to apply preconditioning to the linear system.
!  If the system is A*x = b, the matrix is A = dG/dY + CJ*dG/dYPRIME
!  (of order NEQ).  This system can then be preconditioned in the form
!  (P-inverse)*A*x = (P-inverse)*b, with left preconditioner P.
!  (DDASPK does not allow right preconditioning.)
!  Then the Krylov method is applied to this altered, but equivalent,
!  linear system, hopefully with much better performance than without
!  preconditioning.  (In addition, a diagonal scaling matrix based on
!  the tolerances is also introduced into the altered system.)
!
!  The JAC routine evaluates any data needed for solving systems
!  with coefficient matrix P, and PSOL carries out that solution.
!  In any case, in order to improve convergence, you should try to
!  make P approximate the matrix A as much as possible, while keeping
!  the system P*x = b reasonably easy and inexpensive to solve for x,
!  given a vector b.
!
!
! *Description
!
!------INPUT - WHAT TO DO ON THE FIRST CALL TO DDASPK-------------------
!
!
!  The first call of the code is defined to be the start of each new
!  problem.  Read through the descriptions of all the following items,
!  provide sufficient storage space for designated arrays, set
!  appropriate variables for the initialization of the problem, and
!  give information about how you want the problem to be solved.
!
!
!  RES -- Provide a subroutine of the form
!
!             SUBROUTINE RES (T, Y, YPRIME, CJ, DELTA, IRES, RPAR, IPAR)
!
!         to define the system of differential/algebraic
!         equations which is to be solved. For the given values
!         of T, Y and YPRIME, the subroutine should return
!         the residual of the differential/algebraic system
!             DELTA = G(T,Y,YPRIME)
!         DELTA is a vector of length NEQ which is output from RES.
!
!         Subroutine RES must not alter T, Y, YPRIME, or CJ.
!         You must declare the name RES in an EXTERNAL
!         statement in your program that calls DDASPK.
!         You must dimension Y, YPRIME, and DELTA in RES.
!
!         The input argument CJ can be ignored, or used to rescale
!         constraint equations in the system (see Ref. 2, p. 145).
!         Note: In this respect, DDASPK is not downward-compatible
!         with DDASSL, which does not have the RES argument CJ.
!
!         IRES is an integer flag which is always equal to zero
!         on input.  Subroutine RES should alter IRES only if it
!         encounters an illegal value of Y or a stop condition.
!         Set IRES = -1 if an input value is illegal, and DDASPK
!         will try to solve the problem without getting IRES = -1.
!         If IRES = -2, DDASPK will return control to the calling
!         program with IDID = -11.
!
!         RPAR and IPAR are real and integer parameter arrays which
!         you can use for communication between your calling program
!         and subroutine RES. They are not altered by DDASPK. If you
!         do not need RPAR or IPAR, ignore these parameters by treat-
!         ing them as dummy arguments. If you do choose to use them,
!         dimension them in your calling program and in RES as arrays
!         of appropriate length.
!
!  NEQ -- Set it to the number of equations in the system (NEQ .GE. 1).
!
!  T -- Set it to the initial point of the integration. (T must be
!       a variable.)
!
!  Y(*) -- Set this array to the initial values of the NEQ solution
!          components at the initial point.  You must dimension Y of
!          length at least NEQ in your calling program.
!
!  YPRIME(*) -- Set this array to the initial values of the NEQ first
!               derivatives of the solution components at the initial
!               point.  You must dimension YPRIME at least NEQ in your
!               calling program. 
!
!  TOUT - Set it to the first point at which a solution is desired.
!         You cannot take TOUT = T.  Integration either forward in T
!         (TOUT .GT. T) or backward in T (TOUT .LT. T) is permitted.
!
!         The code advances the solution from T to TOUT using step
!         sizes which are automatically selected so as to achieve the
!         desired accuracy.  If you wish, the code will return with the
!         solution and its derivative at intermediate steps (the
!         intermediate-output mode) so that you can monitor them,
!         but you still must provide TOUT in accord with the basic
!         aim of the code.
!
!         The first step taken by the code is a critical one because
!         it must reflect how fast the solution changes near the
!         initial point.  The code automatically selects an initial
!         step size which is practically always suitable for the
!         problem.  By using the fact that the code will not step past
!         TOUT in the first step, you could, if necessary, restrict the
!         length of the initial step.
!
!         For some problems it may not be permissible to integrate
!         past a point TSTOP, because a discontinuity occurs there
!         or the solution or its derivative is not defined beyond
!         TSTOP.  When you have declared a TSTOP point (see INFO(4)
!         and RWORK(1)), you have told the code not to integrate past
!         TSTOP.  In this case any tout beyond TSTOP is invalid input.
!
!  INFO(*) - Use the INFO array to give the code more details about
!            how you want your problem solved.  This array should be
!            dimensioned of length 20, though DDASPK uses only the 
!            first 15 entries.  You must respond to all of the following
!            items, which are arranged as questions.  The simplest use
!            of DDASPK corresponds to setting all entries of INFO to 0.
!
!       INFO(1) - This parameter enables the code to initialize itself.
!              You must set it to indicate the start of every new 
!              problem.
!
!          **** Is this the first call for this problem ...
!                yes - set INFO(1) = 0
!                 no - not applicable here.
!                      See below for continuation calls.  ****
!
!       INFO(2) - How much accuracy you want of your solution
!              is specified by the error tolerances RTOL and ATOL.
!              The simplest use is to take them both to be scalars.
!              To obtain more flexibility, they can both be arrays.
!              The code must be told your choice.
!
!          **** Are both error tolerances RTOL, ATOL scalars ...
!                yes - set INFO(2) = 0
!                      and input scalars for both RTOL and ATOL
!                 no - set INFO(2) = 1
!                      and input arrays for both RTOL and ATOL ****
!
!       INFO(3) - The code integrates from T in the direction of TOUT
!              by steps.  If you wish, it will return the computed
!              solution and derivative at the next intermediate step
!              (the intermediate-output mode) or TOUT, whichever comes
!              first.  This is a good way to proceed if you want to
!              see the behavior of the solution.  If you must have
!              solutions at a great many specific TOUT points, this
!              code will compute them efficiently.
!
!          **** Do you want the solution only at
!               TOUT (and not at the next intermediate step) ...
!                yes - set INFO(3) = 0
!                 no - set INFO(3) = 1 ****
!
!       INFO(4) - To handle solutions at a great many specific
!              values TOUT efficiently, this code may integrate past
!              TOUT and interpolate to obtain the result at TOUT.
!              Sometimes it is not possible to integrate beyond some
!              point TSTOP because the equation changes there or it is
!              not defined past TSTOP.  Then you must tell the code
!              this stop condition.
!
!           **** Can the integration be carried out without any
!                restrictions on the independent variable T ...
!                 yes - set INFO(4) = 0
!                  no - set INFO(4) = 1
!                       and define the stopping point TSTOP by
!                       setting RWORK(1) = TSTOP ****
!
!       INFO(5) - used only when INFO(12) = 0 (direct methods).
!              To solve differential/algebraic systems you may wish
!              to use a matrix of partial derivatives of the
!              system of differential equations.  If you do not
!              provide a subroutine to evaluate it analytically (see
!              description of the item JAC in the call list), it will
!              be approximated by numerical differencing in this code.
!              Although it is less trouble for you to have the code
!              compute partial derivatives by numerical differencing,
!              the solution will be more reliable if you provide the
!              derivatives via JAC.  Usually numerical differencing is
!              more costly than evaluating derivatives in JAC, but
!              sometimes it is not - this depends on your problem.
!
!           **** Do you want the code to evaluate the partial deriv-
!                atives automatically by numerical differences ...
!                 yes - set INFO(5) = 0
!                  no - set INFO(5) = 1
!                       and provide subroutine JAC for evaluating the
!                       matrix of partial derivatives ****
!
!       INFO(6) - used only when INFO(12) = 0 (direct methods).
!              DDASPK will perform much better if the matrix of
!              partial derivatives, dG/dY + CJ*dG/dYPRIME (here CJ is
!              a scalar determined by DDASPK), is banded and the code
!              is told this.  In this case, the storage needed will be
!              greatly reduced, numerical differencing will be performed
!              much cheaper, and a number of important algorithms will
!              execute much faster.  The differential equation is said 
!              to have half-bandwidths ML (lower) and MU (upper) if 
!              equation i involves only unknowns Y(j) with
!                             i-ML .le. j .le. i+MU .
!              For all i=1,2,...,NEQ.  Thus, ML and MU are the widths
!              of the lower and upper parts of the band, respectively,
!              with the main diagonal being excluded.  If you do not
!              indicate that the equation has a banded matrix of partial
!              derivatives the code works with a full matrix of NEQ**2
!              elements (stored in the conventional way).  Computations
!              with banded matrices cost less time and storage than with
!              full matrices if  2*ML+MU .lt. NEQ.  If you tell the
!              code that the matrix of partial derivatives has a banded
!              structure and you want to provide subroutine JAC to
!              compute the partial derivatives, then you must be careful
!              to store the elements of the matrix in the special form
!              indicated in the description of JAC.
!
!          **** Do you want to solve the problem using a full (dense)
!               matrix (and not a special banded structure) ...
!                yes - set INFO(6) = 0
!                 no - set INFO(6) = 1
!                       and provide the lower (ML) and upper (MU)
!                       bandwidths by setting
!                       IWORK(1)=ML
!                       IWORK(2)=MU ****
!
!       INFO(7) - You can specify a maximum (absolute value of)
!              stepsize, so that the code will avoid passing over very
!              large regions.
!
!          ****  Do you want the code to decide on its own the maximum
!                stepsize ...
!                 yes - set INFO(7) = 0
!                  no - set INFO(7) = 1
!                       and define HMAX by setting
!                       RWORK(2) = HMAX ****
!
!       INFO(8) -  Differential/algebraic problems may occasionally
!              suffer from severe scaling difficulties on the first
!              step.  If you know a great deal about the scaling of 
!              your problem, you can help to alleviate this problem 
!              by specifying an initial stepsize H0.
!
!          ****  Do you want the code to define its own initial
!                stepsize ...
!                 yes - set INFO(8) = 0
!                  no - set INFO(8) = 1
!                       and define H0 by setting
!                       RWORK(3) = H0 ****
!
!       INFO(9) -  If storage is a severe problem, you can save some
!              storage by restricting the maximum method order MAXORD.
!              The default value is 5.  For each order decrease below 5,
!              the code requires NEQ fewer locations, but it is likely 
!              to be slower.  In any case, you must have 
!              1 .le. MAXORD .le. 5.
!          ****  Do you want the maximum order to default to 5 ...
!                 yes - set INFO(9) = 0
!                  no - set INFO(9) = 1
!                       and define MAXORD by setting
!                       IWORK(3) = MAXORD ****
!
!       INFO(10) - If you know that certain components of the
!              solutions to your equations are always nonnegative
!              (or nonpositive), it may help to set this
!              parameter.  There are three options that are
!              available:
!              1.  To have constraint checking only in the initial
!                  condition calculation.
!              2.  To enforce nonnegativity in Y during the integration.
!              3.  To enforce both options 1 and 2.
!
!              When selecting option 2 or 3, it is probably best to try the
!              code without using this option first, and only use
!              this option if that does not work very well.
!
!          ****  Do you want the code to solve the problem without
!                invoking any special inequality constraints ...
!                 yes - set INFO(10) = 0
!                  no - set INFO(10) = 1 to have option 1 enforced 
!                  no - set INFO(10) = 2 to have option 2 enforced
!                  no - set INFO(10) = 3 to have option 3 enforced ****
!
!                  If you have specified INFO(10) = 1 or 3, then you
!                  will also need to identify how each component of Y
!                  in the initial condition calculation is constrained.
!                  You must set:
!                  IWORK(40+I) = +1 if Y(I) must be .GE. 0,
!                  IWORK(40+I) = +2 if Y(I) must be .GT. 0,
!                  IWORK(40+I) = -1 if Y(I) must be .LE. 0, while
!                  IWORK(40+I) = -2 if Y(I) must be .LT. 0, while
!                  IWORK(40+I) =  0 if Y(I) is not constrained.
!
!       INFO(11) - DDASPK normally requires the initial T, Y, and
!              YPRIME to be consistent.  That is, you must have
!              G(T,Y,YPRIME) = 0 at the initial T.  If you do not know
!              the initial conditions precisely, in some cases
!              DDASPK may be able to compute it.
!
!              Denoting the differential variables in Y by Y_d
!              and the algebraic variables by Y_a, DDASPK can solve
!              one of two initialization problems:
!              1.  Given Y_d, calculate Y_a and Y'_d, or
!              2.  Given Y', calculate Y.
!              In either case, initial values for the given
!              components are input, and initial guesses for
!              the unknown components must also be provided as input.
!
!          ****  Are the initial T, Y, YPRIME consistent ...
!
!                 yes - set INFO(11) = 0
!                  no - set INFO(11) = 1 to calculate option 1 above,
!                    or set INFO(11) = 2 to calculate option 2 ****
!
!                  If you have specified INFO(11) = 1, then you
!                  will also need to identify  which are the
!                  differential and which are the algebraic
!                  components (algebraic components are components
!                  whose derivatives do not appear explicitly
!                  in the function G(T,Y,YPRIME)).  You must set:
!                  IWORK(LID+I) = +1 if Y(I) is a differential variable
!                  IWORK(LID+I) = -1 if Y(I) is an algebraic variable,
!                  where LID = 40 if INFO(10) = 0 or 2 and LID = 40+NEQ
!                  if INFO(10) = 1 or 3.
!
!       INFO(12) - Except for the addition of the RES argument CJ,
!              DDASPK by default is downward-compatible with DDASSL,
!              which uses only direct (dense or band) methods to solve 
!              the linear systems involved.  You must set INFO(12) to
!              indicate whether you want the direct methods or the
!              Krylov iterative method.
!          ****   Do you want DDASPK to use standard direct methods
!                 (dense or band) or the Krylov (iterative) method ...
!                   direct methods - set INFO(12) = 0.
!                   Krylov method  - set INFO(12) = 1,
!                       and check the settings of INFO(13) and INFO(15).
!
!       INFO(13) - used when INFO(12) = 1 (Krylov methods).  
!              DDASPK uses scalars MAXL, KMP, NRMAX, and EPLI for the
!              iterative solution of linear systems.  INFO(13) allows 
!              you to override the default values of these parameters.  
!              These parameters and their defaults are as follows:
!              MAXL = maximum number of iterations in the SPIGMR 
!                 algorithm (MAXL .le. NEQ).  The default is 
!                 MAXL = MIN(5,NEQ).
!              KMP = number of vectors on which orthogonalization is 
!                 done in the SPIGMR algorithm.  The default is 
!                 KMP = MAXL, which corresponds to complete GMRES 
!                 iteration, as opposed to the incomplete form.  
!              NRMAX = maximum number of restarts of the SPIGMR 
!                 algorithm per nonlinear iteration.  The default is
!                 NRMAX = 5.
!              EPLI = convergence test constant in SPIGMR algorithm.
!                 The default is EPLI = 0.05.
!              Note that the length of RWORK depends on both MAXL 
!              and KMP.  See the definition of LRW below.
!          ****   Are MAXL, KMP, and EPLI to be given their
!                 default values ...
!                  yes - set INFO(13) = 0
!                   no - set INFO(13) = 1,
!                        and set all of the following:
!                        IWORK(24) = MAXL (1 .le. MAXL .le. NEQ)
!                        IWORK(25) = KMP  (1 .le. KMP .le. MAXL)
!                        IWORK(26) = NRMAX  (NRMAX .ge. 0)
!                        RWORK(10) = EPLI (0 .lt. EPLI .lt. 1.0) ****
!
!        INFO(14) - used with INFO(11) > 0 (initial condition 
!               calculation is requested).  In this case, you may
!               request control to be returned to the calling program
!               immediately after the initial condition calculation,
!               before proceeding to the integration of the system
!               (e.g. to examine the computed Y and YPRIME).
!               If this is done, and if the initialization succeeded
!               (IDID = 4), you should reset INFO(11) to 0 for the
!               next call, to prevent the solver from repeating the 
!               initialization (and to avoid an infinite loop). 
!          ****   Do you want to proceed to the integration after
!                 the initial condition calculation is done ...
!                 yes - set INFO(14) = 0
!                  no - set INFO(14) = 1                        ****
!
!        INFO(15) - used when INFO(12) = 1 (Krylov methods).
!               When using preconditioning in the Krylov method,
!               you must supply a subroutine, PSOL, which solves the
!               associated linear systems using P.
!               The usage of DDASPK is simpler if PSOL can carry out
!               the solution without any prior calculation of data.
!               However, if some partial derivative data is to be
!               calculated in advance and used repeatedly in PSOL,
!               then you must supply a JAC routine to do this,
!               and set INFO(15) to indicate that JAC is to be called
!               for this purpose.  For example, P might be an
!               approximation to a part of the matrix A which can be
!               calculated and LU-factored for repeated solutions of
!               the preconditioner system.  The arrays WP and IWP
!               (described under JAC and PSOL) can be used to
!               communicate data between JAC and PSOL.
!          ****   Does PSOL operate with no prior preparation ...
!                 yes - set INFO(15) = 0 (no JAC routine)
!                  no - set INFO(15) = 1
!                       and supply a JAC routine to evaluate and
!                       preprocess any required Jacobian data.  ****
!
!         INFO(16) - option to exclude algebraic variables from
!               the error test.  
!          ****   Do you wish to control errors locally on
!                 all the variables...
!                 yes - set INFO(16) = 0
!                  no - set INFO(16) = 1
!                       If you have specified INFO(16) = 1, then you
!                       will also need to identify  which are the
!                       differential and which are the algebraic
!                       components (algebraic components are components
!                       whose derivatives do not appear explicitly
!                       in the function G(T,Y,YPRIME)).  You must set:
!                       IWORK(LID+I) = +1 if Y(I) is a differential 
!                                      variable, and
!                       IWORK(LID+I) = -1 if Y(I) is an algebraic
!                                      variable,
!                       where LID = 40 if INFO(10) = 0 or 2 and 
!                       LID = 40 + NEQ if INFO(10) = 1 or 3.
!
!       INFO(17) - used when INFO(11) > 0 (DDASPK is to do an 
!              initial condition calculation).
!              DDASPK uses several heuristic control quantities in the
!              initial condition calculation.  They have default values,
!              but can  also be set by the user using INFO(17).
!              These parameters and their defaults are as follows:
!              MXNIT  = maximum number of Newton iterations
!                 per Jacobian or preconditioner evaluation.
!                 The default is:
!                 MXNIT =  5 in the direct case (INFO(12) = 0), and
!                 MXNIT = 15 in the Krylov case (INFO(12) = 1).
!              MXNJ   = maximum number of Jacobian or preconditioner
!                 evaluations.  The default is:
!                 MXNJ = 6 in the direct case (INFO(12) = 0), and
!                 MXNJ = 2 in the Krylov case (INFO(12) = 1).
!              MXNH   = maximum number of values of the artificial
!                 stepsize parameter H to be tried if INFO(11) = 1.
!                 The default is MXNH = 5.
!                 NOTE: the maximum number of Newton iterations
!                 allowed in all is MXNIT*MXNJ*MXNH if INFO(11) = 1,
!                 and MXNIT*MXNJ if INFO(11) = 2.
!              LSOFF  = flag to turn off the linesearch algorithm
!                 (LSOFF = 0 means linesearch is on, LSOFF = 1 means
!                 it is turned off).  The default is LSOFF = 0.
!              STPTOL = minimum scaled step in linesearch algorithm.
!                 The default is STPTOL = (unit roundoff)**(2/3).
!              EPINIT = swing factor in the Newton iteration convergence
!                 test.  The test is applied to the residual vector,
!                 premultiplied by the approximate Jacobian (in the
!                 direct case) or the preconditioner (in the Krylov
!                 case).  For convergence, the weighted RMS norm of
!                 this vector (scaled by the error weights) must be
!                 less than EPINIT*EPCON, where EPCON = .33 is the
!                 analogous test constant used in the time steps.
!                 The default is EPINIT = .01.
!          ****   Are the initial condition heuristic controls to be 
!                 given their default values...
!                  yes - set INFO(17) = 0
!                   no - set INFO(17) = 1,
!                        and set all of the following:
!                        IWORK(32) = MXNIT (.GT. 0)
!                        IWORK(33) = MXNJ (.GT. 0)
!                        IWORK(34) = MXNH (.GT. 0)
!                        IWORK(35) = LSOFF ( = 0 or 1)
!                        RWORK(14) = STPTOL (.GT. 0.0)
!                        RWORK(15) = EPINIT (.GT. 0.0)  ****
!
!         INFO(18) - option to get extra printing in initial condition 
!                calculation.
!          ****   Do you wish to have extra printing...
!                 no  - set INFO(18) = 0
!                 yes - set INFO(18) = 1 for minimal printing, or
!                       set INFO(18) = 2 for full printing.
!                       If you have specified INFO(18) .ge. 1, data
!                       will be printed with the error handler routines.
!                       To print to a non-default unit number L, include
!                       the line  CALL XSETUN(L)  in your program.  ****
!
!   RTOL, ATOL -- You must assign relative (RTOL) and absolute (ATOL)
!               error tolerances to tell the code how accurately you
!               want the solution to be computed.  They must be defined
!               as variables because the code may change them.
!               you have two choices --
!                     Both RTOL and ATOL are scalars (INFO(2) = 0), or
!                     both RTOL and ATOL are vectors (INFO(2) = 1).
!               In either case all components must be non-negative.
!
!               The tolerances are used by the code in a local error
!               test at each step which requires roughly that
!                        abs(local error in Y(i)) .le. EWT(i) ,
!               where EWT(i) = RTOL*abs(Y(i)) + ATOL is an error weight 
!               quantity, for each vector component.
!               (More specifically, a root-mean-square norm is used to
!               measure the size of vectors, and the error test uses the
!               magnitude of the solution at the beginning of the step.)
!
!               The true (global) error is the difference between the
!               true solution of the initial value problem and the
!               computed approximation.  Practically all present day
!               codes, including this one, control the local error at
!               each step and do not even attempt to control the global
!               error directly.
!
!               Usually, but not always, the true accuracy of
!               the computed Y is comparable to the error tolerances.
!               This code will usually, but not always, deliver a more
!               accurate solution if you reduce the tolerances and
!               integrate again.  By comparing two such solutions you 
!               can get a fairly reliable idea of the true error in the
!               solution at the larger tolerances.
!
!               Setting ATOL = 0. results in a pure relative error test
!               on that component.  Setting RTOL = 0. results in a pure
!               absolute error test on that component.  A mixed test
!               with non-zero RTOL and ATOL corresponds roughly to a
!               relative error test when the solution component is
!               much bigger than ATOL and to an absolute error test
!               when the solution component is smaller than the
!               threshold ATOL.
!
!               The code will not attempt to compute a solution at an
!               accuracy unreasonable for the machine being used.  It
!               will advise you if you ask for too much accuracy and
!               inform you as to the maximum accuracy it believes
!               possible.
!
!  RWORK(*) -- a real work array, which should be dimensioned in your
!               calling program with a length equal to the value of
!               LRW (or greater).
!
!  LRW -- Set it to the declared length of the RWORK array.  The
!               minimum length depends on the options you have selected,
!               given by a base value plus additional storage as described
!               below.
!
!               If INFO(12) = 0 (standard direct method), the base value is
!               base = 50 + max(MAXORD+4,7)*NEQ.
!               The default value is MAXORD = 5 (see INFO(9)).  With the
!               default MAXORD, base = 50 + 9*NEQ.
!               Additional storage must be added to the base value for
!               any or all of the following options:
!                 if INFO(6) = 0 (dense matrix), add NEQ**2
!                 if INFO(6) = 1 (banded matrix), then
!                    if INFO(5) = 0, add (2*ML+MU+1)*NEQ + 2*(NEQ/(ML+MU+1)+1),
!                    if INFO(5) = 1, add (2*ML+MU+1)*NEQ,
!                 if INFO(16) = 1, add NEQ.
!
!              If INFO(12) = 1 (Krylov method), the base value is
!              base = 50 + (MAXORD+5)*NEQ + (MAXL+3+MIN0(1,MAXL-KMP))*NEQ +
!                      + (MAXL+3)*MAXL + 1 + LENWP.
!              See PSOL for description of LENWP.  The default values are:
!              MAXORD = 5 (see INFO(9)), MAXL = min(5,NEQ) and KMP = MAXL 
!              (see INFO(13)).
!              With the default values for MAXORD, MAXL and KMP,
!              base = 91 + 18*NEQ + LENWP.
!              Additional storage must be added to the base value for
!              any or all of the following options:
!                if INFO(16) = 1, add NEQ.
!
!
!  IWORK(*) -- an integer work array, which should be dimensioned in
!              your calling program with a length equal to the value
!              of LIW (or greater).
!
!  LIW -- Set it to the declared length of the IWORK array.  The
!             minimum length depends on the options you have selected,
!             given by a base value plus additional storage as described
!             below.
!
!             If INFO(12) = 0 (standard direct method), the base value is
!             base = 40 + NEQ.
!             IF INFO(10) = 1 or 3, add NEQ to the base value.
!             If INFO(11) = 1 or INFO(16) =1, add NEQ to the base value.
!
!             If INFO(12) = 1 (Krylov method), the base value is
!             base = 40 + LENIWP.
!             See PSOL for description of LENIWP.
!             IF INFO(10) = 1 or 3, add NEQ to the base value.
!             If INFO(11) = 1 or INFO(16) = 1, add NEQ to the base value.
!
!
!  RPAR, IPAR -- These are arrays of double precision and integer type,
!             respectively, which are available for you to use
!             for communication between your program that calls
!             DDASPK and the RES subroutine (and the JAC and PSOL
!             subroutines).  They are not altered by DDASPK.
!             If you do not need RPAR or IPAR, ignore these
!             parameters by treating them as dummy arguments.
!             If you do choose to use them, dimension them in
!             your calling program and in RES (and in JAC and PSOL)
!             as arrays of appropriate length.
!
!  JAC -- This is the name of a routine that you may supply
!         (optionally) that relates to the Jacobian matrix of the
!         nonlinear system that the code must solve at each T step.
!         The role of JAC (and its call sequence) depends on whether
!         a direct (INFO(12) = 0) or Krylov (INFO(12) = 1) method 
!         is selected.
!
!         **** INFO(12) = 0 (direct methods):
!           If you are letting the code generate partial derivatives
!           numerically (INFO(5) = 0), then JAC can be absent
!           (or perhaps a dummy routine to satisfy the loader).
!           Otherwise you must supply a JAC routine to compute
!           the matrix A = dG/dY + CJ*dG/dYPRIME.  It must have
!           the form
!
!           SUBROUTINE JAC (T, Y, YPRIME, PD, CJ, RPAR, IPAR)
!
!           The JAC routine must dimension Y, YPRIME, and PD (and RPAR
!           and IPAR if used).  CJ is a scalar which is input to JAC.
!           For the given values of T, Y, and YPRIME, the JAC routine
!           must evaluate the nonzero elements of the matrix A, and 
!           store these values in the array PD.  The elements of PD are 
!           set to zero before each call to JAC, so that only nonzero
!           elements need to be defined.
!           The way you store the elements into the PD array depends
!           on the structure of the matrix indicated by INFO(6).
!           *** INFO(6) = 0 (full or dense matrix) ***
!               Give PD a first dimension of NEQ.  When you evaluate the
!               nonzero partial derivatives of equation i (i.e. of G(i))
!               with respect to component j (of Y and YPRIME), you must
!               store the element in PD according to
!                  PD(i,j) = dG(i)/dY(j) + CJ*dG(i)/dYPRIME(j).
!           *** INFO(6) = 1 (banded matrix with half-bandwidths ML, MU
!                            as described under INFO(6)) ***
!               Give PD a first dimension of 2*ML+MU+1.  When you 
!               evaluate the nonzero partial derivatives of equation i 
!               (i.e. of G(i)) with respect to component j (of Y and 
!               YPRIME), you must store the element in PD according to 
!                  IROW = i - j + ML + MU + 1
!                  PD(IROW,j) = dG(i)/dY(j) + CJ*dG(i)/dYPRIME(j).
!
!          **** INFO(12) = 1 (Krylov method):
!            If you are not calculating Jacobian data in advance for use
!            in PSOL (INFO(15) = 0), JAC can be absent (or perhaps a
!            dummy routine to satisfy the loader).  Otherwise, you may
!            supply a JAC routine to compute and preprocess any parts of
!            of the Jacobian matrix  A = dG/dY + CJ*dG/dYPRIME that are
!            involved in the preconditioner matrix P.
!            It is to have the form
!
!            SUBROUTINE JAC (RES, IRES, NEQ, T, Y, YPRIME, REWT, SAVR,
!                            WK, H, CJ, WP, IWP, IER, RPAR, IPAR)
!
!           The JAC routine must dimension Y, YPRIME, REWT, SAVR, WK,
!           and (if used) WP, IWP, RPAR, and IPAR.
!           The Y, YPRIME, and SAVR arrays contain the current values
!           of Y, YPRIME, and the residual G, respectively.  
!           The array WK is work space of length NEQ.  
!           H is the step size.  CJ is a scalar, input to JAC, that is
!           normally proportional to 1/H.  REWT is an array of 
!           reciprocal error weights, 1/EWT(i), where EWT(i) is
!           RTOL*abs(Y(i)) + ATOL (unless you supplied routine DDAWTS
!           instead), for use in JAC if needed.  For example, if JAC
!           computes difference quotient approximations to partial
!           derivatives, the REWT array may be useful in setting the
!           increments used.  The JAC routine should do any
!           factorization operations called for, in preparation for
!           solving linear systems in PSOL.  The matrix P should
!           be an approximation to the Jacobian,
!           A = dG/dY + CJ*dG/dYPRIME.
!
!           WP and IWP are real and integer work arrays which you may
!           use for communication between your JAC routine and your
!           PSOL routine.  These may be used to store elements of the 
!           preconditioner P, or related matrix data (such as factored
!           forms).  They are not altered by DDASPK.
!           If you do not need WP or IWP, ignore these parameters by
!           treating them as dummy arguments.  If you do use them,
!           dimension them appropriately in your JAC and PSOL routines.
!           See the PSOL description for instructions on setting 
!           the lengths of WP and IWP.
!
!           On return, JAC should set the error flag IER as follows..
!             IER = 0    if JAC was successful,
!             IER .ne. 0 if JAC was unsuccessful (e.g. if Y or YPRIME
!                        was illegal, or a singular matrix is found).
!           (If IER .ne. 0, a smaller stepsize will be tried.)
!           IER = 0 on entry to JAC, so need be reset only on a failure.
!           If RES is used within JAC, then a nonzero value of IRES will
!           override any nonzero value of IER (see the RES description).
!
!         Regardless of the method type, subroutine JAC must not
!         alter T, Y(*), YPRIME(*), H, CJ, or REWT(*).
!         You must declare the name JAC in an EXTERNAL statement in
!         your program that calls DDASPK.
!
! PSOL --  This is the name of a routine you must supply if you have
!         selected a Krylov method (INFO(12) = 1) with preconditioning.
!         In the direct case (INFO(12) = 0), PSOL can be absent 
!         (a dummy routine may have to be supplied to satisfy the 
!         loader).  Otherwise, you must provide a PSOL routine to 
!         solve linear systems arising from preconditioning.
!         When supplied with INFO(12) = 1, the PSOL routine is to 
!         have the form
!
!         SUBROUTINE PSOL (NEQ, T, Y, YPRIME, SAVR, WK, CJ, WGHT,
!                          WP, IWP, B, EPLIN, IER, RPAR, IPAR)
!
!         The PSOL routine must solve linear systems of the form 
!         P*x = b where P is the left preconditioner matrix.
!
!         The right-hand side vector b is in the B array on input, and
!         PSOL must return the solution vector x in B.
!         The Y, YPRIME, and SAVR arrays contain the current values
!         of Y, YPRIME, and the residual G, respectively.  
!
!         Work space required by JAC and/or PSOL, and space for data to
!         be communicated from JAC to PSOL is made available in the form
!         of arrays WP and IWP, which are parts of the RWORK and IWORK
!         arrays, respectively.  The lengths of these real and integer
!         work spaces WP and IWP must be supplied in LENWP and LENIWP,
!         respectively, as follows..
!           IWORK(27) = LENWP = length of real work space WP
!           IWORK(28) = LENIWP = length of integer work space IWP.
!
!         WK is a work array of length NEQ for use by PSOL.
!         CJ is a scalar, input to PSOL, that is normally proportional
!         to 1/H (H = stepsize).  If the old value of CJ
!         (at the time of the last JAC call) is needed, it must have
!         been saved by JAC in WP.
!
!         WGHT is an array of weights, to be used if PSOL uses an
!         iterative method and performs a convergence test.  (In terms
!         of the argument REWT to JAC, WGHT is REWT/sqrt(NEQ).)
!         If PSOL uses an iterative method, it should use EPLIN
!         (a heuristic parameter) as the bound on the weighted norm of
!         the residual for the computed solution.  Specifically, the
!         residual vector R should satisfy
!              SQRT (SUM ( (R(i)*WGHT(i))**2 ) ) .le. EPLIN
!
!         PSOL must not alter NEQ, T, Y, YPRIME, SAVR, CJ, WGHT, EPLIN.
!
!         On return, PSOL should set the error flag IER as follows..
!           IER = 0 if PSOL was successful,
!           IER .lt. 0 if an unrecoverable error occurred, meaning
!                 control will be passed to the calling routine,
!           IER .gt. 0 if a recoverable error occurred, meaning that
!                 the step will be retried with the same step size
!                 but with a call to JAC to update necessary data,
!                 unless the Jacobian data is current, in which case
!                 the step will be retried with a smaller step size.
!           IER = 0 on entry to PSOL so need be reset only on a failure.
!
!         You must declare the name PSOL in an EXTERNAL statement in
!         your program that calls DDASPK.
!
!
!  OPTIONALLY REPLACEABLE SUBROUTINE:
!
!  DDASPK uses a weighted root-mean-square norm to measure the 
!  size of various error vectors.  The weights used in this norm
!  are set in the following subroutine:
!
!    SUBROUTINE DDAWTS (NEQ, IWT, RTOL, ATOL, Y, EWT, RPAR, IPAR)
!    DIMENSION RTOL(*), ATOL(*), Y(*), EWT(*), RPAR(*), IPAR(*)
!
!  A DDAWTS routine has been included with DDASPK which sets the
!  weights according to
!    EWT(I) = RTOL*ABS(Y(I)) + ATOL
!  in the case of scalar tolerances (IWT = 0) or
!    EWT(I) = RTOL(I)*ABS(Y(I)) + ATOL(I)
!  in the case of array tolerances (IWT = 1).  (IWT is INFO(2).)
!  In some special cases, it may be appropriate for you to define
!  your own error weights by writing a subroutine DDAWTS to be 
!  called instead of the version supplied.  However, this should 
!  be attempted only after careful thought and consideration. 
!  If you supply this routine, you may use the tolerances and Y 
!  as appropriate, but do not overwrite these variables.  You
!  may also use RPAR and IPAR to communicate data as appropriate.
!  ***Note: Aside from the values of the weights, the choice of 
!  norm used in DDASPK (weighted root-mean-square) is not subject
!  to replacement by the user.  In this respect, DDASPK is not
!  downward-compatible with the original DDASSL solver (in which
!  the norm routine was optionally user-replaceable).
!
!
!------OUTPUT - AFTER ANY RETURN FROM DDASPK----------------------------
!
!  The principal aim of the code is to return a computed solution at
!  T = TOUT, although it is also possible to obtain intermediate
!  results along the way.  To find out whether the code achieved its
!  goal or if the integration process was interrupted before the task
!  was completed, you must check the IDID parameter.
!
!
!   T -- The output value of T is the point to which the solution
!        was successfully advanced.
!
!   Y(*) -- contains the computed solution approximation at T.
!
!   YPRIME(*) -- contains the computed derivative approximation at T.
!
!   IDID -- reports what the code did, described as follows:
!
!                     *** TASK COMPLETED ***
!                Reported by positive values of IDID
!
!           IDID = 1 -- a step was successfully taken in the
!                   intermediate-output mode.  The code has not
!                   yet reached TOUT.
!
!           IDID = 2 -- the integration to TSTOP was successfully
!                   completed (T = TSTOP) by stepping exactly to TSTOP.
!
!           IDID = 3 -- the integration to TOUT was successfully
!                   completed (T = TOUT) by stepping past TOUT.
!                   Y(*) and YPRIME(*) are obtained by interpolation.
!
!           IDID = 4 -- the initial condition calculation, with
!                   INFO(11) > 0, was successful, and INFO(14) = 1.
!                   No integration steps were taken, and the solution
!                   is not considered to have been started.
!
!                    *** TASK INTERRUPTED ***
!                Reported by negative values of IDID
!
!           IDID = -1 -- a large amount of work has been expended
!                     (about 500 steps).
!
!           IDID = -2 -- the error tolerances are too stringent.
!
!           IDID = -3 -- the local error test cannot be satisfied
!                     because you specified a zero component in ATOL
!                     and the corresponding computed solution component
!                     is zero.  Thus, a pure relative error test is
!                     impossible for this component.
!
!           IDID = -5 -- there were repeated failures in the evaluation
!                     or processing of the preconditioner (in JAC).
!
!           IDID = -6 -- DDASPK had repeated error test failures on the
!                     last attempted step.
!
!           IDID = -7 -- the nonlinear system solver in the time integration
!                     could not converge.
!
!           IDID = -8 -- the matrix of partial derivatives appears
!                     to be singular (direct method).
!
!           IDID = -9 -- the nonlinear system solver in the time integration
!                     failed to achieve convergence, and there were repeated 
!                     error test failures in this step.
!
!           IDID =-10 -- the nonlinear system solver in the time integration 
!                     failed to achieve convergence because IRES was equal 
!                     to -1.
!
!           IDID =-11 -- IRES = -2 was encountered and control is
!                     being returned to the calling program.
!
!           IDID =-12 -- DDASPK failed to compute the initial Y, YPRIME.
!
!           IDID =-13 -- unrecoverable error encountered inside user's
!                     PSOL routine, and control is being returned to
!                     the calling program.
!
!           IDID =-14 -- the Krylov linear system solver could not 
!                     achieve convergence.
!
!           IDID =-15,..,-32 -- Not applicable for this code.
!
!                    *** TASK TERMINATED ***
!                reported by the value of IDID=-33
!
!           IDID = -33 -- the code has encountered trouble from which
!                   it cannot recover.  A message is printed
!                   explaining the trouble and control is returned
!                   to the calling program.  For example, this occurs
!                   when invalid input is detected.
!
!   RTOL, ATOL -- these quantities remain unchanged except when
!               IDID = -2.  In this case, the error tolerances have been
!               increased by the code to values which are estimated to
!               be appropriate for continuing the integration.  However,
!               the reported solution at T was obtained using the input
!               values of RTOL and ATOL.
!
!   RWORK, IWORK -- contain information which is usually of no interest
!               to the user but necessary for subsequent calls. 
!               However, you may be interested in the performance data
!               listed below.  These quantities are accessed in RWORK 
!               and IWORK but have internal mnemonic names, as follows..
!
!               RWORK(3)--contains H, the step size h to be attempted
!                        on the next step.
!
!               RWORK(4)--contains TN, the current value of the
!                        independent variable, i.e. the farthest point
!                        integration has reached.  This will differ 
!                        from T if interpolation has been performed 
!                        (IDID = 3).
!
!               RWORK(7)--contains HOLD, the stepsize used on the last
!                        successful step.  If INFO(11) = INFO(14) = 1,
!                        this contains the value of H used in the
!                        initial condition calculation.
!
!               IWORK(7)--contains K, the order of the method to be 
!                        attempted on the next step.
!
!               IWORK(8)--contains KOLD, the order of the method used
!                        on the last step.
!
!               IWORK(11)--contains NST, the number of steps (in T) 
!                        taken so far.
!
!               IWORK(12)--contains NRE, the number of calls to RES 
!                        so far.
!
!               IWORK(13)--contains NJE, the number of calls to JAC so
!                        far (Jacobian or preconditioner evaluations).
!
!               IWORK(14)--contains NETF, the total number of error test
!                        failures so far.
!
!               IWORK(15)--contains NCFN, the total number of nonlinear
!                        convergence failures so far (includes counts
!                        of singular iteration matrix or singular
!                        preconditioners).
!
!               IWORK(16)--contains NCFL, the number of convergence
!                        failures of the linear iteration so far.
!
!               IWORK(17)--contains LENIW, the length of IWORK actually
!                        required.  This is defined on normal returns 
!                        and on an illegal input return for
!                        insufficient storage.
!
!               IWORK(18)--contains LENRW, the length of RWORK actually
!                        required.  This is defined on normal returns 
!                        and on an illegal input return for
!                        insufficient storage.
!
!               IWORK(19)--contains NNI, the total number of nonlinear
!                        iterations so far (each of which calls a
!                        linear solver).
!
!               IWORK(20)--contains NLI, the total number of linear
!                        (Krylov) iterations so far.
!
!               IWORK(21)--contains NPS, the number of PSOL calls so
!                        far, for preconditioning solve operations or
!                        for solutions with the user-supplied method.
!
!               Note: The various counters in IWORK do not include 
!               counts during a call made with INFO(11) > 0 and
!               INFO(14) = 1.
!
!
!------INPUT - WHAT TO DO TO CONTINUE THE INTEGRATION  -----------------
!              (CALLS AFTER THE FIRST)
!
!     This code is organized so that subsequent calls to continue the
!     integration involve little (if any) additional effort on your
!     part.  You must monitor the IDID parameter in order to determine
!     what to do next.
!
!     Recalling that the principal task of the code is to integrate
!     from T to TOUT (the interval mode), usually all you will need
!     to do is specify a new TOUT upon reaching the current TOUT.
!
!     Do not alter any quantity not specifically permitted below.  In
!     particular do not alter NEQ, T, Y(*), YPRIME(*), RWORK(*), 
!     IWORK(*), or the differential equation in subroutine RES.  Any 
!     such alteration constitutes a new problem and must be treated 
!     as such, i.e. you must start afresh.
!
!     You cannot change from array to scalar error control or vice
!     versa (INFO(2)), but you can change the size of the entries of
!     RTOL or ATOL.  Increasing a tolerance makes the equation easier
!     to integrate.  Decreasing a tolerance will make the equation
!     harder to integrate and should generally be avoided.
!
!     You can switch from the intermediate-output mode to the
!     interval mode (INFO(3)) or vice versa at any time.
!
!     If it has been necessary to prevent the integration from going
!     past a point TSTOP (INFO(4), RWORK(1)), keep in mind that the
!     code will not integrate to any TOUT beyond the currently
!     specified TSTOP.  Once TSTOP has been reached, you must change
!     the value of TSTOP or set INFO(4) = 0.  You may change INFO(4)
!     or TSTOP at any time but you must supply the value of TSTOP in
!     RWORK(1) whenever you set INFO(4) = 1.
!
!     Do not change INFO(5), INFO(6), INFO(12-17) or their associated
!     IWORK/RWORK locations unless you are going to restart the code.
!
!                    *** FOLLOWING A COMPLETED TASK ***
!
!     If..
!     IDID = 1, call the code again to continue the integration
!                  another step in the direction of TOUT.
!
!     IDID = 2 or 3, define a new TOUT and call the code again.
!                  TOUT must be different from T.  You cannot change
!                  the direction of integration without restarting.
!
!     IDID = 4, reset INFO(11) = 0 and call the code again to begin
!                  the integration.  (If you leave INFO(11) > 0 and
!                  INFO(14) = 1, you may generate an infinite loop.)
!                  In this situation, the next call to DASPK is 
!                  considered to be the first call for the problem,
!                  in that all initializations are done.
!
!                    *** FOLLOWING AN INTERRUPTED TASK ***
!
!     To show the code that you realize the task was interrupted and
!     that you want to continue, you must take appropriate action and
!     set INFO(1) = 1.
!
!     If..
!     IDID = -1, the code has taken about 500 steps.  If you want to
!                  continue, set INFO(1) = 1 and call the code again.
!                  An additional 500 steps will be allowed.
!
!
!     IDID = -2, the error tolerances RTOL, ATOL have been increased
!                  to values the code estimates appropriate for
!                  continuing.  You may want to change them yourself.
!                  If you are sure you want to continue with relaxed
!                  error tolerances, set INFO(1) = 1 and call the code
!                  again.
!
!     IDID = -3, a solution component is zero and you set the
!                  corresponding component of ATOL to zero.  If you
!                  are sure you want to continue, you must first alter
!                  the error criterion to use positive values of ATOL 
!                  for those components corresponding to zero solution
!                  components, then set INFO(1) = 1 and call the code
!                  again.
!
!     IDID = -4  --- cannot occur with this code.
!
!     IDID = -5, your JAC routine failed with the Krylov method.  Check
!                  for errors in JAC and restart the integration.
!
!     IDID = -6, repeated error test failures occurred on the last
!                  attempted step in DDASPK.  A singularity in the
!                  solution may be present.  If you are absolutely
!                  certain you want to continue, you should restart
!                  the integration.  (Provide initial values of Y and
!                  YPRIME which are consistent.)
!
!     IDID = -7, repeated convergence test failures occurred on the last
!                  attempted step in DDASPK.  An inaccurate or ill-
!                  conditioned Jacobian or preconditioner may be the
!                  problem.  If you are absolutely certain you want
!                  to continue, you should restart the integration.
!
!
!     IDID = -8, the matrix of partial derivatives is singular, with
!                  the use of direct methods.  Some of your equations
!                  may be redundant.  DDASPK cannot solve the problem
!                  as stated.  It is possible that the redundant
!                  equations could be removed, and then DDASPK could
!                  solve the problem.  It is also possible that a
!                  solution to your problem either does not exist
!                  or is not unique.
!
!     IDID = -9, DDASPK had multiple convergence test failures, preceded
!                  by multiple error test failures, on the last
!                  attempted step.  It is possible that your problem is
!                  ill-posed and cannot be solved using this code.  Or,
!                  there may be a discontinuity or a singularity in the
!                  solution.  If you are absolutely certain you want to
!                  continue, you should restart the integration.
!
!     IDID = -10, DDASPK had multiple convergence test failures
!                  because IRES was equal to -1.  If you are
!                  absolutely certain you want to continue, you
!                  should restart the integration.
!
!     IDID = -11, there was an unrecoverable error (IRES = -2) from RES
!                  inside the nonlinear system solver.  Determine the
!                  cause before trying again.
!
!     IDID = -12, DDASPK failed to compute the initial Y and YPRIME
!                  vectors.  This could happen because the initial 
!                  approximation to Y or YPRIME was not very good, or
!                  because no consistent values of these vectors exist.
!                  The problem could also be caused by an inaccurate or
!                  singular iteration matrix, or a poor preconditioner.
!
!     IDID = -13, there was an unrecoverable error encountered inside 
!                  your PSOL routine.  Determine the cause before 
!                  trying again.
!
!     IDID = -14, the Krylov linear system solver failed to achieve
!                  convergence.  This may be due to ill-conditioning
!                  in the iteration matrix, or a singularity in the
!                  preconditioner (if one is being used).
!                  Another possibility is that there is a better
!                  choice of Krylov parameters (see INFO(13)).
!                  Possibly the failure is caused by redundant equations
!                  in the system, or by inconsistent equations.
!                  In that case, reformulate the system to make it
!                  consistent and non-redundant.
!
!     IDID = -15,..,-32 --- Cannot occur with this code.
!
!                       *** FOLLOWING A TERMINATED TASK ***
!
!     If IDID = -33, you cannot continue the solution of this problem.
!                  An attempt to do so will result in your run being
!                  terminated.
!
!  ---------------------------------------------------------------------
!
!***REFERENCES
!  1.  L. R. Petzold, A Description of DASSL: A Differential/Algebraic
!      System Solver, in Scientific Computing, R. S. Stepleman et al.
!      (Eds.), North-Holland, Amsterdam, 1983, pp. 65-68.
!  2.  K. E. Brenan, S. L. Campbell, and L. R. Petzold, Numerical 
!      Solution of Initial-Value Problems in Differential-Algebraic
!      Equations, Elsevier, New York, 1989.
!  3.  P. N. Brown and A. C. Hindmarsh, Reduced Storage Matrix Methods
!      in Stiff ODE Systems, J. Applied Mathematics and Computation,
!      31 (1989), pp. 40-91.
!  4.  P. N. Brown, A. C. Hindmarsh, and L. R. Petzold, Using Krylov
!      Methods in the Solution of Large-Scale Differential-Algebraic
!      Systems, SIAM J. Sci. Comp., 15 (1994), pp. 1467-1488.
!  5.  P. N. Brown, A. C. Hindmarsh, and L. R. Petzold, Consistent
!      Initial Condition Calculation for Differential-Algebraic
!      Systems, LLNL Report UCRL-JC-122175, August 1995; submitted to
!      SIAM J. Sci. Comp.
!
!***ROUTINES CALLED
!
!   The following are all the subordinate routines used by DDASPK.
!
!   DDASIC computes consistent initial conditions.
!   DYYPNW updates Y and YPRIME in linesearch for initial condition
!          calculation.
!   DDSTP  carries out one step of the integration.
!   DCNSTR/DCNST0 check the current solution for constraint violations.
!   DDAWTS sets error weight quantities.
!   DINVWT tests and inverts the error weights.
!   DDATRP performs interpolation to get an output solution.
!   DDWNRM computes the weighted root-mean-square norm of a vector.
!   D1MACH provides the unit roundoff of the computer.
!   XERRWD/XSETF/XSETUN/IXSAV is a package to handle error messages. 
!   DDASID nonlinear equation driver to initialize Y and YPRIME using
!          direct linear system solver methods.  Interfaces to Newton
!          solver (direct case).
!   DNSID  solves the nonlinear system for unknown initial values by
!          modified Newton iteration and direct linear system methods.
!   DLINSD carries out linesearch algorithm for initial condition
!          calculation (direct case).
!   DFNRMD calculates weighted norm of preconditioned residual in
!          initial condition calculation (direct case).
!   DNEDD  nonlinear equation driver for direct linear system solver
!          methods.  Interfaces to Newton solver (direct case).
!   DMATD  assembles the iteration matrix (direct case).
!   DNSD   solves the associated nonlinear system by modified
!          Newton iteration and direct linear system methods.
!   DSLVD  interfaces to linear system solver (direct case).
!   DDASIK nonlinear equation driver to initialize Y and YPRIME using
!          Krylov iterative linear system methods.  Interfaces to
!          Newton solver (Krylov case).
!   DNSIK  solves the nonlinear system for unknown initial values by
!          Newton iteration and Krylov iterative linear system methods.
!   DLINSK carries out linesearch algorithm for initial condition
!          calculation (Krylov case).
!   DFNRMK calculates weighted norm of preconditioned residual in
!          initial condition calculation (Krylov case).
!   DNEDK  nonlinear equation driver for iterative linear system solver
!          methods.  Interfaces to Newton solver (Krylov case).
!   DNSK   solves the associated nonlinear system by Inexact Newton
!          iteration and (linear) Krylov iteration.
!   DSLVK  interfaces to linear system solver (Krylov case).
!   DSPIGM solves a linear system by SPIGMR algorithm.
!   DATV   computes matrix-vector product in Krylov algorithm.
!   DORTH  performs orthogonalization of Krylov basis vectors.
!   DHEQR  performs QR factorization of Hessenberg matrix.
!   DHELS  finds least-squares solution of Hessenberg linear system.
!   DGEFA, DGESL, DGBFA, DGBSL are LINPACK routines for solving 
!          linear systems (dense or band direct methods).
!   DAXPY, DCOPY, DDOT, DNRM2, DSCAL are Basic Linear Algebra (BLAS)
!          routines.
!
! The routines called directly by DDASPK are:
!   DCNST0, DDAWTS, DINVWT, D1MACH, DDWNRM, DDASIC, DDATRP, DDSTP,
!   XERRWD
!
!***END PROLOGUE DDASPK
!
!
IMPLICIT DOUBLE PRECISION(A-H,O-Z)
LOGICAL DONE, LAVL, LCFN, LCFL, LWARN
DIMENSION Y(*),YPRIME(*)
DIMENSION INFO(25)        ! Karline: increased from 20 -> 25
INTEGER NIND(3)           !          added
DIMENSION RWORK(LRW),IWORK(LIW)
DIMENSION RTOL(*),ATOL(*)
DIMENSION RPAR(*),IPAR(*)
EXTERNAL  RES, JAC, PSOL, DDASID, DDASIK, DNEDD, DNEDK
!
!     Set pointers into IWORK.
!
INTEGER, PARAMETER :: LML=1, LMU=2, LMTYPE=4, &
& LIWM=1, LMXORD=3, LJCALC=5, LPHASE=6, LK=7, LKOLD=8, &
& LNS=9, LNSTL=10, LNST=11, LNRE=12, LNJE=13, LETF=14, LNCFN=15, &
& LNCFL=16, LNIW=17, LNRW=18, LNNI=19, LNLI=20, LNPS=21, &
& LNPD=22, LMITER=23, LMAXL=24, LKMP=25, LNRMAX=26, LLNWP=27, &
& LLNIWP=28, LLOCWP=29, LLCIWP=30, LKPRIN=31, &
& LMXNIT=32, LMXNJ=33, LMXNH=34, LLSOFF=35, LICNS=41
!
!     Set pointers into RWORK.
!
INTEGER, PARAMETER :: LTSTOP=1, LHMAX=2, LH=3, LTN=4, LCJ=5, LCJOLD=6, &
& LHOLD=7, LS=8, LROUND=9, LEPLI=10, LSQRN=11, LRSQRN=12, &
& LEPCON=13, LSTOL=14, LEPIN=15, &
& LALPHA=21, LBETA=27, LGAMMA=33, LPSI=39, LSIGMA=45, LDELTA=51
!
SAVE LID, LENID, NONNEG
!
!
!***FIRST EXECUTABLE STATEMENT  DDASPK
!
!
! Karline: the index of each variable
DO I = 1, 3
  NIND(I) = INFO(20+I)
ENDDO
! Karline: initialised some variables to avoid compiler warnings - should have no effect
LENRW = 0
LENIW = 0
LENPD = 0
LENWP = 0
! Karline      
IF(INFO(1).NE.0) GO TO 100
!
!-----------------------------------------------------------------------
!     This block is executed for the initial call only.
!     It contains checking of inputs and initializations.
!-----------------------------------------------------------------------
!
!     First check INFO array to make sure all elements of INFO
!     Are within the proper range.  (INFO(1) is checked later, because
!     it must be tested on every call.) ITEMP holds the location
!     within INFO which may be out of range.
!
DO 10 I=2,9
   ITEMP = I
   IF (INFO(I) .NE. 0 .AND. INFO(I) .NE. 1) GO TO 701
10 CONTINUE
ITEMP = 10
IF(INFO(10).LT.0 .OR. INFO(10).GT.3) GO TO 701
ITEMP = 11
IF(INFO(11).LT.0 .OR. INFO(11).GT.2) GO TO 701
DO 15 I=12,17
   ITEMP = I
   IF (INFO(I) .NE. 0 .AND. INFO(I) .NE. 1) GO TO 701
15 CONTINUE
ITEMP = 18
IF(INFO(18).LT.0 .OR. INFO(18).GT.2) GO TO 701

!
!     Check NEQ to see if it is positive.
!
IF (NEQ .LE. 0) GO TO 702
!
!     Check and compute maximum order.
!
MXORD=5
IF (INFO(9) .NE. 0) THEN
   MXORD=IWORK(LMXORD)
   IF (MXORD .LT. 1 .OR. MXORD .GT. 5) GO TO 703
   ENDIF
IWORK(LMXORD)=MXORD
!
!     Set and/or check inputs for constraint checking (INFO(10) .NE. 0).
!     Set values for ICNFLG, NONNEG, and pointer LID.
!
ICNFLG = 0
NONNEG = 0
LID = LICNS
IF (INFO(10) .EQ. 0) GO TO 20
IF (INFO(10) .EQ. 1) THEN
   ICNFLG = 1
   NONNEG = 0
   LID = LICNS + NEQ
ELSEIF (INFO(10) .EQ. 2) THEN
   ICNFLG = 0
   NONNEG = 1
ELSE
   ICNFLG = 1
   NONNEG = 1
   LID = LICNS + NEQ
ENDIF
!
20 CONTINUE
!
!     Set and/or check inputs for Krylov solver (INFO(12) .NE. 0).
!     If indicated, set default values for MAXL, KMP, NRMAX, and EPLI.
!     Otherwise, verify inputs required for iterative solver.
!
IF (INFO(12) .EQ. 0) GO TO 25
!
IWORK(LMITER) = INFO(12)
IF (INFO(13) .EQ. 0) THEN
   IWORK(LMAXL) = MIN(5,NEQ)
   IWORK(LKMP) = IWORK(LMAXL)
   IWORK(LNRMAX) = 5
   RWORK(LEPLI) = 0.05D0
ELSE
   IF(IWORK(LMAXL) .LT. 1 .OR. IWORK(LMAXL) .GT. NEQ) GO TO 720
   IF(IWORK(LKMP) .LT. 1 .OR. IWORK(LKMP) .GT. IWORK(LMAXL)) &
& GO TO 721
   IF(IWORK(LNRMAX) .LT. 0) GO TO 722
   IF(RWORK(LEPLI).LE.0.0D0 .OR. RWORK(LEPLI).GE.1.0D0)GO TO 723
   ENDIF
!
25 CONTINUE
!
!     Set and/or check controls for the initial condition calculation
!     (INFO(11) .GT. 0).  If indicated, set default values.
!     Otherwise, verify inputs required for iterative solver.
!
IF (INFO(11) .EQ. 0) GO TO 30
IF (INFO(17) .EQ. 0) THEN
  IWORK(LMXNIT) = 5
  IF (INFO(12) .GT. 0) IWORK(LMXNIT) = 15
  IWORK(LMXNJ) = 6
  IF (INFO(12) .GT. 0) IWORK(LMXNJ) = 2
  IWORK(LMXNH) = 5
  IWORK(LLSOFF) = 0
  RWORK(LEPIN) = 0.01D0
ELSE
  IF (IWORK(LMXNIT) .LE. 0) GO TO 725
  IF (IWORK(LMXNJ) .LE. 0) GO TO 725
  IF (IWORK(LMXNH) .LE. 0) GO TO 725
  LSOFF = IWORK(LLSOFF)
  IF (LSOFF .LT. 0 .OR. LSOFF .GT. 1) GO TO 725
  IF (RWORK(LEPIN) .LE. 0.0D0) GO TO 725
  ENDIF
!
30 CONTINUE
!
!     Below is the computation and checking of the work array lengths
!     LENIW and LENRW, using direct methods (INFO(12) = 0) or
!     the Krylov methods (INFO(12) = 1).
!
LENIC = 0
IF (INFO(10) .EQ. 1 .OR. INFO(10) .EQ. 3) LENIC = NEQ
LENID = 0
IF (INFO(11) .EQ. 1 .OR. INFO(16) .EQ. 1) LENID = NEQ
IF (INFO(12) .EQ. 0) THEN
!
!        Compute MTYPE, etc.  Check ML and MU.
!
   NCPHI = MAX(MXORD + 1, 4)
   IF(INFO(6).EQ.0) THEN
      LENPD = NEQ**2
      LENRW = 50 + (NCPHI+3)*NEQ + LENPD
      IF(INFO(5).EQ.0) THEN
         IWORK(LMTYPE)=2
      ELSE
         IWORK(LMTYPE)=1
      ENDIF
   ELSE
      IF(IWORK(LML).LT.0.OR.IWORK(LML).GE.NEQ)GO TO 717
      IF(IWORK(LMU).LT.0.OR.IWORK(LMU).GE.NEQ)GO TO 718
      LENPD=(2*IWORK(LML)+IWORK(LMU)+1)*NEQ
      IF(INFO(5).EQ.0) THEN
         IWORK(LMTYPE)=5
         MBAND=IWORK(LML)+IWORK(LMU)+1
         MSAVE=(NEQ/MBAND)+1
         LENRW = 50 + (NCPHI+3)*NEQ + LENPD + 2*MSAVE
      ELSE
         IWORK(LMTYPE)=4
         LENRW = 50 + (NCPHI+3)*NEQ + LENPD
      ENDIF
   ENDIF
!
!        Compute LENIW, LENWP, LENIWP.
!
   LENIW = 40 + LENIC + LENID + NEQ
   LENWP = 0
   LENIWP = 0
!
ELSE IF (INFO(12) .EQ. 1)  THEN
   MAXL = IWORK(LMAXL)
   LENWP = IWORK(LLNWP)
   LENIWP = IWORK(LLNIWP)
   LENPD = (MAXL+3+MIN0(1,MAXL-IWORK(LKMP)))*NEQ &
& + (MAXL+3)*MAXL + 1 + LENWP
   LENRW = 50 + (IWORK(LMXORD)+5)*NEQ + LENPD
   LENIW = 40 + LENIC + LENID + LENIWP
!
ENDIF
IF(INFO(16) .NE. 0) LENRW = LENRW + NEQ
!
!     Check lengths of RWORK and IWORK.
!
IWORK(LNIW)=LENIW
IWORK(LNRW)=LENRW
IWORK(LNPD)=LENPD
IWORK(LLOCWP) = LENPD-LENWP+1
IF(LRW.LT.LENRW)GO TO 704
IF(LIW.LT.LENIW)GO TO 705
!
!     Check ICNSTR for legality.
!
IF (LENIC .GT. 0) THEN
  DO 40 I = 1,NEQ
    ICI = IWORK(LICNS-1+I)
    IF (ICI .LT. -2 .OR. ICI .GT. 2) GO TO 726
40 CONTINUE
  ENDIF
!
!     Check Y for consistency with constraints.
!
IF (LENIC .GT. 0) THEN
  CALL DCNST0(NEQ,Y,IWORK(LICNS),IRET)
  IF (IRET .NE. 0) GO TO 727
  ENDIF
!
!     Check ID for legality.
!
IF (LENID .GT. 0) THEN
  DO 50 I = 1,NEQ
    IDI = IWORK(LID-1+I)
    IF (IDI .NE. 1 .AND. IDI .NE. -1) GO TO 724
50 CONTINUE
  ENDIF
!
!     Check to see that TOUT is different from T.
!
IF(TOUT .EQ. T)GO TO 719
!
!     Check HMAX.
!
IF(INFO(7) .NE. 0) THEN
   HMAX = RWORK(LHMAX)
   IF (HMAX .LE. 0.0D0) GO TO 710
   ENDIF
!
!     Initialize counters and other flags.
!
IWORK(LNST)=0
IWORK(LNRE)=0
IWORK(LNJE)=0
IWORK(LETF)=0
IWORK(LNCFN)=0
IWORK(LNNI)=0
IWORK(LNLI)=0
IWORK(LNPS)=0
IWORK(LNCFL)=0
IWORK(LKPRIN)=INFO(18)
IDID=1
GO TO 200
!
!-----------------------------------------------------------------------
!     This block is for continuation calls only.
!     Here we check INFO(1), and if the last step was interrupted,
!     we check whether appropriate action was taken.
!-----------------------------------------------------------------------
!
100 CONTINUE
IF(INFO(1).EQ.1)GO TO 110
ITEMP = 1
IF(INFO(1).NE.-1)GO TO 701
!
!     If we are here, the last step was interrupted by an error
!     condition from DDSTP, and appropriate action was not taken.
!     This is a fatal error.
!
call rprintf( &
& 'daspk--  warning.. the last step terminated with a negative')
call rprintfi1( &
& 'value of idid and no appropriate action was taken %i' &
& // char(0),idid)
call rexit('- run terminated')

RETURN
110 CONTINUE
!
!-----------------------------------------------------------------------
!     This block is executed on all calls.
!
!     Counters are saved for later checks of performance.
!     Then the error tolerance parameters are checked, and the
!     work array pointers are set.
!-----------------------------------------------------------------------
!
200 CONTINUE
!
!     Save counters for use later.
!
IWORK(LNSTL)=IWORK(LNST)
NLI0 = IWORK(LNLI)
NNI0 = IWORK(LNNI)
NCFN0 = IWORK(LNCFN)
NCFL0 = IWORK(LNCFL)
NWARN = 0
!
!     Check RTOL and ATOL.
!
NZFLG = 0
RTOLI = RTOL(1)
ATOLI = ATOL(1)
DO 210 I=1,NEQ
   IF (INFO(2) .EQ. 1) RTOLI = RTOL(I)
   IF (INFO(2) .EQ. 1) ATOLI = ATOL(I)
   IF (RTOLI .GT. 0.0D0 .OR. ATOLI .GT. 0.0D0) NZFLG = 1
   IF (RTOLI .LT. 0.0D0) GO TO 706
   IF (ATOLI .LT. 0.0D0) GO TO 707
210 CONTINUE
IF (NZFLG .EQ. 0) GO TO 708
!
!     Set pointers to RWORK and IWORK segments.
!     For direct methods, SAVR is not used.
!
IWORK(LLCIWP) = LID + LENID
LSAVR = LDELTA
IF (INFO(12) .NE. 0) LSAVR = LDELTA + NEQ
LE = LSAVR + NEQ
LWT = LE + NEQ
LVT = LWT
IF (INFO(16) .NE. 0) LVT = LWT + NEQ
LPHI = LVT + NEQ
LWM = LPHI + (IWORK(LMXORD)+1)*NEQ
IF (INFO(1) .EQ. 1) GO TO 400
!
!-----------------------------------------------------------------------
!     This block is executed on the initial call only.
!     Set the initial step size, the error weight vector, and PHI.
!     Compute unknown initial components of Y and YPRIME, if requested.
!-----------------------------------------------------------------------
!
CONTINUE
TN=T
IDID=1
!
!     Set error weight array WT and altered weight array VT.
!
CALL DDAWTS(NEQ,INFO(2),RTOL,ATOL,Y,RWORK(LWT),RPAR,IPAR)
! KARLINE      
IF(NIND(1) < NEQ)CALL SCALE(NEQ, NIND, RWORK(LWT), 1.d-1)  !H not known yet
CALL DINVWT(NEQ,RWORK(LWT),IER)
IF (IER .NE. 0) GO TO 713
IF (INFO(16) .NE. 0) THEN
  DO 305 I = 1, NEQ
    RWORK(LVT+I-1) = MAX(IWORK(LID+I-1),0)*RWORK(LWT+I-1)
305 CONTINUE
  ENDIF
!
!     Compute unit roundoff and HMIN.
!
UROUND = D1MACH(4)
RWORK(LROUND) = UROUND
HMIN = 4.0D0*UROUND*MAX(ABS(T),ABS(TOUT))
!
!     Set/check STPTOL control for initial condition calculation.
!     
IF (INFO(11) .NE. 0) THEN
  IF( INFO(17) .EQ. 0) THEN
    RWORK(LSTOL) = UROUND**.6667D0
  ELSE
    IF (RWORK(LSTOL) .LE. 0.0D0) GO TO 725
    ENDIF
  ENDIF
!
!     Compute EPCON and square root of NEQ and its reciprocal, used
!     inside iterative solver.
!
RWORK(LEPCON) = 0.33D0
FLOATN = NEQ
RWORK(LSQRN) = SQRT(FLOATN)
RWORK(LRSQRN) = 1.D0/RWORK(LSQRN)
!
!     Check initial interval to see that it is long enough.
!
TDIST = ABS(TOUT - T)
IF(TDIST .LT. HMIN) GO TO 714
!
!     Check H0, if this was input.
!
IF (INFO(8) .EQ. 0) GO TO 310
   H0 = RWORK(LH)
   IF ((TOUT - T)*H0 .LT. 0.0D0) GO TO 711
   IF (H0 .EQ. 0.0D0) GO TO 712
   GO TO 320
310 CONTINUE
!
!     Compute initial stepsize, to be used by either
!     DDSTP or DDASIC, depending on INFO(11).
!
H0 = 0.001D0*TDIST
YPNORM = DDWNRM(NEQ,YPRIME,RWORK(LVT),RPAR,IPAR)
IF (YPNORM .GT. 0.5D0/H0) H0 = 0.5D0/YPNORM
H0 = SIGN(H0,TOUT-T)
!
!     Adjust H0 if necessary to meet HMAX bound.
!
320 IF (INFO(7) .EQ. 0) GO TO 330
   RH = ABS(H0)/RWORK(LHMAX)
   IF (RH .GT. 1.0D0) H0 = H0/RH
!
!     Check against TSTOP, if applicable.
!
330 IF (INFO(4) .EQ. 0) GO TO 340
   TSTOP = RWORK(LTSTOP)
   IF ((TSTOP - T)*H0 .LT. 0.0D0) GO TO 715
   IF ((T + H0 - TSTOP)*H0 .GT. 0.0D0) H0 = TSTOP - T
   IF ((TSTOP - TOUT)*H0 .LT. 0.0D0) GO TO 709
!
340 IF (INFO(11) .EQ. 0) GO TO 370
!
!     Compute unknown components of initial Y and YPRIME, depending
!     on INFO(11) and INFO(12).  INFO(12) represents the nonlinear
!     solver type (direct/Krylov).  Pass the name of the specific 
!     nonlinear solver, depending on INFO(12).  The location of the work
!     arrays SAVR, YIC, YPIC, PWK also differ in the two cases.
!
NWT = 1
EPCONI = RWORK(LEPIN)*RWORK(LEPCON)
350 IF (INFO(12) .EQ. 0) THEN
   LYIC = LPHI + 2*NEQ
   LYPIC = LYIC + NEQ
   LPWK = LYPIC
   CALL DDASIC(TN,Y,YPRIME,NEQ,INFO(11),IWORK(LID), &
& RES,JAC,PSOL,H0,RWORK(LWT),NWT,IDID,RPAR,IPAR, &
& RWORK(LPHI),RWORK(LSAVR),RWORK(LDELTA),RWORK(LE), &
& RWORK(LYIC),RWORK(LYPIC),RWORK(LPWK),RWORK(LWM),IWORK(LIWM), &
& HMIN,RWORK(LROUND),RWORK(LEPLI),RWORK(LSQRN),RWORK(LRSQRN), &
& EPCONI,RWORK(LSTOL),INFO(15),ICNFLG,IWORK(LICNS),DDASID)
ELSE IF (INFO(12) .EQ. 1) THEN
   LYIC = LWM
   LYPIC = LYIC + NEQ
   LPWK = LYPIC + NEQ
   CALL DDASIC(TN,Y,YPRIME,NEQ,INFO(11),IWORK(LID), &
& RES,JAC,PSOL,H0,RWORK(LWT),NWT,IDID,RPAR,IPAR, &
& RWORK(LPHI),RWORK(LSAVR),RWORK(LDELTA),RWORK(LE), &
& RWORK(LYIC),RWORK(LYPIC),RWORK(LPWK),RWORK(LWM),IWORK(LIWM), &
& HMIN,RWORK(LROUND),RWORK(LEPLI),RWORK(LSQRN),RWORK(LRSQRN), &
& EPCONI,RWORK(LSTOL),INFO(15),ICNFLG,IWORK(LICNS),DDASIK)
ENDIF
!
IF (IDID .LT. 0) GO TO 600
!
!     DDASIC was successful.  If this was the first call to DDASIC,
!     update the WT array (with the current Y) and call it again.
!
IF (NWT .EQ. 2) GO TO 355
NWT = 2
CALL DDAWTS(NEQ,INFO(2),RTOL,ATOL,Y,RWORK(LWT),RPAR,IPAR)
! KARLINE      
IF(NIND(1) < NEQ) CALL SCALE(NEQ, NIND, RWORK(LWT), H0)
CALL DINVWT(NEQ,RWORK(LWT),IER)
IF (IER .NE. 0) GO TO 713
GO TO 350
!
!     If INFO(14) = 1, return now with IDID = 4.
!
355 IF (INFO(14) .EQ. 1) THEN
  IDID = 4
  H = H0
  IF (INFO(11) .EQ. 1) RWORK(LHOLD) = H0
  GO TO 590
ENDIF
!
!     Update the WT and VT arrays one more time, with the new Y.
!
CALL DDAWTS(NEQ,INFO(2),RTOL,ATOL,Y,RWORK(LWT),RPAR,IPAR)
! KARLINE      
IF(NIND(1) < NEQ) CALL SCALE(NEQ, NIND, RWORK(LWT), H0)
CALL DINVWT(NEQ,RWORK(LWT),IER)
IF (IER .NE. 0) GO TO 713
IF (INFO(16) .NE. 0) THEN
  DO 357 I = 1, NEQ
    RWORK(LVT+I-1) = MAX(IWORK(LID+I-1),0)*RWORK(LWT+I-1)
357 CONTINUE
  ENDIF
!
!     Reset the initial stepsize to be used by DDSTP.
!     Use H0, if this was input.  Otherwise, recompute H0,
!     and adjust it if necessary to meet HMAX bound.
!
IF (INFO(8) .NE. 0) THEN
   H0 = RWORK(LH)
   GO TO 360
   ENDIF
!
H0 = 0.001D0*TDIST
YPNORM = DDWNRM(NEQ,YPRIME,RWORK(LVT),RPAR,IPAR)
IF (YPNORM .GT. 0.5D0/H0) H0 = 0.5D0/YPNORM
H0 = SIGN(H0,TOUT-T)
!
360 IF (INFO(7) .NE. 0) THEN
   RH = ABS(H0)/RWORK(LHMAX)
   IF (RH .GT. 1.0D0) H0 = H0/RH
   ENDIF
!
!     Check against TSTOP, if applicable.
!
IF (INFO(4) .NE. 0) THEN
   TSTOP = RWORK(LTSTOP)
   IF ((T + H0 - TSTOP)*H0 .GT. 0.0D0) H0 = TSTOP - T
   ENDIF
!
!     Load H and RWORK(LH) with H0.
!
370 H = H0
RWORK(LH) = H
!
!     Load Y and H*YPRIME into PHI(*,1) and PHI(*,2).
!
ITEMP = LPHI + NEQ
DO 380 I = 1,NEQ
   RWORK(LPHI + I - 1) = Y(I)
   RWORK(ITEMP + I - 1) = H*YPRIME(I)
380 CONTINUE
!
GO TO 500
!
!-----------------------------------------------------------------------
!     This block is for continuation calls only.
!     Its purpose is to check stop conditions before taking a step.
!     Adjust H if necessary to meet HMAX bound.
!-----------------------------------------------------------------------
!
400 CONTINUE
UROUND=RWORK(LROUND)
DONE = .FALSE.
TN=RWORK(LTN)
H=RWORK(LH)
IF(INFO(7) .EQ. 0) GO TO 410
   RH = ABS(H)/RWORK(LHMAX)
   IF(RH .GT. 1.0D0) H = H/RH
410 CONTINUE
IF(T .EQ. TOUT) GO TO 719
IF((T - TOUT)*H .GT. 0.0D0) GO TO 711
IF(INFO(4) .EQ. 1) GO TO 430
IF(INFO(3) .EQ. 1) GO TO 420
IF((TN-TOUT)*H.LT.0.0D0)GO TO 490
CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ,IWORK(LKOLD), &
& RWORK(LPHI),RWORK(LPSI))
T=TOUT
IDID = 3
DONE = .TRUE.
GO TO 490
420 IF((TN-T)*H .LE. 0.0D0) GO TO 490
IF((TN - TOUT)*H .GT. 0.0D0) GO TO 425
CALL DDATRP(TN,TN,Y,YPRIME,NEQ,IWORK(LKOLD), &
& RWORK(LPHI),RWORK(LPSI))
T = TN
IDID = 1
DONE = .TRUE.
GO TO 490
425 CONTINUE
CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ,IWORK(LKOLD), &
& RWORK(LPHI),RWORK(LPSI))
T = TOUT
IDID = 3
DONE = .TRUE.
GO TO 490
430 IF(INFO(3) .EQ. 1) GO TO 440
TSTOP=RWORK(LTSTOP)
IF((TN-TSTOP)*H.GT.0.0D0) GO TO 715
IF((TSTOP-TOUT)*H.LT.0.0D0)GO TO 709
IF((TN-TOUT)*H.LT.0.0D0)GO TO 450
CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ,IWORK(LKOLD), &
& RWORK(LPHI),RWORK(LPSI))
T=TOUT
IDID = 3
DONE = .TRUE.
GO TO 490
440 TSTOP = RWORK(LTSTOP)
IF((TN-TSTOP)*H .GT. 0.0D0) GO TO 715
IF((TSTOP-TOUT)*H .LT. 0.0D0) GO TO 709
IF((TN-T)*H .LE. 0.0D0) GO TO 450
IF((TN - TOUT)*H .GT. 0.0D0) GO TO 445
CALL DDATRP(TN,TN,Y,YPRIME,NEQ,IWORK(LKOLD), &
& RWORK(LPHI),RWORK(LPSI))
T = TN
IDID = 1
DONE = .TRUE.
GO TO 490
445 CONTINUE
CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ,IWORK(LKOLD), &
& RWORK(LPHI),RWORK(LPSI))
T = TOUT
IDID = 3
DONE = .TRUE.
GO TO 490
450 CONTINUE
!
!     Check whether we are within roundoff of TSTOP.
!
IF(ABS(TN-TSTOP).GT.100.0D0*UROUND* &
& (ABS(TN)+ABS(H)))GO TO 460
CALL DDATRP(TN,TSTOP,Y,YPRIME,NEQ,IWORK(LKOLD), &
& RWORK(LPHI),RWORK(LPSI))
IDID=2
T=TSTOP
DONE = .TRUE.
GO TO 490
460 TNEXT=TN+H
IF((TNEXT-TSTOP)*H.LE.0.0D0)GO TO 490
H=TSTOP-TN
RWORK(LH)=H
!
490 IF (DONE) GO TO 590
!
!-----------------------------------------------------------------------
!     The next block contains the call to the one-step integrator DDSTP.
!     This is a looping point for the integration steps.
!     Check for too many steps.
!     Check for poor Newton/Krylov performance.
!     Update WT.  Check for too much accuracy requested.
!     Compute minimum stepsize.
!-----------------------------------------------------------------------
!
500 CONTINUE
!
!     Check for too many steps.
!
IF((IWORK(LNST)-IWORK(LNSTL)).LT.500) GO TO 505
     IDID=-1
     GO TO 527
!
! Check for poor Newton/Krylov performance.
!
505 IF (INFO(12) .EQ. 0) GO TO 510
NSTD = IWORK(LNST) - IWORK(LNSTL)
NNID = IWORK(LNNI) - NNI0
IF (NSTD .LT. 10 .OR. NNID .EQ. 0) GO TO 510
AVLIN = REAL(IWORK(LNLI) - NLI0)/REAL(NNID)
RCFN = REAL(IWORK(LNCFN) - NCFN0)/REAL(NSTD)
RCFL = REAL(IWORK(LNCFL) - NCFL0)/REAL(NNID)
FMAXL = IWORK(LMAXL)
LAVL = AVLIN .GT. FMAXL
LCFN = RCFN .GT. 0.9D0
LCFL = RCFL .GT. 0.9D0
LWARN = LAVL .OR. LCFN .OR. LCFL
IF (.NOT.LWARN) GO TO 510
NWARN = NWARN + 1
IF (NWARN .GT. 10) GO TO 510
IF (LAVL) THEN
call rprintf( &
& 'daspk--  warning.. Poor iterative algorithm performance' &
& // char(0))
call rprintfd2( &
& '      at T = R1. Average no. of linear iterations = R2' &
& // ' %g, %g' // char(0), TN, AVLIN)
  ENDIF
IF (LCFN) THEN
call rprintf( &
& 'daspk--  warning.. Poor iterative algorithm performance ' &
& // char(0))
call rprintfd2( &
& '      at T = R1.  Nonlinear convergence failure rate = R2' &
& // '%g, %g' // char(0), TN, RCFN)

  ENDIF

IF (LCFL) THEN
call rprintf( &
& 'daspk--  warning.. Poor iterative algorithm performance ' &
& // char(0))
call rprintfd2( &
& '      at T = R1. Linear convergence failure rate = R2 ' &
& // char(0), TN, RCFL)

  ENDIF
!
!     Update WT and VT, if this is not the first call.
!
510 CALL DDAWTS(NEQ,INFO(2),RTOL,ATOL,RWORK(LPHI),RWORK(LWT), &
& RPAR,IPAR)
IF(NIND(1) < NEQ) CALL SCALE(NEQ, NIND, RWORK(LWT), H)
CALL DINVWT(NEQ,RWORK(LWT),IER)
IF (IER .NE. 0) THEN
  IDID = -3
  GO TO 527
  ENDIF
IF (INFO(16) .NE. 0) THEN
  DO 515 I = 1, NEQ
    RWORK(LVT+I-1) = MAX(IWORK(LID+I-1),0)*RWORK(LWT+I-1)
515 CONTINUE
ENDIF
!
!     Test for too much accuracy requested.
!
R = DDWNRM(NEQ,RWORK(LPHI),RWORK(LWT),RPAR,IPAR)*100.0D0*UROUND
IF (R .LE. 1.0D0) GO TO 525
!
!     Multiply RTOL and ATOL by R and return.
!
IF(INFO(2).EQ.1)GO TO 523
     RTOL(1)=R*RTOL(1)
     ATOL(1)=R*ATOL(1)
     IDID=-2
     GO TO 527
523 DO 524 I=1,NEQ
     RTOL(I)=R*RTOL(I)
     ATOL(I)=R*ATOL(I)
524 CONTINUE
IDID=-2
GO TO 527
525 CONTINUE
!
!     Compute minimum stepsize.
!
HMIN=4.0D0*UROUND*MAX(ABS(TN),ABS(TOUT))
!
!     Test H vs. HMAX
IF (INFO(7) .NE. 0) THEN
   RH = ABS(H)/RWORK(LHMAX)
   IF (RH .GT. 1.0D0) H = H/RH
   ENDIF
!
!     Call the one-step integrator.
!     Note that INFO(12) represents the nonlinear solver type.
!     Pass the required nonlinear solver, depending upon INFO(12).
!
IF (INFO(12) .EQ. 0) THEN
   CALL DDSTP(TN,Y,YPRIME,NEQ, &
& RES,JAC,PSOL,H,RWORK(LWT),RWORK(LVT),INFO(1),IDID,RPAR,IPAR, &
& RWORK(LPHI),RWORK(LSAVR),RWORK(LDELTA),RWORK(LE), &
& RWORK(LWM),IWORK(LIWM), &
& RWORK(LALPHA),RWORK(LBETA),RWORK(LGAMMA), &
& RWORK(LPSI),RWORK(LSIGMA), &
& RWORK(LCJ),RWORK(LCJOLD),RWORK(LHOLD),RWORK(LS),HMIN, &
& RWORK(LROUND), RWORK(LEPLI),RWORK(LSQRN),RWORK(LRSQRN), &
& RWORK(LEPCON), IWORK(LPHASE),IWORK(LJCALC),INFO(15), &
& IWORK(LK), IWORK(LKOLD),IWORK(LNS),NONNEG,INFO(12), &
& DNEDD)
ELSE IF (INFO(12) .EQ. 1) THEN
   CALL DDSTP(TN,Y,YPRIME,NEQ, &
& RES,JAC,PSOL,H,RWORK(LWT),RWORK(LVT),INFO(1),IDID,RPAR,IPAR, &
& RWORK(LPHI),RWORK(LSAVR),RWORK(LDELTA),RWORK(LE), &
& RWORK(LWM),IWORK(LIWM), &
& RWORK(LALPHA),RWORK(LBETA),RWORK(LGAMMA), &
& RWORK(LPSI),RWORK(LSIGMA), &
& RWORK(LCJ),RWORK(LCJOLD),RWORK(LHOLD),RWORK(LS),HMIN, &
& RWORK(LROUND), RWORK(LEPLI),RWORK(LSQRN),RWORK(LRSQRN), &
& RWORK(LEPCON), IWORK(LPHASE),IWORK(LJCALC),INFO(15), &
& IWORK(LK), IWORK(LKOLD),IWORK(LNS),NONNEG,INFO(12), &
& DNEDK)
ENDIF
!
527 IF(IDID.LT.0)GO TO 600
!
!-----------------------------------------------------------------------
!     This block handles the case of a successful return from DDSTP
!     (IDID=1).  Test for stop conditions.
!-----------------------------------------------------------------------
!
IF(INFO(4).NE.0)GO TO 540
     IF(INFO(3).NE.0)GO TO 530
       IF((TN-TOUT)*H.LT.0.0D0)GO TO 500
       CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ, &
& IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
       IDID=3
       T=TOUT
       GO TO 580
530 IF((TN-TOUT)*H.GE.0.0D0)GO TO 535
       T=TN
       IDID=1
       GO TO 580
535 CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ, &
& IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
       IDID=3
       T=TOUT
       GO TO 580
540 IF(INFO(3).NE.0)GO TO 550
IF((TN-TOUT)*H.LT.0.0D0)GO TO 542
   CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ, &
& IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
   T=TOUT
   IDID=3
   GO TO 580
542 IF(ABS(TN-TSTOP).LE.100.0D0*UROUND* &
& (ABS(TN)+ABS(H)))GO TO 545
TNEXT=TN+H
IF((TNEXT-TSTOP)*H.LE.0.0D0)GO TO 500
H=TSTOP-TN
GO TO 500
545 CALL DDATRP(TN,TSTOP,Y,YPRIME,NEQ, &
& IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
IDID=2
T=TSTOP
GO TO 580
550 IF((TN-TOUT)*H.GE.0.0D0)GO TO 555
IF(ABS(TN-TSTOP).LE.100.0D0*UROUND*(ABS(TN)+ABS(H)))GO TO 552
T=TN
IDID=1
GO TO 580
552 CALL DDATRP(TN,TSTOP,Y,YPRIME,NEQ, &
& IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
IDID=2
T=TSTOP
GO TO 580
555 CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ, &
& IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
T=TOUT
IDID=3
580 CONTINUE
!
!-----------------------------------------------------------------------
!     All successful returns from DDASPK are made from this block.
!-----------------------------------------------------------------------
!
590 CONTINUE
RWORK(LTN)=TN
RWORK(LH)=H
RETURN
!
!-----------------------------------------------------------------------
!     This block handles all unsuccessful returns other than for
!     illegal input.
!-----------------------------------------------------------------------
!
600 CONTINUE
ITEMP = -IDID
IF (ITEMP .EQ. 1) THEN
  GOTO 610
ELSE IF (ITEMP .EQ. 2) THEN
  GOTO 620
ELSE IF (ITEMP .EQ. 3) THEN
  GOTO 630
ELSE IF (ITEMP .EQ. 4) THEN
  GOTO 700
ELSE IF (ITEMP .EQ. 5) THEN
  GOTO 655
ELSE IF (ITEMP .EQ. 6) THEN
  GOTO 640
ELSE IF (ITEMP .EQ. 7) THEN
  GOTO 650
ELSE IF (ITEMP .EQ. 8) THEN
  GOTO 660
ELSE IF (ITEMP .EQ. 9) THEN
  GOTO 670
ELSE IF (ITEMP .EQ. 10) THEN
  GOTO 675
ELSE IF (ITEMP .EQ. 11) THEN
  GOTO 680
ELSE IF (ITEMP .EQ. 12) THEN
  GOTO 685
ELSE IF (ITEMP .EQ. 13) THEN
  GOTO 690
ELSE IF (ITEMP .EQ. 14) THEN
  GOTO 695
ENDIF

!      GO TO (610,620,630,700,655,640,650,660,670,675,
!     *  680,685,690,695), ITEMP
!
!     The maximum number of steps was taken before
!     reaching tout.
!! Karline toggled this off, version > 1.10.3 unless lots of printing requested
610 IF(IWORK(LKPRIN) .GE. 2) THEN
call rprintf( &
& 'daspk--  warning.. At current T (=R1) max number steps' &
& // char(0))
call rprintfd1( &
& '      on this call before reaching tout %g' // char(0), TN)
ENDIF
GO TO 700
!
!     Too much accuracy for machine precision.
!
620 call rprintf( &
& 'daspk--  warning.. At T(=R1) too much accuracy requested' &
& // char(0))
call rprintf( &
& '      for precision of machine. rtol and atol were' // char(0))
call rprintfd1( &
& '      increased to appropriate values %g' &
& // char(0),  TN )

GO TO 700
!
!     WT(I) .LE. 0.0D0 for some I (not at start of problem).
!
630 call rprintf( &
& 'daspk--  warning.. At T(=R1) some element of WT ' // char(0))
call rprintfd1( &
& '      has become less or equal than 0 %g' // char(0), TN )
GO TO 700
!
!     Error test failed repeatedly or with H=HMIN.
!
640 call rprintf( &
& 'daspk--  warning.. At T(=R1) and stepsize H (=R2)' &
& //char(0))
call rprintfd2( &
& '      error test failed repeatedly or with abs(H)=Hmin' &
& // ' %g, %g' // char(0), TN, H )
GO TO 700
!
!     Nonlinear solver failed to converge repeatedly or with H=HMIN.
!
650 call rprintf( &
& 'daspk--  warning.. At T(=R1) and stepsize H (=R2) the' &
& // char(0))
call rprintf( &
& '      nonlinear solver failed to converge ' // char(0))
call rprintfd2( &
& '      repeatedly of with abs (H) = HMIN &g, %g' &
& // char(0),  TN, H)
GO TO 700
!
!     The preconditioner had repeated failures.
!
655 call rprintf( &
& 'daspk--  warning.. At T(=R1) and stepsize H (=R2) the' &
& // char(0))
 call rprintfd2( &
& '      preconditioner had repeated failures %g, %g' &
& // char(0), TN, H )
GO TO 700
!
!     The iteration matrix is singular.
!
660 call rprintf( &
& 'daspk--  warning.. At T(=R1) and stepsize H (=R2) the' &
& // char(0))
call rprintfd2( &
& '      iteration matrix is singular %g, %g' &
& // char(0), TN, H)
GO TO 700
!
!     Nonlinear system failure preceded by error test failures.
!
670 call rprintf( &
& 'daspk--  warning.. At T(=R1) and stepsize H (=R2) the' &
& // char(0))
call rprintf( &
& '      nonlinear solver could not converge ' // char(0))
call rprintfd2( &
& '      Also the error test failed repeatedly %g, %g' &
& // char(0), TN, H )
GO TO 700
!
!     Nonlinear system failure because IRES = -1.
!
675 call rprintf( &
& 'daspk--  warning.. At T(=R1) and stepsize H (=R2) the' &
& // char(0))
call rprintf( &
& '      nonlinear system solver could not converge' &
& // char(0))
call rprintfd2( &
& '      because ires was equal to -1 %g, %g' &
& // char(0), TN, H)
GO TO 700
!
!     Failure because IRES = -2.
!
680 call rprintf( &
& 'daspk--  warning.. At T(=R1) and stepsize H (=R2)' // char(0))
call rprintfd2( &
& '      ires was equal to -2 &g, %g' // char(0), TN, H )
GO TO 700
!
!     Failed to compute initial YPRIME.
!
685 call rprintf( &
& 'daspk--  warning.. At T(=R1) and stepsize H (=R2) the' &
& // char(0))
call rprintfd2( &
& '      initial yprime could not be computed %g, %g' &
& // char(0), TN, H0 )
GO TO 700
!
!     Failure because IER was negative from PSOL.
!
690 call rprintf( &
& 'daspk--  warning.. At T(=R1) and stepsize H (=R2)' // char(0))
call rprintfd2( &
& '      IER was negative from psol %g, %g' // char(0), TN, H)
GO TO 700
!
!     Failure because the linear system solver could not converge.
!
695 call rprintf( &
& 'daspk--  warning.. At T(=R1) and stepsize H (=R2) the' &
& // char(0))
call rprintfd2( &
& '      linear system solver could not converge %g, %g' &
& // char(0), TN,H )
GO TO 700
!
!
700 CONTINUE
INFO(1)=-1
T=TN
RWORK(LTN)=TN
RWORK(LH)=H
RETURN
!
!-----------------------------------------------------------------------
!     This block handles all error returns due to illegal input,
!     as detected before calling DDSTP.
!     First the error message routine is called.  If this happens
!     twice in succession, execution is terminated.
!-----------------------------------------------------------------------
!
701 call rprintfi1( &
& 'daspk-- element (= %i) of info vector is not valid' &
& // char(0),ITEMP)
GO TO 750
702 call rprintfi1( &
& 'daspk-- neq (= %i) < 0' // char(0), NEQ)
GO TO 750
703 call rprintfi1( &
& 'daspk-- maxord (= %i) not in range' // char(0), MXORD)
GO TO 750
704 call rprintfi2( &
& 'daspk-- rwork length needed, LENRW(= %i) exceeds LRW(= %i)' &
& // char(0), LENRW, LRW)
GO TO 750
705 call rprintfi2( &
& 'daspk-- iwork length needed, LENIW(= %i) exceeds LIW(= %i)' &
& // char(0), LENIW, LIW)
GO TO 750
706 call rprintf( &
& 'daspk-- some element of rtol is < 0' // char(0))
GO TO 750
707 call rprintf( &
& 'daspk-- some element of atol is < 0' // char(0))
GO TO 750
708 call rprintf( &
& 'daspk-- all elements of rtol and atol are 0' // char(0))
GO TO 750
709 call rprintfd2( &
& 'daspk-- INFO(4)=1 and TSTOP (= %g) behind TOUT (= %g)' &
& // char(0),TSTOP,TOUT)
GO TO 750
710 call rprintfd1( &
& 'daspk-- HMAX (= %g) < 0' // char(0), HMAX)
GO TO 750
711 call rprintfd2( &
& 'daspk-- TOUT (= %g) behind T (= %g)' // char(0), TOUT, T)
GO TO 750
712 call rprintf( &
& 'daspk-- INFO(8)=1 and H0=0' // char(0))
GO TO 750
713 call rprintf( &
& 'daspk-- some element of WT <= 0 ' // char(0))
GO TO 750
714 call rprintfd2( &
& 'daspk-- TOUT (= %g) too close to T (= %g)' &
& // ' to start integration ' // char(0), TOUT, T)
GO TO 750
715 call rprintfd2( &
& 'daspk-- INFO(4)=1 and TSTOP (= %g) behind T (= %g)' &
& // char(0), TSTOP, T)
GO TO 750
717 call rprintfi1( &
& 'daspk-- ML (= %i) illegal - either < 0 or > neq' &
& // char(0),IWORK(LML))
GO TO 750
718 call rprintfi1( &
& 'daspk-- MU (= %i) illegal - either < 0 or > neq' &
& // char(0),IWORK(LMU))
GO TO 750
719 call rprintfd2( &
& 'daspk-- TOUT (= %g) is equal to T (= %g)' &
& // char(0),TOUT,T)
GO TO 750
720 call rprintfi1( &
& 'daspk-- MAXL (= %i) illegal - either < 1 or > neq' &
& // char(0), IWORK(LMAXL))
GO TO 750
721 call rprintfi1( &
& 'daspk-- KMP (= %i) illegal - either < 1 or > MAXL' &
& // char(0), IWORK(LKMP))
GO TO 750
722 call rprintfi1( &
& 'daspk-- NRMAX (= %i) illegal -  < 0 ' &
& // char(0),IWORK(LNRMAX))
GO TO 750
723 call rprintfd1( &
& 'daspk-- EPLI (= %g) illegal - either <= 0 or >= 1' &
& // char(0),RWORK(LEPLI))
GO TO 750
724 call rprintf( &
& 'daspk-- illegal IWORK value for INFO(11) not equal to 0' &
& // char(0))
GO TO 750
725 call rprintf( &
& 'daspk-- one of the inputs for INFO(17) = 1 is illegal' &
& // char(0))
GO TO 750
726 call rprintf( &
& 'daspk-- illegal IWORK value for INFO(10) not equal to 0' &
& // char(0))
GO TO 750
727 call rprintfi1( &
& 'daspk-- Y(I) and IWORK(40+I) (I= %i) inconsistent' &
& // char(0), IRET )
GO TO 750
750 IF(INFO(1).EQ.-1) GO TO 760
INFO(1)=-1
IDID=-33
RETURN
760 call rprintf( &
& 'daspk-- repeated occurrences of illegal input' // char(0))
call rprintf( &
& 'daspk-- run terminated; apparent infinite loop' // char(0))
RETURN
!
!------END OF SUBROUTINE DDASPK-----------------------------------------
END
! Work performed under the auspices of the U.S. Department of Energy
! by Lawrence Livermore National Laboratory under contract number 
! W-7405-Eng-48.
!
SUBROUTINE DDASIC (X, Y, YPRIME, NEQ, ICOPT, ID, RES, JAC, PSOL, &
& H, WT, NIC, IDID, RPAR, IPAR, PHI, SAVR, DELTA, E, YIC, YPIC, &
& PWK, WM, IWM, HMIN, UROUND, EPLI, SQRTN, RSQRTN, EPCONI, &
& STPTOL, JFLG, ICNFLG, ICNSTR, NLSIC)
!
!***BEGIN PROLOGUE  DDASIC
!***REFER TO  DDASPK
!***DATE WRITTEN   940628   (YYMMDD)
!***REVISION DATE  941206   (YYMMDD)
!***REVISION DATE  950714   (YYMMDD)
!
!-----------------------------------------------------------------------
!***DESCRIPTION
!
!     DDASIC is a driver routine to compute consistent initial values
!     for Y and YPRIME.  There are two different options:  
!     Denoting the differential variables in Y by Y_d, and
!     the algebraic variables by Y_a, the problem solved is either:
!     1.  Given Y_d, calculate Y_a and Y_d', or
!     2.  Given Y', calculate Y.
!     In either case, initial values for the given components
!     are input, and initial guesses for the unknown components
!     must also be provided as input.
!
!     The external routine NLSIC solves the resulting nonlinear system.
!
!     The parameters represent
!
!     X  --        Independent variable.
!     Y  --        Solution vector at X.
!     YPRIME --    Derivative of solution vector.
!     NEQ --       Number of equations to be integrated.
!     ICOPT     -- Flag indicating initial condition option chosen.
!                    ICOPT = 1 for option 1 above.
!                    ICOPT = 2 for option 2.
!     ID        -- Array of dimension NEQ, which must be initialized
!                  if option 1 is chosen.
!                    ID(i) = +1 if Y_i is a differential variable,
!                    ID(i) = -1 if Y_i is an algebraic variable. 
!     RES --       External user-supplied subroutine to evaluate the
!                  residual.  See RES description in DDASPK prologue.
!     JAC --       External user-supplied routine to update Jacobian
!                  or preconditioner information in the nonlinear solver
!                  (optional).  See JAC description in DDASPK prologue.
!     PSOL --      External user-supplied routine to solve
!                  a linear system using preconditioning. 
!                  See PSOL in DDASPK prologue.
!     H --         Scaling factor in iteration matrix.  DDASIC may 
!                  reduce H to achieve convergence.
!     WT --        Vector of weights for error criterion.
!     NIC --       Input number of initial condition calculation call 
!                  (= 1 or 2).
!     IDID --      Completion code.  See IDID in DDASPK prologue.
!     RPAR,IPAR -- Real and integer parameter arrays that
!                  are used for communication between the
!                  calling program and external user routines.
!                  They are not altered by DNSK
!     PHI --       Work space for DDASIC of length at least 2*NEQ.
!     SAVR --      Work vector for DDASIC of length NEQ.
!     DELTA --     Work vector for DDASIC of length NEQ.
!     E --         Work vector for DDASIC of length NEQ.
!     YIC,YPIC --  Work vectors for DDASIC, each of length NEQ.
!     PWK --       Work vector for DDASIC of length NEQ.
!     WM,IWM --    Real and integer arrays storing
!                  information required by the linear solver.
!     EPCONI --    Test constant for Newton iteration convergence.
!     ICNFLG --    Flag showing whether constraints on Y are to apply.
!     ICNSTR --    Integer array of length NEQ with constraint types.
!
!     The other parameters are for use internally by DDASIC.
!
!-----------------------------------------------------------------------
!***ROUTINES CALLED
!   DCOPY, NLSIC
!
!***END PROLOGUE  DDASIC
!
!
IMPLICIT DOUBLE PRECISION(A-H,O-Z)
DIMENSION Y(*),YPRIME(*),ID(*),WT(*),PHI(NEQ,*)
DIMENSION SAVR(*),DELTA(*),E(*),YIC(*),YPIC(*),PWK(*)
DIMENSION WM(*),IWM(*), RPAR(*),IPAR(*), ICNSTR(*)
EXTERNAL RES, JAC, PSOL, NLSIC
!
INTEGER, PARAMETER :: LCFN=15
INTEGER, PARAMETER :: LMXNH=34
!
! The following parameters are data-loaded here:
!     RHCUT  = factor by which H is reduced on retry of Newton solve.
!     RATEMX = maximum convergence rate for which Newton iteration
!              is considered converging.
!
SAVE RHCUT, RATEMX
DATA RHCUT/0.1D0/, RATEMX/0.8D0/
!
!
!-----------------------------------------------------------------------
!     BLOCK 1.
!     Initializations.
!     JSKIP is a flag set to 1 when NIC = 2 and NH = 1, to signal that
!     the initial call to the JAC routine is to be skipped then.
!     Save Y and YPRIME in PHI.  Initialize IDID, NH, and CJ.
!-----------------------------------------------------------------------
!
MXNH = IWM(LMXNH)
IDID = 1
NH = 1
JSKIP = 0
IF (NIC .EQ. 2) JSKIP = 1
CALL DCOPY (NEQ, Y, 1, PHI(1,1), 1)
CALL DCOPY (NEQ, YPRIME, 1, PHI(1,2), 1)
!
IF (ICOPT .EQ. 2) THEN
  CJ = 0.0D0
ELSE
  CJ = 1.0D0/H
ENDIF
!
!-----------------------------------------------------------------------
!     BLOCK 2
!     Call the nonlinear system solver to obtain
!     consistent initial values for Y and YPRIME.
!-----------------------------------------------------------------------
!
200 CONTINUE
CALL NLSIC(X,Y,YPRIME,NEQ,ICOPT,ID,RES,JAC,PSOL,H,WT,JSKIP, &
& RPAR,IPAR,SAVR,DELTA,E,YIC,YPIC,PWK,WM,IWM,CJ,UROUND, &
& EPLI,SQRTN,RSQRTN,EPCONI,RATEMX,STPTOL,JFLG,ICNFLG,ICNSTR, &
& IERNLS)
!
IF (IERNLS .EQ. 0) RETURN
!
!-----------------------------------------------------------------------
!     BLOCK 3
!     The nonlinear solver was unsuccessful.  Increment NCFN.
!     Return with IDID = -12 if either
!       IERNLS = -1: error is considered unrecoverable,
!       ICOPT = 2: we are doing initialization problem type 2, or
!       NH = MXNH: the maximum number of H values has been tried.
!     Otherwise (problem 1 with IERNLS .GE. 1), reduce H and try again.
!     If IERNLS > 1, restore Y and YPRIME to their original values.
!-----------------------------------------------------------------------
!
IWM(LCFN) = IWM(LCFN) + 1
JSKIP = 0
!
IF (IERNLS .EQ. -1) GO TO 350
IF (ICOPT .EQ. 2) GO TO 350
IF (NH .EQ. MXNH) GO TO 350
!
NH = NH + 1
H = H*RHCUT
CJ = 1.0D0/H
!
IF (IERNLS .EQ. 1) GO TO 200
!
CALL DCOPY (NEQ, PHI(1,1), 1, Y, 1)
CALL DCOPY (NEQ, PHI(1,2), 1, YPRIME, 1)
GO TO 200
!
350 IDID = -12
RETURN
!
!------END OF SUBROUTINE DDASIC-----------------------------------------
END
! Work performed under the auspices of the U.S. Department of Energy
! by Lawrence Livermore National Laboratory under contract number 
! W-7405-Eng-48.
!
SUBROUTINE DYYPNW (NEQ, Y, YPRIME, CJ, RL, P, ICOPT, ID, &
& YNEW, YPNEW)
!
!***BEGIN PROLOGUE  DYYPNW
!***REFER TO  DLINSK
!***DATE WRITTEN   940830   (YYMMDD)
!
!
!-----------------------------------------------------------------------
!***DESCRIPTION
!
!     DYYPNW calculates the new (Y,YPRIME) pair needed in the
!     linesearch algorithm based on the current lambda value.  It is
!     called by DLINSK and DLINSD.  Based on the ICOPT and ID values,
!     the corresponding entry in Y or YPRIME is updated.
!
!     In addition to the parameters described in the calling programs,
!     the parameters represent
!
!     P      -- Array of length NEQ that contains the current
!               approximate Newton step.
!     RL     -- Scalar containing the current lambda value.
!     YNEW   -- Array of length NEQ containing the updated Y vector.
!     YPNEW  -- Array of length NEQ containing the updated YPRIME
!               vector.
!-----------------------------------------------------------------------
!
!***ROUTINES CALLED (NONE)
!
!***END PROLOGUE  DYYPNW
!
!
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
DIMENSION Y(*), YPRIME(*), YNEW(*), YPNEW(*), ID(*), P(*)
!
IF (ICOPT .EQ. 1) THEN
   DO 10 I=1,NEQ
      IF(ID(I) .LT. 0) THEN
         YNEW(I) = Y(I) - RL*P(I)
         YPNEW(I) = YPRIME(I)
      ELSE
         YNEW(I) = Y(I)
         YPNEW(I) = YPRIME(I) - RL*CJ*P(I)
      ENDIF
10 CONTINUE
ELSE
   DO 20 I = 1,NEQ
      YNEW(I) = Y(I) - RL*P(I)
      YPNEW(I) = YPRIME(I)
20 CONTINUE
ENDIF
RETURN
!----------------------- END OF SUBROUTINE DYYPNW ----------------------
END
! Work performed under the auspices of the U.S. Department of Energy
! by Lawrence Livermore National Laboratory under contract number 
! W-7405-Eng-48.
!
SUBROUTINE DDSTP(X,Y,YPRIME,NEQ,RES,JAC,PSOL,H,WT,VT, &
& JSTART,IDID,RPAR,IPAR,PHI,SAVR,DELTA,E,WM,IWM, &
& ALPHA,BETA,GAMMA,PSI,SIGMA,CJ,CJOLD,HOLD,S,HMIN,UROUND, &
& EPLI,SQRTN,RSQRTN,EPCON,IPHASE,JCALC,JFLG,K,KOLD,NS,NONNEG, &
& NTYPE,NLS)
!
!***BEGIN PROLOGUE  DDSTP
!***REFER TO  DDASPK
!***DATE WRITTEN   890101   (YYMMDD)
!***REVISION DATE  900926   (YYMMDD)
!***REVISION DATE  940909   (YYMMDD) (Reset PSI(1), PHI(*,2) at 690)
!
!
!-----------------------------------------------------------------------
!***DESCRIPTION
!
!     DDSTP solves a system of differential/algebraic equations of 
!     the form G(X,Y,YPRIME) = 0, for one step (normally from X to X+H).
!
!     The methods used are modified divided difference, fixed leading 
!     coefficient forms of backward differentiation formulas.  
!     The code adjusts the stepsize and order to control the local error
!     per step.
!
!
!     The parameters represent
!     X  --        Independent variable.
!     Y  --        Solution vector at X.
!     YPRIME --    Derivative of solution vector
!                  after successful step.
!     NEQ --       Number of equations to be integrated.
!     RES --       External user-supplied subroutine
!                  to evaluate the residual.  See RES description
!                  in DDASPK prologue.
!     JAC --       External user-supplied routine to update
!                  Jacobian or preconditioner information in the
!                  nonlinear solver.  See JAC description in DDASPK
!                  prologue.
!     PSOL --      External user-supplied routine to solve
!                  a linear system using preconditioning. 
!                  (This is optional).  See PSOL in DDASPK prologue.
!     H --         Appropriate step size for next step.
!                  Normally determined by the code.
!     WT --        Vector of weights for error criterion used in Newton test.
!     VT --        Masked vector of weights used in error test.
!     JSTART --    Integer variable set 0 for
!                  first step, 1 otherwise.
!     IDID --      Completion code returned from the nonlinear solver.
!                  See IDID description in DDASPK prologue.
!     RPAR,IPAR -- Real and integer parameter arrays that
!                  are used for communication between the
!                  calling program and external user routines.
!                  They are not altered by DNSK
!     PHI --       Array of divided differences used by
!                  DDSTP. The length is NEQ*(K+1), where
!                  K is the maximum order.
!     SAVR --      Work vector for DDSTP of length NEQ.
!     DELTA,E --   Work vectors for DDSTP of length NEQ.
!     WM,IWM --    Real and integer arrays storing
!                  information required by the linear solver.
!
!     The other parameters are information
!     which is needed internally by DDSTP to
!     continue from step to step.
!
!-----------------------------------------------------------------------
!***ROUTINES CALLED
!   NLS, DDWNRM, DDATRP
!
!***END PROLOGUE  DDSTP
!
!
IMPLICIT DOUBLE PRECISION(A-H,O-Z)
DIMENSION Y(*),YPRIME(*),WT(*),VT(*)
DIMENSION PHI(NEQ,*),SAVR(*),DELTA(*),E(*)
DIMENSION WM(*),IWM(*)
DIMENSION PSI(*),ALPHA(*),BETA(*),GAMMA(*),SIGMA(*)
DIMENSION RPAR(*),IPAR(*)
EXTERNAL  RES, JAC, PSOL, NLS
!
INTEGER, PARAMETER :: LMXORD=3
INTEGER, PARAMETER :: LNST=11, LETF=14, LCFN=15
!
!
!-----------------------------------------------------------------------
!     BLOCK 1.
!     Initialize.  On the first call, set
!     the order to 1 and initialize
!     other variables.
!-----------------------------------------------------------------------
!
!     Initializations for all calls
!
TERKM1 = 0.D0 ! KARLINE: INITIALISED TO AVOID WARNING
KNEW = 1      ! INITIALISED
ERKM1 = 1.D0  ! INITIALISED
EST = 0.D0    ! INITIALISED - should have no effect
XOLD=X
NCF=0
NEF=0
IF(JSTART .NE. 0) GO TO 120
!
!     If this is the first step, perform
!     other initializations
!
K=1
KOLD=0
HOLD=0.0D0
PSI(1)=H
CJ = 1.D0/H
IPHASE = 0
NS=0
120 CONTINUE
!
!
!
!
!
!-----------------------------------------------------------------------
!     BLOCK 2
!     Compute coefficients of formulas for
!     this step.
!-----------------------------------------------------------------------
200 CONTINUE
KP1=K+1
KP2=K+2
KM1=K-1
IF(H.NE.HOLD.OR.K .NE. KOLD) NS = 0
NS=MIN0(NS+1,KOLD+2)
NSP1=NS+1
IF(KP1 .LT. NS)GO TO 230
!
BETA(1)=1.0D0
ALPHA(1)=1.0D0
TEMP1=H
GAMMA(1)=0.0D0
SIGMA(1)=1.0D0
DO 210 I=2,KP1
   TEMP2=PSI(I-1)
   PSI(I-1)=TEMP1
   BETA(I)=BETA(I-1)*PSI(I-1)/TEMP2
   TEMP1=TEMP2+H
   ALPHA(I)=H/TEMP1
   SIGMA(I)=(I-1)*SIGMA(I-1)*ALPHA(I)
   GAMMA(I)=GAMMA(I-1)+ALPHA(I-1)/H
210 CONTINUE
PSI(KP1)=TEMP1
230 CONTINUE
!
!     Compute ALPHAS, ALPHA0
!
ALPHAS = 0.0D0
ALPHA0 = 0.0D0
DO 240 I = 1,K
  ALPHAS = ALPHAS - 1.0D0/I
  ALPHA0 = ALPHA0 - ALPHA(I)
240 CONTINUE
!
!     Compute leading coefficient CJ
!
CJLAST = CJ
CJ = -ALPHAS/H
!
!     Compute variable stepsize error coefficient CK
!
CK = ABS(ALPHA(KP1) + ALPHAS - ALPHA0)
CK = MAX(CK,ALPHA(KP1))
!
!     Change PHI to PHI STAR
!
IF(KP1 .LT. NSP1) GO TO 280
DO 270 J=NSP1,KP1
   DO 260 I=1,NEQ
      PHI(I,J)=BETA(J)*PHI(I,J)
260 CONTINUE
270 CONTINUE
280 CONTINUE
!
!     Update time
!
X=X+H
!
!     Initialize IDID to 1
!
IDID = 1
!
!
!
!
!
!-----------------------------------------------------------------------
!     BLOCK 3
!     Call the nonlinear system solver to obtain the solution and
!     derivative.
!-----------------------------------------------------------------------
!
CALL NLS(X,Y,YPRIME,NEQ, &
& RES,JAC,PSOL,H,WT,JSTART,IDID,RPAR,IPAR,PHI,GAMMA, &
& SAVR,DELTA,E,WM,IWM,CJ,CJOLD,CJLAST,S, &
& UROUND,EPLI,SQRTN,RSQRTN,EPCON,JCALC,JFLG,KP1, &
& NONNEG,NTYPE,IERNLS)
!
IF(IERNLS .NE. 0)GO TO 600
!
!
!
!
!
!-----------------------------------------------------------------------
!     BLOCK 4
!     Estimate the errors at orders K,K-1,K-2
!     as if constant stepsize was used. Estimate
!     the local error at order K and test
!     whether the current step is successful.
!-----------------------------------------------------------------------
!
!     Estimate errors at orders K,K-1,K-2
!
ENORM = DDWNRM(NEQ,E,VT,RPAR,IPAR)
ERK = SIGMA(K+1)*ENORM
TERK = (K+1)*ERK
EST = ERK
KNEW=K
IF(K .EQ. 1)GO TO 430
DO 405 I = 1,NEQ
  DELTA(I) = PHI(I,KP1) + E(I)
405 CONTINUE
ERKM1=SIGMA(K)*DDWNRM(NEQ,DELTA,VT,RPAR,IPAR)
TERKM1 = K*ERKM1
IF(K .GT. 2)GO TO 410
IF(TERKM1 .LE. 0.5*TERK)GO TO 420
GO TO 430
410 CONTINUE
DO 415 I = 1,NEQ
  DELTA(I) = PHI(I,K) + DELTA(I)
415 CONTINUE
ERKM2=SIGMA(K-1)*DDWNRM(NEQ,DELTA,VT,RPAR,IPAR)
TERKM2 = (K-1)*ERKM2
IF(MAX(TERKM1,TERKM2).GT.TERK)GO TO 430
!
!     Lower the order
!
420 CONTINUE
KNEW=K-1
EST = ERKM1
!
!
!     Calculate the local error for the current step
!     to see if the step was successful
!
430 CONTINUE
ERR = CK * ENORM
IF(ERR .GT. 1.0D0)GO TO 600
!
!
!
!
!
!-----------------------------------------------------------------------
!     BLOCK 5
!     The step is successful. Determine
!     the best order and stepsize for
!     the next step. Update the differences
!     for the next step.
!-----------------------------------------------------------------------
IDID=1
IWM(LNST)=IWM(LNST)+1
KDIFF=K-KOLD
KOLD=K
HOLD=H
!
!
!     Estimate the error at order K+1 unless
!        already decided to lower order, or
!        already using maximum order, or
!        stepsize not constant, or
!        order raised in previous step
!
IF(KNEW.EQ.KM1.OR.K.EQ.IWM(LMXORD))IPHASE=1
IF(IPHASE .EQ. 0)GO TO 545
IF(KNEW.EQ.KM1)GO TO 540
IF(K.EQ.IWM(LMXORD)) GO TO 550
IF(KP1.GE.NS.OR.KDIFF.EQ.1)GO TO 550
DO 510 I=1,NEQ
   DELTA(I)=E(I)-PHI(I,KP2)
510 CONTINUE
ERKP1 = (1.0D0/(K+2))*DDWNRM(NEQ,DELTA,VT,RPAR,IPAR)
TERKP1 = (K+2)*ERKP1
IF(K.GT.1)GO TO 520
IF(TERKP1.GE.0.5D0*TERK)GO TO 550
GO TO 530
520 IF(TERKM1.LE.MIN(TERK,TERKP1))GO TO 540
IF(TERKP1.GE.TERK.OR.K.EQ.IWM(LMXORD))GO TO 550
!
!     Raise order
!
530 K=KP1
EST = ERKP1
GO TO 550
!
!     Lower order
!
540 K=KM1
EST = ERKM1
GO TO 550
!
!     If IPHASE = 0, increase order by one and multiply stepsize by
!     factor two
!
545 K = KP1
HNEW = H*2.0D0
H = HNEW
GO TO 575
!
!
!     Determine the appropriate stepsize for
!     the next step.
!
550 HNEW=H
TEMP2=K+1
R=(2.0D0*EST+0.0001D0)**(-1.0D0/TEMP2)
IF(R .LT. 2.0D0) GO TO 555
HNEW = 2.0D0*H
GO TO 560
555 IF(R .GT. 1.0D0) GO TO 560
R = MAX(0.5D0,MIN(0.9D0,R))
HNEW = H*R
560 H=HNEW
!
!
!     Update differences for next step
!
575 CONTINUE
IF(KOLD.EQ.IWM(LMXORD))GO TO 585
DO 580 I=1,NEQ
   PHI(I,KP2)=E(I)
580 CONTINUE
585 CONTINUE
DO 590 I=1,NEQ
   PHI(I,KP1)=PHI(I,KP1)+E(I)
590 CONTINUE
DO 596 J1=2,KP1
   J=KP1-J1+1
   DO 595 I=1,NEQ
    PHI(I,J)=PHI(I,J)+PHI(I,J+1)
595 CONTINUE
596 CONTINUE
JSTART = 1
RETURN
!
!
!
!
!
!-----------------------------------------------------------------------
!     BLOCK 6
!     The step is unsuccessful. Restore X,PSI,PHI
!     Determine appropriate stepsize for
!     continuing the integration, or exit with
!     an error flag if there have been many
!     failures.
!-----------------------------------------------------------------------
600 IPHASE = 1
!
!     Restore X,PHI,PSI
!
X=XOLD
IF(KP1.LT.NSP1)GO TO 630
DO 620 J=NSP1,KP1
   TEMP1=1.0D0/BETA(J)
   DO 610 I=1,NEQ
      PHI(I,J)=TEMP1*PHI(I,J)
610 CONTINUE
620 CONTINUE
630 CONTINUE
DO 640 I=2,KP1
   PSI(I-1)=PSI(I)-H
640 CONTINUE
!
!
!     Test whether failure is due to nonlinear solver
!     or error test
!
IF(IERNLS .EQ. 0)GO TO 660
IWM(LCFN)=IWM(LCFN)+1
!
!
!     The nonlinear solver failed to converge.
!     Determine the cause of the failure and take appropriate action.
!     If IERNLS .LT. 0, then return.  Otherwise, reduce the stepsize
!     and try again, unless too many failures have occurred.
!
IF (IERNLS .LT. 0) GO TO 675
NCF = NCF + 1
R = 0.25D0
H = H*R
IF (NCF .LT. 10 .AND. ABS(H) .GE. HMIN) GO TO 690
IF (IDID .EQ. 1) IDID = -7
IF (NEF .GE. 3) IDID = -9
GO TO 675
!
!
!     The nonlinear solver converged, and the cause
!     of the failure was the error estimate
!     exceeding the tolerance.
!
660 NEF=NEF+1
IWM(LETF)=IWM(LETF)+1
IF (NEF .GT. 1) GO TO 665
!
!     On first error test failure, keep current order or lower
!     order by one.  Compute new stepsize based on differences
!     of the solution.
!
K = KNEW
TEMP2 = K + 1
R = 0.90D0*(2.0D0*EST+0.0001D0)**(-1.0D0/TEMP2)
R = MAX(0.25D0,MIN(0.9D0,R))
H = H*R
IF (ABS(H) .GE. HMIN) GO TO 690
IDID = -6
GO TO 675
!
!     On second error test failure, use the current order or
!     decrease order by one.  Reduce the stepsize by a factor of
!     one quarter.
!
665 IF (NEF .GT. 2) GO TO 670
K = KNEW
R = 0.25D0
H = R*H
IF (ABS(H) .GE. HMIN) GO TO 690
IDID = -6
GO TO 675
!
!     On third and subsequent error test failures, set the order to
!     one, and reduce the stepsize by a factor of one quarter.
!
670 K = 1
R = 0.25D0
H = R*H
IF (ABS(H) .GE. HMIN) GO TO 690
IDID = -6
GO TO 675
!
!
!
!
!     For all crashes, restore Y to its last value,
!     interpolate to find YPRIME at last X, and return.
!
!     Before returning, verify that the user has not set
!     IDID to a nonnegative value.  If the user has set IDID
!     to a nonnegative value, then reset IDID to be -7, indicating
!     a failure in the nonlinear system solver.
!
675 CONTINUE
CALL DDATRP(X,X,Y,YPRIME,NEQ,K,PHI,PSI)
JSTART = 1
IF (IDID .GE. 0) IDID = -7
RETURN
!
!
!     Go back and try this step again.  
!     If this is the first step, reset PSI(1) and rescale PHI(*,2).
!
690 IF (KOLD .EQ. 0) THEN
  PSI(1) = H
  DO 695 I = 1,NEQ
    PHI(I,2) = R*PHI(I,2)
695 CONTINUE
ENDIF
GO TO 200
!
!------END OF SUBROUTINE DDSTP------------------------------------------
END
! Work performed under the auspices of the U.S. Department of Energy
! by Lawrence Livermore National Laboratory under contract number
! W-7405-Eng-48.
!
SUBROUTINE DCNSTR (NEQ, Y, YNEW, ICNSTR, TAU, RLX, IRET, IVAR)
!
!***BEGIN PROLOGUE  DCNSTR
!***DATE WRITTEN   950808   (YYMMDD)
!***REVISION DATE  950814   (YYMMDD)
!
!
!-----------------------------------------------------------------------
!***DESCRIPTION
!
! This subroutine checks for constraint violations in the proposed 
! new approximate solution YNEW.
! If a constraint violation occurs, then a new step length, TAU,
! is calculated, and this value is to be given to the linesearch routine
! to calculate a new approximate solution YNEW.
!
! On entry:
!
!   NEQ    -- size of the nonlinear system, and the length of arrays
!             Y, YNEW and ICNSTR.
!
!   Y      -- real array containing the current approximate y.
!
!   YNEW   -- real array containing the new approximate y.
!
!   ICNSTR -- INTEGER array of length NEQ containing flags indicating
!             which entries in YNEW are to be constrained.
!             if ICNSTR(I) =  2, then YNEW(I) must be .GT. 0,
!             if ICNSTR(I) =  1, then YNEW(I) must be .GE. 0,
!             if ICNSTR(I) = -1, then YNEW(I) must be .LE. 0, while
!             if ICNSTR(I) = -2, then YNEW(I) must be .LT. 0, while
!             if ICNSTR(I) =  0, then YNEW(I) is not constrained.
!
!   RLX    -- real scalar restricting update, if ICNSTR(I) = 2 or -2,
!             to ABS( (YNEW-Y)/Y ) < FAC2*RLX in component I.
!
!   TAU    -- the current size of the step length for the linesearch.
!
! On return
!
!   TAU    -- the adjusted size of the step length if a constraint
!             violation occurred (otherwise, it is unchanged).  it is
!             the step length to give to the linesearch routine.
!
!   IRET   -- output flag.
!             IRET=0 means that YNEW satisfied all constraints.
!             IRET=1 means that YNEW failed to satisfy all the
!                    constraints, and a new linesearch step
!                    must be computed.
!
!   IVAR   -- index of variable causing constraint to be violated.
!
!-----------------------------------------------------------------------
IMPLICIT DOUBLE PRECISION(A-H,O-Z)
DIMENSION Y(NEQ), YNEW(NEQ), ICNSTR(NEQ)
SAVE FAC, FAC2, ZERO
DATA FAC /0.6D0/, FAC2 /0.9D0/, ZERO/0.0D0/
!-----------------------------------------------------------------------
! Check constraints for proposed new step YNEW.  If a constraint has
! been violated, then calculate a new step length, TAU, to be
! used in the linesearch routine.
!-----------------------------------------------------------------------
IRET = 0
RDYMX = ZERO
IVAR = 0
DO 100 I = 1,NEQ
!
   IF (ICNSTR(I) .EQ. 2) THEN
      RDY = ABS( (YNEW(I)-Y(I))/Y(I) )
      IF (RDY .GT. RDYMX) THEN
         RDYMX = RDY
         IVAR = I
      ENDIF
      IF (YNEW(I) .LE. ZERO) THEN
         TAU = FAC*TAU
         IVAR = I
         IRET = 1
         RETURN
      ENDIF
!
   ELSEIF (ICNSTR(I) .EQ. 1) THEN
      IF (YNEW(I) .LT. ZERO) THEN
         TAU = FAC*TAU
         IVAR = I
         IRET = 1
         RETURN
      ENDIF
!
   ELSEIF (ICNSTR(I) .EQ. -1) THEN
      IF (YNEW(I) .GT. ZERO) THEN
         TAU = FAC*TAU
         IVAR = I
         IRET = 1
         RETURN
      ENDIF
!
   ELSEIF (ICNSTR(I) .EQ. -2) THEN
      RDY = ABS( (YNEW(I)-Y(I))/Y(I) )
      IF (RDY .GT. RDYMX) THEN
         RDYMX = RDY
         IVAR = I
      ENDIF
      IF (YNEW(I) .GE. ZERO) THEN
         TAU = FAC*TAU
         IVAR = I
         IRET = 1
         RETURN
      ENDIF
!
   ENDIF
100 CONTINUE

IF(RDYMX .GE. RLX) THEN
   TAU = FAC2*TAU*RLX/RDYMX
   IRET = 1
ENDIF
!
RETURN
!----------------------- END OF SUBROUTINE DCNSTR ----------------------
END
! Work performed under the auspices of the U.S. Department of Energy
! by Lawrence Livermore National Laboratory under contract number
! W-7405-Eng-48.
!
SUBROUTINE DCNST0 (NEQ, Y, ICNSTR, IRET)
!
!***BEGIN PROLOGUE  DCNST0
!***DATE WRITTEN   950808   (YYMMDD)
!***REVISION DATE  950808   (YYMMDD)
!
!
!-----------------------------------------------------------------------
!***DESCRIPTION
!
! This subroutine checks for constraint violations in the initial 
! approximate solution u.
!
! On entry
!
!   NEQ    -- size of the nonlinear system, and the length of arrays
!             Y and ICNSTR.
!
!   Y      -- real array containing the initial approximate root.
!
!   ICNSTR -- INTEGER array of length NEQ containing flags indicating
!             which entries in Y are to be constrained.
!             if ICNSTR(I) =  2, then Y(I) must be .GT. 0,
!             if ICNSTR(I) =  1, then Y(I) must be .GE. 0,
!             if ICNSTR(I) = -1, then Y(I) must be .LE. 0, while
!             if ICNSTR(I) = -2, then Y(I) must be .LT. 0, while
!             if ICNSTR(I) =  0, then Y(I) is not constrained.
!
! On return
!
!   IRET   -- output flag.
!             IRET=0    means that u satisfied all constraints.
!             IRET.NE.0 means that Y(IRET) failed to satisfy its
!                       constraint.
!
!-----------------------------------------------------------------------
IMPLICIT DOUBLE PRECISION(A-H,O-Z)
DIMENSION Y(NEQ), ICNSTR(NEQ)
SAVE ZERO
DATA ZERO/0.D0/
!-----------------------------------------------------------------------
! Check constraints for initial Y.  If a constraint has been violated,
! set IRET = I to signal an error return to calling routine.
!-----------------------------------------------------------------------
IRET = 0
DO 100 I = 1,NEQ
   IF (ICNSTR(I) .EQ. 2) THEN
      IF (Y(I) .LE. ZERO) THEN
         IRET = I
         RETURN
      ENDIF
   ELSEIF (ICNSTR(I) .EQ. 1) THEN
      IF (Y(I) .LT. ZERO) THEN
         IRET = I
         RETURN
      ENDIF
   ELSEIF (ICNSTR(I) .EQ. -1) THEN
      IF (Y(I) .GT. ZERO) THEN
         IRET = I
         RETURN
      ENDIF
   ELSEIF (ICNSTR(I) .EQ. -2) THEN
      IF (Y(I) .GE. ZERO) THEN
         IRET = I
         RETURN
      ENDIF
  ENDIF
100 CONTINUE
RETURN
!----------------------- END OF SUBROUTINE DCNST0 ----------------------
END
! Work performed under the auspices of the U.S. Department of Energy
! by Lawrence Livermore National Laboratory under contract number 
! W-7405-Eng-48.
!
SUBROUTINE DDAWTS(NEQ,IWT,RTOL,ATOL,Y,WT,RPAR,IPAR)
!
!***BEGIN PROLOGUE  DDAWTS
!***REFER TO  DDASPK
!***ROUTINES CALLED  (NONE)
!***DATE WRITTEN   890101   (YYMMDD)
!***REVISION DATE  900926   (YYMMDD)
!***END PROLOGUE  DDAWTS
!-----------------------------------------------------------------------
!     This subroutine sets the error weight vector,
!     WT, according to WT(I)=RTOL(I)*ABS(Y(I))+ATOL(I),
!     I = 1 to NEQ.
!     RTOL and ATOL are scalars if IWT = 0,
!     and vectors if IWT = 1.
!-----------------------------------------------------------------------
!
IMPLICIT DOUBLE PRECISION(A-H,O-Z)
DIMENSION RTOL(*),ATOL(*),Y(*),WT(*)
DIMENSION RPAR(*),IPAR(*)
RTOLI=RTOL(1)
ATOLI=ATOL(1)
DO 20 I=1,NEQ
   IF (IWT .EQ.0) GO TO 10
     RTOLI=RTOL(I)
     ATOLI=ATOL(I)
10 WT(I)=RTOLI*ABS(Y(I))+ATOLI
20 CONTINUE
RETURN
!
!------END OF SUBROUTINE DDAWTS-----------------------------------------
END
! Work performed under the auspices of the U.S. Department of Energy
! by Lawrence Livermore National Laboratory under contract number 
! W-7405-Eng-48.
!
SUBROUTINE DINVWT(NEQ,WT,IER)
!
!***BEGIN PROLOGUE  DINVWT
!***REFER TO  DDASPK
!***ROUTINES CALLED  (NONE)
!***DATE WRITTEN   950125   (YYMMDD)
!***END PROLOGUE  DINVWT
!-----------------------------------------------------------------------
!     This subroutine checks the error weight vector WT, of length NEQ,
!     for components that are .le. 0, and if none are found, it
!     inverts the WT(I) in place.  This replaces division operations
!     with multiplications in all norm evaluations.
!     IER is returned as 0 if all WT(I) were found positive,
!     and the first I with WT(I) .le. 0.0 otherwise.
!-----------------------------------------------------------------------
!
IMPLICIT DOUBLE PRECISION(A-H,O-Z)
DIMENSION WT(*)
!
DO 10 I = 1,NEQ
  IF (WT(I) .LE. 0.0D0) GO TO 30
10 CONTINUE
DO 20 I = 1,NEQ
  WT(I) = 1.0D0/WT(I)
20 CONTINUE
IER = 0
RETURN
!
30 IER = I
RETURN
!
!------END OF SUBROUTINE DINVWT-----------------------------------------
END
! Work performed under the auspices of the U.S. Department of Energy
! by Lawrence Livermore National Laboratory under contract number 
! W-7405-Eng-48.
!
SUBROUTINE DDATRP(X,XOUT,YOUT,YPOUT,NEQ,KOLD,PHI,PSI)
!
!***BEGIN PROLOGUE  DDATRP
!***REFER TO  DDASPK
!***ROUTINES CALLED  (NONE)
!***DATE WRITTEN   890101   (YYMMDD)
!***REVISION DATE  900926   (YYMMDD)
!***END PROLOGUE  DDATRP
!
!-----------------------------------------------------------------------
!     The methods in subroutine DDSTP use polynomials
!     to approximate the solution.  DDATRP approximates the
!     solution and its derivative at time XOUT by evaluating
!     one of these polynomials, and its derivative, there.
!     Information defining this polynomial is passed from
!     DDSTP, so DDATRP cannot be used alone.
!
!     The parameters are
!
!     X     The current time in the integration.
!     XOUT  The time at which the solution is desired.
!     YOUT  The interpolated approximation to Y at XOUT.
!           (This is output.)
!     YPOUT The interpolated approximation to YPRIME at XOUT.
!           (This is output.)
!     NEQ   Number of equations.
!     KOLD  Order used on last successful step.
!     PHI   Array of scaled divided differences of Y.
!     PSI   Array of past stepsize history.
!-----------------------------------------------------------------------
!
IMPLICIT DOUBLE PRECISION(A-H,O-Z)
DIMENSION YOUT(*),YPOUT(*)
DIMENSION PHI(NEQ,*),PSI(*)
KOLDP1=KOLD+1
TEMP1=XOUT-X
DO 10 I=1,NEQ
   YOUT(I)=PHI(I,1)
   YPOUT(I)=0.0D0
10 CONTINUE
C=1.0D0
D=0.0D0
GAMMA=TEMP1/PSI(1)
DO 30 J=2,KOLDP1
   D=D*GAMMA+C/PSI(J-1)
   C=C*GAMMA
   GAMMA=(TEMP1+PSI(J-1))/PSI(J)
   DO 20 I=1,NEQ
      YOUT(I)=YOUT(I)+C*PHI(I,J)
      YPOUT(I)=YPOUT(I)+D*PHI(I,J)
20 CONTINUE
30 CONTINUE
RETURN
!
!------END OF SUBROUTINE DDATRP-----------------------------------------
END
! Work performed under the auspices of the U.S. Department of Energy
! by Lawrence Livermore National Laboratory under contract number 
! W-7405-Eng-48.
!
DOUBLE PRECISION FUNCTION DDWNRM(NEQ,V,RWT,RPAR,IPAR)
!
!***BEGIN PROLOGUE  DDWNRM
!***ROUTINES CALLED  (NONE)
!***DATE WRITTEN   890101   (YYMMDD)
!***REVISION DATE  900926   (YYMMDD)
!***END PROLOGUE  DDWNRM
!-----------------------------------------------------------------------
!     This function routine computes the weighted
!     root-mean-square norm of the vector of length
!     NEQ contained in the array V, with reciprocal weights
!     contained in the array RWT of length NEQ.
!        DDWNRM=SQRT((1/NEQ)*SUM(V(I)*RWT(I))**2)
!-----------------------------------------------------------------------
!
IMPLICIT DOUBLE PRECISION(A-H,O-Z)
DIMENSION V(*),RWT(*)
DIMENSION RPAR(*),IPAR(*)
DDWNRM = 0.0D0
VMAX = 0.0D0
DO 10 I = 1,NEQ
  IF(ABS(V(I)*RWT(I)) .GT. VMAX) VMAX = ABS(V(I)*RWT(I))
10 CONTINUE
IF(VMAX .LE. 0.0D0) GO TO 30
SUM = 0.0D0
DO 20 I = 1,NEQ
  SUM = SUM + ((V(I)*RWT(I))/VMAX)**2
20 CONTINUE
DDWNRM = VMAX*SQRT(SUM/NEQ)
30 CONTINUE
RETURN
!
!------END OF FUNCTION DDWNRM-------------------------------------------
END
! Work performed under the auspices of the U.S. Department of Energy
! by Lawrence Livermore National Laboratory under contract number 
! W-7405-Eng-48.
!
SUBROUTINE DDASID(X,Y,YPRIME,NEQ,ICOPT,ID,RES,JACD,PDUM,H,WT, &
& JSDUM,RPAR,IPAR,DUMSVR,DELTA,R,YIC,YPIC,DUMPWK,WM,IWM,CJ,UROUND, &
& DUME,DUMS,DUMR,EPCON,RATEMX,STPTOL,JFDUM, &
& ICNFLG,ICNSTR,IERNLS)
!
!***BEGIN PROLOGUE  DDASID
!***REFER TO  DDASPK
!***DATE WRITTEN   940701   (YYMMDD)
!***REVISION DATE  950808   (YYMMDD)
!***REVISION DATE  951110   Removed unreachable block 390.
!
!
!-----------------------------------------------------------------------
!***DESCRIPTION
!
!
!     DDASID solves a nonlinear system of algebraic equations of the
!     form G(X,Y,YPRIME) = 0 for the unknown parts of Y and YPRIME in
!     the initial conditions.
!
!     The method used is a modified Newton scheme.
!
!     The parameters represent
!
!     X         -- Independent variable.
!     Y         -- Solution vector.
!     YPRIME    -- Derivative of solution vector.
!     NEQ       -- Number of unknowns.
!     ICOPT     -- Initial condition option chosen (1 or 2).
!     ID        -- Array of dimension NEQ, which must be initialized
!                  if ICOPT = 1.  See DDASIC.
!     RES       -- External user-supplied subroutine to evaluate the
!                  residual.  See RES description in DDASPK prologue.
!     JACD      -- External user-supplied routine to evaluate the
!                  Jacobian.  See JAC description for the case
!                  INFO(12) = 0 in the DDASPK prologue.
!     PDUM      -- Dummy argument.
!     H         -- Scaling factor for this initial condition calc.
!     WT        -- Vector of weights for error criterion.
!     JSDUM     -- Dummy argument.
!     RPAR,IPAR -- Real and integer arrays used for communication
!                  between the calling program and external user
!                  routines.  They are not altered within DASPK.
!     DUMSVR    -- Dummy argument.
!     DELTA     -- Work vector for NLS of length NEQ.
!     R         -- Work vector for NLS of length NEQ.
!     YIC,YPIC  -- Work vectors for NLS, each of length NEQ.
!     DUMPWK    -- Dummy argument.
!     WM,IWM    -- Real and integer arrays storing matrix information
!                  such as the matrix of partial derivatives,
!                  permutation vector, and various other information.
!     CJ        -- Matrix parameter = 1/H (ICOPT = 1) or 0 (ICOPT = 2).
!     UROUND    -- Unit roundoff.
!     DUME      -- Dummy argument.
!     DUMS      -- Dummy argument.
!     DUMR      -- Dummy argument.
!     EPCON     -- Tolerance to test for convergence of the Newton
!                  iteration.
!     RATEMX    -- Maximum convergence rate for which Newton iteration
!                  is considered converging.
!     JFDUM     -- Dummy argument.
!     STPTOL    -- Tolerance used in calculating the minimum lambda
!                  value allowed.
!     ICNFLG    -- Integer scalar.  If nonzero, then constraint
!                  violations in the proposed new approximate solution
!                  will be checked for, and the maximum step length 
!                  will be adjusted accordingly.
!     ICNSTR    -- Integer array of length NEQ containing flags for
!                  checking constraints.
!     IERNLS    -- Error flag for nonlinear solver.
!                   0   ==> nonlinear solver converged.
!                   1,2 ==> recoverable error inside nonlinear solver.
!                           1 => retry with current Y, YPRIME
!                           2 => retry with original Y, YPRIME
!                  -1   ==> unrecoverable error in nonlinear solver.
!
!     All variables with "DUM" in their names are dummy variables
!     which are not used in this routine.
!
!-----------------------------------------------------------------------
!
!***ROUTINES CALLED
!   RES, DMATD, DNSID
!
!***END PROLOGUE  DDASID
!
!
IMPLICIT DOUBLE PRECISION(A-H,O-Z)
DIMENSION Y(*),YPRIME(*),ID(*),WT(*),ICNSTR(*)
DIMENSION DELTA(*),R(*),YIC(*),YPIC(*)
DIMENSION WM(*),IWM(*), RPAR(*),IPAR(*)
EXTERNAL  RES, JACD
!
INTEGER, PARAMETER :: LNRE=12, LNJE=13, LMXNIT=32, LMXNJ=33
!
!
!     Perform initializations.
!
MXNIT = IWM(LMXNIT)
MXNJ = IWM(LMXNJ)
IERNLS = 0
NJ = 0
!
!     Call RES to initialize DELTA.
!
IRES = 0
IWM(LNRE) = IWM(LNRE) + 1
CALL RES(X,Y,YPRIME,CJ,DELTA,IRES,RPAR,IPAR)
IF (IRES .LT. 0) GO TO 370
!
!     Looping point for updating the Jacobian.
!
300 CONTINUE
!
!     Initialize all error flags to zero.
!
IERJ = 0
IRES = 0
IERNEW = 0
!
!     Reevaluate the iteration matrix, J = dG/dY + CJ*dG/dYPRIME,
!     where G(X,Y,YPRIME) = 0.
!
NJ = NJ + 1
IWM(LNJE)=IWM(LNJE)+1
CALL DMATD(NEQ,X,Y,YPRIME,DELTA,CJ,H,IERJ,WT,R, &
& WM,IWM,RES,IRES,UROUND,JACD,RPAR,IPAR)
IF (IRES .LT. 0 .OR. IERJ .NE. 0) GO TO 370
!
!     Call the nonlinear Newton solver for up to MXNIT iterations.
!
CALL DNSID(X,Y,YPRIME,NEQ,ICOPT,ID,RES,WT,RPAR,IPAR,DELTA,R, &
& YIC,YPIC,WM,IWM,CJ,EPCON,RATEMX,MXNIT,STPTOL, &
& ICNFLG,ICNSTR,IERNEW)
!
IF (IERNEW .EQ. 1 .AND. NJ .LT. MXNJ) THEN
!
!        MXNIT iterations were done, the convergence rate is < 1,
!        and the number of Jacobian evaluations is less than MXNJ.
!        Call RES, reevaluate the Jacobian, and try again.
!
   IWM(LNRE)=IWM(LNRE)+1
   CALL RES(X,Y,YPRIME,CJ,DELTA,IRES,RPAR,IPAR)
   IF (IRES .LT. 0) GO TO 370
   GO TO 300
   ENDIF
!
IF (IERNEW .NE. 0) GO TO 380

RETURN
!
!
!     Unsuccessful exits from nonlinear solver.
!     Compute IERNLS accordingly.
!
370 IERNLS = 2
IF (IRES .LE. -2) IERNLS = -1
RETURN
!
380 IERNLS = MIN(IERNEW,2)
RETURN
!
!------END OF SUBROUTINE DDASID-----------------------------------------
END
! Work performed under the auspices of the U.S. Department of Energy
! by Lawrence Livermore National Laboratory under contract number 
! W-7405-Eng-48.
!
SUBROUTINE DNSID(X,Y,YPRIME,NEQ,ICOPT,ID,RES,WT,RPAR,IPAR, &
& DELTA,R,YIC,YPIC,WM,IWM,CJ,EPCON,RATEMX,MAXIT,STPTOL, &
& ICNFLG,ICNSTR,IERNEW)
!
!***BEGIN PROLOGUE  DNSID
!***REFER TO  DDASPK
!***DATE WRITTEN   940701   (YYMMDD)
!***REVISION DATE  950713   (YYMMDD)
!
!
!-----------------------------------------------------------------------
!***DESCRIPTION
!
!     DNSID solves a nonlinear system of algebraic equations of the
!     form G(X,Y,YPRIME) = 0 for the unknown parts of Y and YPRIME
!     in the initial conditions.
!
!     The method used is a modified Newton scheme.
!
!     The parameters represent
!
!     X         -- Independent variable.
!     Y         -- Solution vector.
!     YPRIME    -- Derivative of solution vector.
!     NEQ       -- Number of unknowns.
!     ICOPT     -- Initial condition option chosen (1 or 2).
!     ID        -- Array of dimension NEQ, which must be initialized
!                  if ICOPT = 1.  See DDASIC.
!     RES       -- External user-supplied subroutine to evaluate the
!                  residual.  See RES description in DDASPK prologue.
!     WT        -- Vector of weights for error criterion.
!     RPAR,IPAR -- Real and integer arrays used for communication
!                  between the calling program and external user
!                  routines.  They are not altered within DASPK.
!     DELTA     -- Residual vector on entry, and work vector of
!                  length NEQ for DNSID.
!     WM,IWM    -- Real and integer arrays storing matrix information
!                  such as the matrix of partial derivatives,
!                  permutation vector, and various other information.
!     CJ        -- Matrix parameter = 1/H (ICOPT = 1) or 0 (ICOPT = 2).
!     R         -- Array of length NEQ used as workspace by the 
!                  linesearch routine DLINSD.
!     YIC,YPIC  -- Work vectors for DLINSD, each of length NEQ.
!     EPCON     -- Tolerance to test for convergence of the Newton
!                  iteration.
!     RATEMX    -- Maximum convergence rate for which Newton iteration
!                  is considered converging.
!     MAXIT     -- Maximum allowed number of Newton iterations.
!     STPTOL    -- Tolerance used in calculating the minimum lambda
!                  value allowed.
!     ICNFLG    -- Integer scalar.  If nonzero, then constraint
!                  violations in the proposed new approximate solution
!                  will be checked for, and the maximum step length 
!                  will be adjusted accordingly.
!     ICNSTR    -- Integer array of length NEQ containing flags for
!                  checking constraints.
!     IERNEW    -- Error flag for Newton iteration.
!                   0  ==> Newton iteration converged.
!                   1  ==> failed to converge, but RATE .le. RATEMX.
!                   2  ==> failed to converge, RATE .gt. RATEMX.
!                   3  ==> other recoverable error (IRES = -1, or
!                          linesearch failed).
!                  -1  ==> unrecoverable error (IRES = -2).
!
!-----------------------------------------------------------------------
!
!***ROUTINES CALLED
!   DSLVD, DDWNRM, DLINSD, DCOPY
!
!***END PROLOGUE  DNSID
!
!
IMPLICIT DOUBLE PRECISION(A-H,O-Z)
DIMENSION Y(*),YPRIME(*),WT(*),R(*)
DIMENSION ID(*),DELTA(*), YIC(*), YPIC(*)
DIMENSION WM(*),IWM(*), RPAR(*),IPAR(*)
DIMENSION ICNSTR(*)
EXTERNAL  RES
!
INTEGER, PARAMETER :: LNNI=19, LLSOFF=35
!
!
!     Initializations.  M is the Newton iteration counter.
!
LSOFF = IWM(LLSOFF)
M = 0
RATE = 1.0D0
RLX = 0.4D0
!
!     Compute a new step vector DELTA by back-substitution.
!
CALL DSLVD (NEQ, DELTA, WM, IWM)
!
!     Get norm of DELTA.  Return now if norm(DELTA) .le. EPCON.
!
DELNRM = DDWNRM(NEQ,DELTA,WT,RPAR,IPAR)
FNRM = DELNRM
IF (FNRM .LE. EPCON) RETURN
!
!     Newton iteration loop.
!
300 CONTINUE
IWM(LNNI) = IWM(LNNI) + 1
!
!     Call linesearch routine for global strategy and set RATE
!
OLDFNM = FNRM
!
CALL DLINSD (NEQ, Y, X, YPRIME, CJ, DELTA, DELNRM, WT, LSOFF, &
& STPTOL, IRET, RES, IRES, WM, IWM, FNRM, ICOPT, ID, &
& R, YIC, YPIC, ICNFLG, ICNSTR, RLX, RPAR, IPAR)
!
RATE = FNRM/OLDFNM
!
!     Check for error condition from linesearch.
IF (IRET .NE. 0) GO TO 390
!
!     Test for convergence of the iteration, and return or loop.
!
IF (FNRM .LE. EPCON) RETURN
!
!     The iteration has not yet converged.  Update M.
!     Test whether the maximum number of iterations have been tried.
!
M = M + 1
IF (M .GE. MAXIT) GO TO 380
!
!     Copy the residual to DELTA and its norm to DELNRM, and loop for
!     another iteration.
!
CALL DCOPY (NEQ, R, 1, DELTA, 1)
DELNRM = FNRM
GO TO 300
!
!     The maximum number of iterations was done.  Set IERNEW and return.
!
380 IF (RATE .LE. RATEMX) THEN
   IERNEW = 1
ELSE
   IERNEW = 2
ENDIF
RETURN
!
390 IF (IRES .LE. -2) THEN
   IERNEW = -1
ELSE
   IERNEW = 3
ENDIF
RETURN
!
!
!------END OF SUBROUTINE DNSID------------------------------------------
END
! Work performed under the auspices of the U.S. Department of Energy
! by Lawrence Livermore National Laboratory under contract number 
! W-7405-Eng-48.
!
SUBROUTINE DLINSD (NEQ, Y, T, YPRIME, CJ, P, PNRM, WT, LSOFF, &
& STPTOL, IRET, RES, IRES, WM, IWM, &
& FNRM, ICOPT, ID, R, YNEW, YPNEW, ICNFLG, &
& ICNSTR, RLX, RPAR, IPAR)
!
!***BEGIN PROLOGUE  DLINSD
!***REFER TO  DNSID
!***DATE WRITTEN   941025   (YYMMDD)
!***REVISION DATE  941215   (YYMMDD)
!***REVISION DATE  960129   Moved line RL = ONE to top block.
!
!
!-----------------------------------------------------------------------
!***DESCRIPTION
!
!     DLINSD uses a linesearch algorithm to calculate a new (Y,YPRIME)
!     pair (YNEW,YPNEW) such that 
!
!     f(YNEW,YPNEW) .le. (1 - 2*ALPHA*RL)*f(Y,YPRIME) ,
!
!     where 0 < RL <= 1.  Here, f(y,y') is defined as
!
!      f(y,y') = (1/2)*norm( (J-inverse)*G(t,y,y') )**2 ,
!
!     where norm() is the weighted RMS vector norm, G is the DAE
!     system residual function, and J is the system iteration matrix
!     (Jacobian).
!
!     In addition to the parameters defined elsewhere, we have
!
!     P       -- Approximate Newton step used in backtracking.
!     PNRM    -- Weighted RMS norm of P.
!     LSOFF   -- Flag showing whether the linesearch algorithm is
!                to be invoked.  0 means do the linesearch, and
!                1 means turn off linesearch.
!     STPTOL  -- Tolerance used in calculating the minimum lambda
!                value allowed.
!     ICNFLG  -- Integer scalar.  If nonzero, then constraint violations
!                in the proposed new approximate solution will be
!                checked for, and the maximum step length will be
!                adjusted accordingly.
!     ICNSTR  -- Integer array of length NEQ containing flags for
!                checking constraints.
!     RLX     -- Real scalar restricting update size in DCNSTR.
!     YNEW    -- Array of length NEQ used to hold the new Y in
!                performing the linesearch.
!     YPNEW   -- Array of length NEQ used to hold the new YPRIME in
!                performing the linesearch.
!     Y       -- Array of length NEQ containing the new Y (i.e.,=YNEW).
!     YPRIME  -- Array of length NEQ containing the new YPRIME 
!                (i.e.,=YPNEW).
!     FNRM    -- Real scalar containing SQRT(2*f(Y,YPRIME)) for the
!                current (Y,YPRIME) on input and output.
!     R       -- Work array of length NEQ, containing the scaled 
!                residual (J-inverse)*G(t,y,y') on return.
!     IRET    -- Return flag.
!                IRET=0 means that a satisfactory (Y,YPRIME) was found.
!                IRET=1 means that the routine failed to find a new
!                       (Y,YPRIME) that was sufficiently distinct from
!                       the current (Y,YPRIME) pair.
!                IRET=2 means IRES .ne. 0 from RES.
!-----------------------------------------------------------------------
!
!***ROUTINES CALLED
!   DFNRMD, DYYPNW, DCOPY
!
!***END PROLOGUE  DLINSD
!
IMPLICIT DOUBLE PRECISION(A-H,O-Z)
EXTERNAL  RES
DIMENSION Y(*), YPRIME(*), WT(*), R(*), ID(*)
DIMENSION WM(*), IWM(*)
DIMENSION YNEW(*), YPNEW(*), P(*), ICNSTR(*)
DIMENSION RPAR(*), IPAR(*)
!
INTEGER, PARAMETER :: LNRE=12, LKPRIN=31
!
SAVE ALPHA, ONE, TWO
DATA ALPHA/1.0D-4/, ONE/1.0D0/, TWO/2.0D0/
!
KPRIN=IWM(LKPRIN)
!
F1NRM = (FNRM*FNRM)/TWO
RATIO = ONE
IF (KPRIN .GE. 2) THEN
call rprintfd1( &
& 'daspk--  in routine dlinsd--PNRM (= %g)' // char(0),  PNRM)

  ENDIF
TAU = PNRM
IVIO = 0
RL = ONE
!-----------------------------------------------------------------------
! Check for violations of the constraints, if any are imposed.
! If any violations are found, the step vector P is rescaled, and the 
! constraint check is repeated, until no violations are found.
!-----------------------------------------------------------------------
IF (ICNFLG .NE. 0) THEN
10 CONTINUE
   CALL DYYPNW (NEQ,Y,YPRIME,CJ,RL,P,ICOPT,ID,YNEW,YPNEW)
   CALL DCNSTR (NEQ, Y, YNEW, ICNSTR, TAU, RLX, IRET, IVAR)
   IF (IRET .EQ. 1) THEN
      IVIO = 1
      RATIO1 = TAU/PNRM
      RATIO = RATIO*RATIO1
      DO 20 I = 1,NEQ
        P(I) = P(I)*RATIO1
20 CONTINUE
      PNRM = TAU
      IF (KPRIN .GE. 2) THEN
call rprintfid( &
& 'daspk--  constraint violation-PNRM (= %g), index =( %i)' &
& // char(0), &
& IVAR,PNRM)
        ENDIF
      IF (PNRM .LE. STPTOL) THEN
        IRET = 1
        RETURN
        ENDIF
      GO TO 10
      ENDIF
   ENDIF
!
SLPI = (-TWO*F1NRM)*RATIO
RLMIN = STPTOL/PNRM
IF (LSOFF .EQ. 0 .AND. KPRIN .GE. 2) THEN
call rprintfd1( &
& 'daspk--  min lambda (= %g)' // char(0), RLMIN)

  ENDIF
!-----------------------------------------------------------------------
! Begin iteration to find RL value satisfying alpha-condition.
! If RL becomes less than RLMIN, then terminate with IRET = 1.
!-----------------------------------------------------------------------
100 CONTINUE
CALL DYYPNW (NEQ,Y,YPRIME,CJ,RL,P,ICOPT,ID,YNEW,YPNEW)
CALL DFNRMD (NEQ, YNEW, T, YPNEW, R, CJ, WT, RES, IRES, &
& FNRMP, WM, IWM, RPAR, IPAR)
IWM(LNRE) = IWM(LNRE) + 1
IF (IRES .NE. 0) THEN
  IRET = 2
  RETURN
  ENDIF
IF (LSOFF .EQ. 1) GO TO 150
!
F1NRMP = FNRMP*FNRMP/TWO
IF (KPRIN .GE. 2) THEN
call rprintfd1( &
& 'daspk--  LAMBDA (= %g)' // char(0), RL)
call rprintfd2( &
& 'daspk--  NORM(F1) = %g, NORM(F1NEW) = %g' &
& // char(0), F1NRM, F1NRMP)

  ENDIF
IF (F1NRMP .GT. F1NRM + ALPHA*SLPI*RL) GO TO 200
!-----------------------------------------------------------------------
! Alpha-condition is satisfied, or linesearch is turned off.
! Copy YNEW,YPNEW to Y,YPRIME and return.
!-----------------------------------------------------------------------
150 IRET = 0
CALL DCOPY (NEQ, YNEW, 1, Y, 1)
CALL DCOPY (NEQ, YPNEW, 1, YPRIME, 1)
FNRM = FNRMP
IF (KPRIN .GE. 1) THEN
call rprintfd1( &
& 'daspk--  leaving routine dlinsd--FNRM (= %g)' &
& // char(0),FNRM)

  ENDIF
RETURN
!-----------------------------------------------------------------------
! Alpha-condition not satisfied.  Perform backtrack to compute new RL
! value.  If no satisfactory YNEW,YPNEW can be found sufficiently 
! distinct from Y,YPRIME, then return IRET = 1.
!-----------------------------------------------------------------------
200 CONTINUE
IF (RL .LT. RLMIN) THEN
  IRET = 1
  RETURN
  ENDIF
!
RL = RL/TWO
GO TO 100
!
!----------------------- END OF SUBROUTINE DLINSD ----------------------
END
! Work performed under the auspices of the U.S. Department of Energy
! by Lawrence Livermore National Laboratory under contract number 
! W-7405-Eng-48.
!
SUBROUTINE DFNRMD (NEQ, Y, T, YPRIME, R, CJ, WT, RES, IRES, &
& FNORM, WM, IWM, RPAR, IPAR)
!
!***BEGIN PROLOGUE  DFNRMD
!***REFER TO  DLINSD
!***DATE WRITTEN   941025   (YYMMDD)
!
!
!-----------------------------------------------------------------------
!***DESCRIPTION
!
!     DFNRMD calculates the scaled preconditioned norm of the nonlinear
!     function used in the nonlinear iteration for obtaining consistent
!     initial conditions.  Specifically, DFNRMD calculates the weighted
!     root-mean-square norm of the vector (J-inverse)*G(T,Y,YPRIME),
!     where J is the Jacobian matrix.
!
!     In addition to the parameters described in the calling program
!     DLINSD, the parameters represent
!
!     R      -- Array of length NEQ that contains
!               (J-inverse)*G(T,Y,YPRIME) on return.
!     FNORM  -- Scalar containing the weighted norm of R on return.
!-----------------------------------------------------------------------
!
!***ROUTINES CALLED
!   RES, DSLVD, DDWNRM
!
!***END PROLOGUE  DFNRMD
!
!
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
EXTERNAL RES
DIMENSION Y(*), YPRIME(*), WT(*), R(*)
DIMENSION WM(*),IWM(*), RPAR(*),IPAR(*)
!-----------------------------------------------------------------------
!     Call RES routine.
!-----------------------------------------------------------------------
IRES = 0
CALL RES(T,Y,YPRIME,CJ,R,IRES,RPAR,IPAR)
IF (IRES .LT. 0) RETURN
!-----------------------------------------------------------------------
!     Apply inverse of Jacobian to vector R.
!-----------------------------------------------------------------------
CALL DSLVD(NEQ,R,WM,IWM)
!-----------------------------------------------------------------------
!     Calculate norm of R.
!-----------------------------------------------------------------------
FNORM = DDWNRM(NEQ,R,WT,RPAR,IPAR)
!
RETURN
!----------------------- END OF SUBROUTINE DFNRMD ----------------------
END
! Work performed under the auspices of the U.S. Department of Energy
! by Lawrence Livermore National Laboratory under contract number 
! W-7405-Eng-48.
!
SUBROUTINE DNEDD(X,Y,YPRIME,NEQ,RES,JACD,PDUM,H,WT, &
& JSTART,IDID,RPAR,IPAR,PHI,GAMMA,DUMSVR,DELTA,E, &
& WM,IWM,CJ,CJOLD,CJLAST,S,UROUND,DUME,DUMS,DUMR, &
& EPCON,JCALC,JFDUM,KP1,NONNEG,NTYPE,IERNLS)
!
!***BEGIN PROLOGUE  DNEDD
!***REFER TO  DDASPK
!***DATE WRITTEN   891219   (YYMMDD)
!***REVISION DATE  900926   (YYMMDD)
!
!
!-----------------------------------------------------------------------
!***DESCRIPTION
!
!     DNEDD solves a nonlinear system of
!     algebraic equations of the form
!     G(X,Y,YPRIME) = 0 for the unknown Y.
!
!     The method used is a modified Newton scheme.
!
!     The parameters represent
!
!     X         -- Independent variable.
!     Y         -- Solution vector.
!     YPRIME    -- Derivative of solution vector.
!     NEQ       -- Number of unknowns.
!     RES       -- External user-supplied subroutine
!                  to evaluate the residual.  See RES description
!                  in DDASPK prologue.
!     JACD      -- External user-supplied routine to evaluate the
!                  Jacobian.  See JAC description for the case
!                  INFO(12) = 0 in the DDASPK prologue.
!     PDUM      -- Dummy argument.
!     H         -- Appropriate step size for next step.
!     WT        -- Vector of weights for error criterion.
!     JSTART    -- Indicates first call to this routine.
!                  If JSTART = 0, then this is the first call,
!                  otherwise it is not.
!     IDID      -- Completion flag, output by DNEDD.
!                  See IDID description in DDASPK prologue.
!     RPAR,IPAR -- Real and integer arrays used for communication
!                  between the calling program and external user
!                  routines.  They are not altered within DASPK.
!     PHI       -- Array of divided differences used by
!                  DNEDD.  The length is NEQ*(K+1),where
!                  K is the maximum order.
!     GAMMA     -- Array used to predict Y and YPRIME.  The length
!                  is MAXORD+1 where MAXORD is the maximum order.
!     DUMSVR    -- Dummy argument.
!     DELTA     -- Work vector for NLS of length NEQ.
!     E         -- Error accumulation vector for NLS of length NEQ.
!     WM,IWM    -- Real and integer arrays storing
!                  matrix information such as the matrix
!                  of partial derivatives, permutation
!                  vector, and various other information.
!     CJ        -- Parameter always proportional to 1/H.
!     CJOLD     -- Saves the value of CJ as of the last call to DMATD.
!                  Accounts for changes in CJ needed to
!                  decide whether to call DMATD.
!     CJLAST    -- Previous value of CJ.
!     S         -- A scalar determined by the approximate rate
!                  of convergence of the Newton iteration and used
!                  in the convergence test for the Newton iteration.
!
!                  If RATE is defined to be an estimate of the
!                  rate of convergence of the Newton iteration,
!                  then S = RATE/(1.D0-RATE).
!
!                  The closer RATE is to 0., the faster the Newton
!                  iteration is converging; the closer RATE is to 1.,
!                  the slower the Newton iteration is converging.
!
!                  On the first Newton iteration with an up-dated
!                  preconditioner S = 100.D0, Thus the initial
!                  RATE of convergence is approximately 1.
!
!                  S is preserved from call to call so that the rate
!                  estimate from a previous step can be applied to
!                  the current step.
!     UROUND    -- Unit roundoff.
!     DUME      -- Dummy argument.
!     DUMS      -- Dummy argument.
!     DUMR      -- Dummy argument.
!     EPCON     -- Tolerance to test for convergence of the Newton
!                  iteration.
!     JCALC     -- Flag used to determine when to update
!                  the Jacobian matrix.  In general:
!
!                  JCALC = -1 ==> Call the DMATD routine to update
!                                 the Jacobian matrix.
!                  JCALC =  0 ==> Jacobian matrix is up-to-date.
!                  JCALC =  1 ==> Jacobian matrix is out-dated,
!                                 but DMATD will not be called unless
!                                 JCALC is set to -1.
!     JFDUM     -- Dummy argument.
!     KP1       -- The current order(K) + 1;  updated across calls.
!     NONNEG    -- Flag to determine nonnegativity constraints.
!     NTYPE     -- Identification code for the NLS routine.
!                   0  ==> modified Newton; direct solver.
!     IERNLS    -- Error flag for nonlinear solver.
!                   0  ==> nonlinear solver converged.
!                   1  ==> recoverable error inside nonlinear solver.
!                  -1  ==> unrecoverable error inside nonlinear solver.
!
!     All variables with "DUM" in their names are dummy variables
!     which are not used in this routine.
!
!     Following is a list and description of local variables which
!     may not have an obvious usage.  They are listed in roughly the
!     order they occur in this subroutine.
!
!     The following group of variables are passed as arguments to
!     the Newton iteration solver.  They are explained in greater detail
!     in DNSD:
!        TOLNEW, MULDEL, MAXIT, IERNEW
!
!     IERTYP -- Flag which tells whether this subroutine is correct.
!               0 ==> correct subroutine.
!               1 ==> incorrect subroutine.
! 
!-----------------------------------------------------------------------
!***ROUTINES CALLED
!   DDWNRM, RES, DMATD, DNSD
!
!***END PROLOGUE  DNEDD
!
!
IMPLICIT DOUBLE PRECISION(A-H,O-Z)
DIMENSION Y(*),YPRIME(*),WT(*)
DIMENSION DELTA(*),E(*)
DIMENSION WM(*),IWM(*), RPAR(*),IPAR(*)
DIMENSION PHI(NEQ,*),GAMMA(*)
EXTERNAL  RES, JACD
!
INTEGER, PARAMETER :: LNRE=12, LNJE=13
!
SAVE MULDEL, MAXIT, XRATE
DATA MULDEL/1/, MAXIT/4/, XRATE/0.25D0/
!
!     Verify that this is the correct subroutine.
!
IERTYP = 0
IF (NTYPE .NE. 0) THEN
   IERTYP = 1
   GO TO 380
   ENDIF
!
!     If this is the first step, perform initializations.
!
IF (JSTART .EQ. 0) THEN
   CJOLD = CJ
   JCALC = -1
   ENDIF
!
!     Perform all other initializations.
!
IERNLS = 0
!
!     Decide whether new Jacobian is needed.
!
TEMP1 = (1.0D0 - XRATE)/(1.0D0 + XRATE)
TEMP2 = 1.0D0/TEMP1
IF (CJ/CJOLD .LT. TEMP1 .OR. CJ/CJOLD .GT. TEMP2) JCALC = -1
IF (CJ .NE. CJLAST) S = 100.D0
!
!-----------------------------------------------------------------------
!     Entry point for updating the Jacobian with current
!     stepsize.
!-----------------------------------------------------------------------
300 CONTINUE
!
!     Initialize all error flags to zero.
!
IERJ = 0
IRES = 0
IERNEW = 0
!
!     Predict the solution and derivative and compute the tolerance
!     for the Newton iteration.
!
DO 310 I=1,NEQ
   Y(I)=PHI(I,1)
   YPRIME(I)=0.0D0
310 CONTINUE
DO 330 J=2,KP1
   DO 320 I=1,NEQ
      Y(I)=Y(I)+PHI(I,J)
      YPRIME(I)=YPRIME(I)+GAMMA(J)*PHI(I,J)
320 CONTINUE
330 CONTINUE
PNORM = DDWNRM (NEQ,Y,WT,RPAR,IPAR)
TOLNEW = 100.D0*UROUND*PNORM
!     
!     Call RES to initialize DELTA.
!
IWM(LNRE)=IWM(LNRE)+1
CALL RES(X,Y,YPRIME,CJ,DELTA,IRES,RPAR,IPAR)
IF (IRES .LT. 0) GO TO 380
!
!     If indicated, reevaluate the iteration matrix 
!     J = dG/dY + CJ*dG/dYPRIME (where G(X,Y,YPRIME)=0).
!     Set JCALC to 0 as an indicator that this has been done.
!
IF(JCALC .EQ. -1) THEN
   IWM(LNJE)=IWM(LNJE)+1
   JCALC=0
   CALL DMATD(NEQ,X,Y,YPRIME,DELTA,CJ,H,IERJ,WT,E,WM,IWM, &
& RES,IRES,UROUND,JACD,RPAR,IPAR)
   CJOLD=CJ
   S = 100.D0
   IF (IRES .LT. 0) GO TO 380
   IF(IERJ .NE. 0)GO TO 380
ENDIF
!
!     Call the nonlinear Newton solver.
!
TEMP1 = 2.0D0/(1.0D0 + CJ/CJOLD)
CALL DNSD(X,Y,YPRIME,NEQ,RES,PDUM,WT,RPAR,IPAR,DUMSVR, &
& DELTA,E,WM,IWM,CJ,DUMS,DUMR,DUME,EPCON,S,TEMP1, &
& TOLNEW,MULDEL,MAXIT,IRES,IDUM,IERNEW)
!
IF (IERNEW .GT. 0 .AND. JCALC .NE. 0) THEN
!
!        The Newton iteration had a recoverable failure with an old
!        iteration matrix.  Retry the step with a new iteration matrix.
!
   JCALC = -1
   GO TO 300
ENDIF
!
IF (IERNEW .NE. 0) GO TO 380
!
!     The Newton iteration has converged.  If nonnegativity of
!     solution is required, set the solution nonnegative, if the
!     perturbation to do it is small enough.  If the change is too
!     large, then consider the corrector iteration to have failed.
!
IF(NONNEG .EQ. 0) GO TO 390
DO 377 I = 1,NEQ
   DELTA(I) = MIN(Y(I),0.0D0)
377 CONTINUE
DELNRM = DDWNRM(NEQ,DELTA,WT,RPAR,IPAR)
IF(DELNRM .GT. EPCON) GO TO 380
DO 378 I = 1,NEQ
   E(I) = E(I) - DELTA(I)
378 CONTINUE
GO TO 390
!
!
!     Exits from nonlinear solver.
!     No convergence with current iteration
!     matrix, or singular iteration matrix.
!     Compute IERNLS and IDID accordingly.
!
380 CONTINUE
IF (IRES .LE. -2 .OR. IERTYP .NE. 0) THEN
   IERNLS = -1
   IF (IRES .LE. -2) IDID = -11
   IF (IERTYP .NE. 0) IDID = -15
ELSE
   IERNLS = 1
   IF (IRES .LT. 0) IDID = -10
   IF (IERJ .NE. 0) IDID = -8
ENDIF
!
390 JCALC = 1
RETURN
!
!------END OF SUBROUTINE DNEDD------------------------------------------
END
! Work performed under the auspices of the U.S. Department of Energy
! by Lawrence Livermore National Laboratory under contract number 
! W-7405-Eng-48.
!
SUBROUTINE DNSD(X,Y,YPRIME,NEQ,RES,PDUM,WT,RPAR,IPAR, &
& DUMSVR,DELTA,E,WM,IWM,CJ,DUMS,DUMR,DUME,EPCON, &
& S,CONFAC,TOLNEW,MULDEL,MAXIT,IRES,IDUM,IERNEW)
!
!***BEGIN PROLOGUE  DNSD
!***REFER TO  DDASPK
!***DATE WRITTEN   891219   (YYMMDD)
!***REVISION DATE  900926   (YYMMDD)
!***REVISION DATE  950126   (YYMMDD)
!
!
!-----------------------------------------------------------------------
!***DESCRIPTION
!
!     DNSD solves a nonlinear system of
!     algebraic equations of the form
!     G(X,Y,YPRIME) = 0 for the unknown Y.
!
!     The method used is a modified Newton scheme.
!
!     The parameters represent
!
!     X         -- Independent variable.
!     Y         -- Solution vector.
!     YPRIME    -- Derivative of solution vector.
!     NEQ       -- Number of unknowns.
!     RES       -- External user-supplied subroutine
!                  to evaluate the residual.  See RES description
!                  in DDASPK prologue.
!     PDUM      -- Dummy argument.
!     WT        -- Vector of weights for error criterion.
!     RPAR,IPAR -- Real and integer arrays used for communication
!                  between the calling program and external user
!                  routines.  They are not altered within DASPK.
!     DUMSVR    -- Dummy argument.
!     DELTA     -- Work vector for DNSD of length NEQ.
!     E         -- Error accumulation vector for DNSD of length NEQ.
!     WM,IWM    -- Real and integer arrays storing
!                  matrix information such as the matrix
!                  of partial derivatives, permutation
!                  vector, and various other information.
!     CJ        -- Parameter always proportional to 1/H (step size).
!     DUMS      -- Dummy argument.
!     DUMR      -- Dummy argument.
!     DUME      -- Dummy argument.
!     EPCON     -- Tolerance to test for convergence of the Newton
!                  iteration.
!     S         -- Used for error convergence tests.
!                  In the Newton iteration: S = RATE/(1 - RATE),
!                  where RATE is the estimated rate of convergence
!                  of the Newton iteration.
!                  The calling routine passes the initial value
!                  of S to the Newton iteration.
!     CONFAC    -- A residual scale factor to improve convergence.
!     TOLNEW    -- Tolerance on the norm of Newton correction in
!                  alternative Newton convergence test.
!     MULDEL    -- A flag indicating whether or not to multiply
!                  DELTA by CONFAC.
!                  0  ==> do not scale DELTA by CONFAC.
!                  1  ==> scale DELTA by CONFAC.
!     MAXIT     -- Maximum allowed number of Newton iterations.
!     IRES      -- Error flag returned from RES.  See RES description
!                  in DDASPK prologue.  If IRES = -1, then IERNEW
!                  will be set to 1.
!                  If IRES < -1, then IERNEW will be set to -1.
!     IDUM      -- Dummy argument.
!     IERNEW    -- Error flag for Newton iteration.
!                   0  ==> Newton iteration converged.
!                   1  ==> recoverable error inside Newton iteration.
!                  -1  ==> unrecoverable error inside Newton iteration.
!
!     All arguments with "DUM" in their names are dummy arguments
!     which are not used in this routine.
!-----------------------------------------------------------------------
!
!***ROUTINES CALLED
!   DSLVD, DDWNRM, RES
!
!***END PROLOGUE  DNSD
!
!
IMPLICIT DOUBLE PRECISION(A-H,O-Z)
DIMENSION Y(*),YPRIME(*),WT(*),DELTA(*),E(*)
DIMENSION WM(*),IWM(*), RPAR(*),IPAR(*)
EXTERNAL  RES
!
INTEGER, PARAMETER :: LNRE=12, LNNI=19
!
!     Initialize Newton counter M and accumulation vector E. 
!
M = 0
DO 100 I=1,NEQ
  E(I)=0.0D0
100 CONTINUE
!
!     Corrector loop.
!
300 CONTINUE
IWM(LNNI) = IWM(LNNI) + 1
!
!     If necessary, multiply residual by convergence factor.
!
IF (MULDEL .EQ. 1) THEN
   DO 320 I = 1,NEQ
     DELTA(I) = DELTA(I) * CONFAC
320 CONTINUE
  ENDIF
!
!     Compute a new iterate (back-substitution).
!     Store the correction in DELTA.
!
CALL DSLVD(NEQ,DELTA,WM,IWM)
!
!     Update Y, E, and YPRIME.
!
DO 340 I=1,NEQ
   Y(I)=Y(I)-DELTA(I)
   E(I)=E(I)-DELTA(I)
   YPRIME(I)=YPRIME(I)-CJ*DELTA(I)
340 CONTINUE
!
!     Test for convergence of the iteration.
!
DELNRM=DDWNRM(NEQ,DELTA,WT,RPAR,IPAR)
IF (DELNRM .LE. TOLNEW) GO TO 370
IF (M .EQ. 0) THEN
  OLDNRM = DELNRM
ELSE
  RATE = (DELNRM/OLDNRM)**(1.0D0/M)
  IF (RATE .GT. 0.9D0) GO TO 380
  S = RATE/(1.0D0 - RATE)
ENDIF
IF (S*DELNRM .LE. EPCON) GO TO 370
!
!     The corrector has not yet converged.
!     Update M and test whether the
!     maximum number of iterations have
!     been tried.
!
M=M+1
IF(M.GE.MAXIT) GO TO 380
!
!     Evaluate the residual,
!     and go back to do another iteration.
!
IWM(LNRE)=IWM(LNRE)+1
CALL RES(X,Y,YPRIME,CJ,DELTA,IRES,RPAR,IPAR)
IF (IRES .LT. 0) GO TO 380
GO TO 300
!
!     The iteration has converged.
!
370 RETURN
!
!     The iteration has not converged.  Set IERNEW appropriately.
!
380 CONTINUE
IF (IRES .LE. -2 ) THEN
   IERNEW = -1
ELSE
   IERNEW = 1
ENDIF
RETURN
!
!
!------END OF SUBROUTINE DNSD-------------------------------------------
END
! Work performed under the auspices of the U.S. Department of Energy
! by Lawrence Livermore National Laboratory under contract number 
! W-7405-Eng-48.
!
SUBROUTINE DMATD(NEQ,X,Y,YPRIME,DELTA,CJ,H,IER,EWT,E, &
& WM,IWM,RES,IRES,UROUND,JACD,RPAR,IPAR)
!
!***BEGIN PROLOGUE  DMATD
!***REFER TO  DDASPK
!***DATE WRITTEN   890101   (YYMMDD)
!***REVISION DATE  900926   (YYMMDD)
!***REVISION DATE  940701   (YYMMDD) (new LIPVT)
!
!-----------------------------------------------------------------------
!***DESCRIPTION
!
!     This routine computes the iteration matrix
!     J = dG/dY+CJ*dG/dYPRIME (where G(X,Y,YPRIME)=0).
!     Here J is computed by:
!       the user-supplied routine JACD if IWM(MTYPE) is 1 or 4, or
!       by numerical difference quotients if IWM(MTYPE) is 2 or 5.
!
!     The parameters have the following meanings.
!     X        = Independent variable.
!     Y        = Array containing predicted values.
!     YPRIME   = Array containing predicted derivatives.
!     DELTA    = Residual evaluated at (X,Y,YPRIME).
!                (Used only if IWM(MTYPE)=2 or 5).
!     CJ       = Scalar parameter defining iteration matrix.
!     H        = Current stepsize in integration.
!     IER      = Variable which is .NE. 0 if iteration matrix
!                is singular, and 0 otherwise.
!     EWT      = Vector of error weights for computing norms.
!     E        = Work space (temporary) of length NEQ.
!     WM       = Real work space for matrices.  On output
!                it contains the LU decomposition
!                of the iteration matrix.
!     IWM      = Integer work space containing
!                matrix information.
!     RES      = External user-supplied subroutine
!                to evaluate the residual.  See RES description
!                in DDASPK prologue.
!     IRES     = Flag which is equal to zero if no illegal values
!                in RES, and less than zero otherwise.  (If IRES
!                is less than zero, the matrix was not completed).
!                In this case (if IRES .LT. 0), then IER = 0.
!     UROUND   = The unit roundoff error of the machine being used.
!     JACD     = Name of the external user-supplied routine
!                to evaluate the iteration matrix.  (This routine
!                is only used if IWM(MTYPE) is 1 or 4)
!                See JAC description for the case INFO(12) = 0
!                in DDASPK prologue.
!     RPAR,IPAR= Real and integer parameter arrays that
!                are used for communication between the
!                calling program and external user routines.
!                They are not altered by DMATD.
!-----------------------------------------------------------------------
!***ROUTINES CALLED
!   JACD, RES, DGEFA, DGBFA
!
!***END PROLOGUE  DMATD
!
!
IMPLICIT DOUBLE PRECISION(A-H,O-Z)
DIMENSION Y(*),YPRIME(*),DELTA(*),EWT(*),E(*)
DIMENSION WM(*),IWM(*), RPAR(*),IPAR(*)
EXTERNAL  RES, JACD
!
INTEGER, PARAMETER :: LML=1, LMU=2, LMTYPE=4, LNRE=12, LNPD=22, LLCIWP=30
!
LIPVT = IWM(LLCIWP)
IER = 0
MTYPE=IWM(LMTYPE)
IF (MTYPE .EQ. 1) THEN
  GOTO 100
ELSE IF (MTYPE .EQ. 2) THEN
  GOTO 200
ELSE IF (MTYPE .EQ. 3) THEN
  GOTO 300
ELSE IF (MTYPE .EQ. 4) THEN
  GOTO 400
ELSE IF (MTYPE .EQ. 5) THEN
  GOTO 500
ENDIF
!      GO TO (100,200,300,400,500),MTYPE
!
!
!     Dense user-supplied matrix.
!
100 LENPD=IWM(LNPD)
DO 110 I=1,LENPD
   WM(I)=0.0D0
110 CONTINUE
CALL JACD(X,Y,YPRIME,WM,CJ,RPAR,IPAR)
GO TO 230
!
!
!     Dense finite-difference-generated matrix.
!
200 IRES=0
NROW=0
SQUR = SQRT(UROUND)
DO 210 I=1,NEQ
   DEL=SQUR*MAX(ABS(Y(I)),ABS(H*YPRIME(I)), &
& ABS(1.D0/EWT(I)))
   DEL=SIGN(DEL,H*YPRIME(I))
   DEL=(Y(I)+DEL)-Y(I)
   YSAVE=Y(I)
   YPSAVE=YPRIME(I)
   Y(I)=Y(I)+DEL
   YPRIME(I)=YPRIME(I)+CJ*DEL
   IWM(LNRE)=IWM(LNRE)+1
   CALL RES(X,Y,YPRIME,CJ,E,IRES,RPAR,IPAR)
   IF (IRES .LT. 0) RETURN
   DELINV=1.0D0/DEL
   DO 220 L=1,NEQ
     WM(NROW+L)=(E(L)-DELTA(L))*DELINV
220 CONTINUE
NROW=NROW+NEQ
Y(I)=YSAVE
YPRIME(I)=YPSAVE
210 CONTINUE
!
!
!     Do dense-matrix LU decomposition on J.
!
230 CALL DGEFA(WM,NEQ,NEQ,IWM(LIPVT),IER)
RETURN
!
!
!     Dummy section for IWM(MTYPE)=3.
!
300 RETURN
!
!
!     Banded user-supplied matrix.
!
400 LENPD=IWM(LNPD)
DO 410 I=1,LENPD
   WM(I)=0.0D0
410 CONTINUE
CALL JACD(X,Y,YPRIME,WM,CJ,RPAR,IPAR)
MEBAND=2*IWM(LML)+IWM(LMU)+1
GO TO 550
!
!
!     Banded finite-difference-generated matrix.
!
500 MBAND=IWM(LML)+IWM(LMU)+1
MBA=MIN0(MBAND,NEQ)
MEBAND=MBAND+IWM(LML)
MEB1=MEBAND-1
MSAVE=(NEQ/MBAND)+1
ISAVE=IWM(LNPD)
IPSAVE=ISAVE+MSAVE
IRES=0
SQUR=SQRT(UROUND)
DO 540 J=1,MBA
  DO 510 N=J,NEQ,MBAND
    K= (N-J)/MBAND + 1
    WM(ISAVE+K)=Y(N)
    WM(IPSAVE+K)=YPRIME(N)
    DEL=SQUR*MAX(ABS(Y(N)),ABS(H*YPRIME(N)), &
& ABS(1.D0/EWT(N)))
    DEL=SIGN(DEL,H*YPRIME(N))
    DEL=(Y(N)+DEL)-Y(N)
    Y(N)=Y(N)+DEL
    YPRIME(N)=YPRIME(N)+CJ*DEL
510 CONTINUE
  IWM(LNRE)=IWM(LNRE)+1
  CALL RES(X,Y,YPRIME,CJ,E,IRES,RPAR,IPAR)
  IF (IRES .LT. 0) RETURN
  DO 530 N=J,NEQ,MBAND
    K= (N-J)/MBAND + 1
    Y(N)=WM(ISAVE+K)
    YPRIME(N)=WM(IPSAVE+K)
    DEL=SQUR*MAX(ABS(Y(N)),ABS(H*YPRIME(N)), &
& ABS(1.D0/EWT(N)))
    DEL=SIGN(DEL,H*YPRIME(N))
    DEL=(Y(N)+DEL)-Y(N)
    DELINV=1.0D0/DEL
    I1=MAX0(1,(N-IWM(LMU)))
    I2=MIN0(NEQ,(N+IWM(LML)))
    II=N*MEB1-IWM(LML)
    DO 520 I=I1,I2
      WM(II+I)=(E(I)-DELTA(I))*DELINV
520 CONTINUE
530 CONTINUE
540 CONTINUE
!
!
!     Do LU decomposition of banded J.
!
550 CALL DGBFA (WM,MEBAND,NEQ,IWM(LML),IWM(LMU),IWM(LIPVT),IER)
RETURN
!
!------END OF SUBROUTINE DMATD------------------------------------------
END
! Work performed under the auspices of the U.S. Department of Energy
! by Lawrence Livermore National Laboratory under contract number 
! W-7405-Eng-48.
!
SUBROUTINE DSLVD(NEQ,DELTA,WM,IWM)
!
!***BEGIN PROLOGUE  DSLVD
!***REFER TO  DDASPK
!***DATE WRITTEN   890101   (YYMMDD)
!***REVISION DATE  900926   (YYMMDD)
!***REVISION DATE  940701   (YYMMDD) (new LIPVT)
!
!-----------------------------------------------------------------------
!***DESCRIPTION
!
!     This routine manages the solution of the linear
!     system arising in the Newton iteration.
!     Real matrix information and real temporary storage
!     is stored in the array WM.
!     Integer matrix information is stored in the array IWM.
!     For a dense matrix, the LINPACK routine DGESL is called.
!     For a banded matrix, the LINPACK routine DGBSL is called.
!-----------------------------------------------------------------------
!***ROUTINES CALLED
!   DGESL, DGBSL
!
!***END PROLOGUE  DSLVD
!
!
IMPLICIT DOUBLE PRECISION(A-H,O-Z)
DIMENSION DELTA(*),WM(*),IWM(*)
!
INTEGER, PARAMETER :: LML=1, LMU=2, LMTYPE=4, LLCIWP=30
!
LIPVT = IWM(LLCIWP)
MTYPE=IWM(LMTYPE)
IF (MTYPE .EQ. 1 .OR. MTYPE .EQ. 2) THEN
  GOTO 100
ELSE IF (MTYPE .EQ. 3) THEN
  GOTO 300
ELSE IF (MTYPE .EQ. 4 .OR. MTYPE .EQ. 5) THEN
  GOTO 400
ENDIF


!      GO TO(100,100,300,400,400),MTYPE
!
!     Dense matrix.
!
100 CALL DGESL(WM,NEQ,NEQ,IWM(LIPVT),DELTA,0)
RETURN
!
!     Dummy section for MTYPE=3.
!
300 CONTINUE
RETURN
!
!     Banded matrix.
!
400 MEBAND=2*IWM(LML)+IWM(LMU)+1
CALL DGBSL(WM,MEBAND,NEQ,IWM(LML), &
& IWM(LMU),IWM(LIPVT),DELTA,0)
RETURN
!
!------END OF SUBROUTINE DSLVD------------------------------------------
END
! Work perfored under the auspices of the U.S. Department of Energy
! by Lawrence Livermore National Laboratory under contract number 
! W-7405-Eng-48.
!
SUBROUTINE DDASIK(X,Y,YPRIME,NEQ,ICOPT,ID,RES,JACK,PSOL,H,WT, &
& JSKIP,RPAR,IPAR,SAVR,DELTA,R,YIC,YPIC,PWK,WM,IWM,CJ,UROUND, &
& EPLI,SQRTN,RSQRTN,EPCON,RATEMX,STPTOL,JFLG, &
& ICNFLG,ICNSTR,IERNLS)
!
!***BEGIN PROLOGUE  DDASIK
!***REFER TO  DDASPK
!***DATE WRITTEN   941026   (YYMMDD)
!***REVISION DATE  950808   (YYMMDD)
!***REVISION DATE  951110   Removed unreachable block 390.
!
!
!-----------------------------------------------------------------------
!***DESCRIPTION
!
!
!     DDASIK solves a nonlinear system of algebraic equations of the
!     form G(X,Y,YPRIME) = 0 for the unknown parts of Y and YPRIME in
!     the initial conditions.
!
!     An initial value for Y and initial guess for YPRIME are input.
!
!     The method used is a Newton scheme with Krylov iteration and a
!     linesearch algorithm.
!
!     The parameters represent
!
!     X         -- Independent variable.
!     Y         -- Solution vector at x.
!     YPRIME    -- Derivative of solution vector.
!     NEQ       -- Number of equations to be integrated.
!     ICOPT     -- Initial condition option chosen (1 or 2).
!     ID        -- Array of dimension NEQ, which must be initialized
!                  if ICOPT = 1.  See DDASIC.
!     RES       -- External user-supplied subroutine
!                  to evaluate the residual.  See RES description
!                  in DDASPK prologue.
!     JACK     --  External user-supplied routine to update
!                  the preconditioner.  (This is optional).
!                  See JAC description for the case
!                  INFO(12) = 1 in the DDASPK prologue.
!     PSOL      -- External user-supplied routine to solve
!                  a linear system using preconditioning.
!                  (This is optional).  See explanation inside DDASPK.
!     H         -- Scaling factor for this initial condition calc.
!     WT        -- Vector of weights for error criterion.
!     JSKIP     -- input flag to signal if initial JAC call is to be
!                  skipped.  1 => skip the call, 0 => do not skip call.
!     RPAR,IPAR -- Real and integer arrays used for communication
!                  between the calling program and external user
!                  routines.  They are not altered within DASPK.
!     SAVR      -- Work vector for DDASIK of length NEQ.
!     DELTA     -- Work vector for DDASIK of length NEQ.
!     R         -- Work vector for DDASIK of length NEQ.
!     YIC,YPIC  -- Work vectors for DDASIK, each of length NEQ.
!     PWK       -- Work vector for DDASIK of length NEQ.
!     WM,IWM    -- Real and integer arrays storing
!                  matrix information for linear system
!                  solvers, and various other information.
!     CJ        -- Matrix parameter = 1/H (ICOPT = 1) or 0 (ICOPT = 2).
!     UROUND    -- Unit roundoff.
!     EPLI      -- convergence test constant.
!                  See DDASPK prologue for more details.
!     SQRTN     -- Square root of NEQ.
!     RSQRTN    -- reciprical of square root of NEQ.
!     EPCON     -- Tolerance to test for convergence of the Newton
!                  iteration.
!     RATEMX    -- Maximum convergence rate for which Newton iteration
!                  is considered converging.
!     JFLG      -- Flag showing whether a Jacobian routine is supplied.
!     ICNFLG    -- Integer scalar.  If nonzero, then constraint
!                  violations in the proposed new approximate solution
!                  will be checked for, and the maximum step length 
!                  will be adjusted accordingly.
!     ICNSTR    -- Integer array of length NEQ containing flags for
!                  checking constraints.
!     IERNLS    -- Error flag for nonlinear solver.
!                   0   ==> nonlinear solver converged.
!                   1,2 ==> recoverable error inside nonlinear solver.
!                           1 => retry with current Y, YPRIME
!                           2 => retry with original Y, YPRIME
!                  -1   ==> unrecoverable error in nonlinear solver.
!
!-----------------------------------------------------------------------
!
!***ROUTINES CALLED
!   RES, JACK, DNSIK, DCOPY
!
!***END PROLOGUE  DDASIK
!
!
IMPLICIT DOUBLE PRECISION(A-H,O-Z)
DIMENSION Y(*),YPRIME(*),ID(*),WT(*),ICNSTR(*)
DIMENSION SAVR(*),DELTA(*),R(*),YIC(*),YPIC(*),PWK(*)
DIMENSION WM(*),IWM(*), RPAR(*),IPAR(*)
EXTERNAL RES, JACK, PSOL
!
INTEGER, PARAMETER :: LNRE=12, LNJE=13, LLOCWP=29, LLCIWP=30
INTEGER, PARAMETER :: LMXNIT=32, LMXNJ=33
!
!
!     Perform initializations.
!
LWP = IWM(LLOCWP)
LIWP = IWM(LLCIWP)
MXNIT = IWM(LMXNIT)
MXNJ = IWM(LMXNJ)
IERNLS = 0
NJ = 0
EPLIN = EPLI*EPCON
!
!     Call RES to initialize DELTA.
!
IRES = 0
IWM(LNRE) = IWM(LNRE) + 1
CALL RES(X,Y,YPRIME,CJ,DELTA,IRES,RPAR,IPAR)
IF (IRES .LT. 0) GO TO 370
!
!     Looping point for updating the preconditioner.
!
300 CONTINUE
!
!     Initialize all error flags to zero.
!
IERPJ = 0
IRES = 0
IERNEW = 0
!
!     If a Jacobian routine was supplied, call it.
!
IF (JFLG .EQ. 1 .AND. JSKIP .EQ. 0) THEN
  NJ = NJ + 1
  IWM(LNJE)=IWM(LNJE)+1
  CALL JACK (RES, IRES, NEQ, X, Y, YPRIME, WT, DELTA, R, H, CJ, &
& WM(LWP), IWM(LIWP), IERPJ, RPAR, IPAR)
  IF (IRES .LT. 0 .OR. IERPJ .NE. 0) GO TO 370
  ENDIF
JSKIP = 0
!
!     Call the nonlinear Newton solver for up to MXNIT iterations.
!
CALL DNSIK(X,Y,YPRIME,NEQ,ICOPT,ID,RES,PSOL,WT,RPAR,IPAR, &
& SAVR,DELTA,R,YIC,YPIC,PWK,WM,IWM,CJ,SQRTN,RSQRTN, &
& EPLIN,EPCON,RATEMX,MXNIT,STPTOL,ICNFLG,ICNSTR,IERNEW)
!
IF (IERNEW .EQ. 1 .AND. NJ .LT. MXNJ .AND. JFLG .EQ. 1) THEN
!
!       Up to MXNIT iterations were done, the convergence rate is < 1,
!       a Jacobian routine is supplied, and the number of JACK calls
!       is less than MXNJ.  
!       Copy the residual SAVR to DELTA, call JACK, and try again.
!
  CALL DCOPY (NEQ,  SAVR, 1, DELTA, 1)
  GO TO 300
  ENDIF
!
IF (IERNEW .NE. 0) GO TO 380
RETURN
!
!
!     Unsuccessful exits from nonlinear solver.
!     Set IERNLS accordingly.
!
370 IERNLS = 2
IF (IRES .LE. -2) IERNLS = -1
RETURN
!
380 IERNLS = MIN(IERNEW,2)
RETURN
!
!----------------------- END OF SUBROUTINE DDASIK-----------------------
END
! Work performed under the auspices of the U.S. Department of Energy
! by Lawrence Livermore National Laboratory under contract number 
! W-7405-Eng-48.
!
SUBROUTINE DNSIK(X,Y,YPRIME,NEQ,ICOPT,ID,RES,PSOL,WT,RPAR,IPAR, &
& SAVR,DELTA,R,YIC,YPIC,PWK,WM,IWM,CJ,SQRTN,RSQRTN,EPLIN,EPCON, &
& RATEMX,MAXIT,STPTOL,ICNFLG,ICNSTR,IERNEW)
!
!***BEGIN PROLOGUE  DNSIK
!***REFER TO  DDASPK
!***DATE WRITTEN   940701   (YYMMDD)
!***REVISION DATE  950714   (YYMMDD)
!
!
!-----------------------------------------------------------------------
!***DESCRIPTION
!
!     DNSIK solves a nonlinear system of algebraic equations of the
!     form G(X,Y,YPRIME) = 0 for the unknown parts of Y and YPRIME in
!     the initial conditions.
!
!     The method used is a Newton scheme combined with a linesearch
!     algorithm, using Krylov iterative linear system methods.
!
!     The parameters represent
!
!     X         -- Independent variable.
!     Y         -- Solution vector.
!     YPRIME    -- Derivative of solution vector.
!     NEQ       -- Number of unknowns.
!     ICOPT     -- Initial condition option chosen (1 or 2).
!     ID        -- Array of dimension NEQ, which must be initialized
!                  if ICOPT = 1.  See DDASIC.
!     RES       -- External user-supplied subroutine
!                  to evaluate the residual.  See RES description
!                  in DDASPK prologue.
!     PSOL      -- External user-supplied routine to solve
!                  a linear system using preconditioning. 
!                  See explanation inside DDASPK.
!     WT        -- Vector of weights for error criterion.
!     RPAR,IPAR -- Real and integer arrays used for communication
!                  between the calling program and external user
!                  routines.  They are not altered within DASPK.
!     SAVR      -- Work vector for DNSIK of length NEQ.
!     DELTA     -- Residual vector on entry, and work vector of
!                  length NEQ for DNSIK.
!     R         -- Work vector for DNSIK of length NEQ.
!     YIC,YPIC  -- Work vectors for DNSIK, each of length NEQ.
!     PWK       -- Work vector for DNSIK of length NEQ.
!     WM,IWM    -- Real and integer arrays storing
!                  matrix information such as the matrix
!                  of partial derivatives, permutation
!                  vector, and various other information.
!     CJ        -- Matrix parameter = 1/H (ICOPT = 1) or 0 (ICOPT = 2).
!     SQRTN     -- Square root of NEQ.
!     RSQRTN    -- reciprical of square root of NEQ.
!     EPLIN     -- Tolerance for linear system solver.
!     EPCON     -- Tolerance to test for convergence of the Newton
!                  iteration.
!     RATEMX    -- Maximum convergence rate for which Newton iteration
!                  is considered converging.
!     MAXIT     -- Maximum allowed number of Newton iterations.
!     STPTOL    -- Tolerance used in calculating the minimum lambda
!                  value allowed.
!     ICNFLG    -- Integer scalar.  If nonzero, then constraint
!                  violations in the proposed new approximate solution
!                  will be checked for, and the maximum step length
!                  will be adjusted accordingly.
!     ICNSTR    -- Integer array of length NEQ containing flags for
!                  checking constraints.
!     IERNEW    -- Error flag for Newton iteration.
!                   0  ==> Newton iteration converged.
!                   1  ==> failed to converge, but RATE .lt. 1.
!                   2  ==> failed to converge, RATE .gt. RATEMX.
!                   3  ==> other recoverable error.
!                  -1  ==> unrecoverable error inside Newton iteration.
!-----------------------------------------------------------------------
!
!***ROUTINES CALLED
!   DFNRMK, DSLVK, DDWNRM, DLINSK, DCOPY
!
!***END PROLOGUE  DNSIK
!
!
IMPLICIT DOUBLE PRECISION(A-H,O-Z)
DIMENSION Y(*),YPRIME(*),WT(*),ID(*),DELTA(*),R(*),SAVR(*)
DIMENSION YIC(*),YPIC(*),PWK(*),WM(*),IWM(*), RPAR(*),IPAR(*)
DIMENSION ICNSTR(*)
EXTERNAL RES, PSOL
!
INTEGER, PARAMETER :: LNNI=19, LNPS=21, LLOCWP=29, LLCIWP=30
INTEGER, PARAMETER :: LLSOFF=35, LSTOL=14
!
!
!     Initializations.  M is the Newton iteration counter.
!
LSOFF = IWM(LLSOFF)
M = 0
RATE = 1.0D0
LWP = IWM(LLOCWP)
LIWP = IWM(LLCIWP)
RLX = 0.4D0
!
!     Save residual in SAVR.
!
CALL DCOPY (NEQ, DELTA, 1, SAVR, 1)
!
!     Compute norm of (P-inverse)*(residual).
!
CALL DFNRMK (NEQ, Y, X, YPRIME, SAVR, R, CJ, WT, SQRTN, RSQRTN, &
& RES, IRES, PSOL, 1, IER, FNRM, EPLIN, WM(LWP), IWM(LIWP), &
& PWK, RPAR, IPAR)
IWM(LNPS) = IWM(LNPS) + 1
IF (IER .NE. 0) THEN
  IERNEW = 3
  RETURN
ENDIF
!
!     Return now if residual norm is .le. EPCON.
!
IF (FNRM .LE. EPCON) RETURN
!
!     Newton iteration loop.
!
300 CONTINUE
IWM(LNNI) = IWM(LNNI) + 1
!
!     Compute a new step vector DELTA.
!
CALL DSLVK (NEQ, Y, X, YPRIME, SAVR, DELTA, WT, WM, IWM, &
& RES, IRES, PSOL, IERSL, CJ, EPLIN, SQRTN, RSQRTN, RHOK, &
& RPAR, IPAR)
IF (IRES .NE. 0 .OR. IERSL .NE. 0) GO TO 390
!
!     Get norm of DELTA.  Return now if DELTA is zero.
!
DELNRM = DDWNRM(NEQ,DELTA,WT,RPAR,IPAR)
IF (DELNRM .EQ. 0.0D0) RETURN
!
!     Call linesearch routine for global strategy and set RATE.
!
OLDFNM = FNRM
!
CALL DLINSK (NEQ, Y, X, YPRIME, SAVR, CJ, DELTA, DELNRM, WT, &
& SQRTN, RSQRTN, LSOFF, STPTOL, IRET, RES, IRES, PSOL, WM, IWM, &
& RHOK, FNRM, ICOPT, ID, WM(LWP), IWM(LIWP), R, EPLIN, YIC, YPIC, &
& PWK, ICNFLG, ICNSTR, RLX, RPAR, IPAR)
!
RATE = FNRM/OLDFNM
!
!     Check for error condition from linesearch.
IF (IRET .NE. 0) GO TO 390
!
!     Test for convergence of the iteration, and return or loop.
!
IF (FNRM .LE. EPCON) RETURN
!
!     The iteration has not yet converged.  Update M.
!     Test whether the maximum number of iterations have been tried.
!
M=M+1
IF(M .GE. MAXIT) GO TO 380
!
!     Copy the residual SAVR to DELTA and loop for another iteration.
!
CALL DCOPY (NEQ,  SAVR, 1, DELTA, 1)
GO TO 300
!
!     The maximum number of iterations was done.  Set IERNEW and return.
!
380 IF (RATE .LE. RATEMX) THEN
   IERNEW = 1
ELSE
   IERNEW = 2
ENDIF
RETURN
!
390 IF (IRES .LE. -2 .OR. IERSL .LT. 0) THEN
   IERNEW = -1
ELSE
   IERNEW = 3
   IF (IRES .EQ. 0 .AND. IERSL .EQ. 1 .AND. M .GE. 2 &
& .AND. RATE .LT. 1.0D0) IERNEW = 1
ENDIF
RETURN
!
!
!----------------------- END OF SUBROUTINE DNSIK------------------------
END
! Work performed under the auspices of the U.S. Department of Energy
! by Lawrence Livermore National Laboratory under contract number 
! W-7405-Eng-48.
!
SUBROUTINE DLINSK (NEQ, Y, T, YPRIME, SAVR, CJ, P, PNRM, WT, &
& SQRTN, RSQRTN, LSOFF, STPTOL, IRET, RES, IRES, PSOL, WM, IWM, &
& RHOK, FNRM, ICOPT, ID, WP, IWP, R, EPLIN, YNEW, YPNEW, PWK, &
& ICNFLG, ICNSTR, RLX, RPAR, IPAR)
!
!***BEGIN PROLOGUE  DLINSK
!***REFER TO  DNSIK
!***DATE WRITTEN   940830   (YYMMDD)
!***REVISION DATE  951006   (Arguments SQRTN, RSQRTN added.)
!***REVISION DATE  960129   Moved line RL = ONE to top block.
!
!
!-----------------------------------------------------------------------
!***DESCRIPTION
!
!     DLINSK uses a linesearch algorithm to calculate a new (Y,YPRIME)
!     pair (YNEW,YPNEW) such that 
!
!     f(YNEW,YPNEW) .le. (1 - 2*ALPHA*RL)*f(Y,YPRIME) + 
!                          ALPHA*RL*RHOK*RHOK ,
!
!     where 0 < RL <= 1, and RHOK is the scaled preconditioned norm of
!     the final residual vector in the Krylov iteration.  
!     Here, f(y,y') is defined as
!
!      f(y,y') = (1/2)*norm( (P-inverse)*G(t,y,y') )**2 ,
!
!     where norm() is the weighted RMS vector norm, G is the DAE
!     system residual function, and P is the preconditioner used
!     in the Krylov iteration.
!
!     In addition to the parameters defined elsewhere, we have
!
!     SAVR    -- Work array of length NEQ, containing the residual
!                vector G(t,y,y') on return.
!     P       -- Approximate Newton step used in backtracking.
!     PNRM    -- Weighted RMS norm of P.
!     LSOFF   -- Flag showing whether the linesearch algorithm is
!                to be invoked.  0 means do the linesearch, 
!                1 means turn off linesearch.
!     STPTOL  -- Tolerance used in calculating the minimum lambda
!                value allowed.
!     ICNFLG  -- Integer scalar.  If nonzero, then constraint violations
!                in the proposed new approximate solution will be
!                checked for, and the maximum step length will be
!                adjusted accordingly.
!     ICNSTR  -- Integer array of length NEQ containing flags for
!                checking constraints.
!     RHOK    -- Weighted norm of preconditioned Krylov residual.
!     RLX     -- Real scalar restricting update size in DCNSTR.
!     YNEW    -- Array of length NEQ used to hold the new Y in
!                performing the linesearch.
!     YPNEW   -- Array of length NEQ used to hold the new YPRIME in
!                performing the linesearch.
!     PWK     -- Work vector of length NEQ for use in PSOL.
!     Y       -- Array of length NEQ containing the new Y (i.e.,=YNEW).
!     YPRIME  -- Array of length NEQ containing the new YPRIME 
!                (i.e.,=YPNEW).
!     FNRM    -- Real scalar containing SQRT(2*f(Y,YPRIME)) for the
!                current (Y,YPRIME) on input and output.
!     R       -- Work space length NEQ for residual vector.
!     IRET    -- Return flag.
!                IRET=0 means that a satisfactory (Y,YPRIME) was found.
!                IRET=1 means that the routine failed to find a new
!                       (Y,YPRIME) that was sufficiently distinct from
!                       the current (Y,YPRIME) pair.
!                IRET=2 means a failure in RES or PSOL.
!-----------------------------------------------------------------------
!
!***ROUTINES CALLED
!   DFNRMK, DYYPNW, DCOPY
!
!***END PROLOGUE  DLINSK
!
IMPLICIT DOUBLE PRECISION(A-H,O-Z)
EXTERNAL  RES, PSOL
DIMENSION Y(*), YPRIME(*), P(*), WT(*), SAVR(*), R(*), ID(*)
DIMENSION WM(*), IWM(*), YNEW(*), YPNEW(*), PWK(*), ICNSTR(*)
DIMENSION WP(*), IWP(*), RPAR(*), IPAR(*)
!
INTEGER, PARAMETER :: LNRE=12, LNPS=21, LKPRIN=31
!
SAVE ALPHA, ONE, TWO
DATA ALPHA/1.0D-4/, ONE/1.0D0/, TWO/2.0D0/
!
KPRIN=IWM(LKPRIN)
F1NRM = (FNRM*FNRM)/TWO
RATIO = ONE
!
IF (KPRIN .GE. 2) THEN
call rprintfd1( &
& 'daspk--  in routine dlinsd--PNRM (= %g)' // char(0),PNRM)

  ENDIF
TAU = PNRM
IVIO = 0
RL = ONE
!-----------------------------------------------------------------------
! Check for violations of the constraints, if any are imposed.
! If any violations are found, the step vector P is rescaled, and the 
! constraint check is repeated, until no violations are found.
!-----------------------------------------------------------------------
IF (ICNFLG .NE. 0) THEN
10 CONTINUE
   CALL DYYPNW (NEQ,Y,YPRIME,CJ,RL,P,ICOPT,ID,YNEW,YPNEW)
   CALL DCNSTR (NEQ, Y, YNEW, ICNSTR, TAU, RLX, IRET, IVAR)
   IF (IRET .EQ. 1) THEN
      IVIO = 1
      RATIO1 = TAU/PNRM
      RATIO = RATIO*RATIO1
      DO 20 I = 1,NEQ
        P(I) = P(I)*RATIO1
20 CONTINUE
      PNRM = TAU
      IF (KPRIN .GE. 2) THEN
call rprintfid( &
& 'daspk--  constraint violation, PNRM(%g), INDEX(%i)' &
& // char(0),IVAR,PNRM)

        ENDIF
      IF (PNRM .LE. STPTOL) THEN
        IRET = 1
        RETURN
        ENDIF
      GO TO 10
      ENDIF
   ENDIF
!
SLPI = (-TWO*F1NRM + RHOK*RHOK)*RATIO
RLMIN = STPTOL/PNRM
IF (LSOFF .EQ. 0 .AND. KPRIN .GE. 2) THEN
call rprintfd1( &
& 'daspk--  Min. LAMBDA &g' // char(0), RLMIN)
  ENDIF
!-----------------------------------------------------------------------
! Begin iteration to find RL value satisfying alpha-condition.
! Update YNEW and YPNEW, then compute norm of new scaled residual and
! perform alpha condition test.
!-----------------------------------------------------------------------
100 CONTINUE
CALL DYYPNW (NEQ,Y,YPRIME,CJ,RL,P,ICOPT,ID,YNEW,YPNEW)
CALL DFNRMK (NEQ, YNEW, T, YPNEW, SAVR, R, CJ, WT, SQRTN, RSQRTN, &
& RES, IRES, PSOL, 0, IER, FNRMP, EPLIN, WP, IWP, PWK, RPAR, IPAR)
IWM(LNRE) = IWM(LNRE) + 1
IF (IRES .GE. 0) IWM(LNPS) = IWM(LNPS) + 1
IF (IRES .NE. 0 .OR. IER .NE. 0) THEN
  IRET = 2
  RETURN
  ENDIF
IF (LSOFF .EQ. 1) GO TO 150
!
F1NRMP = FNRMP*FNRMP/TWO
IF (KPRIN .GE. 2) THEN
call rprintfd1( &
& 'daspk--  LAMBDA (= %g)' // char(0), RL)
call rprintfd2( &
& '     --  NORM(F1) (= %g), NORM(F1NEW) (= %g)' &
& // char(0), F1NRM, F1NRMP)

  ENDIF
IF (F1NRMP .GT. F1NRM + ALPHA*SLPI*RL) GO TO 200
!-----------------------------------------------------------------------
! Alpha-condition is satisfied, or linesearch is turned off.
! Copy YNEW,YPNEW to Y,YPRIME and return.
!-----------------------------------------------------------------------
150 IRET = 0
CALL DCOPY(NEQ, YNEW, 1, Y, 1)
CALL DCOPY(NEQ, YPNEW, 1, YPRIME, 1)
FNRM = FNRMP
IF (KPRIN .GE. 1) THEN
call rprintfd1( &
& 'daspk--  leaving routine dlinsk--FNRM %g' // char(0), FNRM)

  ENDIF
RETURN
!-----------------------------------------------------------------------
! Alpha-condition not satisfied.  Perform backtrack to compute new RL
! value.  If RL is less than RLMIN, i.e. no satisfactory YNEW,YPNEW can
! be found sufficiently distinct from Y,YPRIME, then return IRET = 1.
!-----------------------------------------------------------------------
200 CONTINUE
IF (RL .LT. RLMIN) THEN
  IRET = 1
  RETURN
  ENDIF
!
RL = RL/TWO
GO TO 100
!
!----------------------- END OF SUBROUTINE DLINSK ----------------------
END
! Work performed under the auspices of the U.S. Department of Energy
! by Lawrence Livermore National Laboratory under contract number 
! W-7405-Eng-48.
!
SUBROUTINE DFNRMK (NEQ, Y, T, YPRIME, SAVR, R, CJ, WT, &
& SQRTN, RSQRTN, RES, IRES, PSOL, IRIN, IER, &
& FNORM, EPLIN, WP, IWP, PWK, RPAR, IPAR)
!
!***BEGIN PROLOGUE  DFNRMK
!***REFER TO  DLINSK
!***DATE WRITTEN   940830   (YYMMDD)
!***REVISION DATE  951006   (SQRTN, RSQRTN, and scaling of WT added.)
!
!
!-----------------------------------------------------------------------
!***DESCRIPTION
!
!     DFNRMK calculates the scaled preconditioned norm of the nonlinear
!     function used in the nonlinear iteration for obtaining consistent
!     initial conditions.  Specifically, DFNRMK calculates the weighted
!     root-mean-square norm of the vector (P-inverse)*G(T,Y,YPRIME),
!     where P is the preconditioner matrix.
!
!     In addition to the parameters described in the calling program
!     DLINSK, the parameters represent
!
!     IRIN   -- Flag showing whether the current residual vector is
!               input in SAVR.  1 means it is, 0 means it is not.
!     R      -- Array of length NEQ that contains
!               (P-inverse)*G(T,Y,YPRIME) on return.
!     FNORM  -- Scalar containing the weighted norm of R on return.
!-----------------------------------------------------------------------
!
!***ROUTINES CALLED
!   RES, DCOPY, DSCAL, PSOL, DDWNRM
!
!***END PROLOGUE  DFNRMK
!
!
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
EXTERNAL RES, PSOL
DIMENSION Y(*), YPRIME(*), WT(*), SAVR(*), R(*), PWK(*)
DIMENSION WP(*), IWP(*), RPAR(*), IPAR(*)
!-----------------------------------------------------------------------
!     Call RES routine if IRIN = 0.
!-----------------------------------------------------------------------
IF (IRIN .EQ. 0) THEN
  IRES = 0
  CALL RES (T, Y, YPRIME, CJ, SAVR, IRES, RPAR, IPAR)
  IF (IRES .LT. 0) RETURN
  ENDIF
!-----------------------------------------------------------------------
!     Apply inverse of left preconditioner to vector R.
!     First scale WT array by 1/sqrt(N), and undo scaling afterward.
!-----------------------------------------------------------------------
CALL DCOPY(NEQ, SAVR, 1, R, 1)
CALL DSCAL (NEQ, RSQRTN, WT, 1)
IER = 0
CALL PSOL (NEQ, T, Y, YPRIME, SAVR, PWK, CJ, WT, WP, IWP, &
& R, EPLIN, IER, RPAR, IPAR)
CALL DSCAL (NEQ, SQRTN, WT, 1)
IF (IER .NE. 0) RETURN
!-----------------------------------------------------------------------
!     Calculate norm of R.
!-----------------------------------------------------------------------
FNORM = DDWNRM (NEQ, R, WT, RPAR, IPAR)
!
RETURN
!----------------------- END OF SUBROUTINE DFNRMK ----------------------
END
! Work performed under the auspices of the U.S. Department of Energy
! by Lawrence Livermore National Laboratory under contract number 
! W-7405-Eng-48.
!
SUBROUTINE DNEDK(X,Y,YPRIME,NEQ,RES,JACK,PSOL, &
& H,WT,JSTART,IDID,RPAR,IPAR,PHI,GAMMA,SAVR,DELTA,E, &
& WM,IWM,CJ,CJOLD,CJLAST,S,UROUND,EPLI,SQRTN,RSQRTN, &
& EPCON,JCALC,JFLG,KP1,NONNEG,NTYPE,IERNLS)
!
!***BEGIN PROLOGUE  DNEDK
!***REFER TO  DDASPK
!***DATE WRITTEN   891219   (YYMMDD)
!***REVISION DATE  900926   (YYMMDD)
!***REVISION DATE  940701   (YYMMDD)
!
!
!-----------------------------------------------------------------------
!***DESCRIPTION
!
!     DNEDK solves a nonlinear system of
!     algebraic equations of the form
!     G(X,Y,YPRIME) = 0 for the unknown Y.
!
!     The method used is a matrix-free Newton scheme.
!
!     The parameters represent
!     X         -- Independent variable.
!     Y         -- Solution vector at x.
!     YPRIME    -- Derivative of solution vector
!                  after successful step.
!     NEQ       -- Number of equations to be integrated.
!     RES       -- External user-supplied subroutine
!                  to evaluate the residual.  See RES description
!                  in DDASPK prologue.
!     JACK     --  External user-supplied routine to update
!                  the preconditioner.  (This is optional).
!                  See JAC description for the case
!                  INFO(12) = 1 in the DDASPK prologue.
!     PSOL      -- External user-supplied routine to solve
!                  a linear system using preconditioning. 
!                  (This is optional).  See explanation inside DDASPK.
!     H         -- Appropriate step size for this step.
!     WT        -- Vector of weights for error criterion.
!     JSTART    -- Indicates first call to this routine.
!                  If JSTART = 0, then this is the first call,
!                  otherwise it is not.
!     IDID      -- Completion flag, output by DNEDK.
!                  See IDID description in DDASPK prologue.
!     RPAR,IPAR -- Real and integer arrays used for communication
!                  between the calling program and external user
!                  routines.  They are not altered within DASPK.
!     PHI       -- Array of divided differences used by
!                  DNEDK.  The length is NEQ*(K+1), where
!                  K is the maximum order.
!     GAMMA     -- Array used to predict Y and YPRIME.  The length
!                  is K+1, where K is the maximum order.
!     SAVR      -- Work vector for DNEDK of length NEQ.
!     DELTA     -- Work vector for DNEDK of length NEQ.
!     E         -- Error accumulation vector for DNEDK of length NEQ.
!     WM,IWM    -- Real and integer arrays storing
!                  matrix information for linear system
!                  solvers, and various other information.
!     CJ        -- Parameter always proportional to 1/H.
!     CJOLD     -- Saves the value of CJ as of the last call to DITMD.
!                  Accounts for changes in CJ needed to
!                  decide whether to call DITMD.
!     CJLAST    -- Previous value of CJ.
!     S         -- A scalar determined by the approximate rate
!                  of convergence of the Newton iteration and used
!                  in the convergence test for the Newton iteration.
!
!                  If RATE is defined to be an estimate of the
!                  rate of convergence of the Newton iteration,
!                  then S = RATE/(1.D0-RATE).
!
!                  The closer RATE is to 0., the faster the Newton
!                  iteration is converging; the closer RATE is to 1.,
!                  the slower the Newton iteration is converging.
!
!                  On the first Newton iteration with an up-dated
!                  preconditioner S = 100.D0, Thus the initial
!                  RATE of convergence is approximately 1.
!
!                  S is preserved from call to call so that the rate
!                  estimate from a previous step can be applied to
!                  the current step.
!     UROUND    -- Unit roundoff.
!     EPLI      -- convergence test constant.
!                  See DDASPK prologue for more details.
!     SQRTN     -- Square root of NEQ.
!     RSQRTN    -- reciprical of square root of NEQ.
!     EPCON     -- Tolerance to test for convergence of the Newton
!                  iteration.
!     JCALC     -- Flag used to determine when to update
!                  the Jacobian matrix.  In general:
!
!                  JCALC = -1 ==> Call the DITMD routine to update
!                                 the Jacobian matrix.
!                  JCALC =  0 ==> Jacobian matrix is up-to-date.
!                  JCALC =  1 ==> Jacobian matrix is out-dated,
!                                 but DITMD will not be called unless
!                                 JCALC is set to -1.
!     JFLG      -- Flag showing whether a Jacobian routine is supplied.
!     KP1       -- The current order + 1;  updated across calls.
!     NONNEG    -- Flag to determine nonnegativity constraints.
!     NTYPE     -- Identification code for the DNEDK routine.
!                   1 ==> modified Newton; iterative linear solver.
!                   2 ==> modified Newton; user-supplied linear solver.
!     IERNLS    -- Error flag for nonlinear solver.
!                   0 ==> nonlinear solver converged.
!                   1 ==> recoverable error inside non-linear solver.
!                  -1 ==> unrecoverable error inside non-linear solver.
!
!     The following group of variables are passed as arguments to
!     the Newton iteration solver.  They are explained in greater detail
!     in DNSK:
!        TOLNEW, MULDEL, MAXIT, IERNEW
!
!     IERTYP -- Flag which tells whether this subroutine is correct.
!               0 ==> correct subroutine.
!               1 ==> incorrect subroutine.
!
!-----------------------------------------------------------------------
!***ROUTINES CALLED
!   RES, JACK, DDWNRM, DNSK
!
!***END PROLOGUE  DNEDK
!
!
IMPLICIT DOUBLE PRECISION(A-H,O-Z)
DIMENSION Y(*),YPRIME(*),WT(*)
DIMENSION PHI(NEQ,*),SAVR(*),DELTA(*),E(*)
DIMENSION WM(*),IWM(*)
DIMENSION GAMMA(*),RPAR(*),IPAR(*)
EXTERNAL  RES, JACK, PSOL
!
INTEGER, PARAMETER :: LNRE=12, LNJE=13, LLOCWP=29, LLCIWP=30
!
SAVE MULDEL, MAXIT, XRATE
DATA MULDEL/0/, MAXIT/4/, XRATE/0.25D0/
!
!     Verify that this is the correct subroutine.
!
IRES = 0    ! Karline-added that to avoid warning
IERTYP = 0
IF (NTYPE .NE. 1) THEN
   IERTYP = 1
   GO TO 380
ENDIF
!
!     If this is the first step, perform initializations.
!
IF (JSTART .EQ. 0) THEN
   CJOLD = CJ
   JCALC = -1
   S = 100.D0
   ENDIF
!
!     Perform all other initializations.
!
IERNLS = 0
LWP = IWM(LLOCWP)
LIWP = IWM(LLCIWP)
!
!     Decide whether to update the preconditioner.
!
IF (JFLG .NE. 0) THEN
   TEMP1 = (1.0D0 - XRATE)/(1.0D0 + XRATE)
   TEMP2 = 1.0D0/TEMP1
   IF (CJ/CJOLD .LT. TEMP1 .OR. CJ/CJOLD .GT. TEMP2) JCALC = -1
   IF (CJ .NE. CJLAST) S = 100.D0
ELSE
   JCALC = 0
   ENDIF
!
!     Looping point for updating preconditioner with current stepsize.
!
300 CONTINUE
!
!     Initialize all error flags to zero.
!
IERPJ = 0
IRES = 0
IERSL = 0
IERNEW = 0
!
!     Predict the solution and derivative and compute the tolerance
!     for the Newton iteration.
!
DO 310 I=1,NEQ
   Y(I)=PHI(I,1)
   YPRIME(I)=0.0D0
310 CONTINUE
DO 330 J=2,KP1
   DO 320 I=1,NEQ
      Y(I)=Y(I)+PHI(I,J)
      YPRIME(I)=YPRIME(I)+GAMMA(J)*PHI(I,J)
320 CONTINUE
330 CONTINUE
EPLIN = EPLI*EPCON
TOLNEW = EPLIN
!
!     Call RES to initialize DELTA.
!
IWM(LNRE)=IWM(LNRE)+1
CALL RES(X,Y,YPRIME,CJ,DELTA,IRES,RPAR,IPAR)
IF (IRES .LT. 0) GO TO 380
!
!
!     If indicated, update the preconditioner.
!     Set JCALC to 0 as an indicator that this has been done.
!
IF(JCALC .EQ. -1)THEN
   IWM(LNJE) = IWM(LNJE) + 1
   JCALC=0
   CALL JACK (RES, IRES, NEQ, X, Y, YPRIME, WT, DELTA, E, H, CJ, &
& WM(LWP), IWM(LIWP), IERPJ, RPAR, IPAR)
   CJOLD=CJ
   S = 100.D0
   IF (IRES .LT. 0)  GO TO 380
   IF (IERPJ .NE. 0) GO TO 380
ENDIF
!
!     Call the nonlinear Newton solver.
!
CALL DNSK(X,Y,YPRIME,NEQ,RES,PSOL,WT,RPAR,IPAR,SAVR, &
& DELTA,E,WM,IWM,CJ,SQRTN,RSQRTN,EPLIN,EPCON, &
& S,TEMP1,TOLNEW,MULDEL,MAXIT,IRES,IERSL,IERNEW)
!
IF (IERNEW .GT. 0 .AND. JCALC .NE. 0) THEN
!
!     The Newton iteration had a recoverable failure with an old
!     preconditioner.  Retry the step with a new preconditioner.
!
   JCALC = -1
   GO TO 300
ENDIF
!
IF (IERNEW .NE. 0) GO TO 380
!
!     The Newton iteration has converged.  If nonnegativity of
!     solution is required, set the solution nonnegative, if the
!     perturbation to do it is small enough.  If the change is too
!     large, then consider the corrector iteration to have failed.
!
IF(NONNEG .EQ. 0) GO TO 390
DO 360 I = 1,NEQ
  DELTA(I) = MIN(Y(I),0.0D0)
360 CONTINUE
DELNRM = DDWNRM(NEQ,DELTA,WT,RPAR,IPAR)
IF(DELNRM .GT. EPCON) GO TO 380
DO 370 I = 1,NEQ
  E(I) = E(I) - DELTA(I)
370 CONTINUE
GO TO 390
!
!
!     Exits from nonlinear solver.
!     No convergence with current preconditioner.
!     Compute IERNLS and IDID accordingly.
!
380 CONTINUE
IF (IRES .LE. -2 .OR. IERSL .LT. 0 .OR. IERTYP .NE. 0) THEN
   IERNLS = -1
   IF (IRES .LE. -2) IDID = -11
   IF (IERSL .LT. 0) IDID = -13
   IF (IERTYP .NE. 0) IDID = -15
ELSE
   IERNLS =  1
   IF (IRES .EQ. -1) IDID = -10
   IF (IERPJ .NE. 0) IDID = -5
   IF (IERSL .GT. 0) IDID = -14
ENDIF
!
!
390 JCALC = 1
RETURN
!
!------END OF SUBROUTINE DNEDK------------------------------------------
END
! Work performed under the auspices of the U.S. Department of Energy
! by Lawrence Livermore National Laboratory under contract number 
! W-7405-Eng-48.
!
SUBROUTINE DNSK(X,Y,YPRIME,NEQ,RES,PSOL,WT,RPAR,IPAR, &
& SAVR,DELTA,E,WM,IWM,CJ,SQRTN,RSQRTN,EPLIN,EPCON, &
& S,CONFAC,TOLNEW,MULDEL,MAXIT,IRES,IERSL,IERNEW)
!
!***BEGIN PROLOGUE  DNSK
!***REFER TO  DDASPK
!***DATE WRITTEN   891219   (YYMMDD)
!***REVISION DATE  900926   (YYMMDD)
!***REVISION DATE  950126   (YYMMDD)
!
!
!-----------------------------------------------------------------------
!***DESCRIPTION
!
!     DNSK solves a nonlinear system of
!     algebraic equations of the form
!     G(X,Y,YPRIME) = 0 for the unknown Y.
!
!     The method used is a modified Newton scheme.
!
!     The parameters represent
!
!     X         -- Independent variable.
!     Y         -- Solution vector.
!     YPRIME    -- Derivative of solution vector.
!     NEQ       -- Number of unknowns.
!     RES       -- External user-supplied subroutine
!                  to evaluate the residual.  See RES description
!                  in DDASPK prologue.
!     PSOL      -- External user-supplied routine to solve
!                  a linear system using preconditioning. 
!                  See explanation inside DDASPK.
!     WT        -- Vector of weights for error criterion.
!     RPAR,IPAR -- Real and integer arrays used for communication
!                  between the calling program and external user
!                  routines.  They are not altered within DASPK.
!     SAVR      -- Work vector for DNSK of length NEQ.
!     DELTA     -- Work vector for DNSK of length NEQ.
!     E         -- Error accumulation vector for DNSK of length NEQ.
!     WM,IWM    -- Real and integer arrays storing
!                  matrix information such as the matrix
!                  of partial derivatives, permutation
!                  vector, and various other information.
!     CJ        -- Parameter always proportional to 1/H (step size).
!     SQRTN     -- Square root of NEQ.
!     RSQRTN    -- reciprical of square root of NEQ.
!     EPLIN     -- Tolerance for linear system solver.
!     EPCON     -- Tolerance to test for convergence of the Newton
!                  iteration.
!     S         -- Used for error convergence tests.
!                  In the Newton iteration: S = RATE/(1.D0-RATE),
!                  where RATE is the estimated rate of convergence
!                  of the Newton iteration.
!
!                  The closer RATE is to 0., the faster the Newton
!                  iteration is converging; the closer RATE is to 1.,
!                  the slower the Newton iteration is converging.
!
!                  The calling routine sends the initial value
!                  of S to the Newton iteration.
!     CONFAC    -- A residual scale factor to improve convergence.
!     TOLNEW    -- Tolerance on the norm of Newton correction in
!                  alternative Newton convergence test.
!     MULDEL    -- A flag indicating whether or not to multiply
!                  DELTA by CONFAC.
!                  0  ==> do not scale DELTA by CONFAC.
!                  1  ==> scale DELTA by CONFAC.
!     MAXIT     -- Maximum allowed number of Newton iterations.
!     IRES      -- Error flag returned from RES.  See RES description
!                  in DDASPK prologue.  If IRES = -1, then IERNEW
!                  will be set to 1.
!                  If IRES < -1, then IERNEW will be set to -1.
!     IERSL     -- Error flag for linear system solver.
!                  See IERSL description in subroutine DSLVK.
!                  If IERSL = 1, then IERNEW will be set to 1.
!                  If IERSL < 0, then IERNEW will be set to -1.
!     IERNEW    -- Error flag for Newton iteration.
!                   0  ==> Newton iteration converged.
!                   1  ==> recoverable error inside Newton iteration.
!                  -1  ==> unrecoverable error inside Newton iteration.
!-----------------------------------------------------------------------
!
!***ROUTINES CALLED
!   RES, DSLVK, DDWNRM
!
!***END PROLOGUE  DNSK
!
!
IMPLICIT DOUBLE PRECISION(A-H,O-Z)
DIMENSION Y(*),YPRIME(*),WT(*),DELTA(*),E(*),SAVR(*)
DIMENSION WM(*),IWM(*), RPAR(*),IPAR(*)
EXTERNAL  RES, PSOL
!
INTEGER, PARAMETER :: LNNI=19, LNRE=12
!
!     Initialize Newton counter M and accumulation vector E.
!
M = 0
DO 100 I=1,NEQ
  E(I) = 0.0D0
100 CONTINUE
!
!     Corrector loop.
!
300 CONTINUE
IWM(LNNI) = IWM(LNNI) + 1
!
!     If necessary, multiply residual by convergence factor.
!
IF (MULDEL .EQ. 1) THEN
  DO 320 I = 1,NEQ
    DELTA(I) = DELTA(I) * CONFAC
320 CONTINUE
ENDIF
!
!     Save residual in SAVR.
!
DO 340 I = 1,NEQ
  SAVR(I) = DELTA(I)
340 CONTINUE
!
!     Compute a new iterate.  Store the correction in DELTA.
!
CALL DSLVK (NEQ, Y, X, YPRIME, SAVR, DELTA, WT, WM, IWM, &
& RES, IRES, PSOL, IERSL, CJ, EPLIN, SQRTN, RSQRTN, RHOK, &
& RPAR, IPAR)
IF (IRES .NE. 0 .OR. IERSL .NE. 0) GO TO 380
!
!     Update Y, E, and YPRIME.
!
DO 360 I=1,NEQ
   Y(I) = Y(I) - DELTA(I)
   E(I) = E(I) - DELTA(I)
   YPRIME(I) = YPRIME(I) - CJ*DELTA(I)
360 CONTINUE
!
!     Test for convergence of the iteration.
!
DELNRM = DDWNRM(NEQ,DELTA,WT,RPAR,IPAR)
IF (DELNRM .LE. TOLNEW) GO TO 370
IF (M .EQ. 0) THEN
  OLDNRM = DELNRM
ELSE
  RATE = (DELNRM/OLDNRM)**(1.0D0/M)
  IF (RATE .GT. 0.9D0) GO TO 380
  S = RATE/(1.0D0 - RATE)
ENDIF
IF (S*DELNRM .LE. EPCON) GO TO 370
!
!     The corrector has not yet converged.  Update M and test whether
!     the maximum number of iterations have been tried.
!
M = M + 1
IF (M .GE. MAXIT) GO TO 380
!
!     Evaluate the residual, and go back to do another iteration.
!
IWM(LNRE) = IWM(LNRE) + 1
CALL RES(X,Y,YPRIME,CJ,DELTA,IRES,RPAR,IPAR)
IF (IRES .LT. 0) GO TO 380
GO TO 300
!
!     The iteration has converged.
!
370 RETURN
!
!     The iteration has not converged.  Set IERNEW appropriately.
!
380 CONTINUE
IF (IRES .LE. -2 .OR. IERSL .LT. 0) THEN
   IERNEW = -1
ELSE
   IERNEW = 1
ENDIF
RETURN
!
!
!------END OF SUBROUTINE DNSK-------------------------------------------
END
! Work performed under the auspices of the U.S. Department of Energy
! by Lawrence Livermore National Laboratory under contract number 
! W-7405-Eng-48.
!
SUBROUTINE DSLVK (NEQ, Y, TN, YPRIME, SAVR, X, EWT, WM, IWM, &
& RES, IRES, PSOL, IERSL, CJ, EPLIN, SQRTN, RSQRTN, RHOK, &
& RPAR, IPAR)
!
!***BEGIN PROLOGUE  DSLVK
!***REFER TO  DDASPK
!***DATE WRITTEN   890101   (YYMMDD)
!***REVISION DATE  900926   (YYMMDD)
!***REVISION DATE  940928   Removed MNEWT and added RHOK in call list.
!
!
!-----------------------------------------------------------------------
!***DESCRIPTION
!
! DSLVK uses a restart algorithm and interfaces to DSPIGM for
! the solution of the linear system arising from a Newton iteration.
!
! In addition to variables described elsewhere,
! communication with DSLVK uses the following variables..
! WM    = Real work space containing data for the algorithm
!         (Krylov basis vectors, Hessenberg matrix, etc.).
! IWM   = Integer work space containing data for the algorithm.
! X     = The right-hand side vector on input, and the solution vector
!         on output, of length NEQ.
! IRES  = Error flag from RES.
! IERSL = Output flag ..
!         IERSL =  0 means no trouble occurred (or user RES routine
!                    returned IRES < 0)
!         IERSL =  1 means the iterative method failed to converge
!                    (DSPIGM returned IFLAG > 0.)
!         IERSL = -1 means there was a nonrecoverable error in the
!                    iterative solver, and an error exit will occur.
!-----------------------------------------------------------------------
!***ROUTINES CALLED
!   DSCAL, DCOPY, DSPIGM
!
!***END PROLOGUE  DSLVK
!
INTEGER NEQ, IWM, IRES, IERSL, IPAR
DOUBLE PRECISION Y, TN, YPRIME, SAVR, X, EWT, WM, CJ, EPLIN, &
& SQRTN, RSQRTN, RHOK, RPAR
DIMENSION Y(*), YPRIME(*), SAVR(*), X(*), EWT(*), &
& WM(*), IWM(*), RPAR(*), IPAR(*)
!
INTEGER IFLAG, IRST, NRSTS, NRMAX, LR, LDL, LHES, LGMR, LQ, LV, &
& LWK, LZ, MAXLP1, NPSL
INTEGER NLI, NPS, NCFL, NRE, MAXL, KMP, MITER
EXTERNAL  RES, PSOL
!    
INTEGER, PARAMETER :: LNRE=12, LNCFL=16, LNLI=20, LNPS=21
INTEGER, PARAMETER :: LLOCWP=29, LLCIWP=30
INTEGER, PARAMETER :: LMITER=23, LMAXL=24, LKMP=25, LNRMAX=26
!
!-----------------------------------------------------------------------
! IRST is set to 1, to indicate restarting is in effect.
! NRMAX is the maximum number of restarts.
!-----------------------------------------------------------------------
DATA IRST/1/
!
LIWP = IWM(LLCIWP)
NLI = IWM(LNLI)
NPS = IWM(LNPS)
NCFL = IWM(LNCFL)
NRE = IWM(LNRE)
LWP = IWM(LLOCWP)
MAXL = IWM(LMAXL)
KMP = IWM(LKMP)
NRMAX = IWM(LNRMAX)
MITER = IWM(LMITER)
IERSL = 0
IRES = 0
!-----------------------------------------------------------------------
! Use a restarting strategy to solve the linear system
! P*X = -F.  Parse the work vector, and perform initializations.
! Note that zero is the initial guess for X.
!-----------------------------------------------------------------------
MAXLP1 = MAXL + 1
LV = 1
LR = LV + NEQ*MAXL
LHES = LR + NEQ + 1
LQ = LHES + MAXL*MAXLP1
LWK = LQ + 2*MAXL
LDL = LWK + MIN0(1,MAXL-KMP)*NEQ
LZ = LDL + NEQ
CALL DSCAL (NEQ, RSQRTN, EWT, 1)
CALL DCOPY (NEQ, X, 1, WM(LR), 1)
DO 110 I = 1,NEQ
   X(I) = 0.D0
110 CONTINUE
!-----------------------------------------------------------------------
! Top of loop for the restart algorithm.  Initial pass approximates
! X and sets up a transformed system to perform subsequent restarts
! to update X.  NRSTS is initialized to -1, because restarting
! does not occur until after the first pass.
! Update NRSTS; conditionally copy DL to R; call the DSPIGM
! algorithm to solve A*Z = R;  updated counters;  update X with
! the residual solution.
! Note:  if convergence is not achieved after NRMAX restarts,
! then the linear solver is considered to have failed.
!-----------------------------------------------------------------------
NRSTS = -1
115 CONTINUE
NRSTS = NRSTS + 1
IF (NRSTS .GT. 0) CALL DCOPY (NEQ, WM(LDL), 1, WM(LR),1)
CALL DSPIGM (NEQ, TN, Y, YPRIME, SAVR, WM(LR), EWT, MAXL, MAXLP1, &
& KMP, EPLIN, CJ, RES, IRES, NRES, PSOL, NPSL, WM(LZ), WM(LV), &
& WM(LHES), WM(LQ), LGMR, WM(LWP), IWM(LIWP), WM(LWK), &
& WM(LDL), RHOK, IFLAG, IRST, NRSTS, RPAR, IPAR)
NLI = NLI + LGMR
NPS = NPS + NPSL
NRE = NRE + NRES
DO 120 I = 1,NEQ
 X(I) = X(I) + WM(LZ+I-1)
120 CONTINUE
IF ((IFLAG .EQ. 1) .AND. (NRSTS .LT. NRMAX) .AND. (IRES .EQ. 0)) &
& GO TO 115
!-----------------------------------------------------------------------
! The restart scheme is finished.  Test IRES and IFLAG to see if
! convergence was not achieved, and set flags accordingly.
!-----------------------------------------------------------------------
IF (IRES .LT. 0) THEN
   NCFL = NCFL + 1
ELSE IF (IFLAG .NE. 0) THEN
   NCFL = NCFL + 1
   IF (IFLAG .GT. 0) IERSL = 1
   IF (IFLAG .LT. 0) IERSL = -1
ENDIF
!-----------------------------------------------------------------------
! Update IWM with counters, rescale EWT, and return.
!-----------------------------------------------------------------------
IWM(LNLI)  = NLI
IWM(LNPS)  = NPS
IWM(LNCFL) = NCFL
IWM(LNRE)  = NRE
CALL DSCAL (NEQ, SQRTN, EWT, 1)
RETURN
!
!------END OF SUBROUTINE DSLVK------------------------------------------
END
! Work performed under the auspices of the U.S. Department of Energy
! by Lawrence Livermore National Laboratory under contract number 
! W-7405-Eng-48.
!
SUBROUTINE DSPIGM (NEQ, TN, Y, YPRIME, SAVR, R, WGHT, MAXL, &
& MAXLP1, KMP, EPLIN, CJ, RES, IRES, NRE, PSOL, NPSL, Z, V, &
& HES, Q, LGMR, WP, IWP, WK, DL, RHOK, IFLAG, IRST, NRSTS, &
& RPAR, IPAR)
!
!***BEGIN PROLOGUE  DSPIGM
!***DATE WRITTEN   890101   (YYMMDD)
!***REVISION DATE  900926   (YYMMDD)
!***REVISION DATE  940927   Removed MNEWT and added RHOK in call list.
!
!
!-----------------------------------------------------------------------
!***DESCRIPTION
!
! This routine solves the linear system A * Z = R using a scaled
! preconditioned version of the generalized minimum residual method.
! An initial guess of Z = 0 is assumed.
!
!      On entry
!
!          NEQ = Problem size, passed to PSOL.
!
!           TN = Current Value of T.
!
!            Y = Array Containing current dependent variable vector.
!
!       YPRIME = Array Containing current first derivative of Y.
!
!         SAVR = Array containing current value of G(T,Y,YPRIME).
!
!            R = The right hand side of the system A*Z = R.
!                R is also used as work space when computing
!                the final approximation and will therefore be
!                destroyed.
!                (R is the same as V(*,MAXL+1) in the call to DSPIGM.)
!
!         WGHT = The vector of length NEQ containing the nonzero
!                elements of the diagonal scaling matrix.
!
!         MAXL = The maximum allowable order of the matrix H.
!
!       MAXLP1 = MAXL + 1, used for dynamic dimensioning of HES.
!
!          KMP = The number of previous vectors the new vector, VNEW,
!                must be made orthogonal to.  (KMP .LE. MAXL.)
!
!        EPLIN = Tolerance on residuals R-A*Z in weighted rms norm.
!
!           CJ = Scalar proportional to current value of 
!                1/(step size H).
!
!           WK = Real work array used by routine DATVPK and PSOL.
!
!           DL = Real work array used for calculation of the residual
!                norm RHO when the method is incomplete (KMP.LT.MAXL)
!                and/or when using restarting.
!
!           WP = Real work array used by preconditioner PSOL.
!
!          IWP = Integer work array used by preconditioner PSOL.
!
!         IRST = Method flag indicating if restarting is being
!                performed.  IRST .GT. 0 means restarting is active,
!                while IRST = 0 means restarting is not being used.
!
!        NRSTS = Counter for the number of restarts on the current
!                call to DSPIGM.  If NRSTS .GT. 0, then the residual
!                R is already scaled, and so scaling of R is not
!                necessary.
!
!
!      On Return
!
!         Z    = The final computed approximation to the solution
!                of the system A*Z = R.
!
!         LGMR = The number of iterations performed and
!                the current order of the upper Hessenberg
!                matrix HES.
!
!         NRE  = The number of calls to RES (i.e. DATVPK)
!
!         NPSL = The number of calls to PSOL.
!
!         V    = The neq by (LGMR+1) array containing the LGMR
!                orthogonal vectors V(*,1) to V(*,LGMR).
!
!         HES  = The upper triangular factor of the QR decomposition
!                of the (LGMR+1) by LGMR upper Hessenberg matrix whose
!                entries are the scaled inner-products of A*V(*,I)
!                and V(*,K).
!
!         Q    = Real array of length 2*MAXL containing the components
!                of the givens rotations used in the QR decomposition
!                of HES.  It is loaded in DHEQR and used in DHELS.
!
!         IRES = Error flag from RES.
!
!           DL = Scaled preconditioned residual, 
!                (D-inverse)*(P-inverse)*(R-A*Z). Only loaded when
!                performing restarts of the Krylov iteration.
!
!         RHOK = Weighted norm of final preconditioned residual.
!
!        IFLAG = Integer error flag..
!                0 Means convergence in LGMR iterations, LGMR.LE.MAXL.
!                1 Means the convergence test did not pass in MAXL
!                  iterations, but the new residual norm (RHO) is
!                  .LT. the old residual norm (RNRM), and so Z is
!                  computed.
!                2 Means the convergence test did not pass in MAXL
!                  iterations, new residual norm (RHO) .GE. old residual
!                  norm (RNRM), and the initial guess, Z = 0, is
!                  returned.
!                3 Means there was a recoverable error in PSOL
!                  caused by the preconditioner being out of date.
!               -1 Means there was an unrecoverable error in PSOL.
!
!-----------------------------------------------------------------------
!***ROUTINES CALLED
!   PSOL, DNRM2, DSCAL, DATVPK, DORTH, DHEQR, DCOPY, DHELS, DAXPY
!
!***END PROLOGUE  DSPIGM
!
INTEGER NEQ,MAXL,MAXLP1,KMP,IRES,NRE,NPSL,LGMR,IWP, &
& IFLAG,IRST,NRSTS,IPAR
DOUBLE PRECISION TN,Y,YPRIME,SAVR,R,WGHT,EPLIN,CJ,Z,V,HES,Q,WP,WK, &
& DL,RHOK,RPAR
DIMENSION Y(*), YPRIME(*), SAVR(*), R(*), WGHT(*), Z(*), &
& V(NEQ,*), HES(MAXLP1,*), Q(*), WP(*), IWP(*), WK(*), DL(*), &
& RPAR(*), IPAR(*)
INTEGER I, IER, INFO, IP1, I2, J, K, LL, LLP1
DOUBLE PRECISION RNRM,C,DLNRM,PROD,RHO,S,SNORMW,DNRM2,TEM
EXTERNAL  RES, PSOL
!
! KARLINE: INITIALISED RHO TO AVOID A WARNING - should have no effect
RHO = 0.D0
IER = 0
IFLAG = 0
LGMR = 0
NPSL = 0
NRE = 0
!-----------------------------------------------------------------------
! The initial guess for Z is 0.  The initial residual is therefore
! the vector R.  Initialize Z to 0.
!-----------------------------------------------------------------------
DO 10 I = 1,NEQ
  Z(I) = 0.0D0
10 CONTINUE
!-----------------------------------------------------------------------
! Apply inverse of left preconditioner to vector R if NRSTS .EQ. 0.
! Form V(*,1), the scaled preconditioned right hand side.
!-----------------------------------------------------------------------
IF (NRSTS .EQ. 0) THEN
   CALL PSOL (NEQ, TN, Y, YPRIME, SAVR, WK, CJ, WGHT, WP, IWP, &
& R, EPLIN, IER, RPAR, IPAR)
   NPSL = 1
   IF (IER .NE. 0) GO TO 300
   DO 30 I = 1,NEQ
      V(I,1) = R(I)*WGHT(I)
30 CONTINUE
ELSE
   DO 35 I = 1,NEQ
     V(I,1) = R(I)
35 CONTINUE
ENDIF
!-----------------------------------------------------------------------
! Calculate norm of scaled vector V(*,1) and normalize it
! If, however, the norm of V(*,1) (i.e. the norm of the preconditioned
! residual) is .le. EPLIN, then return with Z=0.
!-----------------------------------------------------------------------
RNRM = DNRM2 (NEQ, V, 1)
IF (RNRM .LE. EPLIN) THEN
  RHOK = RNRM
  RETURN
  ENDIF
TEM = 1.0D0/RNRM
CALL DSCAL (NEQ, TEM, V(1,1), 1)
!-----------------------------------------------------------------------
! Zero out the HES array.
!-----------------------------------------------------------------------
DO 65 J = 1,MAXL
  DO 60 I = 1,MAXLP1
    HES(I,J) = 0.0D0
60 CONTINUE
65 CONTINUE
!-----------------------------------------------------------------------
! Main loop to compute the vectors V(*,2) to V(*,MAXL).
! The running product PROD is needed for the convergence test.
!-----------------------------------------------------------------------
PROD = 1.0D0
DO 90 LL = 1,MAXL
  LGMR = LL
!-----------------------------------------------------------------------
! Call routine DATVPK to compute VNEW = ABAR*V(LL), where ABAR is
! the matrix A with scaling and inverse preconditioner factors applied.
! Call routine DORTH to orthogonalize the new vector VNEW = V(*,LL+1).
! call routine DHEQR to update the factors of HES.
!-----------------------------------------------------------------------
  CALL DATVPK (NEQ, Y, TN, YPRIME, SAVR, V(1,LL), WGHT, Z, &
& RES, IRES, PSOL, V(1,LL+1), WK, WP, IWP, CJ, EPLIN, &
& IER, NRE, NPSL, RPAR, IPAR)
  IF (IRES .LT. 0) RETURN
  IF (IER .NE. 0) GO TO 300
  CALL DORTH (V(1,LL+1), V, HES, NEQ, LL, MAXLP1, KMP, SNORMW)
  HES(LL+1,LL) = SNORMW
  CALL DHEQR (HES, MAXLP1, LL, Q, INFO, LL)
  IF (INFO .EQ. LL) GO TO 120
!-----------------------------------------------------------------------
! Update RHO, the estimate of the norm of the residual R - A*ZL.
! If KMP .LT. MAXL, then the vectors V(*,1),...,V(*,LL+1) are not
! necessarily orthogonal for LL .GT. KMP.  The vector DL must then
! be computed, and its norm used in the calculation of RHO.
!-----------------------------------------------------------------------
  PROD = PROD*Q(2*LL)
  RHO = ABS(PROD*RNRM)
  IF ((LL.GT.KMP) .AND. (KMP.LT.MAXL)) THEN
    IF (LL .EQ. KMP+1) THEN
      CALL DCOPY (NEQ, V(1,1), 1, DL, 1)
      DO 75 I = 1,KMP
        IP1 = I + 1
        I2 = I*2
        S = Q(I2)
        C = Q(I2-1)
        DO 70 K = 1,NEQ
          DL(K) = S*DL(K) + C*V(K,IP1)
70 CONTINUE
75 CONTINUE
      ENDIF
    S = Q(2*LL)
    C = Q(2*LL-1)/SNORMW
    LLP1 = LL + 1
    DO 80 K = 1,NEQ
      DL(K) = S*DL(K) + C*V(K,LLP1)
80 CONTINUE
    DLNRM = DNRM2 (NEQ, DL, 1)
    RHO = RHO*DLNRM
    ENDIF
!-----------------------------------------------------------------------
! Test for convergence.  If passed, compute approximation ZL.
! If failed and LL .LT. MAXL, then continue iterating.
!-----------------------------------------------------------------------
  IF (RHO .LE. EPLIN) GO TO 200
  IF (LL .EQ. MAXL) GO TO 100
!-----------------------------------------------------------------------
! Rescale so that the norm of V(1,LL+1) is one.
!-----------------------------------------------------------------------
  TEM = 1.0D0/SNORMW
  CALL DSCAL (NEQ, TEM, V(1,LL+1), 1)
90 CONTINUE
100 CONTINUE
IF (RHO .LT. RNRM) GO TO 150
120 CONTINUE
IFLAG = 2
DO 130 I = 1,NEQ
   Z(I) = 0.D0
130 CONTINUE
RETURN
150 IFLAG = 1
!-----------------------------------------------------------------------
! The tolerance was not met, but the residual norm was reduced.
! If performing restarting (IRST .gt. 0) calculate the residual vector
! RL and store it in the DL array.  If the incomplete version is 
! being used (KMP .lt. MAXL) then DL has already been calculated.
!-----------------------------------------------------------------------
IF (IRST .GT. 0) THEN
   IF (KMP .EQ. MAXL) THEN
!
!           Calculate DL from the V(I)'s.
!
      CALL DCOPY (NEQ, V(1,1), 1, DL, 1)
      MAXLM1 = MAXL - 1
      DO 175 I = 1,MAXLM1
         IP1 = I + 1
         I2 = I*2
         S = Q(I2)
         C = Q(I2-1)
         DO 170 K = 1,NEQ
            DL(K) = S*DL(K) + C*V(K,IP1)
170 CONTINUE
175 CONTINUE
      S = Q(2*MAXL)
      C = Q(2*MAXL-1)/SNORMW
      DO 180 K = 1,NEQ
         DL(K) = S*DL(K) + C*V(K,MAXLP1)
180 CONTINUE
   ENDIF
!
!        Scale DL by RNRM*PROD to obtain the residual RL.
!
   TEM = RNRM*PROD
   CALL DSCAL(NEQ, TEM, DL, 1)
ENDIF
!-----------------------------------------------------------------------
! Compute the approximation ZL to the solution.
! Since the vector Z was used as work space, and the initial guess
! of the Newton correction is zero, Z must be reset to zero.
!-----------------------------------------------------------------------
200 CONTINUE
LL = LGMR
LLP1 = LL + 1
DO 210 K = 1,LLP1
  R(K) = 0.0D0
210 CONTINUE
R(1) = RNRM
CALL DHELS (HES, MAXLP1, LL, Q, R)
DO 220 K = 1,NEQ
  Z(K) = 0.0D0
220 CONTINUE
DO 230 I = 1,LL
  CALL DAXPY (NEQ, R(I), V(1,I), 1, Z, 1)
230 CONTINUE
DO 240 I = 1,NEQ
  Z(I) = Z(I)/WGHT(I)
240 CONTINUE
! Load RHO into RHOK.
RHOK = RHO
RETURN
!-----------------------------------------------------------------------
! This block handles error returns forced by routine PSOL.
!-----------------------------------------------------------------------
300 CONTINUE
IF (IER .LT. 0) IFLAG = -1
IF (IER .GT. 0) IFLAG = 3
!
RETURN
!
!------END OF SUBROUTINE DSPIGM-----------------------------------------
END
! Work performed under the auspices of the U.S. Department of Energy
! by Lawrence Livermore National Laboratory under contract number 
! W-7405-Eng-48.
!
SUBROUTINE DORTH (VNEW, V, HES, N, LL, LDHES, KMP, SNORMW)
!
!***BEGIN PROLOGUE  DORTH
!***DATE WRITTEN   890101   (YYMMDD)
!***REVISION DATE  900926   (YYMMDD)
!
!
!-----------------------------------------------------------------------
!***DESCRIPTION
!
! This routine orthogonalizes the vector VNEW against the previous
! KMP vectors in the V array.  It uses a modified Gram-Schmidt
! orthogonalization procedure with conditional reorthogonalization.
!
!      On entry
!
!         VNEW = The vector of length N containing a scaled product
!                OF The Jacobian and the vector V(*,LL).
!
!         V    = The N x LL array containing the previous LL
!                orthogonal vectors V(*,1) to V(*,LL).
!
!         HES  = An LL x LL upper Hessenberg matrix containing,
!                in HES(I,K), K.LT.LL, scaled inner products of
!                A*V(*,K) and V(*,I).
!
!        LDHES = The leading dimension of the HES array.
!
!         N    = The order of the matrix A, and the length of VNEW.
!
!         LL   = The current order of the matrix HES.
!
!          KMP = The number of previous vectors the new vector VNEW
!                must be made orthogonal to (KMP .LE. MAXL).
!
!
!      On return
!
!         VNEW = The new vector orthogonal to V(*,I0),
!                where I0 = MAX(1, LL-KMP+1).
!
!         HES  = Upper Hessenberg matrix with column LL filled in with
!                scaled inner products of A*V(*,LL) and V(*,I).
!
!       SNORMW = L-2 norm of VNEW.
!
!-----------------------------------------------------------------------
!***ROUTINES CALLED
!   DDOT, DNRM2, DAXPY 
!
!***END PROLOGUE  DORTH
!
INTEGER N, LL, LDHES, KMP
DOUBLE PRECISION VNEW, V, HES, SNORMW
DIMENSION VNEW(*), V(N,*), HES(LDHES,*)
INTEGER I, I0
DOUBLE PRECISION ARG, DDOT, DNRM2, SUMDSQ, TEM, VNRM
!
!-----------------------------------------------------------------------
! Get norm of unaltered VNEW for later use.
!-----------------------------------------------------------------------
VNRM = DNRM2 (N, VNEW, 1)
!-----------------------------------------------------------------------
! Do Modified Gram-Schmidt on VNEW = A*V(LL).
! Scaled inner products give new column of HES.
! Projections of earlier vectors are subtracted from VNEW.
!-----------------------------------------------------------------------
I0 = MAX0(1,LL-KMP+1)
DO 10 I = I0,LL
  HES(I,LL) = DDOT (N, V(1,I), 1, VNEW, 1)
  TEM = -HES(I,LL)
  CALL DAXPY (N, TEM, V(1,I), 1, VNEW, 1)
10 CONTINUE
!-----------------------------------------------------------------------
! Compute SNORMW = norm of VNEW.
! If VNEW is small compared to its input value (in norm), then
! Reorthogonalize VNEW to V(*,1) through V(*,LL).
! Correct if relative correction exceeds 1000*(unit roundoff).
! Finally, correct SNORMW using the dot products involved.
!-----------------------------------------------------------------------
SNORMW = DNRM2 (N, VNEW, 1)
IF (VNRM + 0.001D0*SNORMW .NE. VNRM) RETURN
SUMDSQ = 0.0D0
DO 30 I = I0,LL
  TEM = -DDOT (N, V(1,I), 1, VNEW, 1)
  IF (HES(I,LL) + 0.001D0*TEM .EQ. HES(I,LL)) GO TO 30
  HES(I,LL) = HES(I,LL) - TEM
  CALL DAXPY (N, TEM, V(1,I), 1, VNEW, 1)
  SUMDSQ = SUMDSQ + TEM**2
30 CONTINUE
IF (SUMDSQ .EQ. 0.0D0) RETURN
ARG = MAX(0.0D0,SNORMW**2 - SUMDSQ)
SNORMW = SQRT(ARG)
RETURN
!
!------END OF SUBROUTINE DORTH------------------------------------------
END

!-----------------------------------------------------------------------
! Karline:
! rescaling of error term according to the index of each variable
! index 2 variables are scaled with 1/H
! index 3 variables are scaled with 1/H^2
!-----------------------------------------------------------------------

SUBROUTINE SCALE(NEQ, NIND, SCAL, H)
INTEGER NEQ, NIND(3) , I
DOUBLE PRECISION SCAL(*), H
IF(NIND(2).NE.0) THEN
       DO I=NIND(1)+1,NIND(1)+NIND(2)
          SCAL(I)=SCAL(I)/min(1.D0,H)
       END DO
ENDIF

IF(NIND(3).NE.0) THEN
      DO I=NIND(1)+NIND(2)+1,NIND(1)+NIND(2)+NIND(3)
          SCAL(I)=SCAL(I)/min(1.D0, H*H)
       END DO
ENDIF
RETURN
END

SUBROUTINE DATVPK(NEQ, Y, TN, YPRIME, SAVR, V, WGHT, YPTEM, RES, &
& IRES, PSOL, Z, VTEM, WP, IWP, CJ, EPLIN, IER, NRE, NPSL, &
& RPAR,IPAR)
!
!***BEGIN PROLOGUE  DATV
!***DATE WRITTEN   890101   (YYMMDD)
!***REVISION DATE  900926   (YYMMDD)
!KS  30-06-2019: name changed to DATVPK 
!
!-----------------------------------------------------------------------
!***DESCRIPTION
!
! This routine computes the product
!
!   Z = (D-inverse)*(P-inverse)*(dF/dY)*(D*V),
!
! where F(Y) = G(T, Y, CJ*(Y-A)), CJ is a scalar proportional to 1/H,
! and A involves the past history of Y.  The quantity CJ*(Y-A) is
! an approximation to the first derivative of Y and is stored
! in the array YPRIME.  Note that dF/dY = dG/dY + CJ*dG/dYPRIME.
!
! D is a diagonal scaling matrix, and P is the left preconditioning
! matrix.  V is assumed to have L2 norm equal to 1.
! The product is stored in Z and is computed by means of a
! difference quotient, a call to RES, and one call to PSOL.
!
!      On entry
!
!          NEQ = Problem size, passed to RES and PSOL.
!
!            Y = Array containing current dependent variable vector.
!
!       YPRIME = Array containing current first derivative of y.
!
!         SAVR = Array containing current value of G(T,Y,YPRIME).
!
!            V = Real array of length NEQ (can be the same array as Z).
!
!         WGHT = Array of length NEQ containing scale factors.
!                1/WGHT(I) are the diagonal elements of the matrix D.
!
!        YPTEM = Work array of length NEQ.
!
!         VTEM = Work array of length NEQ used to store the
!                unscaled version of V.
!
!         WP = Real work array used by preconditioner PSOL.
!
!         IWP = Integer work array used by preconditioner PSOL.
!
!           CJ = Scalar proportional to current value of 
!                1/(step size H).
!
!
!      On return
!
!            Z = Array of length NEQ containing desired scaled
!                matrix-vector product.
!
!         IRES = Error flag from RES.
!
!          IER = Error flag from PSOL.
!
!         NRE  = The number of calls to RES.
!
!         NPSL = The number of calls to PSOL.
!
!-----------------------------------------------------------------------
!***ROUTINES CALLED
!   RES, PSOL
!
!***END PROLOGUE  DATVPK
!
INTEGER NEQ, IRES, IWP, IER, NRE, NPSL, IPAR
DOUBLE PRECISION Y, TN, YPRIME, SAVR, V, WGHT, YPTEM, Z, VTEM, &
& WP, CJ, RPAR
DIMENSION Y(*), YPRIME(*), SAVR(*), V(*), WGHT(*), YPTEM(*), &
& Z(*), VTEM(*), WP(*), IWP(*), RPAR(*), IPAR(*)
INTEGER I
DOUBLE PRECISION EPLIN
EXTERNAL  RES, PSOL
!
IRES = 0
!-----------------------------------------------------------------------
! Set VTEM = D * V.
!-----------------------------------------------------------------------
DO 10 I = 1,NEQ
  VTEM(I) = V(I)/WGHT(I)
10 CONTINUE
IER = 0
!-----------------------------------------------------------------------
! Store Y in Z and increment Z by VTEM.
! Store YPRIME in YPTEM and increment YPTEM by VTEM*CJ.
!-----------------------------------------------------------------------
DO 20 I = 1,NEQ
  YPTEM(I) = YPRIME(I) + VTEM(I)*CJ
  Z(I) = Y(I) + VTEM(I)
20 CONTINUE
!-----------------------------------------------------------------------
! Call RES with incremented Y, YPRIME arguments
! stored in Z, YPTEM.  VTEM is overwritten with new residual.
!-----------------------------------------------------------------------
CONTINUE
CALL RES(TN,Z,YPTEM,CJ,VTEM,IRES,RPAR,IPAR)
NRE = NRE + 1
IF (IRES .LT. 0) RETURN
!-----------------------------------------------------------------------
! Set Z = (dF/dY) * VBAR using difference quotient.
! (VBAR is old value of VTEM before calling RES)
!-----------------------------------------------------------------------
DO 70 I = 1,NEQ
  Z(I) = VTEM(I) - SAVR(I)
70 CONTINUE
!-----------------------------------------------------------------------
! Apply inverse of left preconditioner to Z.
!-----------------------------------------------------------------------
CALL PSOL (NEQ, TN, Y, YPRIME, SAVR, YPTEM, CJ, WGHT, WP, IWP, &
& Z, EPLIN, IER, RPAR, IPAR)
NPSL = NPSL + 1
IF (IER .NE. 0) RETURN
!-----------------------------------------------------------------------
! Apply D-inverse to Z and return.
!-----------------------------------------------------------------------
DO 90 I = 1,NEQ
  Z(I) = Z(I)*WGHT(I)
90 CONTINUE
RETURN
!
!------END OF SUBROUTINE DATVPK-------------------------------------------
END
