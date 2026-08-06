module minqa_module
   use, intrinsic :: iso_fortran_env, only : real64, output_unit
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
   implicit none
   private

   integer, parameter, public :: dp = real64

   abstract interface
      function objective_function(x) result(f)
         import :: dp
         real(dp), intent(in) :: x(:)
         real(dp) :: f
      end function objective_function
   end interface

   type, public :: minqa_control_t
      integer :: npt = 0
      real(dp) :: rhobeg = -1.0_dp
      real(dp) :: rhoend = -1.0_dp
      integer :: iprint = 0
      integer :: maxfun = 10000
      logical :: adjust_start = .false.
   end type minqa_control_t

   type, public :: minqa_result_t
      real(dp), allocatable :: x(:)
      real(dp) :: fval = huge(1.0_dp)
      integer :: evaluations = 0
      integer :: status = 0
      integer :: raw_status = 0
      character(len=:), allocatable :: message
   end type minqa_result_t

   interface bobyqa
      module procedure bobyqa_arrays
      module procedure bobyqa_scalar_bounds
   end interface bobyqa

   public :: objective_function
   public :: bobyqa, newuoa, uobyqa
   public :: minqa_status_message

   procedure(objective_function), pointer :: active_objective => null()
   integer :: active_evaluations = 0

contains

   subroutine bobyqa_arrays(objective, x, result, lower, upper, control)
      procedure(objective_function) :: objective
      real(dp), intent(inout) :: x(:)
      type(minqa_result_t), intent(out) :: result
      real(dp), intent(in), optional :: lower(:), upper(:)
      type(minqa_control_t), intent(in), optional :: control

      type(minqa_control_t) :: ctrl
      real(dp), allocatable :: xl(:), xu(:), work(:)
      real(dp) :: rb, re, min_range
      integer :: n, npt, ierr, lwork

      n = size(x)
      if (present(control)) ctrl = control
      call initialize_result(result, x)

      if (n < 2) then
         call set_input_error(result, 'minqa requires at least two variables')
         return
      end if
      if (ctrl%maxfun <= 0) then
         call set_input_error(result, 'maxfun must be positive')
         return
      end if

      allocate(xl(n), xu(n))
      xl = -huge(1.0_dp) / 4.0_dp
      xu =  huge(1.0_dp) / 4.0_dp

      if (present(lower)) then
         if (size(lower) /= n) then
            call set_input_error(result, 'lower must have the same size as x')
            return
         end if
         xl = lower
      end if
      if (present(upper)) then
         if (size(upper) /= n) then
            call set_input_error(result, 'upper must have the same size as x')
            return
         end if
         xu = upper
      end if

      if (any(.not. ieee_is_finite(x)) .or. any(ieee_is_nan(xl)) .or. any(ieee_is_nan(xu))) then
         call set_input_error(result, 'x must be finite and bounds must not be NaN')
         return
      end if
      where (.not. ieee_is_finite(xl) .and. xl < 0.0_dp)
         xl = -huge(1.0_dp) / 4.0_dp
      end where
      where (.not. ieee_is_finite(xu) .and. xu > 0.0_dp)
         xu = huge(1.0_dp) / 4.0_dp
      end where
      if (any(.not. ieee_is_finite(xl)) .or. any(.not. ieee_is_finite(xu))) then
         call set_input_error(result, 'lower may only use -Inf and upper may only use +Inf')
         return
      end if
      if (any(xl >= xu)) then
         call set_input_error(result, 'each lower bound must be less than its upper bound')
         return
      end if

      if (any(x < xl) .or. any(x > xu)) then
         if (ctrl%adjust_start) then
            x = min(max(x, xl), xu)
         else
            call set_input_error(result, 'starting values violate the bounds')
            return
         end if
      end if

      call resolve_controls(n, x, ctrl, npt, rb, re, result)
      if (result%status /= 0) return

      min_range = minval(xu - xl)
      if (min_range < 2.0_dp * rb) rb = 0.2_dp * min_range
      if (re > rb) re = min(re, 1.0e-6_dp * rb)
      if (rb <= 0.0_dp .or. re <= 0.0_dp) then
         call set_input_error(result, 'bounds do not permit a positive trust-region radius')
         return
      end if

      lwork = (npt + 5) * (npt + n) + (3 * n * (n + 5)) / 2
      allocate(work(lwork))

      active_objective => objective
      active_evaluations = 0
      ierr = 0
      call bobyqa_core(n, npt, x, xl, xu, rb, re, ctrl%iprint, ctrl%maxfun, work, ierr)
      result%evaluations = active_evaluations
      result%raw_status = ierr
      result%status = map_status(ierr)
      result%message = minqa_status_message(result%status, 'bobyqa')
      result%x = x
      result%fval = safe_objective_value(objective, x)
      nullify(active_objective)
   end subroutine bobyqa_arrays

   subroutine bobyqa_scalar_bounds(objective, x, result, lower, upper, control)
      procedure(objective_function) :: objective
      real(dp), intent(inout) :: x(:)
      type(minqa_result_t), intent(out) :: result
      real(dp), intent(in) :: lower, upper
      type(minqa_control_t), intent(in), optional :: control

      real(dp), allocatable :: xl(:), xu(:)

      allocate(xl(size(x)), xu(size(x)))
      xl = lower
      xu = upper
      if (present(control)) then
         call bobyqa_arrays(objective, x, result, xl, xu, control)
      else
         call bobyqa_arrays(objective, x, result, xl, xu)
      end if
   end subroutine bobyqa_scalar_bounds

   subroutine newuoa(objective, x, result, control)
      procedure(objective_function) :: objective
      real(dp), intent(inout) :: x(:)
      type(minqa_result_t), intent(out) :: result
      type(minqa_control_t), intent(in), optional :: control

      type(minqa_control_t) :: ctrl
      real(dp), allocatable :: work(:)
      real(dp) :: rb, re
      integer :: n, npt, ierr, lwork

      n = size(x)
      if (present(control)) ctrl = control
      call initialize_result(result, x)

      if (n < 2) then
         call set_input_error(result, 'minqa requires at least two variables')
         return
      end if
      if (any(.not. ieee_is_finite(x))) then
         call set_input_error(result, 'starting values must be finite')
         return
      end if
      if (ctrl%maxfun <= 0) then
         call set_input_error(result, 'maxfun must be positive')
         return
      end if

      call resolve_controls(n, x, ctrl, npt, rb, re, result)
      if (result%status /= 0) return

      lwork = (npt + 13) * (npt + n) + (3 * n * (n + 3)) / 2
      allocate(work(lwork))

      active_objective => objective
      active_evaluations = 0
      ierr = 0
      call newuoa_core(n, npt, x, rb, re, ctrl%iprint, ctrl%maxfun, work, ierr)
      result%evaluations = active_evaluations
      result%raw_status = ierr
      result%status = map_status(ierr)
      result%message = minqa_status_message(result%status, 'newuoa')
      result%x = x
      result%fval = safe_objective_value(objective, x)
      nullify(active_objective)
   end subroutine newuoa

   subroutine uobyqa(objective, x, result, control)
      procedure(objective_function) :: objective
      real(dp), intent(inout) :: x(:)
      type(minqa_result_t), intent(out) :: result
      type(minqa_control_t), intent(in), optional :: control

      type(minqa_control_t) :: ctrl
      real(dp), allocatable :: work(:)
      real(dp) :: rb, re
      integer :: n, npt, ierr, lwork

      n = size(x)
      if (present(control)) ctrl = control
      call initialize_result(result, x)

      if (n < 2) then
         call set_input_error(result, 'minqa requires at least two variables')
         return
      end if
      if (any(.not. ieee_is_finite(x))) then
         call set_input_error(result, 'starting values must be finite')
         return
      end if
      if (ctrl%maxfun <= 0) then
         call set_input_error(result, 'maxfun must be positive')
         return
      end if

      call resolve_controls(n, x, ctrl, npt, rb, re, result)
      if (result%status /= 0) return

      lwork = (n * (42 + n * (23 + n * (8 + n))) + max(2 * n * n + 4, 18 * n)) / 4
      allocate(work(lwork))

      active_objective => objective
      active_evaluations = 0
      ierr = 0
      call uobyqa_core(n, x, rb, re, ctrl%iprint, ctrl%maxfun, work, ierr)
      result%evaluations = active_evaluations
      result%raw_status = ierr
      result%status = map_status(ierr)
      result%message = minqa_status_message(result%status, 'uobyqa')
      result%x = x
      result%fval = safe_objective_value(objective, x)
      nullify(active_objective)
   end subroutine uobyqa

   subroutine resolve_controls(n, x, ctrl, npt, rhobeg, rhoend, result)
      integer, intent(in) :: n
      real(dp), intent(in) :: x(:)
      type(minqa_control_t), intent(in) :: ctrl
      integer, intent(out) :: npt
      real(dp), intent(out) :: rhobeg, rhoend
      type(minqa_result_t), intent(inout) :: result

      integer :: npt_max

      npt_max = ((n + 1) * (n + 2)) / 2
      if (ctrl%npt > 0) then
         npt = max(n + 2, min(ctrl%npt, npt_max))
      else
         npt = max(n + 2, min(min(n + 2, 2 * n), npt_max))
      end if

      if (ctrl%rhobeg > 0.0_dp) then
         rhobeg = ctrl%rhobeg
      else
         rhobeg = min(0.95_dp, 0.2_dp * maxval(abs(x)))
         if (rhobeg <= 0.0_dp) rhobeg = 0.5_dp
      end if

      if (ctrl%rhoend > 0.0_dp) then
         rhoend = ctrl%rhoend
      else
         rhoend = 1.0e-6_dp * rhobeg
      end if

      if (rhoend <= 0.0_dp .or. rhoend > rhobeg) then
         call set_input_error(result, 'require 0 < rhoend <= rhobeg')
      end if
   end subroutine resolve_controls

   subroutine initialize_result(result, x)
      type(minqa_result_t), intent(out) :: result
      real(dp), intent(in) :: x(:)

      result%x = x
      result%fval = huge(1.0_dp)
      result%evaluations = 0
      result%status = 0
      result%raw_status = 0
      result%message = 'not started'
   end subroutine initialize_result

   subroutine set_input_error(result, message)
      type(minqa_result_t), intent(inout) :: result
      character(len=*), intent(in) :: message

      result%status = 6
      result%raw_status = -1
      result%message = message
   end subroutine set_input_error

   integer function map_status(raw_status) result(status)
      integer, intent(in) :: raw_status

      select case (raw_status)
      case (0)
         status = 0
      case (390)
         status = 1
      case (10, 101)
         status = 2
      case (430, 3701, 2101)
         status = 3
      case (20)
         status = 4
      case (320)
         status = 5
      case default
         status = 7
      end select
   end function map_status

   function minqa_status_message(status, algorithm) result(message)
      integer, intent(in) :: status
      character(len=*), intent(in), optional :: algorithm
      character(len=:), allocatable :: message
      character(len=:), allocatable :: prefix

      if (present(algorithm)) then
         prefix = trim(algorithm) // ': '
      else
         prefix = ''
      end if

      select case (status)
      case (0)
         message = prefix // 'normal exit'
      case (1)
         message = prefix // 'maximum number of objective evaluations exceeded'
      case (2)
         message = prefix // 'npt is outside the required interval'
      case (3)
         message = prefix // 'a trust-region step failed to reduce the model'
      case (4)
         message = prefix // 'a bound interval is smaller than 2*rhobeg'
      case (5)
         message = prefix // 'excessive cancellation in an update denominator'
      case (6)
         message = prefix // 'invalid input'
      case default
         message = prefix // 'unknown optimizer status'
      end select
   end function minqa_status_message

   function safe_objective_value(objective, x) result(f)
      procedure(objective_function) :: objective
      real(dp), intent(in) :: x(:)
      real(dp) :: f

      f = objective(x)
      if (.not. ieee_is_finite(f)) f = huge(1.0_dp)
   end function safe_objective_value

   function calfun(n, x, ip) result(f)
      integer, intent(in) :: n, ip
      real(dp), intent(in) :: x(*)
      real(dp) :: f

      if (.not. associated(active_objective)) error stop 'minqa: no active objective'
      active_evaluations = active_evaluations + 1
      if (any(.not. ieee_is_finite(x(1:n)))) then
         f = huge(1.0_dp)
      else
         f = active_objective(x(1:n))
         if (.not. ieee_is_finite(f)) f = huge(1.0_dp)
      end if

      if (ip == 3) then
         write(output_unit, '(i0,": ",es16.8,":",*(1x,es13.5))') active_evaluations, f, x(1:n)
      else if (ip > 3) then
         if (mod(active_evaluations, ip) == 0) then
            write(output_unit, '(i0,": ",es16.8,":",*(1x,es13.5))') active_evaluations, f, x(1:n)
         end if
      end if
   end function calfun

   subroutine minqit(iprint, rho, nf, fopt, n, xbase, xopt)
      integer, intent(in) :: iprint, nf, n
      real(dp), intent(in) :: rho, fopt
      real(dp), intent(in) :: xbase(*), xopt(*)

      if (iprint >= 2) then
         write(output_unit, '("rho: ",es10.3," eval: ",i0," fn: ",es16.8," par:",*(1x,es13.5))') &
            rho, nf, fopt, xbase(1:n) + xopt(1:n)
      end if
   end subroutine minqit

   subroutine minqir(iprint, f, nf, n, x)
      integer, intent(in) :: iprint, nf, n
      real(dp), intent(in) :: f
      real(dp), intent(in) :: x(*)

      if (iprint > 0) then
         write(output_unit, '("At return: eval=",i0," fn=",es16.8," par:",*(1x,es13.5))') nf, f, x(1:n)
      end if
   end subroutine minqir

SUBROUTINE ALTMOV (N,NPT,XPT,XOPT,BMAT,ZMAT,NDIM,SL,SU,KOPT, KNEW,ADELT,XNEW,XALT,ALPHA,CAUCHY,GLAG,HCOL, &
    & W)
implicit integer (i-n), real(real64) (a-h,o-z)
DIMENSION XPT(NPT,*),XOPT(*),BMAT(NDIM,*),ZMAT(NPT,*),SL(*), SU(*),XNEW(*),XALT(*),GLAG(*),HCOL(*),W(*)
!
!     The arguments N, NPT, XPT, XOPT, BMAT, ZMAT, NDIM, SL and SU all have
!       the same meanings as the corresponding arguments of BOBYQB.
!     KOPT is the index of the optimal interpolation point.
!     KNEW is the index of the interpolation point that is going to be moved.
!     ADELT is the current trust region bound.
!     XNEW will be set to a suitable new position for the interpolation point
!       XPT(KNEW,.). Specifically, it satisfies the SL, SU and trust region
!       bounds and it should provide a large denominator in the next call of
!       UPDATE. The step XNEW-XOPT from XOPT is restricted to moves along the
!       straight lines through XOPT and another interpolation point.
!     XALT also provides a large value of the modulus of the KNEW-th Lagrange
!       function subject to the constraints that have been mentioned, its main
!       difference from XNEW being that XALT-XOPT is a constrained version of
!       the Cauchy step within the trust region. An exception is that XALT is
!       not calculated if all components of GLAG (see below) are zero.
!     ALPHA will be set to the KNEW-th diagonal element of the H matrix.
!     CAUCHY will be set to the square of the KNEW-th Lagrange function at
!       the step XALT-XOPT from XOPT for the vector XALT that is returned,
!       except that CAUCHY is set to zero if XALT is not calculated.
!     GLAG is a working space vector of length N for the gradient of the
!       KNEW-th Lagrange function at XOPT.
!     HCOL is a working space vector of length NPT for the second derivative
!       coefficients of the KNEW-th Lagrange function.
!     W is a working space vector of length 2N that is going to hold the
!       constrained Cauchy step from XOPT of the Lagrange function, followed
!       by the downhill version of XALT when the uphill step is calculated.
!
!     Set the first NPT components of W to the leading elements of the
!     KNEW-th column of the H matrix.
!
HALF=0.5D0
ONE=1.0D0
ZERO=0.0D0
CONST=ONE+DSQRT(2.0D0)
DO K=1,NPT
HCOL(K)=ZERO
END DO
DO J=1,NPT-N-1
TEMP=ZMAT(KNEW,J)
DO K=1,NPT
HCOL(K)=HCOL(K)+TEMP*ZMAT(K,J)
END DO
END DO
ALPHA=HCOL(KNEW)
HA=HALF*ALPHA
!
!     Calculate the gradient of the KNEW-th Lagrange function at XOPT.
!
DO  I=1,N
GLAG(I)=BMAT(KNEW,I)
END DO
DO K=1,NPT
TEMP=ZERO
DO J=1,N
TEMP=TEMP+XPT(K,J)*XOPT(J)
END DO
TEMP=HCOL(K)*TEMP
DO I=1,N
GLAG(I)=GLAG(I)+TEMP*XPT(K,I)
END DO
END DO
!
!     Search for a large denominator along the straight lines through XOPT
!     and another interpolation point. SLBD and SUBD will be lower and upper
!     bounds on the step along each of these lines in turn. PREDSQ will be
!     set to the square of the predicted denominator for each line. PRESAV
!     will be set to the largest admissible value of PREDSQ that occurs.
!
PRESAV=ZERO
DO K=1,NPT
IF (K .EQ. KOPT) GOTO 80
DDERIV=ZERO
DISTSQ=ZERO
DO I=1,N
TEMP=XPT(K,I)-XOPT(I)
DDERIV=DDERIV+GLAG(I)*TEMP
DISTSQ=DISTSQ+TEMP*TEMP
END DO
SUBD=ADELT/DSQRT(DISTSQ)
SLBD=-SUBD
ILBD=0
IUBD=0
SUMIN=DMIN1(ONE,SUBD)
!
!     Revise SLBD and SUBD if necessary because of the bounds in SL and SU.
!
DO I=1,N
TEMP=XPT(K,I)-XOPT(I)
IF (TEMP .GT. ZERO) THEN
IF (SLBD*TEMP .LT. SL(I)-XOPT(I)) THEN
SLBD=(SL(I)-XOPT(I))/TEMP
ILBD=-I
END IF
IF (SUBD*TEMP .GT. SU(I)-XOPT(I)) THEN
SUBD=DMAX1(SUMIN,(SU(I)-XOPT(I))/TEMP)
IUBD=I
END IF
ELSE IF (TEMP .LT. ZERO) THEN
IF (SLBD*TEMP .GT. SU(I)-XOPT(I)) THEN
SLBD=(SU(I)-XOPT(I))/TEMP
ILBD=I
END IF
IF (SUBD*TEMP .LT. SL(I)-XOPT(I)) THEN
SUBD=DMAX1(SUMIN,(SL(I)-XOPT(I))/TEMP)
IUBD=-I
END IF
END IF
END DO
!
!     Seek a large modulus of the KNEW-th Lagrange function when the index
!     of the other interpolation point on the line through XOPT is KNEW.
!
IF (K .EQ. KNEW) THEN
DIFF=DDERIV-ONE
STEP=SLBD
VLAG=SLBD*(DDERIV-SLBD*DIFF)
ISBD=ILBD
TEMP=SUBD*(DDERIV-SUBD*DIFF)
IF (DABS(TEMP) .GT. DABS(VLAG)) THEN
STEP=SUBD
VLAG=TEMP
ISBD=IUBD
END IF
TEMPD=HALF*DDERIV
TEMPA=TEMPD-DIFF*SLBD
TEMPB=TEMPD-DIFF*SUBD
IF (TEMPA*TEMPB .LT. ZERO) THEN
TEMP=TEMPD*TEMPD/DIFF
IF (DABS(TEMP) .GT. DABS(VLAG)) THEN
STEP=TEMPD/DIFF
VLAG=TEMP
ISBD=0
END IF
END IF
!
!     Search along each of the other lines through XOPT and another point.
!
ELSE
STEP=SLBD
VLAG=SLBD*(ONE-SLBD)
ISBD=ILBD
TEMP=SUBD*(ONE-SUBD)
IF (DABS(TEMP) .GT. DABS(VLAG)) THEN
STEP=SUBD
VLAG=TEMP
ISBD=IUBD
END IF
IF (SUBD .GT. HALF) THEN
IF (DABS(VLAG) .LT. 0.25D0) THEN
STEP=HALF
VLAG=0.25D0
ISBD=0
END IF
END IF
VLAG=VLAG*DDERIV
END IF
!
!     Calculate PREDSQ for the current line search and maintain PRESAV.
!
TEMP=STEP*(ONE-STEP)*DISTSQ
PREDSQ=VLAG*VLAG*(VLAG*VLAG+HA*TEMP*TEMP)
IF (PREDSQ .GT. PRESAV) THEN
PRESAV=PREDSQ
KSAV=K
STPSAV=STEP
IBDSAV=ISBD
END IF
80 END DO
!
!     Construct XNEW in a way that satisfies the bound constraints exactly.
!
IBDSAV=0
DO I=1,N
TEMP=XOPT(I)+STPSAV*(XPT(KSAV,I)-XOPT(I))
XNEW(I)=DMAX1(SL(I),DMIN1(SU(I),TEMP))
END DO
IF (IBDSAV .LT. 0) XNEW(-IBDSAV)=SL(-IBDSAV)
IF (IBDSAV .GT. 0) XNEW(IBDSAV)=SU(IBDSAV)
!
!     Prepare for the iterative method that assembles the constrained Cauchy
!     step in W. The sum of squares of the fixed components of W is formed in
!     WFIXSQ, and the free components of W are set to BIGSTP.
!
BIGSTP=ADELT+ADELT
IFLAG=0
100 WFIXSQ=ZERO
GGFREE=ZERO
DO I=1,N
W(I)=ZERO
TEMPA=DMIN1(XOPT(I)-SL(I),GLAG(I))
TEMPB=DMAX1(XOPT(I)-SU(I),GLAG(I))
IF (TEMPA .GT. ZERO .OR. TEMPB .LT. ZERO) THEN
W(I)=BIGSTP
GGFREE=GGFREE+GLAG(I)**2
END IF
END DO
IF (GGFREE .EQ. ZERO) THEN
CAUCHY=ZERO
GOTO 200
END IF
!
!     Investigate whether more components of W can be fixed.
!
120 TEMP=ADELT*ADELT-WFIXSQ
IF (TEMP .GT. ZERO) THEN
WSQSAV=WFIXSQ
STEP=DSQRT(TEMP/GGFREE)
GGFREE=ZERO
DO I=1,N
IF (W(I) .EQ. BIGSTP) THEN
TEMP=XOPT(I)-STEP*GLAG(I)
IF (TEMP .LE. SL(I)) THEN
W(I)=SL(I)-XOPT(I)
WFIXSQ=WFIXSQ+W(I)**2
ELSE IF (TEMP .GE. SU(I)) THEN
W(I)=SU(I)-XOPT(I)
WFIXSQ=WFIXSQ+W(I)**2
ELSE
GGFREE=GGFREE+GLAG(I)**2
END IF
END IF
END DO
IF (WFIXSQ .GT. WSQSAV .AND. GGFREE .GT. ZERO) GOTO 120
END IF
!
!     Set the remaining free components of W and all components of XALT,
!     except that W may be scaled later.
!
GW=ZERO
DO I=1,N
IF (W(I) .EQ. BIGSTP) THEN
W(I)=-STEP*GLAG(I)
XALT(I)=DMAX1(SL(I),DMIN1(SU(I),XOPT(I)+W(I)))
ELSE IF (W(I) .EQ. ZERO) THEN
XALT(I)=XOPT(I)
ELSE IF (GLAG(I) .GT. ZERO) THEN
XALT(I)=SL(I)
ELSE
XALT(I)=SU(I)
END IF
GW=GW+GLAG(I)*W(I)
END DO
!
!     Set CURV to the curvature of the KNEW-th Lagrange function along W.
!     Scale W by a factor less than one if that can reduce the modulus of
!     the Lagrange function at XOPT+W. Set CAUCHY to the final value of
!     the square of this function.
!
CURV=ZERO
DO K=1,NPT
TEMP=ZERO
DO J=1,N
TEMP=TEMP+XPT(K,J)*W(J)
END DO
CURV=CURV+HCOL(K)*TEMP*TEMP
END DO
IF (IFLAG .EQ. 1) CURV=-CURV
IF (CURV .GT. -GW .AND. CURV .LT. -CONST*GW) THEN
SCALE=-GW/CURV
DO I=1,N
TEMP=XOPT(I)+SCALE*W(I)
XALT(I)=DMAX1(SL(I),DMIN1(SU(I),TEMP))
END DO
CAUCHY=(HALF*GW*SCALE)**2
ELSE
CAUCHY=(GW+HALF*CURV)**2
END IF
!
!     If IFLAG is zero, then XALT is calculated as before after reversing
!     the sign of GLAG. Thus two XALT vectors become available. The one that
!     is chosen is the one that gives the larger value of CAUCHY.
!
IF (IFLAG .EQ. 0) THEN
DO I=1,N
GLAG(I)=-GLAG(I)
W(N+I)=XALT(I)
END DO
CSAVE=CAUCHY
IFLAG=1
GOTO 100
END IF
IF (CSAVE .GT. CAUCHY) THEN
DO I=1,N
XALT(I)=W(N+I)
END DO
CAUCHY=CSAVE
END IF
200 RETURN
END


SUBROUTINE BIGDEN (N,NPT,XOPT,XPT,BMAT,ZMAT,IDZ,NDIM,KOPT, KNEW,D,W,VLAG,BETA,S,WVEC,PROD)
implicit integer (i-n), real(real64) (a-h,o-z)
DIMENSION XOPT(*),XPT(NPT,*),BMAT(NDIM,*),ZMAT(NPT,*),D(*), W(*),VLAG(*),S(*),WVEC(NDIM,*),PROD(NDIM,*)
DIMENSION DEN(9),DENEX(9),PAR(9)
!
!     N is the number of variables.
!     NPT is the number of interpolation equations.
!     XOPT is the best interpolation point so far.
!     XPT contains the coordinates of the current interpolation points.
!     BMAT provides the last N columns of H.
!     ZMAT and IDZ give a factorization of the first NPT by NPT submatrix of H.
!     NDIM is the first dimension of BMAT and has the value NPT+N.
!     KOPT is the index of the optimal interpolation point.
!     KNEW is the index of the interpolation point that is going to be moved.
!     D will be set to the step from XOPT to the new point, and on entry it
!       should be the D that was calculated by the last call of BIGLAG. The
!       length of the initial D provides a trust region bound on the final D.
!     W will be set to Wcheck for the final choice of D.
!     VLAG will be set to Theta*Wcheck+e_b for the final choice of D.
!     BETA will be set to the value that will occur in the updating formula
!       when the KNEW-th interpolation point is moved to its new position.
!     S, WVEC, PROD and the private arrays DEN, DENEX and PAR will be used
!       for working space.
!
!     D is calculated in a way that should provide a denominator with a large
!     modulus in the updating formula when the KNEW-th interpolation point is
!     shifted to the new position XOPT+D.
!
!     Set some constants.
!
HALF=0.5D0
ONE=1.0D0
QUART=0.25D0
TWO=2.0D0
ZERO=0.0D0
TWOPI=8.0D0*DATAN(ONE)
NPTM=NPT-N-1
!
!     Store the first NPT elements of the KNEW-th column of H in W(N+1)
!     to W(N+NPT).
!
DO K=1,NPT
W(N+K)=ZERO
END DO
DO J=1,NPTM
TEMP=ZMAT(KNEW,J)
IF (J .LT. IDZ) TEMP=-TEMP
DO K=1,NPT
20 W(N+K)=W(N+K)+TEMP*ZMAT(K,J)
END DO
END DO
ALPHA=W(N+KNEW)
!
!     The initial search direction D is taken from the last call of BIGLAG,
!     and the initial S is set below, usually to the direction from X_OPT
!     to X_KNEW, but a different direction to an interpolation point may
!     be chosen, in order to prevent S from being nearly parallel to D.
!
DD=ZERO
DS=ZERO
SS=ZERO
XOPTSQ=ZERO
DO I=1,N
DD=DD+D(I)**2
S(I)=XPT(KNEW,I)-XOPT(I)
DS=DS+D(I)*S(I)
SS=SS+S(I)**2
XOPTSQ=XOPTSQ+XOPT(I)**2
END DO
IF (DS*DS .GT. 0.99D0*DD*SS) THEN
KSAV=KNEW
DTEST=DS*DS/SS
DO K=1,NPT
IF (K .NE. KOPT) THEN
DSTEMP=ZERO
SSTEMP=ZERO
DO I=1,N
DIFF=XPT(K,I)-XOPT(I)
DSTEMP=DSTEMP+D(I)*DIFF
SSTEMP=SSTEMP+DIFF*DIFF
END DO
IF (DSTEMP*DSTEMP/SSTEMP .LT. DTEST) THEN
KSAV=K
DTEST=DSTEMP*DSTEMP/SSTEMP
DS=DSTEMP
SS=SSTEMP
END IF
END IF
END DO
DO I=1,N
S(I)=XPT(KSAV,I)-XOPT(I)
END DO
END IF
SSDEN=DD*SS-DS*DS
ITERC=0
DENSAV=ZERO
!
!     Begin the iteration by overwriting S with a vector that has the
!     required length and direction.
!
70 ITERC=ITERC+1
TEMP=ONE/DSQRT(SSDEN)
XOPTD=ZERO
XOPTS=ZERO
DO I=1,N
S(I)=TEMP*(DD*S(I)-DS*D(I))
XOPTD=XOPTD+XOPT(I)*D(I)
XOPTS=XOPTS+XOPT(I)*S(I)
END DO
!
!     Set the coefficients of the first two terms of BETA.
!
TEMPA=HALF*XOPTD*XOPTD
TEMPB=HALF*XOPTS*XOPTS
DEN(1)=DD*(XOPTSQ+HALF*DD)+TEMPA+TEMPB
DEN(2)=TWO*XOPTD*DD
DEN(3)=TWO*XOPTS*DD
DEN(4)=TEMPA-TEMPB
DEN(5)=XOPTD*XOPTS
DO I=6,9
DEN(I)=ZERO
END DO
!
!     Put the coefficients of Wcheck in WVEC.
!
DO K=1,NPT
TEMPA=ZERO
TEMPB=ZERO
TEMPC=ZERO
DO I=1,N
TEMPA=TEMPA+XPT(K,I)*D(I)
TEMPB=TEMPB+XPT(K,I)*S(I)
TEMPC=TEMPC+XPT(K,I)*XOPT(I)
END DO
WVEC(K,1)=QUART*(TEMPA*TEMPA+TEMPB*TEMPB)
WVEC(K,2)=TEMPA*TEMPC
WVEC(K,3)=TEMPB*TEMPC
WVEC(K,4)=QUART*(TEMPA*TEMPA-TEMPB*TEMPB)
WVEC(K,5)=HALF*TEMPA*TEMPB
END DO
DO I=1,N
IP=I+NPT
WVEC(IP,1)=ZERO
WVEC(IP,2)=D(I)
WVEC(IP,3)=S(I)
WVEC(IP,4)=ZERO
WVEC(IP,5)=ZERO
END DO
!
!     Put the coefficents of THETA*Wcheck in PROD.
!
DO JC=1,5
NW=NPT
IF (JC .EQ. 2 .OR. JC .EQ. 3) NW=NDIM
DO K=1,NPT
PROD(K,JC)=ZERO
END DO
DO J=1,NPTM
SUM=ZERO
DO K=1,NPT
SUM=SUM+ZMAT(K,J)*WVEC(K,JC)
END DO
IF (J .LT. IDZ) SUM=-SUM
DO K=1,NPT
PROD(K,JC)=PROD(K,JC)+SUM*ZMAT(K,J)
END DO
END DO
IF (NW .EQ. NDIM) THEN
DO K=1,NPT
SUM=ZERO
DO J=1,N
SUM=SUM+BMAT(K,J)*WVEC(NPT+J,JC)
END DO
PROD(K,JC)=PROD(K,JC)+SUM
END DO
END IF
DO J=1,N
SUM=ZERO
DO I=1,NW
SUM=SUM+BMAT(I,J)*WVEC(I,JC)
END DO
PROD(NPT+J,JC)=SUM
END DO
END DO
!
!     Include in DEN the part of BETA that depends on THETA.
!
DO K=1,NDIM
SUM=ZERO
DO I=1,5
PAR(I)=HALF*PROD(K,I)*WVEC(K,I)
SUM=SUM+PAR(I)
END DO
DEN(1)=DEN(1)-PAR(1)-SUM
TEMPA=PROD(K,1)*WVEC(K,2)+PROD(K,2)*WVEC(K,1)
TEMPB=PROD(K,2)*WVEC(K,4)+PROD(K,4)*WVEC(K,2)
TEMPC=PROD(K,3)*WVEC(K,5)+PROD(K,5)*WVEC(K,3)
DEN(2)=DEN(2)-TEMPA-HALF*(TEMPB+TEMPC)
DEN(6)=DEN(6)-HALF*(TEMPB-TEMPC)
TEMPA=PROD(K,1)*WVEC(K,3)+PROD(K,3)*WVEC(K,1)
TEMPB=PROD(K,2)*WVEC(K,5)+PROD(K,5)*WVEC(K,2)
TEMPC=PROD(K,3)*WVEC(K,4)+PROD(K,4)*WVEC(K,3)
DEN(3)=DEN(3)-TEMPA-HALF*(TEMPB-TEMPC)
DEN(7)=DEN(7)-HALF*(TEMPB+TEMPC)
TEMPA=PROD(K,1)*WVEC(K,4)+PROD(K,4)*WVEC(K,1)
DEN(4)=DEN(4)-TEMPA-PAR(2)+PAR(3)
TEMPA=PROD(K,1)*WVEC(K,5)+PROD(K,5)*WVEC(K,1)
TEMPB=PROD(K,2)*WVEC(K,3)+PROD(K,3)*WVEC(K,2)
DEN(5)=DEN(5)-TEMPA-HALF*TEMPB
DEN(8)=DEN(8)-PAR(4)+PAR(5)
TEMPA=PROD(K,4)*WVEC(K,5)+PROD(K,5)*WVEC(K,4)
DEN(9)=DEN(9)-HALF*TEMPA
END DO
!
!     Extend DEN so that it holds all the coefficients of DENOM.
!
SUM=ZERO
DO I=1,5
PAR(I)=HALF*PROD(KNEW,I)**2
SUM=SUM+PAR(I)
END DO
DENEX(1)=ALPHA*DEN(1)+PAR(1)+SUM
TEMPA=TWO*PROD(KNEW,1)*PROD(KNEW,2)
TEMPB=PROD(KNEW,2)*PROD(KNEW,4)
TEMPC=PROD(KNEW,3)*PROD(KNEW,5)
DENEX(2)=ALPHA*DEN(2)+TEMPA+TEMPB+TEMPC
DENEX(6)=ALPHA*DEN(6)+TEMPB-TEMPC
TEMPA=TWO*PROD(KNEW,1)*PROD(KNEW,3)
TEMPB=PROD(KNEW,2)*PROD(KNEW,5)
TEMPC=PROD(KNEW,3)*PROD(KNEW,4)
DENEX(3)=ALPHA*DEN(3)+TEMPA+TEMPB-TEMPC
DENEX(7)=ALPHA*DEN(7)+TEMPB+TEMPC
TEMPA=TWO*PROD(KNEW,1)*PROD(KNEW,4)
DENEX(4)=ALPHA*DEN(4)+TEMPA+PAR(2)-PAR(3)
TEMPA=TWO*PROD(KNEW,1)*PROD(KNEW,5)
DENEX(5)=ALPHA*DEN(5)+TEMPA+PROD(KNEW,2)*PROD(KNEW,3)
DENEX(8)=ALPHA*DEN(8)+PAR(4)-PAR(5)
DENEX(9)=ALPHA*DEN(9)+PROD(KNEW,4)*PROD(KNEW,5)
!
!     Seek the value of the angle that maximizes the modulus of DENOM.
!
SUM=DENEX(1)+DENEX(2)+DENEX(4)+DENEX(6)+DENEX(8)
DENOLD=SUM
DENMAX=SUM
ISAVE=0
IU=49
TEMP=TWOPI/DBLE(IU+1)
PAR(1)=ONE
DO I=1,IU
ANGLE=DBLE(I)*TEMP
PAR(2)=DCOS(ANGLE)
PAR(3)=DSIN(ANGLE)
DO J=4,8,2
PAR(J)=PAR(2)*PAR(J-2)-PAR(3)*PAR(J-1)
PAR(J+1)=PAR(2)*PAR(J-1)+PAR(3)*PAR(J-2)
END DO
SUMOLD=SUM
SUM=ZERO
DO J=1,9
SUM=SUM+DENEX(J)*PAR(J)
END DO
IF (DABS(SUM) .GT. DABS(DENMAX)) THEN
DENMAX=SUM
ISAVE=I
TEMPA=SUMOLD
ELSE IF (I .EQ. ISAVE+1) THEN
TEMPB=SUM
END IF
END DO
IF (ISAVE .EQ. 0) TEMPA=SUM
IF (ISAVE .EQ. IU) TEMPB=DENOLD
STEP=ZERO
IF (TEMPA .NE. TEMPB) THEN
TEMPA=TEMPA-DENMAX
TEMPB=TEMPB-DENMAX
STEP=HALF*(TEMPA-TEMPB)/(TEMPA+TEMPB)
END IF
ANGLE=TEMP*(DBLE(ISAVE)+STEP)
!
!     Calculate the new parameters of the denominator, the new VLAG vector
!     and the new D. Then test for convergence.
!
PAR(2)=DCOS(ANGLE)
PAR(3)=DSIN(ANGLE)
DO J=4,8,2
PAR(J)=PAR(2)*PAR(J-2)-PAR(3)*PAR(J-1)
PAR(J+1)=PAR(2)*PAR(J-1)+PAR(3)*PAR(J-2)
END DO
BETA=ZERO
DENMAX=ZERO
DO J=1,9
BETA=BETA+DEN(J)*PAR(J)
DENMAX=DENMAX+DENEX(J)*PAR(J)
END DO
DO K=1,NDIM
VLAG(K)=ZERO
DO J=1,5
VLAG(K)=VLAG(K)+PROD(K,J)*PAR(J)
END DO
END DO
TAU=VLAG(KNEW)
DD=ZERO
TEMPA=ZERO
TEMPB=ZERO
DO I=1,N
D(I)=PAR(2)*D(I)+PAR(3)*S(I)
W(I)=XOPT(I)+D(I)
DD=DD+D(I)**2
TEMPA=TEMPA+D(I)*W(I)
TEMPB=TEMPB+W(I)*W(I)
END DO
IF (ITERC .GE. N) GOTO 340
IF (ITERC .GT. 1) DENSAV=DMAX1(DENSAV,DENOLD)
IF (DABS(DENMAX) .LE. 1.1D0*DABS(DENSAV)) GOTO 340
DENSAV=DENMAX
!
!     Set S to half the gradient of the denominator with respect to D.
!     Then branch for the next iteration.
!
DO I=1,N
TEMP=TEMPA*XOPT(I)+TEMPB*D(I)-VLAG(NPT+I)
S(I)=TAU*BMAT(KNEW,I)+ALPHA*TEMP
END DO
DO K=1,NPT
SUM=ZERO
DO J=1,N
SUM=SUM+XPT(K,J)*W(J)
END DO
TEMP=(TAU*W(N+K)-ALPHA*VLAG(K))*SUM
DO I=1,N
S(I)=S(I)+TEMP*XPT(K,I)
END DO
END DO
SS=ZERO
DS=ZERO
DO I=1,N
SS=SS+S(I)**2
DS=DS+D(I)*S(I)
END DO
SSDEN=DD*SS-DS*DS
IF (SSDEN .GE. 1.0D-8*DD*SS) GOTO 70
!
!     Set the vector W before the RETURN from the subroutine.
!
340 DO K=1,NDIM
W(K)=ZERO
DO J=1,5
W(K)=W(K)+WVEC(K,J)*PAR(J)
END DO
END DO
VLAG(KOPT)=VLAG(KOPT)+ONE
RETURN
END


SUBROUTINE BIGLAG (N,NPT,XOPT,XPT,BMAT,ZMAT,IDZ,NDIM,KNEW, DELTA,D,ALPHA,HCOL,GC,GD,S,W)
implicit integer (i-n), real(real64) (a-h,o-z)
DIMENSION XOPT(*),XPT(NPT,*),BMAT(NDIM,*),ZMAT(NPT,*),D(*), HCOL(*),GC(*),GD(*),S(*),W(*)
!
!     N is the number of variables.
!     NPT is the number of interpolation equations.
!     XOPT is the best interpolation point so far.
!     XPT contains the coordinates of the current interpolation points.
!     BMAT provides the last N columns of H.
!     ZMAT and IDZ give a factorization of the first NPT by NPT submatrix of H.
!     NDIM is the first dimension of BMAT and has the value NPT+N.
!     KNEW is the index of the interpolation point that is going to be moved.
!     DELTA is the current trust region bound.
!     D will be set to the step from XOPT to the new point.
!     ALPHA will be set to the KNEW-th diagonal element of the H matrix.
!     HCOL, GC, GD, S and W will be used for working space.
!
!     The step D is calculated in a way that attempts to maximize the modulus
!     of LFUNC(XOPT+D), subject to the bound ||D|| .LE. DELTA, where LFUNC is
!     the KNEW-th Lagrange function.
!
!     Set some constants.
!
HALF=0.5D0
ONE=1.0D0
ZERO=0.0D0
TWOPI=8.0D0*DATAN(ONE)
DELSQ=DELTA*DELTA
NPTM=NPT-N-1
!
!     Set the first NPT components of HCOL to the leading elements of the
!     KNEW-th column of H.
!
ITERC=0
DO K=1,NPT
HCOL(K)=ZERO
END DO
DO J=1,NPTM
TEMP=ZMAT(KNEW,J)
IF (J .LT. IDZ) TEMP=-TEMP
DO K=1,NPT
HCOL(K)=HCOL(K)+TEMP*ZMAT(K,J)
END DO
END DO
ALPHA=HCOL(KNEW)
!
!     Set the unscaled initial direction D. Form the gradient of LFUNC at
!     XOPT, and multiply D by the second derivative matrix of LFUNC.
!
DD=ZERO
DO I=1,N
D(I)=XPT(KNEW,I)-XOPT(I)
GC(I)=BMAT(KNEW,I)
GD(I)=ZERO
DD=DD+D(I)**2
END DO
DO K=1,NPT
TEMP=ZERO
SUM=ZERO
DO J=1,N
TEMP=TEMP+XPT(K,J)*XOPT(J)
SUM=SUM+XPT(K,J)*D(J)
END DO
TEMP=HCOL(K)*TEMP
SUM=HCOL(K)*SUM
DO I=1,N
GC(I)=GC(I)+TEMP*XPT(K,I)
GD(I)=GD(I)+SUM*XPT(K,I)
END DO
END DO
!
!     Scale D and GD, with a sign change if required. Set S to another
!     vector in the initial two dimensional subspace.
!
GG=ZERO
SP=ZERO
DHD=ZERO
DO I=1,N
GG=GG+GC(I)**2
SP=SP+D(I)*GC(I)
DHD=DHD+D(I)*GD(I)
END DO
SCALE=DELTA/DSQRT(DD)
IF (SP*DHD .LT. ZERO) SCALE=-SCALE
TEMP=ZERO
IF (SP*SP .GT. 0.99D0*DD*GG) TEMP=ONE
TAU=SCALE*(DABS(SP)+HALF*SCALE*DABS(DHD))
IF (GG*DELSQ .LT. 0.01D0*TAU*TAU) TEMP=ONE
DO I=1,N
D(I)=SCALE*D(I)
GD(I)=SCALE*GD(I)
S(I)=GC(I)+TEMP*GD(I)
END DO
!
!     Begin the iteration by overwriting S with a vector that has the
!     required length and direction, except that termination occurs if
!     the given D and S are nearly parallel.
!
80 ITERC=ITERC+1
DD=ZERO
SP=ZERO
SS=ZERO
DO I=1,N
DD=DD+D(I)**2
SP=SP+D(I)*S(I)
SS=SS+S(I)**2
END DO
TEMP=DD*SS-SP*SP
IF (TEMP .LE. 1.0D-8*DD*SS) GOTO 160
DENOM=DSQRT(TEMP)
DO I=1,N
S(I)=(DD*S(I)-SP*D(I))/DENOM
W(I)=ZERO
END DO
!
!     Calculate the coefficients of the objective function on the circle,
!     beginning with the multiplication of S by the second derivative matrix.
!
DO K=1,NPT
SUM=ZERO
DO J=1,N
SUM=SUM+XPT(K,J)*S(J)
END DO
SUM=HCOL(K)*SUM
DO I=1,N
W(I)=W(I)+SUM*XPT(K,I)
END DO
END DO
CF1=ZERO
CF2=ZERO
CF3=ZERO
CF4=ZERO
CF5=ZERO
DO I=1,N
CF1=CF1+S(I)*W(I)
CF2=CF2+D(I)*GC(I)
CF3=CF3+S(I)*GC(I)
CF4=CF4+D(I)*GD(I)
CF5=CF5+S(I)*GD(I)
END DO
CF1=HALF*CF1
CF4=HALF*CF4-CF1
!
!     Seek the value of the angle that maximizes the modulus of TAU.
!
TAUBEG=CF1+CF2+CF4
TAUMAX=TAUBEG
TAUOLD=TAUBEG
ISAVE=0
IU=49
TEMP=TWOPI/DBLE(IU+1)
DO I=1,IU
ANGLE=DBLE(I)*TEMP
CTH=DCOS(ANGLE)
STH=DSIN(ANGLE)
TAU=CF1+(CF2+CF4*CTH)*CTH+(CF3+CF5*CTH)*STH
IF (DABS(TAU) .GT. DABS(TAUMAX)) THEN
TAUMAX=TAU
ISAVE=I
TEMPA=TAUOLD
ELSE IF (I .EQ. ISAVE+1) THEN
TEMPB=TAU
END IF
TAUOLD=TAU
END DO
IF (ISAVE .EQ. 0) TEMPA=TAU
IF (ISAVE .EQ. IU) TEMPB=TAUBEG
STEP=ZERO
IF (TEMPA .NE. TEMPB) THEN
TEMPA=TEMPA-TAUMAX
TEMPB=TEMPB-TAUMAX
STEP=HALF*(TEMPA-TEMPB)/(TEMPA+TEMPB)
END IF
ANGLE=TEMP*(DBLE(ISAVE)+STEP)
!
!     Calculate the new D and GD. Then test for convergence.
!
CTH=DCOS(ANGLE)
STH=DSIN(ANGLE)
TAU=CF1+(CF2+CF4*CTH)*CTH+(CF3+CF5*CTH)*STH
DO I=1,N
D(I)=CTH*D(I)+STH*S(I)
GD(I)=CTH*GD(I)+STH*W(I)
S(I)=GC(I)+GD(I)
END DO
IF (DABS(TAU) .LE. 1.1D0*DABS(TAUBEG)) GOTO 160
IF (ITERC .LT. N) GOTO 80
160 RETURN
END


subroutine bobyqa_core (N,NPT,X,XL,XU,RHOBEG,RHOEND,IPRINT, MAXFUN,W,IERR)
implicit integer (i-n), real(real64) (a-h,o-z)
DIMENSION X(*),XL(*),XU(*),W(*)
!
!     This subroutine seeks the least value of a function of many variables,
!     by applying a trust region method that forms quadratic models by
!     interpolation. There is usually some freedom in the interpolation
!     conditions, which is taken up by minimizing the Frobenius norm of
!     the change to the second derivative of the model, beginning with the
!     zero matrix. The values of the variables are constrained by upper and
!     lower bounds. The arguments of the subroutine are as follows.
!
!     N must be set to the number of variables and must be at least two.
!     NPT is the number of interpolation conditions. Its value must be in
!       the interval [N+2,(N+1)(N+2)/2]. Choices that exceed 2*N+1 are not
!       recommended.
!     Initial values of the variables must be set in X(1),X(2),...,X(N). They
!       will be changed to the values that give the least calculated F.
!     For I=1,2,...,N, XL(I) and XU(I) must provide the lower and upper
!       bounds, respectively, on X(I). The construction of quadratic models
!       requires XL(I) to be strictly less than XU(I) for each I. Further,
!       the contribution to a model from changes to the I-th variable is
!       damaged severely by rounding errors if XU(I)-XL(I) is too small.
!     RHOBEG and RHOEND must be set to the initial and final values of a trust
!       region radius, so both must be positive with RHOEND no greater than
!       RHOBEG. Typically, RHOBEG should be about one tenth of the greatest
!       expected change to a variable, while RHOEND should indicate the
!       accuracy that is required in the final values of the variables. An
!       error return occurs if any of the differences XU(I)-XL(I), I=1,...,N,
!       is less than 2*RHOBEG.
!     The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
!       amount of printing. Specifically, there is no output if IPRINT=0 and
!       there is output only at the return if IPRINT=1. Otherwise, each new
!       value of RHO is printed, with the best vector of variables so far and
!       the corresponding value of the objective function. Further, each new
!       value of F with its variables are output if IPRINT=3.
!     MAXFUN must be set to an upper bound on the number of calls of CALFUN.
!     The array W will be used for working space. Its length must be at least
!       (NPT+5)*(NPT+N)+3*N*(N+5)/2.
!JN   Add IERR to tell what the error is to minqer.
!     IERR gives an error code that is interpreted by minqer interface routine.
!       and passed back to minqa.R
!
!     DOUBLE PRECISION FUNCTION CALFUN (N,X,IP) has to be provided by
!     the user. It returns the value of the objective function for
!     the current values of the variables X(1),X(2),...,X(N), which are
!     generated automatically in a way that satisfies the bounds given
!     in XL and XU.
!
!     Return if the value of NPT is unacceptable.
!
!  Modified by John Nash to put the setup controls into the R code.
NP=N+1
!J  Comment out to END IF as in R code
IF (NPT .LT. N+2 .OR. NPT .GT. ((N+2)*NP)/2) THEN
!JN         CALL minqer(10)
!$$$          PRINT 10
!$$$   10     FORMAT (/4X,'Return from BOBYQA because NPT is not in',
!$$$     1      ' the required interval')
!$$$          GO TO 40
IERR = 10
!JN    Fail out NPT has unacceptable value
GO TO 40
END IF
!
!     Partition the working space array, so that different parts of it can
!     be treated separately during the calculation of BOBYQB. The partition
!     requires the first (NPT+2)*(NPT+N)+3*N*(N+5)/2 elements of W plus the
!     space that is taken by the last array in the argument list of BOBYQB.
!
NDIM=NPT+N
IXB=1
IXP=IXB+N
IFV=IXP+N*NPT
IXO=IFV+NPT
IGO=IXO+N
IHQ=IGO+N
IPQ=IHQ+(N*NP)/2
IBMAT=IPQ+NPT
IZMAT=IBMAT+NDIM*N
ISL=IZMAT+NPT*(NPT-NP)
ISU=ISL+N
IXN=ISU+N
IXA=IXN+N
ID=IXA+N
IVL=ID+N
IW=IVL+NDIM
!JN 100807 to ensure a defined value for the return code
IERR = 0
!
!     Return if there is insufficient space between the bounds. Modify the
!     initial X if necessary in order to avoid conflicts between the bounds
!     and the construction of the first quadratic model. The lower and upper
!     bounds on moves from the updated X are set now, in the ISL and ISU
!     partitions of W, in order to provide useful and exact information about
!     components of X that become within distance RHOBEG from their bounds.
!
ZERO=0.0D0
DO J=1,N
TEMP=XU(J)-XL(J)
IF (TEMP .LT. RHOBEG+RHOBEG) THEN
!     JN         CALL minqer(20)
!$$$  PRINT 20
!$$$  20     FORMAT (/4X,'Return from BOBYQA because one of the',
!$$$  1      ' differences XU(I)-XL(I)'/6X,' is less than 2*RHOBEG.')
!$$$  GO TO 40
IERR = 20
GOTO 40
END IF
JSL=ISL+J-1
JSU=JSL+N
W(JSL)=XL(J)-X(J)
W(JSU)=XU(J)-X(J)
IF (W(JSL) .GE. -RHOBEG) THEN
IF (W(JSL) .GE. ZERO) THEN
X(J)=XL(J)
W(JSL)=ZERO
W(JSU)=TEMP
ELSE
X(J)=XL(J)+RHOBEG
W(JSL)=-RHOBEG
W(JSU)=DMAX1(XU(J)-X(J),RHOBEG)
END IF
ELSE IF (W(JSU) .LE. RHOBEG) THEN
IF (W(JSU) .LE. ZERO) THEN
X(J)=XU(J)
W(JSL)=-TEMP
W(JSU)=ZERO
ELSE
X(J)=XU(J)-RHOBEG
W(JSL)=DMIN1(XL(J)-X(J),-RHOBEG)
W(JSU)=RHOBEG
END IF
END IF
END DO
!
!     Make the call of BOBYQB.
!
CALL BOBYQB (N,NPT,X,XL,XU,RHOBEG,RHOEND,IPRINT,MAXFUN,W(IXB), W(IXP),W(IFV),W(IXO),W(IGO),W(IHQ),W(IPQ), &
    & W(IBMAT),W(IZMAT), NDIM,W(ISL),W(ISU),W(IXN),W(IXA),W(ID),W(IVL),W(IW), IERR)
40 RETURN
!JN Added 100807 -- put label 40
END


SUBROUTINE BOBYQB (N,NPT,X,XL,XU,RHOBEG,RHOEND,IPRINT, MAXFUN,XBASE,XPT,FVAL,XOPT,GOPT,HQ,PQ,BMAT,ZMAT, &
    & NDIM, SL,SU,XNEW,XALT,D,VLAG,W,IERR)
implicit integer (i-n), real(real64) (a-h,o-z)
DIMENSION X(*),XL(*),XU(*),XBASE(*),XPT(NPT,*),FVAL(*), XOPT(*),GOPT(*),HQ(*),PQ(*),BMAT(NDIM,*), &
    & ZMAT(NPT,*), SL(*),SU(*),XNEW(*),XALT(*),D(*),VLAG(*),W(*)
!
!     The arguments N, NPT, X, XL, XU, RHOBEG, RHOEND, IPRINT and MAXFUN
!       are identical to the corresponding arguments in SUBROUTINE BOBYQA.
!     XBASE holds a shift of origin that should reduce the contributions
!       from rounding errors to values of the model and Lagrange functions.
!     XPT is a two-dimensional array that holds the coordinates of the
!       interpolation points relative to XBASE.
!     FVAL holds the values of F at the interpolation points.
!     XOPT is set to the displacement from XBASE of the trust region centre.
!     GOPT holds the gradient of the quadratic model at XBASE+XOPT.
!     HQ holds the explicit second derivatives of the quadratic model.
!     PQ contains the parameters of the implicit second derivatives of the
!       quadratic model.
!     BMAT holds the last N columns of H.
!     ZMAT holds the factorization of the leading NPT by NPT submatrix of H,
!       this factorization being ZMAT times ZMAT^T, which provides both the
!       correct rank and positive semi-definiteness.
!     NDIM is the first dimension of BMAT and has the value NPT+N.
!     SL and SU hold the differences XL-XBASE and XU-XBASE, respectively.
!       All the components of every XOPT are going to satisfy the bounds
!       SL(I) .LEQ. XOPT(I) .LEQ. SU(I), with appropriate equalities when
!       XOPT is on a constraint boundary.
!     XNEW is chosen by SUBROUTINE TRSBOX or ALTMOV. Usually XBASE+XNEW is the
!       vector of variables for the next call of CALFUN. XNEW also satisfies
!       the SL and SU constraints in the way that has just been mentioned.
!     XALT is an alternative to XNEW, chosen by ALTMOV, that may replace XNEW
!       in order to increase the denominator in the updating of UPDATE.
!     D is reserved for a trial step from XOPT, which is usually XNEW-XOPT.
!     VLAG contains the values of the Lagrange functions at a new point X.
!       They are part of a product that requires VLAG to be of length NDIM.
!     W is a one-dimensional array that is used for working space. Its length
!       must be at least 3*NDIM = 3*(NPT+N).
!JN 100807
!     IERR is an error code to tell calling program WHICH error occurred.
!JN       Note that it is defined in BOBYQA to 0 initially.
!
!     Set some constants.
!
HALF=0.5D0
ONE=1.0D0
TEN=10.0D0
TENTH=0.1D0
TWO=2.0D0
ZERO=0.0D0
NP=N+1
NPTM=NPT-NP
NH=(N*NP)/2
!
!     The call of PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ,
!     BMAT and ZMAT for the first iteration, with the corresponding values of
!     of NF and KOPT, which are the number of calls of CALFUN so far and the
!     index of the interpolation point at the trust region centre. Then the
!     initial XOPT is set too. The branch to label 720 occurs if MAXFUN is
!     less than NPT. GOPT will be updated if KOPT is different from KBASE.
!
CALL PRELIM (N,NPT,X,XL,XU,RHOBEG,IPRINT,MAXFUN,XBASE,XPT, FVAL,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,SL,SU,NF,KOPT)
XOPTSQ=ZERO
DO I=1,N
XOPT(I)=XPT(KOPT,I)
XOPTSQ=XOPTSQ+XOPT(I)**2
END DO
FSAVE=FVAL(1)
IF (NF .LT. NPT) THEN
!JN 100807
IERR=390
GOTO 720
!     JN         CALL minqer(390)
!$$$          IF (IPRINT .GT. 0) PRINT 390
!$$$          GOTO 720
END IF
KBASE=1
!
!     Complete the settings that are required for the iterative procedure.
!
RHO=RHOBEG
DELTA=RHO
NRESC=NF
NTRITS=0
DIFFA=ZERO
DIFFB=ZERO
ITEST=0
NFSAV=NF
!
!     Update GOPT if necessary before the first iteration and after each
!     call of RESCUE that makes a call of CALFUN.
!
20 IF (KOPT .NE. KBASE) THEN
IH=0
DO J=1,N
DO I=1,J
IH=IH+1
IF (I .LT. J) GOPT(J)=GOPT(J)+HQ(IH)*XOPT(I)
GOPT(I)=GOPT(I)+HQ(IH)*XOPT(J)
END DO
END DO
IF (NF .GT. NPT) THEN
DO K=1,NPT
TEMP=ZERO
DO J=1,N
TEMP=TEMP+XPT(K,J)*XOPT(J)
END DO
TEMP=PQ(K)*TEMP
DO I=1,N
GOPT(I)=GOPT(I)+TEMP*XPT(K,I)
END DO
END DO
END IF
END IF
!
!     Generate the next point in the trust region that provides a small value
!     of the quadratic model subject to the constraints on the variables.
!     The integer NTRITS is set to the number "trust region" iterations that
!     have occurred since the last "alternative" iteration. If the length
!     of XNEW-XOPT is less than HALF*RHO, however, then there is a branch to
!     label 650 or 680 with NTRITS=-1, instead of calculating F at XNEW.
!
60 CALL TRSBOX (N,NPT,XPT,XOPT,GOPT,HQ,PQ,SL,SU,DELTA,XNEW,D, W,W(NP),W(NP+N),W(NP+2*N),W(NP+3*N),DSQ, &
    & CRVMIN)
DNORM=DMIN1(DELTA,DSQRT(DSQ))
IF (DNORM .LT. HALF*RHO) THEN
NTRITS=-1
DISTSQ=(TEN*RHO)**2
IF (NF .LE. NFSAV+2) GOTO 650
!
!     The following choice between labels 650 and 680 depends on whether or
!     not our work with the current RHO seems to be complete. Either RHO is
!     decreased or termination occurs if the errors in the quadratic model at
!     the last three interpolation points compare favourably with predictions
!     of likely improvements to the model within distance HALF*RHO of XOPT.
!
ERRBIG=DMAX1(DIFFA,DIFFB,DIFFC)
FRHOSQ=0.125D0*RHO*RHO
IF (CRVMIN .GT. ZERO .AND. ERRBIG .GT. FRHOSQ*CRVMIN) GOTO 650
BDTOL=ERRBIG/RHO
DO J=1,N
BDTEST=BDTOL
IF (XNEW(J) .EQ. SL(J)) BDTEST=W(J)
IF (XNEW(J) .EQ. SU(J)) BDTEST=-W(J)
IF (BDTEST .LT. BDTOL) THEN
CURV=HQ((J+J*J)/2)
DO K=1,NPT
CURV=CURV+PQ(K)*XPT(K,J)**2
END DO
BDTEST=BDTEST+HALF*CURV*RHO
IF (BDTEST .LT. BDTOL) GOTO 650
END IF
END DO
GOTO 680
END IF
NTRITS=NTRITS+1
!
!     Severe cancellation is likely to occur if XOPT is too far from XBASE.
!     If the following test holds, then XBASE is shifted so that XOPT becomes
!     zero. The appropriate changes are made to BMAT and to the second
!     derivatives of the current model, beginning with the changes to BMAT
!     that do not depend on ZMAT. VLAG is used temporarily for working space.
!
90 IF (DSQ .LE. 1.0D-3*XOPTSQ) THEN
FRACSQ=0.25D0*XOPTSQ
SUMPQ=ZERO
DO K=1,NPT
SUMPQ=SUMPQ+PQ(K)
SUM=-HALF*XOPTSQ
DO I=1,N
SUM=SUM+XPT(K,I)*XOPT(I)
END DO
W(NPT+K)=SUM
TEMP=FRACSQ-HALF*SUM
DO I=1,N
W(I)=BMAT(K,I)
VLAG(I)=SUM*XPT(K,I)+TEMP*XOPT(I)
IP=NPT+I
DO J=1,I
BMAT(IP,J)=BMAT(IP,J)+W(I)*VLAG(J)+VLAG(I)*W(J)
END DO
END DO
END DO
!
!     Then the revisions of BMAT that depend on ZMAT are calculated.
!
DO JJ=1,NPTM
SUMZ=ZERO
SUMW=ZERO
DO K=1,NPT
SUMZ=SUMZ+ZMAT(K,JJ)
VLAG(K)=W(NPT+K)*ZMAT(K,JJ)
SUMW=SUMW+VLAG(K)
END DO
DO J=1,N
SUM=(FRACSQ*SUMZ-HALF*SUMW)*XOPT(J)
DO K=1,NPT
SUM=SUM+VLAG(K)*XPT(K,J)
END DO
W(J)=SUM
DO K=1,NPT
BMAT(K,J)=BMAT(K,J)+SUM*ZMAT(K,JJ)
END DO
END DO
DO I=1,N
IP=I+NPT
TEMP=W(I)
DO J=1,I
BMAT(IP,J)=BMAT(IP,J)+TEMP*W(J)
END DO
END DO
END DO
!
!     The following instructions complete the shift, including the changes
!     to the second derivative parameters of the quadratic model.
!
IH=0
DO J=1,N
W(J)=-HALF*SUMPQ*XOPT(J)
DO K=1,NPT
W(J)=W(J)+PQ(K)*XPT(K,J)
XPT(K,J)=XPT(K,J)-XOPT(J)
END DO
DO I=1,J
IH=IH+1
HQ(IH)=HQ(IH)+W(I)*XOPT(J)+XOPT(I)*W(J)
BMAT(NPT+I,J)=BMAT(NPT+J,I)
END DO
END DO
DO I=1,N
XBASE(I)=XBASE(I)+XOPT(I)
XNEW(I)=XNEW(I)-XOPT(I)
SL(I)=SL(I)-XOPT(I)
SU(I)=SU(I)-XOPT(I)
XOPT(I)=ZERO
END DO
XOPTSQ=ZERO
END IF
IF (NTRITS .EQ. 0) GOTO 210
GOTO 230
!
!     XBASE is also moved to XOPT by a call of RESCUE. This calculation is
!     more expensive than the previous shift, because new matrices BMAT and
!     ZMAT are generated from scratch, which may include the replacement of
!     interpolation points whose positions seem to be causing near linear
!     dependence in the interpolation conditions. Therefore RESCUE is called
!     only if rounding errors have reduced by at least a factor of two the
!     denominator of the formula for updating the H matrix. It provides a
!     useful safeguard, but is not invoked in most applications of BOBYQA.
!
190 NFSAV=NF
KBASE=KOPT
CALL RESCUE (N,NPT,XL,XU,IPRINT,MAXFUN,XBASE,XPT,FVAL, XOPT,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,SL,SU,NF,DELTA, &
    & KOPT, VLAG,W,W(N+NP),W(NDIM+NP))
!
!     XOPT is updated now in case the branch below to label 720 is taken.
!     Any updating of GOPT occurs after the branch below to label 20, which
!     leads to a trust region iteration as does the branch to label 60.
!
XOPTSQ=ZERO
IF (KOPT .NE. KBASE) THEN
DO I=1,N
XOPT(I)=XPT(KOPT,I)
XOPTSQ=XOPTSQ+XOPT(I)**2
END DO
END IF
IF (NF .LT. 0) THEN
NF=MAXFUN
!JN          CALL minqer(390)
!$$$          IF (IPRINT .GT. 0) PRINT 390
!$$$          GOTO 720
!JN 100807
IERR=390
GOTO 720
END IF
NRESC=NF
IF (NFSAV .LT. NF) THEN
NFSAV=NF
GOTO 20
END IF
IF (NTRITS .GT. 0) GOTO 60
!
!     Pick two alternative vectors of variables, relative to XBASE, that
!     are suitable as new positions of the KNEW-th interpolation point.
!     Firstly, XNEW is set to the point on a line through XOPT and another
!     interpolation point that minimizes the predicted value of the next
!     denominator, subject to ||XNEW - XOPT|| .LEQ. ADELT and to the SL
!     and SU bounds. Secondly, XALT is set to the best feasible point on
!     a constrained version of the Cauchy step of the KNEW-th Lagrange
!     function, the corresponding value of the square of this function
!     being returned in CAUCHY. The choice between these alternatives is
!     going to be made when the denominator is calculated.
!
210 CALL ALTMOV (N,NPT,XPT,XOPT,BMAT,ZMAT,NDIM,SL,SU,KOPT, KNEW,ADELT,XNEW,XALT,ALPHA,CAUCHY,W,W(NP), &
    & W(NDIM+1))
DO I=1,N
D(I)=XNEW(I)-XOPT(I)
END DO
!
!     Calculate VLAG and BETA for the current choice of D. The scalar
!     product of D with XPT(K,.) is going to be held in W(NPT+K) for
!     use when VQUAD is calculated.
!
230 DO K=1,NPT
SUMA=ZERO
SUMB=ZERO
SUM=ZERO
DO J=1,N
SUMA=SUMA+XPT(K,J)*D(J)
SUMB=SUMB+XPT(K,J)*XOPT(J)
SUM=SUM+BMAT(K,J)*D(J)
END DO
W(K)=SUMA*(HALF*SUMA+SUMB)
VLAG(K)=SUM
W(NPT+K)=SUMA
END DO
BETA=ZERO
DO JJ=1,NPTM
SUM=ZERO
DO K=1,NPT
SUM=SUM+ZMAT(K,JJ)*W(K)
END DO
BETA=BETA-SUM*SUM
DO K=1,NPT
VLAG(K)=VLAG(K)+SUM*ZMAT(K,JJ)
END DO
END DO
DSQ=ZERO
BSUM=ZERO
DX=ZERO
DO J=1,N
DSQ=DSQ+D(J)**2
SUM=ZERO
DO K=1,NPT
SUM=SUM+W(K)*BMAT(K,J)
END DO
BSUM=BSUM+SUM*D(J)
JP=NPT+J
DO I=1,N
SUM=SUM+BMAT(JP,I)*D(I)
END DO
VLAG(JP)=SUM
BSUM=BSUM+SUM*D(J)
DX=DX+D(J)*XOPT(J)
END DO
BETA=DX*DX+DSQ*(XOPTSQ+DX+DX+HALF*DSQ)+BETA-BSUM
VLAG(KOPT)=VLAG(KOPT)+ONE
!
!     If NTRITS is zero, the denominator may be increased by replacing
!     the step D of ALTMOV by a Cauchy step. Then RESCUE may be called if
!     rounding errors have damaged the chosen denominator.
!
IF (NTRITS .EQ. 0) THEN
DENOM=VLAG(KNEW)**2+ALPHA*BETA
IF (DENOM .LT. CAUCHY .AND. CAUCHY .GT. ZERO) THEN
DO I=1,N
XNEW(I)=XALT(I)
D(I)=XNEW(I)-XOPT(I)
END DO
CAUCHY=ZERO
GO TO 230
END IF
IF (DENOM .LE. HALF*VLAG(KNEW)**2) THEN
IF (NF .GT. NRESC) GOTO 190
!     JN              IF (IPRINT .GT. 0) CALL minqer(320)
!$$$              PRINT 320
!$$$  320         FORMAT (/5X,'Return from BOBYQA because of much',
!$$$     1          ' cancellation in a denominator.')
!$$$              GOTO 720
!JN 100807
IERR=320
GOTO 720
END IF
!
!     Alternatively, if NTRITS is positive, then set KNEW to the index of
!     the next interpolation point to be deleted to make room for a trust
!     region step. Again RESCUE may be called if rounding errors have damaged
!     the chosen denominator, which is the reason for attempting to select
!     KNEW before calculating the next value of the objective function.
!
ELSE
DELSQ=DELTA*DELTA
SCADEN=ZERO
BIGLSQ=ZERO
KNEW=0
DO K=1,NPT
IF (K .EQ. KOPT) GOTO 350
HDIAG=ZERO
DO JJ=1,NPTM
HDIAG=HDIAG+ZMAT(K,JJ)**2
END DO
DEN=BETA*HDIAG+VLAG(K)**2
DISTSQ=ZERO
DO J=1,N
DISTSQ=DISTSQ+(XPT(K,J)-XOPT(J))**2
END DO
TEMP=DMAX1(ONE,(DISTSQ/DELSQ)**2)
IF (TEMP*DEN .GT. SCADEN) THEN
SCADEN=TEMP*DEN
KNEW=K
DENOM=DEN
END IF
BIGLSQ=DMAX1(BIGLSQ,TEMP*VLAG(K)**2)
350 END DO
IF (SCADEN .LE. HALF*BIGLSQ) THEN
IF (NF .GT. NRESC) GOTO 190
!JN              IF (IPRINT .GT. 0) CALL minqer(320)
!$$$              PRINT 320
!$$$              GOTO 720
!JN 100807
IERR=320
GOTO 720
END IF
END IF
!
!     Put the variables for the next calculation of the objective function
!       in XNEW, with any adjustments for the bounds.
!
!
!     Calculate the value of the objective function at XBASE+XNEW, unless
!       the limit on the number of calculations of F has been reached.
!
360 DO I=1,N
X(I)=DMIN1(DMAX1(XL(I),XBASE(I)+XNEW(I)),XU(I))
IF (XNEW(I) .EQ. SL(I)) X(I)=XL(I)
IF (XNEW(I) .EQ. SU(I)) X(I)=XU(I)
END DO
IF (NF .GE. MAXFUN) THEN
!JN          IF (IPRINT .GT. 0) CALL minqer(390)
!$$$          PRINT 390
!$$$  390     FORMAT (/4X,'Return from BOBYQA because CALFUN has been',
!$$$     1      ' called MAXFUN times.')
!$$$          GOTO 720
!JN 100807
IERR=390
GOTO 720
END IF
NF=NF+1
F = CALFUN (N,X,IPRINT)
!$$$      CALL minqi3(IPRINT, F, NF, N, X)
!$$$      IF (IPRINT .EQ. 3) THEN
!$$$          PRINT 400, NF,F,(X(I),I=1,N)
!$$$  400      FORMAT (/4X,'Function number',I6,'    F =',1PD18.10,
!$$$     1       '    The corresponding X is:'/(2X,5D15.6))
!$$$      END IF
IF (NTRITS .EQ. -1) THEN
FSAVE=F
GOTO 720
END IF
!
!     Use the quadratic model to predict the change in F due to the step D,
!       and set DIFF to the error of this prediction.
!
FOPT=FVAL(KOPT)
VQUAD=ZERO
IH=0
DO J=1,N
VQUAD=VQUAD+D(J)*GOPT(J)
DO I=1,J
IH=IH+1
TEMP=D(I)*D(J)
IF (I .EQ. J) TEMP=HALF*TEMP
VQUAD=VQUAD+HQ(IH)*TEMP
END DO
END DO
DO K=1,NPT
VQUAD=VQUAD+HALF*PQ(K)*W(NPT+K)**2
END DO
DIFF=F-FOPT-VQUAD
DIFFC=DIFFB
DIFFB=DIFFA
DIFFA=DABS(DIFF)
IF (DNORM .GT. RHO) NFSAV=NF
!
!     Pick the next value of DELTA after a trust region step.
!
IF (NTRITS .GT. 0) THEN
IF (VQUAD .GE. ZERO) THEN
!JN              IF (IPRINT .GT. 0) CALL minqer(430)
!$$$              PRINT 430
!$$$  430         FORMAT (/4X,'Return from BOBYQA because a trust',
!$$$     1          ' region step has failed to reduce Q.')
!$$$              GOTO 720
!JN 100807
IERR=430
GOTO 720
END IF
RATIO=(F-FOPT)/VQUAD
IF (RATIO .LE. TENTH) THEN
DELTA=DMIN1(HALF*DELTA,DNORM)
ELSE IF (RATIO.LE. 0.7D0) THEN
DELTA=DMAX1(HALF*DELTA,DNORM)
ELSE
DELTA=DMAX1(HALF*DELTA,DNORM+DNORM)
END IF
IF (DELTA .LE. 1.5D0*RHO) DELTA=RHO
!
!     Recalculate KNEW and DENOM if the new F is less than FOPT.
!
IF (F .LT. FOPT) THEN
KSAV=KNEW
DENSAV=DENOM
DELSQ=DELTA*DELTA
SCADEN=ZERO
BIGLSQ=ZERO
KNEW=0
DO K=1,NPT
HDIAG=ZERO
DO JJ=1,NPTM
HDIAG=HDIAG+ZMAT(K,JJ)**2
END DO
DEN=BETA*HDIAG+VLAG(K)**2
DISTSQ=ZERO
DO J=1,N
DISTSQ=DISTSQ+(XPT(K,J)-XNEW(J))**2
END DO
TEMP=DMAX1(ONE,(DISTSQ/DELSQ)**2)
IF (TEMP*DEN .GT. SCADEN) THEN
SCADEN=TEMP*DEN
KNEW=K
DENOM=DEN
END IF
BIGLSQ=DMAX1(BIGLSQ,TEMP*VLAG(K)**2)
END DO
IF (SCADEN .LE. HALF*BIGLSQ) THEN
KNEW=KSAV
DENOM=DENSAV
END IF
END IF
END IF
!
!     Update BMAT and ZMAT, so that the KNEW-th interpolation point can be
!     moved. Also update the second derivative terms of the model.
!
CALL UPDATEBOBYQA (N,NPT,BMAT,ZMAT,NDIM,VLAG,BETA,DENOM,KNEW,W)
IH=0
PQOLD=PQ(KNEW)
PQ(KNEW)=ZERO
DO I=1,N
TEMP=PQOLD*XPT(KNEW,I)
DO J=1,I
IH=IH+1
HQ(IH)=HQ(IH)+TEMP*XPT(KNEW,J)
END DO
END DO
DO JJ=1,NPTM
TEMP=DIFF*ZMAT(KNEW,JJ)
DO K=1,NPT
PQ(K)=PQ(K)+TEMP*ZMAT(K,JJ)
END DO
END DO
!
!     Include the new interpolation point, and make the changes to GOPT at
!     the old XOPT that are caused by the updating of the quadratic model.
!
FVAL(KNEW)=F
DO I=1,N
XPT(KNEW,I)=XNEW(I)
W(I)=BMAT(KNEW,I)
END DO
DO K=1,NPT
SUMA=ZERO
DO JJ=1,NPTM
SUMA=SUMA+ZMAT(KNEW,JJ)*ZMAT(K,JJ)
END DO
SUMB=ZERO
DO J=1,N
SUMB=SUMB+XPT(K,J)*XOPT(J)
END DO
TEMP=SUMA*SUMB
DO I=1,N
W(I)=W(I)+TEMP*XPT(K,I)
END DO
END DO
DO I=1,N
GOPT(I)=GOPT(I)+DIFF*W(I)
END DO
!
!     Update XOPT, GOPT and KOPT if the new calculated F is less than FOPT.
!
IF (F .LT. FOPT) THEN
KOPT=KNEW
XOPTSQ=ZERO
IH=0
DO J=1,N
XOPT(J)=XNEW(J)
XOPTSQ=XOPTSQ+XOPT(J)**2
DO I=1,J
IH=IH+1
IF (I .LT. J) GOPT(J)=GOPT(J)+HQ(IH)*D(I)
GOPT(I)=GOPT(I)+HQ(IH)*D(J)
END DO
END DO
DO K=1,NPT
TEMP=ZERO
DO J=1,N
TEMP=TEMP+XPT(K,J)*D(J)
END DO
TEMP=PQ(K)*TEMP
DO I=1,N
GOPT(I)=GOPT(I)+TEMP*XPT(K,I)
END DO
END DO
END IF
!
!     Calculate the parameters of the least Frobenius norm interpolant to
!     the current data, the gradient of this interpolant at XOPT being put
!     into VLAG(NPT+I), I=1,2,...,N.
!
IF (NTRITS .GT. 0) THEN
DO K=1,NPT
VLAG(K)=FVAL(K)-FVAL(KOPT)
W(K)=ZERO
END DO
DO J=1,NPTM
SUM=ZERO
DO  K=1,NPT
SUM=SUM+ZMAT(K,J)*VLAG(K)
END DO
DO K=1,NPT
W(K)=W(K)+SUM*ZMAT(K,J)
END DO
END DO
DO K=1,NPT
SUM=ZERO
DO J=1,N
SUM=SUM+XPT(K,J)*XOPT(J)
END DO
W(K+NPT)=W(K)
W(K)=SUM*W(K)
END DO
GQSQ=ZERO
GISQ=ZERO
DO I=1,N
SUM=ZERO
DO K=1,NPT
SUM=SUM+BMAT(K,I)*VLAG(K)+XPT(K,I)*W(K)
END DO
IF (XOPT(I) .EQ. SL(I)) THEN
GQSQ=GQSQ+DMIN1(ZERO,GOPT(I))**2
GISQ=GISQ+DMIN1(ZERO,SUM)**2
ELSE IF (XOPT(I) .EQ. SU(I)) THEN
GQSQ=GQSQ+DMAX1(ZERO,GOPT(I))**2
GISQ=GISQ+DMAX1(ZERO,SUM)**2
ELSE
GQSQ=GQSQ+GOPT(I)**2
GISQ=GISQ+SUM*SUM
END IF
VLAG(NPT+I)=SUM
END DO
!
!     Test whether to replace the new quadratic model by the least Frobenius
!     norm interpolant, making the replacement if the test is satisfied.
!
ITEST=ITEST+1
IF (GQSQ .LT. TEN*GISQ) ITEST=0
IF (ITEST .GE. 3) THEN
DO I=1,MAX0(NPT,NH)
IF (I .LE. N) GOPT(I)=VLAG(NPT+I)
IF (I .LE. NPT) PQ(I)=W(NPT+I)
IF (I .LE. NH) HQ(I)=ZERO
ITEST=0
END DO
END IF
END IF
!
!     If a trust region step has provided a sufficient decrease in F, then
!     branch for another trust region calculation. The case NTRITS=0 occurs
!     when the new interpolation point was reached by an alternative step.
!
IF (NTRITS .EQ. 0) GOTO 60
IF (F .LE. FOPT+TENTH*VQUAD) GOTO 60
!
!     Alternatively, find out if the interpolation points are close enough
!       to the best point so far.
!
DISTSQ=DMAX1((TWO*DELTA)**2,(TEN*RHO)**2)
650 KNEW=0
DO K=1,NPT
SUM=ZERO
DO J=1,N
SUM=SUM+(XPT(K,J)-XOPT(J))**2
END DO
IF (SUM .GT. DISTSQ) THEN
KNEW=K
DISTSQ=SUM
END IF
END DO
!
!     If KNEW is positive, then ALTMOV finds alternative new positions for
!     the KNEW-th interpolation point within distance ADELT of XOPT. It is
!     reached via label 90. Otherwise, there is a branch to label 60 for
!     another trust region iteration, unless the calculations with the
!     current RHO are complete.
!
IF (KNEW .GT. 0) THEN
DIST=DSQRT(DISTSQ)
IF (NTRITS .EQ. -1) THEN
DELTA=DMIN1(TENTH*DELTA,HALF*DIST)
IF (DELTA .LE. 1.5D0*RHO) DELTA=RHO
END IF
NTRITS=0
ADELT=DMAX1(DMIN1(TENTH*DIST,DELTA),RHO)
DSQ=ADELT*ADELT
GOTO 90
END IF
IF (NTRITS .EQ. -1) GOTO 680
IF (RATIO .GT. ZERO) GOTO 60
IF (DMAX1(DELTA,DNORM) .GT. RHO) GOTO 60
!
!     The calculations with the current value of RHO are complete. Pick the
!       next values of RHO and DELTA.
!
680 IF (RHO .GT. RHOEND) THEN
DELTA=HALF*RHO
RATIO=RHO/RHOEND
IF (RATIO .LE. 16.0D0) THEN
RHO=RHOEND
ELSE IF (RATIO .LE. 250.0D0) THEN
RHO=DSQRT(RATIO)*RHOEND
ELSE
RHO=TENTH*RHO
END IF
DELTA=DMAX1(DELTA,RHO)
CALL minqit(IPRINT, RHO, NF, FVAL(KOPT), N, XBASE, XOPT)
!$$$          IF (IPRINT .GE. 2) THEN
!$$$              IF (IPRINT .GE. 3) PRINT 690
!$$$  690         FORMAT (5X)
!$$$              PRINT 700, RHO,NF
!$$$  700         FORMAT (/4X,'New RHO =',1PD11.4,5X,'Number of',
!$$$     1          ' function values =',I6)
!$$$              PRINT 710, FVAL(KOPT),(XBASE(I)+XOPT(I),I=1,N)
!$$$  710         FORMAT (4X,'Least value of F =',1PD23.15,9X,
!$$$     1          'The corresponding X is:'/(2X,5D15.6))
!$$$          END IF
NTRITS=0
NFSAV=NF
GOTO 60
END IF
!
!     Return from the calculation, after another Newton-Raphson step, if
!       it is too short to have been tried before.
!
IF (NTRITS .EQ. -1) GOTO 360
720 IF (FVAL(KOPT) .LE. FSAVE) THEN
DO I=1,N
X(I)=DMIN1(DMAX1(XL(I),XBASE(I)+XOPT(I)),XU(I))
IF (XOPT(I) .EQ. SL(I)) X(I)=XL(I)
IF (XOPT(I) .EQ. SU(I)) X(I)=XU(I)
END DO
F=FVAL(KOPT)
END IF
!     JN 100807 Do we want to add IERR to minqir as a diagnostic. If zero, not print,
!JN          if not, then use minqer output or similar.
!JN ??      IF (IERR.NE.0) CALL minqer(IERR)
CALL minqir(IPRINT, F, NF, N, X)
!$$$  IF (IPRINT .GE. 1) THEN
!$$$          PRINT 740, NF
!$$$  740     FORMAT (/4X,'At the return from BOBYQA',5X,
!$$$  1      'Number of function values =',I6)
!$$$          PRINT 710, F,(X(I),I=1,N)
!$$$  END IF
RETURN
END


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% lagmax.f %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE LAGMAX (N,G,H,RHO,D,V,VMAX)
implicit integer (i-n), real(real64) (a-h,o-z)
DIMENSION G(*),H(N,*),D(*),V(*)
!
!     N is the number of variables of a quadratic objective function, Q say.
!     G is the gradient of Q at the origin.
!     H is the symmetric Hessian matrix of Q. Only the upper triangular and
!       diagonal parts need be set.
!     RHO is the trust region radius, and has to be positive.
!     D will be set to the calculated vector of variables.
!     The array V will be used for working space.
!     VMAX will be set to |Q(0)-Q(D)|.
!
!     Calculating the D that maximizes |Q(0)-Q(D)| subject to ||D|| .LEQ. RHO
!     requires of order N**3 operations, but sometimes it is adequate if
!     |Q(0)-Q(D)| is within about 0.9 of its greatest possible value. This
!     subroutine provides such a solution in only of order N**2 operations,
!     where the claim of accuracy has been tested by numerical experiments.
!
!     Preliminary calculations.
!
HALF=0.5D0
HALFRT=DSQRT(HALF)
ONE=1.0D0
ZERO=0.0D0
!
!     Pick V such that ||HV|| / ||V|| is large.
!
HMAX=ZERO
DO I=1,N
SUM=ZERO
DO J=1,N
H(J,I)=H(I,J)
SUM=SUM+H(I,J)**2
END DO
IF (SUM .GT. HMAX) THEN
HMAX=SUM
K=I
END IF
END DO
DO J=1,N
V(J)=H(K,J)
END DO
!
!     Set D to a vector in the subspace spanned by V and HV that maximizes
!     |(D,HD)|/(D,D), except that we set D=HV if V and HV are nearly parallel.
!     The vector that has the name D at label 60 used to be the vector W.
!
VSQ=ZERO
VHV=ZERO
DSQ=ZERO
DO I=1,N
VSQ=VSQ+V(I)**2
D(I)=ZERO
DO J=1,N
D(I)=D(I)+H(I,J)*V(J)
END DO
VHV=VHV+V(I)*D(I)
DSQ=DSQ+D(I)**2
END DO
IF (VHV*VHV .LE. 0.9999D0*DSQ*VSQ) THEN
TEMP=VHV/VSQ
WSQ=ZERO
DO I=1,N
D(I)=D(I)-TEMP*V(I)
WSQ=WSQ+D(I)**2
END DO
WHW=ZERO
RATIO=DSQRT(WSQ/VSQ)
DO I=1,N
TEMP=ZERO
DO J=1,N
TEMP=TEMP+H(I,J)*D(J)
END DO
WHW=WHW+TEMP*D(I)
V(I)=RATIO*V(I)
END DO
VHV=RATIO*RATIO*VHV
VHW=RATIO*WSQ
TEMP=HALF*(WHW-VHV)
TEMP=TEMP+DSIGN(DSQRT(TEMP**2+VHW**2),WHW+VHV)
DO I=1,N
D(I)=VHW*V(I)+TEMP*D(I)
END DO
END IF
!
!     We now turn our attention to the subspace spanned by G and D. A multiple
!     of the current D is returned if that choice seems to be adequate.
!
GG=ZERO
GD=ZERO
DD=ZERO
DHD=ZERO
DO I=1,N
GG=GG+G(I)**2
GD=GD+G(I)*D(I)
DD=DD+D(I)**2
SUM=ZERO
DO J=1,N
SUM=SUM+H(I,J)*D(J)
END DO
DHD=DHD+SUM*D(I)
END DO
TEMP=GD/GG
VV=ZERO
SCALE=DSIGN(RHO/DSQRT(DD),GD*DHD)
DO I=1,N
V(I)=D(I)-TEMP*G(I)
VV=VV+V(I)**2
D(I)=SCALE*D(I)
END DO
GNORM=DSQRT(GG)
IF (GNORM*DD .LE. 0.5D-2*RHO*DABS(DHD) .OR. VV/DD .LE. 1.0D-4) THEN
VMAX=DABS(SCALE*(GD+HALF*SCALE*DHD))
GOTO 170
END IF
!
!     G and V are now orthogonal in the subspace spanned by G and D. Hence
!     we generate an orthonormal basis of this subspace such that (D,HV) is
!     negligible or zero, where D and V will be the basis vectors.
!
GHG=ZERO
VHG=ZERO
VHV=ZERO
DO I=1,N
SUM=ZERO
SUMV=ZERO
DO J=1,N
SUM=SUM+H(I,J)*G(J)
SUMV=SUMV+H(I,J)*V(J)
END DO
GHG=GHG+SUM*G(I)
VHG=VHG+SUMV*G(I)
VHV=VHV+SUMV*V(I)
END DO
VNORM=DSQRT(VV)
GHG=GHG/GG
VHG=VHG/(VNORM*GNORM)
VHV=VHV/VV
IF (DABS(VHG) .LE. 0.01D0*DMAX1(DABS(GHG),DABS(VHV))) THEN
VMU=GHG-VHV
WCOS=ONE
WSIN=ZERO
ELSE
TEMP=HALF*(GHG-VHV)
VMU=TEMP+DSIGN(DSQRT(TEMP**2+VHG**2),TEMP)
TEMP=DSQRT(VMU**2+VHG**2)
WCOS=VMU/TEMP
WSIN=VHG/TEMP
END IF
TEMPA=WCOS/GNORM
TEMPB=WSIN/VNORM
TEMPC=WCOS/VNORM
TEMPD=WSIN/GNORM
DO I=1,N
D(I)=TEMPA*G(I)+TEMPB*V(I)
V(I)=TEMPC*V(I)-TEMPD*G(I)
END DO
!
!     The final D is a multiple of the current D, V, D+V or D-V. We make the
!     choice from these possibilities that is optimal.
!
DLIN=WCOS*GNORM/RHO
VLIN=-WSIN*GNORM/RHO
TEMPA=DABS(DLIN)+HALF*DABS(VMU+VHV)
TEMPB=DABS(VLIN)+HALF*DABS(GHG-VMU)
TEMPC=HALFRT*(DABS(DLIN)+DABS(VLIN))+0.25D0*DABS(GHG+VHV)
IF (TEMPA .GE. TEMPB .AND. TEMPA .GE. TEMPC) THEN
TEMPD=DSIGN(RHO,DLIN*(VMU+VHV))
TEMPV=ZERO
ELSE IF (TEMPB .GE. TEMPC) THEN
TEMPD=ZERO
TEMPV=DSIGN(RHO,VLIN*(GHG-VMU))
ELSE
TEMPD=DSIGN(HALFRT*RHO,DLIN*(GHG+VHV))
TEMPV=DSIGN(HALFRT*RHO,VLIN*(GHG+VHV))
END IF
DO I=1,N
D(I)=TEMPD*D(I)+TEMPV*V(I)
END DO
VMAX=RHO*RHO*DMAX1(TEMPA,TEMPB,TEMPC)
170 RETURN
END


subroutine newuoa_core (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W,IERR)
implicit integer (i-n), real(real64) (a-h,o-z)
DIMENSION X(*),W(*)
!JN Declare error flag
INTEGER IERR
!
!     This subroutine seeks the least value of a function of many variables,
!     by a trust region method that forms quadratic models by interpolation.
!     There can be some freedom in the interpolation conditions, which is
!     taken up by minimizing the Frobenius norm of the change to the second
!     derivative of the quadratic model, beginning with a zero matrix. The
!     arguments of the subroutine are as follows.
!
!     N must be set to the number of variables and must be at least two.
!     NPT is the number of interpolation conditions. Its value must be in the
!       interval [N+2,(N+1)(N+2)/2].
!     Initial values of the variables must be set in X(1),X(2),...,X(N). They
!       will be changed to the values that give the least calculated F.
!     RHOBEG and RHOEND must be set to the initial and final values of a trust
!       region radius, so both must be positive with RHOEND<=RHOBEG. Typically
!       RHOBEG should be about one tenth of the greatest expected change to a
!       variable, and RHOEND should indicate the accuracy that is required in
!       the final values of the variables.
!     The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
!       amount of printing. Specifically, there is no output if IPRINT=0 and
!       there is output only at the return if IPRINT=1. Otherwise, each new
!       value of RHO is printed, with the best vector of variables so far and
!       the corresponding value of the objective function. Further, each new
!       value of F with its variables are output if IPRINT=3.
!     MAXFUN must be set to an upper bound on the number of calls of CALFUN.
!     The array W will be used for working space. Its length must be at least
!     (NPT+13)*(NPT+N)+3*N*(N+3)/2.
!JN   Add IERR to tell what the error is to minqer.
!     IERR gives an error code that is interpreted by minqer interface routine.
!       and passed back to minqa.R
!
!     SUBROUTINE CALFUN (N,X,F) must be provided by the user. It must set F to
!     the value of the objective function for the variables X(1),X(2),...,X(N).
!
!     Partition the working space array, so that different parts of it can be
!     treated separately by the subroutine that performs the main calculation.
!
NP=N+1
NPTM=NPT-NP
IF (NPT .LT. N+2 .OR. NPT .GT. ((N+2)*NP)/2) THEN
!JN         CALL minqer(101)
!$$$          PRINT 10
!$$$   10     FORMAT (/4X,'Return from NEWUOA because NPT is not in',
!$$$     1      ' the required interval')
!$$$          GO TO 20
IERR = 10
GO TO 20
END IF
NDIM=NPT+N
IXB=1
IXO=IXB+N
IXN=IXO+N
IXP=IXN+N
IFV=IXP+N*NPT
IGQ=IFV+NPT
IHQ=IGQ+N
IPQ=IHQ+(N*NP)/2
IBMAT=IPQ+NPT
IZMAT=IBMAT+NDIM*N
ID=IZMAT+NPT*NPTM
IVL=ID+N
IW=IVL+NDIM
!
!     The above settings provide a partition of W for subroutine NEWUOB.
!     The partition requires the first NPT*(NPT+N)+5*N*(N+3)/2 elements of
!     W plus the space that is needed by the last array of NEWUOB.
!
CALL NEWUOB (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W(IXB), W(IXO),W(IXN),W(IXP),W(IFV),W(IGQ),W(IHQ), &
    & W(IPQ),W(IBMAT), W(IZMAT),NDIM,W(ID),W(IVL),W(IW), IERR)
20 RETURN
END


SUBROUTINE NEWUOB (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,XBASE, XOPT,XNEW,XPT,FVAL,GQ,HQ,PQ,BMAT,ZMAT,NDIM, &
    & D,VLAG,W,IERR)
!JN Remember IERR here.
implicit integer (i-n), real(real64) (a-h,o-z)
DIMENSION X(*),XBASE(*),XOPT(*),XNEW(*),XPT(NPT,*),FVAL(*), GQ(*),HQ(*),PQ(*),BMAT(NDIM,*),ZMAT(NPT,*), &
    & D(*),VLAG(*),W(*)
!
!     The arguments N, NPT, X, RHOBEG, RHOEND, IPRINT and MAXFUN are identical
!       to the corresponding arguments in SUBROUTINE NEWUOA.
!     XBASE will hold a shift of origin that should reduce the contributions
!       from rounding errors to values of the model and Lagrange functions.
!     XOPT will be set to the displacement from XBASE of the vector of
!       variables that provides the least calculated F so far.
!     XNEW will be set to the displacement from XBASE of the vector of
!       variables for the current calculation of F.
!     XPT will contain the interpolation point coordinates relative to XBASE.
!     FVAL will hold the values of F at the interpolation points.
!     GQ will hold the gradient of the quadratic model at XBASE.
!     HQ will hold the explicit second derivatives of the quadratic model.
!     PQ will contain the parameters of the implicit second derivatives of
!       the quadratic model.
!     BMAT will hold the last N columns of H.
!     ZMAT will hold the factorization of the leading NPT by NPT submatrix of
!       H, this factorization being ZMAT times Diag(DZ) times ZMAT^T, where
!       the elements of DZ are plus or minus one, as specified by IDZ.
!     NDIM is the first dimension of BMAT and has the value NPT+N.
!     D is reserved for trial steps from XOPT.
!     VLAG will contain the values of the Lagrange functions at a new point X.
!       They are part of a product that requires VLAG to be of length NDIM.
!     The array W will be used for working space. Its length must be at least
!       10*NDIM = 10*(NPT+N).
!
!     Set some constants.
!
HALF=0.5D0
ONE=1.0D0
TENTH=0.1D0
ZERO=0.0D0
NP=N+1
NH=(N*NP)/2
NPTM=NPT-NP
NFTEST=MAX0(MAXFUN,1)
!
!     Set the initial elements of XPT, BMAT, HQ, PQ and ZMAT to zero.
!
DO J=1,N
XBASE(J)=X(J)
DO K=1,NPT
XPT(K,J)=ZERO
END DO
DO I=1,NDIM
BMAT(I,J)=ZERO
END DO
END DO
DO IH=1,NH
HQ(IH)=ZERO
END DO
DO K=1,NPT
PQ(K)=ZERO
DO J=1,NPTM
ZMAT(K,J)=ZERO
END DO
END DO
!
!     Begin the initialization procedure. NF becomes one more than the number
!     of function values so far. The coordinates of the displacement of the
!     next initial interpolation point from XBASE are set in XPT(NF,.).
!
RHOSQ=RHOBEG*RHOBEG
RECIP=ONE/RHOSQ
RECIQ=DSQRT(HALF)/RHOSQ
NF=0
50 NFM=NF
NFMM=NF-N
NF=NF+1
IF (NFM .LE. 2*N) THEN
IF (NFM .GE. 1 .AND. NFM .LE. N) THEN
XPT(NF,NFM)=RHOBEG
ELSE IF (NFM .GT. N) THEN
XPT(NF,NFMM)=-RHOBEG
END IF
ELSE
ITEMP=(NFMM-1)/N
JPT=NFM-ITEMP*N-N
IPT=JPT+ITEMP
IF (IPT .GT. N) THEN
ITEMP=JPT
JPT=IPT-N
IPT=ITEMP
END IF
XIPT=RHOBEG
IF (FVAL(IPT+NP) .LT. FVAL(IPT+1)) XIPT=-XIPT
XJPT=RHOBEG
IF (FVAL(JPT+NP) .LT. FVAL(JPT+1)) XJPT=-XJPT
XPT(NF,IPT)=XIPT
XPT(NF,JPT)=XJPT
END IF
!
!     Calculate the next value of F, label 70 being reached immediately
!     after this calculation. The least function value so far and its index
!     are required.
!
DO J=1,N
X(J)=XPT(NF,J)+XBASE(J)
END DO
GOTO 310
70 FVAL(NF)=F
IF (NF .EQ. 1) THEN
FBEG=F
FOPT=F
KOPT=1
ELSE IF (F .LT. FOPT) THEN
FOPT=F
KOPT=NF
END IF
!
!     Set the nonzero initial elements of BMAT and the quadratic model in
!     the cases when NF is at most 2*N+1.
!
IF (NFM .LE. 2*N) THEN
IF (NFM .GE. 1 .AND. NFM .LE. N) THEN
GQ(NFM)=(F-FBEG)/RHOBEG
IF (NPT .LT. NF+N) THEN
BMAT(1,NFM)=-ONE/RHOBEG
BMAT(NF,NFM)=ONE/RHOBEG
BMAT(NPT+NFM,NFM)=-HALF*RHOSQ
END IF
ELSE IF (NFM .GT. N) THEN
BMAT(NF-N,NFMM)=HALF/RHOBEG
BMAT(NF,NFMM)=-HALF/RHOBEG
ZMAT(1,NFMM)=-RECIQ-RECIQ
ZMAT(NF-N,NFMM)=RECIQ
ZMAT(NF,NFMM)=RECIQ
IH=(NFMM*(NFMM+1))/2
TEMP=(FBEG-F)/RHOBEG
HQ(IH)=(GQ(NFMM)-TEMP)/RHOBEG
GQ(NFMM)=HALF*(GQ(NFMM)+TEMP)
END IF
!
!     Set the off-diagonal second derivatives of the Lagrange functions and
!     the initial quadratic model.
!
ELSE
IH=(IPT*(IPT-1))/2+JPT
IF (XIPT .LT. ZERO) IPT=IPT+N
IF (XJPT .LT. ZERO) JPT=JPT+N
ZMAT(1,NFMM)=RECIP
ZMAT(NF,NFMM)=RECIP
ZMAT(IPT+1,NFMM)=-RECIP
ZMAT(JPT+1,NFMM)=-RECIP
HQ(IH)=(FBEG-FVAL(IPT+1)-FVAL(JPT+1)+F)/(XIPT*XJPT)
END IF
IF (NF .LT. NPT) GOTO 50
!
!     Begin the iterative procedure, because the initial model is complete.
!
RHO=RHOBEG
DELTA=RHO
IDZ=1
DIFFA=ZERO
DIFFB=ZERO
ITEST=0
XOPTSQ=ZERO
DO I=1,N
XOPT(I)=XPT(KOPT,I)
XOPTSQ=XOPTSQ+XOPT(I)**2
END DO
90 NFSAV=NF
!
!     Generate the next trust region step and test its length. Set KNEW
!     to -1 if the purpose of the next F will be to improve the model.
!
100 KNEW=0
CALL TRSAPP (N,NPT,XOPT,XPT,GQ,HQ,PQ,DELTA,D,W,W(NP), W(NP+N),W(NP+2*N),CRVMIN)
DSQ=ZERO
DO I=1,N
DSQ=DSQ+D(I)**2
END DO
DNORM=DMIN1(DELTA,DSQRT(DSQ))
IF (DNORM .LT. HALF*RHO) THEN
KNEW=-1
DELTA=TENTH*DELTA
RATIO=-1.0D0
IF (DELTA .LE. 1.5D0*RHO) DELTA=RHO
IF (NF .LE. NFSAV+2) GOTO 460
TEMP=0.125D0*CRVMIN*RHO*RHO
IF (TEMP .LE. DMAX1(DIFFA,DIFFB,DIFFC)) GOTO 460
GOTO 490
END IF
!
!     Shift XBASE if XOPT may be too far from XBASE. First make the changes
!     to BMAT that do not depend on ZMAT.
!
120 IF (DSQ .LE. 1.0D-3*XOPTSQ) THEN
TEMPQ=0.25D0*XOPTSQ
DO K=1,NPT
SUM=ZERO
DO I=1,N
SUM=SUM+XPT(K,I)*XOPT(I)
END DO
TEMP=PQ(K)*SUM
SUM=SUM-HALF*XOPTSQ
W(NPT+K)=SUM
DO I=1,N
GQ(I)=GQ(I)+TEMP*XPT(K,I)
XPT(K,I)=XPT(K,I)-HALF*XOPT(I)
VLAG(I)=BMAT(K,I)
W(I)=SUM*XPT(K,I)+TEMPQ*XOPT(I)
IP=NPT+I
DO J=1,I
BMAT(IP,J)=BMAT(IP,J)+VLAG(I)*W(J)+W(I)*VLAG(J)
END DO
END DO
END DO
!
!     Then the revisions of BMAT that depend on ZMAT are calculated.
!
DO K=1,NPTM
SUMZ=ZERO
DO I=1,NPT
SUMZ=SUMZ+ZMAT(I,K)
W(I)=W(NPT+I)*ZMAT(I,K)
END DO
DO J=1,N
SUM=TEMPQ*SUMZ*XOPT(J)
DO I=1,NPT
SUM=SUM+W(I)*XPT(I,J)
END DO
VLAG(J)=SUM
IF (K .LT. IDZ) SUM=-SUM
DO I=1,NPT
BMAT(I,J)=BMAT(I,J)+SUM*ZMAT(I,K)
END DO
END DO
DO I=1,N
IP=I+NPT
TEMP=VLAG(I)
IF (K .LT. IDZ) TEMP=-TEMP
DO J=1,I
BMAT(IP,J)=BMAT(IP,J)+TEMP*VLAG(J)
END DO
END DO
END DO
!
!     The following instructions complete the shift of XBASE, including
!     the changes to the parameters of the quadratic model.
!
IH=0
DO J=1,N
W(J)=ZERO
DO K=1,NPT
W(J)=W(J)+PQ(K)*XPT(K,J)
XPT(K,J)=XPT(K,J)-HALF*XOPT(J)
END DO
DO I=1,J
IH=IH+1
IF (I .LT. J) GQ(J)=GQ(J)+HQ(IH)*XOPT(I)
GQ(I)=GQ(I)+HQ(IH)*XOPT(J)
HQ(IH)=HQ(IH)+W(I)*XOPT(J)+XOPT(I)*W(J)
BMAT(NPT+I,J)=BMAT(NPT+J,I)
END DO
END DO
DO J=1,N
XBASE(J)=XBASE(J)+XOPT(J)
XOPT(J)=ZERO
END DO
XOPTSQ=ZERO
END IF
!
!     Pick the model step if KNEW is positive. A different choice of D
!     may be made later, if the choice of D by BIGLAG causes substantial
!     cancellation in DENOM.
!
IF (KNEW .GT. 0) THEN
CALL BIGLAG (N,NPT,XOPT,XPT,BMAT,ZMAT,IDZ,NDIM,KNEW,DSTEP, D,ALPHA,VLAG,VLAG(NPT+1),W,W(NP),W(NP+N))
END IF
!
!     Calculate VLAG and BETA for the current choice of D. The first NPT
!     components of W_check will be held in W.
!
DO K=1,NPT
SUMA=ZERO
SUMB=ZERO
SUM=ZERO
DO J=1,N
SUMA=SUMA+XPT(K,J)*D(J)
SUMB=SUMB+XPT(K,J)*XOPT(J)
SUM=SUM+BMAT(K,J)*D(J)
END DO
W(K)=SUMA*(HALF*SUMA+SUMB)
VLAG(K)=SUM
END DO
BETA=ZERO
DO K=1,NPTM
SUM=ZERO
DO I=1,NPT
SUM=SUM+ZMAT(I,K)*W(I)
END DO
IF (K .LT. IDZ) THEN
BETA=BETA+SUM*SUM
SUM=-SUM
ELSE
BETA=BETA-SUM*SUM
END IF
DO I=1,NPT
VLAG(I)=VLAG(I)+SUM*ZMAT(I,K)
END DO
END DO
BSUM=ZERO
DX=ZERO
DO J=1,N
SUM=ZERO
DO I=1,NPT
SUM=SUM+W(I)*BMAT(I,J)
END DO
BSUM=BSUM+SUM*D(J)
JP=NPT+J
DO K=1,N
SUM=SUM+BMAT(JP,K)*D(K)
END DO
VLAG(JP)=SUM
BSUM=BSUM+SUM*D(J)
DX=DX+D(J)*XOPT(J)
END DO
BETA=DX*DX+DSQ*(XOPTSQ+DX+DX+HALF*DSQ)+BETA-BSUM
VLAG(KOPT)=VLAG(KOPT)+ONE
!
!     If KNEW is positive and if the cancellation in DENOM is unacceptable,
!     then BIGDEN calculates an alternative model step, XNEW being used for
!     working space.
!
IF (KNEW .GT. 0) THEN
TEMP=ONE+ALPHA*BETA/VLAG(KNEW)**2
IF (DABS(TEMP) .LE. 0.8D0) THEN
CALL BIGDEN (N,NPT,XOPT,XPT,BMAT,ZMAT,IDZ,NDIM,KOPT, KNEW,D,W,VLAG,BETA,XNEW,W(NDIM+1),W(6*NDIM+1))
END IF
END IF
!
!     Calculate the next value of the objective function.
!
290 DO I=1,N
XNEW(I)=XOPT(I)+D(I)
X(I)=XBASE(I)+XNEW(I)
END DO
NF=NF+1
310 IF (NF .GT. NFTEST) THEN
NF=NF-1
!JN          CALL MINQER (390)
!$$$          IF (IPRINT .GT. 0) PRINT 320
!$$$  320     FORMAT (/4X,'Return from NEWUOA because CALFUN has been',
!$$$     1      ' called MAXFUN times.')
!$$$          GOTO 530
IERR = 390
GO TO 530
END IF
F = CALFUN (N,X,IPRINT)
!$$$      IF (IPRINT .EQ. 3) THEN
!$$$  PRINT 70, NF,F,(X(I),I=1,N)
!$$$  70      FORMAT (/4X,'Function number',I6,'    F =',1PD18.10,
!$$$  1       '    The corresponding X is:'/(2X,5D15.6))
!$$$      END IF
!$$$      CALL minqi3 (IPRINT, F, NF, N, X)
IF (NF .LE. NPT) GOTO 70
IF (KNEW .EQ. -1) GOTO 530
!
!     Use the quadratic model to predict the change in F due to the step D,
!     and set DIFF to the error of this prediction.
!
VQUAD=ZERO
IH=0
DO J=1,N
VQUAD=VQUAD+D(J)*GQ(J)
DO I=1,J
IH=IH+1
TEMP=D(I)*XNEW(J)+D(J)*XOPT(I)
IF (I .EQ. J) TEMP=HALF*TEMP
VQUAD=VQUAD+TEMP*HQ(IH)
END DO
END DO
DO K=1,NPT
VQUAD=VQUAD+PQ(K)*W(K)
END DO
DIFF=F-FOPT-VQUAD
DIFFC=DIFFB
DIFFB=DIFFA
DIFFA=DABS(DIFF)
IF (DNORM .GT. RHO) NFSAV=NF
!
!     Update FOPT and XOPT if the new F is the least value of the objective
!     function so far. The branch when KNEW is positive occurs if D is not
!     a trust region step.
!
FSAVE=FOPT
IF (F .LT. FOPT) THEN
FOPT=F
XOPTSQ=ZERO
DO I=1,N
XOPT(I)=XNEW(I)
XOPTSQ=XOPTSQ+XOPT(I)**2
END DO
END IF
KSAVE=KNEW
IF (KNEW .GT. 0) GOTO 410
!
!     Pick the next value of DELTA after a trust region step.
!
IF (VQUAD .GE. ZERO) THEN
!JN              IF (IPRINT .GT. 0) CALL minqer(3701)
!$$$          IF (IPRINT .GT. 0) PRINT 370
!$$$  370     FORMAT (/4X,'Return from NEWUOA because a trust',
!$$$     1      ' region step has failed to reduce Q.')
IERR = 3701
GOTO 530
END IF
RATIO=(F-FSAVE)/VQUAD
IF (RATIO .LE. TENTH) THEN
DELTA=HALF*DNORM
ELSE IF (RATIO.LE. 0.7D0) THEN
DELTA=DMAX1(HALF*DELTA,DNORM)
ELSE
DELTA=DMAX1(HALF*DELTA,DNORM+DNORM)
END IF
IF (DELTA .LE. 1.5D0*RHO) DELTA=RHO
!
!     Set KNEW to the index of the next interpolation point to be deleted.
!
RHOSQ=DMAX1(TENTH*DELTA,RHO)**2
KTEMP=0
DETRAT=ZERO
IF (F .GE. FSAVE) THEN
KTEMP=KOPT
DETRAT=ONE
END IF
DO K=1,NPT
HDIAG=ZERO
DO J=1,NPTM
TEMP=ONE
IF (J .LT. IDZ) TEMP=-ONE
HDIAG=HDIAG+TEMP*ZMAT(K,J)**2
END DO
TEMP=DABS(BETA*HDIAG+VLAG(K)**2)
DISTSQ=ZERO
DO J=1,N
DISTSQ=DISTSQ+(XPT(K,J)-XOPT(J))**2
END DO
IF (DISTSQ .GT. RHOSQ) TEMP=TEMP*(DISTSQ/RHOSQ)**3
IF (TEMP .GT. DETRAT .AND. K .NE. KTEMP) THEN
DETRAT=TEMP
KNEW=K
END IF
END DO
IF (KNEW .EQ. 0) GOTO 460
!
!     Update BMAT, ZMAT and IDZ, so that the KNEW-th interpolation point
!     can be moved. Begin the updating of the quadratic model, starting
!     with the explicit second derivative term.
!
410 CALL UPDATE (N,NPT,BMAT,ZMAT,IDZ,NDIM,VLAG,BETA,KNEW,W)
FVAL(KNEW)=F
IH=0
DO I=1,N
TEMP=PQ(KNEW)*XPT(KNEW,I)
DO J=1,I
IH=IH+1
HQ(IH)=HQ(IH)+TEMP*XPT(KNEW,J)
END DO
END DO
PQ(KNEW)=ZERO
!
!     Update the other second derivative parameters, and then the gradient
!     vector of the model. Also include the new interpolation point.
!
DO J=1,NPTM
TEMP=DIFF*ZMAT(KNEW,J)
IF (J .LT. IDZ) TEMP=-TEMP
DO K=1,NPT
PQ(K)=PQ(K)+TEMP*ZMAT(K,J)
END DO
END DO
GQSQ=ZERO
DO I=1,N
GQ(I)=GQ(I)+DIFF*BMAT(KNEW,I)
GQSQ=GQSQ+GQ(I)**2
XPT(KNEW,I)=XNEW(I)
END DO
!
!     If a trust region step makes a small change to the objective function,
!     then calculate the gradient of the least Frobenius norm interpolant at
!     XBASE, and store it in W, using VLAG for a vector of right hand sides.
!
IF (KSAVE .EQ. 0 .AND. DELTA .EQ. RHO) THEN
IF (DABS(RATIO) .GT. 1.0D-2) THEN
ITEST=0
ELSE
DO K=1,NPT
VLAG(K)=FVAL(K)-FVAL(KOPT)
END DO
GISQ=ZERO
DO I=1,N
SUM=ZERO
DO K=1,NPT
SUM=SUM+BMAT(K,I)*VLAG(K)
END DO
GISQ=GISQ+SUM*SUM
W(I)=SUM
END DO
!
!     Test whether to replace the new quadratic model by the least Frobenius
!     norm interpolant, making the replacement if the test is satisfied.
!
ITEST=ITEST+1
IF (GQSQ .LT. 1.0D2*GISQ) ITEST=0
IF (ITEST .GE. 3) THEN
DO I=1,N
GQ(I)=W(I)
END DO
DO IH=1,NH
HQ(IH)=ZERO
END DO
DO J=1,NPTM
W(J)=ZERO
DO K=1,NPT
W(J)=W(J)+VLAG(K)*ZMAT(K,J)
END DO
IF (J .LT. IDZ) W(J)=-W(J)
END DO
DO K=1,NPT
PQ(K)=ZERO
DO J=1,NPTM
PQ(K)=PQ(K)+ZMAT(K,J)*W(J)
END DO
END DO
ITEST=0
END IF
END IF
END IF
IF (F .LT. FSAVE) KOPT=KNEW
!
!     If a trust region step has provided a sufficient decrease in F, then
!     branch for another trust region calculation. The case KSAVE>0 occurs
!     when the new function value was calculated by a model step.
!
IF (F .LE. FSAVE+TENTH*VQUAD) GOTO 100
IF (KSAVE .GT. 0) GOTO 100
!
!     Alternatively, find out if the interpolation points are close enough
!     to the best point so far.
!
KNEW=0
460 DISTSQ=4.0D0*DELTA*DELTA
DO K=1,NPT
SUM=ZERO
DO J=1,N
SUM=SUM+(XPT(K,J)-XOPT(J))**2
END DO
IF (SUM .GT. DISTSQ) THEN
KNEW=K
DISTSQ=SUM
END IF
END DO
!
!     If KNEW is positive, then set DSTEP, and branch back for the next
!     iteration, which will generate a "model step".
!
IF (KNEW .GT. 0) THEN
DSTEP=DMAX1(DMIN1(TENTH*DSQRT(DISTSQ),HALF*DELTA),RHO)
DSQ=DSTEP*DSTEP
GOTO 120
END IF
IF (RATIO .GT. ZERO) GOTO 100
IF (DMAX1(DELTA,DNORM) .GT. RHO) GOTO 100
!
!     The calculations with the current value of RHO are complete. Pick the
!     next values of RHO and DELTA.
!
490 IF (RHO .GT. RHOEND) THEN
DELTA=HALF*RHO
RATIO=RHO/RHOEND
IF (RATIO .LE. 16.0D0) THEN
RHO=RHOEND
ELSE IF (RATIO .LE. 250.0D0) THEN
RHO=DSQRT(RATIO)*RHOEND
ELSE
RHO=TENTH*RHO
END IF
DELTA=DMAX1(DELTA,RHO)
IF (IPRINT .GE. 2) THEN
CALL minqit(IPRINT, RHO, NF, FOPT, N, XBASE, XOPT)
!$$$              IF (IPRINT .GE. 3) PRINT 500
!$$$  500         FORMAT (5X)
!$$$              PRINT 510, RHO,NF
!$$$  510         FORMAT (/4X,'New RHO =',1PD11.4,5X,'Number of',
!$$$     1          ' function values =',I6)
!$$$              PRINT 520, FOPT,(XBASE(I)+XOPT(I),I=1,N)
!$$$  520         FORMAT (4X,'Least value of F =',1PD23.15,9X,
!$$$     1          'The corresponding X is:'/(2X,5D15.6))
END IF
GOTO 90
END IF
!
!     Return from the calculation, after another Newton-Raphson step, if
!     it is too short to have been tried before.
!
IF (KNEW .EQ. -1) GOTO 290
530 IF (FOPT .LE. F) THEN
DO I=1,N
X(I)=XBASE(I)+XOPT(I)
END DO
F=FOPT
END IF
IF (IPRINT .GE. 1) THEN
CALL minqir(IPRINT, F, NF, N, X)
!$$$          PRINT 550, NF
!$$$  550     FORMAT (/4X,'At the return from NEWUOA',5X,
!$$$     1      'Number of function values =',I6)
!$$$          PRINT 520, F,(X(I),I=1,N)
END IF
RETURN
END


SUBROUTINE PRELIM (N,NPT,X,XL,XU,RHOBEG,IPRINT,MAXFUN,XBASE, XPT,FVAL,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,SL,SU,NF, &
    & KOPT)
implicit integer (i-n), real(real64) (a-h,o-z)
DIMENSION X(*),XL(*),XU(*),XBASE(*),XPT(NPT,*),FVAL(*),GOPT(*), HQ(*),PQ(*),BMAT(NDIM,*),ZMAT(NPT,*), &
    & SL(*),SU(*)
!
!     The arguments N, NPT, X, XL, XU, RHOBEG, IPRINT and MAXFUN are the
!       same as the corresponding arguments in SUBROUTINE BOBYQA.
!     The arguments XBASE, XPT, FVAL, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU
!       are the same as the corresponding arguments in BOBYQB, the elements
!       of SL and SU being set in BOBYQA.
!     GOPT is usually the gradient of the quadratic model at XOPT+XBASE, but
!       it is set by PRELIM to the gradient of the quadratic model at XBASE.
!       If XOPT is nonzero, BOBYQB will change it to its usual value later.
!     NF is maintaned as the number of calls of CALFUN so far.
!     KOPT will be such that the least calculated value of F so far is at
!       the point XPT(KOPT,.)+XBASE in the space of the variables.
!
!     SUBROUTINE PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ,
!     BMAT and ZMAT for the first iteration, and it maintains the values of
!     NF and KOPT. The vector X is also changed by PRELIM.
!
!     Set some constants.
!
HALF=0.5D0
ONE=1.0D0
TWO=2.0D0
ZERO=0.0D0
RHOSQ=RHOBEG*RHOBEG
RECIP=ONE/RHOSQ
NP=N+1
!
!     Set XBASE to the initial vector of variables, and set the initial
!     elements of XPT, BMAT, HQ, PQ and ZMAT to zero.
!
DO J=1,N
XBASE(J)=X(J)
DO K=1,NPT
XPT(K,J)=ZERO
END DO
DO I=1,NDIM
BMAT(I,J)=ZERO
END DO
END DO
DO IH=1,(N*NP)/2
HQ(IH)=ZERO
END DO
DO K=1,NPT
PQ(K)=ZERO
DO J=1,NPT-NP
ZMAT(K,J)=ZERO
END DO
END DO
!
!     Begin the initialization procedure. NF becomes one more than the number
!     of function values so far. The coordinates of the displacement of the
!     next initial interpolation point from XBASE are set in XPT(NF+1,.).
!
NF=0
50 NFM=NF
NFX=NF-N
NF=NF+1
IF (NFM .LE. 2*N) THEN
IF (NFM .GE. 1 .AND. NFM .LE. N) THEN
STEPA=RHOBEG
IF (SU(NFM) .EQ. ZERO) STEPA=-STEPA
XPT(NF,NFM)=STEPA
ELSE IF (NFM .GT. N) THEN
STEPA=XPT(NF-N,NFX)
STEPB=-RHOBEG
IF (SL(NFX) .EQ. ZERO) STEPB=DMIN1(TWO*RHOBEG,SU(NFX))
IF (SU(NFX) .EQ. ZERO) STEPB=DMAX1(-TWO*RHOBEG,SL(NFX))
XPT(NF,NFX)=STEPB
END IF
ELSE
ITEMP=(NFM-NP)/N
JPT=NFM-ITEMP*N-N
IPT=JPT+ITEMP
IF (IPT .GT. N) THEN
ITEMP=JPT
JPT=IPT-N
IPT=ITEMP
END IF
XPT(NF,IPT)=XPT(IPT+1,IPT)
XPT(NF,JPT)=XPT(JPT+1,JPT)
END IF
!
!     Calculate the next value of F. The least function value so far and
!     its index are required.
!
DO J=1,N
X(J)=DMIN1(DMAX1(XL(J),XBASE(J)+XPT(NF,J)),XU(J))
IF (XPT(NF,J) .EQ. SL(J)) X(J)=XL(J)
IF (XPT(NF,J) .EQ. SU(J)) X(J)=XU(J)
END DO
F = CALFUN (N,X,IPRINT)
!$$$      IF (IPRINT .EQ. 3) THEN
!$$$          PRINT 70, NF,F,(X(I),I=1,N)
!$$$   70      FORMAT (/4X,'Function number',I6,'    F =',1PD18.10,
!$$$     1       '    The corresponding X is:'/(2X,5D15.6))
!$$$      END IF
!$$$      CALL minqi3 (IPRINT, F, NF, N, X)
FVAL(NF)=F
IF (NF .EQ. 1) THEN
FBEG=F
KOPT=1
ELSE IF (F .LT. FVAL(KOPT)) THEN
KOPT=NF
END IF
!
!     Set the nonzero initial elements of BMAT and the quadratic model in the
!     cases when NF is at most 2*N+1. If NF exceeds N+1, then the positions
!     of the NF-th and (NF-N)-th interpolation points may be switched, in
!     order that the function value at the first of them contributes to the
!     off-diagonal second derivative terms of the initial quadratic model.
!
IF (NF .LE. 2*N+1) THEN
IF (NF .GE. 2 .AND. NF .LE. N+1) THEN
GOPT(NFM)=(F-FBEG)/STEPA
IF (NPT .LT. NF+N) THEN
BMAT(1,NFM)=-ONE/STEPA
BMAT(NF,NFM)=ONE/STEPA
BMAT(NPT+NFM,NFM)=-HALF*RHOSQ
END IF
ELSE IF (NF .GE. N+2) THEN
IH=(NFX*(NFX+1))/2
TEMP=(F-FBEG)/STEPB
DIFF=STEPB-STEPA
HQ(IH)=TWO*(TEMP-GOPT(NFX))/DIFF
GOPT(NFX)=(GOPT(NFX)*STEPB-TEMP*STEPA)/DIFF
IF (STEPA*STEPB .LT. ZERO) THEN
IF (F .LT. FVAL(NF-N)) THEN
FVAL(NF)=FVAL(NF-N)
FVAL(NF-N)=F
IF (KOPT .EQ. NF) KOPT=NF-N
XPT(NF-N,NFX)=STEPB
XPT(NF,NFX)=STEPA
END IF
END IF
BMAT(1,NFX)=-(STEPA+STEPB)/(STEPA*STEPB)
BMAT(NF,NFX)=-HALF/XPT(NF-N,NFX)
BMAT(NF-N,NFX)=-BMAT(1,NFX)-BMAT(NF,NFX)
ZMAT(1,NFX)=DSQRT(TWO)/(STEPA*STEPB)
ZMAT(NF,NFX)=DSQRT(HALF)/RHOSQ
ZMAT(NF-N,NFX)=-ZMAT(1,NFX)-ZMAT(NF,NFX)
END IF
!
!     Set the off-diagonal second derivatives of the Lagrange functions and
!     the initial quadratic model.
!
ELSE
IH=(IPT*(IPT-1))/2+JPT
ZMAT(1,NFX)=RECIP
ZMAT(NF,NFX)=RECIP
ZMAT(IPT+1,NFX)=-RECIP
ZMAT(JPT+1,NFX)=-RECIP
TEMP=XPT(NF,IPT)*XPT(NF,JPT)
HQ(IH)=(FBEG-FVAL(IPT+1)-FVAL(JPT+1)+F)/TEMP
END IF
IF (NF .LT. NPT .AND. NF .LT. MAXFUN) GOTO 50
RETURN
END


SUBROUTINE RESCUE (N,NPT,XL,XU,IPRINT,MAXFUN,XBASE,XPT, FVAL,XOPT,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,SL,SU,NF, &
    & DELTA, KOPT,VLAG,PTSAUX,PTSID,W)
implicit integer (i-n), real(real64) (a-h,o-z)
DIMENSION XL(*),XU(*),XBASE(*),XPT(NPT,*),FVAL(*),XOPT(*), GOPT(*),HQ(*),PQ(*),BMAT(NDIM,*),ZMAT(NPT,*), &
    & SL(*),SU(*), VLAG(*),PTSAUX(2,*),PTSID(*),W(*)
!
!     The arguments N, NPT, XL, XU, IPRINT, MAXFUN, XBASE, XPT, FVAL, XOPT,
!       GOPT, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU have the same meanings as
!       the corresponding arguments of BOBYQB on the entry to RESCUE.
!     NF is maintained as the number of calls of CALFUN so far, except that
!       NF is set to -1 if the value of MAXFUN prevents further progress.
!     KOPT is maintained so that FVAL(KOPT) is the least calculated function
!       value. Its correct value must be given on entry. It is updated if a
!       new least function value is found, but the corresponding changes to
!       XOPT and GOPT have to be made later by the calling program.
!     DELTA is the current trust region radius.
!     VLAG is a working space vector that will be used for the values of the
!       provisional Lagrange functions at each of the interpolation points.
!       They are part of a product that requires VLAG to be of length NDIM.
!     PTSAUX is also a working space array. For J=1,2,...,N, PTSAUX(1,J) and
!       PTSAUX(2,J) specify the two positions of provisional interpolation
!       points when a nonzero step is taken along e_J (the J-th coordinate
!       direction) through XBASE+XOPT, as specified below. Usually these
!       steps have length DELTA, but other lengths are chosen if necessary
!       in order to satisfy the given bounds on the variables.
!     PTSID is also a working space array. It has NPT components that denote
!       provisional new positions of the original interpolation points, in
!       case changes are needed to restore the linear independence of the
!       interpolation conditions. The K-th point is a candidate for change
!       if and only if PTSID(K) is nonzero. In this case let p and q be the
!       integer parts of PTSID(K) and (PTSID(K)-p) multiplied by N+1. If p
!       and q are both positive, the step from XBASE+XOPT to the new K-th
!       interpolation point is PTSAUX(1,p)*e_p + PTSAUX(1,q)*e_q. Otherwise
!       the step is PTSAUX(1,p)*e_p or PTSAUX(2,q)*e_q in the cases q=0 or
!       p=0, respectively.
!     The first NDIM+NPT elements of the array W are used for working space.
!     The final elements of BMAT and ZMAT are set in a well-conditioned way
!       to the values that are appropriate for the new interpolation points.
!     The elements of GOPT, HQ and PQ are also revised to the values that are
!       appropriate to the final quadratic model.
!
!     Set some constants.
!
HALF=0.5D0
ONE=1.0D0
ZERO=0.0D0
NP=N+1
SFRAC=HALF/DBLE(NP)
NPTM=NPT-NP
!
!     Shift the interpolation points so that XOPT becomes the origin, and set
!     the elements of ZMAT to zero. The value of SUMPQ is required in the
!     updating of HQ below. The squares of the distances from XOPT to the
!     other interpolation points are set at the end of W. Increments of WINC
!     may be added later to these squares to balance the consideration of
!     the choice of point that is going to become current.
!
SUMPQ=ZERO
WINC=ZERO
DO K=1,NPT
DISTSQ=ZERO
DO J=1,N
XPT(K,J)=XPT(K,J)-XOPT(J)
DISTSQ=DISTSQ+XPT(K,J)**2
END DO
SUMPQ=SUMPQ+PQ(K)
W(NDIM+K)=DISTSQ
WINC=DMAX1(WINC,DISTSQ)
DO J=1,NPTM
ZMAT(K,J)=ZERO
END DO
END DO
!
!     Update HQ so that HQ and PQ define the second derivatives of the model
!     after XBASE has been shifted to the trust region centre.
!
IH=0
DO J=1,N
W(J)=HALF*SUMPQ*XOPT(J)
DO K=1,NPT
W(J)=W(J)+PQ(K)*XPT(K,J)
END DO
DO I=1,J
IH=IH+1
HQ(IH)=HQ(IH)+W(I)*XOPT(J)+W(J)*XOPT(I)
END DO
END DO
!
!     Shift XBASE, SL, SU and XOPT. Set the elements of BMAT to zero, and
!     also set the elements of PTSAUX.
!
DO J=1,N
XBASE(J)=XBASE(J)+XOPT(J)
SL(J)=SL(J)-XOPT(J)
SU(J)=SU(J)-XOPT(J)
XOPT(J)=ZERO
PTSAUX(1,J)=DMIN1(DELTA,SU(J))
PTSAUX(2,J)=DMAX1(-DELTA,SL(J))
IF (PTSAUX(1,J)+PTSAUX(2,J) .LT. ZERO) THEN
TEMP=PTSAUX(1,J)
PTSAUX(1,J)=PTSAUX(2,J)
PTSAUX(2,J)=TEMP
END IF
IF (DABS(PTSAUX(2,J)) .LT. HALF*DABS(PTSAUX(1,J))) THEN
PTSAUX(2,J)=HALF*PTSAUX(1,J)
END IF
DO I=1,NDIM
BMAT(I,J)=ZERO
END DO
END DO
FBASE=FVAL(KOPT)
!
!     Set the identifiers of the artificial interpolation points that are
!     along a coordinate direction from XOPT, and set the corresponding
!     nonzero elements of BMAT and ZMAT.
!
PTSID(1)=SFRAC
DO J=1,N
JP=J+1
JPN=JP+N
PTSID(JP)=DBLE(J)+SFRAC
IF (JPN .LE. NPT) THEN
PTSID(JPN)=DBLE(J)/DBLE(NP)+SFRAC
TEMP=ONE/(PTSAUX(1,J)-PTSAUX(2,J))
BMAT(JP,J)=-TEMP+ONE/PTSAUX(1,J)
BMAT(JPN,J)=TEMP+ONE/PTSAUX(2,J)
BMAT(1,J)=-BMAT(JP,J)-BMAT(JPN,J)
ZMAT(1,J)=DSQRT(2.0D0)/DABS(PTSAUX(1,J)*PTSAUX(2,J))
ZMAT(JP,J)=ZMAT(1,J)*PTSAUX(2,J)*TEMP
ZMAT(JPN,J)=-ZMAT(1,J)*PTSAUX(1,J)*TEMP
ELSE
BMAT(1,J)=-ONE/PTSAUX(1,J)
BMAT(JP,J)=ONE/PTSAUX(1,J)
BMAT(J+NPT,J)=-HALF*PTSAUX(1,J)**2
END IF
END DO
!
!     Set any remaining identifiers with their nonzero elements of ZMAT.
!
IF (NPT .GE. N+NP) THEN
DO K=2*NP,NPT
IW=(DBLE(K-NP)-HALF)/DBLE(N)
IP=K-NP-IW*N
IQ=IP+IW
IF (IQ .GT. N) IQ=IQ-N
PTSID(K)=DBLE(IP)+DBLE(IQ)/DBLE(NP)+SFRAC
TEMP=ONE/(PTSAUX(1,IP)*PTSAUX(1,IQ))
ZMAT(1,K-NP)=TEMP
ZMAT(IP+1,K-NP)=-TEMP
ZMAT(IQ+1,K-NP)=-TEMP
ZMAT(K,K-NP)=TEMP
END DO
END IF
NREM=NPT
KOLD=1
KNEW=KOPT
!
!     Reorder the provisional points in the way that exchanges PTSID(KOLD)
!     with PTSID(KNEW).
!
80 DO J=1,N
TEMP=BMAT(KOLD,J)
BMAT(KOLD,J)=BMAT(KNEW,J)
BMAT(KNEW,J)=TEMP
END DO
DO J=1,NPTM
TEMP=ZMAT(KOLD,J)
ZMAT(KOLD,J)=ZMAT(KNEW,J)
ZMAT(KNEW,J)=TEMP
END DO
PTSID(KOLD)=PTSID(KNEW)
PTSID(KNEW)=ZERO
W(NDIM+KNEW)=ZERO
NREM=NREM-1
IF (KNEW .NE. KOPT) THEN
TEMP=VLAG(KOLD)
VLAG(KOLD)=VLAG(KNEW)
VLAG(KNEW)=TEMP
!
!     Update the BMAT and ZMAT matrices so that the status of the KNEW-th
!     interpolation point can be changed from provisional to original. The
!     branch to label 350 occurs if all the original points are reinstated.
!     The nonnegative values of W(NDIM+K) are required in the search below.
!
CALL UPDATEBOBYQA (N,NPT,BMAT,ZMAT,NDIM,VLAG,BETA,DENOM,KNEW, W)
IF (NREM .EQ. 0) GOTO 350
DO K=1,NPT
W(NDIM+K)=DABS(W(NDIM+K))
END DO
END IF
!
!     Pick the index KNEW of an original interpolation point that has not
!     yet replaced one of the provisional interpolation points, giving
!     attention to the closeness to XOPT and to previous tries with KNEW.
!
120 DSQMIN=ZERO
DO K=1,NPT
IF (W(NDIM+K) .GT. ZERO) THEN
IF (DSQMIN .EQ. ZERO .OR. W(NDIM+K) .LT. DSQMIN) THEN
KNEW=K
DSQMIN=W(NDIM+K)
END IF
END IF
END DO
IF (DSQMIN .EQ. ZERO) GOTO 260
!
!     Form the W-vector of the chosen original interpolation point.
!
DO J=1,N
W(NPT+J)=XPT(KNEW,J)
END DO
DO K=1,NPT
SUM=ZERO
IF (K .EQ. KOPT) THEN
CONTINUE
ELSE IF (PTSID(K) .EQ. ZERO) THEN
DO J=1,N
SUM=SUM+W(NPT+J)*XPT(K,J)
END DO
ELSE
IP=PTSID(K)
IF (IP .GT. 0) SUM=W(NPT+IP)*PTSAUX(1,IP)
IQ=DBLE(NP)*PTSID(K)-DBLE(IP*NP)
IF (IQ .GT. 0) THEN
IW=1
IF (IP .EQ. 0) IW=2
SUM=SUM+W(NPT+IQ)*PTSAUX(IW,IQ)
END IF
END IF
W(K)=HALF*SUM*SUM
END DO
!
!     Calculate VLAG and BETA for the required updating of the H matrix if
!     XPT(KNEW,.) is reinstated in the set of interpolation points.
!
DO K=1,NPT
SUM=ZERO
DO J=1,N
SUM=SUM+BMAT(K,J)*W(NPT+J)
END DO
VLAG(K)=SUM
END DO
BETA=ZERO
DO J=1,NPTM
SUM=ZERO
DO K=1,NPT
SUM=SUM+ZMAT(K,J)*W(K)
END DO
BETA=BETA-SUM*SUM
DO K=1,NPT
VLAG(K)=VLAG(K)+SUM*ZMAT(K,J)
END DO
END DO
BSUM=ZERO
DISTSQ=ZERO
DO J=1,N
SUM=ZERO
DO K=1,NPT
SUM=SUM+BMAT(K,J)*W(K)
END DO
JP=J+NPT
BSUM=BSUM+SUM*W(JP)
DO IP=NPT+1,NDIM
SUM=SUM+BMAT(IP,J)*W(IP)
END DO
BSUM=BSUM+SUM*W(JP)
VLAG(JP)=SUM
DISTSQ=DISTSQ+XPT(KNEW,J)**2
END DO
BETA=HALF*DISTSQ*DISTSQ+BETA-BSUM
VLAG(KOPT)=VLAG(KOPT)+ONE
!
!     KOLD is set to the index of the provisional interpolation point that is
!     going to be deleted to make way for the KNEW-th original interpolation
!     point. The choice of KOLD is governed by the avoidance of a small value
!     of the denominator in the updating calculation of UPDATE.
!
DENOM=ZERO
VLMXSQ=ZERO
DO K=1,NPT
IF (PTSID(K) .NE. ZERO) THEN
HDIAG=ZERO
DO J=1,NPTM
HDIAG=HDIAG+ZMAT(K,J)**2
END DO
DEN=BETA*HDIAG+VLAG(K)**2
IF (DEN .GT. DENOM) THEN
KOLD=K
DENOM=DEN
END IF
END IF
VLMXSQ=DMAX1(VLMXSQ,VLAG(K)**2)
END DO
IF (DENOM .LE. 1.0D-2*VLMXSQ) THEN
W(NDIM+KNEW)=-W(NDIM+KNEW)-WINC
GOTO 120
END IF
GOTO 80
!
!     When label 260 is reached, all the final positions of the interpolation
!     points have been chosen although any changes have not been included yet
!     in XPT. Also the final BMAT and ZMAT matrices are complete, but, apart
!     from the shift of XBASE, the updating of the quadratic model remains to
!     be done. The following cycle through the new interpolation points begins
!     by putting the new point in XPT(KPT,.) and by setting PQ(KPT) to zero,
!     except that a RETURN occurs if MAXFUN prohibits another value of F.
!
260 DO KPT=1,NPT
IF (PTSID(KPT) .EQ. ZERO) GOTO 340
IF (NF .GE. MAXFUN) THEN
NF=-1
GOTO 350
END IF
IH=0
DO J=1,N
W(J)=XPT(KPT,J)
XPT(KPT,J)=ZERO
TEMP=PQ(KPT)*W(J)
DO I=1,J
IH=IH+1
HQ(IH)=HQ(IH)+TEMP*W(I)
END DO
END DO
PQ(KPT)=ZERO
IP=PTSID(KPT)
IQ=DBLE(NP)*PTSID(KPT)-DBLE(IP*NP)
IF (IP .GT. 0) THEN
XP=PTSAUX(1,IP)
XPT(KPT,IP)=XP
END IF
IF (IQ .GT. 0) THEN
XQ=PTSAUX(1,IQ)
IF (IP .EQ. 0) XQ=PTSAUX(2,IQ)
XPT(KPT,IQ)=XQ
END IF
!
!     Set VQUAD to the value of the current model at the new point.
!
VQUAD=FBASE
IF (IP .GT. 0) THEN
IHP=(IP+IP*IP)/2
VQUAD=VQUAD+XP*(GOPT(IP)+HALF*XP*HQ(IHP))
END IF
IF (IQ .GT. 0) THEN
IHQ=(IQ+IQ*IQ)/2
VQUAD=VQUAD+XQ*(GOPT(IQ)+HALF*XQ*HQ(IHQ))
IF (IP .GT. 0) THEN
IW=MAX0(IHP,IHQ)-IABS(IP-IQ)
VQUAD=VQUAD+XP*XQ*HQ(IW)
END IF
END IF
DO K=1,NPT
TEMP=ZERO
IF (IP .GT. 0) TEMP=TEMP+XP*XPT(K,IP)
IF (IQ .GT. 0) TEMP=TEMP+XQ*XPT(K,IQ)
VQUAD=VQUAD+HALF*PQ(K)*TEMP*TEMP
END DO
!
!     Calculate F at the new interpolation point, and set DIFF to the factor
!     that is going to multiply the KPT-th Lagrange function when the model
!     is updated to provide interpolation to the new function value.
!
DO I=1,N
W(I)=DMIN1(DMAX1(XL(I),XBASE(I)+XPT(KPT,I)),XU(I))
IF (XPT(KPT,I) .EQ. SL(I)) W(I)=XL(I)
IF (XPT(KPT,I) .EQ. SU(I)) W(I)=XU(I)
END DO
NF=NF+1
F = CALFUN (N,W,IPRINT)
!$$$  IF (IPRINT .EQ. 3) THEN
!$$$  PRINT 70, NF,F,(X(I),I=1,N)
!$$$   70      FORMAT (/4X,'Function number',I6,'    F =',1PD18.10,
!$$$     1       '    The corresponding X is:'/(2X,5D15.6))
!$$$      END IF
!$$$      CALL minqi3 (IPRINT, F, NF, N, X)
FVAL(KPT)=F
IF (F .LT. FVAL(KOPT)) KOPT=KPT
DIFF=F-VQUAD
!
!     Update the quadratic model. The RETURN from the subroutine occurs when
!     all the new interpolation points are included in the model.
!
DO I=1,N
GOPT(I)=GOPT(I)+DIFF*BMAT(KPT,I)
END DO
DO K=1,NPT
SUM=ZERO
DO J=1,NPTM
SUM=SUM+ZMAT(K,J)*ZMAT(KPT,J)
END DO
TEMP=DIFF*SUM
IF (PTSID(K) .EQ. ZERO) THEN
PQ(K)=PQ(K)+TEMP
ELSE
IP=PTSID(K)
IQ=DBLE(NP)*PTSID(K)-DBLE(IP*NP)
IHQ=(IQ*IQ+IQ)/2
IF (IP .EQ. 0) THEN
HQ(IHQ)=HQ(IHQ)+TEMP*PTSAUX(2,IQ)**2
ELSE
IHP=(IP*IP+IP)/2
HQ(IHP)=HQ(IHP)+TEMP*PTSAUX(1,IP)**2
IF (IQ .GT. 0) THEN
HQ(IHQ)=HQ(IHQ)+TEMP*PTSAUX(1,IQ)**2
IW=MAX0(IHP,IHQ)-IABS(IQ-IP)
HQ(IW)=HQ(IW)+TEMP*PTSAUX(1,IP)*PTSAUX(1,IQ)
END IF
END IF
END IF
END DO
PTSID(KPT)=ZERO
340 END DO
350 RETURN
END


SUBROUTINE TRSAPP (N,NPT,XOPT,XPT,GQ,HQ,PQ,DELTA,STEP, D,G,HD,HS,CRVMIN)
implicit integer (i-n), real(real64) (a-h,o-z)
DIMENSION XOPT(*),XPT(NPT,*),GQ(*),HQ(*),PQ(*),STEP(*), D(*),G(*),HD(*),HS(*)
!
!     N is the number of variables of a quadratic objective function, Q say.
!     The arguments NPT, XOPT, XPT, GQ, HQ and PQ have their usual meanings,
!       in order to define the current quadratic model Q.
!     DELTA is the trust region radius, and has to be positive.
!     STEP will be set to the calculated trial step.
!     The arrays D, G, HD and HS will be used for working space.
!     CRVMIN will be set to the least curvature of H along the conjugate
!       directions that occur, except that it is set to zero if STEP goes
!       all the way to the trust region boundary.
!
!     The calculation of STEP begins with the truncated conjugate gradient
!     method. If the boundary of the trust region is reached, then further
!     changes to STEP may be made, each one being in the 2D space spanned
!     by the current STEP and the corresponding gradient of Q. Thus STEP
!     should provide a substantial reduction to Q within the trust region.
!
!     Initialization, which includes setting HD to H times XOPT.
!
HALF=0.5D0
ZERO=0.0D0
TWOPI=8.0D0*DATAN(1.0D0)
DELSQ=DELTA*DELTA
ITERC=0
ITERMAX=N
ITERSW=ITERMAX
DO I=1,N
D(I)=XOPT(I)
END DO
GOTO 170
!
!     Prepare for the first line search.
!
20 QRED=ZERO
DD=ZERO
DO I=1,N
STEP(I)=ZERO
HS(I)=ZERO
G(I)=GQ(I)+HD(I)
D(I)=-G(I)
DD=DD+D(I)**2
END DO
CRVMIN=ZERO
IF (DD .EQ. ZERO) GOTO 160
DS=ZERO
SS=ZERO
GG=DD
GGBEG=GG
!
!     Calculate the step to the trust region boundary and the product HD.
!
40 ITERC=ITERC+1
TEMP=DELSQ-SS
BSTEP=TEMP/(DS+DSQRT(DS*DS+DD*TEMP))
GOTO 170
50 DHD=ZERO
DO J=1,N
DHD=DHD+D(J)*HD(J)
END DO
!
!     Update CRVMIN and set the step-length ALPHA.
!
ALPHA=BSTEP
IF (DHD .GT. ZERO) THEN
TEMP=DHD/DD
IF (ITERC .EQ. 1) CRVMIN=TEMP
CRVMIN=DMIN1(CRVMIN,TEMP)
ALPHA=DMIN1(ALPHA,GG/DHD)
END IF
QADD=ALPHA*(GG-HALF*ALPHA*DHD)
QRED=QRED+QADD
!
!     Update STEP and HS.
!
GGSAV=GG
GG=ZERO
DO I=1,N
STEP(I)=STEP(I)+ALPHA*D(I)
HS(I)=HS(I)+ALPHA*HD(I)
GG=GG+(G(I)+HS(I))**2
END DO
!
!     Begin another conjugate direction iteration if required.
!
IF (ALPHA .LT. BSTEP) THEN
IF (QADD .LE. 0.01D0*QRED) GOTO 160
IF (GG .LE. 1.0D-4*GGBEG) GOTO 160
IF (ITERC .EQ. ITERMAX) GOTO 160
TEMP=GG/GGSAV
DD=ZERO
DS=ZERO
SS=ZERO
DO I=1,N
D(I)=TEMP*D(I)-G(I)-HS(I)
DD=DD+D(I)**2
DS=DS+D(I)*STEP(I)
SS=SS+STEP(I)**2
END DO
IF (DS .LE. ZERO) GOTO 160
IF (SS .LT. DELSQ) GOTO 40
END IF
CRVMIN=ZERO
ITERSW=ITERC
!
!     Test whether an alternative iteration is required.
!
90 IF (GG .LE. 1.0D-4*GGBEG) GOTO 160
SG=ZERO
SHS=ZERO
DO I=1,N
SG=SG+STEP(I)*G(I)
SHS=SHS+STEP(I)*HS(I)
END DO
SGK=SG+SHS
ANGTEST=SGK/DSQRT(GG*DELSQ)
IF (ANGTEST .LE. -0.99D0) GOTO 160
!
!     Begin the alternative iteration by calculating D and HD and some
!     scalar products.
!
ITERC=ITERC+1
TEMP=DSQRT(DELSQ*GG-SGK*SGK)
TEMPA=DELSQ/TEMP
TEMPB=SGK/TEMP
DO I=1,N
D(I)=TEMPA*(G(I)+HS(I))-TEMPB*STEP(I)
END DO
GOTO 170
120 DG=ZERO
DHD=ZERO
DHS=ZERO
DO I=1,N
DG=DG+D(I)*G(I)
DHD=DHD+HD(I)*D(I)
DHS=DHS+HD(I)*STEP(I)
END DO
!
!     Seek the value of the angle that minimizes Q.
!
CF=HALF*(SHS-DHD)
QBEG=SG+CF
QSAV=QBEG
QMIN=QBEG
ISAVE=0
IU=49
TEMP=TWOPI/DBLE(IU+1)
DO I=1,IU
ANGLE=DBLE(I)*TEMP
CTH=DCOS(ANGLE)
STH=DSIN(ANGLE)
QNEW=(SG+CF*CTH)*CTH+(DG+DHS*CTH)*STH
IF (QNEW .LT. QMIN) THEN
QMIN=QNEW
ISAVE=I
TEMPA=QSAV
ELSE IF (I .EQ. ISAVE+1) THEN
TEMPB=QNEW
END IF
QSAV=QNEW
END DO
IF (ISAVE .EQ. ZERO) TEMPA=QNEW
IF (ISAVE .EQ. IU) TEMPB=QBEG
ANGLE=ZERO
IF (TEMPA .NE. TEMPB) THEN
TEMPA=TEMPA-QMIN
TEMPB=TEMPB-QMIN
ANGLE=HALF*(TEMPA-TEMPB)/(TEMPA+TEMPB)
END IF
ANGLE=TEMP*(DBLE(ISAVE)+ANGLE)
!
!     Calculate the new STEP and HS. Then test for convergence.
!
CTH=DCOS(ANGLE)
STH=DSIN(ANGLE)
REDUC=QBEG-(SG+CF*CTH)*CTH-(DG+DHS*CTH)*STH
GG=ZERO
DO I=1,N
STEP(I)=CTH*STEP(I)+STH*D(I)
HS(I)=CTH*HS(I)+STH*HD(I)
GG=GG+(G(I)+HS(I))**2
END DO
QRED=QRED+REDUC
RATIO=REDUC/QRED
IF (ITERC .LT. ITERMAX .AND. RATIO .GT. 0.01D0) GOTO 90
160 RETURN
!
!     The following instructions act as a subroutine for setting the vector
!     HD to the vector D multiplied by the second derivative matrix of Q.
!     They are called from three different places, which are distinguished
!     by the value of ITERC.
!
170 DO I=1,N
HD(I)=ZERO
END DO
DO K=1,NPT
TEMP=ZERO
DO J=1,N
TEMP=TEMP+XPT(K,J)*D(J)
END DO
TEMP=TEMP*PQ(K)
DO I=1,N
HD(I)=HD(I)+TEMP*XPT(K,I)
END DO
END DO
IH=0
DO J=1,N
DO I=1,J
IH=IH+1
IF (I .LT. J) HD(J)=HD(J)+HQ(IH)*D(I)
HD(I)=HD(I)+HQ(IH)*D(J)
END DO
END DO
IF (ITERC .EQ. 0) GOTO 20
IF (ITERC .LE. ITERSW) GOTO 50
GOTO 120
END


SUBROUTINE TRSBOX (N,NPT,XPT,XOPT,GOPT,HQ,PQ,SL,SU,DELTA, XNEW,D,GNEW,XBDI,S,HS,HRED,DSQ,CRVMIN)
implicit integer (i-n), real(real64) (a-h,o-z)
DIMENSION XPT(NPT,*),XOPT(*),GOPT(*),HQ(*),PQ(*),SL(*),SU(*), XNEW(*),D(*),GNEW(*),XBDI(*),S(*),HS(*), &
    & HRED(*)
!
!     The arguments N, NPT, XPT, XOPT, GOPT, HQ, PQ, SL and SU have the same
!       meanings as the corresponding arguments of BOBYQB.
!     DELTA is the trust region radius for the present calculation, which
!       seeks a small value of the quadratic model within distance DELTA of
!       XOPT subject to the bounds on the variables.
!     XNEW will be set to a new vector of variables that is approximately
!       the one that minimizes the quadratic model within the trust region
!       subject to the SL and SU constraints on the variables. It satisfies
!       as equations the bounds that become active during the calculation.
!     D is the calculated trial step from XOPT, generated iteratively from an
!       initial value of zero. Thus XNEW is XOPT+D after the final iteration.
!     GNEW holds the gradient of the quadratic model at XOPT+D. It is updated
!       when D is updated.
!     XBDI is a working space vector. For I=1,2,...,N, the element XBDI(I) is
!       set to -1.0, 0.0, or 1.0, the value being nonzero if and only if the
!       I-th variable has become fixed at a bound, the bound being SL(I) or
!       SU(I) in the case XBDI(I)=-1.0 or XBDI(I)=1.0, respectively. This
!       information is accumulated during the construction of XNEW.
!     The arrays S, HS and HRED are also used for working space. They hold the
!       current search direction, and the changes in the gradient of Q along S
!       and the reduced D, respectively, where the reduced D is the same as D,
!       except that the components of the fixed variables are zero.
!     DSQ will be set to the square of the length of XNEW-XOPT.
!     CRVMIN is set to zero if D reaches the trust region boundary. Otherwise
!       it is set to the least curvature of H that occurs in the conjugate
!       gradient searches that are not restricted by any constraints. The
!       value CRVMIN=-1.0D0 is set, however, if all of these searches are
!       constrained.
!
!     A version of the truncated conjugate gradient is applied. If a line
!     search is restricted by a constraint, then the procedure is restarted,
!     the values of the variables that are at their bounds being fixed. If
!     the trust region boundary is reached, then further changes may be made
!     to D, each one being in the two dimensional space that is spanned
!     by the current D and the gradient of Q at XOPT+D, staying on the trust
!     region boundary. Termination occurs when the reduction in Q seems to
!     be close to the greatest reduction that can be achieved.
!
!     Set some constants.
!
HALF=0.5D0
ONE=1.0D0
ONEMIN=-1.0D0
ZERO=0.0D0
!
!     The sign of GOPT(I) gives the sign of the change to the I-th variable
!     that will reduce Q from its value at XOPT. Thus XBDI(I) shows whether
!     or not to fix the I-th variable at one of its bounds initially, with
!     NACT being set to the number of fixed variables. D and GNEW are also
!     set for the first iteration. DELSQ is the upper bound on the sum of
!     squares of the free variables. QRED is the reduction in Q so far.
!
ITERC=0
NACT=0
SQSTP=ZERO
DO I=1,N
XBDI(I)=ZERO
IF (XOPT(I) .LE. SL(I)) THEN
IF (GOPT(I) .GE. ZERO) XBDI(I)=ONEMIN
ELSE IF (XOPT(I) .GE. SU(I)) THEN
IF (GOPT(I) .LE. ZERO) XBDI(I)=ONE
END IF
IF (XBDI(I) .NE. ZERO) NACT=NACT+1
D(I)=ZERO
GNEW(I)=GOPT(I)
END DO
DELSQ=DELTA*DELTA
QRED=ZERO
CRVMIN=ONEMIN
!
!     Set the next search direction of the conjugate gradient method. It is
!     the steepest descent direction initially and when the iterations are
!     restarted because a variable has just been fixed by a bound, and of
!     course the components of the fixed variables are zero. ITERMAX is an
!     upper bound on the indices of the conjugate gradient iterations.
!
20 BETA=ZERO
30 STEPSQ=ZERO
DO I=1,N
IF (XBDI(I) .NE. ZERO) THEN
S(I)=ZERO
ELSE IF (BETA .EQ. ZERO) THEN
S(I)=-GNEW(I)
ELSE
S(I)=BETA*S(I)-GNEW(I)
END IF
STEPSQ=STEPSQ+S(I)**2
END DO
IF (STEPSQ .EQ. ZERO) GOTO 190
IF (BETA .EQ. ZERO) THEN
GREDSQ=STEPSQ
ITERMAX=ITERC+N-NACT
END IF
IF (GREDSQ*DELSQ .LE. 1.0D-4*QRED*QRED) GO TO 190
!
!     Multiply the search direction by the second derivative matrix of Q and
!     calculate some scalars for the choice of steplength. Then set BLEN to
!     the length of the the step to the trust region boundary and STPLEN to
!     the steplength, ignoring the simple bounds.
!
GOTO 210
50 RESID=DELSQ
DS=ZERO
SHS=ZERO
DO I=1,N
IF (XBDI(I) .EQ. ZERO) THEN
RESID=RESID-D(I)**2
DS=DS+S(I)*D(I)
SHS=SHS+S(I)*HS(I)
END IF
END DO
IF (RESID .LE. ZERO) GOTO 90
TEMP=DSQRT(STEPSQ*RESID+DS*DS)
IF (DS .LT. ZERO) THEN
BLEN=(TEMP-DS)/STEPSQ
ELSE
BLEN=RESID/(TEMP+DS)
END IF
STPLEN=BLEN
IF (SHS .GT. ZERO) THEN
STPLEN=DMIN1(BLEN,GREDSQ/SHS)
END IF

!
!     Reduce STPLEN if necessary in order to preserve the simple bounds,
!     letting IACT be the index of the new constrained variable.
!
IACT=0
DO I=1,N
IF (S(I) .NE. ZERO) THEN
XSUM=XOPT(I)+D(I)
IF (S(I) .GT. ZERO) THEN
TEMP=(SU(I)-XSUM)/S(I)
ELSE
TEMP=(SL(I)-XSUM)/S(I)
END IF
IF (TEMP .LT. STPLEN) THEN
STPLEN=TEMP
IACT=I
END IF
END IF
END DO
!
!     Update CRVMIN, GNEW and D. Set SDEC to the decrease that occurs in Q.
!
SDEC=ZERO
IF (STPLEN .GT. ZERO) THEN
ITERC=ITERC+1
TEMP=SHS/STEPSQ
IF (IACT .EQ. 0 .AND. TEMP .GT. ZERO) THEN
CRVMIN=DMIN1(CRVMIN,TEMP)
IF (CRVMIN .EQ. ONEMIN) CRVMIN=TEMP
END IF
GGSAV=GREDSQ
GREDSQ=ZERO
DO I=1,N
GNEW(I)=GNEW(I)+STPLEN*HS(I)
IF (XBDI(I) .EQ. ZERO) GREDSQ=GREDSQ+GNEW(I)**2
D(I)=D(I)+STPLEN*S(I)
END DO
SDEC=DMAX1(STPLEN*(GGSAV-HALF*STPLEN*SHS),ZERO)
QRED=QRED+SDEC
END IF
!
!     Restart the conjugate gradient method if it has hit a new bound.
!
IF (IACT .GT. 0) THEN
NACT=NACT+1
XBDI(IACT)=ONE
IF (S(IACT) .LT. ZERO) XBDI(IACT)=ONEMIN
DELSQ=DELSQ-D(IACT)**2
IF (DELSQ .LE. ZERO) GOTO 90
GOTO 20
END IF
!
!     If STPLEN is less than BLEN, then either apply another conjugate
!     gradient iteration or RETURN.
!
IF (STPLEN .LT. BLEN) THEN
IF (ITERC .EQ. ITERMAX) GOTO 190
IF (SDEC .LE. 0.01D0*QRED) GOTO 190
BETA=GREDSQ/GGSAV
GOTO 30
END IF
90 CRVMIN=ZERO
!
!     Prepare for the alternative iteration by calculating some scalars
!     and by multiplying the reduced D by the second derivative matrix of
!     Q, where S holds the reduced D in the call of GGMULT.
!
100 IF (NACT .GE. N-1) GOTO 190
DREDSQ=ZERO
DREDG=ZERO
GREDSQ=ZERO
DO I=1,N
IF (XBDI(I) .EQ. ZERO) THEN
DREDSQ=DREDSQ+D(I)**2
DREDG=DREDG+D(I)*GNEW(I)
GREDSQ=GREDSQ+GNEW(I)**2
S(I)=D(I)
ELSE
S(I)=ZERO
END IF
END DO
ITCSAV=ITERC
GOTO 210
!
!     Let the search direction S be a linear combination of the reduced D
!     and the reduced G that is orthogonal to the reduced D.
!
120 ITERC=ITERC+1
TEMP=GREDSQ*DREDSQ-DREDG*DREDG
IF (TEMP .LE. 1.0D-4*QRED*QRED) GOTO 190
TEMP=DSQRT(TEMP)
DO I=1,N
IF (XBDI(I) .EQ. ZERO) THEN
S(I)=(DREDG*D(I)-DREDSQ*GNEW(I))/TEMP
ELSE
S(I)=ZERO
END IF
END DO
SREDG=-TEMP
!
!     By considering the simple bounds on the variables, calculate an upper
!     bound on the tangent of half the angle of the alternative iteration,
!     namely ANGBD, except that, if already a free variable has reached a
!     bound, there is a branch back to label 100 after fixing that variable.
!
ANGBD=ONE
IACT=0
DO I=1,N
IF (XBDI(I) .EQ. ZERO) THEN
TEMPA=XOPT(I)+D(I)-SL(I)
TEMPB=SU(I)-XOPT(I)-D(I)
IF (TEMPA .LE. ZERO) THEN
NACT=NACT+1
XBDI(I)=ONEMIN
GOTO 100
ELSE IF (TEMPB .LE. ZERO) THEN
NACT=NACT+1
XBDI(I)=ONE
GOTO 100
END IF
RATIO=ONE
SSQ=D(I)**2+S(I)**2
TEMP=SSQ-(XOPT(I)-SL(I))**2
IF (TEMP .GT. ZERO) THEN
TEMP=DSQRT(TEMP)-S(I)
IF (ANGBD*TEMP .GT. TEMPA) THEN
ANGBD=TEMPA/TEMP
IACT=I
XSAV=ONEMIN
END IF
END IF
TEMP=SSQ-(SU(I)-XOPT(I))**2
IF (TEMP .GT. ZERO) THEN
TEMP=DSQRT(TEMP)+S(I)
IF (ANGBD*TEMP .GT. TEMPB) THEN
ANGBD=TEMPB/TEMP
IACT=I
XSAV=ONE
END IF
END IF
END IF
END DO
!
!     Calculate HHD and some curvatures for the alternative iteration.
!
GOTO 210
150 SHS=ZERO
DHS=ZERO
DHD=ZERO
DO I=1,N
IF (XBDI(I) .EQ. ZERO) THEN
SHS=SHS+S(I)*HS(I)
DHS=DHS+D(I)*HS(I)
DHD=DHD+D(I)*HRED(I)
END IF
END DO
!
!     Seek the greatest reduction in Q for a range of equally spaced values
!     of ANGT in [0,ANGBD], where ANGT is the tangent of half the angle of
!     the alternative iteration.
!
REDMAX=ZERO
ISAV=0
REDSAV=ZERO
IU=17.0D0*ANGBD+3.1D0
DO I=1,IU
ANGT=ANGBD*DBLE(I)/DBLE(IU)
STH=(ANGT+ANGT)/(ONE+ANGT*ANGT)
TEMP=SHS+ANGT*(ANGT*DHD-DHS-DHS)
REDNEW=STH*(ANGT*DREDG-SREDG-HALF*STH*TEMP)
IF (REDNEW .GT. REDMAX) THEN
REDMAX=REDNEW
ISAV=I
RDPREV=REDSAV
ELSE IF (I .EQ. ISAV+1) THEN
RDNEXT=REDNEW
END IF
REDSAV=REDNEW
END DO
!
!     Return if the reduction is zero. Otherwise, set the sine and cosine
!     of the angle of the alternative iteration, and calculate SDEC.
!
IF (ISAV .EQ. 0) GOTO 190
IF (ISAV .LT. IU) THEN
TEMP=(RDNEXT-RDPREV)/(REDMAX+REDMAX-RDPREV-RDNEXT)
ANGT=ANGBD*(DBLE(ISAV)+HALF*TEMP)/DBLE(IU)
END IF
CTH=(ONE-ANGT*ANGT)/(ONE+ANGT*ANGT)
STH=(ANGT+ANGT)/(ONE+ANGT*ANGT)
TEMP=SHS+ANGT*(ANGT*DHD-DHS-DHS)
SDEC=STH*(ANGT*DREDG-SREDG-HALF*STH*TEMP)
IF (SDEC .LE. ZERO) GOTO 190
!
!     Update GNEW, D and HRED. If the angle of the alternative iteration
!     is restricted by a bound on a free variable, that variable is fixed
!     at the bound.
!
DREDG=ZERO
GREDSQ=ZERO
DO I=1,N
GNEW(I)=GNEW(I)+(CTH-ONE)*HRED(I)+STH*HS(I)
IF (XBDI(I) .EQ. ZERO) THEN
D(I)=CTH*D(I)+STH*S(I)
DREDG=DREDG+D(I)*GNEW(I)
GREDSQ=GREDSQ+GNEW(I)**2
END IF
HRED(I)=CTH*HRED(I)+STH*HS(I)
END DO
QRED=QRED+SDEC
IF (IACT .GT. 0 .AND. ISAV .EQ. IU) THEN
NACT=NACT+1
XBDI(IACT)=XSAV
GOTO 100
END IF
!
!     If SDEC is sufficiently small, then RETURN after setting XNEW to
!     XOPT+D, giving careful attention to the bounds.
!
IF (SDEC .GT. 0.01D0*QRED) GOTO 120
190 DSQ=ZERO
DO I=1,N
XNEW(I)=DMAX1(DMIN1(XOPT(I)+D(I),SU(I)),SL(I))
IF (XBDI(I) .EQ. ONEMIN) XNEW(I)=SL(I)
IF (XBDI(I) .EQ. ONE) XNEW(I)=SU(I)
D(I)=XNEW(I)-XOPT(I)
DSQ=DSQ+D(I)**2
END DO
RETURN

!     The following instructions multiply the current S-vector by the second
!     derivative matrix of the quadratic model, putting the product in HS.
!     They are reached from three different parts of the software above and
!     they can be regarded as an external subroutine.
!
210 IH=0
DO J=1,N
HS(J)=ZERO
DO I=1,J
IH=IH+1
IF (I .LT. J) HS(J)=HS(J)+HQ(IH)*S(I)
HS(I)=HS(I)+HQ(IH)*S(J)
END DO
END DO
DO K=1,NPT
IF (PQ(K) .NE. ZERO) THEN
TEMP=ZERO
DO J=1,N
TEMP=TEMP+XPT(K,J)*S(J)
END DO
TEMP=TEMP*PQ(K)
DO I=1,N
HS(I)=HS(I)+TEMP*XPT(K,I)
END DO
END IF
END DO
IF (CRVMIN .NE. ZERO) GOTO 50
IF (ITERC .GT. ITCSAV) GOTO 150
DO I=1,N
HRED(I)=HS(I)
END DO
GOTO 120
END


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% trstep.f %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE TRSTEP (N,G,H,DELTA,TOL,D,GG,TD,TN,W,PIV,Z,EVALUE)
implicit integer (i-n), real(real64) (a-h,o-z)
DIMENSION G(*),H(N,*),D(*),GG(*),TD(*),TN(*),W(*),PIV(*),Z(*)
!
!     N is the number of variables of a quadratic objective function, Q say.
!     G is the gradient of Q at the origin.
!     H is the Hessian matrix of Q. Only the upper triangular and diagonal
!       parts need be set. The lower triangular part is used to store the
!       elements of a Householder similarity transformation.
!     DELTA is the trust region radius, and has to be positive.
!     TOL is the value of a tolerance from the open interval (0,1).
!     D will be set to the calculated vector of variables.
!     The arrays GG, TD, TN, W, PIV and Z will be used for working space.
!     EVALUE will be set to the least eigenvalue of H if and only if D is a
!     Newton-Raphson step. Then EVALUE will be positive, but otherwise it
!     will be set to zero.
!
!     Let MAXRED be the maximum of Q(0)-Q(D) subject to ||D|| .LEQ. DELTA,
!     and let ACTRED be the value of Q(0)-Q(D) that is actually calculated.
!     We take the view that any D is acceptable if it has the properties
!
!             ||D|| .LEQ. DELTA  and  ACTRED .LEQ. (1-TOL)*MAXRED.
!
!     The calculation of D is done by the method of Section 2 of the paper
!     by MJDP in the 1997 Dundee Numerical Analysis Conference Proceedings,
!     after transforming H to tridiagonal form.
!
!     Initialization.
!
ONE=1.0D0
TWO=2.0D0
ZERO=0.0D0
DELSQ=DELTA*DELTA
EVALUE=ZERO
NM=N-1
DO I=1,N
D(I)=ZERO
TD(I)=H(I,I)
DO J=1,I
H(I,J)=H(J,I)
END DO
END DO
!
!     Apply Householder transformations to obtain a tridiagonal matrix that
!     is similar to H, and put the elements of the Householder vectors in
!     the lower triangular part of H. Further, TD and TN will contain the
!     diagonal and other nonzero elements of the tridiagonal matrix.
!
DO K=1,NM
KP=K+1
SUM=ZERO
IF (KP .LT. N) THEN
KPP=KP+1
DO I=KPP,N
SUM=SUM+H(I,K)**2
END DO
END IF
IF (SUM .EQ. ZERO) THEN
TN(K)=H(KP,K)
H(KP,K)=ZERO
ELSE
TEMP=H(KP,K)
TN(K)=DSIGN(DSQRT(SUM+TEMP*TEMP),TEMP)
H(KP,K)=-SUM/(TEMP+TN(K))
TEMP=DSQRT(TWO/(SUM+H(KP,K)**2))
DO I=KP,N
W(I)=TEMP*H(I,K)
H(I,K)=W(I)
Z(I)=TD(I)*W(I)
END DO
WZ=ZERO
DO J=KP,NM
JP=J+1
DO I=JP,N
Z(I)=Z(I)+H(I,J)*W(J)
Z(J)=Z(J)+H(I,J)*W(I)
END DO
WZ=WZ+W(J)*Z(J)
END DO
WZ=WZ+W(N)*Z(N)
DO J=KP,N
TD(J)=TD(J)+W(J)*(WZ*W(J)-TWO*Z(J))
IF (J .LT. N) THEN
JP=J+1
DO I=JP,N
H(I,J)=H(I,J)-W(I)*Z(J)-W(J)*(Z(I)-WZ*W(I))
END DO
END IF
END DO
END IF
END DO
!
!     Form GG by applying the similarity transformation to G.
!
GSQ=ZERO
DO I=1,N
GG(I)=G(I)
GSQ=GSQ+G(I)**2
END DO
GNORM=DSQRT(GSQ)
DO K=1,NM
KP=K+1
SUM=ZERO
DO I=KP,N
SUM=SUM+GG(I)*H(I,K)
END DO
DO I=KP,N
GG(I)=GG(I)-SUM*H(I,K)
END DO
END DO
!
!     Begin the trust region calculation with a tridiagonal matrix by
!     calculating the norm of H. Then treat the case when H is zero.
!
HNORM=DABS(TD(1))+DABS(TN(1))
TDMIN=TD(1)
TN(N)=ZERO
DO I=2,N
TEMP=DABS(TN(I-1))+DABS(TD(I))+DABS(TN(I))
HNORM=DMAX1(HNORM,TEMP)
TDMIN=DMIN1(TDMIN,TD(I))
END DO
IF (HNORM .EQ. ZERO) THEN
IF (GNORM .EQ. ZERO) GOTO 400
SCALE=DELTA/GNORM
DO I=1,N
D(I)=-SCALE*GG(I)
END DO
GOTO 370
END IF
!
!     Set the initial values of PAR and its bounds.
!
PARL=DMAX1(ZERO,-TDMIN,GNORM/DELTA-HNORM)
PARLEST=PARL
PAR=PARL
PARU=ZERO
PARUEST=ZERO
POSDEF=ZERO
ITERC=0
!
!     Calculate the pivots of the Cholesky factorization of (H+PAR*I).
!
140 ITERC=ITERC+1
KSAV=0
PIV(1)=TD(1)+PAR
K=1
150 IF (PIV(K) .GT. ZERO) THEN
PIV(K+1)=TD(K+1)+PAR-TN(K)**2/PIV(K)
ELSE
IF (PIV(K) .LT. ZERO .OR. TN(K) .NE. ZERO) GOTO 160
KSAV=K
PIV(K+1)=TD(K+1)+PAR
END IF
K=K+1
IF (K .LT. N) GOTO 150
IF (PIV(K) .LT. ZERO) GOTO 160
IF (PIV(K) .EQ. ZERO) KSAV=K
!
!     Branch if all the pivots are positive, allowing for the case when
!     G is zero.
!
IF (KSAV .EQ. 0 .AND. GSQ .GT. ZERO) GOTO 230
IF (GSQ .EQ. ZERO) THEN
IF (PAR .EQ. ZERO) GOTO 370
PARU=PAR
PARUEST=PAR
IF (KSAV .EQ. 0) GOTO 190
END IF
K=KSAV
!
!     Set D to a direction of nonpositive curvature of the given tridiagonal
!     matrix, and thus revise PARLEST.
!
160 D(K)=ONE
IF (DABS(TN(K)) .LE. DABS(PIV(K))) THEN
DSQ=ONE
DHD=PIV(K)
ELSE
TEMP=TD(K+1)+PAR
IF (TEMP .LE. DABS(PIV(K))) THEN
D(K+1)=DSIGN(ONE,-TN(K))
DHD=PIV(K)+TEMP-TWO*DABS(TN(K))
ELSE
D(K+1)=-TN(K)/TEMP
DHD=PIV(K)+TN(K)*D(K+1)
END IF
DSQ=ONE+D(K+1)**2
END IF
170 IF (K .GT. 1) THEN
K=K-1
IF (TN(K) .NE. ZERO) THEN
D(K)=-TN(K)*D(K+1)/PIV(K)
DSQ=DSQ+D(K)**2
GOTO 170
END IF
DO I=1,K
D(I)=ZERO
END DO
END IF
PARL=PAR
PARLEST=PAR-DHD/DSQ
!
!     Terminate with D set to a multiple of the current D if the following
!     test suggests that it suitable to do so.
!
190 TEMP=PARUEST
IF (GSQ .EQ. ZERO) TEMP=TEMP*(ONE-TOL)
IF (PARUEST .GT. ZERO .AND. PARLEST .GE. TEMP) THEN
DTG=ZERO
DO I=1,N
DTG=DTG+D(I)*GG(I)
END DO
SCALE=-DSIGN(DELTA/DSQRT(DSQ),DTG)
DO I=1,N
D(I)=SCALE*D(I)
END DO
GOTO 370
END IF
!
!     Pick the value of PAR for the next iteration.
!
220 IF (PARU .EQ. ZERO) THEN
PAR=TWO*PARLEST+GNORM/DELTA
ELSE
PAR=0.5D0*(PARL+PARU)
PAR=DMAX1(PAR,PARLEST)
END IF
IF (PARUEST .GT. ZERO) PAR=DMIN1(PAR,PARUEST)
GOTO 140
!
!     Calculate D for the current PAR in the positive definite case.
!
230 W(1)=-GG(1)/PIV(1)
DO I=2,N
W(I)=(-GG(I)-TN(I-1)*W(I-1))/PIV(I)
END DO
D(N)=W(N)
DO I=NM,1,-1
D(I)=W(I)-TN(I)*D(I+1)/PIV(I)
END DO
!
!     Branch if a Newton-Raphson step is acceptable.
!
DSQ=ZERO
WSQ=ZERO
DO I=1,N
DSQ=DSQ+D(I)**2
WSQ=WSQ+PIV(I)*W(I)**2
END DO
IF (PAR .EQ. ZERO .AND. DSQ .LE. DELSQ) GOTO 320
!
!     Make the usual test for acceptability of a full trust region step.
!
DNORM=DSQRT(DSQ)
PHI=ONE/DNORM-ONE/DELTA
TEMP=TOL*(ONE+PAR*DSQ/WSQ)-DSQ*PHI*PHI
IF (TEMP .GE. ZERO) THEN
SCALE=DELTA/DNORM
DO I=1,N
D(I)=SCALE*D(I)
END DO
GOTO 370
END IF
IF (ITERC .GE. 2 .AND. PAR .LE. PARL) GOTO 370
IF (PARU .GT. ZERO .AND. PAR .GE. PARU) GOTO 370
!
!     Complete the iteration when PHI is negative.
!
IF (PHI .LT. ZERO) THEN
PARLEST=PAR
IF (POSDEF.EQ. ONE) THEN
IF (PHI .LE. PHIL) GOTO 370
SLOPE=(PHI-PHIL)/(PAR-PARL)
PARLEST=PAR-PHI/SLOPE
END IF
SLOPE=ONE/GNORM
IF (PARU .GT. ZERO) SLOPE=(PHIU-PHI)/(PARU-PAR)
TEMP=PAR-PHI/SLOPE
IF (PARUEST .GT. ZERO) TEMP=DMIN1(TEMP,PARUEST)
PARUEST=TEMP
POSDEF=ONE
PARL=PAR
PHIL=PHI
GOTO 220
END IF
!
!     If required, calculate Z for the alternative test for convergence.
!
IF (POSDEF .EQ. ZERO) THEN
W(1)=ONE/PIV(1)
DO I=2,N
TEMP=-TN(I-1)*W(I-1)
W(I)=(DSIGN(ONE,TEMP)+TEMP)/PIV(I)
END DO
Z(N)=W(N)
DO I=NM,1,-1
Z(I)=W(I)-TN(I)*Z(I+1)/PIV(I)
END DO
WWSQ=ZERO
ZSQ=ZERO
DTZ=ZERO
DO I=1,N
WWSQ=WWSQ+PIV(I)*W(I)**2
ZSQ=ZSQ+Z(I)**2
DTZ=DTZ+D(I)*Z(I)
END DO
!
!     Apply the alternative test for convergence.
!
TEMPA=DABS(DELSQ-DSQ)
TEMPB=DSQRT(DTZ*DTZ+TEMPA*ZSQ)
GAM=TEMPA/(DSIGN(TEMPB,DTZ)+DTZ)
TEMP=TOL*(WSQ+PAR*DELSQ)-GAM*GAM*WWSQ
IF (TEMP .GE. ZERO) THEN
DO I=1,N
D(I)=D(I)+GAM*Z(I)
END DO
GOTO 370
END IF
PARLEST=DMAX1(PARLEST,PAR-WWSQ/ZSQ)
END IF
!
!     Complete the iteration when PHI is positive.
!
SLOPE=ONE/GNORM
IF (PARU .GT. ZERO) THEN
IF (PHI .GE. PHIU) GOTO 370
SLOPE=(PHIU-PHI)/(PARU-PAR)
END IF
PARLEST=DMAX1(PARLEST,PAR-PHI/SLOPE)
PARUEST=PAR
IF (POSDEF .EQ. ONE) THEN
SLOPE=(PHI-PHIL)/(PAR-PARL)
PARUEST=PAR-PHI/SLOPE
END IF
PARU=PAR
PHIU=PHI
GOTO 220
!
!     Set EVALUE to the least eigenvalue of the second derivative matrix if
!     D is a Newton-Raphson step. SHFMAX will be an upper bound on EVALUE.
!
320 SHFMIN=ZERO
PIVOT=TD(1)
SHFMAX=PIVOT
DO K=2,N
PIVOT=TD(K)-TN(K-1)**2/PIVOT
SHFMAX=DMIN1(SHFMAX,PIVOT)
END DO
!
!     Find EVALUE by a bisection method, but occasionally SHFMAX may be
!     adjusted by the rule of false position.
!
KSAVE=0
340 SHIFT=0.5D0*(SHFMIN+SHFMAX)
K=1
TEMP=TD(1)-SHIFT
350 IF (TEMP .GT. ZERO) THEN
PIV(K)=TEMP
IF (K .LT. N) THEN
TEMP=TD(K+1)-SHIFT-TN(K)**2/TEMP
K=K+1
GOTO 350
END IF
SHFMIN=SHIFT
ELSE
IF (K .LT. KSAVE) GOTO 360
IF (K .EQ. KSAVE) THEN
IF (PIVKSV .EQ. ZERO) GOTO 360
IF (PIV(K)-TEMP .LT. TEMP-PIVKSV) THEN
PIVKSV=TEMP
SHFMAX=SHIFT
ELSE
PIVKSV=ZERO
SHFMAX=(SHIFT*PIV(K)-SHFMIN*TEMP)/(PIV(K)-TEMP)
END IF
ELSE
KSAVE=K
PIVKSV=TEMP
SHFMAX=SHIFT
END IF
END IF
IF (SHFMIN .LE. 0.99D0*SHFMAX) GOTO 340
360 EVALUE=SHFMIN
!
!     Apply the inverse Householder transformations to D.
!
370 NM=N-1
DO K=NM,1,-1
KP=K+1
SUM=ZERO
DO I=KP,N
SUM=SUM+D(I)*H(I,K)
END DO
DO I=KP,N
D(I)=D(I)-SUM*H(I,K)
END DO
END DO
!
!     Return from the subroutine.
!
400 RETURN
END


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% uobyqa.f %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine uobyqa_core (N,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W, IERR)
implicit integer (i-n), real(real64) (a-h,o-z)
!JN Declare IERR
INTEGER IERR
DIMENSION X(*),W(*)
!
!     This subroutine seeks the least value of a function of many variables,
!     by a trust region method that forms quadratic models by interpolation.
!     The algorithm is described in "UOBYQA: unconstrained optimization by
!     quadratic approximation" by M.J.D. Powell, Report DAMTP 2000/NA14,
!     University of Cambridge. The arguments of the subroutine are as follows.
!
!     N must be set to the number of variables and must be at least two.
!     Initial values of the variables must be set in X(1),X(2),...,X(N). They
!       will be changed to the values that give the least calculated F.
!     RHOBEG and RHOEND must be set to the initial and final values of a trust
!       region radius, so both must be positive with RHOEND<=RHOBEG. Typically
!       RHOBEG should be about one tenth of the greatest expected change to a
!       variable, and RHOEND should indicate the accuracy that is required in
!       the final values of the variables.
!     The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
!       amount of printing. Specifically, there is no output if IPRINT=0 and
!       there is output only at the return if IPRINT=1. Otherwise, each new
!       value of RHO is printed, with the best vector of variables so far and
!       the corresponding value of the objective function. Further, each new
!       value of F with its variables are output if IPRINT=3.
!     MAXFUN must be set to an upper bound on the number of calls of CALFUN.
!     The array W will be used for working space. Its length must be at least
!       ( N**4 + 8*N**3 + 23*N**2 + 42*N + max [ 2*N**2 + 4, 18*N ] ) / 4.
!JN   Add IERR to tell what the error is to minqer.
!     IERR gives an error code that is interpreted by minqer interface routine.
!       and passed back to minqa.R
!
!     SUBROUTINE CALFUN (N,X,F) must be provided by the user. It must set F to
!     the value of the objective function for the variables X(1),X(2),...,X(N).
!
!     Partition the working space array, so that different parts of it can be
!     treated separately by the subroutine that performs the main calculation.
!
NPT=(N*N+3*N+2)/2
IXB=1
IXO=IXB+N
IXN=IXO+N
IXP=IXN+N
IPQ=IXP+N*NPT
IPL=IPQ+NPT-1
IH=IPL+(NPT-1)*NPT
IG=IH+N*N
ID=IG+N
IVL=IH
IW=ID+N
!JN Initialize IERR to 0 (normal exit)
IERR=0
CALL UOBYQB (N,X,RHOBEG,RHOEND,IPRINT,MAXFUN,NPT,W(IXB),W(IXO), W(IXN),W(IXP),W(IPQ),W(IPL),W(IH),W(IG), &
    & W(ID),W(IVL),W(IW),IERR)
RETURN
END


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% uobyqb.f %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE UOBYQB (N,X,RHOBEG,RHOEND,IPRINT,MAXFUN,NPT,XBASE, XOPT,XNEW,XPT,PQ,PL,H,G,D,VLAG,W,IERR)
implicit integer (i-n), real(real64) (a-h,o-z)
!JN Declare IERR
INTEGER IERR
DIMENSION X(*),XBASE(*),XOPT(*),XNEW(*),XPT(NPT,*),PQ(*), PL(NPT,*),H(N,*),G(*),D(*),VLAG(*),W(*)
!
!     The arguments N, X, RHOBEG, RHOEND, IPRINT and MAXFUN are identical to
!       the corresponding arguments in SUBROUTINE UOBYQA.
!     NPT is set by UOBYQA to (N*N+3*N+2)/2 for the above dimension statement.
!     XBASE will contain a shift of origin that reduces the contributions from
!       rounding errors to values of the model and Lagrange functions.
!     XOPT will be set to the displacement from XBASE of the vector of
!       variables that provides the least calculated F so far.
!     XNEW will be set to the displacement from XBASE of the vector of
!       variables for the current calculation of F.
!     XPT will contain the interpolation point coordinates relative to XBASE.
!     PQ will contain the parameters of the quadratic model.
!     PL will contain the parameters of the Lagrange functions.
!     H will provide the second derivatives that TRSTEP and LAGMAX require.
!     G will provide the first derivatives that TRSTEP and LAGMAX require.
!     D is reserved for trial steps from XOPT, except that it will contain
!       diagonal second derivatives during the initialization procedure.
!     VLAG will contain the values of the Lagrange functions at a new point X.
!     The array W will be used for working space. Its length must be at least
!     max [ 6*N, ( N**2 + 3*N + 2 ) / 2 ].
!
!     Set some constants.
!
ONE=1.0D0
TWO=2.0D0
ZERO=0.0D0
HALF=0.5D0
TOL=0.01D0
NNP=N+N+1
NPTM=NPT-1
NFTEST=MAX0(MAXFUN,1)
!
!     Initialization. NF is the number of function calculations so far.
!
RHO=RHOBEG
RHOSQ=RHO*RHO
NF=0
DO I=1,N
XBASE(I)=X(I)
DO K=1,NPT
XPT(K,I)=ZERO
END DO
END DO
DO K=1,NPT
DO J=1,NPTM
PL(K,J)=ZERO
END DO
END DO
!
!     The branch to label 120 obtains a new value of the objective function
!     and then there is a branch back to label 50, because the new function
!     value is needed to form the initial quadratic model. The least function
!     value so far and its index are noted below.
!
30 DO I=1,N
X(I)=XBASE(I)+XPT(NF+1,I)
END DO
GOTO 120
50 IF (NF .EQ. 1) THEN
FOPT=F
KOPT=NF
FBASE=F
J=0
JSWITCH=-1
IH=N
ELSE
IF (F .LT. FOPT) THEN
FOPT=F
KOPT=NF
END IF
END IF
!
!     Form the gradient and diagonal second derivatives of the initial
!     quadratic model and Lagrange functions.
!
IF (NF .LE. NNP) THEN
JSWITCH=-JSWITCH
IF (JSWITCH .GT. 0) THEN
IF (J .GE. 1) THEN
IH=IH+J
IF (W(J) .LT. ZERO) THEN
D(J)=(FSAVE+F-TWO*FBASE)/RHOSQ
PQ(J)=(FSAVE-F)/(TWO*RHO)
PL(1,IH)=-TWO/RHOSQ
PL(NF-1,J)=HALF/RHO
PL(NF-1,IH)=ONE/RHOSQ
ELSE
PQ(J)=(4.0D0*FSAVE-3.0D0*FBASE-F)/(TWO*RHO)
D(J)=(FBASE+F-TWO*FSAVE)/RHOSQ
PL(1,J)=-1.5D0/RHO
PL(1,IH)=ONE/RHOSQ
PL(NF-1,J)=TWO/RHO
PL(NF-1,IH)=-TWO/RHOSQ
END IF
PQ(IH)=D(J)
PL(NF,J)=-HALF/RHO
PL(NF,IH)=ONE/RHOSQ
END IF
!
!     Pick the shift from XBASE to the next initial interpolation point
!     that provides diagonal second derivatives.
!
IF (J .LT. N) THEN
J=J+1
XPT(NF+1,J)=RHO
END IF
ELSE
FSAVE=F
IF (F .LT. FBASE) THEN
W(J)=RHO
XPT(NF+1,J)=TWO*RHO
ELSE
W(J)=-RHO
XPT(NF+1,J)=-RHO
END IF
END IF
IF (NF .LT. NNP) GOTO 30
!
!     Form the off-diagonal second derivatives of the initial quadratic model.
!
IH=N
IP=1
IQ=2
END IF
IH=IH+1
IF (NF .GT. NNP) THEN
TEMP=ONE/(W(IP)*W(IQ))
TEMPA=F-FBASE-W(IP)*PQ(IP)-W(IQ)*PQ(IQ)
PQ(IH)=(TEMPA-HALF*RHOSQ*(D(IP)+D(IQ)))*TEMP
PL(1,IH)=TEMP
IW=IP+IP
IF (W(IP) .LT. ZERO) IW=IW+1
PL(IW,IH)=-TEMP
IW=IQ+IQ
IF (W(IQ) .LT. ZERO) IW=IW+1
PL(IW,IH)=-TEMP
PL(NF,IH)=TEMP
!
!     Pick the shift from XBASE to the next initial interpolation point
!     that provides off-diagonal second derivatives.
!
IP=IP+1
END IF
IF (IP .EQ. IQ) THEN
IH=IH+1
IP=1
IQ=IQ+1
END IF
IF (NF .LT. NPT) THEN
XPT(NF+1,IP)=W(IP)
XPT(NF+1,IQ)=W(IQ)
GOTO 30
END IF
!
!     Set parameters to begin the iterations for the current RHO.
!
SIXTHM=ZERO
DELTA=RHO
60 TWORSQ=(TWO*RHO)**2
RHOSQ=RHO*RHO
!
!     Form the gradient of the quadratic model at the trust region centre.
!
70 KNEW=0
IH=N
DO J=1,N
XOPT(J)=XPT(KOPT,J)
G(J)=PQ(J)
DO I=1,J
IH=IH+1
G(I)=G(I)+PQ(IH)*XOPT(J)
IF (I .LT. J) G(J)=G(J)+PQ(IH)*XOPT(I)
H(I,J)=PQ(IH)
END DO
END DO
!
!     Generate the next trust region step and test its length. Set KNEW
!     to -1 if the purpose of the next F will be to improve conditioning,
!     and also calculate a lower bound on the Hessian term of the model Q.
!
CALL TRSTEP (N,G,H,DELTA,TOL,D,W(1),W(N+1),W(2*N+1),W(3*N+1), W(4*N+1),W(5*N+1),EVALUE)
TEMP=ZERO
DO I=1,N
TEMP=TEMP+D(I)**2
END DO
DNORM=DMIN1(DELTA,DSQRT(TEMP))
ERRTOL=-ONE
IF (DNORM .LT. HALF*RHO) THEN
KNEW=-1
ERRTOL=HALF*EVALUE*RHO*RHO
IF (NF .LE. NPT+9) ERRTOL=ZERO
GOTO 290
END IF
!
!     Calculate the next value of the objective function.
!
100 DO I=1,N
XNEW(I)=XOPT(I)+D(I)
X(I)=XBASE(I)+XNEW(I)
END DO
120 IF (NF .GE. NFTEST) THEN
!$$$          IF (IPRINT .GT. 0) PRINT 130
!$$$  130     FORMAT (/4X,'Return from UOBYQA because CALFUN has been',
!$$$     1      ' called MAXFUN times')
!$$$          GOTO 420
!JN Too many function evaluations
IERR=390
GOTO 420
END IF
NF=NF+1
F = CALFUN (N,X,IPRINT)
!$$$      IF (IPRINT .EQ. 3) THEN
!$$$          PRINT 70, NF,F,(X(I),I=1,N)
!$$$   70      FORMAT (/4X,'Function number',I6,'    F =',1PD18.10,
!$$$     1       '    The corresponding X is:'/(2X,5D15.6))
!$$$      END IF
!$$$      CALL minqi3 (IPRINT, F, NF, N, X)
IF (NF .LE. NPT) GOTO 50
IF (KNEW .EQ. -1) GOTO 420
!
!     Use the quadratic model to predict the change in F due to the step D,
!     and find the values of the Lagrange functions at the new point.
!
VQUAD=ZERO
IH=N
DO J=1,N
W(J)=D(J)
VQUAD=VQUAD+W(J)*PQ(J)
DO I=1,J
IH=IH+1
W(IH)=D(I)*XNEW(J)+D(J)*XOPT(I)
IF (I .EQ. J) W(IH)=HALF*W(IH)
VQUAD=VQUAD+W(IH)*PQ(IH)
END DO
END DO
DO K=1,NPT
TEMP=ZERO
DO J=1,NPTM
TEMP=TEMP+W(J)*PL(K,J)
END DO
VLAG(K)=TEMP
END DO
VLAG(KOPT)=VLAG(KOPT)+ONE
!
!     Update SIXTHM, which is a lower bound on one sixth of the greatest
!     third derivative of F.
!
DIFF=F-FOPT-VQUAD
SUM=ZERO
DO K=1,NPT
TEMP=ZERO
DO I=1,N
TEMP=TEMP+(XPT(K,I)-XNEW(I))**2
END DO
TEMP=DSQRT(TEMP)
SUM=SUM+DABS(TEMP*TEMP*TEMP*VLAG(K))
END DO
SIXTHM=DMAX1(SIXTHM,DABS(DIFF)/SUM)
!
!     Update FOPT and XOPT if the new F is the least value of the objective
!     function so far. Then branch if D is not a trust region step.
!
FSAVE=FOPT
IF (F .LT. FOPT) THEN
FOPT=F
DO I=1,N
XOPT(I)=XNEW(I)
END DO
END IF
KSAVE=KNEW
IF (KNEW .GT. 0) GOTO 240
!
!     Pick the next value of DELTA after a trust region step.
!
IF (VQUAD .GE. ZERO) THEN
!JN              IF (IPRINT .GT. 0) CALL minqer(2101)
!$$$          IF (IPRINT .GT. 0) PRINT 210
!$$$  210     FORMAT (/4X,'Return from UOBYQA because a trust',
!$$$     1      ' region step has failed to reduce Q')
IERR=2101
GOTO 420
END IF
RATIO=(F-FSAVE)/VQUAD
IF (RATIO .LE. 0.1D0) THEN
DELTA=HALF*DNORM
ELSE IF (RATIO.LE. 0.7D0) THEN
DELTA=DMAX1(HALF*DELTA,DNORM)
ELSE
DELTA=DMAX1(DELTA,1.25D0*DNORM,DNORM+RHO)
END IF
IF (DELTA .LE. 1.5D0*RHO) DELTA=RHO
!
!     Set KNEW to the index of the next interpolation point to be deleted.
!
KTEMP=0
DETRAT=ZERO
IF (F .GE. FSAVE) THEN
KTEMP=KOPT
DETRAT=ONE
END IF
DO K=1,NPT
SUM=ZERO
DO I=1,N
SUM=SUM+(XPT(K,I)-XOPT(I))**2
END DO
TEMP=DABS(VLAG(K))
IF (SUM .GT. RHOSQ) TEMP=TEMP*(SUM/RHOSQ)**1.5D0
IF (TEMP .GT. DETRAT .AND. K .NE. KTEMP) THEN
DETRAT=TEMP
DDKNEW=SUM
KNEW=K
END IF
END DO
IF (KNEW .EQ. 0) GOTO 290
!
!     Replace the interpolation point that has index KNEW by the point XNEW,
!     and also update the Lagrange functions and the quadratic model.
!
240 DO I=1,N
XPT(KNEW,I)=XNEW(I)
END DO
TEMP=ONE/VLAG(KNEW)
DO J=1,NPTM
PL(KNEW,J)=TEMP*PL(KNEW,J)
PQ(J)=PQ(J)+DIFF*PL(KNEW,J)
END DO
DO K=1,NPT
IF (K .NE. KNEW) THEN
TEMP=VLAG(K)
DO J=1,NPTM
PL(K,J)=PL(K,J)-TEMP*PL(KNEW,J)
END DO
END IF
END DO
!
!     Update KOPT if F is the least calculated value of the objective
!     function. Then branch for another trust region calculation. The
!     case KSAVE>0 indicates that a model step has just been taken.
!
IF (F .LT. FSAVE) THEN
KOPT=KNEW
GOTO 70
END IF
IF (KSAVE .GT. 0) GOTO 70
IF (DNORM .GT. TWO*RHO) GOTO 70
IF (DDKNEW .GT. TWORSQ) GOTO 70
!
!     Alternatively, find out if the interpolation points are close
!     enough to the best point so far.
!
290 DO K=1,NPT
W(K)=ZERO
DO I=1,N
W(K)=W(K)+(XPT(K,I)-XOPT(I))**2
END DO
END DO
310 KNEW=-1
DISTEST=TWORSQ
DO K=1,NPT
IF (W(K) .GT. DISTEST) THEN
KNEW=K
DISTEST=W(K)
END IF
END DO
!
!     If a point is sufficiently far away, then set the gradient and Hessian
!     of its Lagrange function at the centre of the trust region, and find
!     half the sum of squares of components of the Hessian.
!
IF (KNEW .GT. 0) THEN
IH=N
SUMH=ZERO
DO J=1,N
G(J)=PL(KNEW,J)
DO I=1,J
IH=IH+1
TEMP=PL(KNEW,IH)
G(J)=G(J)+TEMP*XOPT(I)
IF (I .LT. J) THEN
G(I)=G(I)+TEMP*XOPT(J)
SUMH=SUMH+TEMP*TEMP
END IF
H(I,J)=TEMP
END DO
SUMH=SUMH+HALF*TEMP*TEMP
END DO
!
!     If ERRTOL is positive, test whether to replace the interpolation point
!     with index KNEW, using a bound on the maximum modulus of its Lagrange
!     function in the trust region.
!
IF (ERRTOL .GT. ZERO) THEN
W(KNEW)=ZERO
SUMG=ZERO
DO I=1,N
SUMG=SUMG+G(I)**2
END DO
ESTIM=RHO*(DSQRT(SUMG)+RHO*DSQRT(HALF*SUMH))
WMULT=SIXTHM*DISTEST**1.5D0
IF (WMULT*ESTIM .LE. ERRTOL) GOTO 310
END IF
!
!     If the KNEW-th point may be replaced, then pick a D that gives a large
!     value of the modulus of its Lagrange function within the trust region.
!     Here the vector XNEW is used as temporary working space.
!
CALL LAGMAX (N,G,H,RHO,D,XNEW,VMAX)
IF (ERRTOL .GT. ZERO) THEN
IF (WMULT*VMAX .LE. ERRTOL) GOTO 310
END IF
GOTO 100
END IF
IF (DNORM .GT. RHO) GOTO 70
!
!     Prepare to reduce RHO by shifting XBASE to the best point so far,
!     and make the corresponding changes to the gradients of the Lagrange
!     functions and the quadratic model.
!
IF (RHO .GT. RHOEND) THEN
IH=N
DO J=1,N
XBASE(J)=XBASE(J)+XOPT(J)
DO K=1,NPT
XPT(K,J)=XPT(K,J)-XOPT(J)
END DO
DO I=1,J
IH=IH+1
PQ(I)=PQ(I)+PQ(IH)*XOPT(J)
IF (I .LT. J) THEN
PQ(J)=PQ(J)+PQ(IH)*XOPT(I)
DO K=1,NPT
PL(K,J)=PL(K,J)+PL(K,IH)*XOPT(I)
END DO
END IF
DO K=1,NPT
PL(K,I)=PL(K,I)+PL(K,IH)*XOPT(J)
END DO
END DO
END DO
!
!     Pick the next values of RHO and DELTA.
!
DELTA=HALF*RHO
RATIO=RHO/RHOEND
IF (RATIO .LE. 16.0D0) THEN
RHO=RHOEND
ELSE IF (RATIO .LE. 250.0D0) THEN
RHO=DSQRT(RATIO)*RHOEND
ELSE
RHO=0.1D0*RHO
END IF
DELTA=DMAX1(DELTA,RHO)
IF (IPRINT .GE. 2) THEN
CALL minqit(IPRINT, RHO, NF, FOPT, N, XBASE, XOPT)
!$$$              IF (IPRINT .GE. 3) PRINT 390
!$$$  390         FORMAT (5X)
!$$$              PRINT 400, RHO,NF
!$$$  400         FORMAT (/4X,'New RHO =',1PD11.4,5X,'Number of',
!$$$     1          ' function values =',I6)
!$$$              PRINT 410, FOPT,(XBASE(I),I=1,N)
!$$$  410         FORMAT (4X,'Least value of F =',1PD23.15,9X,
!$$$     1          'The corresponding X is:'/(2X,5D15.6))
END IF
GOTO 60
END IF
!
!     Return from the calculation, after another Newton-Raphson step, if
!     it is too short to have been tried before.
!
IF (ERRTOL .GE. ZERO) GOTO 100
420 IF (FOPT .LE. F) THEN
DO I=1,N
X(I)=XBASE(I)+XOPT(I)
END DO
F=FOPT
END IF
IF (IPRINT .GE. 1) THEN
CALL minqir(IPRINT, F, NF, N, X)
!$$$          PRINT 440, NF
!$$$  440     FORMAT (/4X,'At the return from UOBYQA',5X,
!$$$     1      'Number of function values =',I6)
!$$$          PRINT 410, F,(X(I),I=1,N)
END IF
RETURN
END


SUBROUTINE UPDATE (N,NPT,BMAT,ZMAT,IDZ,NDIM,VLAG,BETA,KNEW,W)
implicit integer (i-n), real(real64) (a-h,o-z)
DIMENSION BMAT(NDIM,*),ZMAT(NPT,*),VLAG(*),W(*)
!
!     The arrays BMAT and ZMAT with IDZ are updated, in order to shift the
!     interpolation point that has index KNEW. On entry, VLAG contains the
!     components of the vector Theta*Wcheck+e_b of the updating formula
!     (6.11), and BETA holds the value of the parameter that has this name.
!     The vector W is used for working space.
!
!     Set some constants.
!
ONE=1.0D0
ZERO=0.0D0
NPTM=NPT-N-1
!
!     Apply the rotations that put zeros in the KNEW-th row of ZMAT.
!
JL=1
DO J=2,NPTM
IF (J .EQ. IDZ) THEN
JL=IDZ
ELSE IF (ZMAT(KNEW,J) .NE. ZERO) THEN
TEMP=DSQRT(ZMAT(KNEW,JL)**2+ZMAT(KNEW,J)**2)
TEMPA=ZMAT(KNEW,JL)/TEMP
TEMPB=ZMAT(KNEW,J)/TEMP
DO I=1,NPT
TEMP=TEMPA*ZMAT(I,JL)+TEMPB*ZMAT(I,J)
ZMAT(I,J)=TEMPA*ZMAT(I,J)-TEMPB*ZMAT(I,JL)
ZMAT(I,JL)=TEMP
END DO
ZMAT(KNEW,J)=ZERO
END IF
END DO
!
!     Put the first NPT components of the KNEW-th column of HLAG into W,
!     and calculate the parameters of the updating formula.
!
TEMPA=ZMAT(KNEW,1)
IF (IDZ .GE. 2) TEMPA=-TEMPA
IF (JL .GT. 1) TEMPB=ZMAT(KNEW,JL)
DO I=1,NPT
W(I)=TEMPA*ZMAT(I,1)
IF (JL .GT. 1) W(I)=W(I)+TEMPB*ZMAT(I,JL)
END DO
ALPHA=W(KNEW)
TAU=VLAG(KNEW)
TAUSQ=TAU*TAU
DENOM=ALPHA*BETA+TAUSQ
VLAG(KNEW)=VLAG(KNEW)-ONE
!
!     Complete the updating of ZMAT when there is only one nonzero element
!     in the KNEW-th row of the new matrix ZMAT, but, if IFLAG is set to one,
!     then the first column of ZMAT will be exchanged with another one later.
!
IFLAG=0
IF (JL .EQ. 1) THEN
TEMP=DSQRT(DABS(DENOM))
TEMPB=TEMPA/TEMP
TEMPA=TAU/TEMP
DO I=1,NPT
ZMAT(I,1)=TEMPA*ZMAT(I,1)-TEMPB*VLAG(I)
END DO
IF (IDZ .EQ. 1 .AND. TEMP .LT. ZERO) IDZ=2
IF (IDZ .GE. 2 .AND. TEMP .GE. ZERO) IFLAG=1
ELSE
!
!     Complete the updating of ZMAT in the alternative case.
!
JA=1
IF (BETA .GE. ZERO) JA=JL
JB=JL+1-JA
TEMP=ZMAT(KNEW,JB)/DENOM
TEMPA=TEMP*BETA
TEMPB=TEMP*TAU
TEMP=ZMAT(KNEW,JA)
SCALA=ONE/DSQRT(DABS(BETA)*TEMP*TEMP+TAUSQ)
SCALB=SCALA*DSQRT(DABS(DENOM))
DO I=1,NPT
ZMAT(I,JA)=SCALA*(TAU*ZMAT(I,JA)-TEMP*VLAG(I))
ZMAT(I,JB)=SCALB*(ZMAT(I,JB)-TEMPA*W(I)-TEMPB*VLAG(I))
END DO
IF (DENOM .LE. ZERO) THEN
IF (BETA .LT. ZERO) IDZ=IDZ+1
IF (BETA .GE. ZERO) IFLAG=1
END IF
END IF
!
!     IDZ is reduced in the following case, and usually the first column
!     of ZMAT is exchanged with a later one.
!
IF (IFLAG .EQ. 1) THEN
IDZ=IDZ-1
DO I=1,NPT
TEMP=ZMAT(I,1)
ZMAT(I,1)=ZMAT(I,IDZ)
ZMAT(I,IDZ)=TEMP
END DO
END IF
!
!     Finally, update the matrix BMAT.
!
DO J=1,N
JP=NPT+J
W(JP)=BMAT(KNEW,J)
TEMPA=(ALPHA*VLAG(JP)-TAU*W(JP))/DENOM
TEMPB=(-BETA*W(JP)-TAU*VLAG(JP))/DENOM
DO I=1,JP
BMAT(I,J)=BMAT(I,J)+TEMPA*VLAG(I)+TEMPB*W(I)
IF (I .GT. NPT) BMAT(JP,I-NPT)=BMAT(I,J)
END DO
END DO
RETURN
END


SUBROUTINE UPDATEBOBYQA (N,NPT,BMAT,ZMAT,NDIM,VLAG,BETA,DENOM, KNEW,W)
implicit integer (i-n), real(real64) (a-h,o-z)
DIMENSION BMAT(NDIM,*),ZMAT(NPT,*),VLAG(*),W(*)
!
!     The arrays BMAT and ZMAT are updated, as required by the new position
!     of the interpolation point that has the index KNEW. The vector VLAG has
!     N+NPT components, set on entry to the first NPT and last N components
!     of the product Hw in equation (4.11) of the Powell (2006) paper on
!     NEWUOA. Further, BETA is set on entry to the value of the parameter
!     with that name, and DENOM is set to the denominator of the updating
!     formula. Elements of ZMAT may be treated as zero if their moduli are
!     at most ZTEST. The first NDIM elements of W are used for working space.
!
!     Set some constants.
!
ONE=1.0D0
ZERO=0.0D0
NPTM=NPT-N-1
ZTEST=ZERO
DO K=1,NPT
DO J=1,NPTM
ZTEST=DMAX1(ZTEST,DABS(ZMAT(K,J)))
END DO
END DO
ZTEST=1.0D-20*ZTEST
!
!     Apply the rotations that put zeros in the KNEW-th row of ZMAT.
!
JL=1
DO J=2,NPTM
IF (DABS(ZMAT(KNEW,J)) .GT. ZTEST) THEN
TEMP=DSQRT(ZMAT(KNEW,1)**2+ZMAT(KNEW,J)**2)
TEMPA=ZMAT(KNEW,1)/TEMP
TEMPB=ZMAT(KNEW,J)/TEMP
DO I=1,NPT
TEMP=TEMPA*ZMAT(I,1)+TEMPB*ZMAT(I,J)
ZMAT(I,J)=TEMPA*ZMAT(I,J)-TEMPB*ZMAT(I,1)
ZMAT(I,1)=TEMP
END DO
END IF
ZMAT(KNEW,J)=ZERO
END DO
!
!     Put the first NPT components of the KNEW-th column of HLAG into W,
!     and calculate the parameters of the updating formula.
!
DO I=1,NPT
W(I)=ZMAT(KNEW,1)*ZMAT(I,1)
END DO
ALPHA=W(KNEW)
TAU=VLAG(KNEW)
VLAG(KNEW)=VLAG(KNEW)-ONE
!
!     Complete the updating of ZMAT.
!
TEMP=DSQRT(DENOM)
TEMPB=ZMAT(KNEW,1)/TEMP
TEMPA=TAU/TEMP
DO I=1,NPT
ZMAT(I,1)=TEMPA*ZMAT(I,1)-TEMPB*VLAG(I)
END DO
!
!     Finally, update the matrix BMAT.
!
DO J=1,N
JP=NPT+J
W(JP)=BMAT(KNEW,J)
TEMPA=(ALPHA*VLAG(JP)-TAU*W(JP))/DENOM
TEMPB=(-BETA*W(JP)-TAU*VLAG(JP))/DENOM
DO I=1,JP
BMAT(I,J)=BMAT(I,J)+TEMPA*VLAG(I)+TEMPB*W(I)
IF (I .GT. NPT) BMAT(JP,I-NPT)=BMAT(I,J)
END DO
END DO
RETURN
END


end module minqa_module
