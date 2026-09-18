module pbs_api
   use pbs_kinds, only : dp
   use pbs_splines, only : bspline_design, quantile_type7, sort_real
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
   implicit none
   private

   integer, parameter, public :: PBS_OK = 0
   integer, parameter, public :: PBS_INVALID_ARGUMENT = 1
   integer, parameter, public :: PBS_OUTSIDE_BOUNDARY = 2

   type, public :: pbs_basis
      real(dp), allocatable :: values(:, :)
      real(dp), allocatable :: knots(:)
      real(dp) :: boundary_knots(2) = 0.0_dp
      integer :: degree = 3
      logical :: intercept = .false.
      logical :: periodic = .true.
      integer :: status = PBS_OK
      character(len=256) :: message = ''
   end type pbs_basis

   public :: dp
   public :: pbs
   public :: predict_pbs

contains

   pure function pbs(x, df, knots, degree, intercept, boundary_knots, periodic) result(out)
      real(dp), intent(in) :: x(:) !! Predictor values; IEEE NaNs are preserved as missing rows in the returned basis.
      integer, intent(in), optional :: df !! Requested basis degrees of freedom when knots are omitted.
      real(dp), intent(in), optional :: knots(:) !! Internal breakpoints; when present they take precedence over df.
      integer, intent(in), optional :: degree !! Piecewise-polynomial degree, default 3 and required to be at least 1.
      logical, intent(in), optional :: intercept !! Whether to retain the first basis column; default false.
      real(dp), intent(in), optional :: boundary_knots(2) !! Lower and upper spline boundaries; defaults to finite x range.
      logical, intent(in), optional :: periodic !! Whether to construct the periodic basis; default true.
      type(pbs_basis) :: out
      real(dp), allocatable :: active_x(:)
      real(dp), allocatable :: all_knots(:)
      real(dp), allocatable :: basis_full(:)
      real(dp), allocatable :: internal_knots(:)
      real(dp), allocatable :: temp(:, :)
      real(dp) :: nan_value
      real(dp) :: prob
      real(dp) :: scale
      real(dp) :: pivot
      real(dp) :: delta
      integer :: d
      integer :: deg
      integer :: i
      integer :: j
      integer :: n_active
      integer :: n_internal
      integer :: nbase
      integer :: ncol
      integer :: ord
      logical :: has_boundary
      logical :: keep_intercept
      logical :: is_periodic

      nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
      deg = 3
      if (present(degree)) deg = degree
      keep_intercept = .false.
      if (present(intercept)) keep_intercept = intercept
      is_periodic = .true.
      if (present(periodic)) is_periodic = periodic

      out%degree = deg
      out%intercept = keep_intercept
      out%periodic = is_periodic
      out%status = PBS_OK
      out%message = ''

      if (deg < 1) then
         call set_error(out, PBS_INVALID_ARGUMENT, "'degree' must be integer >= 1")
         return
      end if
      ord = deg + 1

      n_active = count(.not. ieee_is_nan(x))
      if (n_active <= 0) then
         call set_error(out, PBS_INVALID_ARGUMENT, 'x contains no finite non-missing values')
         return
      end if
      allocate(active_x(n_active))
      j = 0
      do i = 1, size(x)
         if (.not. ieee_is_nan(x(i))) then
            j = j + 1
            active_x(j) = x(i)
         end if
      end do

      has_boundary = present(boundary_knots)
      if (has_boundary) then
         out%boundary_knots = boundary_knots
         if (out%boundary_knots(1) > out%boundary_knots(2)) then
            pivot = out%boundary_knots(1)
            out%boundary_knots(1) = out%boundary_knots(2)
            out%boundary_knots(2) = pivot
         end if
      else
         out%boundary_knots = [minval(active_x), maxval(active_x)]
      end if
      if (out%boundary_knots(1) >= out%boundary_knots(2)) then
         call set_error(out, PBS_INVALID_ARGUMENT, 'boundary knots must be distinct')
         return
      end if

      if (present(knots)) then
         allocate(internal_knots(size(knots)))
         internal_knots = knots
         out%knots = knots
         if (size(internal_knots) > 1) call sort_real(internal_knots)
      else if (present(df)) then
         if (is_periodic) then
            n_internal = df - 1 + merge(0, 1, keep_intercept)
         else
            n_internal = df - ord + merge(0, 1, keep_intercept)
         end if
         n_internal = max(0, n_internal)
         allocate(internal_knots(n_internal))
         do i = 1, n_internal
            prob = real(i, dp) / real(n_internal + 1, dp)
            internal_knots(i) = quantile_inside(active_x, out%boundary_knots, has_boundary, prob)
         end do
      else
         allocate(internal_knots(0))
      end if
      if (.not. allocated(out%knots)) out%knots = internal_knots

      if (is_periodic) then
         if (size(internal_knots) < deg) then
            call set_error(out, PBS_INVALID_ARGUMENT, &
               'periodic basis requires at least degree internal knots')
            return
         end if
         if (any(active_x < out%boundary_knots(1)) .or. any(active_x > out%boundary_knots(2))) then
            call set_error(out, PBS_OUTSIDE_BOUNDARY, &
               "some 'x' values beyond boundary knots may cause ill-conditioned bases")
            return
         end if
         call periodic_knot_sequence(out%boundary_knots, internal_knots, deg, all_knots)
         nbase = size(all_knots) - ord
         ncol = nbase - deg
         if (.not. keep_intercept) ncol = ncol - 1
         allocate(out%values(size(x), ncol))
         out%values = nan_value
         allocate(basis_full(nbase))
         allocate(temp(1, nbase - deg))
         do i = 1, size(x)
            if (ieee_is_nan(x(i))) cycle
            call bspline_design(all_knots, x(i), deg, 0, basis_full)
            temp(1, :) = 0.0_dp
            if (nbase > 2 * deg) temp(1, 1:nbase - 2 * deg) = basis_full(deg + 1:nbase - deg)
            temp(1, nbase - 2 * deg + 1:nbase - deg) = basis_full(1:deg) + basis_full(nbase - deg + 1:nbase)
            if (keep_intercept) then
               out%values(i, :) = temp(1, :)
            else
               out%values(i, :) = temp(1, 2:)
            end if
         end do
      else
         call ordinary_knot_sequence(out%boundary_knots, internal_knots, ord, all_knots)
         nbase = size(all_knots) - ord
         ncol = nbase
         if (.not. keep_intercept) ncol = ncol - 1
         allocate(out%values(size(x), ncol))
         out%values = nan_value
         allocate(basis_full(nbase))
         allocate(temp(1, nbase))

         do i = 1, size(x)
            if (ieee_is_nan(x(i))) cycle
            if (x(i) < out%boundary_knots(1) .or. x(i) > out%boundary_knots(2)) then
               if (x(i) < out%boundary_knots(1)) then
                  pivot = out%boundary_knots(1)
               else
                  pivot = out%boundary_knots(2)
               end if
               delta = x(i) - pivot
               temp(1, :) = 0.0_dp
               scale = 1.0_dp
               do d = 0, deg
                  if (d > 1) scale = scale * real(d, dp)
                  call bspline_design(all_knots, pivot, deg, d, basis_full)
                  temp(1, :) = temp(1, :) + delta**d * basis_full / scale
               end do
            else
               call bspline_design(all_knots, x(i), deg, 0, basis_full)
               temp(1, :) = basis_full
            end if
            if (keep_intercept) then
               out%values(i, :) = temp(1, :)
            else
               out%values(i, :) = temp(1, 2:)
            end if
         end do
      end if
   end function pbs

   pure function predict_pbs(object, newx) result(out)
      type(pbs_basis), intent(in) :: object !! Previously constructed basis carrying spline knots and options.
      real(dp), intent(in) :: newx(:) !! New predictor values evaluated with the stored spline specification.
      type(pbs_basis) :: out

      out = pbs(newx, knots=object%knots, degree=object%degree, intercept=object%intercept, &
         boundary_knots=object%boundary_knots, periodic=object%periodic)
   end function predict_pbs

   pure subroutine set_error(out, status, message)
      type(pbs_basis), intent(inout) :: out !! Result object whose status and diagnostic message are set.
      integer, intent(in) :: status !! Nonzero package error code describing the validation failure.
      character(len=*), intent(in) :: message !! Human-readable explanation of the validation failure.

      out%status = status
      out%message = message
      if (.not. allocated(out%values)) allocate(out%values(0, 0))
      if (.not. allocated(out%knots)) allocate(out%knots(0))
   end subroutine set_error

   pure function quantile_inside(x, boundaries, restrict_boundary, prob) result(q)
      real(dp), intent(in) :: x(:) !! Non-missing predictor values used for knot selection.
      real(dp), intent(in) :: boundaries(2) !! Sorted lower and upper boundary knots.
      logical, intent(in) :: restrict_boundary !! Whether values outside explicit boundaries are excluded.
      real(dp), intent(in) :: prob !! Type-7 quantile probability in [0, 1].
      real(dp) :: q
      real(dp), allocatable :: inside(:)
      integer :: i
      integer :: n

      if (.not. restrict_boundary) then
         q = quantile_type7(x, prob)
         return
      end if
      n = count(x >= boundaries(1) .and. x <= boundaries(2))
      if (n <= 0) then
         q = 0.5_dp * (boundaries(1) + boundaries(2))
         return
      end if
      allocate(inside(n))
      n = 0
      do i = 1, size(x)
         if (x(i) >= boundaries(1) .and. x(i) <= boundaries(2)) then
            n = n + 1
            inside(n) = x(i)
         end if
      end do
      q = quantile_type7(inside, prob)
   end function quantile_inside

   pure subroutine periodic_knot_sequence(boundaries, internal_knots, degree, all_knots)
      real(dp), intent(in) :: boundaries(2) !! Sorted lower and upper period boundaries.
      real(dp), intent(in) :: internal_knots(:) !! Sorted internal breakpoints within the period.
      integer, intent(in) :: degree !! Polynomial degree controlling the number of wrapped knots.
      real(dp), allocatable, intent(out) :: all_knots(:) !! Extended periodic knot sequence for B-spline evaluation.
      real(dp), allocatable :: base(:)
      integer :: i
      integer :: n

      allocate(base(size(internal_knots) + 2))
      base = [boundaries(1), internal_knots, boundaries(2)]
      call sort_real(base)
      n = size(base)
      allocate(all_knots(n + 2 * degree))
      all_knots(degree + 1:degree + n) = base
      do i = 1, degree
         all_knots(degree + n + i) = base(n) + base(i + 1) - base(1)
         all_knots(degree + 1 - i) = base(1) - base(n) + base(n - i)
      end do
   end subroutine periodic_knot_sequence

   pure subroutine ordinary_knot_sequence(boundaries, internal_knots, ord, all_knots)
      real(dp), intent(in) :: boundaries(2) !! Sorted lower and upper spline boundaries.
      real(dp), intent(in) :: internal_knots(:) !! Sorted internal spline breakpoints.
      integer, intent(in) :: ord !! Spline order equal to degree plus one.
      real(dp), allocatable, intent(out) :: all_knots(:) !! Open knot sequence with repeated boundary knots.
      integer :: n

      n = size(internal_knots)
      allocate(all_knots(2 * ord + n))
      all_knots(1:ord) = boundaries(1)
      if (n > 0) all_knots(ord + 1:ord + n) = internal_knots
      all_knots(ord + n + 1:) = boundaries(2)
      call sort_real(all_knots)
   end subroutine ordinary_knot_sequence

end module pbs_api
