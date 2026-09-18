! SPDX-License-Identifier: MIT
! Modern Fortran translation of the computational API of R package kde1d.
! Upstream copyright: Thomas Nagler and Thibault Vatter.
! Numerical algorithms are derived from the MIT-licensed kde1d/kde1d-cpp sources.
module kde1d_api
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
   use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
   use r_kinds, only : dp
   use r_distributions, only : r_dnorm, r_pnorm, r_qnorm
   use r_quantiles, only : r_median
   use kde1d_bandwidth, only : gaussian_derivative, kde_drv_binned, linbin
   use kde1d_bandwidth, only : plugin_bandwidth, weighted_quantile_cpp
   use kde1d_interpolation, only : interpolation_grid
   implicit none
   private

   integer, parameter, public :: kde1d_continuous = 1
   integer, parameter, public :: kde1d_discrete = 2
   integer, parameter, public :: kde1d_zero_inflated = 3
   real(dp), parameter :: k0 = 0.3989425_dp

   type, public :: kde1d_model
      type(interpolation_grid) :: grid
      real(dp) :: xmin = 0.0_dp
      real(dp) :: xmax = 0.0_dp
      logical :: has_xmin = .false.
      logical :: has_xmax = .false.
      integer :: var_type = kde1d_continuous
      real(dp) :: mult = 1.0_dp
      real(dp) :: bw_spec = 0.0_dp
      logical :: bw_automatic = .true.
      real(dp) :: bw = 0.0_dp
      integer :: deg = 2
      integer :: grid_size = 400
      logical :: boundary_repair = .true.
      real(dp) :: prob0 = 0.0_dp
      real(dp) :: edf = 0.0_dp
      real(dp) :: loglik = 0.0_dp
      integer :: nobs = 0
      real(dp) :: boundary_offset = 0.0_dp
      real(dp) :: boundary_scale = 1.0_dp
      logical :: fitted = .false.
   end type kde1d_model

   interface equi_jitter
      module procedure equi_jitter_real
      module procedure equi_jitter_integer
   end interface equi_jitter

   public :: dp
   public :: dkde1d
   public :: equi_jitter
   public :: kde1d_fit
   public :: kde1d_loglik
   public :: kde1d_summary
   public :: pkde1d
   public :: qkde1d
   public :: rkde1d

contains

   pure function equi_jitter_real(x) result(y)
      real(dp), intent(in) :: x(:) !! Numeric R-style input, returned unchanged because R equi_jitter only jitters factors.
      real(dp) :: y(size(x))

      y = x
   end function equi_jitter_real

   pure function equi_jitter_integer(x) result(y)
      integer, intent(in) :: x(:) !! Ordered integer category codes to spread deterministically inside unit cells.
      real(dp) :: y(size(x))

      real(dp) :: xr(size(x))

      xr = real(x, dp)
      y = equi_jitter_codes(xr)
   end function equi_jitter_integer

   pure subroutine kde1d_fit(x, model, ierr, type_name, xmin, xmax, mult, bw, deg, weights, boundary_repair, grid_size)
      real(dp), intent(in) :: x(:) !! Observations; NaNs are omitted, and discrete observations must be integer-valued.
      type(kde1d_model), intent(out) :: model !! Fitted KDE object containing grid, controls, and diagnostics.
      integer, intent(out) :: ierr !! Zero on success; positive values identify invalid input or numerical failure.
      character(len=*), intent(in), optional :: type_name !! Variable type alias: continuous, discrete, or zero-inflated.
      real(dp), intent(in), optional :: xmin !! Finite lower support bound; absent means unbounded below.
      real(dp), intent(in), optional :: xmax !! Finite upper support bound; absent means unbounded above.
      real(dp), intent(in), optional :: mult !! Positive bandwidth multiplier; default 1.
      real(dp), intent(in), optional :: bw !! Positive fixed bandwidth; absent or NaN requests plug-in selection.
      integer, intent(in), optional :: deg !! Local log-polynomial degree 0, 1, or 2; default 2.
      real(dp), intent(in), optional :: weights(:) !! Nonnegative observation weights; absent means equal weights.
      logical, intent(in), optional :: boundary_repair !! Enable data-adaptive local-linear endpoint repair; default true.
      integer, intent(in), optional :: grid_size !! Number of regular fit cells; default 400 and at least 4.

      logical :: had_weights
      logical :: use_boundary_repair
      real(dp), allocatable :: boundary_x(:)
      real(dp), allocatable :: fitted(:, :)
      real(dp), allocatable :: grid_points(:)
      real(dp), allocatable :: influences(:)
      real(dp), allocatable :: values(:)
      real(dp), allocatable :: w(:)
      real(dp), allocatable :: work_x(:)
      real(dp) :: qnan_value
      integer :: interp_err

      ierr = 0
      qnan_value = quiet_nan()
      model%xmin = 0.0_dp
      model%xmax = 0.0_dp
      model%has_xmin = present(xmin)
      model%has_xmax = present(xmax)
      if (present(xmin)) model%xmin = xmin
      if (present(xmax)) model%xmax = xmax
      model%mult = 1.0_dp
      if (present(mult)) model%mult = mult
      model%deg = 2
      if (present(deg)) model%deg = deg
      model%grid_size = 400
      if (present(grid_size)) model%grid_size = grid_size
      model%boundary_repair = .true.
      if (present(boundary_repair)) model%boundary_repair = boundary_repair
      model%bw_automatic = .true.
      model%bw_spec = qnan_value
      if (present(bw)) then
         if (.not. ieee_is_nan(bw)) then
            model%bw_automatic = .false.
            model%bw_spec = bw
         end if
      end if
      model%bw = model%bw_spec
      model%var_type = parse_type(type_name)
      if (model%var_type == 0) then
         ierr = 1
         return
      end if

      if (size(x) == 0) then
         ierr = 10
         return
      end if
      if (model%mult <= 0.0_dp) then
         ierr = 11
         return
      end if
      if (.not. model%bw_automatic .and. model%bw_spec <= 0.0_dp) then
         ierr = 12
         return
      end if
      if (model%deg < 0 .or. model%deg > 2) then
         ierr = 13
         return
      end if
      if (model%grid_size < 4) then
         ierr = 14
         return
      end if
      if (model%has_xmin .and. model%has_xmax) then
         if (model%xmin > model%xmax) then
            ierr = 15
            return
         end if
      end if
      if (model%var_type == kde1d_discrete) then
         if (model%has_xmin .and. model%xmin /= floor(model%xmin)) then
            ierr = 16
            return
         end if
         if (model%has_xmax .and. model%xmax /= floor(model%xmax)) then
            ierr = 16
            return
         end if
      end if
      had_weights = present(weights)
      if (had_weights) then
         if (size(weights) /= size(x)) then
            ierr = 17
            return
         end if
         if (any((.not. ieee_is_nan(weights)) .and. weights < 0.0_dp)) then
            ierr = 18
            return
         end if
         if (any((.not. ieee_is_nan(weights)) .and. (.not. ieee_is_finite(weights)))) then
            ierr = 18
            return
         end if
      end if
      if (model%var_type == kde1d_discrete) then
         if (any((.not. ieee_is_nan(x)) .and. x /= floor(x))) then
            ierr = 19
            return
         end if
      end if
      if (.not. check_boundaries(x, model)) then
         ierr = 20
         return
      end if

      call preprocess_data(x, weights, work_x, w, had_weights, ierr)
      if (ierr /= 0) return
      model%nobs = count(.not. ieee_is_nan(x))

      if (model%var_type == kde1d_zero_inflated) then
         call preprocess_zero_inflated(work_x, w, model%prob0, had_weights, ierr)
         if (ierr /= 0) return
         if (size(work_x) == 0) then
            grid_points = [-2.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp]
            values = 0.0_dp * grid_points
            call model%grid%initialize(grid_points, values, 0, interp_err)
            if (interp_err /= 0) then
               ierr = 31
               return
            end if
            model%bw = qnan_value
            model%loglik = 0.0_dp
            model%edf = 1.0_dp
            model%fitted = .true.
            return
         end if
      else if (model%var_type == kde1d_discrete) then
         work_x = equi_jitter_codes(work_x)
      end if

      use_boundary_repair = model%boundary_repair .and. size(work_x) >= 16 &
         .and. (model%has_xmin .or. model%has_xmax)
      if (use_boundary_repair) boundary_x = work_x

      call set_one_sided_scale(model, work_x, w, had_weights)
      work_x = boundary_transform(model, work_x, .false.)

      if (model%bw_automatic) then
         model%bw = plugin_bandwidth(work_x, model%deg, merge_weights(w, had_weights))
      else
         model%bw = model%bw_spec
      end if
      model%bw = model%bw * model%mult
      if (model%var_type == kde1d_discrete .and. .not. model%has_xmin .and. .not. model%has_xmax) then
         model%bw = max(model%bw, 0.1_dp)
      end if
      if (.not. ieee_is_finite(model%bw) .or. model%bw <= 0.0_dp) then
         ierr = 21
         return
      end if

      grid_points = construct_grid_points(model, work_x)
      call fit_lp(model, work_x, boundary_transform(model, grid_points, .false.), &
         merge_weights(w, had_weights), fitted, ierr)
      if (ierr /= 0) return
      values = boundary_correct(model, grid_points, fitted(:, 1))
      influences = fitted(:, 2)

      if (.not. model%has_xmin .and. model%has_xmax) then
         grid_points = grid_points(size(grid_points):1:-1)
         influences = influences(size(influences):1:-1)
      end if
      call model%grid%initialize(grid_points, max(values, 0.0_dp), 3, interp_err)
      if (interp_err /= 0) then
         ierr = 22
         return
      end if

      if (use_boundary_repair) then
         call repair_boundaries(model, boundary_x, grid_points, influences, merge_weights(w, had_weights), ierr)
         if (ierr /= 0) return
      end if

      work_x = boundary_transform(model, work_x, .true.)
      if (model%var_type == kde1d_discrete) work_x = anint(work_x)
      if (.not. had_weights) then
         if (allocated(w)) deallocate(w)
         allocate(w(size(work_x)))
         w = 1.0_dp
         had_weights = .true.
      end if
      model%loglik = weighted_loglik(work_x, w, model)
      if (model%prob0 > 0.0_dp) then
         model%loglik = model%loglik + real(size(x), dp) * model%prob0 * log(model%prob0)
      end if

      influences = max(0.0_dp, min(3.0_dp, influences))
      call calculate_edf(model, grid_points, influences, work_x, ierr)
      if (ierr /= 0) return
      model%edf = model%edf + merge(1.0_dp, 0.0_dp, model%prob0 > 0.0_dp)
      model%bw = model%bw / model%mult
      model%fitted = .true.
   end subroutine kde1d_fit

   pure function dkde1d(x, model) result(values)
      real(dp), intent(in) :: x(:) !! Evaluation points for density or probability-mass values.
      type(kde1d_model), intent(in) :: model !! Fitted kde1d model.
      real(dp) :: values(size(x))

      select case (model%var_type)
      case (kde1d_discrete)
         values = pdf_discrete(x, model)
      case (kde1d_zero_inflated)
         values = pdf_zero_inflated(x, model)
      case default
         values = pdf_continuous(x, model)
      end select
   end function dkde1d

   pure function pkde1d(q, model) result(values)
      real(dp), intent(in) :: q(:) !! Quantile points at which the fitted CDF is evaluated.
      type(kde1d_model), intent(in) :: model !! Fitted kde1d model.
      real(dp) :: values(size(q))

      select case (model%var_type)
      case (kde1d_discrete)
         values = cdf_discrete(q, model)
      case (kde1d_zero_inflated)
         values = cdf_zero_inflated(q, model)
      case default
         values = model%grid%integrate(q, .true.)
      end select
   end function pkde1d

   pure function qkde1d(p, model) result(values)
      real(dp), intent(in) :: p(:) !! Probabilities in [0,1]; invalid finite probabilities return NaN.
      type(kde1d_model), intent(in) :: model !! Fitted kde1d model.
      real(dp) :: values(size(p))

      integer :: i

      select case (model%var_type)
      case (kde1d_discrete)
         values = quantile_discrete(p, model)
      case (kde1d_zero_inflated)
         values = quantile_zero_inflated(p, model)
      case default
         values = model%grid%quantile(max(0.0_dp, min(1.0_dp, p)))
      end select
      do i = 1, size(p)
         if ((.not. ieee_is_nan(p(i))) .and. (p(i) < 0.0_dp .or. p(i) > 1.0_dp)) values(i) = quiet_nan()
         if (ieee_is_nan(p(i))) values(i) = quiet_nan()
      end do
   end function qkde1d

   subroutine rkde1d(n, model, values, quasi, seed)
      integer, intent(in) :: n !! Number of deviates to generate; nonpositive values return an empty/unchanged prefix.
      type(kde1d_model), intent(in) :: model !! Fitted kde1d model used for the quantile transform.
      real(dp), intent(out) :: values(:) !! Output sample; caller supplies at least max(n,0) elements.
      logical, intent(in), optional :: quasi !! If true, use a one-dimensional Sobol/van-der-Corput sequence.
      integer, intent(in), optional :: seed !! Optional deterministic seed or digital shift for generated uniforms.

      logical :: use_quasi
      real(dp), allocatable :: u(:)
      integer :: i

      if (n <= 0) return
      allocate(u(n))
      use_quasi = .false.
      if (present(quasi)) use_quasi = quasi
      if (use_quasi) then
         do i = 1, n
            u(i) = sobol_1d(i, seed)
         end do
      else
         call seed_fortran_rng(seed)
         call random_number(u)
      end if
      values(:n) = qkde1d(u, model)
   end subroutine rkde1d

   pure real(dp) function kde1d_loglik(model) result(value)
      type(kde1d_model), intent(in) :: model !! Fitted kde1d model whose weighted log-likelihood is requested.

      value = model%loglik
   end function kde1d_loglik

   pure function kde1d_summary(model) result(values)
      type(kde1d_model), intent(in) :: model !! Fitted kde1d model summarized as nobs,bw,mult,loglik,edf.
      real(dp) :: values(5)

      values = [real(model%nobs, dp), model%bw, model%mult, model%loglik, model%edf]
   end function kde1d_summary

   pure integer function parse_type(type_name) result(var_type)
      character(len=*), intent(in), optional :: type_name !! Optional R-compatible variable-type alias.

      character(len=:), allocatable :: name

      if (.not. present(type_name)) then
         var_type = kde1d_continuous
         return
      end if
      name = trim(adjustl(type_name))
      select case (name)
      case ('c', 'cont', 'continuous')
         var_type = kde1d_continuous
      case ('d', 'disc', 'discrete')
         var_type = kde1d_discrete
      case ('zi', 'zinfl', 'zero-inflated', 'zero_inflated')
         var_type = kde1d_zero_inflated
      case default
         var_type = 0
      end select
   end function parse_type

   pure logical function check_boundaries(x, model) result(ok)
      real(dp), intent(in) :: x(:) !! Raw observations checked against model support restrictions.
      type(kde1d_model), intent(in) :: model !! Model carrying support bounds and zero-inflated semantics.

      integer :: i

      ok = .true.
      do i = 1, size(x)
         if (ieee_is_nan(x(i))) cycle
         if (model%var_type == kde1d_zero_inflated .and. x(i) == 0.0_dp) cycle
         if (model%has_xmin .and. x(i) < model%xmin) ok = .false.
         if (model%has_xmax .and. x(i) > model%xmax) ok = .false.
         if (.not. ok) return
      end do
   end function check_boundaries

   pure subroutine preprocess_data(x, weights, work_x, w, had_weights, ierr)
      real(dp), intent(in) :: x(:) !! Raw observations before removal of NaN/zero-weight rows.
      real(dp), intent(in), optional :: weights(:) !! Optional weights aligned with x.
      real(dp), allocatable, intent(out) :: work_x(:) !! Retained observations after drop-marker filtering.
      real(dp), allocatable, intent(out) :: w(:) !! Retained weights normalized to mean one, or empty if absent.
      logical, intent(in) :: had_weights !! Whether weights were supplied by the caller.
      integer, intent(out) :: ierr !! Zero on success; nonzero if no observation remains.

      logical, allocatable :: keep(:)
      integer :: n

      allocate(keep(size(x)))
      keep = .not. ieee_is_nan(x)
      if (had_weights) then
         keep = keep .and. (.not. ieee_is_nan(weights)) .and. weights /= 0.0_dp
      end if
      n = count(keep)
      if (n == 0) then
         ierr = 1
         allocate(work_x(0), w(0))
         return
      end if
      allocate(work_x(n))
      work_x = pack(x, keep)
      if (had_weights) then
         allocate(w(n))
         w = pack(weights, keep)
         w = w / (sum(w) / real(n, dp))
      else
         allocate(w(0))
      end if
      ierr = 0
   end subroutine preprocess_data

   pure subroutine preprocess_zero_inflated(x, w, prob0, had_weights, ierr)
      real(dp), allocatable, intent(inout) :: x(:) !! Retained observations, with zeros removed on output.
      real(dp), allocatable, intent(inout) :: w(:) !! Normalized weights, retained for nonzero observations on output.
      real(dp), intent(out) :: prob0 !! Estimated weighted point mass at zero.
      logical, intent(inout) :: had_weights !! Updated true because zero-inflation estimation requires explicit weights.
      integer, intent(out) :: ierr !! Zero on success; reserved for allocation/numerical failures.

      logical, allocatable :: keep(:)
      real(dp), allocatable :: w_all(:)

      ierr = 0
      allocate(w_all(size(x)))
      if (size(w) == size(x)) then
         w_all = w
      else
         w_all = 1.0_dp
      end if
      where (x == 0.0_dp) w_all = 0.0_dp
      prob0 = 1.0_dp - sum(w_all) / real(size(w_all), dp)
      allocate(keep(size(x)))
      keep = w_all /= 0.0_dp
      x = pack(x, keep)
      w = pack(w_all, keep)
      had_weights = .true.
   end subroutine preprocess_zero_inflated

   pure function merge_weights(weights, use_weights) result(w)
      real(dp), intent(in) :: weights(:) !! Candidate observation weights.
      logical, intent(in) :: use_weights !! True to pass the weights, false to return a zero-length vector.
      real(dp), allocatable :: w(:)

      if (use_weights) then
         w = weights
      else
         allocate(w(0))
      end if
   end function merge_weights

   pure real(dp) function effective_bound(model, lower) result(value)
      type(kde1d_model), intent(in) :: model !! Model carrying support bounds and variable type.
      logical, intent(in) :: lower !! True requests the lower effective bound; false requests upper.

      if (lower) then
         if (.not. model%has_xmin) then
            value = quiet_nan()
         else if (model%var_type == kde1d_discrete) then
            value = model%xmin - 0.5_dp
         else
            value = model%xmin
         end if
      else
         if (.not. model%has_xmax) then
            value = quiet_nan()
         else if (model%var_type == kde1d_discrete) then
            value = model%xmax + 0.5_dp
         else
            value = model%xmax
         end if
      end if
   end function effective_bound

   pure subroutine set_one_sided_scale(model, x, weights, has_weights)
      type(kde1d_model), intent(inout) :: model !! Model whose one-sided transform scale and offset are initialized.
      real(dp), intent(in) :: x(:) !! Preprocessed observations on the original scale.
      real(dp), intent(in) :: weights(:) !! Optional normalized case weights.
      logical, intent(in) :: has_weights !! Whether weights should be used for the robust median distance.

      real(dp), allocatable :: distances(:)
      real(dp) :: bound

      if (model%has_xmin .eqv. model%has_xmax) return
      allocate(distances(size(x)))
      if (model%has_xmin) then
         bound = effective_bound(model, .true.)
         distances = x - bound
      else
         bound = effective_bound(model, .false.)
         distances = bound - x
      end if
      if (has_weights .and. size(weights) == size(x)) then
         model%boundary_scale = weighted_quantile_cpp(distances, 0.5_dp, weights)
      else
         model%boundary_scale = r_median(distances)
      end if
      if (.not. (model%boundary_scale > 0.0_dp)) model%boundary_scale = maxval(distances)
      if (.not. (model%boundary_scale > 0.0_dp)) model%boundary_scale = 1.0_dp
      model%boundary_offset = 1.0e-5_dp * model%boundary_scale
   end subroutine set_one_sided_scale

   pure function boundary_transform(model, x, inverse) result(y)
      type(kde1d_model), intent(in) :: model !! Model defining finite support and one-sided transform scaling.
      real(dp), intent(in) :: x(:) !! Values on the original or transformed scale.
      logical, intent(in) :: inverse !! False applies the support transform; true applies its inverse.
      real(dp) :: y(size(x))

      real(dp) :: bound
      real(dp) :: rng

      y = x
      if (.not. inverse) then
         if (model%has_xmin .and. model%has_xmax) then
            rng = effective_bound(model, .false.) - effective_bound(model, .true.)
            y = (x - effective_bound(model, .true.) + 5.0e-5_dp * rng) / (1.0001_dp * rng)
            y = r_qnorm(y)
         else if (model%has_xmin) then
            bound = effective_bound(model, .true.)
            y = 4.0_dp * (((model%boundary_offset + x - bound) / model%boundary_scale)**0.25_dp &
               - (model%boundary_offset / model%boundary_scale)**0.25_dp)
         else if (model%has_xmax) then
            bound = effective_bound(model, .false.)
            y = 4.0_dp * (((model%boundary_offset + bound - x) / model%boundary_scale)**0.25_dp &
               - (model%boundary_offset / model%boundary_scale)**0.25_dp)
         end if
      else
         if (model%has_xmin .and. model%has_xmax) then
            rng = effective_bound(model, .false.) - effective_bound(model, .true.)
            y = r_pnorm(x) * 1.0001_dp * rng + effective_bound(model, .true.) - 5.0e-5_dp * rng
         else if (model%has_xmin) then
            bound = effective_bound(model, .true.)
            y = model%boundary_scale * (x / 4.0_dp &
               + (model%boundary_offset / model%boundary_scale)**0.25_dp)**4 &
               + bound - model%boundary_offset
         else if (model%has_xmax) then
            bound = effective_bound(model, .false.)
            y = bound + model%boundary_offset - model%boundary_scale * (x / 4.0_dp &
               + (model%boundary_offset / model%boundary_scale)**0.25_dp)**4
         end if
      end if
   end function boundary_transform

   pure function boundary_correct(model, x, fhat) result(values)
      type(kde1d_model), intent(in) :: model !! Model defining the support transform whose Jacobian is undone.
      real(dp), intent(in) :: x(:) !! Original-scale grid points.
      real(dp), intent(in) :: fhat(:) !! Density estimate evaluated on the transformed scale.
      real(dp) :: values(size(fhat))

      real(dp) :: corr(size(fhat))
      real(dp) :: rng

      if (model%has_xmin .and. model%has_xmax) then
         rng = effective_bound(model, .false.) - effective_bound(model, .true.)
         corr = (x - effective_bound(model, .true.) + 5.0e-5_dp * rng) &
            / (effective_bound(model, .false.) - effective_bound(model, .true.) + 1.0e-4_dp * rng)
         corr = r_dnorm(r_qnorm(corr))
         corr = 1.0_dp / (corr * (effective_bound(model, .false.) &
            - effective_bound(model, .true.) + 1.0e-4_dp * rng))
      else if (model%has_xmin) then
         corr = ((model%boundary_offset + x - effective_bound(model, .true.)) / model%boundary_scale)**(-0.75_dp) &
            / model%boundary_scale
      else if (model%has_xmax) then
         corr = ((model%boundary_offset + effective_bound(model, .false.) - x) / model%boundary_scale)**(-0.75_dp) &
            / model%boundary_scale
      else
         corr = 1.0_dp
      end if
      values = fhat * corr
      if (.not. model%has_xmin .and. model%has_xmax) values = values(size(values):1:-1)
   end function boundary_correct

   pure function construct_grid_points(model, x_transformed) result(grid_points)
      type(kde1d_model), intent(in) :: model !! Model defining bandwidth, support transform, and requested grid size.
      real(dp), intent(in) :: x_transformed(:) !! Preprocessed observations already on transformed scale.
      real(dp), allocatable :: grid_points(:)

      real(dp), allocatable :: zgrid(:)
      real(dp) :: lo
      real(dp) :: hi
      real(dp) :: boundaries(2)
      integer :: i

      lo = minval(x_transformed)
      hi = maxval(x_transformed)
      if (model%var_type == kde1d_discrete .and. .not. model%has_xmin .and. .not. model%has_xmax) then
         continue
      else if (.not. model%has_xmin .and. .not. model%has_xmax) then
         lo = lo - 4.0_dp * model%bw
         hi = hi + 4.0_dp * model%bw
      else if (model%has_xmin .and. model%has_xmax) then
         boundaries = [effective_bound(model, .true.), effective_bound(model, .false.)]
         boundaries = boundary_transform(model, boundaries, .false.)
         lo = boundaries(1)
         hi = boundaries(2)
      else
         if (model%has_xmin) then
            boundaries(1) = effective_bound(model, .true.)
         else
            boundaries(1) = effective_bound(model, .false.)
         end if
         boundaries(1:1) = boundary_transform(model, boundaries(1:1), .false.)
         lo = boundaries(1)
         hi = hi + 4.0_dp * model%bw
      end if

      allocate(zgrid(model%grid_size + 1))
      do i = 1, model%grid_size + 1
         zgrid(i) = lo + real(i - 1, dp) * (hi - lo) / real(model%grid_size, dp)
      end do
      grid_points = boundary_transform(model, zgrid, .true.)
      if (model%has_xmin) grid_points(1) = effective_bound(model, .true.)
      if (model%has_xmax) then
         if (model%has_xmin) then
            grid_points(size(grid_points)) = effective_bound(model, .false.)
         else
            grid_points(1) = effective_bound(model, .false.)
         end if
      end if
   end function construct_grid_points

   pure subroutine fit_lp(model, x, grid, weights, result, ierr)
      type(kde1d_model), intent(in) :: model !! Model carrying bandwidth and local-polynomial degree.
      real(dp), intent(in) :: x(:) !! Transformed observations.
      real(dp), intent(in) :: grid(:) !! Regular transformed-domain evaluation knots.
      real(dp), intent(in) :: weights(:) !! Optional normalized case weights; zero length means equal weights.
      real(dp), allocatable, intent(out) :: result(:, :) !! Two columns: fitted density and influence values.
      integer, intent(out) :: ierr !! Zero on success; nonzero if the transformed grid is invalid.

      real(dp), allocatable :: count(:)
      real(dp), allocatable :: f0(:)
      real(dp), allocatable :: f1(:)
      real(dp), allocatable :: f2(:)
      real(dp), allocatable :: wcount(:)
      real(dp), allocatable :: wbin(:)
      real(dp), allocatable :: b(:)
      real(dp), allocatable :: s(:)
      real(dp) :: density_floor
      real(dp) :: d
      real(dp) :: r
      integer :: k
      integer :: m

      ierr = 0
      m = size(grid)
      allocate(result(m, 2), f0(m), f1(m), f2(m), wbin(m), s(m), b(m), wcount(m))
      call kde_drv_binned(x, model%bw, grid(1), grid(m), m - 1, weights, 0, f0, wcount)
      density_floor = epsilon(1.0_dp) * maxval(abs(f0))
      if (size(weights) == size(x)) then
         allocate(count(m))
         call linbin(x, grid(1), grid(m), m - 1, [(1.0_dp, k = 1, size(x))], count)
         wbin = 0.0_dp
         do k = 1, m
            if (count(k) > 0.0_dp) wbin(k) = wcount(k) / count(k)
         end do
      else
         wbin = 1.0_dp
      end if

      result(:, 1) = f0
      result(:, 2) = 0.0_dp
      do k = 1, m
         if (f0(k) > density_floor) then
            result(k, 2) = k0 / (real(size(x), dp) * model%bw) * wbin(k) / f0(k)
         end if
      end do
      if (model%deg == 0) return

      call kde_drv_binned(x, model%bw, grid(1), grid(m), m - 1, weights, 1, f1)
      s = model%bw
      b = 0.0_dp
      where (f0 > density_floor) b = f1 / f0
      if (model%deg == 2) then
         call kde_drv_binned(x, model%bw, grid(1), grid(m), m - 1, weights, 2, f2)
         do k = 1, m
            if (f0(k) <= density_floor) cycle
            d = f2(k) / f0(k) - b(k)**2
            if (1.0_dp + model%bw**2 * d <= 0.0_dp) then
               result(k, :) = 0.0_dp
               cycle
            end if
            r = 1.0_dp / sqrt(1.0_dp + model%bw**2 * d)
            s(k) = (r / model%bw)**2
            b(k) = b(k) * model%bw**2
            result(k, 1) = model%bw * sqrt(s(k)) * result(k, 1)
         end do
      end if
      result(:, 1) = result(:, 1) * exp(-0.5_dp * b * b * s)

      do k = 1, m
         if (.not. ieee_is_finite(f0(k)) .or. f0(k) <= density_floor) then
            result(k, :) = 0.0_dp
            cycle
         end if
         result(k, 2) = calculate_influence(model%deg, size(x), f0(k), f1(k), f2(k), &
            model%bw, s(k), wbin(k))
         if (.not. ieee_is_finite(result(k, 1)) .or. result(k, 1) < 0.0_dp) result(k, :) = 0.0_dp
      end do
   end subroutine fit_lp

   pure real(dp) function calculate_influence(degree, n, f0, f1, f2, bandwidth, s, weight) result(value)
      integer, intent(in) :: degree !! Local-polynomial degree 0, 1, or 2.
      integer, intent(in) :: n !! Number of retained observations.
      real(dp), intent(in) :: f0 !! Zeroth derivative pilot density at the grid knot.
      real(dp), intent(in) :: f1 !! First derivative pilot density at the grid knot.
      real(dp), intent(in) :: f2 !! Second derivative pilot density at the grid knot.
      real(dp), intent(in) :: bandwidth !! Positive fitting bandwidth.
      real(dp), intent(in) :: s !! Local precision term used by log-quadratic fitting.
      real(dp), intent(in) :: weight !! Average normalized case weight in the grid bin.

      real(dp) :: b2
      real(dp) :: det
      real(dp) :: minv00
      real(dp) :: m(3, 3)
      real(dp) :: s2

      b2 = bandwidth**2
      select case (degree)
      case (0)
         minv00 = 1.0_dp / f0
      case (1)
         m = 0.0_dp
         m(1, 1) = f0
         m(1, 2) = b2 * f1
         m(2, 1) = m(1, 2)
         m(2, 2) = f0 * b2 + b2 * f1 * f1 * b2 / f0
         det = m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1)
         if (abs(det) <= tiny(1.0_dp)) then
            value = 0.0_dp
            return
         end if
         minv00 = m(2, 2) / det
      case default
         m = 0.0_dp
         m(1, 1) = f0
         m(1, 2) = b2 * f1
         m(2, 1) = m(1, 2)
         m(2, 2) = b2 * f2 * b2 + b2 * f0
         m(3, 1) = m(2, 2) / 2.0_dp
         m(1, 3) = m(3, 1)
         s2 = b2 * f1 / f0
         m(2, 3) = f0 / 2.0_dp * (3.0_dp / s * s2 + s2**3)
         m(3, 2) = m(2, 3)
         m(3, 3) = f0 / 4.0_dp * (3.0_dp / s**2 + 6.0_dp / s * s2**2 + s2**4)
         det = determinant3(m)
         if (abs(det) <= tiny(1.0_dp)) then
            value = 0.0_dp
            return
         end if
         minv00 = (m(2, 2) * m(3, 3) - m(2, 3) * m(3, 2)) / det
      end select
      value = k0 * weight / (real(n, dp) * bandwidth) * minv00
   end function calculate_influence

   pure real(dp) function determinant3(a) result(value)
      real(dp), intent(in) :: a(3, 3) !! Three-by-three matrix whose determinant is required.

      value = a(1, 1) * (a(2, 2) * a(3, 3) - a(2, 3) * a(3, 2)) &
         - a(1, 2) * (a(2, 1) * a(3, 3) - a(2, 3) * a(3, 1)) &
         + a(1, 3) * (a(2, 1) * a(3, 2) - a(2, 2) * a(3, 1))
   end function determinant3

   pure function pdf_continuous(x, model) result(values)
      real(dp), intent(in) :: x(:) !! Continuous-scale density evaluation points.
      type(kde1d_model), intent(in) :: model !! Fitted model supplying interpolation grid and support bounds.
      real(dp) :: values(size(x))

      values = max(model%grid%interpolate(x), 0.0_dp)
      if (model%has_xmin) where (x < model%xmin) values = 0.0_dp
      if (model%has_xmax) where (x > model%xmax) values = 0.0_dp
   end function pdf_continuous

   pure function pdf_discrete(x, model) result(values)
      real(dp), intent(in) :: x(:) !! Candidate integer support levels; nonintegers receive zero mass.
      type(kde1d_model), intent(in) :: model !! Fitted discrete kde1d model.
      real(dp) :: values(size(x))

      real(dp), allocatable :: levels(:)
      real(dp) :: lb
      real(dp) :: normalizer
      real(dp) :: ub
      integer :: i
      integer :: nlevels

      lb = discrete_support_min(model)
      ub = discrete_support_max(model)
      nlevels = max(1, nint(ub - lb) + 1)
      allocate(levels(nlevels))
      do i = 1, nlevels
         levels(i) = lb + real(i - 1, dp)
      end do
      values = pdf_continuous(x, model)
      where (x < lb .or. x > ub .or. x /= anint(x)) values = 0.0_dp
      normalizer = sum(pdf_continuous(levels, model))
      if (normalizer > 0.0_dp .and. ieee_is_finite(normalizer)) then
         values = values / normalizer
      else
         values = quiet_nan()
      end if
   end function pdf_discrete

   pure function pdf_zero_inflated(x, model) result(values)
      real(dp), intent(in) :: x(:) !! Evaluation points for the hurdle density/point mass.
      type(kde1d_model), intent(in) :: model !! Fitted zero-inflated kde1d model.
      real(dp) :: values(size(x))

      values = (1.0_dp - model%prob0) * pdf_continuous(x, model)
      where (x == 0.0_dp) values = model%prob0
   end function pdf_zero_inflated

   pure function cdf_discrete(x, model) result(values)
      real(dp), intent(in) :: x(:) !! Points at which the discrete fitted CDF is evaluated.
      type(kde1d_model), intent(in) :: model !! Fitted discrete kde1d model.
      real(dp) :: values(size(x))

      real(dp), allocatable :: cumulative(:)
      real(dp), allocatable :: levels(:)
      real(dp) :: lb
      real(dp) :: ub
      integer :: i
      integer :: idx
      integer :: nlevels

      lb = discrete_support_min(model)
      ub = discrete_support_max(model)
      nlevels = max(1, nint(ub - lb) + 1)
      allocate(levels(nlevels), cumulative(nlevels))
      do i = 1, nlevels
         levels(i) = lb + real(i - 1, dp)
      end do
      cumulative = pdf_discrete(levels, model)
      do i = 2, nlevels
         cumulative(i) = cumulative(i) + cumulative(i - 1)
      end do
      do i = 1, size(x)
         if (ieee_is_nan(x(i))) then
            values(i) = quiet_nan()
         else if (x(i) < lb) then
            values(i) = 0.0_dp
         else if (x(i) >= ub) then
            values(i) = 1.0_dp
         else
            idx = int(x(i) - lb) + 1
            values(i) = min(1.0_dp, max(0.0_dp, cumulative(idx)))
         end if
      end do
   end function cdf_discrete

   pure function cdf_zero_inflated(x, model) result(values)
      real(dp), intent(in) :: x(:) !! Points at which the zero-inflated fitted CDF is evaluated.
      type(kde1d_model), intent(in) :: model !! Fitted zero-inflated kde1d model.
      real(dp) :: values(size(x))

      real(dp) :: continuous(size(x))

      if (model%prob0 < 1.0_dp) then
         continuous = model%grid%integrate(x, .true.)
      else
         continuous = 0.0_dp
      end if
      values = (1.0_dp - model%prob0) * continuous
      where (x >= 0.0_dp) values = values + model%prob0
   end function cdf_zero_inflated

   pure function quantile_discrete(p, model) result(values)
      real(dp), intent(in) :: p(:) !! Probabilities whose fitted discrete quantiles are requested.
      type(kde1d_model), intent(in) :: model !! Fitted discrete kde1d model.
      real(dp) :: values(size(p))

      real(dp), allocatable :: cumulative(:)
      real(dp), allocatable :: levels(:)
      real(dp) :: lb
      real(dp) :: ub
      integer :: i
      integer :: j
      integer :: nlevels

      lb = discrete_support_min(model)
      ub = discrete_support_max(model)
      nlevels = max(1, nint(ub - lb) + 1)
      allocate(levels(nlevels), cumulative(nlevels))
      do i = 1, nlevels
         levels(i) = lb + real(i - 1, dp)
      end do
      cumulative = cdf_discrete(levels, model)
      do i = 1, size(p)
         if (ieee_is_nan(p(i)) .or. p(i) < 0.0_dp .or. p(i) > 1.0_dp) then
            values(i) = quiet_nan()
            cycle
         end if
         values(i) = levels(nlevels)
         do j = 1, nlevels
            if (cumulative(j) >= p(i)) then
               values(i) = levels(j)
               exit
            end if
         end do
      end do
   end function quantile_discrete

   pure function quantile_zero_inflated(p, model) result(values)
      real(dp), intent(in) :: p(:) !! Probabilities whose hurdle-model quantiles are requested.
      type(kde1d_model), intent(in) :: model !! Fitted zero-inflated kde1d model.
      real(dp) :: values(size(p))

      real(dp) :: newp(size(p))
      real(dp) :: p0
      integer :: i

      if (model%prob0 >= 1.0_dp) then
         values = 0.0_dp
         do i = 1, size(p)
            if (ieee_is_nan(p(i)) .or. p(i) < 0.0_dp .or. p(i) > 1.0_dp) values(i) = quiet_nan()
         end do
         return
      end if
      block
         real(dp) :: c0(1)
         c0 = cdf_zero_inflated([0.0_dp], model)
         p0 = c0(1)
      end block
      do i = 1, size(p)
         if (p(i) <= p0 - model%prob0) then
            newp(i) = p(i) / (1.0_dp - model%prob0)
         else
            newp(i) = max(p(i) - model%prob0, 0.0_dp) / (1.0_dp - model%prob0)
         end if
      end do
      values = model%grid%quantile(max(0.0_dp, min(1.0_dp, newp)))
      do i = 1, size(p)
         if (p(i) > p0 - model%prob0 .and. p(i) <= p0) values(i) = 0.0_dp
         if (ieee_is_nan(p(i)) .or. p(i) < 0.0_dp .or. p(i) > 1.0_dp) values(i) = quiet_nan()
      end do
   end function quantile_zero_inflated

   pure real(dp) function discrete_support_min(model) result(value)
      type(kde1d_model), intent(in) :: model !! Fitted discrete model whose integer support lower end is requested.

      if (model%has_xmin) then
         value = model%xmin
      else
         value = floor(model%grid%grid_points(1))
      end if
   end function discrete_support_min

   pure real(dp) function discrete_support_max(model) result(value)
      type(kde1d_model), intent(in) :: model !! Fitted discrete model whose integer support upper end is requested.

      if (model%has_xmax) then
         value = model%xmax
      else
         value = max(discrete_support_min(model), &
            real(ceiling(model%grid%grid_points(size(model%grid%grid_points))), dp))
      end if
   end function discrete_support_max

   pure real(dp) function weighted_loglik(x, weights, model) result(value)
      real(dp), intent(in) :: x(:) !! Retained original-scale observations for likelihood evaluation.
      real(dp), intent(in) :: weights(:) !! Mean-one likelihood weights aligned with x.
      type(kde1d_model), intent(in) :: model !! Fitted model used to evaluate observation densities.

      real(dp) :: dens(size(x))

      dens = dkde1d(x, model)
      if (any(dens <= 0.0_dp)) then
         value = -huge(1.0_dp)
      else
         value = sum(weights * log(dens))
      end if
   end function weighted_loglik

   pure subroutine calculate_edf(model, grid_points, influences, x, ierr)
      type(kde1d_model), intent(inout) :: model !! Fitted model receiving the effective-degrees-of-freedom estimate.
      real(dp), intent(in) :: grid_points(:) !! Original-scale fit grid carrying influence ordinates.
      real(dp), intent(in) :: influences(:) !! Truncated diagonal-influence estimates on grid_points.
      real(dp), intent(in) :: x(:) !! Retained original-scale observations at which influences are summed.
      integer, intent(out) :: ierr !! Zero on success; nonzero if influence interpolation cannot be initialized.

      type(interpolation_grid) :: infl_grid
      real(dp) :: infl(size(x))

      call infl_grid%initialize(grid_points, influences, 0, ierr)
      if (ierr /= 0) return
      infl = infl_grid%interpolate(x)
      model%edf = sum(infl)
   end subroutine calculate_edf

   pure subroutine repair_boundaries(model, x, grid, influences, weights, ierr)
      type(kde1d_model), intent(inout) :: model !! Fitted model whose transformed bulk fit may be repaired at endpoints.
      real(dp), intent(in) :: x(:) !! Pre-transform observations on the original jitter/continuous scale.
      real(dp), intent(in) :: grid(:) !! Increasing original-scale evaluation grid.
      real(dp), intent(inout) :: influences(:) !! Bulk influences replaced by fused influences on output.
      real(dp), intent(in) :: weights(:) !! Mean-one case weights or zero length for equal weighting.
      integer, intent(out) :: ierr !! Zero on success; nonzero only if the fused spline cannot be normalized.

      real(dp), allocatable :: case_weights(:)
      real(dp), allocatable :: lower_dist(:)
      real(dp), allocatable :: lower_weights(:)
      real(dp), allocatable :: upper_dist(:)
      real(dp), allocatable :: upper_weights(:)
      real(dp), allocatable :: w_lower(:)
      real(dp), allocatable :: w_upper(:)
      real(dp), allocatable :: lower_density(:)
      real(dp), allocatable :: lower_infl(:)
      real(dp), allocatable :: upper_density(:)
      real(dp), allocatable :: upper_infl(:)
      real(dp), allocatable :: bulk_cdf(:)
      real(dp) :: bandwidth_fraction
      real(dp) :: h
      real(dp) :: n_eff
      real(dp) :: q
      logical :: lower_finite
      logical :: upper_finite
      integer :: i
      integer :: n_lower
      integer :: n_upper

      ierr = 0
      if (size(weights) == size(x)) then
         case_weights = weights
      else
         allocate(case_weights(size(x)))
         case_weights = 1.0_dp
      end if
      n_eff = sum(case_weights)**2 / sum(case_weights**2)
      if (n_eff < 16.0_dp) return

      lower_finite = .false.
      upper_finite = .false.
      if (model%has_xmin) then
         call prepare_endpoint(model, x, case_weights, .true., lower_dist, lower_weights, lower_finite)
      end if
      if (model%has_xmax) then
         call prepare_endpoint(model, x, case_weights, .false., upper_dist, upper_weights, upper_finite)
      end if
      if (.not. lower_finite .and. .not. upper_finite) return

      q = min(0.25_dp, 1.0_dp / sqrt(n_eff))
      bulk_cdf = model%grid%integrate(grid, .true.)
      allocate(w_lower(size(grid)), w_upper(size(grid)))
      w_lower = 0.0_dp
      w_upper = 0.0_dp
      do i = 1, size(grid)
         if (lower_finite) w_lower(i) = endpoint_weight(bulk_cdf(i), q)
         if (upper_finite) w_upper(i) = endpoint_weight(1.0_dp - bulk_cdf(i), q)
      end do
      n_lower = count(w_lower > 0.0_dp)
      n_upper = count(w_upper > 0.0_dp)

      bandwidth_fraction = merge(1.0_dp, 0.75_dp, model%has_xmin .and. model%has_xmax)
      if (lower_finite) then
         h = select_boundary_bandwidth(lower_dist, lower_weights, bandwidth_fraction, model%mult)
      else
         h = select_boundary_bandwidth(upper_dist, upper_weights, bandwidth_fraction, model%mult)
      end if
      if (.not. (h > 0.0_dp) .or. .not. ieee_is_finite(h)) return

      allocate(lower_density(0), lower_infl(0), upper_density(0), upper_infl(0))
      if (lower_finite .and. n_lower > 0) then
         call fit_boundary_component(lower_dist, grid(:n_lower) - effective_bound(model, .true.), h, lower_weights, &
            lower_density, lower_infl)
      end if
      if (upper_finite .and. n_upper > 0) then
         block
            real(dp), allocatable :: upper_eval(:)
            upper_eval = effective_bound(model, .false.) - grid(size(grid) - n_upper + 1:size(grid))
            upper_eval = upper_eval(size(upper_eval):1:-1)
            call fit_boundary_component(upper_dist, upper_eval, h, upper_weights, upper_density, upper_infl)
            upper_density = upper_density(size(upper_density):1:-1)
            upper_infl = upper_infl(size(upper_infl):1:-1)
         end block
      end if
      call fuse_boundary_components(model, grid, w_lower, w_upper, lower_density, lower_infl, &
         upper_density, upper_infl, influences, ierr)
   end subroutine repair_boundaries

   pure real(dp) function endpoint_weight(probability, q) result(value)
      real(dp), intent(in) :: probability !! Bulk-tail probability determining endpoint fusion strength.
      real(dp), intent(in) :: q !! Shrinking endpoint probability width.

      real(dp) :: z

      z = min(1.0_dp, probability / q)
      value = (1.0_dp - z * z)**2
   end function endpoint_weight

   pure subroutine prepare_endpoint(model, x, weights, lower, dist, sorted_weights, finite_endpoint)
      type(kde1d_model), intent(in) :: model !! Model supplying the chosen finite support endpoint.
      real(dp), intent(in) :: x(:) !! Original-scale observations measured from the endpoint.
      real(dp), intent(in) :: weights(:) !! Mean-one case weights aligned with x.
      logical, intent(in) :: lower !! True measures distance from lower endpoint; false from upper endpoint.
      real(dp), allocatable, intent(out) :: dist(:) !! Ascending endpoint distances.
      real(dp), allocatable, intent(out) :: sorted_weights(:) !! Weights permuted into endpoint-distance order.
      logical, intent(out) :: finite_endpoint !! Classification that the endpoint has finite nonzero limiting density.

      if (lower) then
         dist = x - effective_bound(model, .true.)
      else
         dist = effective_bound(model, .false.) - x
      end if
      sorted_weights = weights
      call sort_pairs_local(dist, sorted_weights)
      finite_endpoint = is_finite_endpoint(dist, sorted_weights)
   end subroutine prepare_endpoint

   pure logical function is_finite_endpoint(dist, weights) result(finite_endpoint)
      real(dp), intent(in) :: dist(:) !! Ascending nonnegative distances from a finite support endpoint.
      real(dp), intent(in) :: weights(:) !! Correspondingly ordered mean-one case weights.

      real(dp), allocatable :: counts(:)
      real(dp) :: beta
      real(dp) :: beta_min
      real(dp) :: cumulative
      real(dp) :: denom
      real(dp) :: dist_k1
      real(dp) :: dist_min
      real(dp) :: n_eff
      real(dp) :: tail_mass
      real(dp) :: target
      integer :: i
      integer :: k1

      finite_endpoint = .false.
      n_eff = sum(weights)**2 / sum(weights**2)
      if (.not. (n_eff > 1.0_dp) .or. .not. ieee_is_finite(n_eff)) return
      counts = weights * n_eff / sum(weights)
      target = min(nearest(n_eff, -1.0_dp), real(ceiling(2.0_dp * sqrt(n_eff)), dp))
      k1 = 1
      cumulative = counts(1)
      do while (cumulative <= target .and. k1 < size(dist))
         k1 = k1 + 1
         cumulative = cumulative + counts(k1)
      end do
      dist_k1 = dist(k1)
      if (.not. (dist_k1 > 0.0_dp) .or. .not. ieee_is_finite(dist_k1)) return

      tail_mass = 0.0_dp
      denom = 0.0_dp
      dist_min = epsilon(1.0_dp) * dist_k1
      do i = 1, k1 - 1
         tail_mass = tail_mass + counts(i)
         denom = denom + counts(i) * log(dist_k1 / max(dist(i), dist_min))
      end do
      if (.not. (tail_mass > 0.0_dp) .or. .not. (denom > 0.0_dp) .or. .not. ieee_is_finite(denom)) return
      beta = tail_mass / denom
      beta_min = merge(0.9_dp, 0.975_dp, maxval(weights) == minval(weights))
      finite_endpoint = beta >= beta_min &
         .and. beta * (1.0_dp - 1.6448536269514722_dp / sqrt(tail_mass)) <= 1.0_dp
   end function is_finite_endpoint

   pure real(dp) function select_boundary_bandwidth(dist, weights, fraction, multiplier) result(h)
      real(dp), intent(in) :: dist(:) !! Ascending endpoint distances.
      real(dp), intent(in) :: weights(:) !! Mean-one weights aligned with endpoint distances.
      real(dp), intent(in) :: fraction !! Fraction of positive-weight rows used for bandwidth selection.
      real(dp), intent(in) :: multiplier !! Public bandwidth multiplier applied to the selected boundary bandwidth.

      real(dp), allocatable :: d(:)
      real(dp), allocatable :: w(:)
      integer :: i
      integer :: j
      integer :: n_bw
      integer :: n_positive

      n_positive = count(weights > 0.0_dp)
      if (n_positive == 0) then
         h = quiet_nan()
         return
      end if
      allocate(d(n_positive), w(n_positive))
      j = 0
      do i = 1, size(dist)
         if (weights(i) > 0.0_dp) then
            j = j + 1
            d(j) = dist(i)
            w(j) = weights(i)
         end if
      end do
      n_bw = min(n_positive, max(4, ceiling(fraction * real(n_positive, dp))))
      h = plugin_bandwidth(d(:n_bw), 2, w(:n_bw)) * multiplier
   end function select_boundary_bandwidth

   pure subroutine fit_boundary_component(dist, eval_dist, bandwidth, weights, density, influence_num)
      real(dp), intent(in) :: dist(:) !! Ascending endpoint distances for retained observations.
      real(dp), intent(in) :: eval_dist(:) !! Nonnegative endpoint distances at which the component is evaluated.
      real(dp), intent(in) :: bandwidth !! Positive local-linear boundary bandwidth.
      real(dp), intent(in) :: weights(:) !! Mean-one case weights aligned with dist.
      real(dp), allocatable, intent(out) :: density(:) !! Nonnegative endpoint-component density estimates.
      real(dp), allocatable, intent(out) :: influence_num(:) !! Endpoint influence numerators on eval_dist.

      integer, parameter :: nb = 256
      real(dp), allocatable :: count(:)
      real(dp), allocatable :: wcount(:)
      real(dp) :: c0
      real(dp) :: c1
      real(dp) :: f0
      real(dp) :: f1
      real(dp) :: infl_w
      real(dp) :: phi
      real(dp) :: sumw
      real(dp) :: u
      real(dp) :: upper
      real(dp) :: unit_weights(size(dist))
      integer :: i
      integer :: j

      allocate(density(size(eval_dist)), influence_num(size(eval_dist)))
      density = 0.0_dp
      influence_num = 0.0_dp
      if (size(eval_dist) == 0 .or. bandwidth <= 0.0_dp) return
      upper = maxval(eval_dist) + 6.0_dp * bandwidth
      sumw = sum(weights)
      unit_weights = 1.0_dp
      allocate(count(nb + 1), wcount(nb + 1))
      call linbin(dist, 0.0_dp, upper, nb, unit_weights, count)
      call linbin(dist, 0.0_dp, upper, nb, weights, wcount)

      do j = 1, size(eval_dist)
         call boundary_kernel_coefficients(eval_dist(j) / bandwidth, c0, c1)
         f0 = 0.0_dp
         f1 = 0.0_dp
         do i = 1, size(dist)
            if (dist(i) > upper) cycle
            u = (eval_dist(j) - dist(i)) / bandwidth
            phi = r_dnorm(u)
            if (abs(u) <= 4.0_dp) f0 = f0 + weights(i) * phi / bandwidth
            if (abs(u) <= 5.0_dp) f1 = f1 - weights(i) * u * phi / bandwidth**2
         end do
         f0 = f0 / sumw
         f1 = f1 / sumw
         density(j) = c0 * f0 - bandwidth * c1 * f1
         infl_w = interpolated_bin_weight(eval_dist(j), upper, count, wcount, weights)
         influence_num(j) = k0 * c0 * infl_w / (real(size(dist), dp) * bandwidth)
         if (.not. ieee_is_finite(density(j)) .or. density(j) <= 0.0_dp &
            .or. .not. ieee_is_finite(influence_num(j))) then
            density(j) = 0.0_dp
            influence_num(j) = 0.0_dp
         end if
      end do
   end subroutine fit_boundary_component

   pure subroutine boundary_kernel_coefficients(a, c0, c1)
      real(dp), intent(in) :: a !! Standardized endpoint distance t/h.
      real(dp), intent(out) :: c0 !! Constant equivalent-kernel coefficient.
      real(dp), intent(out) :: c1 !! Linear equivalent-kernel coefficient.

      real(dp) :: denom
      real(dp) :: mu0
      real(dp) :: mu1
      real(dp) :: mu2
      real(dp) :: phi

      phi = k0 * exp(-0.5_dp * a * a)
      mu0 = 0.5_dp * erfc(-a / sqrt(2.0_dp))
      mu1 = -phi
      mu2 = mu0 - a * phi
      denom = mu0 * mu2 - mu1 * mu1
      c0 = mu2 / denom
      c1 = -mu1 / denom
   end subroutine boundary_kernel_coefficients

   pure real(dp) function interpolated_bin_weight(x, upper, count, wcount, weights) result(value)
      real(dp), intent(in) :: x !! Endpoint distance at which the binned average weight is interpolated.
      real(dp), intent(in) :: upper !! Upper edge of the regular 256-cell endpoint grid.
      real(dp), intent(in) :: count(:) !! Unweighted linearly binned row counts.
      real(dp), intent(in) :: wcount(:) !! Weighted linearly binned counts using mean-one case weights.
      real(dp), intent(in) :: weights(:) !! Original endpoint weights used to detect the equal-weight shortcut.

      real(dp) :: avg(size(count))
      real(dp) :: fraction
      real(dp) :: position
      integer :: bin
      integer :: nb

      if (maxval(weights) == minval(weights)) then
         value = 1.0_dp
         return
      end if
      avg = 0.0_dp
      where (count > 0.0_dp) avg = wcount / count
      nb = size(count) - 1
      position = max(0.0_dp, min(real(nb, dp), x * real(nb, dp) / upper))
      bin = min(int(floor(position)) + 1, nb)
      fraction = position - real(bin - 1, dp)
      value = (1.0_dp - fraction) * avg(bin) + fraction * avg(bin + 1)
   end function interpolated_bin_weight

   pure subroutine fuse_boundary_components(model, grid, w_lower, w_upper, lower_density, lower_infl, &
      upper_density, upper_infl, influences, ierr)
      type(kde1d_model), intent(inout) :: model !! Model whose interpolation grid is replaced by the fused endpoint estimate.
      real(dp), intent(in) :: grid(:) !! Increasing original-scale grid shared by bulk and fused estimates.
      real(dp), intent(in) :: w_lower(:) !! Lower-endpoint fusion weights on grid.
      real(dp), intent(in) :: w_upper(:) !! Upper-endpoint fusion weights on grid.
      real(dp), intent(in) :: lower_density(:) !! Lower boundary density component, possibly zero length.
      real(dp), intent(in) :: lower_infl(:) !! Lower boundary influence numerators, possibly zero length.
      real(dp), intent(in) :: upper_density(:) !! Upper boundary density in increasing grid order, possibly zero length.
      real(dp), intent(in) :: upper_infl(:) !! Upper boundary influence numerators in increasing grid order.
      real(dp), intent(inout) :: influences(:) !! Bulk influences replaced by fused influences.
      integer, intent(out) :: ierr !! Zero on success; nonzero if fused density normalization fails.

      real(dp) :: f(size(grid))
      real(dp) :: f_bulk(size(grid))
      real(dp) :: influence_num(size(grid))
      real(dp) :: wbulk
      integer :: i
      integer :: upper_start

      f_bulk = model%grid%values
      f = f_bulk
      influence_num = f_bulk * influences
      upper_start = size(grid) - size(upper_density) + 1
      do i = 1, size(grid)
         wbulk = 1.0_dp
         if (i <= size(lower_density) .and. w_lower(i) > 0.0_dp) wbulk = wbulk - w_lower(i)
         if (i >= upper_start .and. size(upper_density) > 0 .and. w_upper(i) > 0.0_dp) wbulk = wbulk - w_upper(i)
         f(i) = wbulk * f_bulk(i)
         influence_num(i) = wbulk * f_bulk(i) * influences(i)
         if (i <= size(lower_density) .and. w_lower(i) > 0.0_dp) then
            f(i) = f(i) + w_lower(i) * lower_density(i)
            influence_num(i) = influence_num(i) + w_lower(i) * lower_infl(i)
         end if
         if (i >= upper_start .and. size(upper_density) > 0 .and. w_upper(i) > 0.0_dp) then
            f(i) = f(i) + w_upper(i) * upper_density(i - upper_start + 1)
            influence_num(i) = influence_num(i) + w_upper(i) * upper_infl(i - upper_start + 1)
         end if
         if (f(i) > 0.0_dp) then
            influences(i) = influence_num(i) / f(i)
         else
            influences(i) = 0.0_dp
         end if
      end do
      call model%grid%initialize(grid, max(f, 0.0_dp), 3, ierr)
   end subroutine fuse_boundary_components

   pure function equi_jitter_codes(x) result(y)
      real(dp), intent(in) :: x(:) !! Integer-valued category codes jittered conditionally and equidistantly.
      real(dp) :: y(size(x))

      integer, allocatable :: order(:)
      real(dp), allocatable :: sorted(:)
      integer :: first
      integer :: i
      integer :: j
      integer :: last
      integer :: n_tie

      allocate(order(size(x)), sorted(size(x)))
      do i = 1, size(x)
         order(i) = i
         sorted(i) = x(i)
      end do
      call sort_values_with_order(sorted, order)
      y = x
      first = 1
      do while (first <= size(x))
         last = first
         do while (last < size(x))
            if (sorted(last + 1) /= sorted(first)) exit
            last = last + 1
         end do
         n_tie = last - first + 1
         do j = first, last
            y(order(j)) = sorted(j) - 0.5_dp + real(j - first + 1, dp) / real(n_tie + 1, dp)
         end do
         first = last + 1
      end do
   end function equi_jitter_codes

   pure subroutine sort_values_with_order(x, order)
      real(dp), intent(inout) :: x(:) !! Values sorted ascending in place.
      integer, intent(inout) :: order(:) !! Original one-based positions permuted stably with x.

      integer :: i
      integer :: j
      integer :: key_order
      real(dp) :: key

      do i = 2, size(x)
         key = x(i)
         key_order = order(i)
         j = i - 1
         do while (j >= 1)
            if (x(j) <= key) exit
            x(j + 1) = x(j)
            order(j + 1) = order(j)
            j = j - 1
         end do
         x(j + 1) = key
         order(j + 1) = key_order
      end do
   end subroutine sort_values_with_order

   pure subroutine sort_pairs_local(x, weights)
      real(dp), intent(inout) :: x(:) !! Values sorted ascending in place.
      real(dp), intent(inout) :: weights(:) !! Companion weights stably permuted with x.

      integer :: i
      integer :: j
      real(dp) :: key
      real(dp) :: wkey

      do i = 2, size(x)
         key = x(i)
         wkey = weights(i)
         j = i - 1
         do while (j >= 1)
            if (x(j) <= key) exit
            x(j + 1) = x(j)
            weights(j + 1) = weights(j)
            j = j - 1
         end do
         x(j + 1) = key
         weights(j + 1) = wkey
      end do
   end subroutine sort_pairs_local

   subroutine seed_fortran_rng(seed)
      integer, intent(in), optional :: seed !! Scalar seed expanded deterministically to the compiler RNG state.

      integer, allocatable :: put(:)
      integer :: i
      integer :: n
      integer :: base

      if (.not. present(seed)) return
      call random_seed(size=n)
      allocate(put(n))
      base = abs(seed)
      do i = 1, n
         put(i) = modulo(base + 104729 * i + 7919 * i * i, huge(1) - 1)
         if (put(i) == 0) put(i) = i
      end do
      call random_seed(put=put)
   end subroutine seed_fortran_rng

   pure real(dp) function sobol_1d(index, seed) result(value)
      integer, intent(in) :: index !! One-based quasi-random sequence index.
      integer, intent(in), optional :: seed !! Optional digital XOR shift; differs from randtoolbox scrambling streams.

      integer :: bits
      integer :: gray
      integer :: j
      integer :: shift

      gray = ieor(index - 1, shiftr(index - 1, 1))
      shift = 0
      if (present(seed)) shift = seed
      bits = ieor(gray, shift)
      value = 0.0_dp
      do j = 1, bit_size(bits) - 1
         if (btest(bits, j - 1)) value = value + 2.0_dp**(-j)
      end do
      value = max(0.0_dp, min(nearest(1.0_dp, -1.0_dp), value))
   end function sobol_1d

   pure real(dp) function quiet_nan() result(value)
      value = ieee_value(0.0_dp, ieee_quiet_nan)
   end function quiet_nan

end module kde1d_api
