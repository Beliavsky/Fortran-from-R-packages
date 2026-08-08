module nfcp_types
  use, intrinsic :: iso_fortran_env, only : real64, int64
  implicit none
  private

  integer, parameter, public :: dp = real64
  integer, parameter, public :: i8 = int64

  integer, parameter, public :: nfcp_ok = 0
  integer, parameter, public :: nfcp_invalid_input = 1
  integer, parameter, public :: nfcp_singular = 2
  integer, parameter, public :: nfcp_not_converged = 3
  integer, parameter, public :: nfcp_nonfinite = 4

  type, public :: nfcp_model_t
    integer :: n_factors = 0
    integer :: n_season = 0
    integer :: n_me = 0
    logical :: gbm = .true.
    real(dp) :: mu = 0.0_dp
    real(dp) :: mu_rn = 0.0_dp
    real(dp) :: equilibrium = 0.0_dp
    real(dp), allocatable :: lambda(:)
    real(dp), allocatable :: kappa(:)
    real(dp), allocatable :: sigma(:)
    real(dp), allocatable :: rho(:,:)
    real(dp), allocatable :: season_cos(:)
    real(dp), allocatable :: season_sin(:)
    real(dp), allocatable :: measurement_error(:)
  contains
    procedure :: valid => nfcp_model_valid
  end type nfcp_model_t

  type, public :: nfcp_parameterization_t
    integer :: n_factors = 0
    integer :: n_season = 0
    integer :: n_me = 0
    logical :: gbm = .true.
    logical :: estimate_initial_state = .false.
  contains
    procedure :: count => parameter_count
    procedure :: names => parameter_names
  end type nfcp_parameterization_t

  type, public :: nfcp_filter_result_t
    integer :: status = nfcp_invalid_input
    character(len=:), allocatable :: message
    real(dp) :: log_likelihood = -huge(1.0_dp)
    real(dp) :: aic = huge(1.0_dp)
    real(dp) :: bic = huge(1.0_dp)
    real(dp), allocatable :: final_state(:)
    real(dp), allocatable :: final_covariance(:,:)
    real(dp), allocatable :: states(:,:)
    real(dp), allocatable :: state_variances(:,:)
    real(dp), allocatable :: fitted_log_futures(:,:)
    real(dp), allocatable :: residuals(:,:)
    real(dp), allocatable :: log_likelihood_path(:)
  end type nfcp_filter_result_t

  type, public :: nfcp_simulation_result_t
    integer :: status = nfcp_invalid_input
    character(len=:), allocatable :: message
    real(dp), allocatable :: times(:)
    real(dp), allocatable :: states(:,:,:)
    real(dp), allocatable :: spot_prices(:,:)
  end type nfcp_simulation_result_t

  type, public :: nfcp_futures_simulation_result_t
    integer :: status = nfcp_invalid_input
    character(len=:), allocatable :: message
    real(dp), allocatable :: states(:,:)
    real(dp), allocatable :: futures_prices(:,:)
    real(dp), allocatable :: spot_prices(:)
  end type nfcp_futures_simulation_result_t

  type, public :: nfcp_option_result_t
    integer :: status = nfcp_invalid_input
    character(len=:), allocatable :: message
    real(dp) :: value = 0.0_dp
    real(dp) :: standard_error = 0.0_dp
    real(dp) :: annualized_volatility = 0.0_dp
    real(dp) :: delta = 0.0_dp
    real(dp) :: gamma = 0.0_dp
    real(dp) :: vega = 0.0_dp
    real(dp) :: theta = 0.0_dp
    real(dp) :: rho = 0.0_dp
    real(dp), allocatable :: exercise_probability(:)
  end type nfcp_option_result_t

  type, public :: nfcp_mle_control_t
    integer :: population_size = 0
    integer :: generations = 50
    integer :: wait_generations = 12
    integer :: local_max_evaluations = 3000
    integer :: seed = 12345
    real(dp) :: differential_weight = 0.8_dp
    real(dp) :: crossover_probability = 0.9_dp
    real(dp) :: solution_tolerance = 1.0e-6_dp
    logical :: local_refinement = .true.
    integer :: trace = 0
  end type nfcp_mle_control_t

  type, public :: nfcp_mle_result_t
    integer :: status = nfcp_invalid_input
    character(len=:), allocatable :: message
    logical :: converged = .false.
    real(dp) :: log_likelihood = -huge(1.0_dp)
    integer :: evaluations = 0
    integer :: generations = 0
    real(dp), allocatable :: parameters(:)
    real(dp), allocatable :: standard_errors(:)
    real(dp), allocatable :: hessian(:,:)
    type(nfcp_filter_result_t) :: filter
  end type nfcp_mle_result_t

  public :: initialize_model

contains

  subroutine initialize_model(model, n_factors, gbm, n_me, n_season)
    type(nfcp_model_t), intent(out) :: model
    integer, intent(in) :: n_factors
    logical, intent(in), optional :: gbm
    integer, intent(in), optional :: n_me, n_season
    integer :: i

    model%n_factors = max(0, n_factors)
    model%gbm = .true.
    if (present(gbm)) model%gbm = gbm
    model%n_me = 1
    if (present(n_me)) model%n_me = max(0, n_me)
    model%n_season = 0
    if (present(n_season)) model%n_season = max(0, n_season)

    allocate(model%lambda(model%n_factors), model%kappa(model%n_factors), &
             model%sigma(model%n_factors), model%rho(model%n_factors, model%n_factors))
    allocate(model%season_cos(model%n_season), model%season_sin(model%n_season))
    allocate(model%measurement_error(model%n_me))

    model%lambda = 0.0_dp
    model%kappa = 1.0_dp
    model%sigma = 0.2_dp
    model%rho = 0.0_dp
    do i = 1, model%n_factors
      model%rho(i,i) = 1.0_dp
    end do
    model%season_cos = 0.0_dp
    model%season_sin = 0.0_dp
    model%measurement_error = 0.01_dp
    if (model%gbm .and. model%n_factors > 0) model%kappa(1) = 0.0_dp
  end subroutine initialize_model

  logical function nfcp_model_valid(self, message) result(ok)
    class(nfcp_model_t), intent(in) :: self
    character(len=:), allocatable, intent(out), optional :: message
    integer :: i

    ok = .false.
    if (self%n_factors < 1) then
      if (present(message)) message = 'n_factors must be at least one'
      return
    end if
    if (.not. allocated(self%lambda) .or. size(self%lambda) /= self%n_factors .or. &
        .not. allocated(self%kappa) .or. size(self%kappa) /= self%n_factors .or. &
        .not. allocated(self%sigma) .or. size(self%sigma) /= self%n_factors) then
      if (present(message)) message = 'factor arrays have inconsistent dimensions'
      return
    end if
    if (.not. allocated(self%rho) .or. any(shape(self%rho) /= [self%n_factors, self%n_factors])) then
      if (present(message)) message = 'correlation matrix has inconsistent dimensions'
      return
    end if
    if (any(self%sigma < 0.0_dp) .or. any(self%kappa < 0.0_dp)) then
      if (present(message)) message = 'sigma and kappa must be nonnegative'
      return
    end if
    if (maxval(abs(self%rho - transpose(self%rho))) > 1.0e-10_dp .or. &
        any(abs(self%rho) > 1.0_dp + 1.0e-12_dp)) then
      if (present(message)) message = 'rho must be a symmetric correlation matrix'
      return
    end if
    do i = 1, self%n_factors
      if (abs(self%rho(i,i) - 1.0_dp) > 1.0e-10_dp) then
        if (present(message)) message = 'rho diagonal elements must equal one'
        return
      end if
    end do
    if (self%gbm .and. abs(self%kappa(1)) > 1.0e-12_dp) then
      if (present(message)) message = 'kappa(1) must be zero for a GBM first factor'
      return
    end if
    if (self%n_me > 0) then
      if (.not. allocated(self%measurement_error) .or. size(self%measurement_error) /= self%n_me) then
        if (present(message)) message = 'measurement-error array has inconsistent dimensions'
        return
      end if
      if (any(self%measurement_error < 0.0_dp)) then
        if (present(message)) message = 'measurement errors must be nonnegative'
        return
      end if
    end if
    if (self%n_season > 0) then
      if (.not. allocated(self%season_cos) .or. .not. allocated(self%season_sin) .or. &
          size(self%season_cos) /= self%n_season .or. size(self%season_sin) /= self%n_season) then
        if (present(message)) message = 'seasonal arrays have inconsistent dimensions'
        return
      end if
    end if
    ok = .true.
    if (present(message)) message = 'ok'
  end function nfcp_model_valid

  integer function parameter_count(self) result(n)
    class(nfcp_parameterization_t), intent(in) :: self
    n = merge(self%n_factors, 0, self%estimate_initial_state)
    n = n + 1
    n = n + self%n_factors
    n = n + self%n_factors - merge(1, 0, self%gbm)
    n = n + self%n_factors
    n = n + self%n_factors * (self%n_factors - 1) / 2
    n = n + 2 * self%n_season + self%n_me
  end function parameter_count

  function parameter_names(self) result(names)
    class(nfcp_parameterization_t), intent(in) :: self
    character(len=32), allocatable :: names(:)
    integer :: i, j, k

    allocate(names(self%count()))
    names = ''
    k = 0
    if (self%estimate_initial_state) then
      do i = 1, self%n_factors
        k = k + 1
        write(names(k),'("x_0_",i0)') i
      end do
    end if
    k = k + 1
    names(k) = merge('mu', 'E ', self%gbm)
    do i = 1, self%n_factors
      k = k + 1
      if (self%gbm .and. i == 1) then
        names(k) = 'mu_rn'
      else
        write(names(k),'("lambda_",i0)') i
      end if
    end do
    do i = merge(2, 1, self%gbm), self%n_factors
      k = k + 1
      write(names(k),'("kappa_",i0)') i
    end do
    do i = 1, self%n_factors
      k = k + 1
      write(names(k),'("sigma_",i0)') i
    end do
    do i = 1, self%n_factors - 1
      do j = i + 1, self%n_factors
        k = k + 1
        write(names(k),'("rho_",i0,"_",i0)') i, j
      end do
    end do
    do i = 1, self%n_season
      k = k + 1
      write(names(k),'("season_",i0,"_1")') i
      k = k + 1
      write(names(k),'("season_",i0,"_2")') i
    end do
    do i = 1, self%n_me
      k = k + 1
      write(names(k),'("ME_",i0)') i
    end do
  end function parameter_names

end module nfcp_types
