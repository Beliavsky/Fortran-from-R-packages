module nfcp_parameters
  use nfcp_types, only : dp, nfcp_model_t, nfcp_parameterization_t, initialize_model, &
                         nfcp_ok, nfcp_invalid_input
  implicit none
  private
  public :: pack_parameters, unpack_parameters, default_parameter_bounds

contains

  subroutine pack_parameters(model, spec, initial_state, theta, status)
    type(nfcp_model_t), intent(in) :: model
    type(nfcp_parameterization_t), intent(in) :: spec
    real(dp), intent(in), optional :: initial_state(:)
    real(dp), allocatable, intent(out) :: theta(:)
    integer, intent(out) :: status
    integer :: i, j, k

    status = nfcp_invalid_input
    if (model%n_factors /= spec%n_factors .or. model%n_season /= spec%n_season .or. &
        model%n_me /= spec%n_me .or. model%gbm .neqv. spec%gbm) return
    if (spec%estimate_initial_state) then
      if (.not. present(initial_state)) return
      if (size(initial_state) /= spec%n_factors) return
    end if

    allocate(theta(spec%count()))
    k = 0
    if (spec%estimate_initial_state) then
      theta(1:spec%n_factors) = initial_state
      k = spec%n_factors
    end if
    k = k + 1
    theta(k) = merge(model%mu, model%equilibrium, spec%gbm)
    do i = 1, spec%n_factors
      k = k + 1
      if (spec%gbm .and. i == 1) then
        theta(k) = model%mu_rn
      else
        theta(k) = model%lambda(i)
      end if
    end do
    do i = merge(2,1,spec%gbm), spec%n_factors
      k = k + 1
      theta(k) = model%kappa(i)
    end do
    do i = 1, spec%n_factors
      k = k + 1
      theta(k) = model%sigma(i)
    end do
    do i = 1, spec%n_factors-1
      do j = i+1, spec%n_factors
        k = k + 1
        theta(k) = model%rho(i,j)
      end do
    end do
    do i = 1, spec%n_season
      k = k + 1; theta(k) = model%season_cos(i)
      k = k + 1; theta(k) = model%season_sin(i)
    end do
    do i = 1, spec%n_me
      k = k + 1; theta(k) = model%measurement_error(i)
    end do
    status = nfcp_ok
  end subroutine pack_parameters

  subroutine unpack_parameters(theta, spec, model, initial_state, status)
    real(dp), intent(in) :: theta(:)
    type(nfcp_parameterization_t), intent(in) :: spec
    type(nfcp_model_t), intent(out) :: model
    real(dp), allocatable, intent(out) :: initial_state(:)
    integer, intent(out) :: status
    integer :: i, j, k

    status = nfcp_invalid_input
    if (size(theta) /= spec%count()) return
    call initialize_model(model, spec%n_factors, spec%gbm, spec%n_me, spec%n_season)
    allocate(initial_state(spec%n_factors))
    initial_state = 0.0_dp
    k = 0
    if (spec%estimate_initial_state) then
      initial_state = theta(1:spec%n_factors)
      k = spec%n_factors
    end if
    k = k + 1
    if (spec%gbm) then
      model%mu = theta(k)
      model%equilibrium = 0.0_dp
    else
      model%equilibrium = theta(k)
      model%mu = 0.0_dp
    end if
    do i = 1, spec%n_factors
      k = k + 1
      if (spec%gbm .and. i == 1) then
        model%mu_rn = theta(k)
        model%lambda(1) = model%mu-model%mu_rn
      else
        model%lambda(i) = theta(k)
      end if
    end do
    if (spec%gbm) model%kappa(1) = 0.0_dp
    do i = merge(2,1,spec%gbm), spec%n_factors
      k = k + 1
      model%kappa(i) = theta(k)
    end do
    do i = 1, spec%n_factors
      k = k + 1
      model%sigma(i) = theta(k)
    end do
    model%rho = 0.0_dp
    do i = 1, spec%n_factors
      model%rho(i,i) = 1.0_dp
    end do
    do i = 1, spec%n_factors-1
      do j = i+1, spec%n_factors
        k = k + 1
        model%rho(i,j) = theta(k)
        model%rho(j,i) = theta(k)
      end do
    end do
    do i = 1, spec%n_season
      k = k + 1; model%season_cos(i) = theta(k)
      k = k + 1; model%season_sin(i) = theta(k)
    end do
    do i = 1, spec%n_me
      k = k + 1; model%measurement_error(i) = theta(k)
    end do
    status = nfcp_ok
  end subroutine unpack_parameters

  subroutine default_parameter_bounds(spec, lower, upper)
    type(nfcp_parameterization_t), intent(in) :: spec
    real(dp), allocatable, intent(out) :: lower(:), upper(:)
    character(len=32), allocatable :: names(:)
    integer :: i

    names = spec%names()
    allocate(lower(size(names)), upper(size(names)))
    do i = 1, size(names)
      if (index(names(i),'x_0_') == 1) then
        lower(i) = -10.0_dp; upper(i) = 10.0_dp
      else if (trim(names(i)) == 'E') then
        lower(i) = -10.0_dp; upper(i) = 10.0_dp
      else if (trim(names(i)) == 'mu' .or. trim(names(i)) == 'mu_rn') then
        lower(i) = -10.0_dp; upper(i) = 10.0_dp
      else if (index(names(i),'lambda_') == 1) then
        lower(i) = -10.0_dp; upper(i) = 10.0_dp
      else if (index(names(i),'kappa_') == 1) then
        lower(i) = 1.0e-5_dp; upper(i) = 50.0_dp
      else if (index(names(i),'sigma_') == 1) then
        lower(i) = 1.0e-8_dp; upper(i) = 10.0_dp
      else if (index(names(i),'rho_') == 1) then
        lower(i) = -0.999_dp; upper(i) = 0.999_dp
      else if (index(names(i),'season_') == 1) then
        lower(i) = -1.0_dp; upper(i) = 1.0_dp
      else if (index(names(i),'ME_') == 1) then
        lower(i) = 1.0e-5_dp; upper(i) = 1.0_dp
      else
        lower(i) = -10.0_dp; upper(i) = 10.0_dp
      end if
    end do
  end subroutine default_parameter_bounds

end module nfcp_parameters
