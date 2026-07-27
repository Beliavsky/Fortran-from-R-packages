! SPDX-License-Identifier: GPL-3.0-only
module svdnf_models
  use svdnf_kinds, only : dp
  use svdnf_stats, only : poisson_pmf, binomial_pmf, random_poisson, random_bernoulli
  use svdnf_types
  implicit none
  private
  public :: dynamics_svm, set_custom_dynamics, validate_dynamics
  public :: evaluate_mu_y, evaluate_sigma_y, evaluate_mu_x, evaluate_sigma_x
  public :: jump_probability, draw_jump_count, parameter_vector, set_parameter_vector
  public :: parameter_count, model_parameter_names

contains

  function dynamics_svm(model, mu, kappa, theta, sigma, rho, omega, delta, alpha, &
      rho_z, nu, p, phi, h, coefs) result(dynamics)
    character(len=*), intent(in), optional :: model
    real(dp), intent(in), optional :: mu, kappa, theta, sigma, rho, omega, delta
    real(dp), intent(in), optional :: alpha, rho_z, nu, p, phi, h
    real(dp), intent(in), optional :: coefs(:)
    type(svm_dynamics) :: dynamics
    character(len=32) :: selected

    if (present(mu)) dynamics%mu = mu
    if (present(kappa)) dynamics%kappa = kappa
    if (present(theta)) dynamics%theta = theta
    if (present(sigma)) dynamics%sigma = sigma
    if (present(rho)) dynamics%rho = rho
    if (present(omega)) dynamics%omega = omega
    if (present(delta)) dynamics%delta = delta
    if (present(alpha)) dynamics%alpha = alpha
    if (present(rho_z)) dynamics%rho_z = rho_z
    if (present(nu)) dynamics%nu = nu
    if (present(p)) dynamics%p = p
    if (present(phi)) dynamics%phi = phi
    if (present(h)) dynamics%h = h
    if (present(coefs)) then
      allocate(dynamics%coefs(size(coefs)))
      dynamics%coefs = coefs
    end if

    selected = 'Heston'
    if (present(model)) selected = adjustl(model)
    select case (trim(lowercase(selected)))
    case ('duffiepansingleton','duffie_pan_singleton','dps')
      dynamics%model_id = model_duffie_pan_singleton
      dynamics%model = 'DuffiePanSingleton'
    case ('bates')
      dynamics%model_id = model_bates
      dynamics%model = 'Bates'
      dynamics%rho_z = 0.0_dp
      dynamics%nu = 0.0_dp
    case ('heston')
      dynamics%model_id = model_heston
      dynamics%model = 'Heston'
      dynamics%rho_z = 0.0_dp
      dynamics%nu = 0.0_dp
      dynamics%alpha = 0.0_dp
      dynamics%delta = 0.0_dp
    case ('pittmalikdoucet','pitt_malik_doucet','pmd')
      dynamics%model_id = model_pitt_malik_doucet
      dynamics%model = 'PittMalikDoucet'
      dynamics%rho_z = 0.0_dp
      dynamics%nu = 0.0_dp
    case ('taylorwithleverage','taylor_with_leverage')
      dynamics%model_id = model_taylor_leverage
      dynamics%model = 'TaylorWithLeverage'
      dynamics%rho_z = 0.0_dp
      dynamics%nu = 0.0_dp
      dynamics%alpha = 0.0_dp
      dynamics%delta = 0.0_dp
    case ('taylor')
      dynamics%model_id = model_taylor
      dynamics%model = 'Taylor'
      dynamics%rho = 0.0_dp
      dynamics%rho_z = 0.0_dp
      dynamics%nu = 0.0_dp
      dynamics%alpha = 0.0_dp
      dynamics%delta = 0.0_dp
    case ('capm_sv','capmsv')
      dynamics%model_id = model_capm_sv
      dynamics%model = 'CAPM_SV'
      dynamics%rho = 0.0_dp
      dynamics%rho_z = 0.0_dp
      dynamics%nu = 0.0_dp
      dynamics%alpha = 0.0_dp
      dynamics%delta = 0.0_dp
      if (.not. allocated(dynamics%coefs)) dynamics%coefs = [0.0_dp,1.0_dp]
    case ('custom')
      dynamics%model_id = model_custom
      dynamics%model = 'Custom'
    case default
      dynamics%model_id = -1
      dynamics%model = selected
    end select
  end function dynamics_svm

  subroutine set_custom_dynamics(dynamics, mu_y, sigma_y, mu_x, sigma_x, &
      mu_y_parameters, sigma_y_parameters, mu_x_parameters, sigma_x_parameters, &
      jump_probability_function_value, jump_count_function_value, jump_parameters)
    type(svm_dynamics), intent(inout) :: dynamics
    procedure(state_function) :: mu_y, sigma_y, mu_x, sigma_x
    real(dp), intent(in), optional :: mu_y_parameters(:), sigma_y_parameters(:)
    real(dp), intent(in), optional :: mu_x_parameters(:), sigma_x_parameters(:)
    procedure(jump_probability_function), optional :: jump_probability_function_value
    procedure(jump_count_function), optional :: jump_count_function_value
    real(dp), intent(in), optional :: jump_parameters(:)
    dynamics%model_id = model_custom
    dynamics%model = 'Custom'
    dynamics%custom_mu_y => mu_y
    dynamics%custom_sigma_y => sigma_y
    dynamics%custom_mu_x => mu_x
    dynamics%custom_sigma_x => sigma_x
    if (present(mu_y_parameters)) dynamics%mu_y_parameters = mu_y_parameters
    if (present(sigma_y_parameters)) dynamics%sigma_y_parameters = sigma_y_parameters
    if (present(mu_x_parameters)) dynamics%mu_x_parameters = mu_x_parameters
    if (present(sigma_x_parameters)) dynamics%sigma_x_parameters = sigma_x_parameters
    if (present(jump_probability_function_value)) then
      dynamics%custom_jump_probability => jump_probability_function_value
    end if
    if (present(jump_count_function_value)) dynamics%custom_jump_count => jump_count_function_value
    if (present(jump_parameters)) dynamics%jump_parameters = jump_parameters
  end subroutine set_custom_dynamics

  real(dp) function evaluate_mu_y(dynamics, x) result(value)
    type(svm_dynamics), intent(in) :: dynamics
    real(dp), intent(in) :: x
    real(dp) :: alpha_bar
    select case (dynamics%model_id)
    case (model_duffie_pan_singleton)
      alpha_bar = exp(dynamics%alpha + 0.5_dp*dynamics%delta**2) / &
        max(1.0_dp-dynamics%rho_z*dynamics%nu,tiny(1.0_dp)) - 1.0_dp
      value = dynamics%h*(dynamics%mu - 0.5_dp*x - alpha_bar*dynamics%omega)
    case (model_bates)
      alpha_bar = exp(dynamics%alpha + 0.5_dp*dynamics%delta**2) - 1.0_dp
      value = dynamics%h*(dynamics%mu - 0.5_dp*x - alpha_bar*dynamics%omega)
    case (model_heston)
      value = dynamics%h*(dynamics%mu - 0.5_dp*x)
    case (model_taylor, model_taylor_leverage, model_pitt_malik_doucet, model_capm_sv)
      value = 0.0_dp
    case (model_custom)
      if (associated(dynamics%custom_mu_y)) then
        if (allocated(dynamics%mu_y_parameters)) then
          value = dynamics%custom_mu_y(x,dynamics%mu_y_parameters)
        else
          value = dynamics%custom_mu_y(x,[real(dp) ::])
        end if
      else
        value = 0.0_dp
      end if
    case default
      value = 0.0_dp
    end select
  end function evaluate_mu_y

  real(dp) function evaluate_sigma_y(dynamics, x) result(value)
    type(svm_dynamics), intent(in) :: dynamics
    real(dp), intent(in) :: x
    select case (dynamics%model_id)
    case (model_duffie_pan_singleton, model_bates, model_heston)
      value = sqrt(dynamics%h*max(x,0.0_dp))
    case (model_taylor, model_taylor_leverage, model_pitt_malik_doucet, model_capm_sv)
      value = exp(0.5_dp*x)
    case (model_custom)
      if (associated(dynamics%custom_sigma_y)) then
        if (allocated(dynamics%sigma_y_parameters)) then
          value = dynamics%custom_sigma_y(x,dynamics%sigma_y_parameters)
        else
          value = dynamics%custom_sigma_y(x,[real(dp) ::])
        end if
      else
        value = 0.0_dp
      end if
    case default
      value = 0.0_dp
    end select
  end function evaluate_sigma_y

  real(dp) function evaluate_mu_x(dynamics, x) result(value)
    type(svm_dynamics), intent(in) :: dynamics
    real(dp), intent(in) :: x
    select case (dynamics%model_id)
    case (model_duffie_pan_singleton, model_bates, model_heston)
      value = x + dynamics%h*dynamics%kappa*(dynamics%theta-max(x,0.0_dp))
    case (model_taylor, model_taylor_leverage, model_pitt_malik_doucet, model_capm_sv)
      value = dynamics%theta + dynamics%phi*(x-dynamics%theta)
    case (model_custom)
      if (associated(dynamics%custom_mu_x)) then
        if (allocated(dynamics%mu_x_parameters)) then
          value = dynamics%custom_mu_x(x,dynamics%mu_x_parameters)
        else
          value = dynamics%custom_mu_x(x,[real(dp) ::])
        end if
      else
        value = x
      end if
    case default
      value = x
    end select
  end function evaluate_mu_x

  real(dp) function evaluate_sigma_x(dynamics, x) result(value)
    type(svm_dynamics), intent(in) :: dynamics
    real(dp), intent(in) :: x
    select case (dynamics%model_id)
    case (model_duffie_pan_singleton, model_bates, model_heston)
      value = dynamics%sigma*sqrt(dynamics%h*max(x,0.0_dp))
    case (model_taylor, model_taylor_leverage, model_pitt_malik_doucet, model_capm_sv)
      value = dynamics%sigma
    case (model_custom)
      if (associated(dynamics%custom_sigma_x)) then
        if (allocated(dynamics%sigma_x_parameters)) then
          value = dynamics%custom_sigma_x(x,dynamics%sigma_x_parameters)
        else
          value = dynamics%custom_sigma_x(x,[real(dp) ::])
        end if
      else
        value = 0.0_dp
      end if
    case default
      value = 0.0_dp
    end select
  end function evaluate_sigma_x

  real(dp) function jump_probability(dynamics, n) result(value)
    type(svm_dynamics), intent(in) :: dynamics
    integer, intent(in) :: n
    select case (dynamics%model_id)
    case (model_duffie_pan_singleton, model_bates)
      value = poisson_pmf(n,dynamics%h*dynamics%omega)
    case (model_pitt_malik_doucet)
      value = binomial_pmf(n,1,dynamics%p)
    case (model_custom)
      if (associated(dynamics%custom_jump_probability)) then
        if (allocated(dynamics%jump_parameters)) then
          value = dynamics%custom_jump_probability(n,dynamics%jump_parameters)
        else
          value = dynamics%custom_jump_probability(n,[real(dp) ::])
        end if
      else
        value = merge(1.0_dp,0.0_dp,n==0)
      end if
    case default
      value = merge(1.0_dp,0.0_dp,n==0)
    end select
  end function jump_probability

  integer function draw_jump_count(dynamics) result(value)
    type(svm_dynamics), intent(in) :: dynamics
    select case (dynamics%model_id)
    case (model_duffie_pan_singleton, model_bates)
      value = random_poisson(dynamics%h*dynamics%omega)
    case (model_pitt_malik_doucet)
      value = random_bernoulli(dynamics%p)
    case (model_custom)
      if (associated(dynamics%custom_jump_count)) then
        if (allocated(dynamics%jump_parameters)) then
          value = dynamics%custom_jump_count(dynamics%jump_parameters)
        else
          value = dynamics%custom_jump_count([real(dp) ::])
        end if
      else
        value = 0
      end if
    case default
      value = 0
    end select
  end function draw_jump_count

  subroutine validate_dynamics(dynamics, ok, message)
    type(svm_dynamics), intent(in) :: dynamics
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message
    ok = .false.
    message = ''
    if (dynamics%model_id < model_custom .or. dynamics%model_id > model_capm_sv) then
      message = 'Unknown stochastic-volatility model.'
      return
    end if
    if (dynamics%sigma <= 0.0_dp) then
      message = 'sigma must be positive.'
      return
    end if
    if (dynamics%rho <= -1.0_dp .or. dynamics%rho >= 1.0_dp) then
      message = 'rho must be strictly between -1 and 1.'
      return
    end if
    if (dynamics%delta < 0.0_dp .or. dynamics%nu < 0.0_dp) then
      message = 'delta and nu must be nonnegative.'
      return
    end if
    select case (dynamics%model_id)
    case (model_duffie_pan_singleton,model_bates,model_heston)
      if (dynamics%kappa <= 0.0_dp .or. dynamics%theta <= 0.0_dp .or. dynamics%h <= 0.0_dp) then
        message = 'kappa, theta, and h must be positive.'
        return
      end if
      if (dynamics%model_id /= model_heston .and. dynamics%omega < 0.0_dp) then
        message = 'omega must be nonnegative.'
        return
      end if
    case (model_taylor,model_taylor_leverage,model_pitt_malik_doucet,model_capm_sv)
      if (abs(dynamics%phi) >= 1.0_dp) then
        message = 'The absolute value of phi must be less than 1.'
        return
      end if
      if (dynamics%model_id == model_pitt_malik_doucet) then
        if (dynamics%p < 0.0_dp .or. dynamics%p > 1.0_dp) then
          message = 'p must be between 0 and 1.'
          return
        end if
      end if
    case (model_custom)
      if (.not. associated(dynamics%custom_mu_y) .or. .not. associated(dynamics%custom_sigma_y) .or. &
          .not. associated(dynamics%custom_mu_x) .or. .not. associated(dynamics%custom_sigma_x)) then
        message = 'Custom models require four state-function callbacks.'
        return
      end if
    end select
    ok = .true.
  end subroutine validate_dynamics

  integer function parameter_count(dynamics) result(n)
    type(svm_dynamics), intent(in) :: dynamics
    select case (dynamics%model_id)
    case (model_duffie_pan_singleton); n = 10
    case (model_bates); n = 8
    case (model_heston); n = 5
    case (model_pitt_malik_doucet); n = 7
    case (model_taylor_leverage); n = 4
    case (model_taylor); n = 3
    case (model_capm_sv)
      if (allocated(dynamics%coefs)) then
        n = size(dynamics%coefs) + 3
      else
        n = 5
      end if
    case default; n = 0
    end select
  end function parameter_count

  function parameter_vector(dynamics) result(parameters)
    type(svm_dynamics), intent(in) :: dynamics
    real(dp), allocatable :: parameters(:)
    select case (dynamics%model_id)
    case (model_duffie_pan_singleton)
      parameters = [dynamics%mu,dynamics%alpha,dynamics%delta,dynamics%rho_z,dynamics%nu, &
        dynamics%omega,dynamics%kappa,dynamics%theta,dynamics%sigma,dynamics%rho]
    case (model_bates)
      parameters = [dynamics%mu,dynamics%alpha,dynamics%delta,dynamics%omega,dynamics%kappa, &
        dynamics%theta,dynamics%sigma,dynamics%rho]
    case (model_heston)
      parameters = [dynamics%mu,dynamics%kappa,dynamics%theta,dynamics%sigma,dynamics%rho]
    case (model_pitt_malik_doucet)
      parameters = [dynamics%phi,dynamics%theta,dynamics%sigma,dynamics%rho,dynamics%delta, &
        dynamics%alpha,dynamics%p]
    case (model_taylor_leverage)
      parameters = [dynamics%phi,dynamics%theta,dynamics%sigma,dynamics%rho]
    case (model_taylor)
      parameters = [dynamics%phi,dynamics%theta,dynamics%sigma]
    case (model_capm_sv)
      if (allocated(dynamics%coefs)) then
        parameters = [dynamics%coefs,dynamics%phi,dynamics%theta,dynamics%sigma]
      else
        parameters = [0.0_dp,1.0_dp,dynamics%phi,dynamics%theta,dynamics%sigma]
      end if
    case default
      allocate(parameters(0))
    end select
  end function parameter_vector

  subroutine set_parameter_vector(dynamics, parameters, ok)
    type(svm_dynamics), intent(inout) :: dynamics
    real(dp), intent(in) :: parameters(:)
    logical, intent(out), optional :: ok
    integer :: k
    logical :: valid
    valid = .true.
    if (size(parameters) /= parameter_count(dynamics)) then
      valid = .false.
    else
      select case (dynamics%model_id)
      case (model_duffie_pan_singleton)
        dynamics%mu=parameters(1); dynamics%alpha=parameters(2); dynamics%delta=parameters(3)
        dynamics%rho_z=parameters(4); dynamics%nu=parameters(5); dynamics%omega=parameters(6)
        dynamics%kappa=parameters(7); dynamics%theta=parameters(8); dynamics%sigma=parameters(9)
        dynamics%rho=parameters(10)
      case (model_bates)
        dynamics%mu=parameters(1); dynamics%alpha=parameters(2); dynamics%delta=parameters(3)
        dynamics%omega=parameters(4); dynamics%kappa=parameters(5); dynamics%theta=parameters(6)
        dynamics%sigma=parameters(7); dynamics%rho=parameters(8)
      case (model_heston)
        dynamics%mu=parameters(1); dynamics%kappa=parameters(2); dynamics%theta=parameters(3)
        dynamics%sigma=parameters(4); dynamics%rho=parameters(5)
      case (model_pitt_malik_doucet)
        dynamics%phi=parameters(1); dynamics%theta=parameters(2); dynamics%sigma=parameters(3)
        dynamics%rho=parameters(4); dynamics%delta=parameters(5); dynamics%alpha=parameters(6)
        dynamics%p=parameters(7)
      case (model_taylor_leverage)
        dynamics%phi=parameters(1); dynamics%theta=parameters(2); dynamics%sigma=parameters(3)
        dynamics%rho=parameters(4)
      case (model_taylor)
        dynamics%phi=parameters(1); dynamics%theta=parameters(2); dynamics%sigma=parameters(3)
      case (model_capm_sv)
        k = size(parameters)-3
        dynamics%coefs = parameters(1:k)
        dynamics%phi=parameters(k+1); dynamics%theta=parameters(k+2); dynamics%sigma=parameters(k+3)
      case default
        valid = .false.
      end select
    end if
    if (present(ok)) ok = valid
  end subroutine set_parameter_vector

  function model_parameter_names(dynamics) result(names)
    type(svm_dynamics), intent(in) :: dynamics
    character(len=24), allocatable :: names(:)
    integer :: i, k
    select case (dynamics%model_id)
    case (model_duffie_pan_singleton)
      names = [character(len=24) :: 'mu','alpha','delta','rho_z','nu','omega','kappa','theta','sigma','rho']
    case (model_bates)
      names = [character(len=24) :: 'mu','alpha','delta','omega','kappa','theta','sigma','rho']
    case (model_heston)
      names = [character(len=24) :: 'mu','kappa','theta','sigma','rho']
    case (model_pitt_malik_doucet)
      names = [character(len=24) :: 'phi','theta','sigma','rho','delta','alpha','p']
    case (model_taylor_leverage)
      names = [character(len=24) :: 'phi','theta','sigma','rho']
    case (model_taylor)
      names = [character(len=24) :: 'phi','theta','sigma']
    case (model_capm_sv)
      k = max(0,parameter_count(dynamics)-3)
      allocate(names(k+3))
      do i = 1, k
        write(names(i),'("c_",i0)') i-1
      end do
      names(k+1:k+3) = [character(len=24) :: 'phi','theta','sigma']
    case default
      allocate(names(0))
    end select
  end function model_parameter_names

  pure function lowercase(text) result(value)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: value
    integer :: i, code
    value = text
    do i = 1, len(text)
      code = iachar(value(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) value(i:i)=achar(code+32)
    end do
  end function lowercase

end module svdnf_models
