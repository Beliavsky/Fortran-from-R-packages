! SPDX-License-Identifier: GPL-3.0-only
module svdnf_simulation
  use svdnf_kinds, only : dp
  use svdnf_types, only : svm_dynamics, simulation_result, filter_result, forecast_result
  use svdnf_models, only : validate_dynamics, evaluate_mu_y, evaluate_sigma_y, &
    evaluate_mu_x, evaluate_sigma_x, draw_jump_count
  use svdnf_stats, only : seed_random, random_normal, random_gamma, sample_discrete, &
    normal_quantile, mean_value, standard_deviation
  implicit none
  private
  public :: model_simulate, model_sim, predict_filter, predict_svdnf

contains

  function model_simulate(dynamics, n_steps, initial_volatility, seed) result(output)
    type(svm_dynamics), intent(in) :: dynamics
    integer, intent(in) :: n_steps
    real(dp), intent(in), optional :: initial_volatility
    integer, intent(in), optional :: seed
    type(simulation_result) :: output
    real(dp) :: x_previous, mu_y, sigma_y, mu_x, sigma_x, epsilon_x, epsilon_y
    integer :: i, count
    logical :: ok
    character(len=160) :: message

    call validate_dynamics(dynamics,ok,message)
    if (.not. ok) then
      output%ok = .false.
      output%message = message
      return
    end if
    if (n_steps < 1) then
      output%ok = .false.
      output%message = 'n_steps must be positive.'
      return
    end if
    if (present(seed)) call seed_random(seed)
    allocate(output%volatility_factor(n_steps),output%returns(n_steps))
    allocate(output%jump_counts(n_steps),output%volatility_jumps(n_steps),output%return_jumps(n_steps))
    x_previous = dynamics%theta
    if (present(initial_volatility)) x_previous = initial_volatility
    do i = 1, n_steps
      count = draw_jump_count(dynamics)
      output%jump_counts(i) = count
      if (count > 0 .and. dynamics%nu > 0.0_dp) then
        output%volatility_jumps(i) = random_gamma(real(count,dp),dynamics%nu)
      else
        output%volatility_jumps(i) = 0.0_dp
      end if
      if (count > 0 .and. dynamics%delta > 0.0_dp) then
        output%return_jumps(i) = real(count,dp)*dynamics%alpha + &
          sqrt(real(count,dp))*dynamics%delta*random_normal() + &
          dynamics%rho_z*output%volatility_jumps(i)
      else
        output%return_jumps(i) = real(count,dp)*dynamics%alpha + &
          dynamics%rho_z*output%volatility_jumps(i)
      end if
      epsilon_x = random_normal()
      epsilon_y = random_normal()
      mu_y = evaluate_mu_y(dynamics,x_previous)
      sigma_y = evaluate_sigma_y(dynamics,x_previous)
      mu_x = evaluate_mu_x(dynamics,x_previous)
      sigma_x = evaluate_sigma_x(dynamics,x_previous)
      output%returns(i) = mu_y + sigma_y*(dynamics%rho*epsilon_x + &
        sqrt(max(1.0_dp-dynamics%rho**2,0.0_dp))*epsilon_y) + output%return_jumps(i)
      output%volatility_factor(i) = mu_x + sigma_x*epsilon_x + output%volatility_jumps(i)
      x_previous = output%volatility_factor(i)
    end do
    output%ok = .true.
  end function model_simulate

  function model_sim(dynamics, t, init_vol, seed) result(output)
    type(svm_dynamics), intent(in) :: dynamics
    integer, intent(in) :: t
    real(dp), intent(in), optional :: init_vol
    integer, intent(in), optional :: seed
    type(simulation_result) :: output
    output = model_simulate(dynamics,t,init_vol,seed)
  end function model_sim

  function predict_filter(filtered, n_ahead, n_sim, confidence, new_factors, seed, &
      upstream_filter_index) result(output)
    type(filter_result), intent(in) :: filtered
    integer, intent(in), optional :: n_ahead, n_sim
    real(dp), intent(in), optional :: confidence
    real(dp), intent(in), optional :: new_factors(:,:)
    integer, intent(in), optional :: seed
    logical, intent(in), optional :: upstream_filter_index
    type(forecast_result) :: output
    integer :: horizon, simulations, i, j, index_value, filter_column
    real(dp) :: conf, z
    real(dp), allocatable :: vol_paths(:,:), return_paths(:,:), probabilities(:)
    type(simulation_result) :: path
    logical :: upstream_index

    if (.not. filtered%ok) then
      output%message = 'The supplied filtering result is invalid.'
      return
    end if
    horizon = 15
    simulations = 1000
    conf = 0.95_dp
    if (present(n_ahead)) horizon = n_ahead
    if (present(n_sim)) simulations = n_sim
    if (present(confidence)) conf = confidence
    if (present(new_factors)) horizon = size(new_factors,1)
    if (horizon < 1 .or. simulations < 2) then
      output%message = 'The forecast horizon must be positive and n_sim at least two.'
      return
    end if
    if (conf <= 0.0_dp .or. conf >= 1.0_dp) then
      output%message = 'confidence must be between zero and one.'
      return
    end if
    if (present(new_factors)) then
      if (.not. allocated(filtered%dynamics%coefs)) then
        output%message = 'New factors were supplied but the model has no coefficients.'
        return
      end if
      if (size(new_factors,2) /= size(filtered%dynamics%coefs)) then
        output%message = 'new_factors must be horizon by number of coefficients.'
        return
      end if
    end if
    if (present(seed)) call seed_random(seed)
    upstream_index = .false.
    if (present(upstream_filter_index)) upstream_index = upstream_filter_index
    filter_column = size(filtered%filter_grid,2)
    if (upstream_index) filter_column = max(1,filter_column-1)
    probabilities = filtered%filter_grid(:,filter_column)
    allocate(vol_paths(horizon,simulations),return_paths(horizon,simulations))
    do j = 1, simulations
      index_value = sample_discrete(probabilities)
      path = model_simulate(filtered%dynamics,horizon, &
        filtered%grids%var_mid_points(index_value))
      if (.not. path%ok) then
        output%message = path%message
        return
      end if
      vol_paths(:,j) = path%volatility_factor
      return_paths(:,j) = path%returns
      if (present(new_factors)) then
        return_paths(:,j) = return_paths(:,j) + matmul(new_factors,filtered%dynamics%coefs)
      end if
    end do
    allocate(output%mean_volatility(horizon),output%lower_volatility(horizon), &
      output%upper_volatility(horizon),output%mean_return(horizon), &
      output%lower_return(horizon),output%upper_return(horizon))
    z = normal_quantile(0.5_dp+0.5_dp*conf)
    do i = 1, horizon
      output%mean_volatility(i)=mean_value(vol_paths(i,:))
      output%lower_volatility(i)=output%mean_volatility(i)-z*standard_deviation(vol_paths(i,:))
      output%upper_volatility(i)=output%mean_volatility(i)+z*standard_deviation(vol_paths(i,:))
      output%mean_return(i)=mean_value(return_paths(i,:))
      output%lower_return(i)=output%mean_return(i)-z*standard_deviation(return_paths(i,:))
      output%upper_return(i)=output%mean_return(i)+z*standard_deviation(return_paths(i,:))
    end do
    output%confidence = conf
    output%ok = .true.
  end function predict_filter

  function predict_svdnf(filtered, n_ahead, n_sim, confidence, new_data, seed) result(output)
    type(filter_result), intent(in) :: filtered
    integer, intent(in), optional :: n_ahead, n_sim
    real(dp), intent(in), optional :: confidence
    real(dp), intent(in), optional :: new_data(:,:)
    integer, intent(in), optional :: seed
    type(forecast_result) :: output
    if (present(new_data)) then
      output = predict_filter(filtered,n_ahead,n_sim,confidence,new_data,seed)
    else
      output = predict_filter(filtered,n_ahead,n_sim,confidence,seed=seed)
    end if
  end function predict_svdnf

end module svdnf_simulation
