module xva_exposure
  use trading, only : dp, trade_t, csa_t, quantile_type7, set_random_seed, &
    random_normal
  use xva_types, only : simulation_data_t, exposure_profile_t
  use xva_math, only : same_text
  implicit none
  private

  public :: calc_simulated_exposure

contains

  subroutine calc_simulated_exposure(discount_factors, time_points, spot_curve, &
      csa, trades, sim_data, framework, exposure_profile, seed)
    real(dp), intent(in) :: discount_factors(:)
    real(dp), intent(in) :: time_points(:)
    real(dp), intent(in) :: spot_curve(:)
    type(csa_t), intent(in) :: csa
    type(trade_t), intent(inout) :: trades(:)
    type(simulation_data_t), intent(in) :: sim_data
    character(len=*), intent(in) :: framework
    type(exposure_profile_t), intent(out) :: exposure_profile
    integer, intent(in), optional :: seed
    real(dp), allocatable :: forward_curve(:)
    real(dp), allocatable :: theta(:)
    real(dp), allocatable :: time_differences(:)
    real(dp), allocatable :: normal_draws(:,:)
    real(dp), allocatable :: all_draws(:)
    real(dp), allocatable :: swap_mtms(:,:)
    real(dp), allocatable :: collateralized_swap_mtms(:,:)
    real(dp) :: spot_interest_rate
    integer :: i
    integer :: j
    integer :: n
    integer :: number_of_draws

    n = size(time_points)
    if (n < 2) error stop "calc_simulated_exposure: at least two time points are required"
    if (size(discount_factors) /= n .or. size(spot_curve) /= n) then
      error stop "calc_simulated_exposure: curve and time-grid sizes differ"
    end if
    if (any(time_points(2:) <= time_points(:n - 1))) then
      error stop "calc_simulated_exposure: time_points must be strictly increasing"
    end if
    if (sim_data%num_simulations <= 0) then
      error stop "calc_simulated_exposure: num_simulations must be positive"
    end if
    if (sim_data%mean_reversion_a <= 0.0_dp) then
      error stop "calc_simulated_exposure: mean_reversion_a must be positive"
    end if
    if (sim_data%volatility < 0.0_dp) then
      error stop "calc_simulated_exposure: volatility must be nonnegative"
    end if
    if (sim_data%pfe_percentile < 0.0_dp .or. &
        sim_data%pfe_percentile > 1.0_dp) then
      error stop "calc_simulated_exposure: pfe_percentile must be in [0,1]"
    end if
    if (any(discount_factors <= 0.0_dp)) then
      error stop "calc_simulated_exposure: discount factors must be positive"
    end if

    allocate(time_differences(n - 1), forward_curve(n), theta(n - 1))
    time_differences = time_points(2:) - time_points(:n - 1)
    spot_interest_rate = spot_curve(1)
    forward_curve(1) = spot_interest_rate
    forward_curve(2:) = (discount_factors(:n - 1) / discount_factors(2:) - &
      1.0_dp) / time_differences
    theta = forward_curve(2:) - forward_curve(:n - 1) + &
      sim_data%mean_reversion_a * forward_curve(2:)

    if (present(seed)) call set_random_seed(seed)
    number_of_draws = sim_data%num_simulations * (n - 1)
    allocate(all_draws(number_of_draws), &
      normal_draws(sim_data%num_simulations, n - 1))
    call random_normal(all_draws)
    normal_draws = reshape(all_draws, shape(normal_draws))

    allocate(swap_mtms(sim_data%num_simulations, n), &
      collateralized_swap_mtms(sim_data%num_simulations, n))
    swap_mtms = 0.0_dp
    collateralized_swap_mtms = 0.0_dp

    do i = 1, size(trades)
      if (.not. same_text(trades(i)%class_name, "IRDSwap")) then
        error stop "calc_simulated_exposure: only IRDSwap trades are supported"
      end if
      call add_swap_paths(trades(i), discount_factors, time_points, forward_curve, &
        theta, time_differences, normal_draws, sim_data, csa, swap_mtms, &
        collateralized_swap_mtms, same_text(framework, "IMM"))
    end do

    call summarize_paths(swap_mtms, sim_data%pfe_percentile, &
      exposure_profile%ee_uncollateralized, &
      exposure_profile%nee_uncollateralized, &
      exposure_profile%pfe_uncollateralized)
    call summarize_paths(collateralized_swap_mtms, sim_data%pfe_percentile, &
      exposure_profile%ee, exposure_profile%nee, exposure_profile%pfe)

    allocate(exposure_profile%eee(n))
    exposure_profile%eee(1) = exposure_profile%ee(1)
    do j = 2, n
      exposure_profile%eee(j) = max(exposure_profile%eee(j - 1), &
        exposure_profile%ee(j))
    end do
  end subroutine calc_simulated_exposure

  subroutine add_swap_paths(trade, discount_factors, time_points, forward_curve, &
      theta, time_differences, normal_draws, sim_data, csa, swap_mtms, &
      collateralized_swap_mtms, update_mtm)
    type(trade_t), intent(inout) :: trade
    real(dp), intent(in) :: discount_factors(:)
    real(dp), intent(in) :: time_points(:)
    real(dp), intent(in) :: forward_curve(:)
    real(dp), intent(in) :: theta(:)
    real(dp), intent(in) :: time_differences(:)
    real(dp), intent(in) :: normal_draws(:,:)
    type(simulation_data_t), intent(in) :: sim_data
    type(csa_t), intent(in) :: csa
    real(dp), intent(inout) :: swap_mtms(:,:)
    real(dp), intent(inout) :: collateralized_swap_mtms(:,:)
    logical, intent(in) :: update_mtm
    real(dp), allocatable :: a_factor(:)
    real(dp), allocatable :: b_factor(:)
    real(dp), allocatable :: disc_factors(:,:)
    real(dp), allocatable :: fixed_leg(:)
    real(dp), allocatable :: floating_leg(:)
    real(dp), allocatable :: interest_rates(:)
    real(dp), allocatable :: swap_mtm(:)
    real(dp), allocatable :: collateralized_swap_mtm(:)
    real(dp) :: buy_sell_sign
    real(dp) :: dt
    real(dp) :: maturity
    real(dp) :: variance_adjustment
    integer :: i
    integer :: j
    integer :: number_of_points

    maturity = trade%ei
    number_of_points = count(time_points <= maturity + &
      100.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(maturity)))
    if (number_of_points == 0) return

    allocate(a_factor(number_of_points), b_factor(number_of_points), &
      disc_factors(number_of_points, number_of_points), &
      fixed_leg(number_of_points), floating_leg(number_of_points), &
      interest_rates(number_of_points), swap_mtm(number_of_points), &
      collateralized_swap_mtm(number_of_points))

    b_factor = (1.0_dp - exp(-sim_data%mean_reversion_a * &
      (maturity - time_points(:number_of_points)))) / sim_data%mean_reversion_a
    variance_adjustment = sim_data%volatility**2 / &
      (4.0_dp * sim_data%mean_reversion_a)
    dt = maturity / real(number_of_points, dp)
    if (same_text(trade%buy_sell, "Buy")) then
      buy_sell_sign = 1.0_dp
    else
      buy_sell_sign = -1.0_dp
    end if

    do j = 1, sim_data%num_simulations
      interest_rates = 0.0_dp
      interest_rates(1) = forward_curve(1)
      do i = 2, number_of_points
        interest_rates(i) = interest_rates(i - 1) + &
          (theta(i - 1) - sim_data%mean_reversion_a * interest_rates(i - 1)) * &
          time_differences(i - 1) + sim_data%volatility * &
          normal_draws(j, i - 1) * sqrt(time_differences(i - 1))
      end do

      do i = 1, number_of_points
        a_factor(i) = discount_factors(number_of_points) / discount_factors(i) * &
          exp(b_factor(i) * forward_curve(i) - variance_adjustment * &
          (1.0_dp - exp(-2.0_dp * sim_data%mean_reversion_a * &
          time_points(i))) * b_factor(i)**2)
      end do

      disc_factors = 0.0_dp
      do i = 1, number_of_points
        disc_factors(i, :number_of_points - i + 1) = &
          a_factor(number_of_points - i + 1) * &
          exp(-b_factor(number_of_points - i + 1) * &
          interest_rates(:number_of_points - i + 1))
      end do

      do i = 1, number_of_points
        floating_leg(i) = 1.0_dp - &
          disc_factors(number_of_points - i + 1, i)
        fixed_leg(i) = dt * trade%pay_leg_rate * sum(disc_factors(:, i))
      end do
      swap_mtm = (floating_leg - fixed_leg) * buy_sell_sign
      swap_mtms(j, :number_of_points) = swap_mtms(j, :number_of_points) + swap_mtm
      call csa%apply_threshold(swap_mtm, collateralized_swap_mtm)
      collateralized_swap_mtms(j, :number_of_points) = &
        collateralized_swap_mtms(j, :number_of_points) + collateralized_swap_mtm

      if (update_mtm .and. j == 1) trade%mtm = swap_mtm(1)
    end do
  end subroutine add_swap_paths

  subroutine summarize_paths(paths, percentile, ee, nee, pfe)
    real(dp), intent(in) :: paths(:,:)
    real(dp), intent(in) :: percentile
    real(dp), allocatable, intent(out) :: ee(:)
    real(dp), allocatable, intent(out) :: nee(:)
    real(dp), allocatable, intent(out) :: pfe(:)
    real(dp), allocatable :: selected(:)
    integer :: count_selected
    integer :: i
    integer :: j
    integer :: n_times

    n_times = size(paths, 2)
    allocate(ee(n_times), nee(n_times), pfe(n_times), selected(size(paths, 1)))
    ee = 0.0_dp
    nee = 0.0_dp
    pfe = 0.0_dp

    do j = 1, n_times
      count_selected = 0
      do i = 1, size(paths, 1)
        if (paths(i, j) >= 0.0_dp) then
          count_selected = count_selected + 1
          selected(count_selected) = paths(i, j)
        end if
      end do
      if (count_selected > 0) ee(j) = sum(selected(:count_selected)) / &
        real(count_selected, dp)

      count_selected = 0
      do i = 1, size(paths, 1)
        if (paths(i, j) < 0.0_dp) then
          count_selected = count_selected + 1
          selected(count_selected) = paths(i, j)
        end if
      end do
      if (count_selected > 0) nee(j) = sum(selected(:count_selected)) / &
        real(count_selected, dp)

      count_selected = 0
      do i = 1, size(paths, 1)
        if (paths(i, j) > 0.0_dp) then
          count_selected = count_selected + 1
          selected(count_selected) = paths(i, j)
        end if
      end do
      if (count_selected > 0) then
        pfe(j) = quantile_type7(selected(:count_selected), percentile)
      end if
    end do
  end subroutine summarize_paths

end module xva_exposure
