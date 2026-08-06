program test_trading
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use trading
  implicit none

  integer :: failures

  failures = 0
  call test_statistics(failures)
  call test_curve_and_climate(failures)
  call test_trades_and_io(failures)
  call test_lottery(failures)
  call test_betting(failures)
  call test_dynamic_beta(failures)

  if (failures > 0) then
    write(*, '(a,i0)') "FAILED tests: ", failures
    error stop 1
  end if
  write(*, '(a)') "All Trading tests passed."

contains

  subroutine test_statistics(failures)
    integer, intent(inout) :: failures
    real(dp) :: data(6, 2)
    real(dp) :: distance(2, 2)
    real(dp) :: entropy

    data(:, 1) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]
    data(:, 2) = 2.0_dp * data(:, 1)
    call angular_distance(data, distance)
    call assert_close(distance(1, 2), 0.0_dp, 1.0e-12_dp, &
      "perfect angular distance", failures)
    call assert_close(chebyshev_distance([1.0_dp, 3.0_dp], &
      [2.0_dp, 1.0_dp]), 2.0_dp, 1.0e-12_dp, &
      "chebyshev distance", failures)

    entropy = sample_entropy([0.0_dp, 0.1_dp, -0.1_dp, 0.2_dp, -0.2_dp, &
      0.15_dp, -0.05_dp, 0.12_dp, -0.08_dp, 0.04_dp], tolerance=1.0_dp)
    call assert_true(ieee_is_finite(entropy), "sample entropy finite", failures)
    call assert_close(quantile_type7([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], &
      0.25_dp), 1.75_dp, 1.0e-12_dp, "type 7 quantile", failures)
  end subroutine test_statistics

  subroutine test_curve_and_climate(failures)
    integer, intent(inout) :: failures
    type(curve_t) :: curve
    real(dp) :: output(3)

    allocate(curve%tenors(3), curve%rates(3))
    curve%tenors = [0.0_dp, 1.0_dp, 2.0_dp]
    curve%rates = [0.0_dp, 2.0_dp, 4.0_dp]
    call curve%interpolate([0.5_dp, 1.5_dp, 3.0_dp], output, "linear")
    call assert_close(output(1), 1.0_dp, 1.0e-12_dp, &
      "linear interpolation", failures)
    call assert_close(output(3), 6.0_dp, 1.0e-12_dp, &
      "linear extrapolation", failures)

    call assert_close(total_carbon_emissions(&
      [100.0_dp, 200.0_dp, 50.0_dp], [1000.0_dp, 5000.0_dp, 6000.0_dp], &
      [20000.0_dp, 10000.0_dp, 30000.0_dp]), 115.0_dp, 1.0e-12_dp, &
      "total carbon emissions", failures)
    call assert_close(carbon_footprint(&
      [100.0_dp, 200.0_dp, 50.0_dp], [1000.0_dp, 5000.0_dp, 6000.0_dp], &
      [20000.0_dp, 10000.0_dp, 30000.0_dp]), 115.0_dp / 350.0_dp, &
      1.0e-12_dp, "carbon footprint", failures)
  end subroutine test_curve_and_climate

  subroutine test_trades_and_io(failures)
    integer, intent(inout) :: failures
    type(trade_t) :: swap
    type(trade_t) :: cdo
    type(trade_t), allocatable :: trades(:)
    type(trade_t), allocatable :: derivatives(:)
    type(curve_t) :: loaded_curve
    type(csa_t), allocatable :: agreements(:)
    type(collateral_t), allocatable :: collateral(:)
    type(hash_table_t) :: ratings
    real(dp), allocatable :: fund(:)
    real(dp), allocatable :: benchmark(:)
    logical :: found

    call swap%configure_class("IRDSwap")
    swap%notional = 10000.0_dp
    swap%si = 0.0_dp
    swap%ei = 10.0_dp
    swap%buy_sell = "Buy"
    call assert_close(swap%calc_supervisory_duration(), &
      (1.0_dp - exp(-0.5_dp)) / 0.05_dp, 1.0e-12_dp, &
      "supervisory duration", failures)
    call assert_true(swap%set_time_bucket() == 3, "time bucket", failures)

    call cdo%configure_class("CDOTranche")
    cdo%buy_sell = "Sell"
    cdo%cdo_attach_point = 0.3_dp
    cdo%cdo_detach_point = 0.5_dp
    call assert_true(cdo%cdo_supervisory_delta() < 0.0_dp, &
      "cdo sell delta", failures)

    call parse_trades_csv("data/example_trades.csv", trades)
    call assert_true(size(trades) == 29, "trade parser and digital expansion", failures)
    call select_derivatives(trades, derivatives)
    call assert_true(size(derivatives) > 0, "select derivatives", failures)

    call load_curve_csv("data/spot_rates.csv", loaded_curve)
    call assert_true(size(loaded_curve%tenors) > 10, "load curve csv", failures)
    call load_csa_csv("data/CSA.csv", agreements)
    call assert_true(size(agreements) == 2, "load csa csv", failures)
    call assert_true(size(agreements(1)%currencies) == 2, &
      "split csa currencies", failures)
    call load_collateral_csv("data/coll.csv", collateral)
    call assert_true(size(collateral) == 2, "load collateral csv", failures)
    call ratings%load_csv("data/RatingsMapping.csv")
    call assert_close(ratings%find_value_character("AAA", found), 0.007_dp, &
      1.0e-12_dp, "hash table percentage", failures)
    call assert_true(found, "hash table lookup", failures)
    call load_track_record_csv("data/example_track_record.csv", fund, benchmark)
    call assert_true(size(fund) == size(benchmark) .and. size(fund) > 50, &
      "load track record csv", failures)
  end subroutine test_trades_and_io

  subroutine test_lottery(failures)
    integer, intent(inout) :: failures
    type(lottery_draw_t) :: draw(2)
    type(lottery_draw_t), allocatable :: loaded(:)
    type(lottery_pnl_result_t) :: pnl
    type(euro_combination_iterator_t) :: iterator
    integer :: combination(7)
    logical :: has_value

    draw(1)%date_yyyymmdd = 20240101
    draw(1)%n_main = 5
    draw(1)%n_bonus = 2
    draw(1)%main_numbers(:5) = [1, 2, 3, 4, 5]
    draw(1)%bonus_numbers = [6, 7]
    draw(2) = draw(1)
    draw(2)%date_yyyymmdd = 20240102
    draw(2)%main_numbers(:5) = [10, 11, 12, 13, 14]

    call euro_lottery_backtest(draw, [1, 2, 3, 4, 5, 6, 7])
    call assert_true(draw(1)%winning_scale == 13, "euro jackpot scale", failures)
    call calculate_euro_lottery_pnl(draw, pnl)
    call assert_close(pnl%payout(1), 10000000.0_dp, 1.0e-6_dp, &
      "euro jackpot payout", failures)

    call load_euro_lottery_results("data/euromillions_results.csv", loaded)
    call assert_true(size(loaded) > 1000, "load euro results", failures)
    call assert_true(loaded(1)%date_yyyymmdd > 0, "parse textual date", failures)

    call iterator%next(combination, has_value)
    call assert_true(has_value .and. all(combination == [1, 2, 3, 4, 5, 1, 2]), &
      "first euro combination", failures)
    call iterator%next(combination, has_value)
    call assert_true(has_value .and. all(combination == [1, 2, 3, 4, 5, 1, 3]), &
      "second euro combination", failures)
    call assert_true(euro_lottery_combination_count() == 139838160_8, &
      "euro combination count", failures)
  end subroutine test_lottery

  subroutine test_betting(failures)
    integer, intent(inout) :: failures
    type(betting_result_t) :: betting
    type(repetitions_result_t) :: repetitions
    integer, allocatable :: fibonacci(:)

    call capped_fibonacci_sequence(20, fibonacci)
    call assert_true(all(fibonacci == [0, 1, 1, 2, 3, 5, 8, 13]), &
      "capped fibonacci", failures)

    call martingale_strategy_repetitions(4, 0.5_dp, 20, 100, repetitions, &
      quantile_probability=0.5_dp, seed=123)
    call assert_true(size(repetitions%number_of_trials_needed) == 20, &
      "martingale repetitions shape", failures)
    call assert_true(repetitions%has_quantile, "martingale quantile", failures)

    call roulette_dalembert(1.0_dp, 32.0_dp, 100.0_dp, 10, 50, &
      betting, seed=42)
    call assert_true(size(betting%final_capital) == 10, &
      "roulette result shape", failures)
    call assert_true(all(ieee_is_finite(betting%final_capital)), &
      "roulette finite", failures)
  end subroutine test_betting

  subroutine test_dynamic_beta(failures)
    integer, intent(inout) :: failures
    type(dynamic_beta_result) :: result
    real(dp) :: benchmark(80)
    real(dp) :: fund(80)
    real(dp) :: beta_mean
    integer :: i

    do i = 1, size(benchmark)
      benchmark(i) = 0.01_dp * sin(0.23_dp * real(i, dp)) + &
        0.004_dp * cos(0.11_dp * real(i, dp))
      fund(i) = 0.002_dp + 1.5_dp * benchmark(i) + &
        0.0002_dp * sin(0.71_dp * real(i, dp))
    end do

    call dynamic_beta(fund, benchmark, result, em_iterations=50)
    beta_mean = sum(result%smoothed_state(1, 21:80)) / 60.0_dp
    call assert_close(beta_mean, 1.5_dp, 0.15_dp, &
      "dynamic beta recovery", failures)
    call assert_true(result%observation_variance > 0.0_dp, &
      "dynamic beta observation variance", failures)
  end subroutine test_dynamic_beta

  subroutine assert_close(actual, expected, tolerance, name, failures)
    real(dp), intent(in) :: actual
    real(dp), intent(in) :: expected
    real(dp), intent(in) :: tolerance
    character(len=*), intent(in) :: name
    integer, intent(inout) :: failures

    if (.not. ieee_is_finite(actual) .or. abs(actual - expected) > tolerance) then
      failures = failures + 1
      write(*, '(a,2(1x,es16.8))') "FAIL " // trim(name) // ":", actual, expected
    end if
  end subroutine assert_close

  subroutine assert_true(condition, name, failures)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: name
    integer, intent(inout) :: failures

    if (.not. condition) then
      failures = failures + 1
      write(*, '(a)') "FAIL " // trim(name)
    end if
  end subroutine assert_true

end program test_trading
