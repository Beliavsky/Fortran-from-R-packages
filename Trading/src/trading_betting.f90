module trading_betting
  use trading_kinds, only : dp
  use trading_stats, only : quantile_type7, set_random_seed
  implicit none
  private

  type, public :: betting_result_t
    real(dp), allocatable :: maximum_capital(:)
    real(dp), allocatable :: minimum_capital(:)
    real(dp), allocatable :: final_capital(:)
    real(dp), allocatable :: last_path(:)
  end type betting_result_t

  type, public :: repetitions_result_t
    integer, allocatable :: number_of_trials_needed(:)
    real(dp) :: relevant_quantile = 0.0_dp
    logical :: has_quantile = .false.
  end type repetitions_result_t

  public :: capped_fibonacci_sequence
  public :: martingale_strategy_repetitions
  public :: roulette_dalembert
  public :: roulette_fibonacci
  public :: roulette_labouchere
  public :: roulette_martingale
  public :: roulette_specific_number

contains

  subroutine capped_fibonacci_sequence(maximum_number, sequence)
    integer, intent(in) :: maximum_number
    integer, allocatable, intent(out) :: sequence(:)
    integer, allocatable :: work(:)
    integer :: count
    integer :: next_value

    if (maximum_number <= 1) then
      allocate(sequence(2))
      sequence = [0, 1]
      return
    end if

    allocate(work(128))
    work = 0
    work(1:2) = [0, 1]
    count = 2
    do
      next_value = work(count) + work(count - 1)
      if (next_value >= maximum_number) exit
      count = count + 1
      if (count > size(work)) error stop "capped_fibonacci_sequence: limit exceeded"
      work(count) = next_value
    end do

    allocate(sequence(count))
    sequence = work(:count)
  end subroutine capped_fibonacci_sequence

  subroutine martingale_strategy_repetitions(length_of_targeted_sequence, &
      probability_of_success, simulations_number, trials_per_simulation, &
      result, quantile_probability, seed)
    integer, intent(in) :: length_of_targeted_sequence
    real(dp), intent(in) :: probability_of_success
    integer, intent(in) :: simulations_number
    integer, intent(in) :: trials_per_simulation
    type(repetitions_result_t), intent(out) :: result
    real(dp), intent(in), optional :: quantile_probability
    integer, intent(in), optional :: seed
    logical, allocatable :: outcomes(:)
    real(dp), allocatable :: quantile_data(:)
    real(dp), allocatable :: random_values(:)
    integer :: i
    integer :: j
    integer :: run_length

    call validate_simulation_sizes(simulations_number, trials_per_simulation)
    if (length_of_targeted_sequence < 1) then
      error stop "martingale_strategy_repetitions: target length must be positive"
    end if
    if (probability_of_success < 0.0_dp .or. probability_of_success > 1.0_dp) then
      error stop "martingale_strategy_repetitions: probability must be in [0,1]"
    end if
    if (present(seed)) call set_random_seed(seed)

    allocate(result%number_of_trials_needed(simulations_number))
    allocate(outcomes(trials_per_simulation), random_values(trials_per_simulation))
    result%number_of_trials_needed = 0

    do j = 1, simulations_number
      call random_number(random_values)
      outcomes = random_values <= probability_of_success
      run_length = 1
      if (length_of_targeted_sequence == 1) then
        result%number_of_trials_needed(j) = 1
        cycle
      end if
      do i = 1, trials_per_simulation - 1
        if (outcomes(i) .eqv. outcomes(i + 1)) then
          run_length = run_length + 1
          if (run_length == length_of_targeted_sequence) then
            result%number_of_trials_needed(j) = i + 1
            exit
          end if
        else
          run_length = 1
        end if
      end do
    end do

    if (present(quantile_probability)) then
      allocate(quantile_data(simulations_number))
      quantile_data = real(result%number_of_trials_needed, dp)
      do i = 1, simulations_number
        if (result%number_of_trials_needed(i) == 0) then
          quantile_data(i) = real(trials_per_simulation, dp)
        end if
      end do
      result%relevant_quantile = quantile_type7(quantile_data, quantile_probability)
      result%has_quantile = .true.
    end if
  end subroutine martingale_strategy_repetitions

  subroutine roulette_dalembert(bet_minimum, bet_maximum, initial_capital, &
      simulations_number, trials_per_simulation, result, seed)
    real(dp), intent(in) :: bet_minimum
    real(dp), intent(in) :: bet_maximum
    real(dp), intent(in) :: initial_capital
    integer, intent(in) :: simulations_number
    integer, intent(in) :: trials_per_simulation
    type(betting_result_t), intent(out) :: result
    integer, intent(in), optional :: seed
    integer, allocatable :: spins(:)
    real(dp) :: bet_amount
    real(dp) :: capital
    integer :: i
    integer :: j
    integer :: path_length

    call initialize_betting_result(simulations_number, trials_per_simulation - 1, result)
    call validate_bets(bet_minimum, bet_maximum)
    if (present(seed)) call set_random_seed(seed)
    allocate(spins(trials_per_simulation))

    do j = 1, simulations_number
      call generate_spins(spins)
      bet_amount = bet_minimum
      capital = initial_capital
      path_length = 0

      do i = 1, trials_per_simulation - 1
        if (spins(i + 1) == 0) then
          capital = capital - bet_amount
          bet_amount = min(bet_amount + bet_minimum, bet_maximum)
        else if (mod(spins(i), 2) == mod(spins(i + 1), 2)) then
          capital = capital + bet_amount
          if (bet_amount - bet_minimum >= bet_minimum) then
            bet_amount = bet_amount - bet_minimum
          end if
        else
          capital = capital - bet_amount
          if (bet_amount + bet_minimum > bet_maximum) then
            bet_amount = bet_minimum
          else
            bet_amount = bet_amount + bet_minimum
          end if
        end if
        path_length = i
        result%last_path(i) = capital
        if (capital < 0.0_dp) exit
      end do
      call store_summary(result, j, initial_capital, capital, path_length)
    end do
  end subroutine roulette_dalembert

  subroutine roulette_fibonacci(bet_minimum, bet_maximum, initial_capital, &
      simulations_number, trials_per_simulation, result, seed)
    real(dp), intent(in) :: bet_minimum
    real(dp), intent(in) :: bet_maximum
    real(dp), intent(in) :: initial_capital
    integer, intent(in) :: simulations_number
    integer, intent(in) :: trials_per_simulation
    type(betting_result_t), intent(out) :: result
    integer, intent(in), optional :: seed
    integer, allocatable :: fibonacci(:)
    integer, allocatable :: spins(:)
    real(dp) :: capital
    real(dp) :: wager
    integer :: counter
    integer :: i
    integer :: j
    integer :: path_length

    call initialize_betting_result(simulations_number, trials_per_simulation - 1, result)
    call validate_bets(bet_minimum, bet_maximum)
    if (present(seed)) call set_random_seed(seed)
    call capped_fibonacci_sequence(int(bet_maximum / bet_minimum), fibonacci)
    allocate(spins(trials_per_simulation))

    do j = 1, simulations_number
      call generate_spins(spins)
      capital = initial_capital
      counter = min(2, size(fibonacci))
      path_length = 0

      do i = 1, trials_per_simulation - 1
        wager = bet_minimum * real(fibonacci(counter), dp)
        if (spins(i + 1) == 0) then
          capital = capital - wager
          counter = counter + 1
          if (counter > size(fibonacci)) counter = min(2, size(fibonacci))
        else if (mod(spins(i), 2) == mod(spins(i + 1), 2)) then
          capital = capital + wager
          counter = max(counter - 2, min(2, size(fibonacci)))
        else
          capital = capital - wager
          counter = counter + 1
          if (counter > size(fibonacci)) counter = min(2, size(fibonacci))
        end if
        path_length = i
        result%last_path(i) = capital
        if (capital < 0.0_dp) exit
      end do
      call store_summary(result, j, initial_capital, capital, path_length)
    end do
  end subroutine roulette_fibonacci

  subroutine roulette_labouchere(bet_minimum, bet_maximum, initial_capital, &
      profit_target, simulations_number, trials_per_simulation, result, &
      profit_sequence, seed)
    real(dp), intent(in) :: bet_minimum
    real(dp), intent(in) :: bet_maximum
    real(dp), intent(in) :: initial_capital
    real(dp), intent(in) :: profit_target
    integer, intent(in) :: simulations_number
    integer, intent(in) :: trials_per_simulation
    type(betting_result_t), intent(out) :: result
    real(dp), intent(in), optional :: profit_sequence(:)
    integer, intent(in), optional :: seed
    integer, allocatable :: spins(:)
    real(dp), allocatable :: base_sequence(:)
    real(dp), allocatable :: work_sequence(:)
    real(dp) :: bet_amount
    real(dp) :: capital
    integer :: i
    integer :: j
    integer :: path_length
    integer :: sequence_length

    call initialize_betting_result(simulations_number, trials_per_simulation - 1, result)
    call validate_bets(bet_minimum, bet_maximum)
    if (present(seed)) call set_random_seed(seed)

    if (present(profit_sequence)) then
      allocate(base_sequence(size(profit_sequence)))
      base_sequence = profit_sequence
    else
      sequence_length = max(1, nint(profit_target / bet_minimum))
      allocate(base_sequence(sequence_length))
      base_sequence = bet_minimum
    end if
    if (abs(sum(base_sequence) - profit_target) > &
        1.0e-10_dp * max(1.0_dp, abs(profit_target))) then
      error stop "roulette_labouchere: sequence must sum to the profit target"
    end if

    allocate(work_sequence(size(base_sequence) + trials_per_simulation))
    allocate(spins(trials_per_simulation))

    do j = 1, simulations_number
      call generate_spins(spins)
      work_sequence = 0.0_dp
      work_sequence(:size(base_sequence)) = base_sequence
      sequence_length = size(base_sequence)
      capital = initial_capital
      path_length = 0
      bet_amount = labouchere_bet(work_sequence, sequence_length, bet_maximum)

      do i = 1, trials_per_simulation - 1
        if (spins(i + 1) == 0 .or. &
            mod(spins(i), 2) /= mod(spins(i + 1), 2)) then
          capital = capital - bet_amount
          sequence_length = sequence_length + 1
          work_sequence(sequence_length) = bet_amount
        else
          capital = capital + bet_amount
          if (sequence_length <= 2) then
            path_length = i
            result%last_path(i) = capital
            exit
          end if
          work_sequence(1:sequence_length - 2) = work_sequence(2:sequence_length - 1)
          sequence_length = sequence_length - 2
        end if

        path_length = i
        result%last_path(i) = capital
        if (capital < 0.0_dp .or. sequence_length == 0) exit
        bet_amount = labouchere_bet(work_sequence, sequence_length, bet_maximum)
      end do
      call store_summary(result, j, initial_capital, capital, path_length)
    end do
  end subroutine roulette_labouchere

  subroutine roulette_martingale(bet_minimum, bet_maximum, initial_capital, &
      simulations_number, trials_per_simulation, result, seed)
    real(dp), intent(in) :: bet_minimum
    real(dp), intent(in) :: bet_maximum
    real(dp), intent(in) :: initial_capital
    integer, intent(in) :: simulations_number
    integer, intent(in) :: trials_per_simulation
    type(betting_result_t), intent(out) :: result
    integer, intent(in), optional :: seed
    integer, allocatable :: spins(:)
    real(dp) :: bet_amount
    real(dp) :: capital
    logical :: first_loss
    logical :: win
    integer :: i
    integer :: j
    integer :: path_length
    integer :: previous_nonzero

    call initialize_betting_result(simulations_number, trials_per_simulation - 1, result)
    call validate_bets(bet_minimum, bet_maximum)
    if (present(seed)) call set_random_seed(seed)
    allocate(spins(trials_per_simulation))

    do j = 1, simulations_number
      call generate_spins(spins)
      bet_amount = bet_minimum
      capital = initial_capital
      first_loss = .true.
      path_length = 0

      do i = 1, trials_per_simulation - 1
        if (spins(i) == 0 .and. first_loss) then
          capital = capital - bet_amount
          path_length = i
          result%last_path(i) = capital
          cycle
        end if

        if (spins(i + 1) == 0) then
          win = .false.
        else if (spins(i) == 0) then
          previous_nonzero = find_previous_nonzero(spins, i - 1)
          if (previous_nonzero == 0) then
            win = .false.
          else
            win = mod(spins(previous_nonzero), 2) /= mod(spins(i + 1), 2)
          end if
        else
          win = mod(spins(i), 2) /= mod(spins(i + 1), 2)
        end if

        if (win) then
          capital = capital + bet_amount
          bet_amount = bet_minimum
          first_loss = .true.
        else
          capital = capital - bet_amount
          if (bet_amount >= bet_maximum .or. &
              (abs(bet_amount - bet_minimum) <= epsilon(bet_amount) .and. first_loss)) then
            bet_amount = bet_minimum / 2.0_dp
            first_loss = .false.
          end if
          bet_amount = min(2.0_dp * bet_amount, bet_maximum)
        end if

        path_length = i
        result%last_path(i) = capital
        if (capital < 0.0_dp) exit
      end do
      call store_summary(result, j, initial_capital, capital, path_length)
    end do
  end subroutine roulette_martingale

  subroutine roulette_specific_number(bet_minimum, bet_maximum, initial_capital, &
      targeted_number, simulations_number, trials_per_simulation, result, &
      stop_loss, seed)
    real(dp), intent(in) :: bet_minimum
    real(dp), intent(in) :: bet_maximum
    real(dp), intent(in) :: initial_capital
    integer, intent(in) :: targeted_number
    integer, intent(in) :: simulations_number
    integer, intent(in) :: trials_per_simulation
    type(betting_result_t), intent(out) :: result
    integer, intent(in), optional :: stop_loss
    integer, intent(in), optional :: seed
    integer, allocatable :: spins(:)
    real(dp) :: bet_amount
    real(dp) :: capital
    integer :: i
    integer :: j
    integer :: loss_counter
    integer :: path_length

    call initialize_betting_result(simulations_number, trials_per_simulation, result)
    call validate_bets(bet_minimum, bet_maximum)
    if (targeted_number < 0 .or. targeted_number > 37) then
      error stop "roulette_specific_number: target must be in 0:37"
    end if
    if (present(seed)) call set_random_seed(seed)
    allocate(spins(trials_per_simulation))

    do j = 1, simulations_number
      call generate_spins(spins)
      bet_amount = bet_minimum
      capital = initial_capital
      loss_counter = 0
      path_length = 0

      do i = 1, trials_per_simulation
        if (spins(i) /= targeted_number) then
          capital = capital - bet_amount
          loss_counter = loss_counter + 1
          if (mod(loss_counter, 18) == 0) then
            bet_amount = min(2.0_dp * bet_amount, bet_maximum)
          end if
          if (present(stop_loss)) then
            if (loss_counter > stop_loss) bet_amount = bet_minimum
          end if
        else
          capital = capital + 35.0_dp * bet_amount
          bet_amount = bet_minimum
          loss_counter = 0
        end if
        path_length = i
        result%last_path(i) = capital
        if (capital < 0.0_dp) exit
      end do
      call store_summary(result, j, initial_capital, capital, path_length)
    end do
  end subroutine roulette_specific_number

  subroutine initialize_betting_result(simulations_number, path_size, result)
    integer, intent(in) :: simulations_number
    integer, intent(in) :: path_size
    type(betting_result_t), intent(out) :: result

    call validate_simulation_sizes(simulations_number, max(path_size, 1))
    allocate(result%maximum_capital(simulations_number))
    allocate(result%minimum_capital(simulations_number))
    allocate(result%final_capital(simulations_number))
    allocate(result%last_path(max(path_size, 1)))
    result%maximum_capital = 0.0_dp
    result%minimum_capital = 0.0_dp
    result%final_capital = 0.0_dp
    result%last_path = 0.0_dp
  end subroutine initialize_betting_result

  subroutine store_summary(result, simulation_index, initial_capital, capital, path_length)
    type(betting_result_t), intent(inout) :: result
    integer, intent(in) :: simulation_index
    real(dp), intent(in) :: initial_capital
    real(dp), intent(in) :: capital
    integer, intent(in) :: path_length

    result%final_capital(simulation_index) = capital
    if (path_length > 0) then
      result%minimum_capital(simulation_index) = &
        min(initial_capital, minval(result%last_path(:path_length)))
      result%maximum_capital(simulation_index) = &
        max(initial_capital, maxval(result%last_path(:path_length)))
    else
      result%minimum_capital(simulation_index) = initial_capital
      result%maximum_capital(simulation_index) = initial_capital
    end if
  end subroutine store_summary

  subroutine generate_spins(spins)
    integer, intent(out) :: spins(:)
    real(dp), allocatable :: random_values(:)

    allocate(random_values(size(spins)))
    call random_number(random_values)
    spins = floor(38.0_dp * random_values)
  end subroutine generate_spins

  pure real(dp) function labouchere_bet(sequence, sequence_length, bet_maximum) &
      result(value)
    real(dp), intent(in) :: sequence(:)
    integer, intent(in) :: sequence_length
    real(dp), intent(in) :: bet_maximum

    if (sequence_length <= 0) then
      value = 0.0_dp
    else if (sequence_length == 1) then
      value = min(sequence(1), bet_maximum)
    else
      value = min(sequence(1) + sequence(sequence_length), bet_maximum)
    end if
  end function labouchere_bet

  pure integer function find_previous_nonzero(spins, last_index) result(index_value)
    integer, intent(in) :: spins(:)
    integer, intent(in) :: last_index
    integer :: i

    index_value = 0
    do i = last_index, 1, -1
      if (spins(i) /= 0) then
        index_value = i
        return
      end if
    end do
  end function find_previous_nonzero

  subroutine validate_bets(bet_minimum, bet_maximum)
    real(dp), intent(in) :: bet_minimum
    real(dp), intent(in) :: bet_maximum

    if (bet_minimum <= 0.0_dp .or. bet_maximum < bet_minimum) then
      error stop "roulette simulation: invalid betting limits"
    end if
  end subroutine validate_bets

  subroutine validate_simulation_sizes(simulations_number, trials_per_simulation)
    integer, intent(in) :: simulations_number
    integer, intent(in) :: trials_per_simulation

    if (simulations_number < 1 .or. trials_per_simulation < 1) then
      error stop "simulation dimensions must be positive"
    end if
  end subroutine validate_simulation_sizes

end module trading_betting
