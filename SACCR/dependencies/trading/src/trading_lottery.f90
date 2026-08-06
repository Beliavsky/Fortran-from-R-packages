module trading_lottery
  use trading_kinds, only : dp, str_len
  use trading_strings, only : lowercase, parse_integer_or_zero, split_delimited
  implicit none
  private

  type, public :: lottery_draw_t
    integer :: date_yyyymmdd = 0
    integer :: main_numbers(6) = 0
    integer :: bonus_numbers(2) = 0
    integer :: n_main = 0
    integer :: n_bonus = 0
    integer :: matched_main = 0
    integer :: matched_bonus = 0
    integer :: winning_scale = 0
  end type lottery_draw_t

  type, public :: lottery_pnl_result_t
    integer, allocatable :: date_yyyymmdd(:)
    real(dp), allocatable :: payout(:)
    real(dp), allocatable :: cumulative_profit(:)
    real(dp), allocatable :: cumulative_cost(:)
    real(dp), allocatable :: final_pnl(:)
  end type lottery_pnl_result_t

  type, public :: euro_combination_iterator_t
    integer :: main_numbers(5) = [1, 2, 3, 4, 5]
    integer :: bonus_numbers(2) = [1, 2]
    logical :: started = .false.
    logical :: finished = .false.
  contains
    procedure :: next => next_euro_combination
    procedure :: reset => reset_euro_combination
  end type euro_combination_iterator_t

  public :: load_euro_lottery_results
  public :: load_set_for_life_results
  public :: load_uk_lottery_results
  public :: load_uk_thunderball_results
  public :: euro_lottery_backtest
  public :: set_for_life_backtest
  public :: uk_lottery_backtest
  public :: uk_thunderball_backtest
  public :: calculate_euro_lottery_pnl
  public :: calculate_set_for_life_pnl
  public :: calculate_uk_lottery_pnl
  public :: calculate_uk_thunderball_pnl
  public :: top_five_numbers
  public :: euro_lottery_combination_count
  public :: outer_join_merge_integer
  public :: parse_date

contains

  subroutine load_euro_lottery_results(filename, draws)
    character(len=*), intent(in) :: filename
    type(lottery_draw_t), allocatable, intent(out) :: draws(:)

    call load_draw_file(filename, 5, 2, .true., draws)
  end subroutine load_euro_lottery_results

  subroutine load_set_for_life_results(filename, draws)
    character(len=*), intent(in) :: filename
    type(lottery_draw_t), allocatable, intent(out) :: draws(:)

    call load_draw_file(filename, 5, 1, .false., draws)
  end subroutine load_set_for_life_results

  subroutine load_uk_lottery_results(filename, draws)
    character(len=*), intent(in) :: filename
    type(lottery_draw_t), allocatable, intent(out) :: draws(:)

    call load_draw_file(filename, 6, 1, .false., draws)
  end subroutine load_uk_lottery_results

  subroutine load_uk_thunderball_results(filename, draws)
    character(len=*), intent(in) :: filename
    type(lottery_draw_t), allocatable, intent(out) :: draws(:)

    call load_draw_file(filename, 5, 1, .true., draws)
  end subroutine load_uk_thunderball_results

  subroutine load_draw_file(filename, n_main, n_bonus, textual_date, draws)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: n_main
    integer, intent(in) :: n_bonus
    logical, intent(in) :: textual_date
    type(lottery_draw_t), allocatable, intent(out) :: draws(:)
    character(len=2048) :: line
    character(len=str_len) :: fields(16)
    integer :: count
    integer :: field_count
    integer :: i
    integer :: ios
    integer :: row
    integer :: unit

    open(newunit=unit, file=filename, status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "load_draw_file: unable to open file"

    read(unit, '(a)', iostat=ios) line
    count = 0
    do
      read(unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) > 0) count = count + 1
    end do
    rewind(unit)
    read(unit, '(a)', iostat=ios) line

    allocate(draws(count))
    row = 0
    do
      read(unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) cycle
      call split_delimited(line, ',', fields, field_count)
      if (field_count < 1 + n_main + n_bonus) cycle

      row = row + 1
      draws(row)%n_main = n_main
      draws(row)%n_bonus = n_bonus
      draws(row)%date_yyyymmdd = parse_date(fields(1), .true., textual_date)
      do i = 1, n_main
        draws(row)%main_numbers(i) = parse_integer_or_zero(fields(i + 1))
      end do
      do i = 1, n_bonus
        draws(row)%bonus_numbers(i) = &
          parse_integer_or_zero(fields(1 + n_main + i))
      end do
    end do
    close(unit)

    if (row < count) draws = draws(:row)
  end subroutine load_draw_file

  integer function parse_date(text, day_first, textual) result(date_value)
    character(len=*), intent(in) :: text
    logical, intent(in), optional :: day_first
    logical, intent(in), optional :: textual
    character(len=str_len) :: fields(8)
    character(len=str_len) :: normalized
    logical :: is_day_first
    logical :: is_textual
    integer :: count
    integer :: day
    integer :: first
    integer :: month
    integer :: second
    integer :: year

    is_day_first = .true.
    if (present(day_first)) is_day_first = day_first
    is_textual = index(lowercase(text), "january") > 0 .or. &
      index(lowercase(text), "february") > 0 .or. &
      index(lowercase(text), "march") > 0 .or. &
      index(lowercase(text), "april") > 0 .or. &
      index(lowercase(text), "may") > 0 .or. &
      index(lowercase(text), "june") > 0 .or. &
      index(lowercase(text), "july") > 0 .or. &
      index(lowercase(text), "august") > 0 .or. &
      index(lowercase(text), "september") > 0 .or. &
      index(lowercase(text), "october") > 0 .or. &
      index(lowercase(text), "november") > 0 .or. &
      index(lowercase(text), "december") > 0
    if (present(textual)) is_textual = textual

    if (is_textual) then
      normalized = adjustl(trim(text))
      call replace_character(normalized, '-', ' ')
      call split_delimited(normalized, ' ', fields, count)
      call remove_empty_fields(fields, count)
      if (count < 3) then
        date_value = 0
        return
      end if

      year = parse_integer_or_zero(fields(count))
      month = month_number(fields(count - 1))
      day = parse_ordinal_day(fields(count - 2))
      date_value = 10000 * year + 100 * month + day
    else
      normalized = adjustl(trim(text))
      call replace_character(normalized, '-', '/')
      call split_delimited(normalized, '/', fields, count)
      if (count < 3) then
        date_value = 0
        return
      end if
      first = parse_integer_or_zero(fields(1))
      second = parse_integer_or_zero(fields(2))
      year = parse_integer_or_zero(fields(3))
      if (year < 100) year = year + 2000
      if (is_day_first) then
        day = first
        month = second
      else
        month = first
        day = second
      end if
      date_value = 10000 * year + 100 * month + day
    end if
  end function parse_date

  subroutine euro_lottery_backtest(draws, user_numbers, date_since)
    type(lottery_draw_t), intent(inout) :: draws(:)
    integer, intent(in) :: user_numbers(7)
    integer, intent(in), optional :: date_since
    integer :: i

    do i = 1, size(draws)
      if (present(date_since)) then
        if (draws(i)%date_yyyymmdd <= date_since) then
          draws(i)%winning_scale = 0
          cycle
        end if
      end if
      draws(i)%matched_main = match_count(user_numbers(:5), &
        draws(i)%main_numbers(:5))
      draws(i)%matched_bonus = match_count(user_numbers(6:7), &
        draws(i)%bonus_numbers(:2))
      draws(i)%winning_scale = euro_scale(draws(i)%matched_main, &
        draws(i)%matched_bonus)
    end do
  end subroutine euro_lottery_backtest

  subroutine set_for_life_backtest(draws, user_numbers, date_since)
    type(lottery_draw_t), intent(inout) :: draws(:)
    integer, intent(in) :: user_numbers(6)
    integer, intent(in), optional :: date_since
    integer :: i

    do i = 1, size(draws)
      if (present(date_since)) then
        if (draws(i)%date_yyyymmdd <= date_since) then
          draws(i)%winning_scale = 0
          cycle
        end if
      end if
      draws(i)%matched_main = match_count(user_numbers(:5), &
        draws(i)%main_numbers(:5))
      draws(i)%matched_bonus = match_count(user_numbers(6:6), &
        draws(i)%bonus_numbers(:1))
      draws(i)%winning_scale = set_for_life_scale(draws(i)%matched_main, &
        draws(i)%matched_bonus)
    end do
  end subroutine set_for_life_backtest

  subroutine uk_lottery_backtest(draws, user_numbers, date_since)
    type(lottery_draw_t), intent(inout) :: draws(:)
    integer, intent(in) :: user_numbers(6)
    integer, intent(in), optional :: date_since
    integer :: i

    do i = 1, size(draws)
      if (present(date_since)) then
        if (draws(i)%date_yyyymmdd <= date_since) then
          draws(i)%winning_scale = 0
          cycle
        end if
      end if
      draws(i)%matched_main = match_count(user_numbers, &
        draws(i)%main_numbers(:6))
      draws(i)%matched_bonus = match_count(user_numbers, &
        draws(i)%bonus_numbers(:1))
      draws(i)%winning_scale = uk_lottery_scale(draws(i)%matched_main, &
        draws(i)%matched_bonus)
    end do
  end subroutine uk_lottery_backtest

  subroutine uk_thunderball_backtest(draws, user_numbers, date_since)
    type(lottery_draw_t), intent(inout) :: draws(:)
    integer, intent(in) :: user_numbers(6)
    integer, intent(in), optional :: date_since
    integer :: i

    do i = 1, size(draws)
      if (present(date_since)) then
        if (draws(i)%date_yyyymmdd <= date_since) then
          draws(i)%winning_scale = 0
          cycle
        end if
      end if
      draws(i)%matched_main = match_count(user_numbers(:5), &
        draws(i)%main_numbers(:5))
      draws(i)%matched_bonus = match_count(user_numbers(6:6), &
        draws(i)%bonus_numbers(:1))
      draws(i)%winning_scale = thunderball_scale(draws(i)%matched_main, &
        draws(i)%matched_bonus)
    end do
  end subroutine uk_thunderball_backtest

  subroutine calculate_euro_lottery_pnl(draws, result)
    type(lottery_draw_t), intent(in) :: draws(:)
    type(lottery_pnl_result_t), intent(out) :: result
    real(dp), parameter :: prizes(0:13) = &
      [0.0_dp, 4.0_dp, 6.0_dp, 7.0_dp, 9.0_dp, 11.0_dp, 14.0_dp, &
      39.0_dp, 57.0_dp, 120.0_dp, 1299.0_dp, 20851.0_dp, &
      200738.0_dp, 10000000.0_dp]

    call calculate_pnl(draws, prizes, 2.5_dp, result)
  end subroutine calculate_euro_lottery_pnl

  subroutine calculate_set_for_life_pnl(draws, result)
    type(lottery_draw_t), intent(in) :: draws(:)
    type(lottery_pnl_result_t), intent(out) :: result
    real(dp), parameter :: prizes(0:8) = &
      [0.0_dp, 5.0_dp, 10.0_dp, 20.0_dp, 30.0_dp, 50.0_dp, &
      250.0_dp, 120000.0_dp, 3600000.0_dp]

    call calculate_pnl(draws, prizes, 1.5_dp, result)
  end subroutine calculate_set_for_life_pnl

  subroutine calculate_uk_lottery_pnl(draws, result)
    type(lottery_draw_t), intent(in) :: draws(:)
    type(lottery_pnl_result_t), intent(out) :: result
    real(dp), parameter :: prizes(0:6) = &
      [0.0_dp, 0.0_dp, 30.0_dp, 140.0_dp, 1750.0_dp, &
      1000000.0_dp, 100000000.0_dp]

    call calculate_pnl(draws, prizes, 2.0_dp, result)
  end subroutine calculate_uk_lottery_pnl

  subroutine calculate_uk_thunderball_pnl(draws, result)
    type(lottery_draw_t), intent(in) :: draws(:)
    type(lottery_pnl_result_t), intent(out) :: result
    real(dp), parameter :: prizes(0:9) = &
      [0.0_dp, 3.0_dp, 5.0_dp, 10.0_dp, 10.0_dp, 20.0_dp, &
      100.0_dp, 250.0_dp, 5000.0_dp, 500000.0_dp]

    call calculate_pnl(draws, prizes, 1.0_dp, result)
  end subroutine calculate_uk_thunderball_pnl

  subroutine calculate_pnl(draws, prizes, ticket_cost, result)
    type(lottery_draw_t), intent(in) :: draws(:)
    real(dp), intent(in) :: prizes(0:)
    real(dp), intent(in) :: ticket_cost
    type(lottery_pnl_result_t), intent(out) :: result
    integer, allocatable :: order(:)
    integer :: i
    integer :: n
    integer :: scale

    n = size(draws)
    allocate(order(n))
    order = [(i, i = 1, n)]
    call sort_draw_order(draws, order)

    allocate(result%date_yyyymmdd(n), result%payout(n))
    allocate(result%cumulative_profit(n), result%cumulative_cost(n))
    allocate(result%final_pnl(n))

    do i = 1, n
      result%date_yyyymmdd(i) = draws(order(i))%date_yyyymmdd
      scale = draws(order(i))%winning_scale
      if (scale >= lbound(prizes, 1) .and. scale <= ubound(prizes, 1)) then
        result%payout(i) = prizes(scale)
      else
        result%payout(i) = 0.0_dp
      end if
      if (i == 1) then
        result%cumulative_profit(i) = result%payout(i)
      else
        result%cumulative_profit(i) = result%cumulative_profit(i - 1) + &
          result%payout(i)
      end if
      result%cumulative_cost(i) = real(i, dp) * ticket_cost
      result%final_pnl(i) = result%cumulative_profit(i) - &
        result%cumulative_cost(i)
    end do
  end subroutine calculate_pnl

  subroutine top_five_numbers(draws, numbers, least_lucky, max_number, date_since)
    type(lottery_draw_t), intent(in) :: draws(:)
    integer, intent(out) :: numbers(5)
    logical, intent(in), optional :: least_lucky
    integer, intent(in), optional :: max_number
    integer, intent(in), optional :: date_since
    integer, allocatable :: counts(:)
    integer, allocatable :: order(:)
    integer :: i
    integer :: j
    integer :: maximum
    logical :: choose_least

    maximum = 50
    if (present(max_number)) maximum = max_number
    choose_least = .false.
    if (present(least_lucky)) choose_least = least_lucky

    allocate(counts(maximum), order(maximum))
    counts = 0
    order = [(i, i = 1, maximum)]

    do i = 1, size(draws)
      if (present(date_since)) then
        if (draws(i)%date_yyyymmdd <= date_since) cycle
      end if
      do j = 1, draws(i)%n_main
        if (draws(i)%main_numbers(j) >= 1 .and. &
            draws(i)%main_numbers(j) <= maximum) then
          counts(draws(i)%main_numbers(j)) = &
            counts(draws(i)%main_numbers(j)) + 1
        end if
      end do
    end do

    call sort_number_order(counts, order)
    if (choose_least) then
      numbers = order(:5)
    else
      numbers = order(maximum - 4:maximum)
    end if
  end subroutine top_five_numbers

  pure integer(kind=8) function euro_lottery_combination_count() result(value)
    value = 139838160_8
  end function euro_lottery_combination_count

  subroutine reset_euro_combination(self)
    class(euro_combination_iterator_t), intent(inout) :: self

    self%main_numbers = [1, 2, 3, 4, 5]
    self%bonus_numbers = [1, 2]
    self%started = .false.
    self%finished = .false.
  end subroutine reset_euro_combination

  subroutine next_euro_combination(self, numbers, has_value)
    class(euro_combination_iterator_t), intent(inout) :: self
    integer, intent(out) :: numbers(7)
    logical, intent(out) :: has_value

    if (self%finished) then
      numbers = 0
      has_value = .false.
      return
    end if

    if (.not. self%started) then
      self%started = .true.
    else
      if (.not. advance_combination(self%bonus_numbers, 12)) then
        self%bonus_numbers = [1, 2]
        if (.not. advance_combination(self%main_numbers, 50)) then
          self%finished = .true.
          numbers = 0
          has_value = .false.
          return
        end if
      end if
    end if

    numbers(:5) = self%main_numbers
    numbers(6:7) = self%bonus_numbers
    has_value = .true.
  end subroutine next_euro_combination

  subroutine outer_join_merge_integer(a, b, result)
    integer, intent(in) :: a(:, :)
    integer, intent(in) :: b(:, :)
    integer, allocatable, intent(out) :: result(:, :)
    integer :: i
    integer :: j
    integer :: row

    allocate(result(size(a, 1) * size(b, 1), size(a, 2) + size(b, 2)))
    row = 0
    do j = 1, size(b, 1)
      do i = 1, size(a, 1)
        row = row + 1
        result(row, :size(a, 2)) = a(i, :)
        result(row, size(a, 2) + 1:) = b(j, :)
      end do
    end do
  end subroutine outer_join_merge_integer

  pure integer function match_count(selected, drawn) result(count_value)
    integer, intent(in) :: selected(:)
    integer, intent(in) :: drawn(:)
    integer :: i

    count_value = 0
    do i = 1, size(selected)
      if (any(selected(i) == drawn)) count_value = count_value + 1
    end do
  end function match_count

  pure integer function euro_scale(main_matches, bonus_matches) result(scale)
    integer, intent(in) :: main_matches
    integer, intent(in) :: bonus_matches

    scale = 0
    select case (main_matches)
    case (1)
      if (bonus_matches == 2) scale = 3
    case (2)
      if (bonus_matches == 0) scale = 1
      if (bonus_matches == 1) scale = 2
      if (bonus_matches == 2) scale = 6
    case (3)
      if (bonus_matches == 0) scale = 4
      if (bonus_matches == 1) scale = 5
      if (bonus_matches == 2) scale = 8
    case (4)
      if (bonus_matches == 0) scale = 7
      if (bonus_matches == 1) scale = 9
      if (bonus_matches == 2) scale = 10
    case (5)
      if (bonus_matches == 0) scale = 11
      if (bonus_matches == 1) scale = 12
      if (bonus_matches == 2) scale = 13
    end select
  end function euro_scale

  pure integer function set_for_life_scale(main_matches, bonus_matches) result(scale)
    integer, intent(in) :: main_matches
    integer, intent(in) :: bonus_matches

    scale = 0
    if (main_matches >= 2 .and. main_matches <= 5) then
      scale = 2 * (main_matches - 2) + 1 + bonus_matches
    end if
  end function set_for_life_scale

  pure integer function uk_lottery_scale(main_matches, bonus_matches) result(scale)
    integer, intent(in) :: main_matches
    integer, intent(in) :: bonus_matches

    scale = 0
    select case (main_matches)
    case (2)
      scale = 1
    case (3)
      scale = 2
    case (4)
      scale = 3
    case (5)
      scale = 4 + min(bonus_matches, 1)
    case (6)
      scale = 6
    end select
  end function uk_lottery_scale

  pure integer function thunderball_scale(main_matches, bonus_matches) result(scale)
    integer, intent(in) :: main_matches
    integer, intent(in) :: bonus_matches

    scale = 0
    if (bonus_matches == 1 .and. main_matches <= 2) then
      scale = main_matches + 1
    else if (main_matches == 3) then
      scale = 4 + bonus_matches
    else if (main_matches == 4) then
      scale = 6 + bonus_matches
    else if (main_matches == 5) then
      scale = 8 + bonus_matches
    end if
  end function thunderball_scale

  logical function advance_combination(values, maximum) result(advanced)
    integer, intent(inout) :: values(:)
    integer, intent(in) :: maximum
    integer :: i
    integer :: j
    integer :: k

    k = size(values)
    do i = k, 1, -1
      if (values(i) < maximum - (k - i)) then
        values(i) = values(i) + 1
        do j = i + 1, k
          values(j) = values(j - 1) + 1
        end do
        advanced = .true.
        return
      end if
    end do
    advanced = .false.
  end function advance_combination

  subroutine sort_draw_order(draws, order)
    type(lottery_draw_t), intent(in) :: draws(:)
    integer, intent(inout) :: order(:)
    integer :: i
    integer :: j
    integer :: key

    do i = 2, size(order)
      key = order(i)
      j = i - 1
      do while (j >= 1)
        if (draws(order(j))%date_yyyymmdd <= draws(key)%date_yyyymmdd) exit
        order(j + 1) = order(j)
        j = j - 1
      end do
      order(j + 1) = key
    end do
  end subroutine sort_draw_order

  subroutine sort_number_order(counts, order)
    integer, intent(in) :: counts(:)
    integer, intent(inout) :: order(:)
    integer :: i
    integer :: j
    integer :: key

    do i = 2, size(order)
      key = order(i)
      j = i - 1
      do while (j >= 1)
        if (counts(order(j)) < counts(key)) exit
        if (counts(order(j)) == counts(key) .and. order(j) < key) exit
        order(j + 1) = order(j)
        j = j - 1
      end do
      order(j + 1) = key
    end do
  end subroutine sort_number_order

  integer function month_number(text) result(month)
    character(len=*), intent(in) :: text
    character(len=:), allocatable :: name

    name = lowercase(trim(text))
    select case (name)
    case ("january")
      month = 1
    case ("february")
      month = 2
    case ("march")
      month = 3
    case ("april")
      month = 4
    case ("may")
      month = 5
    case ("june")
      month = 6
    case ("july")
      month = 7
    case ("august", "augu")
      month = 8
    case ("september")
      month = 9
    case ("october")
      month = 10
    case ("november")
      month = 11
    case ("december")
      month = 12
    case default
      month = 0
    end select
  end function month_number

  integer function parse_ordinal_day(text) result(day)
    character(len=*), intent(in) :: text
    character(len=str_len) :: work
    integer :: i

    work = trim(text)
    do i = 1, len_trim(work)
      if (work(i:i) < '0' .or. work(i:i) > '9') then
        work(i:) = ""
        exit
      end if
    end do
    day = parse_integer_or_zero(work)
  end function parse_ordinal_day

  subroutine replace_character(text, old_character, new_character)
    character(len=*), intent(inout) :: text
    character(len=1), intent(in) :: old_character
    character(len=1), intent(in) :: new_character
    integer :: i

    do i = 1, len(text)
      if (text(i:i) == old_character) text(i:i) = new_character
    end do
  end subroutine replace_character

  subroutine remove_empty_fields(fields, count)
    character(len=str_len), intent(inout) :: fields(:)
    integer, intent(inout) :: count
    integer :: i
    integer :: output_index

    output_index = 0
    do i = 1, count
      if (len_trim(fields(i)) > 0) then
        output_index = output_index + 1
        fields(output_index) = fields(i)
      end if
    end do
    if (output_index < size(fields)) fields(output_index + 1:) = ""
    count = output_index
  end subroutine remove_empty_fields

end module trading_lottery
