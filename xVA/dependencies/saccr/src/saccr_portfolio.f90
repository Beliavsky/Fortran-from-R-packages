module saccr_portfolio
  use trading, only : dp, str_len, trade_t, csa_t, collateral_t
  use saccr_types, only : supervisory_record_t, exposure_result_t, portfolio_result_t, &
    addon_result_t, replacement_cost_t
  use saccr_addon, only : calc_addon
  use saccr_core, only : calc_rc, calc_pfe, calc_ead
  implicit none
  private

  public :: calculate_exposure
  public :: calculate_portfolio

contains

  subroutine calculate_exposure(trades, result, csa, collaterals, simplified, &
      oem, ignore_margin, records)
    type(trade_t), intent(in) :: trades(:)
    type(exposure_result_t), intent(out) :: result
    type(csa_t), intent(in), optional :: csa
    type(collateral_t), intent(in), optional :: collaterals(:)
    logical, intent(in), optional :: simplified
    logical, intent(in), optional :: oem
    logical, intent(in), optional :: ignore_margin
    type(supervisory_record_t), intent(in), optional :: records(:)
    type(addon_result_t) :: unmargined_addon
    type(replacement_cost_t) :: unmargined_rc
    real(dp) :: maturity_factor
    real(dp) :: unmargined_pfe
    logical :: margin_ignored
    logical :: use_oem
    logical :: use_simplified

    result = exposure_result_t()
    if (size(trades) == 0) return
    result%counterparty = trades(1)%counterparty
    if (len_trim(result%counterparty) == 0) result%counterparty = "Counterparty"

    use_simplified = .false.
    if (present(simplified)) use_simplified = simplified
    use_oem = .false.
    if (present(oem)) use_oem = oem
    if (use_oem) use_simplified = .true.
    margin_ignored = .false.
    if (present(ignore_margin)) margin_ignored = ignore_margin

    if (present(csa)) then
      result%margined = .not. margin_ignored
      result%csa_id = csa%id
      if (margin_ignored) then
        maturity_factor = 1.0_dp
      else
        maturity_factor = csa%maturity_factor(use_simplified)
      end if
      result%maturity_factor = maturity_factor
      if (present(records)) then
        call calc_addon(trades, result%addon, maturity_factor=maturity_factor, &
          simplified=use_simplified, oem=use_oem, records=records)
      else
        call calc_addon(trades, result%addon, maturity_factor=maturity_factor, &
          simplified=use_simplified, oem=use_oem)
      end if
      if (present(collaterals)) then
        call calc_rc(trades, result%replacement_cost, csa=csa, &
          collaterals=collaterals, simplified=use_simplified, &
          ignore_margin=margin_ignored)
      else
        call calc_rc(trades, result%replacement_cost, csa=csa, &
          simplified=use_simplified, ignore_margin=margin_ignored)
      end if
    else
      if (present(records)) then
        call calc_addon(trades, result%addon, simplified=use_simplified, &
          oem=use_oem, records=records)
      else
        call calc_addon(trades, result%addon, simplified=use_simplified, oem=use_oem)
      end if
      call calc_rc(trades, result%replacement_cost, simplified=use_simplified, &
        ignore_margin=margin_ignored)
    end if

    result%pfe = calc_pfe(result%replacement_cost%v_c, result%addon%addon, &
      v=result%replacement_cost%v, simplified=use_simplified)
    result%ead = calc_ead(result%replacement_cost%rc, result%pfe)

    if (present(records)) then
      call calc_addon(trades, unmargined_addon, simplified=use_simplified, &
        oem=use_oem, records=records)
    else
      call calc_addon(trades, unmargined_addon, simplified=use_simplified, oem=use_oem)
    end if
    call calc_rc(trades, unmargined_rc, simplified=use_simplified)
    unmargined_pfe = calc_pfe(unmargined_rc%v_c, unmargined_addon%addon, &
      v=unmargined_rc%v, simplified=use_simplified)
    result%unmargined_ead = calc_ead(unmargined_rc%rc, unmargined_pfe)
    result%ead = min(result%ead, result%unmargined_ead)
  end subroutine calculate_exposure

  subroutine calculate_portfolio(trades, result, csas, collaterals, simplified, &
      oem, ignore_margin, records)
    type(trade_t), intent(in) :: trades(:)
    type(portfolio_result_t), intent(out) :: result
    type(csa_t), intent(in), optional :: csas(:)
    type(collateral_t), intent(in), optional :: collaterals(:)
    logical, intent(in), optional :: simplified
    logical, intent(in), optional :: oem
    logical, intent(in), optional :: ignore_margin
    type(supervisory_record_t), intent(in), optional :: records(:)
    character(len=str_len), allocatable :: counterparties(:)
    type(trade_t), allocatable :: subset(:)
    type(exposure_result_t) :: exposure
    logical, allocatable :: assigned(:)
    logical :: margin_ignored
    logical :: use_oem
    logical :: use_simplified
    integer :: c
    integer :: i

    result = portfolio_result_t()
    allocate(result%exposures(0))
    if (size(trades) == 0) return
    if (all_sold_options(trades)) return

    use_simplified = .false.
    if (present(simplified)) use_simplified = simplified
    use_oem = .false.
    if (present(oem)) use_oem = oem
    if (use_oem) use_simplified = .true.
    margin_ignored = .false.
    if (present(ignore_margin)) margin_ignored = ignore_margin

    call unique_counterparties(trades, counterparties)
    allocate(assigned(size(trades)))

    do c = 1, size(counterparties)
      assigned = .false.
      do i = 1, size(trades)
        if (.not. same_text(trades(i)%counterparty, counterparties(c))) assigned(i) = .true.
      end do

      if (present(csas)) then
        do i = 1, size(csas)
          if (.not. same_text(csas(i)%counterparty, counterparties(c))) cycle
          call collect_csa_trades(trades, assigned, csas(i), subset)
          if (size(subset) == 0) cycle
          if (present(collaterals)) then
            if (present(records)) then
              call calculate_exposure(subset, exposure, csa=csas(i), &
                collaterals=collaterals, simplified=use_simplified, oem=use_oem, &
                ignore_margin=margin_ignored, records=records)
            else
              call calculate_exposure(subset, exposure, csa=csas(i), &
                collaterals=collaterals, simplified=use_simplified, oem=use_oem, &
                ignore_margin=margin_ignored)
            end if
          else
            if (present(records)) then
              call calculate_exposure(subset, exposure, csa=csas(i), &
                simplified=use_simplified, oem=use_oem, &
                ignore_margin=margin_ignored, records=records)
            else
              call calculate_exposure(subset, exposure, csa=csas(i), &
                simplified=use_simplified, oem=use_oem, &
                ignore_margin=margin_ignored)
            end if
          end if
          call append_exposure(result%exposures, exposure)
        end do
      end if

      call collect_unassigned_trades(trades, assigned, subset)
      if (size(subset) > 0) then
        if (present(records)) then
          call calculate_exposure(subset, exposure, simplified=use_simplified, &
            oem=use_oem, ignore_margin=margin_ignored, records=records)
        else
          call calculate_exposure(subset, exposure, simplified=use_simplified, &
            oem=use_oem, ignore_margin=margin_ignored)
        end if
        call append_exposure(result%exposures, exposure)
      end if
    end do

    result%total_ead = 0.0_dp
    do i = 1, size(result%exposures)
      result%total_ead = result%total_ead + result%exposures(i)%ead
    end do
  end subroutine calculate_portfolio

  subroutine collect_csa_trades(trades, assigned, csa, subset)
    type(trade_t), intent(in) :: trades(:)
    logical, intent(inout) :: assigned(:)
    type(csa_t), intent(in) :: csa
    type(trade_t), allocatable, intent(out) :: subset(:)
    logical, allocatable :: selected(:)
    integer :: count
    integer :: i
    integer :: output_index

    allocate(selected(size(trades)))
    selected = .false.
    count = 0
    do i = 1, size(trades)
      if (assigned(i)) cycle
      if (.not. same_text(trades(i)%counterparty, csa%counterparty)) cycle
      if (.not. list_contains(csa%trade_groups, trades(i)%trade_group)) cycle
      if (.not. list_contains(csa%currencies, trades(i)%currency)) cycle
      selected(i) = .true.
      count = count + 1
    end do

    allocate(subset(count))
    output_index = 0
    do i = 1, size(trades)
      if (.not. selected(i)) cycle
      output_index = output_index + 1
      subset(output_index) = trades(i)
      assigned(i) = .true.
    end do
  end subroutine collect_csa_trades

  subroutine collect_unassigned_trades(trades, assigned, subset)
    type(trade_t), intent(in) :: trades(:)
    logical, intent(inout) :: assigned(:)
    type(trade_t), allocatable, intent(out) :: subset(:)
    integer :: n_selected
    integer :: i
    integer :: output_index

    n_selected = count(.not. assigned)
    allocate(subset(n_selected))
    output_index = 0
    do i = 1, size(trades)
      if (assigned(i)) cycle
      output_index = output_index + 1
      subset(output_index) = trades(i)
      assigned(i) = .true.
    end do
  end subroutine collect_unassigned_trades

  subroutine unique_counterparties(trades, values)
    type(trade_t), intent(in) :: trades(:)
    character(len=str_len), allocatable, intent(out) :: values(:)
    character(len=str_len), allocatable :: work(:)
    character(len=str_len) :: value
    integer :: i
    integer :: j
    integer :: n
    logical :: found

    allocate(values(0))
    do i = 1, size(trades)
      value = trades(i)%counterparty
      found = .false.
      do j = 1, size(values)
        if (same_text(values(j), value)) then
          found = .true.
          exit
        end if
      end do
      if (found) cycle
      n = size(values)
      allocate(work(n + 1))
      if (n > 0) work(:n) = values
      work(n + 1) = value
      call move_alloc(work, values)
    end do
  end subroutine unique_counterparties

  subroutine append_exposure(values, value)
    type(exposure_result_t), allocatable, intent(inout) :: values(:)
    type(exposure_result_t), intent(in) :: value
    type(exposure_result_t), allocatable :: work(:)
    integer :: n

    n = size(values)
    allocate(work(n + 1))
    if (n > 0) work(:n) = values
    work(n + 1) = value
    call move_alloc(work, values)
  end subroutine append_exposure

  logical function list_contains(values, value) result(found)
    character(len=str_len), allocatable, intent(in) :: values(:)
    character(len=*), intent(in) :: value
    integer :: i

    if (.not. allocated(values)) then
      found = .true.
      return
    end if
    found = .false.
    do i = 1, size(values)
      if (same_text(values(i), value)) then
        found = .true.
        return
      end if
    end do
  end function list_contains

  logical function all_sold_options(trades) result(value)
    type(trade_t), intent(in) :: trades(:)
    integer :: i

    value = size(trades) > 0
    do i = 1, size(trades)
      if (.not. same_text(trades(i)%trade_type, "Option") .or. &
          .not. same_text(trades(i)%buy_sell, "Sell")) then
        value = .false.
        return
      end if
    end do
  end function all_sold_options

  pure logical function same_text(left, right) result(equal)
    character(len=*), intent(in) :: left
    character(len=*), intent(in) :: right

    equal = uppercase(trim(adjustl(left))) == uppercase(trim(adjustl(right)))
  end function same_text

  pure function uppercase(text) result(value)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: value
    integer :: code
    integer :: i

    value = text
    do i = 1, len(text)
      code = iachar(value(i:i))
      if (code >= iachar('a') .and. code <= iachar('z')) then
        value(i:i) = achar(code - iachar('a') + iachar('A'))
      end if
    end do
  end function uppercase

end module saccr_portfolio
