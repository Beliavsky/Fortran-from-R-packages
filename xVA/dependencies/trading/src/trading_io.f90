module trading_io
  use trading_kinds, only : dp, str_len
  use trading_strings, only : parse_real_or_zero, split_delimited
  use trading_trades, only : trade_t
  use trading_curve, only : curve_t
  use trading_csa, only : csa_t, collateral_t
  implicit none
  private

  public :: parse_trades_csv
  public :: load_curve_csv
  public :: load_csa_csv
  public :: load_collateral_csv
  public :: load_track_record_csv
  public :: write_trade_details

contains

  subroutine parse_trades_csv(filename, trades)
    character(len=*), intent(in) :: filename
    type(trade_t), allocatable, intent(out) :: trades(:)
    character(len=4096) :: line
    character(len=str_len) :: headers(96)
    character(len=str_len) :: fields(96)
    type(trade_t), allocatable :: work(:)
    type(trade_t) :: current
    integer :: count
    integer :: field_count
    integer :: header_count
    integer :: ios
    integer :: output_count
    integer :: unit

    open(newunit=unit, file=filename, status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "parse_trades_csv: unable to open file"
    read(unit, '(a)', iostat=ios) line
    if (ios /= 0) error stop "parse_trades_csv: missing header"
    call split_delimited(line, ',', headers, header_count)

    count = 0
    do
      read(unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) > 0) count = count + 1
    end do
    rewind(unit)
    read(unit, '(a)', iostat=ios) line

    allocate(work(2 * max(count, 1)))
    output_count = 0
    do
      read(unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) cycle
      call split_delimited(line, ',', fields, field_count)
      call populate_trade(headers, header_count, fields, field_count, current)

      if (trim(current%exotic_type) == "Digital") then
        output_count = output_count + 1
        work(output_count) = current
        work(output_count)%buy_sell = "Buy"
        work(output_count)%strike_price = 0.95_dp * current%strike_price
        work(output_count)%external_id = trim(current%external_id) // "L"
        write(work(output_count)%exotic_type, '(a,i0)') "Digital", output_count

        output_count = output_count + 1
        work(output_count) = current
        work(output_count)%buy_sell = "Sell"
        work(output_count)%strike_price = 1.05_dp * current%strike_price
        work(output_count)%external_id = trim(current%external_id) // "S"
        write(work(output_count)%exotic_type, '(a,i0)') "Digital", output_count
      else
        output_count = output_count + 1
        work(output_count) = current
      end if
    end do
    close(unit)

    allocate(trades(output_count))
    trades = work(:output_count)
  end subroutine parse_trades_csv

  subroutine load_curve_csv(filename, curve)
    character(len=*), intent(in) :: filename
    type(curve_t), intent(out) :: curve
    character(len=1024) :: line
    character(len=str_len) :: fields(8)
    real(dp), allocatable :: tenors(:)
    real(dp), allocatable :: rates(:)
    integer :: count
    integer :: field_count
    integer :: ios
    integer :: row
    integer :: unit

    open(newunit=unit, file=filename, status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "load_curve_csv: unable to open file"
    read(unit, '(a)', iostat=ios) line
    count = 0
    do
      read(unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) > 0) count = count + 1
    end do
    rewind(unit)
    read(unit, '(a)', iostat=ios) line

    allocate(tenors(count), rates(count))
    row = 0
    do
      read(unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) cycle
      call split_delimited(line, ',', fields, field_count)
      if (field_count < 2) cycle
      row = row + 1
      tenors(row) = parse_real_or_zero(fields(1))
      rates(row) = parse_real_or_zero(fields(2))
    end do
    close(unit)

    allocate(curve%tenors(row), curve%rates(row))
    curve%tenors = tenors(:row)
    curve%rates = rates(:row)
  end subroutine load_curve_csv

  subroutine load_csa_csv(filename, agreements)
    character(len=*), intent(in) :: filename
    type(csa_t), allocatable, intent(out) :: agreements(:)
    character(len=2048) :: line
    character(len=str_len) :: headers(32)
    character(len=str_len) :: fields(32)
    character(len=str_len) :: list_fields(32)
    integer :: count
    integer :: field_count
    integer :: header_count
    integer :: ios
    integer :: list_count
    integer :: row
    integer :: unit

    open(newunit=unit, file=filename, status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "load_csa_csv: unable to open file"
    read(unit, '(a)', iostat=ios) line
    call split_delimited(line, ',', headers, header_count)
    count = 0
    do
      read(unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) > 0) count = count + 1
    end do
    rewind(unit)
    read(unit, '(a)', iostat=ios) line

    allocate(agreements(count))
    row = 0
    do
      read(unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) cycle
      call split_delimited(line, ',', fields, field_count)
      row = row + 1
      agreements(row)%id = field_value(headers, fields, header_count, "ID")
      agreements(row)%threshold_counterparty = real_field(headers, fields, header_count, "thres_cpty")
      agreements(row)%threshold_processing_organization = &
        real_field(headers, fields, header_count, "thres_PO")
      agreements(row)%initial_margin_counterparty = &
        real_field(headers, fields, header_count, "IM_cpty")
      agreements(row)%initial_margin_processing_organization = &
        real_field(headers, fields, header_count, "IM_PO")
      agreements(row)%minimum_transfer_counterparty = &
        real_field(headers, fields, header_count, "MTA_cpty")
      agreements(row)%minimum_transfer_processing_organization = &
        real_field(headers, fields, header_count, "MTA_PO")
      agreements(row)%mpor_days = real_field(headers, fields, header_count, "mpor_days")
      agreements(row)%remargin_frequency = &
        real_field(headers, fields, header_count, "remargin_freq")
      agreements(row)%rounding = real_field(headers, fields, header_count, "rounding")
      agreements(row)%counterparty = &
        field_value(headers, fields, header_count, "Counterparty")
      agreements(row)%values_type = &
        field_value(headers, fields, header_count, "Values_type")

      call split_delimited(field_value(headers, fields, header_count, "Currency"), &
        ';', list_fields, list_count)
      allocate(agreements(row)%currencies(list_count))
      agreements(row)%currencies = list_fields(:list_count)

      call split_delimited(field_value(headers, fields, header_count, "TradeGroups"), &
        ';', list_fields, list_count)
      allocate(agreements(row)%trade_groups(list_count))
      agreements(row)%trade_groups = list_fields(:list_count)
    end do
    close(unit)
  end subroutine load_csa_csv

  subroutine load_collateral_csv(filename, collateral)
    character(len=*), intent(in) :: filename
    type(collateral_t), allocatable, intent(out) :: collateral(:)
    character(len=1024) :: line
    character(len=str_len) :: headers(16)
    character(len=str_len) :: fields(16)
    integer :: count
    integer :: field_count
    integer :: header_count
    integer :: ios
    integer :: row
    integer :: unit

    open(newunit=unit, file=filename, status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "load_collateral_csv: unable to open file"
    read(unit, '(a)', iostat=ios) line
    call split_delimited(line, ',', headers, header_count)
    count = 0
    do
      read(unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) > 0) count = count + 1
    end do
    rewind(unit)
    read(unit, '(a)', iostat=ios) line

    allocate(collateral(count))
    row = 0
    do
      read(unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) cycle
      call split_delimited(line, ',', fields, field_count)
      row = row + 1
      collateral(row)%id = field_value(headers, fields, header_count, "ID")
      collateral(row)%amount = real_field(headers, fields, header_count, "Amount")
      collateral(row)%csa_id = field_value(headers, fields, header_count, "csa_id")
      collateral(row)%collateral_type = field_value(headers, fields, header_count, "type")
    end do
    close(unit)
  end subroutine load_collateral_csv

  subroutine load_track_record_csv(filename, fund, benchmark)
    character(len=*), intent(in) :: filename
    real(dp), allocatable, intent(out) :: fund(:)
    real(dp), allocatable, intent(out) :: benchmark(:)
    character(len=1024) :: line
    character(len=str_len) :: fields(8)
    real(dp), allocatable :: temp_fund(:)
    real(dp), allocatable :: temp_benchmark(:)
    integer :: count
    integer :: field_count
    integer :: ios
    integer :: row
    integer :: unit

    open(newunit=unit, file=filename, status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "load_track_record_csv: unable to open file"
    read(unit, '(a)', iostat=ios) line
    count = 0
    do
      read(unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) > 0) count = count + 1
    end do
    rewind(unit)
    read(unit, '(a)', iostat=ios) line

    allocate(temp_fund(count), temp_benchmark(count))
    row = 0
    do
      read(unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) cycle
      call split_delimited(line, ',', fields, field_count)
      if (field_count < 3) cycle
      row = row + 1
      temp_fund(row) = parse_real_or_zero(fields(2))
      temp_benchmark(row) = parse_real_or_zero(fields(3))
    end do
    close(unit)

    allocate(fund(row), benchmark(row))
    fund = temp_fund(:row)
    benchmark = temp_benchmark(:row)
  end subroutine load_track_record_csv

  subroutine write_trade_details(trade, unit)
    type(trade_t), intent(in) :: trade
    integer, intent(in), optional :: unit
    integer :: output_unit

    output_unit = 6
    if (present(unit)) output_unit = unit
    write(output_unit, '(a)') "class_name=" // trim(trade%class_name)
    if (len_trim(trade%trade_group) > 0) then
      write(output_unit, '(a)') "trade_group=" // trim(trade%trade_group)
    end if
    if (len_trim(trade%trade_type) > 0) then
      write(output_unit, '(a)') "trade_type=" // trim(trade%trade_type)
    end if
    if (abs(trade%notional) > 0.0_dp) write(output_unit, '(a,es24.16)') "notional=", trade%notional
    if (abs(trade%mtm) > 0.0_dp) write(output_unit, '(a,es24.16)') "mtm=", trade%mtm
    if (len_trim(trade%currency) > 0) then
      write(output_unit, '(a)') "currency=" // trim(trade%currency)
    end if
    if (len_trim(trade%buy_sell) > 0) then
      write(output_unit, '(a)') "buy_sell=" // trim(trade%buy_sell)
    end if
    if (len_trim(trade%external_id) > 0) then
      write(output_unit, '(a)') "external_id=" // trim(trade%external_id)
    end if
  end subroutine write_trade_details

  subroutine populate_trade(headers, header_count, fields, field_count, trade)
    character(len=str_len), intent(in) :: headers(:)
    integer, intent(in) :: header_count
    character(len=str_len), intent(in) :: fields(:)
    integer, intent(in) :: field_count
    type(trade_t), intent(out) :: trade
    character(len=str_len) :: class_name

    class_name = field_value(headers, fields, min(header_count, field_count), "TradeObj")
    trade = trade_t()
    call trade%configure_class(trim(class_name))
    trade%subclass = field_value(headers, fields, header_count, "SubClass")
    trade%notional = real_field(headers, fields, header_count, "Notional")
    trade%mtm = real_field(headers, fields, header_count, "MtM")
    trade%currency = field_value(headers, fields, header_count, "Currency")
    trade%si = real_field(headers, fields, header_count, "Si")
    trade%ei = real_field(headers, fields, header_count, "Ei")
    trade%buy_sell = field_value(headers, fields, header_count, "BuySell")
    trade%external_id = field_value(headers, fields, header_count, "external_id")
    trade%counterparty = field_value(headers, fields, header_count, "Counterparty")
    trade%option_type = field_value(headers, fields, header_count, "OptionType")
    trade%underlying_price = real_field(headers, fields, header_count, "UnderlyingPrice")
    trade%strike_price = real_field(headers, fields, header_count, "StrikePrice")
    trade%commodity_type = field_value(headers, fields, header_count, "commodity_type")
    trade%reference_entity = field_value(headers, fields, header_count, "RefEntity")
    trade%ccy_pair = field_value(headers, fields, header_count, "ccyPair")
    trade%pay_leg_type = field_value(headers, fields, header_count, "pay_leg_type")
    trade%pay_leg_ref = field_value(headers, fields, header_count, "pay_leg_ref")
    trade%pay_leg_tenor = field_value(headers, fields, header_count, "pay_leg_tenor")
    trade%rec_leg_type = field_value(headers, fields, header_count, "rec_leg_type")
    trade%rec_leg_ref = field_value(headers, fields, header_count, "rec_leg_ref")
    trade%rec_leg_tenor = field_value(headers, fields, header_count, "rec_leg_tenor")
    trade%vol_strike = real_field(headers, fields, header_count, "vol_strike")
    trade%trade_price = real_field(headers, fields, header_count, "traded_price")
    trade%delivery_type = field_value(headers, fields, header_count, "del_type")
    trade%underlying_instrument = &
      field_value(headers, fields, header_count, "Underlying_Instrument")
    trade%exotic_type = field_value(headers, fields, header_count, "Exotic_Type")
    trade%isin = field_value(headers, fields, header_count, "ISIN")
    trade%coupon_type = field_value(headers, fields, header_count, "coupon_type")
    trade%issuer = field_value(headers, fields, header_count, "Issuer")
    trade%maturity_date = field_value(headers, fields, header_count, "maturity_date")
    trade%payment_frequency = &
      field_value(headers, fields, header_count, "payment_frequency")
    trade%credit_risk_weight = &
      real_field(headers, fields, header_count, "credit_risk_weight")
    trade%cdo_attach_point = &
      real_field(headers, fields, header_count, "cdo_attach_point")
    trade%cdo_detach_point = &
      real_field(headers, fields, header_count, "cdo_detach_point")
    trade%fx_near_leg_fields = &
      field_value(headers, fields, header_count, "fx_near_leg_fields")
    trade%ccy_paying = field_value(headers, fields, header_count, "ccy_paying")
    trade%amount_paying = real_field(headers, fields, header_count, "amount_paying")
    trade%ccy_receiving = field_value(headers, fields, header_count, "ccy_receiving")
    trade%amount_receiving = &
      real_field(headers, fields, header_count, "amount_receiving")
  end subroutine populate_trade

  function field_value(headers, fields, count, name) result(value)
    character(len=str_len), intent(in) :: headers(:)
    character(len=str_len), intent(in) :: fields(:)
    integer, intent(in) :: count
    character(len=*), intent(in) :: name
    character(len=str_len) :: value
    integer :: i

    value = ""
    do i = 1, min(count, min(size(headers), size(fields)))
      if (trim(adjustl(headers(i))) == trim(name)) then
        value = adjustl(trim(fields(i)))
        return
      end if
    end do
  end function field_value

  real(dp) function real_field(headers, fields, count, name) result(value)
    character(len=str_len), intent(in) :: headers(:)
    character(len=str_len), intent(in) :: fields(:)
    integer, intent(in) :: count
    character(len=*), intent(in) :: name

    value = parse_real_or_zero(field_value(headers, fields, count, name))
  end function real_field

end module trading_io
