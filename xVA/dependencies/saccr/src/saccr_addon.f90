module saccr_addon
  use trading, only : dp, str_len, trade_t
  use saccr_types, only : supervisory_record_t, single_trade_addon_t, &
    hedging_set_result_t, asset_class_result_t, addon_result_t
  use saccr_supervisory, only : default_supervisory_data, &
    supervisory_factor, supervisory_correlation
  use saccr_core, only : calculate_factor_multiplier, single_trade_addon
  implicit none
  private

  public :: calc_addon
  public :: hedging_set_name

contains

  subroutine calc_addon(trades, result, maturity_factor, simplified, oem, records)
    type(trade_t), intent(in) :: trades(:)
    type(addon_result_t), intent(out) :: result
    real(dp), intent(in), optional :: maturity_factor
    logical, intent(in), optional :: simplified
    logical, intent(in), optional :: oem
    type(supervisory_record_t), intent(in), optional :: records(:)
    type(supervisory_record_t), allocatable :: supervisory(:)
    type(single_trade_addon_t), allocatable :: details(:)
    type(trade_t), allocatable :: work(:)
    character(len=str_len), allocatable :: asset_names(:)
    logical :: use_oem
    logical :: use_simplified
    integer :: i

    result = addon_result_t()
    allocate(result%trades(0), result%hedging_sets(0), result%asset_classes(0))
    if (size(trades) == 0) return

    use_simplified = .false.
    if (present(simplified)) use_simplified = simplified
    use_oem = .false.
    if (present(oem)) use_oem = oem
    if (use_oem) use_simplified = .true.

    if (present(records)) then
      allocate(supervisory(size(records)))
      supervisory = records
    else
      call default_supervisory_data(supervisory)
    end if

    allocate(work(size(trades)), details(size(trades)))
    work = trades
    do i = 1, size(work)
      if (same_text(work(i)%trade_group, "FX")) then
        if (len_trim(work(i)%base_ccy) == 0) work(i)%base_ccy = "EUR"
        call work(i)%set_fx_dynamic()
      end if
      if (present(maturity_factor)) then
        call single_trade_addon(work(i), supervisory, details(i), &
          maturity_factor=maturity_factor, simplified=use_simplified, &
          hedging_set=hedging_set_name(work(i)))
      else
        call single_trade_addon(work(i), supervisory, details(i), &
          simplified=use_simplified, hedging_set=hedging_set_name(work(i)))
      end if
    end do
    call move_alloc(details, result%trades)

    call unique_asset_classes(work, asset_names)
    do i = 1, size(asset_names)
      select case (uppercase(trim(asset_names(i))))
      case ("FX")
        call aggregate_fx(work, result%trades, supervisory, result, use_oem)
      case ("IRD")
        call aggregate_ird(work, result%trades, supervisory, result, &
          use_simplified, use_oem)
      case ("CREDIT")
        call aggregate_credit(work, result%trades, supervisory, result, &
          use_simplified, use_oem)
      case ("COMMODITY")
        call aggregate_commodity(work, result%trades, result, use_simplified, use_oem)
      case ("EQ")
        call aggregate_equity(work, result%trades, supervisory, result, &
          use_simplified, use_oem)
      case ("OTHEREXPOSURE")
        call aggregate_other(work, result%trades, result)
      end select
    end do

    result%addon = 0.0_dp
    do i = 1, size(result%asset_classes)
      result%addon = result%addon + result%asset_classes(i)%addon
    end do
    if (present(maturity_factor)) then
      result%maturity_factor = maturity_factor
      result%has_maturity_factor = .true.
    end if
  end subroutine calc_addon

  function hedging_set_name(trade) result(name)
    type(trade_t), intent(in) :: trade
    character(len=str_len) :: name
    character(len=str_len) :: inner
    character(len=str_len) :: prefix

    prefix = super_set_prefix(trade)
    select case (uppercase(trim(trade%trade_group)))
    case ("FX")
      inner = trade%ccy_pair
      if (len_trim(inner) == 0) inner = trade%currency
    case ("IRD")
      inner = trade%currency
    case ("CREDIT")
      inner = trade%reference_entity
    case ("COMMODITY")
      inner = trade%subclass
    case ("EQ")
      if (len_trim(trade%issuer) > 0) then
        inner = trade%issuer
      else
        inner = trade%underlying_instrument
      end if
    case ("OTHEREXPOSURE")
      inner = trade%subclass
    case default
      inner = trade%subclass
    end select

    if (len_trim(prefix) > 0 .and. len_trim(inner) > 0) then
      name = trim(prefix) // "_" // trim(inner)
    else if (len_trim(prefix) > 0) then
      name = prefix
    else
      name = inner
    end if
  end function hedging_set_name

  subroutine aggregate_fx(trades, details, records, result, oem)
    type(trade_t), intent(in) :: trades(:)
    type(single_trade_addon_t), intent(in) :: details(:)
    type(supervisory_record_t), intent(in) :: records(:)
    type(addon_result_t), intent(inout) :: result
    logical, intent(in) :: oem
    character(len=str_len), allocatable :: groups(:)
    type(asset_class_result_t) :: asset
    type(hedging_set_result_t) :: item
    real(dp) :: raw
    real(dp) :: factor
    integer :: i
    integer :: j

    asset%name = "FX"
    call unique_groups(trades, "FX", groups)
    do i = 1, size(groups)
      raw = 0.0_dp
      do j = 1, size(trades)
        if (.not. same_text(trades(j)%trade_group, "FX")) cycle
        if (.not. same_text(hedging_set_name(trades(j)), groups(i))) cycle
        if (oem) then
          raw = raw + trades(j)%notional * trades(j)%ei
        else
          raw = raw + details(j)%effective_notional
        end if
      end do
      if (oem) then
        factor = 0.04_dp
      else
        factor = calculate_factor_multiplier(groups(i)) * &
          supervisory_factor(records, "FX", "")
      end if
      item = hedging_set_result_t()
      item%asset_class = "FX"
      item%name = groups(i)
      item%effective_notional = raw
      item%supervisory_factor = factor
      item%addon = abs(factor * raw)
      asset%addon = asset%addon + item%addon
      call append_hedging_set(result%hedging_sets, item)
    end do
    call append_asset_class(result%asset_classes, asset)
  end subroutine aggregate_fx

  subroutine aggregate_other(trades, details, result)
    type(trade_t), intent(in) :: trades(:)
    type(single_trade_addon_t), intent(in) :: details(:)
    type(addon_result_t), intent(inout) :: result
    character(len=str_len), allocatable :: groups(:)
    type(asset_class_result_t) :: asset
    type(hedging_set_result_t) :: item
    real(dp) :: raw
    integer :: i
    integer :: j

    asset%name = "OtherExposure"
    call unique_groups(trades, "OtherExposure", groups)
    do i = 1, size(groups)
      raw = 0.0_dp
      do j = 1, size(trades)
        if (.not. same_text(trades(j)%trade_group, "OtherExposure")) cycle
        if (same_text(hedging_set_name(trades(j)), groups(i))) then
          raw = raw + details(j)%effective_notional
        end if
      end do
      item = hedging_set_result_t()
      item%asset_class = "OtherExposure"
      item%name = groups(i)
      item%effective_notional = raw
      item%supervisory_factor = 0.08_dp * calculate_factor_multiplier(groups(i))
      item%addon = abs(item%supervisory_factor * raw)
      asset%addon = asset%addon + item%addon
      call append_hedging_set(result%hedging_sets, item)
    end do
    call append_asset_class(result%asset_classes, asset)
  end subroutine aggregate_other

  subroutine aggregate_ird(trades, details, records, result, simplified, oem)
    type(trade_t), intent(in) :: trades(:)
    type(single_trade_addon_t), intent(in) :: details(:)
    type(supervisory_record_t), intent(in) :: records(:)
    type(addon_result_t), intent(inout) :: result
    logical, intent(in) :: simplified
    logical, intent(in) :: oem
    character(len=str_len), allocatable :: groups(:)
    type(asset_class_result_t) :: asset
    type(hedging_set_result_t) :: item
    real(dp) :: bucket(3)
    real(dp) :: factor
    real(dp) :: effective
    integer :: b
    integer :: i
    integer :: j

    asset%name = "IRD"
    call unique_groups(trades, "IRD", groups)
    do i = 1, size(groups)
      bucket = 0.0_dp
      do j = 1, size(trades)
        if (.not. same_text(trades(j)%trade_group, "IRD")) cycle
        if (.not. same_text(hedging_set_name(trades(j)), groups(i))) cycle
        b = trades(j)%set_time_bucket()
        if (oem) then
          bucket(b) = bucket(b) + trades(j)%notional * trades(j)%ei
        else
          bucket(b) = bucket(b) + details(j)%effective_notional
        end if
      end do

      if (.not. oem) call apply_digital_caps(trades, details, groups(i), bucket, records)
      if (simplified) then
        effective = sum(abs(bucket))
      else
        effective = sqrt(max(0.0_dp, bucket(1)**2 + bucket(2)**2 + bucket(3)**2 + &
          1.4_dp * bucket(1) * bucket(2) + 1.4_dp * bucket(2) * bucket(3) + &
          0.6_dp * bucket(1) * bucket(3)))
      end if
      if (oem) then
        factor = 0.005_dp
      else
        factor = calculate_factor_multiplier(groups(i)) * &
          supervisory_factor(records, "IRD", "")
      end if

      item = hedging_set_result_t()
      item%asset_class = "IRD"
      item%name = groups(i)
      item%effective_notional = effective
      item%supervisory_factor = factor
      item%addon = factor * effective
      asset%addon = asset%addon + item%addon
      call append_hedging_set(result%hedging_sets, item)
    end do
    call append_asset_class(result%asset_classes, asset)
  end subroutine aggregate_ird

  subroutine aggregate_credit(trades, details, records, result, simplified, oem)
    type(trade_t), intent(in) :: trades(:)
    type(single_trade_addon_t), intent(in) :: details(:)
    type(supervisory_record_t), intent(in) :: records(:)
    type(addon_result_t), intent(inout) :: result
    logical, intent(in) :: simplified
    logical, intent(in) :: oem
    character(len=str_len), allocatable :: groups(:)
    character(len=str_len) :: lookup_asset
    character(len=str_len) :: subclass
    type(asset_class_result_t) :: asset
    type(hedging_set_result_t) :: item
    real(dp) :: raw
    real(dp) :: rho
    real(dp) :: factor
    integer :: first
    integer :: i
    integer :: j

    asset%name = "Credit"
    call unique_groups(trades, "Credit", groups)
    do i = 1, size(groups)
      raw = 0.0_dp
      first = 0
      do j = 1, size(trades)
        if (.not. same_text(trades(j)%trade_group, "Credit")) cycle
        if (.not. same_text(hedging_set_name(trades(j)), groups(i))) cycle
        if (first == 0) first = j
        if (oem) then
          raw = raw + trades(j)%notional * trades(j)%ei
        else
          raw = raw + details(j)%effective_notional
        end if
      end do
      if (first == 0) cycle

      if (same_text(trades(first)%class_name, "CDS")) then
        lookup_asset = "CreditSingle"
      else
        lookup_asset = "CreditIndex"
      end if
      subclass = trades(first)%subclass
      if (oem) then
        factor = 0.06_dp
      else
        factor = calculate_factor_multiplier(groups(i)) * &
          supervisory_factor(records, lookup_asset, subclass)
      end if
      if (simplified) then
        rho = 1.0_dp
      else
        rho = supervisory_correlation(records, lookup_asset, subclass)
      end if

      item = hedging_set_result_t()
      item%asset_class = "Credit"
      item%name = groups(i)
      item%effective_notional = raw
      item%supervisory_factor = factor
      item%correlation = rho
      item%addon = factor * raw
      asset%systematic_component = asset%systematic_component + item%addon * rho
      asset%idiosyncratic_component = asset%idiosyncratic_component + &
        (1.0_dp - rho**2) * item%addon**2
      call append_hedging_set(result%hedging_sets, item)
    end do
    asset%systematic_component = asset%systematic_component**2
    asset%addon = sqrt(max(0.0_dp, asset%systematic_component + &
      asset%idiosyncratic_component))
    call append_asset_class(result%asset_classes, asset)
  end subroutine aggregate_credit

  subroutine aggregate_equity(trades, details, records, result, simplified, oem)
    type(trade_t), intent(in) :: trades(:)
    type(single_trade_addon_t), intent(in) :: details(:)
    type(supervisory_record_t), intent(in) :: records(:)
    type(addon_result_t), intent(inout) :: result
    logical, intent(in) :: simplified
    logical, intent(in) :: oem
    character(len=str_len), allocatable :: groups(:)
    character(len=str_len) :: subclass
    type(asset_class_result_t) :: asset
    type(hedging_set_result_t) :: item
    real(dp) :: raw
    real(dp) :: rho
    real(dp) :: factor
    integer :: first
    integer :: i
    integer :: j

    asset%name = "EQ"
    call unique_groups(trades, "EQ", groups)
    do i = 1, size(groups)
      raw = 0.0_dp
      first = 0
      do j = 1, size(trades)
        if (.not. same_text(trades(j)%trade_group, "EQ")) cycle
        if (.not. same_text(hedging_set_name(trades(j)), groups(i))) cycle
        if (first == 0) first = j
        raw = raw + details(j)%effective_notional
      end do
      if (first == 0) cycle

      if (same_text(trades(first)%subclass, "Index")) then
        subclass = "Index"
      else
        subclass = ""
      end if
      if (oem) then
        factor = 0.32_dp
      else
        factor = calculate_factor_multiplier(groups(i)) * &
          supervisory_factor(records, "EQ", subclass)
      end if
      if (simplified) then
        rho = 1.0_dp
      else
        rho = supervisory_correlation(records, "EQ", subclass)
      end if

      item = hedging_set_result_t()
      item%asset_class = "EQ"
      item%name = groups(i)
      item%effective_notional = raw
      item%supervisory_factor = factor
      item%correlation = rho
      item%addon = factor * raw
      asset%systematic_component = asset%systematic_component + item%addon * rho
      asset%idiosyncratic_component = asset%idiosyncratic_component + &
        (1.0_dp - rho**2) * item%addon**2
      call append_hedging_set(result%hedging_sets, item)
    end do
    asset%systematic_component = asset%systematic_component**2
    asset%addon = sqrt(max(0.0_dp, asset%systematic_component + &
      asset%idiosyncratic_component))
    call append_asset_class(result%asset_classes, asset)
  end subroutine aggregate_equity

  subroutine aggregate_commodity(trades, details, result, simplified, oem)
    type(trade_t), intent(in) :: trades(:)
    type(single_trade_addon_t), intent(in) :: details(:)
    type(addon_result_t), intent(inout) :: result
    logical, intent(in) :: simplified
    logical, intent(in) :: oem
    character(len=str_len), allocatable :: groups(:)
    character(len=str_len), allocatable :: commodities(:)
    type(asset_class_result_t) :: asset
    type(hedging_set_result_t) :: item
    real(dp) :: addon_sum
    real(dp) :: addon_square_sum
    real(dp) :: group_addon
    real(dp) :: raw
    real(dp) :: rho
    real(dp) :: factor
    integer :: first
    integer :: i
    integer :: j
    integer :: k

    asset%name = "Commodity"
    rho = 0.4_dp
    if (simplified) rho = 1.0_dp
    call unique_groups(trades, "Commodity", groups)
    do i = 1, size(groups)
      call unique_commodity_types(trades, groups(i), commodities)
      addon_sum = 0.0_dp
      addon_square_sum = 0.0_dp
      do j = 1, size(commodities)
        raw = 0.0_dp
        first = 0
        do k = 1, size(trades)
          if (.not. same_text(trades(k)%trade_group, "Commodity")) cycle
          if (.not. same_text(hedging_set_name(trades(k)), groups(i))) cycle
          if (.not. same_text(trades(k)%commodity_type, commodities(j))) cycle
          if (first == 0) first = k
          if (oem) then
            raw = raw + trades(k)%notional * trades(k)%ei
          else
            raw = raw + details(k)%effective_notional
          end if
        end do
        if (first == 0) cycle

        factor = calculate_factor_multiplier(groups(i))
        if (same_text(trades(first)%commodity_type, "Electricity")) then
          factor = factor * 0.40_dp
        else
          factor = factor * 0.18_dp
        end if
        item = hedging_set_result_t()
        item%asset_class = "Commodity"
        item%name = groups(i)
        item%subgroup = commodities(j)
        item%effective_notional = raw
        item%supervisory_factor = factor
        item%correlation = rho
        item%addon = abs(raw) * factor
        addon_sum = addon_sum + item%addon
        addon_square_sum = addon_square_sum + item%addon**2
        call append_hedging_set(result%hedging_sets, item)
      end do

      group_addon = sqrt(max(0.0_dp, (rho * addon_sum)**2 + &
        (1.0_dp - rho**2) * addon_square_sum))
      asset%addon = asset%addon + group_addon
    end do
    call append_asset_class(result%asset_classes, asset)
  end subroutine aggregate_commodity

  subroutine apply_digital_caps(trades, details, group, bucket, records)
    type(trade_t), intent(in) :: trades(:)
    type(single_trade_addon_t), intent(in) :: details(:)
    character(len=*), intent(in) :: group
    real(dp), intent(inout) :: bucket(3)
    type(supervisory_record_t), intent(in) :: records(:)
    character(len=str_len), allocatable :: identifiers(:)
    real(dp) :: cap
    real(dp) :: group_effective
    real(dp) :: factor
    integer :: b
    integer :: first
    integer :: i
    integer :: j

    call unique_digital_identifiers(trades, group, identifiers)
    factor = calculate_factor_multiplier(group) * supervisory_factor(records, "IRD", "")
    if (factor <= 0.0_dp) return

    do i = 1, size(identifiers)
      group_effective = 0.0_dp
      first = 0
      b = 0
      do j = 1, size(trades)
        if (.not. same_text(trades(j)%trade_group, "IRD")) cycle
        if (.not. same_text(hedging_set_name(trades(j)), group)) cycle
        if (.not. same_text(trades(j)%exotic_type, identifiers(i))) cycle
        if (first == 0) then
          first = j
          b = trades(j)%set_time_bucket()
        end if
        group_effective = group_effective + details(j)%effective_notional
      end do
      if (first == 0 .or. b == 0) cycle
      cap = abs(trades(first)%notional / factor)
      if (abs(group_effective) > cap) then
        bucket(b) = bucket(b) - sign(abs(group_effective) - cap, group_effective)
      end if
    end do
  end subroutine apply_digital_caps

  subroutine unique_asset_classes(trades, values)
    type(trade_t), intent(in) :: trades(:)
    character(len=str_len), allocatable, intent(out) :: values(:)
    integer :: i

    allocate(values(0))
    do i = 1, size(trades)
      if (len_trim(trades(i)%trade_group) == 0) cycle
      call append_unique(values, trades(i)%trade_group)
    end do
  end subroutine unique_asset_classes

  subroutine unique_groups(trades, asset_class, values)
    type(trade_t), intent(in) :: trades(:)
    character(len=*), intent(in) :: asset_class
    character(len=str_len), allocatable, intent(out) :: values(:)
    integer :: i

    allocate(values(0))
    do i = 1, size(trades)
      if (.not. same_text(trades(i)%trade_group, asset_class)) cycle
      call append_unique(values, hedging_set_name(trades(i)))
    end do
  end subroutine unique_groups

  subroutine unique_commodity_types(trades, group, values)
    type(trade_t), intent(in) :: trades(:)
    character(len=*), intent(in) :: group
    character(len=str_len), allocatable, intent(out) :: values(:)
    integer :: i

    allocate(values(0))
    do i = 1, size(trades)
      if (.not. same_text(trades(i)%trade_group, "Commodity")) cycle
      if (.not. same_text(hedging_set_name(trades(i)), group)) cycle
      call append_unique(values, trades(i)%commodity_type)
    end do
  end subroutine unique_commodity_types

  subroutine unique_digital_identifiers(trades, group, values)
    type(trade_t), intent(in) :: trades(:)
    character(len=*), intent(in) :: group
    character(len=str_len), allocatable, intent(out) :: values(:)
    integer :: i

    allocate(values(0))
    do i = 1, size(trades)
      if (.not. same_text(trades(i)%trade_group, "IRD")) cycle
      if (.not. same_text(hedging_set_name(trades(i)), group)) cycle
      if (len_trim(trades(i)%exotic_type) == 0) cycle
      call append_unique(values, trades(i)%exotic_type)
    end do
  end subroutine unique_digital_identifiers

  function super_set_prefix(trade) result(prefix)
    type(trade_t), intent(in) :: trade
    character(len=str_len) :: prefix

    prefix = ""
    if (index(uppercase(trim(trade%class_name)), "VOL") > 0) then
      prefix = "Vol_" // trim(trade%underlying_instrument)
    else if (trade%is_basis_swap()) then
      prefix = "Basis_" // trim(trade%pay_leg_ref) // "_" // trim(trade%rec_leg_ref)
    end if
  end function super_set_prefix

  subroutine append_unique(values, value)
    character(len=str_len), allocatable, intent(inout) :: values(:)
    character(len=*), intent(in) :: value
    character(len=str_len), allocatable :: work(:)
    integer :: i
    integer :: n

    do i = 1, size(values)
      if (same_text(values(i), value)) return
    end do
    n = size(values)
    allocate(work(n + 1))
    if (n > 0) work(:n) = values
    work(n + 1) = value
    call move_alloc(work, values)
  end subroutine append_unique

  subroutine append_hedging_set(values, value)
    type(hedging_set_result_t), allocatable, intent(inout) :: values(:)
    type(hedging_set_result_t), intent(in) :: value
    type(hedging_set_result_t), allocatable :: work(:)
    integer :: n

    n = size(values)
    allocate(work(n + 1))
    if (n > 0) work(:n) = values
    work(n + 1) = value
    call move_alloc(work, values)
  end subroutine append_hedging_set

  subroutine append_asset_class(values, value)
    type(asset_class_result_t), allocatable, intent(inout) :: values(:)
    type(asset_class_result_t), intent(in) :: value
    type(asset_class_result_t), allocatable :: work(:)
    integer :: n

    n = size(values)
    allocate(work(n + 1))
    if (n > 0) work(:n) = values
    work(n + 1) = value
    call move_alloc(work, values)
  end subroutine append_asset_class

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

end module saccr_addon
