module saccr_io
  use trading, only : dp, trade_t, csa_t, collateral_t, parse_trades_csv, &
    load_csa_csv, load_collateral_csv, select_derivatives
  use saccr_types, only : portfolio_result_t
  use saccr_portfolio, only : calculate_portfolio
  use saccr_core, only : determine_ccr_methodology
  implicit none
  private

  public :: saccr_calculator
  public :: determine_ccr_methodology_csv
  public :: write_exposure_csv
  public :: write_portfolio_summary

contains

  subroutine saccr_calculator(trades_filename, result, csa_filename, &
      collateral_filename, simplified, oem, ignore_margin)
    character(len=*), intent(in) :: trades_filename
    type(portfolio_result_t), intent(out) :: result
    character(len=*), intent(in), optional :: csa_filename
    character(len=*), intent(in), optional :: collateral_filename
    logical, intent(in), optional :: simplified
    logical, intent(in), optional :: oem
    logical, intent(in), optional :: ignore_margin
    type(trade_t), allocatable :: all_trades(:)
    type(trade_t), allocatable :: trades(:)
    type(csa_t), allocatable :: agreements(:)
    type(collateral_t), allocatable :: collateral(:)
    logical :: margin_ignored
    logical :: use_oem
    logical :: use_simplified

    use_simplified = .false.
    if (present(simplified)) use_simplified = simplified
    use_oem = .false.
    if (present(oem)) use_oem = oem
    if (use_oem) use_simplified = .true.
    margin_ignored = .false.
    if (present(ignore_margin)) margin_ignored = ignore_margin

    call parse_trades_csv(trades_filename, all_trades)
    call select_derivatives(all_trades, trades)

    if (present(csa_filename)) then
      call load_csa_csv(csa_filename, agreements)
    else
      allocate(agreements(0))
    end if
    if (present(collateral_filename)) then
      call load_collateral_csv(collateral_filename, collateral)
    else
      allocate(collateral(0))
    end if

    call calculate_portfolio(trades, result, csas=agreements, &
      collaterals=collateral, simplified=use_simplified, oem=use_oem, &
      ignore_margin=margin_ignored)
  end subroutine saccr_calculator


  function determine_ccr_methodology_csv(trades_filename, total_assets) result(methodology)
    character(len=*), intent(in) :: trades_filename
    real(dp), intent(in) :: total_assets
    character(len=256) :: methodology
    type(trade_t), allocatable :: all_trades(:)
    type(trade_t), allocatable :: trades(:)

    call parse_trades_csv(trades_filename, all_trades)
    call select_derivatives(all_trades, trades)
    methodology = determine_ccr_methodology(trades, total_assets)
  end function determine_ccr_methodology_csv

  subroutine write_exposure_csv(filename, result)
    character(len=*), intent(in) :: filename
    type(portfolio_result_t), intent(in) :: result
    integer :: i
    integer :: ios
    integer :: unit

    open(newunit=unit, file=filename, status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "write_exposure_csv: unable to open file"
    write(unit, '(a)') "Counterparty,CSA_ID,Margined,Addon,EAD,PFE," // &
      "Replacement_Cost_VC,Replacement_Cost"
    do i = 1, size(result%exposures)
      write(unit, '(a,",",a,",",l1,5(",",es24.16))') &
        trim(result%exposures(i)%counterparty), &
        trim(result%exposures(i)%csa_id), &
        result%exposures(i)%margined, &
        result%exposures(i)%addon%addon, &
        result%exposures(i)%ead, &
        result%exposures(i)%pfe, &
        result%exposures(i)%replacement_cost%v_c, &
        result%exposures(i)%replacement_cost%rc
    end do
    close(unit)
  end subroutine write_exposure_csv

  subroutine write_portfolio_summary(result, unit)
    type(portfolio_result_t), intent(in) :: result
    integer, intent(in), optional :: unit
    integer :: i
    integer :: output_unit

    output_unit = 6
    if (present(unit)) output_unit = unit
    write(output_unit, '(a,i0)') "exposure_sets=", size(result%exposures)
    do i = 1, size(result%exposures)
      write(output_unit, '(a,i0)') "exposure=", i
      write(output_unit, '(a)') "counterparty=" // &
        trim(result%exposures(i)%counterparty)
      if (len_trim(result%exposures(i)%csa_id) > 0) then
        write(output_unit, '(a)') "csa_id=" // trim(result%exposures(i)%csa_id)
      end if
      write(output_unit, '(a,es16.8)') "addon=", result%exposures(i)%addon%addon
      write(output_unit, '(a,es16.8)') "replacement_cost=", &
        result%exposures(i)%replacement_cost%rc
      write(output_unit, '(a,es16.8)') "pfe=", result%exposures(i)%pfe
      write(output_unit, '(a,es16.8)') "ead=", result%exposures(i)%ead
    end do
    write(output_unit, '(a,es16.8)') "total_ead=", result%total_ead
  end subroutine write_portfolio_summary

end module saccr_io
