! SPDX-License-Identifier: GPL-2.0-only
program test_vol_money
   use fmbasics
   implicit none
   type(vol_surface_t) :: surface
   type(single_currency_money_t) :: a, b, c
   type(multi_currency_money_t) :: multi, aggregated
   type(cash_flow_t) :: flow
   real(dp), allocatable :: vol(:)
   real(dp), parameter :: expected(6) = [ &
      0.6076543447950257_dp, 0.26853916752886564_dp, &
      0.19909034016558932_dp, 0.25769535624031686_dp, &
      0.25855784359768552_dp, 0.26647898600000003_dp]
   integer :: maturity(6), status

   surface = build_vol_surface('data/volsurface.csv', status)
   call check(status == FM_OK, 'load vol surface')
   maturity = [make_date(2023,8,15), make_date(2023,10,10), make_date(2020,2,29), &
      make_date(2021,4,15), make_date(2022,6,10), make_date(2025,6,10)]
   vol = interpolate_vol(surface, maturity, [3.0_dp, 96.0_dp, 150.0_dp, &
      80.0_dp, 90.0_dp, 300.0_dp], status)
   call check(status == FM_OK, 'vol interpolation status')
   call check(maxval(abs(vol-expected)) < 5.0e-11_dp, 'vol interpolation values')

   a = single_currency_money(1.0_dp, aud())
   b = single_currency_money(2.0_dp, usd())
   c = single_currency_money(3.0_dp, aud())
   multi = multi_currency_money([a, b, c], status)
   call check(status == FM_OK .and. multi%size() == 3, 'multi money')
   aggregated = aggregate_by_currency(multi)
   call check(aggregated%size() == 2, 'money aggregation count')
   call check(abs(aggregated%value(1)-4.0_dp) < 1.0e-15_dp, 'money aggregation value')
   flow = cash_flow([make_date(2017,11,15), make_date(2017,11,16), &
      make_date(2017,11,17)], multi, status)
   call check(status == FM_OK .and. flow%size() == 3, 'cash flow')

   print '(a)', 'test_vol_money: PASS'

contains

   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write(*,'(a,1x,a)') 'FAIL:', trim(label)
         error stop 1
      end if
   end subroutine check

end program test_vol_money
