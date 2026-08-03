! SPDX-License-Identifier: GPL-2.0-only
program test_curves_credit
   use fmbasics
   implicit none
   type(zero_curve_t) :: zc
   type(interest_rate_t) :: zeros, fwds
   type(discount_factor_t) :: dfs
   type(cds_spec_t) :: spec
   type(zero_hazard_rate_t) :: hazard, forward_hazard
   type(survival_probabilities_t) :: survival, boot
   type(credit_curve_t) :: cc
   type(cds_curve_t) :: cds
   real(dp), allocatable :: v(:)
   integer :: status
   integer :: dates(2), froms(1), tos(1)

   zc = build_zero_curve(logdf_interpolation(), 'data/zerocurve.csv', status)
   call check(status == FM_OK .and. zc%size() == 27, 'load zero curve')
   v = interpolate_zero(zc, [0.0_dp, 4.0_dp, 6.0_dp, 50.0_dp], status)
   call check(status == FM_OK .and. all(v > 0.0_dp), 'zero interpolation')

   dates = [make_date(2016,12,31), make_date(2017,12,31)]
   zeros = interpolate_zeros(zc, dates, status=status)
   call check(status == FM_OK .and. zeros%size() == 2, 'date zero interpolation')
   froms = [make_date(2016,12,31)]
   tos = [make_date(2017,12,31)]
   dfs = interpolate_dfs(zc, froms, tos, status)
   fwds = interpolate_fwds(zc, froms, tos, status)
   call check(status == FM_OK .and. dfs%value(1) > 0.0_dp, 'forward df')
   call check(fwds%value(1) > -1.0_dp, 'forward rate')

   spec = cds_spec('Empty')
   hazard = zero_hazard_rate(0.04_dp, 0.0_dp, 'act/360', spec)
   survival = as_survival_probabilities(hazard, [make_date(2010,1,1)], [make_date(2015,1,1)])
   call check(abs(survival%value(1)-0.831331978570109_dp) < 1.0e-14_dp, 'hazard survival')
   hazard = as_zero_hazard_rate(survival, 0.0_dp, 'act/360')
   call check(abs(hazard%value(1)-0.04_dp) < 1.0e-13_dp, 'survival hazard')

   cds = cds_curve(zc%reference_date, [1.0_dp, 3.0_dp, 5.0_dp, 7.0_dp], &
      [0.005_dp, 0.007_dp, 0.009_dp, 0.011_dp], 0.6_dp, 4, spec, status)
   boot = bootstrap_cds_survival(cds, zc, status=status)
   call check(status == FM_OK, 'CDS bootstrap')
   call check(all(boot%value > 0.0_dp) .and. all(boot%value <= 1.0_dp), 'bootstrap range')
   call check(all(boot%value(2:) < boot%value(:boot%size()-1)), 'bootstrap monotonic')

   cc = credit_curve(boot, zc%reference_date, cubic_interpolation(), status=status)
   v = interpolate_credit(cc, [0.5_dp, 2.0_dp, 6.0_dp], status)
   call check(status == FM_OK .and. all(v >= 0.0_dp), 'credit interpolation')
   forward_hazard = interpolate_credit_fwds(cc, froms, tos, status)
   call check(status == FM_OK .and. forward_hazard%value(1) >= 0.0_dp, 'forward hazard')

   print '(a)', 'test_curves_credit: PASS'

contains

   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write(*,'(a,1x,a)') 'FAIL:', trim(label)
         error stop 1
      end if
   end subroutine check

end program test_curves_credit
