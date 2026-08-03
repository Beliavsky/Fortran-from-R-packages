! SPDX-License-Identifier: GPL-2.0-only
program rates_curve_example
   use fmbasics
   implicit none
   type(interest_rate_t) :: rate, zero
   type(discount_factor_t) :: df
   type(zero_curve_t) :: curve
   integer :: start_date, end_date, status

   start_date = make_date(2020, 1, 1)
   end_date = make_date(2025, 1, 1)
   rate = interest_rate(0.04_dp, 2.0_dp, 'act/365', status)
   df = as_discount_factor(rate, start_date, end_date, status)
   write(*,'(a,f12.8)') 'Five-year discount factor: ', df%value(1)

   curve = build_zero_curve(logdf_interpolation(), 'data/zerocurve.csv', status)
   if (status /= FM_OK) error stop 'could not load zero curve'
   zero = interpolate_zeros(curve, [make_date(2018, 12, 31)], status=status)
   write(*,'(a,f10.6)') 'Interpolated continuous zero rate: ', zero%value(1)
end program rates_curve_example
