! SPDX-License-Identifier: GPL-2.0-only
program credit_vol_example
   use fmbasics
   implicit none
   type(zero_curve_t) :: zero_curve_data
   type(cds_curve_t) :: cds
   type(cds_spec_t) :: spec
   type(survival_probabilities_t) :: survival
   type(vol_surface_t) :: surface
   real(dp), allocatable :: vol(:)
   integer :: status

   zero_curve_data = build_zero_curve(logdf_interpolation(), 'data/zerocurve.csv', status)
   if (status /= FM_OK) error stop 'could not load zero curve'

   spec = cds_spec('Example')
   cds = cds_curve(zero_curve_data%reference_date, [1.0_dp, 3.0_dp, 5.0_dp], &
      [0.005_dp, 0.008_dp, 0.011_dp], 0.6_dp, 4, spec, status)
   survival = bootstrap_cds_survival(cds, zero_curve_data, status=status)
   if (status /= FM_OK) error stop 'CDS bootstrap failed'
   write(*,'(a,3(1x,f10.7))') 'Survival probabilities:', survival%value

   surface = build_vol_surface('data/volsurface.csv', status)
   vol = interpolate_vol(surface, [make_date(2022, 6, 10)], [90.0_dp], status)
   if (status /= FM_OK) error stop 'volatility interpolation failed'
   write(*,'(a,f10.6)') 'Interpolated volatility: ', vol(1)
end program credit_vol_example
