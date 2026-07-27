! SPDX-License-Identifier: GPL-2.0-only
! Copyright (C) 2010-2023 Yang Lu and David Kane
! Copyright (C) 2026 Modern Fortran translation contributors
! This program is free software under GNU GPL version 2 only.
program demo_pa
  use pa
  implicit none
  integer :: category(12), period(12), i
  real(dp) :: wb(12), wp(12), ret(12)
  type(brinson_multi_result) :: model
  type(attribution_summary) :: summary

  category = [1,1,2,2,3,3, 1,1,2,2,3,3]
  period = [1,1,1,1,1,1, 2,2,2,2,2,2]
  wb = [0.15_dp,0.15_dp,0.20_dp,0.20_dp,0.15_dp,0.15_dp, &
        0.10_dp,0.20_dp,0.15_dp,0.25_dp,0.10_dp,0.20_dp]
  wp = [0.20_dp,0.10_dp,0.25_dp,0.15_dp,0.10_dp,0.20_dp, &
        0.15_dp,0.20_dp,0.20_dp,0.20_dp,0.15_dp,0.10_dp]
  ret = [0.02_dp,0.01_dp,-0.01_dp,0.03_dp,0.04_dp,0.00_dp, &
         0.01_dp,0.025_dp,0.02_dp,-0.005_dp,0.03_dp,0.015_dp]

  call fit_brinson_multi(period, category, wb, wp, ret, model)
  call summarize_brinson_multi(model, 'geometric', summary)

  print '(a)', 'Brinson raw attribution by period'
  print '(a)', 'period       allocation       selection     interaction          active'
  do i = 1, size(model%period)
    print '(i6,4f16.8)', model%period(i), model%raw(:,i)
  end do
  print '(a,4f16.8)', 'geometric ', summary%aggregate
end program demo_pa
