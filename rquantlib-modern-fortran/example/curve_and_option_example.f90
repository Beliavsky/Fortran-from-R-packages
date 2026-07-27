! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 OpenAI
! This file is part of a modern Fortran translation of RQuantLib.
! It is free software under GNU GPL version 2 or any later version.
program curve_and_option_example
  use rq_kinds, only: dp
  use rq_curves
  use rq_options
  implicit none
  type(discount_curve_t) :: curve
  type(option_result) :: option
  real(dp) :: rate
  call make_flat_curve(0.04_dp,10.0_dp,0.25_dp,curve)
  rate=curve%zero_rate(2.0_dp)
  option=european_option('call',100.0_dp,105.0_dp,0.01_dp,rate,2.0_dp,0.22_dp)
  write(*,'(a,f10.6)') 'two-year zero rate: ',rate
  write(*,'(a,f10.6)') 'call value:         ',option%value
end program curve_and_option_example
