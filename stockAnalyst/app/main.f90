! stockAnalyst-fortran
! Copyright (C) 2022 MaheshP Kumar (original R package)
! Copyright (C) 2026 Fortran port contributors
! SPDX-License-Identifier: GPL-3.0-only
!
! This file is part of stockAnalyst-fortran.
! It is free software: you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software
! Foundation, version 3 of the License.

program stock_analyst_demo
  use stock_analyst, only : dp, share_value_ggm_constant_growth, computing_r_with_capm, &
    computing_wacc, share_value_ri
  implicit none

  print '(a,f10.2)', 'Gordon-growth share value: ', &
    share_value_ggm_constant_growth(1.1024_dp, 0.101_dp, 0.06_dp, 1)
  print '(a,f10.4)', 'CAPM required return:      ', &
    computing_r_with_capm(0.049_dp, 0.74_dp, 0.045_dp)
  print '(a,f10.5)', 'Weighted cost of capital: ', &
    computing_wacc(35.0_dp, 65.0_dp, 0.056_dp, 0.127_dp, 0.29_dp)
  print '(a,f10.2)', 'Residual-income value:    ', &
    share_value_ri(6.0_dp, [1.4_dp, 1.8_dp, 3.175_dp], 0.10_dp, [1.0_dp, 2.0_dp, 3.0_dp])
end program stock_analyst_demo
