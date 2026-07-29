! stockAnalyst-fortran
! Copyright (C) 2022 MaheshP Kumar (original R package)
! Copyright (C) 2026 Fortran port contributors
! SPDX-License-Identifier: GPL-3.0-only
!
! This file is part of stockAnalyst-fortran.
! It is free software: you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software
! Foundation, version 3 of the License.

program reference_examples
  use stock_analyst, only : dp, share_value_using_ddm_n_years, terminal_value_using_pe, &
    share_price_using_past_pe, computing_ev_multiple
  implicit none

  print '(a,f10.2)', 'Multi-period DDM: ', share_value_using_ddm_n_years( &
    [3.0_dp, 3.15_dp], 40.0_dp, [1.0_dp, 2.0_dp], 2.0_dp, 0.08_dp)
  print '(a,f10.2)', 'GGM terminal value: ', terminal_value_using_pe( &
    'GGM', 14.3_dp, 3.0_dp, 0.45_dp, 0.0715_dp, 0.10_dp)
  print '(a,f10.0)', 'Median historical P/E value: ', share_price_using_past_pe( &
    'median', [15.8_dp, 23.1_dp, 10.0_dp, 19.8_dp, 35.8_dp], 203.71_dp)
  print '(a,f10.1)', 'EV/EBITDA: ', computing_ev_multiple('EBITDA', 14411.0_dp, 3320.0_dp, 18962.0_dp)
end program reference_examples
