! bondAnalyst modern Fortran port
! Copyright (C) 2022 MaheshP Kumar
! Copyright (C) 2026 Fortran port contributors
! SPDX-License-Identifier: GPL-3.0-only
module bondanalyst
   use bondanalyst_kinds, only : dp
   use bondanalyst_support, only : ba_success, ba_invalid_argument, &
      ba_size_mismatch, ba_no_root, ba_out_of_range
   use bondanalyst_valuation
   use bondanalyst_rates
   use bondanalyst_money_market
   use bondanalyst_duration
   implicit none
   public

end module bondanalyst
