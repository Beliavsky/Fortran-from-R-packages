! SPDX-License-Identifier: GPL-2.0-or-later
! Based on BCC1997 0.1.1, Copyright (C) 2017 Haoran Zhang.
module bcc1997
   use bcc1997_kinds, only : dp
   use bcc1997_types, only : bcc_parameters, integration_settings, bcc_result
   use bcc1997_model, only : bcc, bcc_price, bcc_price_strikes, &
      bcc_characteristic_1, bcc_characteristic_2, black_scholes_price, &
      validate_parameters
   implicit none
   public
end module bcc1997
