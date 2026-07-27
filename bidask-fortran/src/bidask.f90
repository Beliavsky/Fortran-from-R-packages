! SPDX-License-Identifier: MIT
! Based on bidask 2.1.5, Copyright (c) 2024 Emanuele Guidotti.
module bidask
  use bidask_kinds, only: dp
  use bidask_types, only: ohlc_data, spread_result, spread_series_result
  use bidask_estimators, only: edge, edge_estimate, ar_estimate, cs_estimate, &
    roll_estimate, ohlc_estimate, estimate_method
  use bidask_windows, only: edge_rolling, edge_expanding, spread, &
    spread_expanding, spread_endpoints
  use bidask_simulation, only: simulate_ohlc, sim
  implicit none
  public
end module bidask
