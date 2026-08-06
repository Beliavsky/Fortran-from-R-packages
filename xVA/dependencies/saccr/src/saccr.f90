module saccr
  use trading, only : dp, str_len, trade_t, csa_t, collateral_t, &
    parse_trades_csv, load_csa_csv, load_collateral_csv, select_derivatives
  use saccr_types
  use saccr_supervisory
  use saccr_core
  use saccr_addon
  use saccr_portfolio
  use saccr_io
  use saccr_examples
  implicit none
  public
end module saccr
