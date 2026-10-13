! SPDX-License-Identifier: GPL-2.0-only
module propcis
  !! Public API for all 21 computational exports of PropCIs 0.3-0.
  use propcis_common, only: dp, ci_result, ci_ok, ci_invalid, ci_numerical
  use propcis_frequentist
  use propcis_bayesian
  implicit none
  public
end module propcis
