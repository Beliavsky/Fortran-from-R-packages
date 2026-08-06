! SPDX-License-Identifier: GPL-3.0-only
module spantest
  use spantest_kinds, only : dp, pi_dp
  use spantest_types
  use spantest_classical, only : span_bj, span_f1, span_f2, span_grs, span_hk, span_km, span_py
  use spantest_gl, only : span_gl_a, span_gl_ad
  use spantest_as, only : span_as, cauchy_pvalue
  use spantest_simulation, only : span_simulate, garch_filter, standardized_skew_t
  implicit none
  public
end module spantest
