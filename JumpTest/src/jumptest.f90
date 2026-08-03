! SPDX-License-Identifier: MIT
module jumptest
  use jumptest_kinds, only : dp, i8, pi
  use jumptest_status, only : JT_SUCCESS, JT_INVALID_ARGUMENT, JT_INVALID_DIMENSION, &
    JT_NONFINITE_INPUT, JT_DEGENERATE_SAMPLE, JT_NUMERICAL_FAILURE, status_message
  use jumptest_statistics, only : statp_result, adjp_result, pcombine_result, &
    METHOD_BNS, METHOD_AMED, METHOD_AMIN, bns_statistic, amin_statistic, &
    amed_statistic, jumptestday, jumptestperiod, pcombine, ppool, bh_adjust
  use jumptest_simulation, only : simulation_result, sv, svj, sv1f, sv1fj, sv2f, &
    lp_path, pvc_path, pv2_path
  implicit none
  private

  public :: dp, i8, pi
  public :: JT_SUCCESS, JT_INVALID_ARGUMENT, JT_INVALID_DIMENSION
  public :: JT_NONFINITE_INPUT, JT_DEGENERATE_SAMPLE, JT_NUMERICAL_FAILURE
  public :: status_message
  public :: statp_result, adjp_result, pcombine_result, simulation_result
  public :: METHOD_BNS, METHOD_AMED, METHOD_AMIN
  public :: bns_statistic, amin_statistic, amed_statistic
  public :: jumptestday, jumptestperiod, pcombine, ppool, bh_adjust
  public :: sv, svj, sv1f, sv1fj, sv2f
  public :: lp_path, pvc_path, pv2_path

end module jumptest
