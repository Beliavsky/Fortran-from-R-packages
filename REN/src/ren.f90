! SPDX-License-Identifier: AGPL-3.0-or-later
! Derived from REN 0.1.0 computational code; see NOTICE.md.
module ren
  use ren_kinds, only : dp, i8
  use ren_types
  use ren_portfolio, only : insert_at, po_avg, po_gross_exp, po_cov_shrink, po_cols, po_jm, &
    buh_clust, po_bhu, po_tzt, po_sw, po_sw_lasso
  use ren_analysis, only : prepare_data, perform_analysis, ren_run
  implicit none
  private
  public :: dp, i8
  public :: asset_group, cluster_result, prepared_data_type, analysis_options, analysis_result
  public :: ren_success, ren_invalid_argument, ren_dimension_error, ren_numerical_error, ren_dependency_error
  public :: ren_method_count, ren_method_names, ren_status_message
  public :: insert_at, po_avg, po_gross_exp, po_cov_shrink, po_cols, po_jm
  public :: buh_clust, po_bhu, po_tzt, po_sw, po_sw_lasso
  public :: prepare_data, perform_analysis, ren_run
end module ren
