! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Modern Fortran translation of R package classInt computational code.
module classint
    use classint_kinds, only: dp
    use classint_types, only: classint_options, class_intervals, jenks_indices, hclust_1d_model
    use classint_core, only: classint_fit, get_hclust_class_intervals, get_bclust_class_intervals
    use classint_metrics, only: find_cols, class_counts, jenks_tests, gvf, tai, oai, &
                                classint_loglik, classint_aic, classint_n_partitions
    use classint_fisher, only: fisher_exact, jenks_breaks
    use classint_dpih, only: dpih_bandwidth
    use classint_pretty, only: pretty_breaks
    use classint_api, only: classify_intervals
    implicit none
    private

    public :: dp
    public :: classint_options, class_intervals, jenks_indices, hclust_1d_model
    public :: classint_fit, get_hclust_class_intervals, get_bclust_class_intervals
    public :: find_cols, classify_intervals, class_counts, jenks_tests
    public :: gvf, tai, oai, classint_loglik, classint_aic, classint_n_partitions
    public :: fisher_exact, jenks_breaks, dpih_bandwidth, pretty_breaks
end module classint
