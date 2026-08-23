! SPDX-License-Identifier: GPL-3.0-only
module bigstatsr
    use bigstatsr_kinds, only: dp
    use bigstatsr_fbm, only: fbm_real, fbm_code256, create_fbm, attach_fbm, &
        create_code256, attach_code256, fbm_from_array, fbm_copy, fbm_increment, fbm_transpose
    use bigstatsr_matrix, only: colstats_result, scaling_result, big_colstats, big_scale, &
        big_prod_vec, big_cprod_vec, big_prod_mat, big_cprod_mat, big_crossprod_self, &
        big_tcrossprod_self, big_cor, big_counts_rows, big_counts_cols
    use bigstatsr_stats, only: auc_boot_result, pcor_result, auc, auc_sorted, &
        auc_sorted_weighted, auc_bootstrap, pcor
    use bigstatsr_regression, only: univ_reg_result, enet_path_result, summary_result, &
        big_univ_linreg, big_univ_logreg, elastic_net_gaussian_path, &
        elastic_net_logistic_path, predict_enet, big_summaries
    use bigstatsr_svd, only: big_svd_result, big_svd, big_random_svd, svd_predict
    use bigstatsr_misc, only: block_size, rows_along, cols_along, get_beta
    implicit none
    private
    public :: dp
    public :: fbm_real, fbm_code256, create_fbm, attach_fbm, create_code256, attach_code256
    public :: fbm_from_array, fbm_copy, fbm_increment, fbm_transpose
    public :: colstats_result, scaling_result, big_colstats, big_scale
    public :: big_prod_vec, big_cprod_vec, big_prod_mat, big_cprod_mat
    public :: big_crossprod_self, big_tcrossprod_self, big_cor
    public :: big_counts_rows, big_counts_cols
    public :: auc_boot_result, pcor_result, auc, auc_sorted, auc_sorted_weighted, auc_bootstrap, pcor
    public :: univ_reg_result, enet_path_result, summary_result
    public :: big_univ_linreg, big_univ_logreg, elastic_net_gaussian_path
    public :: elastic_net_logistic_path, predict_enet, big_summaries
    public :: big_svd_result, big_svd, big_random_svd, svd_predict
    public :: block_size, rows_along, cols_along, get_beta
end module bigstatsr
