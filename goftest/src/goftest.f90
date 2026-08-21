! SPDX-License-Identifier: GPL-2.0-or-later
module goftest
    use goftest_kinds, only : dp
    use goftest_ad, only : ad_cdf, ad_quantile, ad_statistic, ad_test_uniform, ad_inf_fast, ad_inf_exact
    use goftest_cvm, only : cvm_cdf, cvm_quantile, cvm_statistic, cvm_test_uniform
    use goftest_api, only : p_ad, q_ad, p_cvm, q_cvm
    use goftest_testing, only : gof_result, ad_test, cvm_test, ad_test_values, cvm_test_values
    use goftest_names, only : recognise_cdf
    implicit none
    public
end module goftest
