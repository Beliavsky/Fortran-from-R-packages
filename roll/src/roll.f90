! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from the R package roll by Jason Foster; see PROVENANCE.md and NOTICE.md.
module roll
   use r_kinds, only : dp
   use roll_univariate, only : roll_na_logical
   use roll_univariate, only : roll_any_vec, roll_all_vec, roll_sum_vec, roll_prod_vec, roll_mean_vec
   use roll_univariate, only : roll_min_vec, roll_max_vec, roll_idxmin_vec, roll_idxmax_vec
   use roll_univariate, only : roll_median_vec, roll_quantile_vec, roll_var_vec, roll_sd_vec, roll_scale_vec
   use roll_matrix, only : roll_any_mat, roll_all_mat, roll_sum_mat, roll_prod_mat, roll_mean_mat
   use roll_matrix, only : roll_min_mat, roll_max_mat, roll_idxmin_mat, roll_idxmax_mat
   use roll_matrix, only : roll_median_mat, roll_quantile_mat, roll_var_mat, roll_sd_mat, roll_scale_mat
   use roll_multivariate, only : roll_lm_result
   use roll_multivariate, only : roll_cov_vec, roll_cor_vec, roll_crossprod_vec, roll_lm_vec
   use roll_multivariate, only : roll_cov_mat, roll_cor_mat, roll_crossprod_mat, roll_lm_mat
   implicit none
   private

   public :: dp, roll_na_logical, roll_lm_result
   public :: roll_any, roll_all, roll_sum, roll_prod, roll_mean
   public :: roll_min, roll_max, roll_idxmin, roll_idxmax, roll_median, roll_quantile
   public :: roll_var, roll_sd, roll_scale, roll_cov, roll_cor, roll_crossprod, roll_lm

   interface roll_any
      module procedure roll_any_vec
      module procedure roll_any_mat
   end interface roll_any

   interface roll_all
      module procedure roll_all_vec
      module procedure roll_all_mat
   end interface roll_all

   interface roll_sum
      module procedure roll_sum_vec
      module procedure roll_sum_mat
   end interface roll_sum

   interface roll_prod
      module procedure roll_prod_vec
      module procedure roll_prod_mat
   end interface roll_prod

   interface roll_mean
      module procedure roll_mean_vec
      module procedure roll_mean_mat
   end interface roll_mean

   interface roll_min
      module procedure roll_min_vec
      module procedure roll_min_mat
   end interface roll_min

   interface roll_max
      module procedure roll_max_vec
      module procedure roll_max_mat
   end interface roll_max

   interface roll_idxmin
      module procedure roll_idxmin_vec
      module procedure roll_idxmin_mat
   end interface roll_idxmin

   interface roll_idxmax
      module procedure roll_idxmax_vec
      module procedure roll_idxmax_mat
   end interface roll_idxmax

   interface roll_median
      module procedure roll_median_vec
      module procedure roll_median_mat
   end interface roll_median

   interface roll_quantile
      module procedure roll_quantile_vec
      module procedure roll_quantile_mat
   end interface roll_quantile

   interface roll_var
      module procedure roll_var_vec
      module procedure roll_var_mat
   end interface roll_var

   interface roll_sd
      module procedure roll_sd_vec
      module procedure roll_sd_mat
   end interface roll_sd

   interface roll_scale
      module procedure roll_scale_vec
      module procedure roll_scale_mat
   end interface roll_scale

   interface roll_cov
      module procedure roll_cov_vec
      module procedure roll_cov_mat
   end interface roll_cov

   interface roll_cor
      module procedure roll_cor_vec
      module procedure roll_cor_mat
   end interface roll_cor

   interface roll_crossprod
      module procedure roll_crossprod_vec
      module procedure roll_crossprod_mat
   end interface roll_crossprod

   interface roll_lm
      module procedure roll_lm_vec
      module procedure roll_lm_mat
   end interface roll_lm

end module roll
