! helper functions for R-to-Fortran transpiler
module r_mod
use, intrinsic :: iso_fortran_env, only: real64, int64
use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf, &
   & ieee_is_finite
#ifdef XR2F_USE_R_RNG
use, intrinsic :: iso_c_binding, only: c_double, c_int
#endif
implicit none
private
public :: dp, runif1, runif_vec, rnorm1, rnorm_vec, rexp, rgamma, rbeta, rchisq, rt_vec, rf_rng, rlogis, rlnorm, &
   & rweibull, rcauchy, rgeom, rnbinom, rhyper, rwilcox, rsignrank, rmultinom, rnorm_mat, rbinom, rpois, &
   & random_choice2_prob, &
   & randint_range, sample_int, sample_int1, sample_value_int, quantile, median, summary, dnorm, tail, cbind2, &
   & cbind, numeric, &
   & pmax, r_round, sd, r_sd, var, r_format_vec, colMeans, apply_col_cumsum, apply_col_sd, apply_row_sd, &
   & count_ws_tokens, &
   & besselJ, besselY, besselI, besselK, &
   & read_real_vector, read_table_real_matrix, read_csv_real_matrix, read_csv_header_names, &
   & write_table_real_matrix, write_table_real_vector, lm_fit_t, glm_fit_t, prcomp_fit_t, eigen_result_t, &
   & optim_result_t, optim_bfgs, optim_lbfgsb, optim_cg, optim_sann, optim_nelder_mead, constr_optim_bfgs, &
   & constr_optim_nelder_mead, nlm_result_t, nlm_stub, nlm_optimize_scalar, nlm_optimize_vec, set_nlm_method, &
   & print_nlm_result, integrate_result_t, integrate, print_integrate_result, hist_result_t, hist, print_hist, &
   & decompose_result_t, ks_test_result_t, lm_fit_general, lm_r_squared_general, lm_predict_general, step_lm, &
   & lm_predict_interval, print_lm_prediction_interval, lm_confint, lm_cooks_distance, print_lm_cooks_top, &
   & lm_coef, print_lm_summary, print_lm_coef_rstyle, print_lm_confint, print_lm_anova, pchisq, normal_cdf, qnorm, &
   & ppois, qpois, &
   & dunif, punif, qunif, dexp, pexp, qexp, dgamma, pgamma, qgamma, dbeta, pbeta, qbeta, dchisq, qchisq, &
   & dt, pt, qt, df, pf, qf, dlogis, plogis, qlogis, dlnorm, plnorm, qlnorm, dweibull, pweibull, qweibull, &
   & dcauchy, pcauchy, qcauchy, dbinom, pbinom, qbinom, dpois, dgeom, pgeom, qgeom, dnbinom, pnbinom, qnbinom, &
   & dhyper, phyper, qhyper, dwilcox, pwilcox, qwilcox, dsignrank, psignrank, qsignrank, cov, cor, cov2cor, &
   & mahalanobis, isSymmetric, scale, all_equal, r_log, r_seq_int, r_seq_len, &
   & glm_binomial_fit, glm_poisson_fit, glm_predict_response, glm_pearson_resid, print_glm_summary, &
   & prcomp, print_prcomp_summary, eigen, print_eigen, arima_fit_t, arima_predict_result_t, arima_sim, arima_fit, &
   & arima_predict, arima_predict_result, print_arima_fit, &
   & acf_fit_t, acf, r_acf, r_acf_values, r_ccf, print_acf, ar_fit_t, ar_fit, ARMAacf, &
   & r_seq_int_by, r_seq_int_length, r_seq_real_by, r_seq_real_length, &
   & r_paste0_real, r_paste0_int, r_index_real, r_index_scalar_real, r_matrix_col, r_matrix_row, r_matrix_rows, &
   & r_matrix_row_filter, r_matrix_col_filter, &
   & r_rep_real, r_rep_char, r_rep_int, r_drop_index, r_drop_indices, r_matrix_index, matrix_elem, r_head, rev_int, &
   & rev_real, r_array_real, r_array_int, r_array_char, matrix, &
   & r_matmul, r_add, r_sub, r_mul, r_div, print_matrix, &
   & print_matrix_rstyle, print_matrix_rstyle_named, print_real_scalar, &
   & print_real_vector, print_complex_vector, print_integer_vector, print_char_vector, &
   & display, print_named_real_vector, print_named_real_row, print_table1, print_table2, print_summary, &
   & print_factor, r_factor_labels, r_as_real, r_ifelse_real, r_na_real, r_inf, r_is_nan, set_print_int_like, &
   & set_print_int_like_tol, set_recycle_warn, set_recycle_stop, set_seed_int, &
   & kmeans_result_t, kmeans, rbind, max_col, tabulate, table2, prop_table, ave, ave_group_key, aggregate, &
   & aggregate_result_t, print_aggregate_result, r_by, by_matrix_result_t, print_by_matrix_result, match, r_in, &
   & unique, duplicated, anyDuplicated, &
   & union, intersect, setdiff, setequal, combn, findInterval, cut, cut_n, outer, &
   & outer_plus, outer_minus, outer_divide, outer_power, &
   & cumsum, cumprod, cummax, diff, diag, toeplitz, chol, chol2inv, forwardsolve, backsolve, sort, sort_list, &
   & polyroot, decompose, ecdf_eval, &
   & nchar, char_join, char_ends_with, int_to_string, real_to_string_f, real_to_string_g, r_to_string_real, &
   & r_substr, r_substr_replace, getwd, tempfile, file_path, file_exists, file_create, file_remove, file_info_t, &
   & file_info, file_isdir, file_size, file_path_value, file_mode, file_mtime, file_ctime, file_atime, file_exe, &
   & file_extension, print_file_info, dir_exists, dir_create, list_files, scan_real, grep_value_char, &
   & r_command_args, strsplit_fixed, toupper, tolower, casefold, trimws, replace_first_fixed, replace_all_fixed, &
   & chartr, urldecode, nextn, kronecker, fft, ar_coef_names, lag_names, lower_tri, upper_tri, row_index_mat, &
   & col_index_mat, matrix_set_grow_real, is_na, which, which_first, which_last, which_arr_ind, replace, rle, &
   & inverse_rle, print_rle, r_typeof, r_character, order_real, rank_average, &
   & rank_first, det_real, kappa_real, eigen_sym_values, solve_real, qr_fit_t, qr, qr_Q, qr_R, qr_coef, qr_rank, &
   & qr_pivot, qr_fitted, qr_resid, qr_qty, qr_qy, print_qr, &
   & rle_real_t, rle_int_t, rle_char_t, rle_logical_t, &
   & nested_matrix_list_len, r_beta, r_lbeta, r_choose, r_lchoose, r_gamma, r_lgamma, r_psigamma, r_digamma, &
   & r_trigamma, &
   & r_factorial, r_lfactorial, r_digit_power_sum, t_test_result_t, t_test, t_test_p_value, print_t_test, &
   & chisq_test_result_t, chisq_test, print_chisq_test, prop_test_result_t, &
   & prop_test, print_prop_test, cor_test_result_t, cor_test, print_cor_test, &
   & fisher_test_result_t, fisher_test, print_fisher_test, wilcox_test_result_t, &
   & wilcox_test, print_wilcox_test, kruskal_test_result_t, kruskal_test, ks_test, print_ks_test, print_factanal, &
   & r_filter_linear, smooth_xy_t, loess_fit_t, smooth_spline_fit_t, smooth, &
   & runmed, ksmooth, lowess, loess_fit, predict_loess, smooth_spline, &
   & predict_smooth_spline, dist, hclust_result_t, hclust, cutree, &
   & print_kruskal_test
public :: date_from_iso, date_from_iso_vec, date_from_yyyymmdd_vec, date_to_char, date_to_char_vec, &
   & date_format, date_format_vec, print_date, print_date_vector, r_elapsed, &
   & date_seq_day, date_seq_length, date_range, sys_time, sys_date, sys_date_string, &
   & sys_timezone, sys_time_format, sys_sleep, proc_time_vec, &
   & sys_getenv, file_rename, unlink_recursive
public :: as_octmode, as_hexmode, as_roman, inttobits, str_to_int, str_to_real, fivenum
public :: r_dataframe_t, table_char_t, table_char, data_frame_real, dataframe_real_col, print_dataframe, &
& print_dataframe_head
public :: r_tibble_real_t, r_tibble_integer_t, tibble_real, tibble_integer, &
   & tibble_nrow, tibble_ncol, tibble_real_col, tibble_integer_col, &
   & tibble_real_filter, tibble_real_select, tibble_real_drop, tibble_real_mutate, &
   & tibble_integer_filter, tibble_integer_select, tibble_integer_drop, &
   & tibble_integer_mutate, tibble_integer_to_real, &
   & tibble_real_log_returns, tibble_real_stats, read_csv_tibble_real, &
   & read_csv_tibble_integer, print_tibble
integer, parameter :: dp = real64
logical :: print_int_like_default = .true.
real(kind=dp) :: print_int_like_tol = 1000.0_dp * epsilon(1.0_dp)
logical :: recycle_warn_default = .false.
logical :: recycle_stop_default = .false.
integer :: nlm_method_default = 2
#ifdef XR2F_USE_R_RNG
interface
   subroutine xr2f_r_set_seed(seed) bind(C, name="xr2f_r_set_seed")
! Runtime helper for R-compatible xr2f r set seed.
      import :: c_int
      integer(c_int), value :: seed
   end subroutine xr2f_r_set_seed

   function xr2f_r_unif_rand() bind(C, name="xr2f_r_unif_rand") result(x)
! Runtime helper for R-compatible xr2f r unif rand.
      import :: c_double
      real(c_double) :: x
   end function xr2f_r_unif_rand

   function xr2f_r_norm_rand() bind(C, name="xr2f_r_norm_rand") result(x)
! Runtime helper for R-compatible xr2f r norm rand.
      import :: c_double
      real(c_double) :: x
   end function xr2f_r_norm_rand
end interface
#endif
type :: lm_fit_t
! Container for fitted lm fit model state.
   real(kind=dp), allocatable :: coef(:), fitted(:), resid(:), cov_unscaled(:,:), y(:), xpred(:,:)
   real(kind=dp) :: sigma, r_squared, adj_r_squared
   integer :: df = 0
   logical :: has_intercept = .true.
end type lm_fit_t

type :: r_dataframe_t
   character(len=:), allocatable :: names(:)
   real(kind=dp), allocatable :: real_cols(:,:)
end type r_dataframe_t

type :: r_tibble_real_t
! Real-only columnar table used by the initial restricted tibble API.
   character(len=:), allocatable :: names(:)
   character(len=:), allocatable :: row_labels(:)
   character(len=:), allocatable :: row_label_name
   real(kind=dp), allocatable :: real_cols(:,:)
end type r_tibble_real_t

type :: r_tibble_integer_t
! Integer-only columnar table used by the restricted tibble API.
   character(len=:), allocatable :: names(:)
   character(len=:), allocatable :: row_labels(:)
   character(len=:), allocatable :: row_label_name
   integer, allocatable :: integer_cols(:,:)
end type r_tibble_integer_t

interface tibble_nrow
   module procedure tibble_real_nrow
   module procedure tibble_integer_nrow
end interface tibble_nrow

interface tibble_ncol
   module procedure tibble_real_ncol
   module procedure tibble_integer_ncol
end interface tibble_ncol

interface tibble_real_mutate
   module procedure tibble_real_mutate_vector
   module procedure tibble_real_mutate_scalar
end interface tibble_real_mutate

interface tibble_integer_mutate
   module procedure tibble_integer_mutate_vector
   module procedure tibble_integer_mutate_scalar
   module procedure tibble_integer_mutate_real_vector
   module procedure tibble_integer_mutate_real_scalar
end interface tibble_integer_mutate

interface print_tibble
   module procedure print_tibble_real
   module procedure print_tibble_integer
end interface print_tibble

type :: table_char_t
   character(len=:), allocatable :: gene(:)
   integer, allocatable :: Freq(:)
end type table_char_t

type :: aggregate_result_t
   character(len=:), allocatable :: group_name, value_name
   character(len=:), allocatable :: labels(:)
   real(kind=dp), allocatable :: values(:)
end type aggregate_result_t

type :: by_matrix_result_t
   character(len=:), allocatable :: labels(:)
   real(kind=dp), allocatable :: values(:,:)
end type by_matrix_result_t

type :: glm_fit_t
! Container for fitted glm fit model state.
   real(kind=dp), allocatable :: coef(:), se(:), z_value(:), p_value(:), fitted(:), resid(:), y(:), xpred(:,:)
   real(kind=dp), allocatable :: offset(:)
   integer :: df = 0
   integer :: convergence = 1
   integer :: iter = 0
   integer :: family = 1
end type glm_fit_t

type :: prcomp_fit_t
! Container for fitted prcomp fit model state.
   real(kind=dp), allocatable :: sdev(:), rotation(:,:), x(:,:), center(:), scale(:)
end type prcomp_fit_t

type :: hist_result_t
! Container for R-like hist result values.
   real(kind=dp), allocatable :: breaks(:), mids(:), density(:)
   integer, allocatable :: counts(:)
end type hist_result_t

type :: file_info_t
! Compact scalar subset of R file.info() result fields.
   real(kind=dp) :: size = 0.0_dp
   logical :: isdir = .false.
   character(len=:), allocatable :: path, mode, mtime, ctime, atime, exe
end type file_info_t

interface hist
   module procedure hist_nbreaks
   module procedure hist_breaks_real
end interface hist

interface r_seq_int
   module procedure r_seq_int_ii, r_seq_int_ir, r_seq_int_ri, r_seq_int_rr
end interface r_seq_int

interface file_info
   module procedure file_info_scalar
   module procedure file_info_vector
end interface file_info

interface print_file_info
   module procedure print_file_info_scalar
   module procedure print_file_info_vector
end interface print_file_info

interface real_to_string_g
   module procedure real_to_string_g_scalar
   module procedure real_to_string_g_vector
end interface real_to_string_g

interface dir_exists
   module procedure dir_exists_scalar
   module procedure dir_exists_vector
end interface dir_exists

type :: eigen_result_t
! Container for R-like eigen result results.
   complex(kind=dp), allocatable :: values(:), vectors(:,:)
end type eigen_result_t

type :: arima_fit_t
! Container for fitted arima fit model state.
   real(kind=dp), allocatable :: coef(:), resid(:)
   real(kind=dp) :: mean = 0.0_dp
   real(kind=dp) :: aic = 0.0_dp
   real(kind=dp) :: sigma2 = 0.0_dp
   real(kind=dp) :: last_x = 0.0_dp
   real(kind=dp) :: last_resid = 0.0_dp
   integer :: p = 0
   integer :: d = 0
   integer :: q = 0
end type arima_fit_t

type :: arima_predict_result_t
! Container for R-like arima predict result results.
   real(kind=dp), allocatable :: pred(:), se(:)
end type arima_predict_result_t

type :: ar_fit_t
! Container for fitted ar fit model state.
   integer :: order = 0
   real(kind=dp), allocatable :: ar(:), aic(:)
   real(kind=dp) :: var_pred = 0.0_dp
end type ar_fit_t

type :: smooth_xy_t
! Container for R-like smooth xy data.
   real(kind=dp), allocatable :: x(:), y(:)
   real(kind=dp) :: df = 0.0_dp
end type smooth_xy_t

type :: loess_fit_t
! Container for fitted loess fit model state.
   real(kind=dp), allocatable :: x(:), y(:)
   real(kind=dp) :: span = 0.75_dp
   integer :: degree = 2
end type loess_fit_t

type :: smooth_spline_fit_t
! Container for fitted smooth spline fit model state.
   real(kind=dp), allocatable :: x(:), y(:)
   real(kind=dp) :: df = 0.0_dp
end type smooth_spline_fit_t

type :: acf_fit_t
! Container for fitted acf fit model state.
   real(kind=dp), allocatable :: acf(:,:,:)
   real(kind=dp), allocatable :: lag(:)
   integer :: n_used = 0
   integer :: type_code = 1
end type acf_fit_t

type :: qr_fit_t
! Container for fitted qr fit model state.
   real(kind=dp), allocatable :: qr(:,:), q(:,:), r(:,:), qraux(:)
   integer, allocatable :: pivot(:)
   integer :: rank = 0
end type qr_fit_t

type :: rle_real_t
! Container for run-length encoded real values.
   integer, allocatable :: lengths(:)
   real(kind=dp), allocatable :: values(:)
end type rle_real_t

type :: rle_int_t
! Container for run-length encoded int values.
   integer, allocatable :: lengths(:)
   integer, allocatable :: values(:)
end type rle_int_t

type :: rle_char_t
! Container for run-length encoded char values.
   integer, allocatable :: lengths(:)
   character(len=:), allocatable :: values(:)
end type rle_char_t

type :: rle_logical_t
! Container for run-length encoded logical values.
   integer, allocatable :: lengths(:)
   logical, allocatable :: values(:)
end type rle_logical_t

type :: kmeans_result_t
! Container for R-like kmeans result results.
   real(kind=dp), allocatable :: centers(:,:)
   integer, allocatable :: cluster(:)
   integer, allocatable :: size(:)
   real(kind=dp), allocatable :: withinss(:)
   real(kind=dp) :: totss = 0.0_dp
   real(kind=dp) :: tot_withinss = 0.0_dp
   real(kind=dp) :: betweenss = 0.0_dp
   integer :: iter = 0
end type kmeans_result_t

type :: hclust_result_t
! Container for R-like hclust result results.
   integer, allocatable :: merge(:,:)
   real(kind=dp), allocatable :: height(:)
   integer, allocatable :: order(:)
   character(len=:), allocatable :: labels(:)
   integer :: method = 1
end type hclust_result_t

type :: optim_result_t
! Container for R-like optim result results.
   real(kind=dp), allocatable :: par(:)
   real(kind=dp) :: value
   integer :: convergence
   integer :: counts(2) = 0
   real(kind=dp), allocatable :: hessian(:,:)
   character(len=:), allocatable :: message
end type optim_result_t

abstract interface
   pure function optim_vec_objective(par) result(value)
      import :: dp
      real(kind=dp), intent(in) :: par(:)
      real(kind=dp) :: value
   end function optim_vec_objective
   pure function optim_vec_gradient(par) result(value)
      import :: dp
      real(kind=dp), intent(in) :: par(:)
      real(kind=dp), allocatable :: value(:)
   end function optim_vec_gradient
end interface

type :: nlm_result_t
! Container for R-like nlm result results.
   real(kind=dp) :: minimum = 0.0_dp
   real(kind=dp), allocatable :: estimate(:), gradient(:), hessian(:,:)
   integer :: code = 1
   integer :: iterations = 0
end type nlm_result_t

type :: integrate_result_t
! Container for R-like integrate result values.
   real(kind=dp) :: value = 0.0_dp
   real(kind=dp) :: abs_error = 0.0_dp
   integer :: subdivisions = 0
   integer :: message = 0
end type integrate_result_t

abstract interface

   function nlm_objective_scalar(x) result(v)
! Support nlm-style optimization for objective scalar.
      import :: dp
      real(kind=dp), intent(in) :: x
      real(kind=dp) :: v
   end function nlm_objective_scalar

   function nlm_objective_vec(p) result(v)
! Support nlm-style optimization for objective vec.
      import :: dp
      real(kind=dp), intent(in) :: p(:)
      real(kind=dp) :: v
   end function nlm_objective_vec

   function integrate_objective(x) result(v)
! Support integrate-style scalar objective functions.
      import :: dp
      real(kind=dp), intent(in) :: x
      real(kind=dp) :: v
   end function integrate_objective

end interface

type :: decompose_result_t
! Container for R-like decompose result results.
   real(kind=dp), allocatable :: trend(:), seasonal(:), random(:), figure(:)
end type decompose_result_t

type :: t_test_result_t
! Container for R-like t test result results.
   real(kind=dp) :: statistic = 0.0_dp
   real(kind=dp) :: parameter = 0.0_dp
   real(kind=dp) :: p_value = 1.0_dp
   real(kind=dp) :: estimate = 0.0_dp
   real(kind=dp) :: estimate2 = 0.0_dp
   real(kind=dp) :: null_value = 0.0_dp
   real(kind=dp) :: conf_low = 0.0_dp
   real(kind=dp) :: conf_high = 0.0_dp
   real(kind=dp) :: stderr = 0.0_dp
   integer :: method = 1
end type t_test_result_t

type :: chisq_test_result_t
! Container for R-like chisq test result results.
   real(kind=dp) :: statistic = 0.0_dp
   integer :: parameter = 0
   real(kind=dp) :: p_value = 1.0_dp
   integer :: method = 1
end type chisq_test_result_t

type :: prop_test_result_t
! Container for R-like prop test result results.
   real(kind=dp) :: statistic = 0.0_dp
   integer :: parameter = 0
   real(kind=dp) :: p_value = 1.0_dp
   real(kind=dp) :: estimate = 0.0_dp
   real(kind=dp) :: estimate2 = 0.0_dp
   real(kind=dp) :: null_value = 0.0_dp
   integer :: method = 1
end type prop_test_result_t

type :: cor_test_result_t
! Container for R-like cor test result results.
   real(kind=dp) :: statistic = 0.0_dp
   integer :: parameter = 0
   real(kind=dp) :: p_value = 1.0_dp
   real(kind=dp) :: estimate = 0.0_dp
   integer :: method = 1
end type cor_test_result_t

type :: fisher_test_result_t
! Container for R-like fisher test result results.
   real(kind=dp) :: p_value = 1.0_dp
   real(kind=dp) :: estimate = 0.0_dp
   integer :: method = 1
end type fisher_test_result_t

type :: wilcox_test_result_t
! Container for R-like wilcox test result results.
   real(kind=dp) :: statistic = 0.0_dp
   real(kind=dp) :: p_value = 1.0_dp
   integer :: method = 1
end type wilcox_test_result_t

type :: kruskal_test_result_t
! Container for R-like kruskal test result results.
   real(kind=dp) :: statistic = 0.0_dp
   integer :: parameter = 0
   real(kind=dp) :: p_value = 1.0_dp
end type kruskal_test_result_t

type :: ks_test_result_t
! Container for R-like ks test result results.
   real(kind=dp) :: statistic = 0.0_dp
   real(kind=dp) :: p_value = 1.0_dp
   integer :: n = 0
end type ks_test_result_t

interface cov
   module procedure cov_vec
   module procedure cov_mat
end interface cov

interface var
   module procedure var_vec
   module procedure var_mat
end interface var

interface sd
   module procedure sd_vec
   module procedure sd_mat
end interface sd

interface cor
   module procedure cor_vec
   module procedure cor_mat
   module procedure cor_mat_pair
end interface cor

interface isSymmetric
   module procedure isSymmetric_real
   module procedure isSymmetric_int
end interface isSymmetric

interface scale
   module procedure scale_vec
   module procedure scale_mat
end interface scale

interface all_equal
   module procedure all_equal_real_scalar
   module procedure all_equal_real_vec
   module procedure all_equal_real_mat
   module procedure all_equal_int_scalar
   module procedure all_equal_int_vec
   module procedure all_equal_int_mat
   module procedure all_equal_complex_scalar
   module procedure all_equal_complex_vec
   module procedure all_equal_complex_mat
   module procedure all_equal_logical_scalar
   module procedure all_equal_logical_vec
   module procedure all_equal_logical_mat
   module procedure all_equal_char_scalar
   module procedure all_equal_char_vec
   module procedure all_equal_char_mat
end interface all_equal

interface r_log
   module procedure r_log_scalar
   module procedure r_log_vec
   module procedure r_log_mat
end interface r_log

interface det_real
   module procedure det_real_mat
   module procedure det_real_int
end interface det_real

interface kmeans
   module procedure kmeans_vec
   module procedure kmeans_mat
end interface kmeans

interface rbind
   module procedure rbind_vec
   module procedure rbind_mat
   module procedure rbind_vec_mat
   module procedure rbind_mat_vec
   module procedure rbind_int_vec
   module procedure rbind_int_mat
   module procedure rbind_int_vec_mat
   module procedure rbind_int_mat_vec
end interface rbind

interface char_join
   module procedure char_join_char
   module procedure char_join_int
end interface char_join

interface t_test
   module procedure t_test_one
   module procedure t_test_two
end interface t_test

interface t_test_p_value
   module procedure t_test_p_value_one
   module procedure t_test_p_value_two
end interface t_test_p_value

interface chisq_test
   module procedure chisq_test_int_vec
   module procedure chisq_test_real_vec
   module procedure chisq_test_int_mat
   module procedure chisq_test_real_mat
end interface chisq_test

interface prop_test
   module procedure prop_test_int_scalar
   module procedure prop_test_real_scalar
   module procedure prop_test_int_vec
   module procedure prop_test_real_vec
end interface prop_test

interface r_head
   module procedure r_head_real
   module procedure r_head_int
   module procedure r_head_real_mat
   module procedure r_head_int_mat
end interface r_head

interface cor_test
   module procedure cor_test_real_vec
   module procedure cor_test_int_vec
end interface cor_test

interface fisher_test
   module procedure fisher_test_int_mat
   module procedure fisher_test_real_mat
end interface fisher_test

interface wilcox_test
   module procedure wilcox_test_two_sample
end interface wilcox_test

interface rbinom
   module procedure rbinom_scalar
   module procedure rbinom_vector_prob
end interface rbinom

interface rpois
   module procedure rpois_scalar
   module procedure rpois_vector
end interface rpois

interface rexp
   module procedure rexp_rng
end interface rexp

interface rgamma
   module procedure rgamma_rng
end interface rgamma

interface rbeta
   module procedure rbeta_rng
end interface rbeta

interface rchisq
   module procedure rchisq_rng
end interface rchisq

interface rf_rng
   module procedure rf_rng_vec
end interface rf_rng

interface rlogis
   module procedure rlogis_rng
end interface rlogis

interface rlnorm
   module procedure rlnorm_rng
end interface rlnorm

interface rweibull
   module procedure rweibull_rng
end interface rweibull

interface rcauchy
   module procedure rcauchy_rng
end interface rcauchy

interface rgeom
   module procedure rgeom_rng
end interface rgeom

interface rnbinom
   module procedure rnbinom_rng
end interface rnbinom

interface rhyper
   module procedure rhyper_rng
end interface rhyper

interface rwilcox
   module procedure rwilcox_rng
end interface rwilcox

interface rsignrank
   module procedure rsignrank_rng
end interface rsignrank

interface rmultinom
   module procedure rmultinom_rng
end interface rmultinom

interface dnorm
   module procedure dnorm_vec
   module procedure dnorm_scalar
   module procedure dnorm_int_scalar
end interface dnorm

interface pchisq
   module procedure pchisq_scalar
   module procedure pchisq_vec
end interface pchisq

interface qchisq
   module procedure qchisq_scalar
   module procedure qchisq_scalar_i
   module procedure qchisq_vec
end interface qchisq

interface qnorm
   module procedure qnorm_scalar
   module procedure qnorm_vec
end interface qnorm

interface ppois
   module procedure ppois_scalar
   module procedure ppois_vec
end interface ppois

interface qpois
   module procedure qpois_scalar
   module procedure qpois_vec
end interface qpois

interface dunif
   module procedure dunif_scalar
   module procedure dunif_vec
end interface dunif

interface punif
   module procedure punif_scalar
   module procedure punif_vec
end interface punif

interface qunif
   module procedure qunif_scalar
   module procedure qunif_vec
end interface qunif

interface dexp
   module procedure dexp_scalar
   module procedure dexp_vec
end interface dexp

interface pexp
   module procedure pexp_scalar
   module procedure pexp_vec
end interface pexp

interface qexp
   module procedure qexp_scalar
   module procedure qexp_vec
end interface qexp

interface dgamma
   module procedure dgamma_scalar
   module procedure dgamma_vec
end interface dgamma

interface pgamma
   module procedure pgamma_scalar
   module procedure pgamma_vec
end interface pgamma

interface qgamma
   module procedure qgamma_scalar
   module procedure qgamma_vec
end interface qgamma

interface dbeta
   module procedure dbeta_scalar
   module procedure dbeta_vec
end interface dbeta

interface pbeta
   module procedure pbeta_scalar
   module procedure pbeta_vec
end interface pbeta

interface qbeta
   module procedure qbeta_scalar
   module procedure qbeta_vec
end interface qbeta

interface dchisq
   module procedure dchisq_scalar
   module procedure dchisq_vec
end interface dchisq

interface dt
   module procedure dt_scalar
   module procedure dt_vec
end interface dt

interface pt
   module procedure pt_scalar
   module procedure pt_vec
end interface pt

interface qt
   module procedure qt_scalar
   module procedure qt_vec
end interface qt

interface df
   module procedure df_scalar
   module procedure df_vec
end interface df

interface pf
   module procedure pf_scalar
   module procedure pf_vec
end interface pf

interface qf
   module procedure qf_scalar
   module procedure qf_vec
end interface qf

interface dlogis
   module procedure dlogis_scalar
   module procedure dlogis_vec
end interface dlogis

interface plogis
   module procedure plogis_scalar
   module procedure plogis_vec
end interface plogis

interface qlogis
   module procedure qlogis_scalar
   module procedure qlogis_vec
end interface qlogis

interface dlnorm
   module procedure dlnorm_scalar
   module procedure dlnorm_vec
end interface dlnorm

interface plnorm
   module procedure plnorm_scalar
   module procedure plnorm_vec
end interface plnorm

interface qlnorm
   module procedure qlnorm_scalar
   module procedure qlnorm_vec
end interface qlnorm

interface dweibull
   module procedure dweibull_scalar
   module procedure dweibull_vec
end interface dweibull

interface pweibull
   module procedure pweibull_scalar
   module procedure pweibull_vec
end interface pweibull

interface qweibull
   module procedure qweibull_scalar
   module procedure qweibull_vec
end interface qweibull

interface dcauchy
   module procedure dcauchy_scalar
   module procedure dcauchy_vec
end interface dcauchy

interface pcauchy
   module procedure pcauchy_scalar
   module procedure pcauchy_vec
end interface pcauchy

interface qcauchy
   module procedure qcauchy_scalar
   module procedure qcauchy_vec
end interface qcauchy

interface dbinom
   module procedure dbinom_scalar
   module procedure dbinom_vec
end interface dbinom

interface pbinom
   module procedure pbinom_scalar
   module procedure pbinom_vec
end interface pbinom

interface qbinom
   module procedure qbinom_scalar
   module procedure qbinom_vec
end interface qbinom

interface dpois
   module procedure dpois_scalar
   module procedure dpois_vec
end interface dpois

interface dgeom
   module procedure dgeom_scalar
   module procedure dgeom_vec
end interface dgeom

interface pgeom
   module procedure pgeom_scalar
   module procedure pgeom_vec
end interface pgeom

interface qgeom
   module procedure qgeom_scalar
   module procedure qgeom_vec
end interface qgeom

interface dnbinom
   module procedure dnbinom_scalar
   module procedure dnbinom_vec
end interface dnbinom

interface pnbinom
   module procedure pnbinom_scalar
   module procedure pnbinom_vec
end interface pnbinom

interface qnbinom
   module procedure qnbinom_scalar
   module procedure qnbinom_vec
end interface qnbinom

interface dhyper
   module procedure dhyper_scalar
   module procedure dhyper_vec
end interface dhyper

interface phyper
   module procedure phyper_scalar
   module procedure phyper_vec
end interface phyper

interface qhyper
   module procedure qhyper_scalar
   module procedure qhyper_vec
end interface qhyper

interface dwilcox
   module procedure dwilcox_scalar
   module procedure dwilcox_vec
end interface dwilcox

interface pwilcox
   module procedure pwilcox_scalar
   module procedure pwilcox_vec
end interface pwilcox

interface qwilcox
   module procedure qwilcox_scalar
   module procedure qwilcox_vec
end interface qwilcox

interface dsignrank
   module procedure dsignrank_scalar
   module procedure dsignrank_vec
end interface dsignrank

interface psignrank
   module procedure psignrank_scalar
   module procedure psignrank_vec
end interface psignrank

interface qsignrank
   module procedure qsignrank_scalar
   module procedure qsignrank_vec
end interface qsignrank


interface glm_binomial_fit
   module procedure glm_binomial_fit_real
   module procedure glm_binomial_fit_int
end interface glm_binomial_fit

interface glm_poisson_fit
   module procedure glm_poisson_fit_real
   module procedure glm_poisson_fit_int
end interface glm_poisson_fit

interface matrix
   module procedure matrix_real
   module procedure matrix_int
end interface matrix

interface summary
   module procedure summary_vec
   module procedure summary_mat
end interface summary

interface print_summary
   module procedure print_summary_vec
   module procedure print_summary_mat
end interface print_summary

interface r_matmul
   module procedure r_matmul_vv_real
   module procedure r_matmul_vv_int
   module procedure r_matmul_vv_real_int
   module procedure r_matmul_vv_int_real
   module procedure r_matmul_mv_real
   module procedure r_matmul_mv_int
   module procedure r_matmul_mv_real_int
   module procedure r_matmul_mv_int_real
   module procedure r_matmul_vm_real
   module procedure r_matmul_vm_int
   module procedure r_matmul_vm_real_int
   module procedure r_matmul_vm_int_real
   module procedure r_matmul_mm_real
   module procedure r_matmul_mm_int
   module procedure r_matmul_mm_real_int
   module procedure r_matmul_mm_int_real
   module procedure r_matmul_mv_complex
   module procedure r_matmul_mv_real_complex
   module procedure r_matmul_mv_int_complex
   module procedure r_matmul_mm_complex
   module procedure r_matmul_mm_complex_real
   module procedure r_matmul_mm_real_complex
end interface r_matmul

interface r_add
   module procedure r_add_vv
   module procedure r_add_vs
   module procedure r_add_sv
   module procedure r_add_mv
   module procedure r_add_vm
end interface r_add

interface r_sub
   module procedure r_sub_vv
   module procedure r_sub_vs
   module procedure r_sub_sv
   module procedure r_sub_mv
   module procedure r_sub_vm
end interface r_sub

interface r_mul
   module procedure r_mul_vv
   module procedure r_mul_vs
   module procedure r_mul_sv
   module procedure r_mul_mv
   module procedure r_mul_vm
end interface r_mul

interface r_div
   module procedure r_div_vv
   module procedure r_div_vs
   module procedure r_div_sv
   module procedure r_div_mv
   module procedure r_div_vm
end interface r_div

interface r_array
   module procedure r_array_real
   module procedure r_array_int
   module procedure r_array_char
end interface r_array

interface tabulate
   module procedure tabulate_int
   module procedure tabulate_real
end interface tabulate

interface table2
   module procedure table2_int
end interface table2

interface prop_table
   module procedure prop_table_int_vec
   module procedure prop_table_int_mat
end interface prop_table

interface ave
   module procedure ave_real_char
   module procedure ave_real_int
end interface ave

interface ave_group_key
   module procedure ave_group_key_char_char
   module procedure ave_group_key_char_int
   module procedure ave_group_key_int_char
   module procedure ave_group_key_int_int
end interface ave_group_key

interface aggregate
   module procedure aggregate_real_char
   module procedure aggregate_real_int
end interface aggregate

interface r_by
   module procedure by_real_vec_char
   module procedure by_real_vec_int
   module procedure by_real_mat_char
   module procedure by_real_mat_int
end interface r_by

interface print_table2
   module procedure print_table2_int
   module procedure print_table2_real
end interface print_table2

interface match
   module procedure match_int
   module procedure match_real
   module procedure match_char
end interface match

interface r_in
   module procedure r_in_int
   module procedure r_in_real
   module procedure r_in_int_real
   module procedure r_in_real_int
   module procedure r_in_int_scalar
   module procedure r_in_real_scalar
   module procedure r_in_int_scalar_real
   module procedure r_in_real_scalar_int
   module procedure r_in_char
   module procedure r_in_logical
end interface r_in

interface r_drop_index
   module procedure r_drop_index_real
   module procedure r_drop_index_int
end interface r_drop_index

interface r_drop_indices
   module procedure r_drop_indices_real
   module procedure r_drop_indices_int
end interface r_drop_indices

interface r_matrix_index
   module procedure r_vector_index_real_int
   module procedure r_vector_index_int_int
   module procedure r_vector_index_real_logical
   module procedure r_vector_index_int_logical
   module procedure r_matrix_index_real
   module procedure r_matrix_index_int
   module procedure r_matrix_index_real_logical_vec
   module procedure r_matrix_index_int_logical_vec
   module procedure r_matrix_index_real_logical
   module procedure r_matrix_index_int_logical
end interface r_matrix_index

interface r_matrix_col
   module procedure r_matrix_col_real
   module procedure r_matrix_col_int
end interface r_matrix_col

interface r_matrix_row
   module procedure r_matrix_row_real
   module procedure r_matrix_row_int
end interface r_matrix_row

interface r_matrix_rows
   module procedure r_matrix_rows_real
   module procedure r_matrix_rows_int
end interface r_matrix_rows

interface r_matrix_row_filter
   module procedure r_matrix_row_filter_real
   module procedure r_matrix_row_filter_int
end interface r_matrix_row_filter

interface r_matrix_col_filter
   module procedure r_matrix_col_filter_real
   module procedure r_matrix_col_filter_int
end interface r_matrix_col_filter

interface unique
   module procedure unique_int
   module procedure unique_real
   module procedure unique_char
   module procedure unique_logical
end interface unique

interface duplicated
   module procedure duplicated_int
   module procedure duplicated_real
   module procedure duplicated_char
   module procedure duplicated_logical
end interface duplicated

interface anyDuplicated
   module procedure anyDuplicated_int
   module procedure anyDuplicated_real
   module procedure anyDuplicated_char
   module procedure anyDuplicated_logical
end interface anyDuplicated

interface union
   module procedure union_int
   module procedure union_real
   module procedure union_char
end interface union

interface intersect
   module procedure intersect_int
   module procedure intersect_real
   module procedure intersect_char
end interface intersect

interface setdiff
   module procedure setdiff_int
   module procedure setdiff_real
   module procedure setdiff_char
end interface setdiff

interface setequal
   module procedure setequal_int
   module procedure setequal_real
   module procedure setequal_char
end interface setequal

interface combn
   module procedure combn_int
   module procedure combn_real
   module procedure combn_char
end interface combn

interface print_matrix
   module procedure print_matrix_real
   module procedure print_matrix_int
   module procedure print_matrix_logical
   module procedure print_matrix_complex
end interface print_matrix

interface print_matrix_rstyle
   module procedure print_matrix_rstyle_real
   module procedure print_matrix_rstyle_int
   module procedure print_matrix_rstyle_logical
   module procedure print_matrix_complex
end interface print_matrix_rstyle

interface print_matrix_rstyle_named
   module procedure print_matrix_rstyle_named_real
   module procedure print_matrix_rstyle_named_int
end interface print_matrix_rstyle_named

interface display
   module procedure print_real_scalar
   module procedure display_integer_scalar, display_logical_scalar, display_complex_scalar, display_char_scalar
   module procedure print_real_vector, print_integer_vector, display_logical_vector, print_complex_vector, &
   & print_char_vector
   module procedure print_matrix_rstyle_real, print_matrix_rstyle_int, print_matrix_rstyle_logical
   module procedure print_matrix_complex, display_char_matrix
   module procedure display_real_array3, display_integer_array3, display_logical_array3
   module procedure display_complex_array3, display_char_array3
end interface display

interface r_as_real
   module procedure r_as_real_char
end interface r_as_real

interface cumsum
   module procedure cumsum_real
   module procedure cumsum_int
end interface cumsum

interface cumprod
   module procedure cumprod_real
   module procedure cumprod_int
end interface cumprod

interface cummax
   module procedure cummax_real
   module procedure cummax_int
end interface cummax

interface diff
   module procedure diff_real
   module procedure diff_mat_real
   module procedure diff_int
end interface diff

interface diag
   module procedure diag_mat_real
   module procedure diag_vec_real
   module procedure diag_vec_real_n
   module procedure diag_mat_complex
   module procedure diag_vec_complex
   module procedure diag_vec_complex_n
   module procedure diag_mat_int
   module procedure diag_vec_int
   module procedure diag_vec_int_n
   module procedure diag_scalar_int
   module procedure diag_scalar_real_n
end interface diag

interface eigen
   module procedure eigen_real
   module procedure eigen_int
end interface eigen

interface sort
   module procedure sort_real
   module procedure sort_int
   module procedure sort_char
end interface sort

interface sort_list
   module procedure sort_list_real
   module procedure sort_list_int
   module procedure sort_list_char
end interface sort_list

interface solve_real
   module procedure solve_real_vec
   module procedure solve_real_vec_i_r
   module procedure solve_real_vec_i_i
   module procedure solve_real_mat
   module procedure solve_real_mat_r_i
   module procedure solve_real_mat_i_r
   module procedure solve_real_mat_i_i
   module procedure solve_real_vec_r_c
   module procedure solve_real_vec_i_c
   module procedure solve_complex_vec
   module procedure solve_complex_mat
end interface solve_real

interface besselJ
   module procedure besselJ_scalar_i
   module procedure besselJ_scalar_r
   module procedure besselJ_vec_i
   module procedure besselJ_vec_r
end interface besselJ

interface besselY
   module procedure besselY_scalar_i
   module procedure besselY_scalar_r
   module procedure besselY_vec_i
   module procedure besselY_vec_r
end interface besselY

interface besselI
   module procedure besselI_scalar_i
   module procedure besselI_scalar_r
   module procedure besselI_vec_i
   module procedure besselI_vec_r
end interface besselI

interface besselK
   module procedure besselK_scalar_i
   module procedure besselK_scalar_r
   module procedure besselK_vec_i
   module procedure besselK_vec_r
end interface besselK

interface chol
   module procedure chol_real
   module procedure chol_int
end interface chol

interface chol2inv
   module procedure chol2inv_real
   module procedure chol2inv_int
end interface chol2inv

interface forwardsolve
   module procedure forwardsolve_vec
   module procedure forwardsolve_vec_i_r
   module procedure forwardsolve_vec_i_i
   module procedure forwardsolve_mat
   module procedure forwardsolve_mat_r_i
   module procedure forwardsolve_mat_i_r
   module procedure forwardsolve_mat_i_i
end interface forwardsolve

interface backsolve
   module procedure backsolve_vec
   module procedure backsolve_vec_i_r
   module procedure backsolve_vec_i_i
   module procedure backsolve_mat
   module procedure backsolve_mat_r_i
   module procedure backsolve_mat_i_r
   module procedure backsolve_mat_i_i
end interface backsolve

interface qr_coef
   module procedure qr_coef_vec
   module procedure qr_coef_mat
end interface qr_coef

interface qr_fitted
   module procedure qr_fitted_vec
   module procedure qr_fitted_mat
end interface qr_fitted

interface qr_resid
   module procedure qr_resid_vec
   module procedure qr_resid_mat
end interface qr_resid

interface qr_qty
   module procedure qr_qty_vec
   module procedure qr_qty_mat
end interface qr_qty

interface qr_qy
   module procedure qr_qy_vec
   module procedure qr_qy_mat
end interface qr_qy

interface hclust
   module procedure hclust_complete
end interface hclust

interface cutree
   module procedure cutree_f90
end interface cutree

interface dist
   module procedure dist_mat
end interface dist

interface arima_sim
   module procedure arima_sim_scalar
   module procedure arima_sim_vector
   module procedure arima_sim_vector_vector
end interface arima_sim

interface acf
   module procedure acf_vec
   module procedure acf_mat
end interface acf

interface r_acf
   module procedure acf_vec
   module procedure acf_mat
end interface r_acf

interface r_acf_values
   module procedure acf_values_vec
   module procedure acf_values_mat
end interface r_acf_values

interface r_ccf
   module procedure ccf_vec
end interface r_ccf

interface ARMAacf
   module procedure ARMAacf_vecma
   module procedure ARMAacf_scalarma
end interface ARMAacf

interface is_na
   module procedure is_na_real_scalar
   module procedure is_na_real_vec
   module procedure is_na_real_mat
   module procedure is_na_int_scalar
   module procedure is_na_int_vec
   module procedure is_na_int_mat
   module procedure is_na_logical_scalar
   module procedure is_na_logical_vec
   module procedure is_na_logical_mat
   module procedure is_na_complex_scalar
   module procedure is_na_complex_vec
   module procedure is_na_complex_mat
   module procedure is_na_char_scalar
   module procedure is_na_char_vec
   module procedure is_na_char_mat
end interface is_na

interface which
   module procedure which_logical
   module procedure which_logical_mat
end interface which

interface replace
   module procedure replace_real_idx_scalar
   module procedure replace_real_idx_vec
   module procedure replace_real_mask_scalar
   module procedure replace_real_mask_vec
   module procedure replace_int_idx_scalar
   module procedure replace_int_idx_vec
   module procedure replace_int_idx_scalar_real
   module procedure replace_int_idx_vec_real
   module procedure replace_int_mask_scalar
   module procedure replace_int_mask_vec
   module procedure replace_int_mask_scalar_real
   module procedure replace_int_mask_vec_real
   module procedure replace_char_mask_scalar
   module procedure replace_char_mask_vec
   module procedure replace_logical_mask_scalar
end interface replace

interface rle
   module procedure rle_real
   module procedure rle_int
   module procedure rle_char
   module procedure rle_logical
end interface rle

interface inverse_rle
   module procedure inverse_rle_real
   module procedure inverse_rle_int
   module procedure inverse_rle_char
   module procedure inverse_rle_logical
end interface inverse_rle

interface print_rle
   module procedure print_rle_real
   module procedure print_rle_int
   module procedure print_rle_char
   module procedure print_rle_logical
end interface print_rle

interface r_typeof
   module procedure r_typeof_real_scalar
   module procedure r_typeof_real_vec
   module procedure r_typeof_real_mat
   module procedure r_typeof_complex_scalar
   module procedure r_typeof_complex_vec
   module procedure r_typeof_complex_mat
   module procedure r_typeof_int_scalar
   module procedure r_typeof_int_vec
   module procedure r_typeof_int_mat
   module procedure r_typeof_char_scalar
   module procedure r_typeof_char_vec
   module procedure r_typeof_char_mat
   module procedure r_typeof_logical_scalar
   module procedure r_typeof_logical_vec
   module procedure r_typeof_logical_mat
end interface r_typeof

contains

function optim_bfgs(fn, par, maxit, reltol, ndeps, ndeps_vec, fnscale, parscale, gr, hessian) result(out)
! Quasi-Newton optimizer for vector-valued parameter objectives.
procedure(optim_vec_objective) :: fn
procedure(optim_vec_gradient), optional :: gr
real(kind=dp), intent(in) :: par(:)
integer, intent(in), optional :: maxit
real(kind=dp), intent(in), optional :: reltol, ndeps, ndeps_vec(:), fnscale, parscale(:)
logical, intent(in), optional :: hessian
type(optim_result_t) :: out
integer :: n, max_iter, n_iter, i, j, iter, fn_count, gr_count, conv_code
logical :: converged, grad_failed
character(len=:), allocatable :: msg, grad_msg
real(kind=dp) :: f, f_new, step_eps, gtol, fscale
real(kind=dp) :: alpha, slope, sy, rho, shift
real(kind=dp), allocatable :: p(:), p_new(:), g(:), g_new(:), pscale(:), g_raw(:)
real(kind=dp), allocatable :: h(:,:), d(:), s(:), y(:), a(:,:), tmp(:,:)
n = size(par)
max_iter = 100
if (present(maxit)) max_iter = maxit
step_eps = 1.0e-3_dp
if (present(ndeps)) step_eps = ndeps
fscale = 1.0_dp
if (present(fnscale)) then
   if (abs(fnscale) > tiny(1.0_dp)) fscale = fnscale
end if
gtol = 1.0e-8_dp
if (present(reltol)) gtol = reltol
gtol = max(gtol, sqrt(epsilon(1.0_dp)))
allocate(p(n), p_new(n), g(n), g_new(n), pscale(n), h(n,n), d(n), s(n), y(n), a(n,n), tmp(n,n))
fn_count = 0
gr_count = 0
conv_code = 1
msg = "maximum iterations reached"
grad_failed = .false.
grad_msg = ""
p = par
pscale = 1.0_dp
if (present(parscale)) pscale(1:min(n, size(parscale))) = max(abs(parscale(1:min(n, size(parscale)))), tiny(1.0_dp))
h = 0.0_dp
do i = 1, n
   h(i,i) = 1.0_dp
end do
f = scaled_fn(p)
fn_count = fn_count + 1
call eval_gradient(p, g)
if (grad_failed) then
   out%par = p
   out%value = f * fscale
   out%convergence = 51
   out%counts = [fn_count, gr_count]
   out%message = grad_msg
   return
end if
if (.not. ieee_is_finite(f) .or. .not. finite_vec(g)) then
   out%par = p
   out%value = f * fscale
   out%convergence = 10
   out%counts = [fn_count, gr_count]
   out%message = "non-finite initial objective or gradient"
   return
end if
converged = .false.
n_iter = 0
do iter = 1, max_iter
   n_iter = iter
   if (sqrt(sum(g**2)) < gtol) then
      converged = .true.
      conv_code = 0
      msg = "converged"
      exit
   end if
   do i = 1, n
      d(i) = -sum(h(i,:) * g)
   end do
   if (dot_product(g, d) >= 0.0_dp) then
      h = 0.0_dp
      do i = 1, n
         h(i,i) = 1.0_dp
      end do
      d = -g
   end if
   alpha = 1.0_dp
   slope = dot_product(g, d)
   f_new = huge(1.0_dp)
   do j = 1, 60
      p_new = p + alpha * d
      f_new = scaled_fn(p_new)
      fn_count = fn_count + 1
      if (ieee_is_finite(f_new) .and. f_new <= f + 1.0e-4_dp * alpha * slope) exit
      if (alpha < 1.0e-12_dp) exit
      alpha = 0.5_dp * alpha
   end do
   if (.not. ieee_is_finite(f_new)) then
      conv_code = 51
      msg = "line search failed"
      exit
   end if
   call eval_gradient(p_new, g_new)
   if (grad_failed) then
      conv_code = 51
      msg = grad_msg
      exit
   end if
   if (.not. finite_vec(g_new)) then
      conv_code = 51
      msg = "non-finite gradient"
      exit
   end if
   s = p_new - p
   y = g_new - g
   sy = dot_product(s, y)
   shift = abs(f - f_new)
   p = p_new
   f = f_new
   g = g_new
   if (ieee_is_finite(sy) .and. finite_vec(s) .and. finite_vec(y) .and. &
      & sy > 1.0e-10_dp * sqrt(sum(s**2)) * sqrt(sum(y**2))) then
      rho = 1.0_dp / sy
      do i = 1, n
         do j = 1, n
            a(i,j) = -rho * s(i) * y(j)
         end do
         a(i,i) = a(i,i) + 1.0_dp
      end do
      do i = 1, n
         do j = 1, n
            tmp(i,j) = sum(h(i,:) * a(j,:))
         end do
      end do
      do i = 1, n
         do j = 1, n
            h(i,j) = sum(a(i,:) * tmp(:,j)) + rho * s(i) * s(j)
         end do
      end do
   end if
   if (shift <= gtol * (1.0_dp + abs(f))) then
      converged = .true.
      conv_code = 0
      msg = "converged"
      exit
   end if
end do
out%par = p
out%value = f * fscale
out%convergence = merge(0, conv_code, converged)
out%counts = [fn_count, gr_count]
out%message = msg
if (present(hessian)) then
   if (hessian) then
      call optim_fd_hessian(scaled_fn, p, step_eps, out%hessian, pscale, ndeps_vec)
      out%hessian = out%hessian * fscale
   end if
end if
contains
pure function scaled_fn(x) result(value)
real(kind=dp), intent(in) :: x(:)
real(kind=dp) :: value
value = fn(x) / fscale
end function scaled_fn

subroutine eval_gradient(x, gout)
real(kind=dp), intent(in) :: x(:)
real(kind=dp), intent(out) :: gout(:)
grad_failed = .false.
grad_msg = ""
if (present(gr)) then
   gr_count = gr_count + 1
   g_raw = gr(x)
   if (size(g_raw) /= size(gout)) then
      grad_failed = .true.
      grad_msg = "analytic gradient has wrong length"
      gout = ieee_value(0.0_dp, ieee_quiet_nan)
      return
   end if
   if (.not. finite_vec(g_raw)) then
      grad_failed = .true.
      grad_msg = "analytic gradient is non-finite"
      gout = ieee_value(0.0_dp, ieee_quiet_nan)
      return
   end if
   gout = g_raw / fscale
   return
end if
gr_count = gr_count + 1
fn_count = fn_count + 2 * size(x)
call optim_fd_gradient(scaled_fn, x, step_eps, gout, pscale, ndeps_vec)
end subroutine eval_gradient

pure function finite_vec(x) result(ok)
real(kind=dp), intent(in) :: x(:)
logical :: ok
ok = all(ieee_is_finite(x))
end function finite_vec
end function optim_bfgs

function optim_lbfgsb(fn, par, lower, upper, maxit, reltol, ndeps, ndeps_vec, fnscale, parscale, gr, &
& hessian) result(out)
! Bound-constrained optimizer for vector-valued parameter objectives.
procedure(optim_vec_objective) :: fn
procedure(optim_vec_gradient), optional :: gr
real(kind=dp), intent(in) :: par(:), lower(:), upper(:)
integer, intent(in), optional :: maxit
real(kind=dp), intent(in), optional :: reltol, ndeps, ndeps_vec(:), fnscale, parscale(:)
logical, intent(in), optional :: hessian
type(optim_result_t) :: out
integer :: n, max_iter, n_iter, i, j, iter, fn_count, gr_count, conv_code
logical :: converged, grad_failed
character(len=:), allocatable :: msg, grad_msg
real(kind=dp) :: f, f_new, step_eps, gtol, fscale
real(kind=dp) :: alpha, slope, sy, rho, shift
real(kind=dp), allocatable :: p(:), p_new(:), g(:), g_new(:), pscale(:), lb(:), ub(:), g_raw(:)
real(kind=dp), allocatable :: h(:,:), d(:), s(:), y(:), a(:,:), tmp(:,:)
n = size(par)
max_iter = 100
if (present(maxit)) max_iter = maxit
step_eps = 1.0e-3_dp
if (present(ndeps)) step_eps = ndeps
fscale = 1.0_dp
if (present(fnscale)) then
   if (abs(fnscale) > tiny(1.0_dp)) fscale = fnscale
end if
gtol = 1.0e-8_dp
if (present(reltol)) gtol = reltol
gtol = max(gtol, sqrt(epsilon(1.0_dp)))
allocate(p(n), p_new(n), g(n), g_new(n), pscale(n), lb(n), ub(n), h(n,n), d(n), s(n), y(n), a(n,n), tmp(n,n))
fn_count = 0
gr_count = 0
conv_code = 1
msg = "maximum iterations reached"
grad_failed = .false.
grad_msg = ""
do i = 1, n
   lb(i) = lower(min(i, size(lower)))
   ub(i) = upper(min(i, size(upper)))
   if (lb(i) > ub(i)) then
      out%par = par
      out%value = huge(1.0_dp)
      out%convergence = 51
      out%counts = [fn_count, gr_count]
      out%message = "invalid bounds"
      return
   end if
end do
p = project(par)
pscale = 1.0_dp
if (present(parscale)) pscale(1:min(n, size(parscale))) = max(abs(parscale(1:min(n, size(parscale)))), tiny(1.0_dp))
h = 0.0_dp
do i = 1, n
   h(i,i) = 1.0_dp
end do
f = scaled_fn(p)
fn_count = fn_count + 1
call eval_gradient(p, g)
if (grad_failed) then
   out%par = p
   out%value = f * fscale
   out%convergence = 51
   out%counts = [fn_count, gr_count]
   out%message = grad_msg
   return
end if
if (.not. ieee_is_finite(f) .or. .not. finite_vec(g)) then
   out%par = p
   out%value = f * fscale
   out%convergence = 10
   out%counts = [fn_count, gr_count]
   out%message = "non-finite initial objective or gradient"
   return
end if
converged = .false.
n_iter = 0
do iter = 1, max_iter
   n_iter = iter
   do i = 1, n
      if ((p(i) <= lb(i) + sqrt(epsilon(1.0_dp)) .and. g(i) > 0.0_dp) .or. &
         & (p(i) >= ub(i) - sqrt(epsilon(1.0_dp)) .and. g(i) < 0.0_dp)) g(i) = 0.0_dp
   end do
   if (sqrt(sum(g**2)) < gtol) then
      converged = .true.
      conv_code = 0
      msg = "converged"
      exit
   end if
   do i = 1, n
      d(i) = -sum(h(i,:) * g)
   end do
   if (dot_product(g, d) >= 0.0_dp) then
      h = 0.0_dp
      do i = 1, n
         h(i,i) = 1.0_dp
      end do
      d = -g
   end if
   alpha = 1.0_dp
   slope = dot_product(g, d)
   f_new = huge(1.0_dp)
   do j = 1, 60
      p_new = project(p + alpha * d)
      f_new = scaled_fn(p_new)
      fn_count = fn_count + 1
      if (ieee_is_finite(f_new) .and. f_new <= f + 1.0e-4_dp * alpha * slope) exit
      if (alpha < 1.0e-12_dp) exit
      alpha = 0.5_dp * alpha
   end do
   if (.not. ieee_is_finite(f_new)) then
      conv_code = 51
      msg = "line search failed"
      exit
   end if
   call eval_gradient(p_new, g_new)
   if (grad_failed) then
      conv_code = 51
      msg = grad_msg
      exit
   end if
   if (.not. finite_vec(g_new)) then
      conv_code = 51
      msg = "non-finite gradient"
      exit
   end if
   s = p_new - p
   y = g_new - g
   sy = dot_product(s, y)
   shift = abs(f - f_new)
   p = p_new
   f = f_new
   g = g_new
   if (ieee_is_finite(sy) .and. finite_vec(s) .and. finite_vec(y) .and. &
      & sy > 1.0e-10_dp * sqrt(sum(s**2)) * sqrt(sum(y**2))) then
      rho = 1.0_dp / sy
      do i = 1, n
         do j = 1, n
            a(i,j) = -rho * s(i) * y(j)
         end do
         a(i,i) = a(i,i) + 1.0_dp
      end do
      do i = 1, n
         do j = 1, n
            tmp(i,j) = sum(h(i,:) * a(j,:))
         end do
      end do
      do i = 1, n
         do j = 1, n
            h(i,j) = sum(a(i,:) * tmp(:,j)) + rho * s(i) * s(j)
         end do
      end do
   end if
   if (shift <= gtol * (1.0_dp + abs(f))) then
      h = 0.0_dp
      do i = 1, n
         h(i,i) = 1.0_dp
      end do
   end if
end do
out%par = p
out%value = f * fscale
out%convergence = merge(0, conv_code, converged)
out%counts = [fn_count, gr_count]
out%message = msg
if (present(hessian)) then
   if (hessian) then
      call optim_fd_hessian(scaled_fn, p, step_eps, out%hessian, pscale, ndeps_vec)
      out%hessian = out%hessian * fscale
   end if
end if
contains
pure function project(x) result(px)
real(kind=dp), intent(in) :: x(:)
real(kind=dp) :: px(size(x))
integer :: k
do k = 1, size(x)
   px(k) = min(max(x(k), lb(k)), ub(k))
end do
end function project

pure function scaled_fn(x) result(value)
real(kind=dp), intent(in) :: x(:)
real(kind=dp) :: value
value = fn(project(x)) / fscale
end function scaled_fn

subroutine eval_gradient(x, gout)
real(kind=dp), intent(in) :: x(:)
real(kind=dp), intent(out) :: gout(:)
grad_failed = .false.
grad_msg = ""
if (present(gr)) then
   gr_count = gr_count + 1
   g_raw = gr(project(x))
   if (size(g_raw) /= size(gout)) then
      grad_failed = .true.
      grad_msg = "analytic gradient has wrong length"
      gout = ieee_value(0.0_dp, ieee_quiet_nan)
      return
   end if
   if (.not. finite_vec(g_raw)) then
      grad_failed = .true.
      grad_msg = "analytic gradient is non-finite"
      gout = ieee_value(0.0_dp, ieee_quiet_nan)
      return
   end if
   gout = g_raw / fscale
   return
end if
gr_count = gr_count + 1
call bounded_fd_gradient(x, gout)
end subroutine eval_gradient

subroutine bounded_fd_gradient(x, gout)
real(kind=dp), intent(in) :: x(:)
real(kind=dp), intent(out) :: gout(:)
real(kind=dp), allocatable :: x0(:), x_plus(:), x_minus(:)
real(kind=dp) :: eps, step_i, scale_i, f_base, f_plus, f_minus
integer :: k
logical :: have_base
allocate(x0(size(x)), x_plus(size(x)), x_minus(size(x)))
x0 = project(x)
have_base = .false.
f_base = 0.0_dp
do k = 1, size(x)
   scale_i = 1.0_dp
   if (k <= size(pscale)) scale_i = max(abs(pscale(k)), tiny(1.0_dp))
   step_i = step_eps
   if (present(ndeps_vec)) then
      if (k <= size(ndeps_vec)) step_i = ndeps_vec(k)
   end if
   eps = step_i * (abs(x0(k)) + scale_i)
   if (abs(eps) <= tiny(1.0_dp)) eps = step_eps * (abs(x0(k)) + scale_i)
   if (abs(eps) <= tiny(1.0_dp)) eps = sqrt(epsilon(1.0_dp))
   if (x0(k) + eps <= ub(k) .and. x0(k) - eps >= lb(k)) then
      x_plus = x0
      x_minus = x0
      x_plus(k) = x0(k) + eps
      x_minus(k) = x0(k) - eps
      f_plus = scaled_fn(x_plus)
      f_minus = scaled_fn(x_minus)
      fn_count = fn_count + 2
      gout(k) = (f_plus - f_minus) / (2.0_dp * eps)
   else if (x0(k) + eps <= ub(k)) then
      if (.not. have_base) then
         f_base = scaled_fn(x0)
         fn_count = fn_count + 1
         have_base = .true.
      end if
      x_plus = x0
      x_plus(k) = x0(k) + eps
      f_plus = scaled_fn(x_plus)
      fn_count = fn_count + 1
      gout(k) = (f_plus - f_base) / eps
   else if (x0(k) - eps >= lb(k)) then
      if (.not. have_base) then
         f_base = scaled_fn(x0)
         fn_count = fn_count + 1
         have_base = .true.
      end if
      x_minus = x0
      x_minus(k) = x0(k) - eps
      f_minus = scaled_fn(x_minus)
      fn_count = fn_count + 1
      gout(k) = (f_base - f_minus) / eps
   else
      gout(k) = 0.0_dp
   end if
end do
end subroutine bounded_fd_gradient

pure function finite_vec(x) result(ok)
real(kind=dp), intent(in) :: x(:)
logical :: ok
ok = all(ieee_is_finite(x))
end function finite_vec
end function optim_lbfgsb

function optim_cg(fn, par, maxit, reltol, ndeps, ndeps_vec, fnscale, parscale) result(out)
! Nonlinear conjugate-gradient optimizer for vector-valued parameters.
procedure(optim_vec_objective) :: fn
real(kind=dp), intent(in) :: par(:)
integer, intent(in), optional :: maxit
real(kind=dp), intent(in), optional :: reltol, ndeps, ndeps_vec(:), fnscale, parscale(:)
type(optim_result_t) :: out
integer :: n, max_iter, iter, j
logical :: converged
real(kind=dp) :: f, f_new, step_eps, gtol, alpha, slope, beta, shift, fscale
real(kind=dp), allocatable :: p(:), p_new(:), g(:), g_new(:), d(:), y(:), pscale(:)
n = size(par)
max_iter = 100
if (present(maxit)) max_iter = maxit
step_eps = 1.0e-3_dp
if (present(ndeps)) step_eps = ndeps
fscale = 1.0_dp
if (present(fnscale)) then
   if (abs(fnscale) > tiny(1.0_dp)) fscale = fnscale
end if
gtol = 1.0e-8_dp
if (present(reltol)) gtol = reltol
gtol = max(gtol, sqrt(epsilon(1.0_dp)))
allocate(p(n), p_new(n), g(n), g_new(n), d(n), y(n), pscale(n))
p = par
pscale = 1.0_dp
if (present(parscale)) pscale(1:min(n, size(parscale))) = max(abs(parscale(1:min(n, size(parscale)))), tiny(1.0_dp))
f = scaled_fn(p)
call optim_fd_gradient(scaled_fn, p, step_eps, g, pscale, ndeps_vec)
d = -g
converged = .false.
do iter = 1, max_iter
   if (sqrt(sum(g**2)) < gtol) then
      converged = .true.
      exit
   end if
   slope = dot_product(g, d)
   if (slope >= 0.0_dp) then
      d = -g
      slope = -dot_product(g, g)
   end if
   alpha = 1.0_dp
   do j = 1, 60
      p_new = p + alpha * d
      f_new = scaled_fn(p_new)
      if (f_new <= f + 1.0e-4_dp * alpha * slope) exit
      if (alpha < 1.0e-12_dp) exit
      alpha = 0.5_dp * alpha
   end do
   call optim_fd_gradient(scaled_fn, p_new, step_eps, g_new, pscale, ndeps_vec)
   shift = abs(f - f_new)
   y = g_new - g
   beta = max(0.0_dp, dot_product(g_new, y) / max(dot_product(g, g), tiny(1.0_dp)))
   d = -g_new + beta * d
   p = p_new
   f = f_new
   g = g_new
   if (shift <= gtol * (1.0_dp + abs(f))) then
      converged = .true.
      exit
   end if
end do
out%par = p
out%value = f * fscale
out%convergence = merge(0, 1, converged)
if (converged) then
   out%message = "converged"
else
   out%message = "maximum iterations reached"
end if
contains
pure function scaled_fn(x) result(value)
real(kind=dp), intent(in) :: x(:)
real(kind=dp) :: value
value = fn(x) / fscale
end function scaled_fn
end function optim_cg

function optim_sann(fn, par, maxit, reltol, ndeps, ndeps_vec, fnscale, parscale) result(out)
! Simulated annealing optimizer for vector-valued parameters.
procedure(optim_vec_objective) :: fn
real(kind=dp), intent(in) :: par(:)
integer, intent(in), optional :: maxit
real(kind=dp), intent(in), optional :: reltol, ndeps, ndeps_vec(:), fnscale, parscale(:)
type(optim_result_t) :: out
integer :: n, max_iter, iter, tmax
real(kind=dp) :: f, f_new, best_f, temp, prob, u, fscale
real(kind=dp), allocatable :: p(:), p_new(:), z(:), best_p(:), pscale(:)
n = size(par)
max_iter = 10000
if (present(maxit)) max_iter = maxit
fscale = 1.0_dp
if (present(fnscale)) then
   if (abs(fnscale) > tiny(1.0_dp)) fscale = fnscale
end if
tmax = 10
allocate(p(n), p_new(n), z(n), best_p(n), pscale(n))
p = par
pscale = 1.0_dp
if (present(parscale)) pscale(1:min(n, size(parscale))) = max(abs(parscale(1:min(n, size(parscale)))), tiny(1.0_dp))
f = scaled_fn(p)
best_p = p
best_f = f
do iter = 1, max_iter
   temp = 10.0_dp / log(real(((iter - 1) / tmax) * tmax, kind=dp) + exp(1.0_dp))
   z = rnorm_vec(n)
   p_new = p + 0.1_dp * temp * z * pscale
   f_new = scaled_fn(p_new)
   if (f_new < f) then
      p = p_new
      f = f_new
   else
      prob = exp(min(0.0_dp, (f - f_new) / max(temp, tiny(1.0_dp))))
      call random_number(u)
      if (u < prob) then
         p = p_new
         f = f_new
      end if
   end if
   if (f < best_f) then
      best_p = p
      best_f = f
   end if
end do
out%par = best_p
out%value = best_f * fscale
out%convergence = 0
out%message = "converged"
contains
pure function scaled_fn(x) result(value)
real(kind=dp), intent(in) :: x(:)
real(kind=dp) :: value
value = fn(x) / fscale
end function scaled_fn
end function optim_sann

function optim_nelder_mead(fn, par, maxit, reltol, ndeps, ndeps_vec, fnscale, parscale) result(out)
! Nelder-Mead simplex optimizer for vector-valued parameters.
procedure(optim_vec_objective) :: fn
real(kind=dp), intent(in) :: par(:)
integer, intent(in), optional :: maxit
real(kind=dp), intent(in), optional :: reltol, ndeps, ndeps_vec(:), fnscale, parscale(:)
type(optim_result_t) :: out
integer :: n, max_iter, iter, i, j, best, worst, second
logical :: converged
real(kind=dp) :: gtol, step, fr, fe, fc, spread, fscale
real(kind=dp), allocatable :: simplex(:,:), fvals(:), centroid(:), xr(:), xe(:), xc(:), pscale(:)
n = size(par)
max_iter = 100
if (present(maxit)) max_iter = maxit
fscale = 1.0_dp
if (present(fnscale)) then
   if (abs(fnscale) > tiny(1.0_dp)) fscale = fnscale
end if
gtol = 1.0e-8_dp
if (present(reltol)) gtol = reltol
gtol = max(gtol, sqrt(epsilon(1.0_dp)))
allocate(simplex(n,n+1), fvals(n+1), centroid(n), xr(n), xe(n), xc(n), pscale(n))
pscale = 1.0_dp
if (present(parscale)) pscale(1:min(n, size(parscale))) = max(abs(parscale(1:min(n, size(parscale)))), tiny(1.0_dp))
simplex(:,1) = par
do i = 1, n
   simplex(:,i+1) = par
   step = 0.05_dp * (abs(par(i)) + pscale(i))
   simplex(i,i+1) = simplex(i,i+1) + step
end do
do j = 1, n + 1
   fvals(j) = scaled_fn(simplex(:,j))
end do
converged = .false.
do iter = 1, max_iter
   best = 1
   worst = 1
   do j = 2, n + 1
      if (fvals(j) < fvals(best)) best = j
      if (fvals(j) > fvals(worst)) worst = j
   end do
   second = merge(2, 1, worst == 1)
   do j = 1, n + 1
      if (j /= worst .and. fvals(j) > fvals(second)) second = j
   end do
   spread = maxval(abs(fvals - fvals(best)))
   if (spread <= gtol * (1.0_dp + abs(fvals(best)))) then
      converged = .true.
      exit
   end if
   centroid = 0.0_dp
   do j = 1, n + 1
      if (j /= worst) centroid = centroid + simplex(:,j)
   end do
   centroid = centroid / real(n, kind=dp)
   xr = centroid + (centroid - simplex(:,worst))
   fr = scaled_fn(xr)
   if (fr < fvals(best)) then
      xe = centroid + 2.0_dp * (xr - centroid)
      fe = scaled_fn(xe)
      if (fe < fr) then
         simplex(:,worst) = xe
         fvals(worst) = fe
      else
         simplex(:,worst) = xr
         fvals(worst) = fr
      end if
   else if (fr < fvals(second)) then
      simplex(:,worst) = xr
      fvals(worst) = fr
   else
      xc = centroid + 0.5_dp * (simplex(:,worst) - centroid)
      fc = scaled_fn(xc)
      if (fc < fvals(worst)) then
         simplex(:,worst) = xc
         fvals(worst) = fc
      else
         do j = 1, n + 1
            if (j /= best) then
               simplex(:,j) = simplex(:,best) + 0.5_dp * (simplex(:,j) - simplex(:,best))
               fvals(j) = scaled_fn(simplex(:,j))
            end if
         end do
      end if
   end if
end do
best = 1
do j = 2, n + 1
   if (fvals(j) < fvals(best)) best = j
end do
out%par = simplex(:,best)
out%value = fvals(best) * fscale
out%convergence = merge(0, 1, converged)
if (converged) then
   out%message = "converged"
else
   out%message = "maximum iterations reached"
end if
contains
pure function scaled_fn(x) result(value)
real(kind=dp), intent(in) :: x(:)
real(kind=dp) :: value
value = fn(x) / fscale
end function scaled_fn
end function optim_nelder_mead

function constr_optim_bfgs(fn, theta, ui, ci, maxit, reltol, ndeps) result(out)
! Approximate R constrOptim() with a logarithmic barrier and BFGS inner solves.
procedure(optim_vec_objective) :: fn
real(kind=dp), intent(in) :: theta(:), ui(:,:), ci(:)
integer, intent(in), optional :: maxit
real(kind=dp), intent(in), optional :: reltol, ndeps
type(optim_result_t) :: out, inner
integer :: outer, max_outer, inner_maxit
real(kind=dp) :: mu, tol, step_eps
real(kind=dp), allocatable :: p(:), slack(:)
max_outer = 8
inner_maxit = 100
if (present(maxit)) inner_maxit = max(10, maxit / max_outer)
tol = 1.0e-8_dp
if (present(reltol)) tol = reltol
step_eps = 1.0e-3_dp
if (present(ndeps)) step_eps = ndeps
p = theta
slack = matmul(ui, p) - ci
if (any(slack <= 0.0_dp)) then
   out%par = p
   out%value = huge(1.0_dp)
   out%convergence = 51
   out%message = "initial value is not feasible"
   return
end if
mu = 1.0_dp
out%par = p
out%value = fn(p)
out%convergence = 1
out%message = "maximum iterations reached"
do outer = 1, max_outer
   inner = optim_bfgs(barrier_obj, p, maxit=inner_maxit, reltol=tol, ndeps=step_eps)
   p = inner%par
   slack = matmul(ui, p) - ci
   if (any(slack <= 0.0_dp)) exit
   out%par = p
   out%value = fn(p)
   out%convergence = inner%convergence
   out%message = inner%message
   if (mu < tol) exit
   mu = 0.2_dp * mu
end do
contains
pure function barrier_obj(x) result(value)
real(kind=dp), intent(in) :: x(:)
real(kind=dp) :: value
real(kind=dp), allocatable :: s(:)
s = matmul(ui, x) - ci
if (any(s <= 0.0_dp)) then
   value = huge(1.0_dp) / 100.0_dp + sum(min(0.0_dp, s)**2) * 1.0e20_dp
else
   value = fn(x) - mu * sum(log(s))
end if
end function barrier_obj
end function constr_optim_bfgs

function constr_optim_nelder_mead(fn, theta, ui, ci, maxit, reltol, ndeps) result(out)
! Approximate R constrOptim() with a logarithmic barrier and Nelder-Mead inner solves.
procedure(optim_vec_objective) :: fn
real(kind=dp), intent(in) :: theta(:), ui(:,:), ci(:)
integer, intent(in), optional :: maxit
real(kind=dp), intent(in), optional :: reltol, ndeps
type(optim_result_t) :: out, inner
integer :: outer, max_outer, inner_maxit
real(kind=dp) :: mu, tol, step_eps
real(kind=dp), allocatable :: p(:), slack(:)
max_outer = 8
inner_maxit = 100
if (present(maxit)) inner_maxit = max(10, maxit / max_outer)
tol = 1.0e-8_dp
if (present(reltol)) tol = reltol
step_eps = 1.0e-3_dp
if (present(ndeps)) step_eps = ndeps
p = theta
slack = matmul(ui, p) - ci
if (any(slack <= 0.0_dp)) then
   out%par = p
   out%value = huge(1.0_dp)
   out%convergence = 51
   out%message = "initial value is not feasible"
   return
end if
mu = 1.0_dp
out%par = p
out%value = fn(p)
out%convergence = 1
out%message = "maximum iterations reached"
do outer = 1, max_outer
   inner = optim_nelder_mead(barrier_obj, p, maxit=inner_maxit, reltol=tol, ndeps=step_eps)
   p = inner%par
   slack = matmul(ui, p) - ci
   if (any(slack <= 0.0_dp)) exit
   out%par = p
   out%value = fn(p)
   out%convergence = inner%convergence
   out%message = inner%message
   if (mu < tol) exit
   mu = 0.2_dp * mu
end do
contains
pure function barrier_obj(x) result(value)
real(kind=dp), intent(in) :: x(:)
real(kind=dp) :: value
real(kind=dp), allocatable :: s(:)
s = matmul(ui, x) - ci
if (any(s <= 0.0_dp)) then
   value = huge(1.0_dp) / 100.0_dp + sum(min(0.0_dp, s)**2) * 1.0e20_dp
else
   value = fn(x) - mu * sum(log(s))
end if
end function barrier_obj
end function constr_optim_nelder_mead

subroutine optim_fd_gradient(fn, p, step_eps, g, parscale, ndeps_vec)
procedure(optim_vec_objective) :: fn
real(kind=dp), intent(in) :: p(:), step_eps
real(kind=dp), intent(in), optional :: parscale(:)
real(kind=dp), intent(in), optional :: ndeps_vec(:)
real(kind=dp), intent(out) :: g(:)
real(kind=dp), allocatable :: p_tmp(:)
real(kind=dp) :: eps, f_plus, f_minus, scale_i, step_i
integer :: i
allocate(p_tmp(size(p)))
do i = 1, size(p)
   scale_i = 1.0_dp
   if (present(parscale)) then
      if (i <= size(parscale)) scale_i = max(abs(parscale(i)), tiny(1.0_dp))
   end if
   step_i = step_eps
   if (present(ndeps_vec)) then
      if (i <= size(ndeps_vec)) step_i = ndeps_vec(i)
   end if
   eps = step_i * (abs(p(i)) + scale_i)
   if (abs(eps) <= tiny(1.0_dp)) eps = step_eps * (abs(p(i)) + scale_i)
   if (abs(eps) <= tiny(1.0_dp)) eps = sqrt(epsilon(1.0_dp))
   p_tmp = p
   p_tmp(i) = p_tmp(i) + eps
   f_plus = fn(p_tmp)
   p_tmp = p
   p_tmp(i) = p_tmp(i) - eps
   f_minus = fn(p_tmp)
   g(i) = (f_plus - f_minus) / (2.0_dp * eps)
end do
end subroutine optim_fd_gradient

subroutine optim_fd_hessian(fn, p, step_eps, hess, parscale, ndeps_vec)
procedure(optim_vec_objective) :: fn
real(kind=dp), intent(in) :: p(:), step_eps
real(kind=dp), intent(in), optional :: parscale(:)
real(kind=dp), intent(in), optional :: ndeps_vec(:)
real(kind=dp), allocatable, intent(out) :: hess(:,:)
real(kind=dp), allocatable :: p_pp(:), p_pm(:), p_mp(:), p_mm(:), step(:)
real(kind=dp) :: f0, f_plus, f_minus, f_pp, f_pm, f_mp, f_mm, scale_i, step_i
integer :: i, j, n
n = size(p)
allocate(hess(n,n), p_pp(n), p_pm(n), p_mp(n), p_mm(n), step(n))
do i = 1, n
   scale_i = 1.0_dp
   if (present(parscale)) then
      if (i <= size(parscale)) scale_i = max(abs(parscale(i)), tiny(1.0_dp))
   end if
   step_i = step_eps
   if (present(ndeps_vec)) then
      if (i <= size(ndeps_vec)) step_i = ndeps_vec(i)
   end if
   step(i) = step_i * (abs(p(i)) + scale_i)
   if (abs(step(i)) <= tiny(1.0_dp)) step(i) = step_eps * (abs(p(i)) + scale_i)
   if (abs(step(i)) <= tiny(1.0_dp)) step(i) = sqrt(epsilon(1.0_dp))
end do
f0 = fn(p)
do i = 1, n
   p_pp = p
   p_mm = p
   p_pp(i) = p_pp(i) + step(i)
   p_mm(i) = p_mm(i) - step(i)
   f_plus = fn(p_pp)
   f_minus = fn(p_mm)
   hess(i,i) = (f_plus - 2.0_dp * f0 + f_minus) / (step(i)**2)
end do
do i = 1, n
   do j = i + 1, n
      p_pp = p
      p_pm = p
      p_mp = p
      p_mm = p
      p_pp(i) = p_pp(i) + step(i)
      p_pp(j) = p_pp(j) + step(j)
      p_pm(i) = p_pm(i) + step(i)
      p_pm(j) = p_pm(j) - step(j)
      p_mp(i) = p_mp(i) - step(i)
      p_mp(j) = p_mp(j) + step(j)
      p_mm(i) = p_mm(i) - step(i)
      p_mm(j) = p_mm(j) - step(j)
      f_pp = fn(p_pp)
      f_pm = fn(p_pm)
      f_mp = fn(p_mp)
      f_mm = fn(p_mm)
      hess(i,j) = (f_pp - f_pm - f_mp + f_mm) / (4.0_dp * step(i) * step(j))
      hess(j,i) = hess(i,j)
   end do
end do
end subroutine optim_fd_hessian

function integrate(fn, lower, upper, rel_tol, subdivisions) result(out)
! Approximate R integrate() for smooth scalar functions using composite Simpson rules.
procedure(integrate_objective) :: fn
real(kind=dp), intent(in) :: lower
real(kind=dp), intent(in) :: upper
real(kind=dp), intent(in), optional :: rel_tol
integer, intent(in), optional :: subdivisions
type(integrate_result_t) :: out
integer :: n, nmax
real(kind=dp) :: prev, curr, tol, scale, inf_threshold
integer :: mode

if (present(rel_tol)) then
   tol = rel_tol
else
   tol = sqrt(epsilon(1.0_dp))
end if
if (present(subdivisions)) then
   nmax = max(2, subdivisions)
else
   nmax = 100
end if
if (mod(nmax, 2) /= 0) nmax = nmax + 1
inf_threshold = sqrt(huge(1.0_dp))
if (abs(lower) > inf_threshold .and. abs(upper) > inf_threshold) then
   mode = 3
else if (abs(upper) > inf_threshold) then
   mode = 1
else if (abs(lower) > inf_threshold) then
   mode = 2
else
   mode = 0
end if
if (mode /= 0) nmax = max(nmax, 16 * nmax)
n = 2
prev = simpson_integral_mapped(fn, lower, upper, n, mode)
do
   n = min(2 * n, nmax)
   curr = simpson_integral_mapped(fn, lower, upper, n, mode)
   scale = max(1.0_dp, abs(curr))
   if (abs(curr - prev) <= tol * scale .or. n >= nmax) exit
   prev = curr
end do
out%value = curr
out%abs_error = abs(curr - prev) / 15.0_dp
out%subdivisions = n
out%message = 0
end function integrate

function simpson_integral_mapped(fn, a, b, n_in, mode) result(out)
! Composite Simpson integral, applying finite transforms for improper bounds.
procedure(integrate_objective) :: fn
real(kind=dp), intent(in) :: a
real(kind=dp), intent(in) :: b
integer, intent(in) :: n_in
integer, intent(in) :: mode
real(kind=dp) :: out
integer :: i, n
real(kind=dp) :: h, t, s
n = max(2, n_in)
if (mod(n, 2) /= 0) n = n + 1
if (mode == 0) then
   out = simpson_integral(fn, a, b, n)
   return
end if
h = 1.0_dp / real(n, kind=dp)
s = integrate_mapped_value(fn, a, b, mode, 0.0_dp) + integrate_mapped_value(fn, a, b, mode, 1.0_dp)
do i = 1, n - 1
   t = h * real(i, kind=dp)
   if (mod(i, 2) == 0) then
      s = s + 2.0_dp * integrate_mapped_value(fn, a, b, mode, t)
   else
      s = s + 4.0_dp * integrate_mapped_value(fn, a, b, mode, t)
   end if
end do
out = s * h / 3.0_dp
end function simpson_integral_mapped

function integrate_mapped_value(fn, a, b, mode, t) result(out)
! Mapped integrand value on t in [0, 1] for improper intervals.
procedure(integrate_objective) :: fn
real(kind=dp), intent(in) :: a
real(kind=dp), intent(in) :: b
integer, intent(in) :: mode
real(kind=dp), intent(in) :: t
real(kind=dp) :: out
real(kind=dp) :: x, jac, tt, pi
pi = acos(-1.0_dp)
select case (mode)
case (1)
   if (t <= 0.0_dp) then
      tt = 0.0_dp
   else if (t >= 1.0_dp) then
      tt = 1.0_dp - epsilon(1.0_dp)
   else
      tt = t
   end if
   x = a + tt / (1.0_dp - tt)
   jac = 1.0_dp / (1.0_dp - tt)**2
case (2)
   if (t <= 0.0_dp) then
      tt = epsilon(1.0_dp)
   else if (t >= 1.0_dp) then
      tt = 1.0_dp
   else
      tt = t
   end if
   x = b - (1.0_dp - tt) / tt
   jac = 1.0_dp / tt**2
case default
   if (t <= 0.0_dp) then
      tt = epsilon(1.0_dp)
   else if (t >= 1.0_dp) then
      tt = 1.0_dp - epsilon(1.0_dp)
   else
      tt = t
   end if
   x = tan(pi * (tt - 0.5_dp))
   jac = pi / cos(pi * (tt - 0.5_dp))**2
end select
out = fn(x) * jac
if (.not. ieee_is_finite(out)) out = 0.0_dp
end function integrate_mapped_value

function simpson_integral(fn, a, b, n_in) result(out)
! Composite Simpson integral with an even number of panels.
procedure(integrate_objective) :: fn
real(kind=dp), intent(in) :: a
real(kind=dp), intent(in) :: b
integer, intent(in) :: n_in
real(kind=dp) :: out
integer :: i, n
real(kind=dp) :: h, x, s, y
n = max(2, n_in)
if (mod(n, 2) /= 0) n = n + 1
h = (b - a) / real(n, kind=dp)
s = 0.0_dp
y = fn(a)
if (ieee_is_finite(y)) s = s + y
y = fn(b)
if (ieee_is_finite(y)) s = s + y
do i = 1, n - 1
   x = a + h * real(i, kind=dp)
   y = fn(x)
   if (.not. ieee_is_finite(y)) cycle
   if (mod(i, 2) == 0) then
      s = s + 2.0_dp * y
   else
      s = s + 4.0_dp * y
   end if
end do
out = s * h / 3.0_dp
end function simpson_integral

subroutine print_integrate_result(fit)
! Print an integrate result in a compact R-like form.
type(integrate_result_t), intent(in) :: fit
write(*,"(g0,a,g0,a,i0)") fit%value, " with absolute error < ", fit%abs_error, ", subdivisions = ", fit%subdivisions
end subroutine print_integrate_result

real(kind=dp) function r_elapsed() result(out)
! Wall-clock elapsed time in seconds.
integer :: count, rate
call system_clock(count=count, count_rate=rate)
if (rate > 0) then
   out = real(count, kind=dp) / real(rate, kind=dp)
else
   out = 0.0_dp
end if
end function r_elapsed

real(kind=dp) function sys_time() result(out)
! Current date-time as seconds since 1970-01-01, local calendar fields.
integer :: vals(8)
integer :: days
call date_and_time(values=vals)
days = date_days_from_civil(vals(1), vals(2), vals(3))
out = real(days, kind=dp) * 86400.0_dp + real(vals(5) * 3600 + vals(6) * 60 + vals(7), kind=dp) &
   & + real(vals(8), kind=dp) / 1000.0_dp
end function sys_time

integer function sys_date() result(out)
! Current date as days since 1970-01-01.
integer :: vals(8)
call date_and_time(values=vals)
out = date_days_from_civil(vals(1), vals(2), vals(3))
end function sys_date

function sys_timezone() result(out)
! Best-effort Tier-1 timezone name.
character(len=:), allocatable :: out
out = "UTC"
end function sys_timezone

function sys_date_string() result(out)
! Current date-time as a character string.
character(len=:), allocatable :: out
out = sys_time_format(sys_time(), "%Y-%m-%d %H:%M:%S")
end function sys_date_string

function proc_time_vec() result(out)
! R-like proc.time vector with elapsed time in the third slot.
real(kind=dp) :: out(5)
out = 0.0_dp
out(3) = r_elapsed()
end function proc_time_vec

function sys_getenv(name) result(out)
! R-like Sys.getenv: value of an environment variable ("" when unset).
character(len=*), intent(in) :: name
character(len=:), allocatable :: out
integer :: len_env, stat_env
call get_environment_variable(name, length=len_env, status=stat_env)
if (stat_env /= 0 .or. len_env <= 0) then
   out = ""
   return
end if
allocate(character(len=len_env) :: out)
call get_environment_variable(name, value=out)
end function sys_getenv

function file_rename(from, to) result(ok)
! R-like file.rename via copy-and-delete; returns .true. on success.
character(len=*), intent(in) :: from, to
logical :: ok
integer :: u_in, u_out, ios
integer(kind=int64) :: nbytes, i
character(len=1) :: byte
ok = .false.
open(newunit=u_in, file=from, access="stream", form="unformatted", &
   & status="old", action="read", iostat=ios)
if (ios /= 0) return
inquire(unit=u_in, size=nbytes)
open(newunit=u_out, file=to, access="stream", form="unformatted", &
   & status="replace", action="write", iostat=ios)
if (ios /= 0) then
   close(u_in)
   return
end if
do i = 1, nbytes
   read(u_in, iostat=ios) byte
   if (ios /= 0) exit
   write(u_out) byte
end do
close(u_out)
close(u_in, status="delete")
ok = .true.
end function file_rename

pure integer function hex_digit_value(ch) result(out)
character(len=1), intent(in) :: ch
select case (ch)
case ("0":"9")
   out = iachar(ch) - iachar("0")
case ("a":"f")
   out = iachar(ch) - iachar("a") + 10
case ("A":"F")
   out = iachar(ch) - iachar("A") + 10
case default
   out = -1
end select
end function hex_digit_value

pure subroutine parse_hex_real(text, out, ok)
character(len=*), intent(in) :: text
real(kind=dp), intent(out) :: out
logical, intent(out) :: ok
integer :: i, n, digit, exponent_value, exponent_sign, signed_exponent
real(kind=dp) :: value, fraction, sign_value
logical :: have_digit, have_exponent_digit
n = len(text)
i = 1
sign_value = 1.0_dp
if (i <= n) then
   if (text(i:i) == "+") then
      i = i + 1
   else if (text(i:i) == "-") then
      sign_value = -1.0_dp
      i = i + 1
   end if
end if
ok = .false.
out = ieee_value(0.0_dp, ieee_quiet_nan)
if (i + 1 > n) return
if (text(i:i) /= "0" .or. (text(i + 1:i + 1) /= "x" .and. text(i + 1:i + 1) /= "X")) return
i = i + 2
value = 0.0_dp
have_digit = .false.
do while (i <= n)
   digit = hex_digit_value(text(i:i))
   if (digit < 0) exit
   value = 16.0_dp * value + real(digit, kind=dp)
   have_digit = .true.
   i = i + 1
end do
if (i <= n) then
   if (text(i:i) == ".") then
      i = i + 1
      fraction = 1.0_dp / 16.0_dp
      do while (i <= n)
         digit = hex_digit_value(text(i:i))
         if (digit < 0) exit
         value = value + fraction * real(digit, kind=dp)
         fraction = fraction / 16.0_dp
         have_digit = .true.
         i = i + 1
      end do
   end if
end if
if (.not. have_digit) return
signed_exponent = 0
if (i <= n) then
   if (text(i:i) == "p" .or. text(i:i) == "P") then
      i = i + 1
      exponent_sign = 1
      if (i <= n) then
         if (text(i:i) == "+") then
            i = i + 1
         else if (text(i:i) == "-") then
            exponent_sign = -1
            i = i + 1
         end if
      end if
      exponent_value = 0
      have_exponent_digit = .false.
      do while (i <= n)
         if (text(i:i) < "0" .or. text(i:i) > "9") exit
         exponent_value = min(100000, 10 * exponent_value + iachar(text(i:i)) - iachar("0"))
         have_exponent_digit = .true.
         i = i + 1
      end do
      if (.not. have_exponent_digit) return
      signed_exponent = exponent_sign * exponent_value
   end if
end if
if (i <= n) return
if (signed_exponent > maxexponent(1.0_dp) + digits(1.0_dp)) then
   if (value == 0.0_dp) then
      out = sign_value * 0.0_dp
   else
      out = sign_value * ieee_value(0.0_dp, ieee_positive_inf)
   end if
else if (signed_exponent < minexponent(1.0_dp) - digits(1.0_dp)) then
   out = sign_value * 0.0_dp
else
   out = sign_value * scale(value, signed_exponent)
end if
ok = .true.
end subroutine parse_hex_real

pure elemental function str_to_real(s) result(out)
! R-like as.numeric() of a character value: parse to real, NA(nan) on failure.
character(len=*), intent(in) :: s
real(kind=dp) :: out
integer :: ios, i, n, sign_offset
logical :: have_digit, parsed_hex
character(len=:), allocatable :: text, lower
text = trim(adjustl(s))
lower = tolower(text)
if (lower == "inf" .or. lower == "+inf" .or. lower == "infinity" .or. &
   lower == "+infinity") then
   out = ieee_value(0.0_dp, ieee_positive_inf)
   return
else if (lower == "-inf" .or. lower == "-infinity") then
   out = ieee_value(0.0_dp, ieee_negative_inf)
   return
else if (lower == "nan" .or. lower == "+nan" .or. lower == "-nan") then
   out = ieee_value(0.0_dp, ieee_quiet_nan)
   return
end if
n = len(text)
sign_offset = 0
if (n > 0) then
   if (text(1:1) == "+" .or. text(1:1) == "-") sign_offset = 1
end if
if (n >= sign_offset + 2) then
   if (text(sign_offset + 1:sign_offset + 1) == "0" .and. &
      (text(sign_offset + 2:sign_offset + 2) == "x" .or. &
      text(sign_offset + 2:sign_offset + 2) == "X")) then
      call parse_hex_real(text, out, parsed_hex)
      if (.not. parsed_hex) out = ieee_value(0.0_dp, ieee_quiet_nan)
      return
   end if
end if
i = 1
if (n > 0) then
   if (text(i:i) == "+" .or. text(i:i) == "-") i = i + 1
end if
have_digit = .false.
do while (i <= n)
   if (text(i:i) < "0" .or. text(i:i) > "9") exit
   have_digit = .true.
   i = i + 1
end do
if (i <= n) then
   if (text(i:i) == ".") then
      i = i + 1
      do while (i <= n)
         if (text(i:i) < "0" .or. text(i:i) > "9") exit
         have_digit = .true.
         i = i + 1
      end do
   end if
end if
if (have_digit .and. i <= n) then
   if (text(i:i) == "e" .or. text(i:i) == "E") then
      i = i + 1
      if (i <= n) then
         if (text(i:i) == "+" .or. text(i:i) == "-") i = i + 1
      end if
      have_digit = .false.
      do while (i <= n)
         if (text(i:i) < "0" .or. text(i:i) > "9") exit
         have_digit = .true.
         i = i + 1
      end do
   end if
end if
if (.not. have_digit .or. i <= n) then
   out = ieee_value(0.0_dp, ieee_quiet_nan)
   return
end if
read(text, *, iostat=ios) out
if (ios /= 0) out = ieee_value(0.0_dp, ieee_quiet_nan)
end function str_to_real

pure elemental integer function str_to_int(s) result(out)
! R-like as.integer() of a character value: parse to integer, NA on failure.
character(len=*), intent(in) :: s
real(kind=dp) :: parsed, truncated
parsed = str_to_real(s)
if (.not. ieee_is_finite(parsed)) then
   out = -huge(0)
   return
end if
truncated = aint(parsed)
if (truncated > real(huge(0), kind=dp) .or. truncated < -real(huge(0), kind=dp)) then
   out = -huge(0)
   return
end if
out = int(truncated)
end function str_to_int

pure function inttobits(n) result(out)
! R-like intToBits(n): the 32 bits of n, least-significant first, as 0/1.
integer, intent(in) :: n
integer :: out(32)
integer :: k
do k = 1, 32
   out(k) = iand(ishft(n, 1 - k), 1)
end do
end function inttobits

function unlink_recursive(path) result(out)
! R-like unlink(path, recursive=TRUE): delete a file or directory tree.
! Returns 0 on success, 1 on failure (R convention).
character(len=*), intent(in) :: path
integer :: out
integer :: stat
if (file_exists(path) .and. .not. dir_exists_scalar(path)) then
   if (file_remove(path)) then
      out = 0
   else
      out = 1
   end if
   return
end if
if (is_windows_path_env()) then
   call execute_command_line('cmd /c rmdir /s /q "' // trim(path) // '" >nul 2>nul', &
      & wait=.true., exitstat=stat)
else
   call execute_command_line('rm -rf "' // trim(path) // '" >/dev/null 2>&1', &
      & wait=.true., exitstat=stat)
end if
out = merge(0, 1, stat == 0)
end function unlink_recursive

function as_octmode(x) result(out)
! R-like as.octmode: integer rendered as octal digits.
integer, intent(in) :: x
character(len=:), allocatable :: out
character(len=32) :: buf
write(buf, "(o0)") x
out = trim(buf)
end function as_octmode

function as_hexmode(x) result(out)
! R-like as.hexmode: integer rendered as lowercase hex digits.
integer, intent(in) :: x
character(len=:), allocatable :: out
character(len=32) :: buf
write(buf, "(z0)") x
out = tolower(trim(buf))
end function as_hexmode

function as_roman(x) result(out)
! R-like as.roman: Roman-numeral string for 1..3899.
integer, intent(in) :: x
character(len=:), allocatable :: out
integer, parameter :: vals(13) = [1000, 900, 500, 400, 100, 90, 50, 40, 10, 9, 5, 4, 1]
character(len=2), parameter :: syms(13) = ["M ", "CM", "D ", "CD", "C ", "XC", &
   & "L ", "XL", "X ", "IX", "V ", "IV", "I "]
integer :: n, k
if (x < 1 .or. x > 3899) then
   out = "NA"
   return
end if
n = x
out = ""
do k = 1, 13
   do while (n >= vals(k))
      out = out // trim(syms(k))
      n = n - vals(k)
   end do
end do
end function as_roman

subroutine sys_sleep(seconds)
! Busy-wait sleep for Tier-1 Sys.sleep support.
real(kind=dp), intent(in) :: seconds
real(kind=dp) :: t0
t0 = r_elapsed()
do while (r_elapsed() - t0 < max(0.0_dp, seconds))
end do
end subroutine sys_sleep

pure integer function date_digit(ch) result(out)
! Convert one decimal character to an integer digit.
character(len=1), intent(in) :: ch ! character to decode
integer :: k
k = iachar(ch) - iachar("0")
if (k < 0 .or. k > 9) then
   out = 0
else
   out = k
end if
end function date_digit

pure integer function date_int_slice(s, i1, i2) result(out)
! Decode a decimal substring.
character(len=*), intent(in) :: s ! source string
integer, intent(in) :: i1 ! first character position
integer, intent(in) :: i2 ! last character position
integer :: i
out = 0
do i = i1, i2
   out = 10 * out + date_digit(s(i:i))
end do
end function date_int_slice

pure integer function date_days_from_civil(y, m, d) result(z)
! Convert a Gregorian date to days since 1970-01-01.
integer, intent(in) :: y ! calendar year
integer, intent(in) :: m ! calendar month
integer, intent(in) :: d ! calendar day
integer :: yy, era, yoe, doy, doe, mp
yy = y
if (m <= 2) yy = yy - 1
if (yy >= 0) then
   era = yy / 400
else
   era = (yy - 399) / 400
end if
yoe = yy - era * 400
mp = m
if (mp > 2) then
   mp = mp - 3
else
   mp = mp + 9
end if
doy = (153 * mp + 2) / 5 + d - 1
doe = yoe * 365 + yoe / 4 - yoe / 100 + doy
z = era * 146097 + doe - 719468
end function date_days_from_civil

pure subroutine date_civil_from_days(z, y, m, d)
! Convert days since 1970-01-01 to Gregorian date fields.
integer, intent(in) :: z ! days since 1970-01-01
integer, intent(out) :: y ! calendar year
integer, intent(out) :: m ! calendar month
integer, intent(out) :: d ! calendar day
integer :: zz, era, doe, yoe, doy, mp
zz = z + 719468
if (zz >= 0) then
   era = zz / 146097
else
   era = (zz - 146096) / 146097
end if
doe = zz - era * 146097
yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365
y = yoe + era * 400
doy = doe - (365 * yoe + yoe / 4 - yoe / 100)
mp = (5 * doy + 2) / 153
d = doy - (153 * mp + 2) / 5 + 1
if (mp < 10) then
   m = mp + 3
else
   m = mp - 9
end if
if (m <= 2) y = y + 1
end subroutine date_civil_from_days

pure function date_from_iso(s) result(out)
! Parse an ISO yyyy-mm-dd date as days since 1970-01-01.
character(len=*), intent(in) :: s ! ISO date string
integer :: out
integer :: y, m, d
y = date_int_slice(s, 1, 4)
m = date_int_slice(s, 6, 7)
d = date_int_slice(s, 9, 10)
out = date_days_from_civil(y, m, d)
end function date_from_iso

pure function date_from_iso_vec(s) result(out)
! Parse a vector of ISO yyyy-mm-dd strings as day counts.
character(len=*), intent(in) :: s(:) ! ISO date strings
integer :: out(size(s))
integer :: i
do i = 1, size(s)
   out(i) = date_from_iso(s(i))
end do
end function date_from_iso_vec

pure function date_from_yyyymmdd_vec(x) result(out)
! Parse numeric yyyymmdd values as day counts.
real(kind=dp), intent(in) :: x(:) ! dates encoded as yyyymmdd
integer :: out(size(x))
integer :: i, v, y, m, d
do i = 1, size(x)
   v = nint(x(i))
   y = v / 10000
   m = mod(v / 100, 100)
   d = mod(v, 100)
   out(i) = date_days_from_civil(y, m, d)
end do
end function date_from_yyyymmdd_vec

pure function date_to_char(x) result(out)
! Format a day count as yyyy-mm-dd.
integer, intent(in) :: x ! days since 1970-01-01
character(len=10) :: out
integer :: y, m, d
call date_civil_from_days(x, y, m, d)
write(out, "(i4.4,a,i2.2,a,i2.2)") y, "-", m, "-", d
end function date_to_char

function sys_time_format(x, fmt) result(out)
! Format a POSIXct-like second count for common Tier-1 strftime patterns.
real(kind=dp), intent(in) :: x
character(len=*), intent(in) :: fmt
character(len=:), allocatable :: out
character(len=64) :: buf
integer :: days, rem, y, m, d, hh, mm, ss
days = floor(x / 86400.0_dp)
rem = nint(x - real(days, kind=dp) * 86400.0_dp)
if (rem < 0) then
   rem = rem + 86400
   days = days - 1
end if
if (rem >= 86400) then
   rem = rem - 86400
   days = days + 1
end if
call date_civil_from_days(days, y, m, d)
hh = rem / 3600
mm = mod(rem / 60, 60)
ss = mod(rem, 60)
select case (trim(fmt))
case ("%Y-%m-%d")
   write(buf, "(i4.4,a,i2.2,a,i2.2)") y, "-", m, "-", d
case ("%H:%M:%S")
   write(buf, "(i2.2,a,i2.2,a,i2.2)") hh, ":", mm, ":", ss
case ("%Y-%m-%d %H:%M:%S", "")
   write(buf, "(i4.4,a,i2.2,a,i2.2,1x,i2.2,a,i2.2,a,i2.2)") y, "-", m, "-", d, hh, ":", mm, ":", ss
case ("%Y-%m-%d %H:%M:%S %Z")
   write(buf, "(i4.4,a,i2.2,a,i2.2,1x,i2.2,a,i2.2,a,i2.2,a)") y, "-", m, "-", d, hh, ":", mm, ":", ss, " UTC"
case default
   write(buf, "(i4.4,a,i2.2,a,i2.2,1x,i2.2,a,i2.2,a,i2.2)") y, "-", m, "-", d, hh, ":", mm, ":", ss
end select
out = trim(buf)
end function sys_time_format

pure function date_to_char_vec(x) result(out)
! Format day counts as yyyy-mm-dd strings.
integer, intent(in) :: x(:) ! days since 1970-01-01
character(len=10) :: out(size(x))
integer :: i
do i = 1, size(x)
   out(i) = date_to_char(x(i))
end do
end function date_to_char_vec

pure function date_format(x, fmt) result(out)
! Format a day count for a small subset of R date formats.
integer, intent(in) :: x ! days since 1970-01-01
character(len=*), intent(in) :: fmt ! one of %Y, %m, %d, or %F
character(len=10) :: out
integer :: y, m, d
call date_civil_from_days(x, y, m, d)
select case (trim(fmt))
case ("%Y")
   write(out, "(i4.4)") y
case ("%Y-%m")
   write(out, "(i4.4, '-', i2.2)") y, m
case ("%m")
   write(out, "(i2.2)") m
case ("%d")
   write(out, "(i2.2)") d
case ("%d/%m/%Y")
   write(out, "(i2.2, '/', i2.2, '/', i4.4)") d, m, y
case default
   out = date_to_char(x)
end select
end function date_format

pure function date_format_vec(x, fmt) result(out)
! Format day counts for a small subset of R date formats.
integer, intent(in) :: x(:) ! days since 1970-01-01
character(len=*), intent(in) :: fmt ! one of %Y, %m, %d, or %F
character(len=10) :: out(size(x))
integer :: i
do i = 1, size(x)
   out(i) = date_format(x(i), fmt)
end do
end function date_format_vec

subroutine print_date(x)
! Print one R-like Date value.
integer, intent(in) :: x ! days since 1970-01-01
write(*, "(a)") trim(date_to_char(x))
end subroutine print_date

subroutine print_date_vector(x)
! Print an R-like Date vector.
integer, intent(in) :: x(:) ! days since 1970-01-01
call print_char_vector(date_to_char_vec(x))
end subroutine print_date_vector

pure function date_seq_day(from, to, by) result(out)
! Create a daily Date sequence between two endpoints.
integer, intent(in) :: from ! first date day count
integer, intent(in) :: to ! final date day count
integer, intent(in), optional :: by ! step in days
integer, allocatable :: out(:)
integer :: step, n, i
step = 1
if (present(by)) step = by
if (step == 0) then
   allocate(out(0))
   return
end if
if ((step > 0 .and. from > to) .or. (step < 0 .and. from < to)) then
   allocate(out(0))
   return
end if
n = abs((to - from) / step) + 1
allocate(out(n))
do i = 1, n
   out(i) = from + (i - 1) * step
end do
end function date_seq_day

pure function date_seq_length(from, by, n) result(out)
! Create a Date sequence with a requested length.
integer, intent(in) :: from ! first date day count
integer, intent(in) :: by ! step in days
integer, intent(in) :: n ! requested length
integer, allocatable :: out(:)
integer :: i
allocate(out(max(0, n)))
do i = 1, size(out)
   out(i) = from + (i - 1) * by
end do
end function date_seq_length

pure function date_range(x) result(out)
! Return minimum and maximum Date day counts.
integer, intent(in) :: x(:) ! Date day counts
integer :: out(2)
out = [minval(x), maxval(x)]
end function date_range

pure function r_character(n) result(out)
! Allocate an R-like character vector initialized to empty strings.
integer, intent(in) :: n ! requested vector length
character(len=:), allocatable :: out(:)
! Deferred-length elements cannot grow independently after allocation.
! Reserve a practical width so later x[i] <- value assignments are retained.
allocate(character(len=256) :: out(max(0, n)))
out = ""
end function r_character

pure function r_drop_index_real(x, k) result(out)
! Return a vector with one or more positions removed.
real(kind=dp), intent(in) :: x(:) ! source vector
integer, intent(in) :: k ! position to remove
real(kind=dp), allocatable :: out(:)
logical, allocatable :: keep(:)
integer :: n, m
n = size(x)
if (n <= 0) then
   allocate(out(0))
   return
end if
allocate(keep(n))
keep = .true.
if (k >= 1 .and. k <= n) keep(k) = .false.
m = count(keep)
allocate(out(m))
if (m > 0) out = pack(x, keep)
end function r_drop_index_real

pure function r_drop_index_int(x, k) result(out)
! Return a vector with one or more positions removed.
integer, intent(in) :: x(:) ! source vector
integer, intent(in) :: k ! position to remove
integer, allocatable :: out(:)
logical, allocatable :: keep(:)
integer :: n, m
n = size(x)
if (n <= 0) then
   allocate(out(0))
   return
end if
allocate(keep(n))
keep = .true.
if (k >= 1 .and. k <= n) keep(k) = .false.
m = count(keep)
allocate(out(m))
if (m > 0) out = pack(x, keep)
end function r_drop_index_int

pure function rev_real(x) result(out)
! Return a real vector in reverse order.
real(kind=dp), intent(in) :: x(:) ! source vector
real(kind=dp), allocatable :: out(:)
integer :: i, n
n = size(x)
allocate(out(n))
do i = 1, n
   out(i) = x(n - i + 1)
end do
end function rev_real

pure function rev_int(x) result(out)
! Return an integer vector in reverse order.
integer, intent(in) :: x(:) ! source vector
integer, allocatable :: out(:)
integer :: i, n
n = size(x)
allocate(out(n))
do i = 1, n
   out(i) = x(n - i + 1)
end do
end function rev_int

pure function r_drop_indices_real(x, drop) result(out)
! Return a vector with one or more positions removed.
real(kind=dp), intent(in) :: x(:) ! source vector
integer, intent(in) :: drop(:) ! positions to remove
real(kind=dp), allocatable :: out(:)
logical, allocatable :: keep(:)
integer :: i, n, m
n = size(x)
if (n <= 0) then
   allocate(out(0))
   return
end if
allocate(keep(n))
keep = .true.
do i = 1, size(drop)
   if (drop(i) >= 1 .and. drop(i) <= n) keep(drop(i)) = .false.
end do
m = count(keep)
allocate(out(m))
if (m > 0) out = pack(x, keep)
end function r_drop_indices_real

pure function r_drop_indices_int(x, drop) result(out)
! Return a vector with one or more positions removed.
integer, intent(in) :: x(:) ! input values
integer, intent(in) :: drop(:)
integer, allocatable :: out(:)
logical, allocatable :: keep(:)
integer :: i, n, m
n = size(x)
if (n <= 0) then
   allocate(out(0))
   return
end if
allocate(keep(n))
keep = .true.
do i = 1, size(drop)
   if (drop(i) >= 1 .and. drop(i) <= n) keep(drop(i)) = .false.
end do
m = count(keep)
allocate(out(m))
if (m > 0) out = pack(x, keep)
end function r_drop_indices_int

pure function matrix_elem(x, i, j) result(out)
real(kind=dp), intent(in) :: x(:,:)
integer, intent(in) :: i, j
real(kind=dp) :: out
out = x(i, j)
end function matrix_elem

pure function r_matrix_index_real(x, idx) result(out)
! Return R-style linear indexing of a matrix in column-major order.
real(kind=dp), intent(in) :: x(:,:) ! source matrix
integer, intent(in) :: idx(:) ! one-dimensional R subscript
real(kind=dp), allocatable :: out(:)
logical, allocatable :: keep(:)
integer :: i, k, n, m
n = size(x)
if (size(idx) == 0 .or. n <= 0) then
   allocate(out(0))
   return
end if
if (all(idx <= 0) .and. any(idx < 0)) then
   allocate(keep(n))
   keep = .true.
   do i = 1, size(idx)
      k = abs(idx(i))
      if (k >= 1 .and. k <= n) keep(k) = .false.
   end do
   m = count(keep)
   allocate(out(m))
   if (m > 0) out = pack(reshape(x, [n]), keep)
else
   m = count(idx /= 0)
   allocate(out(m))
   k = 0
   do i = 1, size(idx)
      if (idx(i) == 0) cycle
      k = k + 1
      if (idx(i) >= 1 .and. idx(i) <= n) then
         out(k) = x(mod(idx(i) - 1, size(x, 1)) + 1, ((idx(i) - 1) / size(x, 1)) + 1)
      else
         out(k) = r_na_real()
      end if
   end do
end if
end function r_matrix_index_real

pure function r_matrix_index_int(x, idx) result(out)
! Return R-style linear indexing of an integer matrix in column-major order.
integer, intent(in) :: x(:,:) ! source matrix
integer, intent(in) :: idx(:) ! one-dimensional R subscript
integer, allocatable :: out(:)
logical, allocatable :: keep(:)
integer :: i, k, n, m
n = size(x)
if (size(idx) == 0 .or. n <= 0) then
   allocate(out(0))
   return
end if
if (all(idx <= 0) .and. any(idx < 0)) then
   allocate(keep(n))
   keep = .true.
   do i = 1, size(idx)
      k = abs(idx(i))
      if (k >= 1 .and. k <= n) keep(k) = .false.
   end do
   m = count(keep)
   allocate(out(m))
   if (m > 0) out = pack(reshape(x, [n]), keep)
else
   m = count(idx /= 0)
   allocate(out(m))
   k = 0
   do i = 1, size(idx)
      if (idx(i) == 0) cycle
      k = k + 1
      if (idx(i) >= 1 .and. idx(i) <= n) then
         out(k) = x(mod(idx(i) - 1, size(x, 1)) + 1, ((idx(i) - 1) / size(x, 1)) + 1)
      else
         out(k) = -huge(0)
      end if
   end do
end if
end function r_matrix_index_int

pure function r_vector_index_real_int(x, idx) result(out)
! Return R-style integer indexing of a real vector, omitting zero indices.
real(kind=dp), intent(in) :: x(:)
integer, intent(in) :: idx(:)
real(kind=dp), allocatable :: out(:)
logical, allocatable :: keep(:)
integer :: i, k
if (all(idx <= 0) .and. any(idx < 0)) then
   allocate(keep(size(x)))
   keep = .true.
   do i = 1, size(idx)
      k = abs(idx(i))
      if (k >= 1 .and. k <= size(x)) keep(k) = .false.
   end do
   out = pack(x, keep)
   return
end if
allocate(out(count(idx /= 0)))
k = 0
do i = 1, size(idx)
   if (idx(i) == 0) cycle
   k = k + 1
   if (idx(i) >= 1 .and. idx(i) <= size(x)) then
      out(k) = x(idx(i))
   else
      out(k) = r_na_real()
   end if
end do
end function r_vector_index_real_int

pure function r_vector_index_int_int(x, idx) result(out)
! Return R-style integer indexing of an integer vector, omitting zero indices.
integer, intent(in) :: x(:)
integer, intent(in) :: idx(:)
integer, allocatable :: out(:)
logical, allocatable :: keep(:)
integer :: i, k, m
if (all(idx <= 0) .and. any(idx < 0)) then
   allocate(keep(size(x)))
   keep = .true.
   do i = 1, size(idx)
      k = abs(idx(i))
      if (k >= 1 .and. k <= size(x)) keep(k) = .false.
   end do
   out = pack(x, keep)
   return
end if
allocate(out(count(idx /= 0)))
k = 0
do i = 1, size(idx)
   if (idx(i) == 0) cycle
   k = k + 1
   if (idx(i) >= 1 .and. idx(i) <= size(x)) then
      out(k) = x(idx(i))
   else
      out(k) = -huge(0)
   end if
end do
end function r_vector_index_int_int

pure function r_vector_index_real_logical(x, mask) result(out)
! Return R-style vector indexing with a recycled logical mask.
real(kind=dp), intent(in) :: x(:)
logical, intent(in) :: mask(:)
real(kind=dp), allocatable :: out(:)
integer :: i, k, m, n
if (size(mask) == 0) then
   allocate(out(0))
   return
end if
n = max(size(x), size(mask))
m = 0
do i = 1, n
   if (mask(mod(i - 1, size(mask)) + 1)) m = m + 1
end do
allocate(out(m))
k = 0
do i = 1, n
   if (.not. mask(mod(i - 1, size(mask)) + 1)) cycle
   k = k + 1
   if (i <= size(x)) then
      out(k) = x(i)
   else
      out(k) = r_na_real()
   end if
end do
end function r_vector_index_real_logical

pure function r_vector_index_int_logical(x, mask) result(out)
! Return R-style integer-vector indexing with a recycled logical mask.
integer, intent(in) :: x(:)
logical, intent(in) :: mask(:)
integer, allocatable :: out(:)
integer :: i, k, m, n
if (size(mask) == 0) then
   allocate(out(0))
   return
end if
n = max(size(x), size(mask))
m = 0
do i = 1, n
   if (mask(mod(i - 1, size(mask)) + 1)) m = m + 1
end do
allocate(out(m))
k = 0
do i = 1, n
   if (.not. mask(mod(i - 1, size(mask)) + 1)) cycle
   k = k + 1
   if (i <= size(x)) then
      out(k) = x(i)
   else
      out(k) = -huge(0)
   end if
end do
end function r_vector_index_int_logical

pure function r_matrix_index_real_logical_vec(x, mask) result(out)
! Return R-style linear matrix indexing with a recycled logical vector.
real(kind=dp), intent(in) :: x(:,:)
logical, intent(in) :: mask(:)
real(kind=dp), allocatable :: out(:)
real(kind=dp), allocatable :: flat(:)
allocate(flat(size(x)))
flat = reshape(x, [size(x)])
out = r_vector_index_real_logical(flat, mask)
end function r_matrix_index_real_logical_vec

pure function r_matrix_index_int_logical_vec(x, mask) result(out)
! Return R-style linear integer-matrix indexing with a recycled logical vector.
integer, intent(in) :: x(:,:)
logical, intent(in) :: mask(:)
integer, allocatable :: out(:)
integer, allocatable :: flat(:)
allocate(flat(size(x)))
flat = reshape(x, [size(x)])
out = r_vector_index_int_logical(flat, mask)
end function r_matrix_index_int_logical_vec

pure function r_matrix_index_real_logical(x, mask) result(out)
! Return R-style logical linear indexing of a matrix in column-major order.
real(kind=dp), intent(in) :: x(:,:) ! source matrix
logical, intent(in) :: mask(:,:) ! same-shape selection mask
real(kind=dp), allocatable :: out(:)
integer :: m
m = count(mask)
allocate(out(m))
if (m > 0) out = pack(x, mask)
end function r_matrix_index_real_logical

pure function r_matrix_index_int_logical(x, mask) result(out)
! Return R-style logical linear indexing of an integer matrix in column-major order.
integer, intent(in) :: x(:,:) ! source matrix
logical, intent(in) :: mask(:,:) ! same-shape selection mask
integer, allocatable :: out(:)
integer :: m
m = count(mask)
allocate(out(m))
if (m > 0) out = pack(x, mask)
end function r_matrix_index_int_logical

pure function r_index_scalar_real(x, idx) result(out)
! Return one real-vector element or R's numeric NA when out of range.
real(kind=dp), intent(in) :: x(:)
integer, intent(in) :: idx
real(kind=dp) :: out
if (idx >= 1 .and. idx <= size(x)) then
   out = x(idx)
else
   out = r_na_real()
end if
end function r_index_scalar_real

pure function r_matrix_col_real(x, j) result(out)
! Return one matrix column as a vector.
real(kind=dp), intent(in) :: x(:,:) ! source matrix
integer, intent(in) :: j ! one-based column index
real(kind=dp), allocatable :: out(:)
allocate(out(size(x, 1)))
if (j >= 1 .and. j <= size(x, 2)) then
   out = x(:, j)
else
   out = 0.0_dp
end if
end function r_matrix_col_real

pure function r_matrix_col_int(x, j) result(out)
! Return one integer matrix column as a vector.
integer, intent(in) :: x(:,:) ! source matrix
integer, intent(in) :: j ! one-based column index
integer, allocatable :: out(:)
allocate(out(size(x, 1)))
if (j >= 1 .and. j <= size(x, 2)) then
   out = x(:, j)
else
   out = 0
end if
end function r_matrix_col_int

pure function r_matrix_row_real(x, i) result(out)
! Return one matrix row as a vector.
real(kind=dp), intent(in) :: x(:,:) ! source matrix
integer, intent(in) :: i ! one-based row index
real(kind=dp), allocatable :: out(:)
allocate(out(size(x, 2)))
if (i >= 1 .and. i <= size(x, 1)) then
   out = x(i, :)
else
   out = 0.0_dp
end if
end function r_matrix_row_real

pure function r_matrix_row_int(x, i) result(out)
! Return one integer matrix row as a vector.
integer, intent(in) :: x(:,:) ! source matrix
integer, intent(in) :: i ! one-based row index
integer, allocatable :: out(:)
allocate(out(size(x, 2)))
if (i >= 1 .and. i <= size(x, 1)) then
   out = x(i, :)
else
   out = 0
end if
end function r_matrix_row_int

pure function r_matrix_rows_real(x, idx) result(out)
! Return matrix rows selected by integer indices.
real(kind=dp), intent(in) :: x(:,:) ! source matrix
integer, intent(in) :: idx(:) ! one-based row indices
real(kind=dp), allocatable :: out(:,:)
integer :: i
allocate(out(size(idx), size(x, 2)))
do i = 1, size(idx)
   if (idx(i) >= 1 .and. idx(i) <= size(x, 1)) then
      out(i, :) = x(idx(i), :)
   else
      out(i, :) = r_na_real()
   end if
end do
end function r_matrix_rows_real

pure function r_matrix_rows_int(x, idx) result(out)
! Return integer matrix rows selected by integer indices.
integer, intent(in) :: x(:,:) ! source matrix
integer, intent(in) :: idx(:) ! one-based row indices
integer, allocatable :: out(:,:)
integer :: i
allocate(out(size(idx), size(x, 2)))
do i = 1, size(idx)
   if (idx(i) >= 1 .and. idx(i) <= size(x, 1)) then
      out(i, :) = x(idx(i), :)
   else
      out(i, :) = -huge(0)
   end if
end do
end function r_matrix_rows_int

pure function r_matrix_row_filter_real(x, mask) result(out)
! Return matrix rows selected by a recycled logical mask.
real(kind=dp), intent(in) :: x(:,:) ! source matrix
logical, intent(in) :: mask(:) ! row mask, recycled to nrow(x)
real(kind=dp), allocatable :: out(:,:)
integer :: i, k, nsel
nsel = 0
if (size(mask) > 0) then
   do i = 1, size(x, 1)
      if (mask(mod(i - 1, size(mask)) + 1)) nsel = nsel + 1
   end do
end if
allocate(out(nsel, size(x, 2)))
k = 0
if (size(mask) > 0) then
   do i = 1, size(x, 1)
      if (mask(mod(i - 1, size(mask)) + 1)) then
         k = k + 1
         out(k, :) = x(i, :)
      end if
   end do
end if
end function r_matrix_row_filter_real

pure function r_matrix_row_filter_int(x, mask) result(out)
! Return integer matrix rows selected by a recycled logical mask.
integer, intent(in) :: x(:,:) ! source matrix
logical, intent(in) :: mask(:) ! row mask, recycled to nrow(x)
integer, allocatable :: out(:,:)
integer :: i, k, nsel
nsel = 0
if (size(mask) > 0) then
   do i = 1, size(x, 1)
      if (mask(mod(i - 1, size(mask)) + 1)) nsel = nsel + 1
   end do
end if
allocate(out(nsel, size(x, 2)))
k = 0
if (size(mask) > 0) then
   do i = 1, size(x, 1)
      if (mask(mod(i - 1, size(mask)) + 1)) then
         k = k + 1
         out(k, :) = x(i, :)
      end if
   end do
end if
end function r_matrix_row_filter_int

pure function r_matrix_col_filter_real(x, mask) result(out)
! Return matrix columns selected by a recycled logical mask.
real(kind=dp), intent(in) :: x(:,:) ! source matrix
logical, intent(in) :: mask(:) ! column mask, recycled to ncol(x)
real(kind=dp), allocatable :: out(:,:)
integer :: j, k, nsel
nsel = 0
if (size(mask) > 0) then
   do j = 1, size(x, 2)
      if (mask(mod(j - 1, size(mask)) + 1)) nsel = nsel + 1
   end do
end if
allocate(out(size(x, 1), nsel))
k = 0
if (size(mask) > 0) then
   do j = 1, size(x, 2)
      if (mask(mod(j - 1, size(mask)) + 1)) then
         k = k + 1
         out(:, k) = x(:, j)
      end if
   end do
end if
end function r_matrix_col_filter_real

pure function r_matrix_col_filter_int(x, mask) result(out)
! Return integer matrix columns selected by a recycled logical mask.
integer, intent(in) :: x(:,:) ! source matrix
logical, intent(in) :: mask(:) ! column mask, recycled to ncol(x)
integer, allocatable :: out(:,:)
integer :: j, k, nsel
nsel = 0
if (size(mask) > 0) then
   do j = 1, size(x, 2)
      if (mask(mod(j - 1, size(mask)) + 1)) nsel = nsel + 1
   end do
end if
allocate(out(size(x, 1), nsel))
k = 0
if (size(mask) > 0) then
   do j = 1, size(x, 2)
      if (mask(mod(j - 1, size(mask)) + 1)) then
         k = k + 1
         out(:, k) = x(:, j)
      end if
   end do
end if
end function r_matrix_col_filter_int

pure function r_index_real(x, idx) result(out)
! Return R-style vector indexing with real indices so NaN indices become real NA.
integer, intent(in) :: x(:) ! source integer vector
real(kind=dp), intent(in) :: idx(:) ! one-based real/NA indices
real(kind=dp), allocatable :: out(:)
logical, allocatable :: keep(:)
logical :: has_negative, has_nonfinite, has_positive
integer :: i, k, m
has_negative = .false.
has_nonfinite = .false.
has_positive = .false.
do i = 1, size(idx)
   if (.not. ieee_is_finite(idx(i))) then
      has_nonfinite = .true.
   else
      k = int(idx(i))
      if (k < 0) has_negative = .true.
      if (k > 0) has_positive = .true.
   end if
end do
if (has_negative .and. .not. has_positive .and. .not. has_nonfinite) then
   allocate(keep(size(x)))
   keep = .true.
   do i = 1, size(idx)
      k = abs(int(idx(i)))
      if (k >= 1 .and. k <= size(x)) keep(k) = .false.
   end do
   out = pack(real(x, kind=dp), keep)
   return
end if
m = 0
do i = 1, size(idx)
   if (.not. ieee_is_finite(idx(i))) then
      m = m + 1
   else if (int(idx(i)) /= 0) then
      m = m + 1
   end if
end do
allocate(out(m))
m = 0
do i = 1, size(idx)
   if (ieee_is_finite(idx(i))) then
      if (int(idx(i)) == 0) cycle
   end if
   m = m + 1
   if (.not. ieee_is_finite(idx(i))) then
      out(m) = r_na_real()
   else
      k = int(idx(i))
      if (k >= 1 .and. k <= size(x)) then
         out(m) = real(x(k), kind=dp)
      else
         out(m) = r_na_real()
      end if
   end if
end do
end function r_index_real


subroutine set_print_int_like(flag)
! Enable/disable integer-like rendering for real matrix printing.
logical, intent(in) :: flag ! logical flag
print_int_like_default = flag
end subroutine set_print_int_like

subroutine set_print_int_like_tol(tol)
! Set tolerance used for integer-like real rendering.
real(kind=dp), intent(in) :: tol ! convergence tolerance
if (tol > 0.0_dp) print_int_like_tol = tol
end subroutine set_print_int_like_tol

subroutine set_recycle_warn(flag)
! Enable/disable warnings for non-multiple recycling lengths.
logical, intent(in) :: flag ! logical flag
recycle_warn_default = flag
end subroutine set_recycle_warn

subroutine set_recycle_stop(flag)
! Enable/disable error stop for non-multiple recycling lengths.
logical, intent(in) :: flag ! logical flag
recycle_stop_default = flag
end subroutine set_recycle_stop

subroutine set_nlm_method(method)
! Select nlm optimizer: "legacy" preserves the historical gradient-descent helper.
character(len=*), intent(in) :: method
character(len=:), allocatable :: m
m = trim(adjustl(method))
select case (m)
case ("newton", "Newton", "NEWTON")
   nlm_method_default = 2
case default
   nlm_method_default = 1
end select
end subroutine set_nlm_method

subroutine set_seed_int(seed)
! Set Fortran RNG seed deterministically from a single integer.
integer, intent(in) :: seed ! random seed
#ifdef XR2F_USE_R_RNG
call xr2f_r_set_seed(int(seed, kind=c_int))
#else
integer :: n, i
integer, allocatable :: put(:)
integer(kind=int64) :: s, m
call random_seed(size=n)
allocate(put(n))
s = int(abs(seed), kind=int64)
if (s == 0_int64) s = 104729_int64
m = int(huge(0), kind=int64) - 1_int64
do i = 1, n
   put(i) = int(modulo(s + 104729_int64 * int(i, kind=int64), m) + 1_int64)
end do
call random_seed(put=put)
deallocate(put)
#endif
end subroutine set_seed_int

function kmeans_vec(x, centers, nstart, iter_max, algorithm) result(out)
! Minimal 1D k-means helper: returns centers and 1-based cluster ids.
real(kind=dp), intent(in) :: x(:) ! input values
integer, intent(in) :: centers ! cluster centers
integer, intent(in), optional :: nstart
integer, intent(in), optional :: iter_max
character(len=*), intent(in), optional :: algorithm
type(kmeans_result_t) :: out
real(kind=dp), allocatable :: c(:), c_new(:), sums(:), best_withinss(:), withinss(:), best_centers(:)
integer, allocatable :: cnt(:), cl(:), cl_best(:), best_size(:)
integer, allocatable :: size_tot(:)
integer, allocatable :: order_idx(:), remap(:)
integer :: i, j, k, n, it, jbest, nstart_loc, iter_max_loc, start, idx, empty_j, far_i, far_j
integer :: t, best_iter
real(kind=dp) :: xmin, xmax, scale, d, dbest, u, best_score, score, far_d, xmean, totss
real(kind=dp) :: delta, best_delta
logical :: changed, use_macqueen, use_hartigan
character(len=32) :: alg
n = size(x)
k = max(1, centers)
nstart_loc = 1
if (present(nstart)) nstart_loc = max(1, nstart)
iter_max_loc = 50
if (present(iter_max)) iter_max_loc = max(1, iter_max)
alg = "Hartigan-Wong"
if (present(algorithm)) alg = trim(adjustl(algorithm))
use_macqueen = alg == "MacQueen" .or. alg == "macqueen" .or. alg == "MACQUEEN"
use_hartigan = alg == "Hartigan-Wong" .or. alg == "hartigan-wong" .or. alg == "HARTIGAN-WONG"
allocate(c(k), c_new(k), sums(k), cnt(k), cl(n), cl_best(n))
allocate(withinss(k), best_withinss(k), best_centers(k), best_size(k), size_tot(k))
if (n <= 0) then
   c = 0.0_dp
   sums = 0.0_dp
   cnt = 0
   withinss = 0.0_dp
   out%centers = reshape(c, [k, 1])
   out%cluster = cl
   out%size = cnt
   out%withinss = withinss
   return
end if
xmin = minval(x)
xmax = maxval(x)
best_score = huge(1.0_dp)
best_withinss = 0.0_dp
best_size = 0
cl_best = 1
best_centers = 0.0_dp
best_iter = 0
if (k > 0) cl_best(1) = 1
xmean = sum(x) / real(n, kind=dp)
totss = sum((x - xmean)**2)
do start = 1, nstart_loc
   if (k == 1) then
      if (nstart_loc > 1 .and. start > 1) then
   u = runif1()
         idx = 1 + int(real(n, kind=dp) * u)
         if (idx > n) idx = n
         c(1) = x(idx)
      else
         c(1) = sum(x) / real(n, kind=dp)
      end if
   else
      if (nstart_loc > 1 .and. start > 1) then
         do j = 1, k
   u = runif1()
            idx = 1 + int(real(n, kind=dp) * u)
            if (idx > n) idx = n
            c(j) = x(idx)
         end do
      else
         c(1) = sum(x) / real(n, kind=dp)
         do j = 2, k
            far_i = 1
            far_d = -1.0_dp
            do i = 1, n
               dbest = abs(x(i) - c(1))
               do t = 2, j - 1
                  dbest = min(dbest, abs(x(i) - c(t)))
               end do
               if (dbest > far_d) then
                  far_d = dbest
                  far_i = i
               end if
            end do
            c(j) = x(far_i)
         end do
      end if
   end if
   cl = 0
   do it = 1, iter_max_loc
      changed = .false.
      if (use_macqueen) then
         cnt = 0
         do i = 1, n
            jbest = 1
            dbest = abs(x(i) - c(1))
            do j = 2, k
               d = abs(x(i) - c(j))
               if (d < dbest) then
                  dbest = d
                  jbest = j
               end if
            end do
            if (cl(i) /= jbest) changed = .true.
            cl(i) = jbest
            cnt(jbest) = cnt(jbest) + 1
            c(jbest) = c(jbest) + (x(i) - c(jbest)) / real(cnt(jbest), kind=dp)
         end do
         if (.not. changed) exit
         cycle
      end if
      do i = 1, n
         jbest = 1
         dbest = abs(x(i) - c(1))
         do j = 2, k
            d = abs(x(i) - c(j))
            if (d < dbest) then
               dbest = d
               jbest = j
            end if
         end do
         if (cl(i) /= jbest) changed = .true.
         cl(i) = jbest
      end do
      sums = 0.0_dp
      cnt = 0
      do i = 1, n
         j = cl(i)
         sums(j) = sums(j) + x(i)
         cnt(j) = cnt(j) + 1
      end do
      do empty_j = 1, k
         if (cnt(empty_j) == 0) then
            far_i = 1
            far_j = cl(1)
            far_d = -1.0_dp
            do i = 1, n
               j = cl(i)
               if (cnt(j) > 1) then
                  d = (x(i) - sums(j) / real(cnt(j), kind=dp))**2
                  if (d > far_d) then
                     far_d = d
                     far_i = i
                     far_j = j
                  end if
               end if
            end do
            if (far_d >= 0.0_dp) then
               sums(far_j) = sums(far_j) - x(far_i)
               cnt(far_j) = cnt(far_j) - 1
               cl(far_i) = empty_j
               sums(empty_j) = x(far_i)
               cnt(empty_j) = 1
               changed = .true.
            end if
         end if
      end do
      c_new = c
      do j = 1, k
         if (cnt(j) > 0) c_new(j) = sums(j) / real(cnt(j), kind=dp)
      end do
      scale = maxval(abs(c_new - c))
      c = c_new
      if ((.not. changed) .or. scale <= 1.0e-12_dp * max(1.0_dp, maxval(abs(c)))) exit
   end do

   if (use_hartigan .and. k > 1) then
      do t = 1, iter_max_loc
         changed = .false.
         sums = 0.0_dp
         cnt = 0
         do i = 1, n
            sums(cl(i)) = sums(cl(i)) + x(i)
            cnt(cl(i)) = cnt(cl(i)) + 1
         end do
         do i = 1, n
            far_j = cl(i)
            if (cnt(far_j) <= 1) cycle
            best_delta = 0.0_dp
            jbest = far_j
            do j = 1, k
               if (j == far_j) cycle
               delta = real(cnt(j), kind=dp) / real(cnt(j) + 1, kind=dp) * (x(i) - c(j))**2 &
                  - real(cnt(far_j), kind=dp) / real(cnt(far_j) - 1, kind=dp) * (x(i) - c(far_j))**2
               if (delta < best_delta) then
                  best_delta = delta
                  jbest = j
               end if
            end do
            if (jbest /= far_j) then
               sums(far_j) = sums(far_j) - x(i)
               cnt(far_j) = cnt(far_j) - 1
               sums(jbest) = sums(jbest) + x(i)
               cnt(jbest) = cnt(jbest) + 1
               c(far_j) = sums(far_j) / real(cnt(far_j), kind=dp)
               c(jbest) = sums(jbest) / real(cnt(jbest), kind=dp)
               cl(i) = jbest
               changed = .true.
            end if
         end do
         if (.not. changed) exit
      end do
      it = min(iter_max_loc, it + t)
   end if

   size_tot = 0
   withinss = 0.0_dp
   do i = 1, n
      j = cl(i)
      size_tot(j) = size_tot(j) + 1
      d = x(i) - c(j)
      withinss(j) = withinss(j) + d * d
   end do
   score = sum(withinss)
   if (score < best_score) then
      best_score = score
      best_centers = c
      cl_best = cl
      best_withinss = withinss
      best_size = size_tot
      best_iter = it
   end if
end do
out%centers = reshape(best_centers, [k, 1])
out%cluster = cl_best
out%size = best_size
out%withinss = best_withinss
out%totss = totss
out%tot_withinss = sum(best_withinss)
out%betweenss = out%totss - out%tot_withinss
out%iter = best_iter

allocate(order_idx(k), remap(k))
do i = 1, k
   order_idx(i) = i
end do
if (k > 1) then
   do i = 1, k - 1
      do j = i + 1, k
         if (best_size(j) > best_size(order_idx(i))) then
            t = order_idx(i)
            order_idx(i) = order_idx(j)
            order_idx(j) = t
         end if
      end do
   end do
end if
if (all(best_size == 0)) order_idx = [(i, i=1,k)]

withinss = best_withinss
sums = best_centers
cnt = best_size
do i = 1, k
   best_size(i) = cnt(order_idx(i))
   best_withinss(i) = withinss(order_idx(i))
   best_centers(i) = sums(order_idx(i))
   remap(order_idx(i)) = i
end do

do i = 1, n
   cl_best(i) = remap(cl_best(i))
end do

out%centers = reshape(best_centers, [k, 1])
out%cluster = cl_best
out%size = best_size
out%withinss = best_withinss
out%totss = totss
out%tot_withinss = sum(best_withinss)
out%betweenss = out%totss - out%tot_withinss
out%iter = best_iter
end function kmeans_vec

function kmeans_mat(x, centers, nstart, iter_max, algorithm) result(out)
! Minimal row-wise k-means helper for matrix observations.
real(kind=dp), intent(in) :: x(:,:) ! input data matrix
integer, intent(in) :: centers ! cluster centers
integer, intent(in), optional :: nstart
integer, intent(in), optional :: iter_max
character(len=*), intent(in), optional :: algorithm
type(kmeans_result_t) :: out
real(kind=dp), allocatable :: c(:,:), c_new(:,:), sums(:,:), best_centers(:,:), best_withinss(:), withinss(:)
integer, allocatable :: cnt(:), cl(:), cl_best(:), best_size(:), size_tot(:)
integer, allocatable :: order_idx(:), remap(:), ci(:)
integer :: i, j, k, n, p, it, jbest, nstart_loc, iter_max_loc, start, empty_j, far_i, far_j
integer :: idx, best_iter
integer :: t
real(kind=dp) :: d, dbest, shift, u, best_score, score, far_d, totss, delta, best_delta
real(kind=dp), allocatable :: xmean(:)
logical :: changed, use_macqueen, use_hartigan
character(len=32) :: alg
n = size(x, 1)
p = size(x, 2)
k = max(1, centers)
nstart_loc = 1
if (present(nstart)) nstart_loc = max(1, nstart)
iter_max_loc = 50
if (present(iter_max)) iter_max_loc = max(1, iter_max)
alg = "Hartigan-Wong"
if (present(algorithm)) alg = trim(adjustl(algorithm))
use_macqueen = alg == "MacQueen" .or. alg == "macqueen" .or. alg == "MACQUEEN"
use_hartigan = alg == "Hartigan-Wong" .or. alg == "hartigan-wong" .or. alg == "HARTIGAN-WONG"
allocate(c(k, p), c_new(k, p), sums(k, p), best_centers(k, p), withinss(k), best_withinss(k))
allocate(cnt(k), size_tot(k), cl(n), cl_best(n), best_size(k), out%size(k), out%withinss(k))
if (n <= 0 .or. p <= 0) then
   allocate(out%centers(k, p), source=0.0_dp)
   cl = 0
   out%cluster = cl
   out%size = 0
   out%withinss = 0.0_dp
   out%totss = 0.0_dp
   out%tot_withinss = 0.0_dp
   out%betweenss = 0.0_dp
   out%iter = 0
   return
end if
best_score = huge(1.0_dp)
best_withinss = 0.0_dp
best_size = 0
best_centers = 0.0_dp
cl_best = 1
best_iter = 0
allocate(xmean(p))
xmean = sum(x, dim=1) / real(n, kind=dp)
totss = 0.0_dp
do i = 1, n
   totss = totss + sum((x(i, :) - xmean)**2)
end do
do start = 1, nstart_loc
   if (k > 1) then
      if (nstart_loc > 1 .and. start > 1) then
         do j = 1, k
   u = runif1()
            idx = 1 + int(real(n, kind=dp) * u)
            if (idx > n) idx = n
            c(j, :) = x(idx, :)
         end do
      else
         c(1, :) = sum(x, dim=1) / real(n, kind=dp)
         do j = 2, k
            far_i = 1
            far_d = -1.0_dp
            do i = 1, n
               dbest = sum((x(i, :) - c(1, :))**2)
               do t = 2, j - 1
                  dbest = min(dbest, sum((x(i, :) - c(t, :))**2))
               end do
               if (dbest > far_d) then
                  far_d = dbest
                  far_i = i
               end if
            end do
            c(j, :) = x(far_i, :)
         end do
      end if
   else
      c(1, :) = sum(x, dim=1) / real(n, kind=dp)
   end if

   cl = 0
   do it = 1, iter_max_loc
      changed = .false.
      if (use_macqueen) then
         cnt = 0
         do i = 1, n
            jbest = 1
            dbest = sum((x(i, :) - c(1, :))**2)
            do j = 2, k
               d = sum((x(i, :) - c(j, :))**2)
               if (d < dbest) then
                  dbest = d
                  jbest = j
               end if
            end do
            if (cl(i) /= jbest) changed = .true.
            cl(i) = jbest
            cnt(jbest) = cnt(jbest) + 1
            c(jbest, :) = c(jbest, :) + (x(i, :) - c(jbest, :)) / real(cnt(jbest), kind=dp)
         end do
         if (.not. changed) exit
         cycle
      end if
      do i = 1, n
         jbest = 1
         dbest = sum((x(i, :) - c(1, :))**2)
         do j = 2, k
            d = sum((x(i, :) - c(j, :))**2)
            if (d < dbest) then
               dbest = d
               jbest = j
            end if
         end do
         if (cl(i) /= jbest) changed = .true.
         cl(i) = jbest
      end do
      sums = 0.0_dp
      cnt = 0
      do i = 1, n
         j = cl(i)
         sums(j, :) = sums(j, :) + x(i, :)
         cnt(j) = cnt(j) + 1
      end do
      do empty_j = 1, k
         if (cnt(empty_j) == 0) then
            far_i = 1
            far_j = cl(1)
            far_d = -1.0_dp
            do i = 1, n
               j = cl(i)
               if (cnt(j) > 1) then
                  d = sum((x(i, :) - sums(j, :) / real(cnt(j), kind=dp))**2)
                  if (d > far_d) then
                     far_d = d
                     far_i = i
                     far_j = j
                  end if
               end if
            end do
            if (far_d >= 0.0_dp) then
               sums(far_j, :) = sums(far_j, :) - x(far_i, :)
               cnt(far_j) = cnt(far_j) - 1
               cl(far_i) = empty_j
               sums(empty_j, :) = x(far_i, :)
               cnt(empty_j) = 1
               changed = .true.
            end if
         end if
      end do
      c_new = c
      do j = 1, k
         if (cnt(j) > 0) c_new(j, :) = sums(j, :) / real(cnt(j), kind=dp)
      end do
      shift = maxval(abs(c_new - c))
      c = c_new
      if ((.not. changed) .or. shift <= 1.0e-12_dp * max(1.0_dp, maxval(abs(c)))) exit
   end do

   if (use_hartigan .and. k > 1) then
      do t = 1, iter_max_loc
         changed = .false.
         sums = 0.0_dp
         cnt = 0
         do i = 1, n
            sums(cl(i), :) = sums(cl(i), :) + x(i, :)
            cnt(cl(i)) = cnt(cl(i)) + 1
         end do
         do i = 1, n
            far_j = cl(i)
            if (cnt(far_j) <= 1) cycle
            best_delta = 0.0_dp
            jbest = far_j
            do j = 1, k
               if (j == far_j) cycle
               delta = real(cnt(j), kind=dp) / real(cnt(j) + 1, kind=dp) * sum((x(i, :) - c(j, :))**2) &
                  - real(cnt(far_j), kind=dp) / real(cnt(far_j) - 1, kind=dp) * sum((x(i, :) - c(far_j, :))**2)
               if (delta < best_delta) then
                  best_delta = delta
                  jbest = j
               end if
            end do
            if (jbest /= far_j) then
               sums(far_j, :) = sums(far_j, :) - x(i, :)
               cnt(far_j) = cnt(far_j) - 1
               sums(jbest, :) = sums(jbest, :) + x(i, :)
               cnt(jbest) = cnt(jbest) + 1
               c(far_j, :) = sums(far_j, :) / real(cnt(far_j), kind=dp)
               c(jbest, :) = sums(jbest, :) / real(cnt(jbest), kind=dp)
               cl(i) = jbest
               changed = .true.
            end if
         end do
         if (.not. changed) exit
      end do
      it = min(iter_max_loc, it + t)
   end if

   size_tot = 0
   withinss = 0.0_dp
   do i = 1, n
      j = cl(i)
      size_tot(j) = size_tot(j) + 1
      d = sum((x(i, :) - c(j, :))**2)
      withinss(j) = withinss(j) + d
   end do
   score = sum(withinss)
   if (score < best_score) then
      best_score = score
      best_centers = c
      cl_best = cl
      best_size = size_tot
      best_withinss = withinss
      best_iter = it
   end if
end do
out%centers = best_centers
out%cluster = cl_best
out%size = best_size
out%withinss = best_withinss
out%totss = totss
out%tot_withinss = sum(best_withinss)
out%betweenss = out%totss - out%tot_withinss
out%iter = best_iter

allocate(order_idx(k), remap(k), ci(n))
do i = 1, k
   order_idx(i) = i
end do
if (k > 1) then
   do i = 1, k - 1
      do j = i + 1, k
         if (best_size(j) > best_size(order_idx(i))) then
            t = order_idx(i)
            order_idx(i) = order_idx(j)
            order_idx(j) = t
         end if
      end do
   end do
end if
if (all(best_size == 0)) order_idx = [(i, i=1,k)]

c_new = best_centers
cnt = best_size
withinss = best_withinss
do i = 1, k
   best_size(i) = cnt(order_idx(i))
   best_withinss(i) = withinss(order_idx(i))
   best_centers(i, :) = c_new(order_idx(i), :)
   remap(order_idx(i)) = i
end do
do i = 1, n
   ci(i) = remap(cl_best(i))
end do
cl_best = ci

out%centers = best_centers
out%cluster = cl_best
out%size = best_size
out%withinss = best_withinss
out%totss = totss
out%tot_withinss = sum(best_withinss)
out%betweenss = out%totss - out%tot_withinss
out%iter = best_iter
end function kmeans_mat

pure function dist_mat(x, method) result(out)
! Compute pairwise distance matrix for observations in rows.
! Valid method values: "euclidean" (default), "manhattan", "maximum", "canberra".
real(kind=dp), intent(in) :: x(:,:) ! observations by rows, variables by columns
character(len=*), intent(in), optional :: method ! distance metric name
real(kind=dp), allocatable :: out(:,:)
integer :: i, j, k, n, p
real(kind=dp) :: acc
character(len=32) :: meth
n = size(x, 1)
p = size(x, 2)
allocate(out(max(0, n), max(0, n)))
if (n <= 0) return
if (present(method)) then
   meth = trim(adjustl(method))
else
   meth = "euclidean"
end if
out = 0.0_dp
do i = 1, n - 1
   out(i, i) = 0.0_dp
   do j = i + 1, n
      acc = 0.0_dp
      if (meth == "manhattan") then
         do k = 1, p
            acc = acc + abs(x(i, k) - x(j, k))
         end do
      else if (meth == "maximum") then
         acc = abs(x(i, 1) - x(j, 1))
         do k = 2, p
            acc = max(acc, abs(x(i, k) - x(j, k)))
         end do
      else if (meth == "canberra") then
         do k = 1, p
            if (abs(x(i, k)) + abs(x(j, k)) /= 0.0_dp) then
               acc = acc + abs(x(i, k) - x(j, k)) / (abs(x(i, k)) + abs(x(j, k)))
            else
               acc = acc + 0.0_dp
            end if
         end do
      else
         do k = 1, p
            acc = acc + (x(i, k) - x(j, k))**2
         end do
         acc = sqrt(acc)
      end if
      out(i, j) = acc
      out(j, i) = acc
   end do
end do
if (n >= 1) out(n, n) = 0.0_dp
end function dist_mat

pure function hclust_complete(d, method, labels) result(out)
! Minimal hierarchical clustering helper on a distance matrix.
! Supported method values: complete, single, average, mcquitty, centroid, median, ward.D, ward.D2.
real(kind=dp), intent(in) :: d(:,:) ! square distance matrix
character(len=*), intent(in), optional :: method ! linkage method name
character(len=*), intent(in), optional :: labels(:) ! optional observation labels
type(hclust_result_t) :: out
logical, allocatable :: alive(:)
real(kind=dp), allocatable :: cdist(:,:)
integer, allocatable :: cluster_rep(:), node_size(:)
integer, allocatable :: stack(:)
integer :: n, step, i, j, k
integer :: node_count, new_node, stack_top, entry, order_pos
real(kind=dp) :: best_d, cand_d, da, db, dab, sa, sb, sk, denom
integer :: best_a, best_b
character(len=16) :: meth
integer :: method_code
n = size(d, 1)
if (n <= 1) then
   allocate(out%merge(0, 2))
   allocate(out%height(0))
   allocate(out%order(max(0, n)))
   allocate(character(len=32) :: out%labels(max(0, n)))
   if (n == 1) then
      out%labels(1) = "1"
      if (present(labels)) then
         if (size(labels) >= 1) then
            deallocate(out%labels)
            allocate(character(len=max(1, len_trim(labels(1)))) :: out%labels(1))
            out%labels(1) = trim(labels(1))
         end if
      end if
      out%order = [1]
   end if
   out%method = 1
   return
end if
if (size(d, 2) /= n) then
   allocate(out%merge(0, 2))
   allocate(out%height(0))
   allocate(out%order(0))
   allocate(character(len=32) :: out%labels(0))
   out%method = 1
   return
end if
if (present(method)) then
   meth = trim(adjustl(method))
else
   meth = "complete"
end if
select case (meth)
case ("single", "SINGLE")
   method_code = 2
case ("average", "AVERAGE", "UPGMA")
   method_code = 3
case ("mcquitty", "MCQUITTY", "WPGMA")
   method_code = 4
case ("centroid", "CENTROID", "UPGMC")
   method_code = 5
case ("median", "MEDIAN", "WPGMC")
   method_code = 6
case ("ward.D", "ward.d", "WARD.D")
   method_code = 7
case ("ward.D2", "ward.d2", "WARD.D2")
   method_code = 8
case default
   meth = "complete"
   method_code = 1
end select
allocate(out%merge(max(0, n - 1), 2))
allocate(out%height(max(0, n - 1)))
allocate(out%order(n))
allocate(character(len=32) :: out%labels(n))
do i = 1, n
   out%labels(i) = int_to_string(i)
end do
if (present(labels)) then
   if (size(labels) >= n) then
      deallocate(out%labels)
      allocate(character(len=max(1, maxval(len_trim(labels(1:n))))) :: out%labels(n))
      do i = 1, n
         out%labels(i) = trim(labels(i))
      end do
   end if
end if
out%order = [(i, i = 1, n)]
out%method = method_code
allocate(alive(2 * n - 1))
allocate(cdist(2 * n - 1, 2 * n - 1))
allocate(cluster_rep(max(0, n - 1)))
allocate(node_size(2 * n - 1))
allocate(stack(max(1, 2 * n)))
cdist = 0.0_dp
alive = .false.
node_size = 0
cdist(1:n, 1:n) = d
alive(1:n) = .true.
node_size(1:n) = 1
do step = 1, n - 1
   best_d = huge(1.0_dp)
   best_a = -1
   best_b = -1
   node_count = n + step - 1
   do i = 1, node_count
      if (.not. alive(i)) cycle
      do j = i + 1, node_count
         if (.not. alive(j)) cycle
         if (cdist(i, j) < best_d) then
            best_d = cdist(i, j)
            best_a = i
            best_b = j
         end if
      end do
   end do
   if (best_a <= 0 .or. best_b <= 0) exit
   if (best_a <= n) then
      out%merge(step, 1) = -best_a
   else
      out%merge(step, 1) = best_a - n
   end if
   if (best_b <= n) then
      out%merge(step, 2) = -best_b
   else
      out%merge(step, 2) = best_b - n
   end if
   out%height(step) = best_d
   new_node = n + step
   alive(best_a) = .false.
   alive(best_b) = .false.
   alive(new_node) = .true.
   node_size(new_node) = node_size(best_a) + node_size(best_b)
   cluster_rep(step) = new_node
   do k = 1, node_count
      if (.not. alive(k) .or. k == new_node) cycle
      da = cdist(best_a, k)
      db = cdist(best_b, k)
      dab = cdist(best_a, best_b)
      sa = real(node_size(best_a), kind=dp)
      sb = real(node_size(best_b), kind=dp)
      sk = real(node_size(k), kind=dp)
      if (method_code == 2) then
         cand_d = min(da, db)
      else if (method_code == 3) then
         cand_d = (sa * da + sb * db) / (sa + sb)
      else if (method_code == 4) then
         cand_d = 0.5_dp * (da + db)
      else if (method_code == 5) then
         denom = sa + sb
         cand_d = (sa * da + sb * db) / denom - (sa * sb * dab) / (denom * denom)
         cand_d = max(0.0_dp, cand_d)
      else if (method_code == 6) then
         cand_d = 0.5_dp * (da + db) - 0.25_dp * dab
         cand_d = max(0.0_dp, cand_d)
      else if (method_code == 7) then
         denom = sa + sb + sk
         cand_d = ((sa + sk) * da + (sb + sk) * db - sk * dab) / denom
         cand_d = max(0.0_dp, cand_d)
      else if (method_code == 8) then
         denom = sa + sb + sk
         cand_d = ((sa + sk) * da * da + (sb + sk) * db * db - sk * dab * dab) / denom
         cand_d = sqrt(max(0.0_dp, cand_d))
      else
         cand_d = max(cdist(best_a, k), cdist(best_b, k))
      end if
      cdist(new_node, k) = cand_d
      cdist(k, new_node) = cand_d
   end do
   cdist(new_node, new_node) = 0.0_dp
   cdist(:, best_a) = 0.0_dp
   cdist(best_a, :) = 0.0_dp
   cdist(:, best_b) = 0.0_dp
   cdist(best_b, :) = 0.0_dp
end do
out%order = 0
order_pos = 0
stack_top = 1
stack(stack_top) = n - 1
do while (stack_top > 0)
   entry = stack(stack_top)
   stack_top = stack_top - 1
   if (entry < 0) then
      order_pos = order_pos + 1
      if (order_pos <= n) out%order(order_pos) = -entry
   else if (entry >= 1 .and. entry <= n - 1) then
      stack_top = stack_top + 1
      stack(stack_top) = out%merge(entry, 2)
      stack_top = stack_top + 1
      stack(stack_top) = out%merge(entry, 1)
   end if
end do
if (any(out%order == 0)) out%order = [(i, i = 1, n)]
end function hclust_complete

pure function cutree_f90(fit, k, h) result(group)
! Cut dendrogram at a target number of groups or height.
! If k and h are absent, return one group; provided k is clamped to the range 1:n.
type(hclust_result_t), intent(in) :: fit ! hierarchical clustering result
integer, intent(in), optional :: k ! requested number of groups
real(kind=dp), intent(in), optional :: h ! requested cut height
integer, allocatable :: group(:)
integer, allocatable :: rep_map(:), parent(:), merge_rep(:)
integer :: n, nmerge, target_groups, nmerge_apply
integer :: i, a, b, root_a, root_b, ra, rb, r, next_group
nmerge = size(fit%merge, 1)
n = max(0, nmerge) + 1
if (n > 0 .and. size(fit%labels) >= n) n = size(fit%labels)
allocate(group(max(0, n)))
if (n <= 0) then
   return
end if
if (size(fit%merge, 2) /= 2) then
   do i = 1, n
      group(i) = i
   end do
   return
end if
if (nmerge <= 0) then
   do i = 1, n
      group(i) = i
   end do
   return
end if
if (present(h)) then
   nmerge_apply = 0
   do i = 1, nmerge
      if (fit%height(i) <= h) then
         nmerge_apply = i
      else
         exit
      end if
   end do
   target_groups = max(1, n - nmerge_apply)
else if (present(k)) then
   target_groups = max(1, min(k, n))
   nmerge_apply = min(nmerge, n - target_groups)
else
   target_groups = 1
   nmerge_apply = min(nmerge, n - target_groups)
end if
if (target_groups >= n) then
   do i = 1, n
      group(i) = i
   end do
   return
end if
allocate(parent(n))
allocate(rep_map(n))
allocate(merge_rep(max(0, nmerge)))
do i = 1, n
   parent(i) = i
   rep_map(i) = 0
   group(i) = 0
end do
merge_rep = 0

do i = 1, nmerge_apply
   a = fit%merge(i, 1)
   b = fit%merge(i, 2)
   if (a < 0) then
      root_a = -a
   else if (a >= 1 .and. a <= nmerge) then
      root_a = merge_rep(a)
   else
      root_a = a
   end if
   if (b < 0) then
      root_b = -b
   else if (b >= 1 .and. b <= nmerge) then
      root_b = merge_rep(b)
   else
      root_b = b
   end if
   root_a = max(1, min(n, root_a))
   root_b = max(1, min(n, root_b))
   ra = root_a
   do while (parent(ra) /= ra)
      parent(ra) = parent(parent(ra))
      ra = parent(ra)
   end do
   rb = root_b
   do while (parent(rb) /= rb)
      parent(rb) = parent(parent(rb))
      rb = parent(rb)
   end do
   if (ra /= rb) parent(rb) = ra
   if (i <= nmerge) merge_rep(i) = ra
end do
next_group = 0
do i = 1, n
   r = i
   do while (parent(r) /= r)
      parent(r) = parent(parent(r))
      r = parent(r)
   end do
   if (rep_map(r) == 0) then
      next_group = next_group + 1
      rep_map(r) = next_group
   end if
   group(i) = rep_map(r)
end do
end function cutree_f90

pure function max_col(x, ties_method) result(idx)
! Return 1-based column index of row-wise maxima (ties -> first).
real(kind=dp), intent(in) :: x(:,:) ! input matrix
character(len=*), intent(in), optional :: ties_method
integer, allocatable :: idx(:)
integer :: i, j, n, k, jbest
real(kind=dp) :: vbest
n = size(x, 1)
k = size(x, 2)
allocate(idx(n))
if (present(ties_method)) then
   ! Only "first" tie handling is supported in this subset.
end if
if (k <= 0) then
   idx = 1
   return
end if
do i = 1, n
   jbest = 1
   vbest = x(i, 1)
   do j = 2, k
      if (x(i, j) > vbest) then
         vbest = x(i, j)
         jbest = j
      end if
   end do
   idx(i) = jbest
end do
end function max_col

pure function tabulate_int(x, nbins) result(out)
! Count occurrences of integer labels 1..nbins.
integer, intent(in) :: x(:) ! input vector
integer, intent(in) :: nbins
integer, allocatable :: out(:)
integer :: i, b
allocate(out(max(0, nbins)))
if (size(out) > 0) out = 0
do i = 1, size(x)
   b = x(i)
   if (b >= 1 .and. b <= size(out)) out(b) = out(b) + 1
end do
end function tabulate_int

pure function tabulate_real(x, nbins) result(out)
! Count occurrences after integer-coding real labels.
real(kind=dp), intent(in) :: x(:) ! input vector
integer, intent(in) :: nbins
integer, allocatable :: out(:)
integer :: i, b
allocate(out(max(0, nbins)))
if (size(out) > 0) out = 0
if (size(out) <= 0) return
do i = 1, size(x)
   if (x(i) /= x(i)) cycle
   if (.not. ieee_is_finite(x(i))) cycle
   if (x(i) < 1.0_dp .or. x(i) >= real(size(out) + 1, kind=dp)) cycle
   b = int(x(i))
   out(b) = out(b) + 1
end do
end function tabulate_real

pure function table_char(x) result(out)
character(len=*), intent(in) :: x(:)
type(table_char_t) :: out
character(len=max(1, len(x))) :: labels(max(1, size(x))), key
integer :: i, j, n
n = 0
do i = 1, size(x)
   if (n == 0 .or. .not. any(labels(1:n) == x(i))) then
      n = n + 1
      labels(n) = x(i)
   end if
end do
do i = 2, n
   key = labels(i)
   j = i - 1
   do while (j >= 1)
      if (labels(j) <= key) exit
      labels(j + 1) = labels(j)
      j = j - 1
   end do
   labels(j + 1) = key
end do
allocate(character(len=max(1, len(x))) :: out%gene(n))
allocate(out%Freq(n))
if (n > 0) out%gene = labels(1:n)
do i = 1, n
   out%Freq(i) = count(x == out%gene(i))
end do
end function table_char

pure function table2_int(x, y, nx, ny) result(out)
! Count paired integer labels into an nx-by-ny contingency table.
integer, intent(in) :: x(:) ! input vector
integer, intent(in) :: y(:) ! response values
integer, intent(in) :: nx, ny
integer, allocatable :: out(:,:)
integer :: i, a, b, n
allocate(out(max(0, nx), max(0, ny)))
if (size(out) > 0) out = 0
n = min(size(x), size(y))
do i = 1, n
   a = x(i)
   b = y(i)
   if (a >= 1 .and. a <= size(out, 1) .and. b >= 1 .and. b <= size(out, 2)) then
      out(a, b) = out(a, b) + 1
   end if
end do
end function table2_int

pure function prop_table_int_vec(x, margin) result(out)
! Convert integer counts to proportions.
integer, intent(in) :: x(:) ! input vector
integer, intent(in), optional :: margin
real(kind=dp), allocatable :: out(:)
integer :: s
if (present(margin)) continue
allocate(out(size(x)))
s = sum(x)
if (s > 0) then
   out = real(x, kind=dp) / real(s, kind=dp)
else
   out = ieee_value(0.0_dp, ieee_quiet_nan)
end if
end function prop_table_int_vec

pure function prop_table_int_mat(x, margin) result(out)
! Convert integer contingency tables to overall, row, or column proportions.
integer, intent(in) :: x(:,:) ! input matrix
integer, intent(in), optional :: margin
real(kind=dp), allocatable :: out(:,:)
integer :: i, j, s
allocate(out(size(x, 1), size(x, 2)))
if (present(margin)) then
   if (margin == 1) then
      do i = 1, size(x, 1)
         s = sum(x(i, :))
         if (s > 0) then
            out(i, :) = real(x(i, :), kind=dp) / real(s, kind=dp)
         else
            out(i, :) = ieee_value(0.0_dp, ieee_quiet_nan)
         end if
      end do
      return
   else if (margin == 2) then
      do j = 1, size(x, 2)
         s = sum(x(:, j))
         if (s > 0) then
            out(:, j) = real(x(:, j), kind=dp) / real(s, kind=dp)
         else
            out(:, j) = ieee_value(0.0_dp, ieee_quiet_nan)
         end if
      end do
      return
   end if
end if
s = sum(x)
if (s > 0) then
   out = real(x, kind=dp) / real(s, kind=dp)
else
   out = ieee_value(0.0_dp, ieee_quiet_nan)
end if
end function prop_table_int_mat

pure function ave_real_char(x, g, fun) result(out)
! Apply a simple group aggregate and return one value per input observation.
real(kind=dp), intent(in) :: x(:)
character(len=*), intent(in) :: g(:)
character(len=*), intent(in) :: fun
real(kind=dp), allocatable :: out(:)
integer :: i, j, n, cnt
real(kind=dp) :: total, vmin, vmax
n = min(size(x), size(g))
allocate(out(size(x)))
out = ieee_value(0.0_dp, ieee_quiet_nan)
do i = 1, n
   total = 0.0_dp
   vmin = huge(0.0_dp)
   vmax = -huge(0.0_dp)
   cnt = 0
   do j = 1, n
      if (trim(g(j)) == trim(g(i))) then
         total = total + x(j)
         vmin = min(vmin, x(j))
         vmax = max(vmax, x(j))
         cnt = cnt + 1
      end if
   end do
   select case (trim(fun))
   case ("mean")
      if (cnt > 0) out(i) = total / real(cnt, kind=dp)
   case ("sum")
      out(i) = total
   case ("length")
      out(i) = real(cnt, kind=dp)
   case ("min")
      if (cnt > 0) out(i) = vmin
   case ("max")
      if (cnt > 0) out(i) = vmax
   case default
      out(i) = ieee_value(0.0_dp, ieee_quiet_nan)
   end select
end do
end function ave_real_char

pure function ave_real_int(x, g, fun) result(out)
! Apply a simple group aggregate and return one value per input observation.
real(kind=dp), intent(in) :: x(:)
integer, intent(in) :: g(:)
character(len=*), intent(in) :: fun
real(kind=dp), allocatable :: out(:)
integer :: i, j, n, cnt
real(kind=dp) :: total, vmin, vmax
n = min(size(x), size(g))
allocate(out(size(x)))
out = ieee_value(0.0_dp, ieee_quiet_nan)
do i = 1, n
   total = 0.0_dp
   vmin = huge(0.0_dp)
   vmax = -huge(0.0_dp)
   cnt = 0
   do j = 1, n
      if (g(j) == g(i)) then
         total = total + x(j)
         vmin = min(vmin, x(j))
         vmax = max(vmax, x(j))
         cnt = cnt + 1
      end if
   end do
   select case (trim(fun))
   case ("mean")
      if (cnt > 0) out(i) = total / real(cnt, kind=dp)
   case ("sum")
      out(i) = total
   case ("length")
      out(i) = real(cnt, kind=dp)
   case ("min")
      if (cnt > 0) out(i) = vmin
   case ("max")
      if (cnt > 0) out(i) = vmax
   case default
      out(i) = ieee_value(0.0_dp, ieee_quiet_nan)
   end select
end do
end function ave_real_int

pure function ave_group_key_char_char(a, b) result(out)
! Combine two grouping vectors into one character key vector.
character(len=*), intent(in) :: a(:), b(:)
character(len=:), allocatable :: out(:)
integer :: i, n, key_len
n = min(size(a), size(b))
key_len = len(a) + len(b) + 1
allocate(character(len=max(1, key_len)) :: out(n))
do i = 1, n
   out(i) = trim(a(i)) // achar(31) // trim(b(i))
end do
end function ave_group_key_char_char

pure function ave_group_key_char_int(a, b) result(out)
! Combine character and integer grouping vectors into one character key vector.
character(len=*), intent(in) :: a(:)
integer, intent(in) :: b(:)
character(len=:), allocatable :: out(:)
integer :: i, n, key_len
n = min(size(a), size(b))
key_len = len(a) + 32
allocate(character(len=max(1, key_len)) :: out(n))
do i = 1, n
   out(i) = trim(a(i)) // achar(31) // trim(int_to_string(b(i)))
end do
end function ave_group_key_char_int

pure function ave_group_key_int_char(a, b) result(out)
! Combine integer and character grouping vectors into one character key vector.
integer, intent(in) :: a(:)
character(len=*), intent(in) :: b(:)
character(len=:), allocatable :: out(:)
integer :: i, n, key_len
n = min(size(a), size(b))
key_len = len(b) + 32
allocate(character(len=max(1, key_len)) :: out(n))
do i = 1, n
   out(i) = trim(int_to_string(a(i))) // achar(31) // trim(b(i))
end do
end function ave_group_key_int_char

pure function ave_group_key_int_int(a, b) result(out)
! Combine two integer grouping vectors into one character key vector.
integer, intent(in) :: a(:), b(:)
character(len=:), allocatable :: out(:)
integer :: i, n
n = min(size(a), size(b))
allocate(character(len=65) :: out(n))
do i = 1, n
   out(i) = trim(int_to_string(a(i))) // achar(31) // trim(int_to_string(b(i)))
end do
end function ave_group_key_int_int

pure function aggregate_real_char(x, g, group_name, value_name, fun) result(out)
! Aggregate a numeric vector by one character grouping vector.
real(kind=dp), intent(in) :: x(:)
character(len=*), intent(in) :: g(:)
character(len=*), intent(in) :: group_name, value_name, fun
type(aggregate_result_t) :: out
integer :: i, j, n, ngrp, cnt
logical :: seen
real(kind=dp) :: total, vmin, vmax
n = min(size(x), size(g))
ngrp = 0
do i = 1, n
   seen = .false.
   do j = 1, i - 1
      if (trim(g(j)) == trim(g(i))) then
         seen = .true.
         exit
      end if
   end do
   if (.not. seen) ngrp = ngrp + 1
end do
out%group_name = trim(group_name)
out%value_name = trim(value_name)
allocate(character(len=max(1, len(g))) :: out%labels(ngrp))
allocate(out%values(ngrp))
ngrp = 0
do i = 1, n
   seen = .false.
   do j = 1, ngrp
      if (trim(out%labels(j)) == trim(g(i))) then
         seen = .true.
         exit
      end if
   end do
   if (seen) cycle
   ngrp = ngrp + 1
   out%labels(ngrp) = trim(g(i))
   total = 0.0_dp
   vmin = huge(0.0_dp)
   vmax = -huge(0.0_dp)
   cnt = 0
   do j = 1, n
      if (trim(g(j)) == trim(g(i))) then
         total = total + x(j)
         vmin = min(vmin, x(j))
         vmax = max(vmax, x(j))
         cnt = cnt + 1
      end if
   end do
   select case (trim(fun))
   case ("mean")
      if (cnt > 0) out%values(ngrp) = total / real(cnt, kind=dp)
   case ("sum")
      out%values(ngrp) = total
   case ("length")
      out%values(ngrp) = real(cnt, kind=dp)
   case ("min")
      out%values(ngrp) = vmin
   case ("max")
      out%values(ngrp) = vmax
   case default
      out%values(ngrp) = ieee_value(0.0_dp, ieee_quiet_nan)
   end select
end do
end function aggregate_real_char

pure function aggregate_real_int(x, g, group_name, value_name, fun) result(out)
! Aggregate a numeric vector by one integer grouping vector.
real(kind=dp), intent(in) :: x(:)
integer, intent(in) :: g(:)
character(len=*), intent(in) :: group_name, value_name, fun
type(aggregate_result_t) :: out
integer :: i, j, n, ngrp, cnt
logical :: seen
real(kind=dp) :: total, vmin, vmax
n = min(size(x), size(g))
ngrp = 0
do i = 1, n
   seen = .false.
   do j = 1, i - 1
      if (g(j) == g(i)) then
         seen = .true.
         exit
      end if
   end do
   if (.not. seen) ngrp = ngrp + 1
end do
out%group_name = trim(group_name)
out%value_name = trim(value_name)
allocate(character(len=32) :: out%labels(ngrp))
allocate(out%values(ngrp))
ngrp = 0
do i = 1, n
   seen = .false.
   do j = 1, ngrp
      if (trim(out%labels(j)) == trim(int_to_string(g(i)))) then
         seen = .true.
         exit
      end if
   end do
   if (seen) cycle
   ngrp = ngrp + 1
   out%labels(ngrp) = trim(int_to_string(g(i)))
   total = 0.0_dp
   vmin = huge(0.0_dp)
   vmax = -huge(0.0_dp)
   cnt = 0
   do j = 1, n
      if (g(j) == g(i)) then
         total = total + x(j)
         vmin = min(vmin, x(j))
         vmax = max(vmax, x(j))
         cnt = cnt + 1
      end if
   end do
   select case (trim(fun))
   case ("mean")
      if (cnt > 0) out%values(ngrp) = total / real(cnt, kind=dp)
   case ("sum")
      out%values(ngrp) = total
   case ("length")
      out%values(ngrp) = real(cnt, kind=dp)
   case ("min")
      out%values(ngrp) = vmin
   case ("max")
      out%values(ngrp) = vmax
   case default
      out%values(ngrp) = ieee_value(0.0_dp, ieee_quiet_nan)
   end select
end do
end function aggregate_real_int

pure function by_real_vec_char(x, g, fun) result(out)
! Apply a scalar aggregate to a vector by one character grouping vector.
real(kind=dp), intent(in) :: x(:)
character(len=*), intent(in) :: g(:)
character(len=*), intent(in) :: fun
real(kind=dp), allocatable :: out(:)
type(aggregate_result_t) :: tmp
tmp = aggregate_real_char(x, g, "", "", fun)
out = tmp%values
end function by_real_vec_char

pure function by_real_vec_int(x, g, fun) result(out)
! Apply a scalar aggregate to a vector by one integer grouping vector.
real(kind=dp), intent(in) :: x(:)
integer, intent(in) :: g(:)
character(len=*), intent(in) :: fun
real(kind=dp), allocatable :: out(:)
type(aggregate_result_t) :: tmp
tmp = aggregate_real_int(x, g, "", "", fun)
out = tmp%values
end function by_real_vec_int

pure function by_real_mat_char(x, g, fun) result(out)
! Apply colMeans/colSums to row groups of a matrix.
real(kind=dp), intent(in) :: x(:,:)
character(len=*), intent(in) :: g(:)
character(len=*), intent(in) :: fun
type(by_matrix_result_t) :: out
integer :: i, j, k, n, ngrp, cnt
logical :: seen
n = min(size(x, 1), size(g))
ngrp = 0
do i = 1, n
   seen = .false.
   do j = 1, i - 1
      if (trim(g(j)) == trim(g(i))) then
         seen = .true.
         exit
      end if
   end do
   if (.not. seen) ngrp = ngrp + 1
end do
allocate(character(len=max(1, len(g))) :: out%labels(ngrp))
allocate(out%values(ngrp, size(x, 2)))
ngrp = 0
do i = 1, n
   seen = .false.
   do j = 1, ngrp
      if (trim(out%labels(j)) == trim(g(i))) then
         seen = .true.
         exit
      end if
   end do
   if (seen) cycle
   ngrp = ngrp + 1
   out%labels(ngrp) = trim(g(i))
   out%values(ngrp, :) = 0.0_dp
   cnt = 0
   do j = 1, n
      if (trim(g(j)) == trim(g(i))) then
         cnt = cnt + 1
         do k = 1, size(x, 2)
            out%values(ngrp, k) = out%values(ngrp, k) + x(j, k)
         end do
      end if
   end do
   if (trim(fun) == "colmeans" .and. cnt > 0) out%values(ngrp, :) = out%values(ngrp, :) / real(cnt, kind=dp)
end do
end function by_real_mat_char

pure function by_real_mat_int(x, g, fun) result(out)
! Apply colMeans/colSums to row groups of a matrix.
real(kind=dp), intent(in) :: x(:,:)
integer, intent(in) :: g(:)
character(len=*), intent(in) :: fun
type(by_matrix_result_t) :: out
integer :: i, j, k, n, ngrp, cnt
logical :: seen
n = min(size(x, 1), size(g))
ngrp = 0
do i = 1, n
   seen = .false.
   do j = 1, i - 1
      if (g(j) == g(i)) then
         seen = .true.
         exit
      end if
   end do
   if (.not. seen) ngrp = ngrp + 1
end do
allocate(character(len=32) :: out%labels(ngrp))
allocate(out%values(ngrp, size(x, 2)))
ngrp = 0
do i = 1, n
   seen = .false.
   do j = 1, ngrp
      if (trim(out%labels(j)) == trim(int_to_string(g(i)))) then
         seen = .true.
         exit
      end if
   end do
   if (seen) cycle
   ngrp = ngrp + 1
   out%labels(ngrp) = trim(int_to_string(g(i)))
   out%values(ngrp, :) = 0.0_dp
   cnt = 0
   do j = 1, n
      if (g(j) == g(i)) then
         cnt = cnt + 1
         do k = 1, size(x, 2)
            out%values(ngrp, k) = out%values(ngrp, k) + x(j, k)
         end do
      end if
   end do
   if (trim(fun) == "colmeans" .and. cnt > 0) out%values(ngrp, :) = out%values(ngrp, :) / real(cnt, kind=dp)
end do
end function by_real_mat_int

pure function nested_matrix_list_len(x) result(n)
! Count non-padding matrix slices in a ragged nested list lowered to rank 3.
real(kind=dp), intent(in) :: x(:,:,:)
integer :: n
integer :: j
n = 0
do j = 1, size(x, 3)
   if (any(ieee_is_finite(x(:,:,j)))) n = j
end do
end function nested_matrix_list_len

pure function match_int(x, table) result(out)
! Return first 1-based match positions, or a sentinel for NA.
integer, intent(in) :: x(:) ! values to match
integer, intent(in) :: table(:) ! lookup table
integer, allocatable :: out(:)
integer :: i, j
allocate(out(size(x)))
out = -huge(0)
do i = 1, size(x)
   do j = 1, size(table)
      if (x(i) == table(j)) then
         out(i) = j
         exit
      end if
   end do
end do
end function match_int

pure function match_real(x, table) result(out)
! Return first 1-based match positions, or a sentinel for NA.
real(kind=dp), intent(in) :: x(:) ! values to match
real(kind=dp), intent(in) :: table(:) ! lookup table
integer, allocatable :: out(:)
integer :: i, j
allocate(out(size(x)))
out = -huge(0)
do i = 1, size(x)
   do j = 1, size(table)
      if (real_values_equal(x(i), table(j))) then
         out(i) = j
         exit
      end if
   end do
end do
end function match_real

pure function match_char(x, table) result(out)
! Return first 1-based match positions, or a sentinel for NA.
character(len=*), intent(in) :: x(:) ! strings to match
character(len=*), intent(in) :: table(:) ! lookup table
integer, allocatable :: out(:)
integer :: i, j
allocate(out(size(x)))
out = -huge(0)
do i = 1, size(x)
   do j = 1, size(table)
      if (trim(x(i)) == trim(table(j))) then
         out(i) = j
         exit
      end if
   end do
end do
end function match_char

pure function r_in_int(x, table) result(out)
! Test membership for int values.
integer, intent(in) :: x(:) ! values to test
integer, intent(in) :: table(:) ! candidate set
logical, allocatable :: out(:)
integer :: i
allocate(out(size(x)))
do i = 1, size(x)
   out(i) = has_int_value(table, size(table), x(i))
end do
end function r_in_int

pure function r_in_real(x, table) result(out)
! Test membership for real values.
real(kind=dp), intent(in) :: x(:) ! values to test
real(kind=dp), intent(in) :: table(:) ! candidate set
logical, allocatable :: out(:)
integer :: i
allocate(out(size(x)))
do i = 1, size(x)
   out(i) = has_real_value(table, size(table), x(i))
end do
end function r_in_real

pure function r_in_int_real(x, table) result(out)
! Test membership for int real values.
integer, intent(in) :: x(:) ! integer values to test
real(kind=dp), intent(in) :: table(:) ! real candidate set
logical, allocatable :: out(:)
integer :: i
allocate(out(size(x)))
do i = 1, size(x)
   out(i) = r_in_int_scalar_real(x(i), table)
end do
end function r_in_int_real

pure function r_in_real_int(x, table) result(out)
! Test membership for real int values.
real(kind=dp), intent(in) :: x(:) ! real values to test
integer, intent(in) :: table(:) ! integer candidate set
logical, allocatable :: out(:)
integer :: i
allocate(out(size(x)))
do i = 1, size(x)
   out(i) = r_in_real_scalar_int(x(i), table)
end do
end function r_in_real_int

pure function r_in_int_scalar(x, table) result(out)
! Test membership for int scalar values.
integer, intent(in) :: x ! value to test
integer, intent(in) :: table(:) ! candidate set
logical :: out
out = has_int_value(table, size(table), x)
end function r_in_int_scalar

pure function r_in_real_scalar(x, table) result(out)
! Test membership for real scalar values.
real(kind=dp), intent(in) :: x ! value to test
real(kind=dp), intent(in) :: table(:) ! candidate set
logical :: out
out = has_real_value(table, size(table), x)
end function r_in_real_scalar

pure function r_in_int_scalar_real(x, table) result(out)
! Test membership for int scalar real values.
integer, intent(in) :: x ! integer value to test
real(kind=dp), intent(in) :: table(:) ! real candidate set
logical :: out
if (x == -huge(0)) then
   out = has_real_value(table, size(table), r_na_real())
else
   out = has_real_value(table, size(table), real(x, kind=dp))
end if
end function r_in_int_scalar_real

pure function r_in_real_scalar_int(x, table) result(out)
! Test membership for real scalar int values.
real(kind=dp), intent(in) :: x ! real value to test
integer, intent(in) :: table(:) ! integer candidate set
logical :: out
integer :: value
if (r_is_na_payload(x)) then
   out = has_int_value(table, size(table), -huge(0))
else if (x /= x .or. .not. ieee_is_finite(x)) then
   out = .false.
else if (x < real(-huge(0) + 1, kind=dp) .or. x > real(huge(0), kind=dp)) then
   out = .false.
else
   value = nint(x)
   out = value /= -huge(0) .and. x == real(value, kind=dp)
   if (out) out = has_int_value(table, size(table), value)
end if
end function r_in_real_scalar_int

pure function r_in_char(x, table) result(out)
! Test membership for char values.
character(len=*), intent(in) :: x(:) ! strings to test
character(len=*), intent(in) :: table(:) ! candidate set
logical, allocatable :: out(:)
integer :: i
allocate(out(size(x)))
do i = 1, size(x)
   out(i) = has_char_value(table, size(table), x(i))
end do
end function r_in_char

pure function r_in_logical(x, table) result(out)
! Test membership for logical values.
logical, intent(in) :: x(:) ! values to test
logical, intent(in) :: table(:) ! candidate set
logical, allocatable :: out(:)
integer :: i
allocate(out(size(x)))
do i = 1, size(x)
   out(i) = has_logical_value(table, size(table), x(i))
end do
end function r_in_logical

pure logical function has_int_value(x, n, value) result(out)
! Test whether a vector contains the requested int value.
integer, intent(in) :: x(:) ! candidate values
integer, intent(in) :: n ! active candidate count
integer, intent(in) :: value ! value to find
integer :: i
out = .false.
do i = 1, min(n, size(x))
   if (x(i) == value) then
      out = .true.
      return
   end if
end do
end function has_int_value

pure logical function has_real_value(x, n, value) result(out)
! Test whether a vector contains the requested real value.
real(kind=dp), intent(in) :: x(:) ! candidate values
integer, intent(in) :: n ! active candidate count
real(kind=dp), intent(in) :: value ! value to find
integer :: i
out = .false.
do i = 1, min(n, size(x))
   if (real_values_equal(x(i), value)) then
      out = .true.
      return
   end if
end do
end function has_real_value

pure elemental logical function real_values_equal(a, b) result(out)
! R equality for set membership: NA and NaN are distinct repeatable values.
real(kind=dp), intent(in) :: a, b
if (r_is_na_payload(a) .or. r_is_na_payload(b)) then
   out = r_is_na_payload(a) .and. r_is_na_payload(b)
else if (a /= a .or. b /= b) then
   out = (a /= a) .and. (b /= b)
else
   out = a == b
end if
end function real_values_equal

pure logical function has_char_value(x, n, value) result(out)
! Test whether a vector contains the requested char value.
character(len=*), intent(in) :: x(:) ! candidate strings
integer, intent(in) :: n ! active candidate count
character(len=*), intent(in) :: value ! string to find
integer :: i
out = .false.
do i = 1, min(n, size(x))
   if (x(i) == value) then
      out = .true.
      return
   end if
end do
end function has_char_value

pure logical function has_logical_value(x, n, value) result(out)
! Test whether a vector contains the requested logical value.
logical, intent(in) :: x(:) ! candidate values
integer, intent(in) :: n ! active candidate count
logical, intent(in) :: value ! value to find
integer :: i
out = .false.
do i = 1, min(n, size(x))
   if (x(i) .eqv. value) then
      out = .true.
      return
   end if
end do
end function has_logical_value

pure function unique_int(x) result(out)
! Return unique int values preserving first occurrence order.
integer, intent(in) :: x(:)
integer, allocatable :: out(:)
integer, allocatable :: tmp(:)
integer :: i, n
allocate(tmp(size(x)))
n = 0
do i = 1, size(x)
   if (.not. has_int_value(tmp, n, x(i))) then
      n = n + 1
      tmp(n) = x(i)
   end if
end do
allocate(out(n))
if (n > 0) out = tmp(1:n)
end function unique_int

pure function unique_real(x) result(out)
! Return unique real values preserving first occurrence order.
real(kind=dp), intent(in) :: x(:)
real(kind=dp), allocatable :: out(:)
real(kind=dp), allocatable :: tmp(:)
integer :: i, n
allocate(tmp(size(x)))
n = 0
do i = 1, size(x)
   if (.not. has_real_value(tmp, n, x(i))) then
      n = n + 1
      tmp(n) = x(i)
   end if
end do
allocate(out(n))
if (n > 0) out = tmp(1:n)
end function unique_real

pure function unique_char(x) result(out)
! Return unique char values preserving first occurrence order.
character(len=*), intent(in) :: x(:)
character(len=:), allocatable :: out(:)
character(len=:), allocatable :: tmp(:)
integer :: i, n, lch
lch = max(1, len(x))
allocate(character(len=lch) :: tmp(size(x)))
n = 0
do i = 1, size(x)
   if (.not. has_char_value(tmp, n, x(i))) then
      n = n + 1
      tmp(n) = x(i)
   end if
end do
allocate(character(len=lch) :: out(n))
if (n > 0) out = tmp(1:n)
end function unique_char

pure function unique_logical(x) result(out)
! Return unique logical values preserving first occurrence order.
logical, intent(in) :: x(:)
logical, allocatable :: out(:)
logical, allocatable :: tmp(:)
integer :: i, n
allocate(tmp(size(x)))
n = 0
do i = 1, size(x)
   if (.not. has_logical_value(tmp, n, x(i))) then
      n = n + 1
      tmp(n) = x(i)
   end if
end do
allocate(out(n))
if (n > 0) out = tmp(1:n)
end function unique_logical

pure function duplicated_int(x, fromLast) result(out)
! Return duplicate flags for int values.
integer, intent(in) :: x(:) ! input vector
logical, intent(in), optional :: fromLast
logical, allocatable :: out(:)
integer, allocatable :: seen(:)
integer :: i, n
logical :: rev
allocate(out(size(x)), seen(size(x)))
out = .false.
n = 0
rev = .false.
if (present(fromLast)) rev = fromLast
if (rev) then
   do i = size(x), 1, -1
      out(i) = has_int_value(seen, n, x(i))
      if (.not. out(i)) then
         n = n + 1
         seen(n) = x(i)
      end if
   end do
else
   do i = 1, size(x)
      out(i) = has_int_value(seen, n, x(i))
      if (.not. out(i)) then
         n = n + 1
         seen(n) = x(i)
      end if
   end do
end if
end function duplicated_int

pure function duplicated_real(x, fromLast) result(out)
! Return duplicate flags for real values.
real(kind=dp), intent(in) :: x(:) ! input vector
logical, intent(in), optional :: fromLast
logical, allocatable :: out(:)
real(kind=dp), allocatable :: seen(:)
integer :: i, n
logical :: rev
allocate(out(size(x)), seen(size(x)))
out = .false.
n = 0
rev = .false.
if (present(fromLast)) rev = fromLast
if (rev) then
   do i = size(x), 1, -1
      out(i) = has_real_value(seen, n, x(i))
      if (.not. out(i)) then
         n = n + 1
         seen(n) = x(i)
      end if
   end do
else
   do i = 1, size(x)
      out(i) = has_real_value(seen, n, x(i))
      if (.not. out(i)) then
         n = n + 1
         seen(n) = x(i)
      end if
   end do
end if
end function duplicated_real

pure function duplicated_char(x, fromLast) result(out)
! Return duplicate flags for char values.
character(len=*), intent(in) :: x(:) ! input vector
logical, intent(in), optional :: fromLast
logical, allocatable :: out(:)
character(len=:), allocatable :: seen(:)
integer :: i, n, lch
logical :: rev
lch = max(1, len(x))
allocate(out(size(x)))
allocate(character(len=lch) :: seen(size(x)))
out = .false.
n = 0
rev = .false.
if (present(fromLast)) rev = fromLast
if (rev) then
   do i = size(x), 1, -1
      out(i) = has_char_value(seen, n, x(i))
      if (.not. out(i)) then
         n = n + 1
         seen(n) = x(i)
      end if
   end do
else
   do i = 1, size(x)
      out(i) = has_char_value(seen, n, x(i))
      if (.not. out(i)) then
         n = n + 1
         seen(n) = x(i)
      end if
   end do
end if
end function duplicated_char

pure function duplicated_logical(x, fromLast) result(out)
! Return duplicate flags for logical values.
logical, intent(in) :: x(:) ! input vector
logical, intent(in), optional :: fromLast
logical, allocatable :: out(:)
logical, allocatable :: seen(:)
integer :: i, n
logical :: rev
allocate(out(size(x)), seen(size(x)))
out = .false.
n = 0
rev = .false.
if (present(fromLast)) rev = fromLast
if (rev) then
   do i = size(x), 1, -1
      out(i) = has_logical_value(seen, n, x(i))
      if (.not. out(i)) then
         n = n + 1
         seen(n) = x(i)
      end if
   end do
else
   do i = 1, size(x)
      out(i) = has_logical_value(seen, n, x(i))
      if (.not. out(i)) then
         n = n + 1
         seen(n) = x(i)
      end if
   end do
end if
end function duplicated_logical

pure function anyDuplicated_int(x, fromLast) result(out)
! Return the first duplicated position for int values.
integer, intent(in) :: x(:) ! input vector
logical, intent(in), optional :: fromLast
integer :: out, i
logical, allocatable :: dup(:)
dup = duplicated_int(x, fromLast)
out = 0
do i = 1, size(dup)
   if (dup(i)) then
      out = i
      return
   end if
end do
end function anyDuplicated_int

pure function anyDuplicated_real(x, fromLast) result(out)
! Return the first duplicated position for real values.
real(kind=dp), intent(in) :: x(:) ! input vector
logical, intent(in), optional :: fromLast
integer :: out, i
logical, allocatable :: dup(:)
dup = duplicated_real(x, fromLast)
out = 0
do i = 1, size(dup)
   if (dup(i)) then
      out = i
      return
   end if
end do
end function anyDuplicated_real

pure function anyDuplicated_char(x, fromLast) result(out)
! Return the first duplicated position for char values.
character(len=*), intent(in) :: x(:) ! input vector
logical, intent(in), optional :: fromLast
integer :: out, i
logical, allocatable :: dup(:)
dup = duplicated_char(x, fromLast)
out = 0
do i = 1, size(dup)
   if (dup(i)) then
      out = i
      return
   end if
end do
end function anyDuplicated_char

pure function anyDuplicated_logical(x, fromLast) result(out)
! Return the first duplicated position for logical values.
logical, intent(in) :: x(:) ! input vector
logical, intent(in), optional :: fromLast
integer :: out, i
logical, allocatable :: dup(:)
dup = duplicated_logical(x, fromLast)
out = 0
do i = 1, size(dup)
   if (dup(i)) then
      out = i
      return
   end if
end do
end function anyDuplicated_logical

pure function union_int(x, y) result(out)
! Return the set union for int vectors.
integer, intent(in) :: x(:), y(:)
integer, allocatable :: out(:)
integer, allocatable :: tmp(:)
integer :: i, n
allocate(tmp(size(x) + size(y)))
n = 0
do i = 1, size(x)
   if (.not. has_int_value(tmp, n, x(i))) then
      n = n + 1
      tmp(n) = x(i)
   end if
end do
do i = 1, size(y)
   if (.not. has_int_value(tmp, n, y(i))) then
      n = n + 1
      tmp(n) = y(i)
   end if
end do
allocate(out(n))
if (n > 0) out = tmp(1:n)
end function union_int

pure function intersect_int(x, y) result(out)
! Return the set intersection for int vectors.
integer, intent(in) :: x(:), y(:)
integer, allocatable :: out(:)
integer, allocatable :: tmp(:)
integer :: i, n
allocate(tmp(size(x)))
n = 0
do i = 1, size(x)
   if ((.not. has_int_value(tmp, n, x(i))) .and. has_int_value(y, size(y), x(i))) then
      n = n + 1
      tmp(n) = x(i)
   end if
end do
allocate(out(n))
if (n > 0) out = tmp(1:n)
end function intersect_int

pure function setdiff_int(x, y) result(out)
! Return values in the first int vector but not the second.
integer, intent(in) :: x(:), y(:)
integer, allocatable :: out(:)
integer, allocatable :: tmp(:)
integer :: i, n
allocate(tmp(size(x)))
n = 0
do i = 1, size(x)
   if ((.not. has_int_value(tmp, n, x(i))) .and. (.not. has_int_value(y, size(y), x(i)))) then
      n = n + 1
      tmp(n) = x(i)
   end if
end do
allocate(out(n))
if (n > 0) out = tmp(1:n)
end function setdiff_int

pure function setequal_int(x, y) result(out)
! Test set equality for int vectors.
integer, intent(in) :: x(:), y(:)
logical :: out
out = size(setdiff_int(x, y)) == 0 .and. size(setdiff_int(y, x)) == 0
end function setequal_int

pure function union_real(x, y) result(out)
! Return the set union for real vectors.
real(kind=dp), intent(in) :: x(:), y(:)
real(kind=dp), allocatable :: out(:)
real(kind=dp), allocatable :: tmp(:)
integer :: i, n
allocate(tmp(size(x) + size(y)))
n = 0
do i = 1, size(x)
   if (.not. has_real_value(tmp, n, x(i))) then
      n = n + 1
      tmp(n) = x(i)
   end if
end do
do i = 1, size(y)
   if (.not. has_real_value(tmp, n, y(i))) then
      n = n + 1
      tmp(n) = y(i)
   end if
end do
allocate(out(n))
if (n > 0) out = tmp(1:n)
end function union_real

pure function intersect_real(x, y) result(out)
! Return the set intersection for real vectors.
real(kind=dp), intent(in) :: x(:), y(:)
real(kind=dp), allocatable :: out(:)
real(kind=dp), allocatable :: tmp(:)
integer :: i, n
allocate(tmp(size(x)))
n = 0
do i = 1, size(x)
   if ((.not. has_real_value(tmp, n, x(i))) .and. has_real_value(y, size(y), x(i))) then
      n = n + 1
      tmp(n) = x(i)
   end if
end do
allocate(out(n))
if (n > 0) out = tmp(1:n)
end function intersect_real

pure function setdiff_real(x, y) result(out)
! Return values in the first real vector but not the second.
real(kind=dp), intent(in) :: x(:), y(:)
real(kind=dp), allocatable :: out(:)
real(kind=dp), allocatable :: tmp(:)
integer :: i, n
allocate(tmp(size(x)))
n = 0
do i = 1, size(x)
   if ((.not. has_real_value(tmp, n, x(i))) .and. (.not. has_real_value(y, size(y), x(i)))) then
      n = n + 1
      tmp(n) = x(i)
   end if
end do
allocate(out(n))
if (n > 0) out = tmp(1:n)
end function setdiff_real

pure function setequal_real(x, y) result(out)
! Test set equality for real vectors.
real(kind=dp), intent(in) :: x(:), y(:)
logical :: out
out = size(setdiff_real(x, y)) == 0 .and. size(setdiff_real(y, x)) == 0
end function setequal_real

pure function union_char(x, y) result(out)
! Return the set union for char vectors.
character(len=*), intent(in) :: x(:), y(:)
character(len=:), allocatable :: out(:)
character(len=:), allocatable :: tmp(:)
integer :: i, n, lch
lch = max(1, len(x), len(y))
allocate(character(len=lch) :: tmp(size(x) + size(y)))
n = 0
do i = 1, size(x)
   if (.not. has_char_value(tmp, n, x(i))) then
      n = n + 1
      tmp(n) = x(i)
   end if
end do
do i = 1, size(y)
   if (.not. has_char_value(tmp, n, y(i))) then
      n = n + 1
      tmp(n) = y(i)
   end if
end do
allocate(character(len=lch) :: out(n))
if (n > 0) out = tmp(1:n)
end function union_char

pure function intersect_char(x, y) result(out)
! Return the set intersection for char vectors.
character(len=*), intent(in) :: x(:), y(:)
character(len=:), allocatable :: out(:)
character(len=:), allocatable :: tmp(:)
integer :: i, n, lch
lch = max(1, len(x), len(y))
allocate(character(len=lch) :: tmp(size(x)))
n = 0
do i = 1, size(x)
   if ((.not. has_char_value(tmp, n, x(i))) .and. has_char_value(y, size(y), x(i))) then
      n = n + 1
      tmp(n) = x(i)
   end if
end do
allocate(character(len=lch) :: out(n))
if (n > 0) out = tmp(1:n)
end function intersect_char

pure function setdiff_char(x, y) result(out)
! Return values in the first char vector but not the second.
character(len=*), intent(in) :: x(:), y(:)
character(len=:), allocatable :: out(:)
character(len=:), allocatable :: tmp(:)
integer :: i, n, lch
lch = max(1, len(x), len(y))
allocate(character(len=lch) :: tmp(size(x)))
n = 0
do i = 1, size(x)
   if ((.not. has_char_value(tmp, n, x(i))) .and. (.not. has_char_value(y, size(y), x(i)))) then
      n = n + 1
      tmp(n) = x(i)
   end if
end do
allocate(character(len=lch) :: out(n))
if (n > 0) out = tmp(1:n)
end function setdiff_char

pure function setequal_char(x, y) result(out)
! Test set equality for char vectors.
character(len=*), intent(in) :: x(:), y(:)
logical :: out
out = size(setdiff_char(x, y)) == 0 .and. size(setdiff_char(y, x)) == 0
end function setequal_char

pure function r_format_vec(x, digits, sep) result(out)
! Format a real vector like paste(sprintf("%.<digits>f", x), collapse=sep).
real(kind=dp), intent(in) :: x(:) ! input vector
integer, intent(in) :: digits ! number of digits
character(len=*), intent(in), optional :: sep
character(len=:), allocatable :: out
character(len=64) :: fmt, buf
integer :: i, d, nbuf
d = max(0, min(30, digits))
write(fmt, '("(f0.", i0, ")")') d
out = ""
do i = 1, size(x)
   write(buf, fmt) x(i)
   buf = adjustl(buf)
   if (d == 0) then
      nbuf = len_trim(buf)
      if (nbuf > 0) then
         if (buf(nbuf:nbuf) == ".") buf(nbuf:nbuf) = " "
      end if
   end if
   if (i > 1) then
      if (present(sep)) then
         out = out // sep
      else
         out = out // " "
      end if
   end if
   out = out // trim(buf)
end do
end function r_format_vec

subroutine maybe_warn_recycle(op, na, nb)
! Warn/stop when any vector recycling occurs (lengths differ).
character(len=*), intent(in) :: op ! input string
integer, intent(in) :: na ! integer argument
integer, intent(in) :: nb ! integer argument
integer :: nmax, nmin
if ((.not. recycle_warn_default) .and. (.not. recycle_stop_default)) return
nmax = max(na, nb)
nmin = min(na, nb)
if (nmin <= 0) return
if (nmax /= nmin) then
   if (recycle_stop_default) then
      error stop "recycle-stop: vector recycling occurred (lengths differ)"
   end if
   if (recycle_warn_default) then
      write(*,'(a)') "Warning message:"
      write(*,'(3a)') "In ", trim(op), " : vector recycling occurred (lengths differ)"
   end if
end if
end subroutine maybe_warn_recycle

subroutine print_real_scalar(x, int_like)
! Print one real value using integer format when integer-like.
real(kind=dp), intent(in) :: x ! value to print
logical, intent(in), optional :: int_like ! force integer-like formatting
logical :: use_int_like, as_int
integer(kind=int64) :: k
real(kind=dp) :: tol
if (r_is_na_payload(x)) then
   write(*,"(a)") "NA"
   return
end if
use_int_like = print_int_like_default
if (present(int_like)) use_int_like = int_like
as_int = .false.
if (use_int_like) then
   if (ieee_is_finite(x) .and. abs(x) <= real(huge(0_int64), kind=dp)) then
      k = nint(x, kind=int64)
      tol = print_int_like_tol * max(1.0_dp, abs(x))
      as_int = abs(x - real(k, kind=dp)) <= tol
   end if
end if
if (as_int) then
   write(*,"(i0)") k
else
   write(*,"(g0)") x
end if
end subroutine print_real_scalar

subroutine display_integer_scalar(x)
! Display one default-kind integer.
integer, intent(in) :: x
write(*,"(g0)") x
end subroutine display_integer_scalar

subroutine display_logical_scalar(x)
! Display one logical value.
logical, intent(in) :: x
write(*,"(g0)") x
end subroutine display_logical_scalar

subroutine display_complex_scalar(x)
! Display one double-precision complex value.
complex(kind=dp), intent(in) :: x
call print_complex_vector([x])
end subroutine display_complex_scalar

subroutine display_char_scalar(x)
! Display one character value without trailing blanks.
character(len=*), intent(in) :: x
write(*,"(a)") trim(x)
end subroutine display_char_scalar

subroutine print_real_vector(x, int_like, digits)
! Print one real vector; use integer format when all values are integer-like.
real(kind=dp), intent(in) :: x(:) ! values to print
logical, intent(in), optional :: int_like ! force integer-like formatting
integer, intent(in), optional :: digits ! digits after decimal point for display
logical :: use_int_like, all_int
integer :: i
integer(kind=int64) :: k
real(kind=dp) :: r, tol
character(len=32) :: fmt
if (size(x) == 0) then
   write(*,"(a)") "numeric(0)"
   return
end if
if (any(r_is_na_payload(x))) then
   do i = 1, size(x)
      if (r_is_na_payload(x(i))) then
         write(*,"(a)", advance="no") "NA"
      else
         write(*,"(g0)", advance="no") x(i)
      end if
      if (i < size(x)) write(*,"(a)", advance="no") " "
   end do
   write(*,*)
   return
end if
if (present(digits)) then
   write(fmt, '("(20(f0.",i0,",1x,:))")') max(0, digits)
   write(*,fmt) (x(i), i = 1, size(x))
   return
end if
use_int_like = print_int_like_default
if (present(int_like)) use_int_like = int_like
all_int = .false.
if (use_int_like) then
   all_int = .true.
   do i = 1, size(x)
      r = x(i)
      if (.not. ieee_is_finite(r)) then
         all_int = .false.
         exit
      end if
      if (abs(r) > real(huge(0_int64), kind=dp)) then
         all_int = .false.
         exit
      end if
      k = nint(r, kind=int64)
      tol = print_int_like_tol * max(1.0_dp, abs(r))
      if (abs(r - real(k, kind=dp)) > tol) then
         all_int = .false.
         exit
      end if
   end do
end if
if (all_int) then
   write(*,"(20(i0,1x,:))") (nint(x(i), kind=int64), i = 1, size(x))
else
   write(*,"(20(g0,1x,:))") (x(i), i = 1, size(x))
end if
end subroutine print_real_vector

subroutine print_complex_vector(x)
! Print one complex vector using a compact R-like complex format.
complex(kind=dp), intent(in) :: x(:) ! values to print
integer :: i
do i = 1, size(x)
   if (abs(aimag(x(i))) <= 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(real(x(i), kind=dp)))) then
      write(*,"(g0)", advance="no") real(x(i), kind=dp)
   else if (abs(real(x(i), kind=dp)) <= 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(aimag(x(i))))) then
      write(*,"(g0,sp,g0,ss,a)", advance="no") 0.0_dp, aimag(x(i)), "i"
   else
      write(*,"(g0,sp,g0,ss,a)", advance="no") real(x(i), kind=dp), aimag(x(i)), "i"
   end if
   if (i < size(x)) write(*,"(a)", advance="no") " "
end do
write(*,*)
end subroutine print_complex_vector

subroutine print_integer_vector(x)
! Print one integer vector with R-like line wrapping.
integer, intent(in) :: x(:) ! values to print
write (*,"(20(g0,1x,:))") x
write(*,*)
end subroutine print_integer_vector

subroutine display_logical_vector(x)
! Display one logical vector.
logical, intent(in) :: x(:)
write(*,"(20(g0,1x,:))") x
end subroutine display_logical_vector

subroutine print_char_vector(x)
! Print char vector values in an R-like format.
character(len=*), intent(in) :: x(:) ! input vector
integer :: i
do i = 1, size(x)
   write(*,"(a)", advance="no") trim(x(i))
   if (i < size(x)) write(*,"(a)", advance="no") " "
end do
write(*,*)
end subroutine print_char_vector

pure function r_factor_labels(codes, levels) result(out)
! Convert one-based factor codes to their character level labels.
integer, intent(in) :: codes(:)
character(len=*), intent(in) :: levels(:)
character(len=:), allocatable :: out(:)
integer :: i
allocate(character(len=len(levels)) :: out(size(codes)))
do i = 1, size(codes)
   if (codes(i) >= 1 .and. codes(i) <= size(levels)) then
      out(i) = levels(codes(i))
   else
      out(i) = "NA"
   end if
end do
end function r_factor_labels

subroutine print_factor(codes, levels)
! Display factor labels followed by the factor's levels.
integer, intent(in) :: codes(:)
character(len=*), intent(in) :: levels(:)
call print_char_vector(r_factor_labels(codes, levels))
write(*,"(a)", advance="no") "Levels: "
call print_char_vector(levels)
end subroutine print_factor

pure elemental function r_as_real_char(x) result(out)
! Convert a character value to double precision, returning NaN on failure.
character(len=*), intent(in) :: x
real(kind=dp) :: out
out = str_to_real(x)
end function r_as_real_char

pure function r_ifelse_real(test, yes, no) result(out)
! Apply scalar real branches to a logical vector encoded as 0, 1, or NaN.
real(kind=dp), intent(in) :: test(:), yes, no
real(kind=dp), allocatable :: out(:)
allocate(out(size(test)))
out = merge(yes, no, test /= 0.0_dp)
where (.not. ieee_is_finite(test))
   out = ieee_value(0.0_dp, ieee_quiet_nan)
end where
end function r_ifelse_real

pure function r_na_real() result(out)
! Return a quiet NaN carrying R's NA payload.
real(kind=dp) :: out
out = transfer(int(z'7ff00000000007a2', kind=int64), out)
end function r_na_real

pure elemental logical function r_is_na_payload(x) result(out)
! Distinguish the runtime's R NA payload from an ordinary IEEE NaN.
real(kind=dp), intent(in) :: x
out = transfer(x, 0_int64) == int(z'7ff00000000007a2', kind=int64)
end function r_is_na_payload

pure function r_inf() result(out)
! Return positive IEEE infinity.
real(kind=dp) :: out
out = ieee_value(0.0_dp, ieee_positive_inf)
end function r_inf

pure elemental logical function r_is_nan(x) result(out)
! True for ordinary NaN but false for the distinct R NA payload.
real(kind=dp), intent(in) :: x
out = (x /= x) .and. transfer(x, 0_int64) /= int(z'7ff00000000007a2', kind=int64)
end function r_is_nan

subroutine display_char_matrix(x)
! Display a rank-2 character array one row per line.
character(len=*), intent(in) :: x(:,:)
integer :: i
do i = 1, size(x, 1)
   call print_char_vector(x(i, :))
end do
end subroutine display_char_matrix

subroutine display_real_array3(x)
! Display a rank-3 real array as consecutive matrix slices.
real(kind=dp), intent(in) :: x(:,:,:)
integer :: k
do k = 1, size(x, 3)
   call print_matrix_rstyle_real(x(:, :, k))
end do
end subroutine display_real_array3

subroutine display_integer_array3(x)
! Display a rank-3 integer array as consecutive matrix slices.
integer, intent(in) :: x(:,:,:)
integer :: k
do k = 1, size(x, 3)
   call print_matrix_rstyle_int(x(:, :, k))
end do
end subroutine display_integer_array3

subroutine display_logical_array3(x)
! Display a rank-3 logical array as consecutive matrix slices.
logical, intent(in) :: x(:,:,:)
integer :: k
do k = 1, size(x, 3)
   call print_matrix_rstyle_logical(x(:, :, k))
end do
end subroutine display_logical_array3

subroutine display_complex_array3(x)
! Display a rank-3 complex array as consecutive matrix slices.
complex(kind=dp), intent(in) :: x(:,:,:)
integer :: k
do k = 1, size(x, 3)
   call print_matrix_complex(x(:, :, k))
end do
end subroutine display_complex_array3

subroutine display_char_array3(x)
! Display a rank-3 character array as consecutive matrix slices.
character(len=*), intent(in) :: x(:,:,:)
integer :: k
do k = 1, size(x, 3)
   call display_char_matrix(x(:, :, k))
end do
end subroutine display_char_array3

subroutine print_named_real_vector(x, names, digits)
! Print named real vector values in an R-like format.
real(kind=dp), intent(in) :: x(:) ! input vector
character(len=*), intent(in) :: names(:) ! display names
integer, intent(in), optional :: digits ! digits after decimal point for display
call print_char_vector(names)
call print_real_vector(x, digits=digits)
end subroutine print_named_real_vector

subroutine print_named_real_row(x, names, digits, row_name)
! Print a named real vector as a one-row table with aligned labels.
real(kind=dp), intent(in) :: x(:) ! input vector
character(len=*), intent(in) :: names(:) ! column names
integer, intent(in), optional :: digits ! number of digits after decimal
character(len=*), intent(in), optional :: row_name ! optional row name
integer :: j, row_w, col_w, digits_use
character(len=32) :: row_fmt, col_fmt, real_fmt
digits_use = 4
if (present(digits)) digits_use = max(0, min(15, digits))
row_w = 12
if (present(row_name)) row_w = max(row_w, len_trim(row_name))
col_w = max(12, digits_use + 8)
do j = 1, size(names)
   col_w = max(col_w, len_trim(names(j)))
end do
write(row_fmt, '("(a", i0, ",1x)")') row_w
write(col_fmt, '("(a", i0, ",1x)")') col_w
write(real_fmt, '("(f", i0, ".", i0, ",1x)")') col_w, digits_use
write(*,'(a)', advance='no') repeat(" ", row_w + 1)
do j = 1, size(x)
   if (j <= size(names)) then
      write(*,col_fmt, advance='no') trim(names(j))
   else
      write(*,'(i0,1x)', advance='no') j
   end if
end do
write(*,*)
if (present(row_name)) then
   write(*,row_fmt, advance='no') trim(row_name)
else
   write(*,row_fmt, advance='no') ""
end if
do j = 1, size(x)
   if (ieee_is_finite(x(j))) then
      write(*,real_fmt, advance='no') x(j)
   else
      write(*,col_fmt, advance='no') "NA"
   end if
end do
write(*,*)
end subroutine print_named_real_row

pure function nlm_stub(p, hessian) result(out)
! Support nlm-style optimization for stub.
real(kind=dp), intent(in) :: p(:) ! dimension count
logical, intent(in), optional :: hessian
type(nlm_result_t) :: out
integer :: n
if (present(hessian)) continue
n = size(p)
allocate(out%estimate(n))
allocate(out%gradient(n))
allocate(out%hessian(n, n))
out%estimate = p
out%gradient = 0.0_dp
out%hessian = 0.0_dp
out%minimum = 0.0_dp
out%code = 1
out%iterations = 0
end function nlm_stub

function nlm_optimize_scalar(fn, p, hessian, stepmax) result(out)
! Support nlm-style optimization for optimize scalar.
procedure(nlm_objective_scalar) :: fn ! callback procedure
real(kind=dp), intent(in) :: p ! dimension count
logical, intent(in), optional :: hessian ! logical flag
real(kind=dp), intent(in), optional :: stepmax
type(nlm_result_t) :: out
real(kind=dp), allocatable :: pv(:)
allocate(pv(1))
pv(1) = p
if (nlm_method_default == 2) then
   out = nlm_optimize_scalar_newton_impl(fn, pv, hessian, stepmax)
else
   out = nlm_optimize_scalar_legacy_impl(fn, pv, hessian, stepmax)
end if
end function nlm_optimize_scalar

function nlm_optimize_vec(fn, p, hessian, stepmax) result(out)
! Support nlm-style optimization for optimize vec.
procedure(nlm_objective_vec) :: fn ! callback procedure
real(kind=dp), intent(in) :: p(:) ! dimension count
logical, intent(in), optional :: hessian ! logical flag
real(kind=dp), intent(in), optional :: stepmax
type(nlm_result_t) :: out
if (nlm_method_default == 2) then
   out = nlm_optimize_vec_newton(fn, p, hessian, stepmax)
else
   out = nlm_optimize_vec_legacy(fn, p, hessian, stepmax)
end if
end function nlm_optimize_vec

function nlm_optimize_vec_legacy(fn, p, hessian, stepmax) result(out)
! Support nlm-style optimization for optimize vec.
procedure(nlm_objective_vec) :: fn ! callback procedure
real(kind=dp), intent(in) :: p(:) ! dimension count
logical, intent(in), optional :: hessian ! logical flag
real(kind=dp), intent(in), optional :: stepmax
type(nlm_result_t) :: out
integer :: n, iter
real(kind=dp), allocatable :: x(:), g(:), trial(:)
real(kind=dp) :: f, f_trial, gnorm, alpha, max_step, step_norm, tol
n = size(p)
allocate(x(n), g(n), trial(n))
x = p
tol = 1.0e-7_dp
max_step = 100.0_dp
if (present(stepmax)) max_step = max(stepmax, 1.0e-12_dp)
f = fn(x)
do iter = 1, 500
   call nlm_fd_grad_vec(fn, x, g)
   gnorm = sqrt(sum(g*g))
   if (gnorm <= tol * max(1.0_dp, sqrt(sum(x*x)))) exit
   alpha = 1.0_dp
   do
      trial = x - alpha * g
      step_norm = sqrt(sum((trial - x) * (trial - x)))
      if (step_norm > max_step) trial = x + (trial - x) * (max_step / step_norm)
      f_trial = fn(trial)
      if (ieee_is_finite(f_trial) .and. f_trial <= f - 1.0e-4_dp * alpha * gnorm * gnorm) exit
      alpha = alpha * 0.5_dp
      if (alpha < 1.0e-10_dp) exit
   end do
   if (.not. ieee_is_finite(f_trial) .or. f_trial >= f) exit
   if (maxval(abs(trial - x)) <= tol * max(1.0_dp, maxval(abs(x)))) then
      x = trial
      f = f_trial
      exit
   end if
   x = trial
   f = f_trial
end do
allocate(out%estimate(n), out%gradient(n), out%hessian(n, n))
out%estimate = x
out%minimum = f
call nlm_fd_grad_vec(fn, x, out%gradient)
out%hessian = 0.0_dp
if (present(hessian)) then
   if (hessian) call nlm_fd_hessian_vec(fn, x, out%hessian)
end if
out%iterations = min(iter, 500)
out%code = merge(1, 4, sqrt(sum(out%gradient*out%gradient)) <= 1.0e-4_dp)
end function nlm_optimize_vec_legacy

function nlm_optimize_vec_newton(fn, p, hessian, stepmax) result(out)
! Newton-style nlm helper using finite-difference gradient and Hessian.
procedure(nlm_objective_vec) :: fn ! callback procedure
real(kind=dp), intent(in) :: p(:) ! starting parameters
logical, intent(in), optional :: hessian ! logical flag
real(kind=dp), intent(in), optional :: stepmax
type(nlm_result_t) :: out
integer :: n, iter, damp_iter, i
real(kind=dp), allocatable :: x(:), g(:), hess(:,:), hwork(:,:), step(:), trial(:)
real(kind=dp) :: f, f_trial, gnorm, alpha, max_step, step_norm, tol, lambda, slope
logical :: have_step
n = size(p)
allocate(x(n), g(n), hess(n,n), hwork(n,n), step(n), trial(n))
x = p
tol = 1.0e-7_dp
max_step = 100.0_dp
if (present(stepmax)) max_step = max(stepmax, 1.0e-12_dp)
f = fn(x)
do iter = 1, 200
   call nlm_fd_grad_vec(fn, x, g)
   gnorm = sqrt(sum(g*g))
   if (gnorm <= tol * max(1.0_dp, sqrt(sum(x*x)))) exit
   call nlm_fd_hessian_vec(fn, x, hess)
   lambda = 1.0e-8_dp
   have_step = .false.
   do damp_iter = 1, 8
      hwork = hess
      do i = 1, n
         hwork(i,i) = hwork(i,i) + lambda
      end do
      step = solve_real(hwork, -g)
      if (all(ieee_is_finite(step)) .and. dot_product(g, step) < 0.0_dp) then
         have_step = .true.
         exit
      end if
      lambda = lambda * 10.0_dp
   end do
   if (.not. have_step) step = -g
   step_norm = sqrt(sum(step*step))
   if (step_norm > max_step) step = step * (max_step / step_norm)
   slope = dot_product(g, step)
   alpha = 1.0_dp
   do
      trial = x + alpha * step
      f_trial = fn(trial)
      if (ieee_is_finite(f_trial) .and. f_trial <= f + 1.0e-4_dp * alpha * slope) exit
      alpha = alpha * 0.5_dp
      if (alpha < 1.0e-10_dp) exit
   end do
   if (.not. ieee_is_finite(f_trial) .or. f_trial >= f) exit
   if (maxval(abs(trial - x)) <= tol * max(1.0_dp, maxval(abs(x)))) then
      x = trial
      f = f_trial
      exit
   end if
   x = trial
   f = f_trial
end do
allocate(out%estimate(n), out%gradient(n), out%hessian(n, n))
out%estimate = x
out%minimum = f
call nlm_fd_grad_vec(fn, x, out%gradient)
out%hessian = 0.0_dp
if (present(hessian)) then
   if (hessian) call nlm_fd_hessian_vec(fn, x, out%hessian)
end if
out%iterations = min(iter, 200)
out%code = merge(1, 4, sqrt(sum(out%gradient*out%gradient)) <= 1.0e-4_dp)
end function nlm_optimize_vec_newton

function nlm_optimize_scalar_legacy_impl(fn, p, hessian, stepmax) result(out)
! Support nlm-style optimization for optimize scalar impl.
procedure(nlm_objective_scalar) :: fn ! callback procedure
real(kind=dp), intent(in) :: p(:) ! dimension count
logical, intent(in), optional :: hessian ! logical flag
real(kind=dp), intent(in), optional :: stepmax
type(nlm_result_t) :: out
integer :: iter
real(kind=dp) :: x, g, f, f_trial, trial, alpha, max_step, tol, h
x = p(1)
tol = 1.0e-7_dp
max_step = 100.0_dp
if (present(stepmax)) max_step = max(stepmax, 1.0e-12_dp)
f = fn(x)
do iter = 1, 500
   g = nlm_fd_grad_scalar(fn, x)
   if (abs(g) <= tol * max(1.0_dp, abs(x))) exit
   alpha = 1.0_dp
   do
      trial = x - alpha * g
      if (abs(trial - x) > max_step) trial = x - sign(max_step, g)
      f_trial = fn(trial)
      if (ieee_is_finite(f_trial) .and. f_trial <= f - 1.0e-4_dp * alpha * g * g) exit
      alpha = alpha * 0.5_dp
      if (alpha < 1.0e-10_dp) exit
   end do
   if (.not. ieee_is_finite(f_trial) .or. f_trial >= f) exit
   if (abs(trial - x) <= tol * max(1.0_dp, abs(x))) then
      x = trial
      f = f_trial
      exit
   end if
   x = trial
   f = f_trial
end do
allocate(out%estimate(1), out%gradient(1), out%hessian(1, 1))
out%estimate(1) = x
out%minimum = f
out%gradient(1) = nlm_fd_grad_scalar(fn, x)
out%hessian = 0.0_dp
if (present(hessian)) then
   if (hessian) then
      h = sqrt(epsilon(1.0_dp)) * max(1.0_dp, abs(x))
      out%hessian(1, 1) = (fn(x + h) - 2.0_dp * fn(x) + fn(x - h)) / (h * h)
   end if
end if
out%iterations = min(iter, 500)
out%code = merge(1, 4, abs(out%gradient(1)) <= 1.0e-4_dp)
end function nlm_optimize_scalar_legacy_impl

function nlm_optimize_scalar_newton_impl(fn, p, hessian, stepmax) result(out)
! Newton-style scalar nlm helper using finite-difference derivative and Hessian.
procedure(nlm_objective_scalar) :: fn ! callback procedure
real(kind=dp), intent(in) :: p(:) ! dimension count
logical, intent(in), optional :: hessian ! logical flag
real(kind=dp), intent(in), optional :: stepmax
type(nlm_result_t) :: out
integer :: iter
real(kind=dp) :: x, g, f, f_trial, trial, alpha, max_step, tol, h, hess, step, slope
x = p(1)
tol = 1.0e-7_dp
max_step = 100.0_dp
if (present(stepmax)) max_step = max(stepmax, 1.0e-12_dp)
f = fn(x)
do iter = 1, 200
   g = nlm_fd_grad_scalar(fn, x)
   if (abs(g) <= tol * max(1.0_dp, abs(x))) exit
   h = sqrt(epsilon(1.0_dp)) * max(1.0_dp, abs(x))
   hess = (fn(x + h) - 2.0_dp * fn(x) + fn(x - h)) / (h * h)
   if (ieee_is_finite(hess) .and. hess > sqrt(epsilon(1.0_dp))) then
      step = -g / hess
   else
      step = -g
   end if
   if (abs(step) > max_step) step = sign(max_step, step)
   slope = g * step
   if (slope >= 0.0_dp) then
      step = -g
      if (abs(step) > max_step) step = sign(max_step, step)
      slope = g * step
   end if
   alpha = 1.0_dp
   do
      trial = x + alpha * step
      f_trial = fn(trial)
      if (ieee_is_finite(f_trial) .and. f_trial <= f + 1.0e-4_dp * alpha * slope) exit
      alpha = alpha * 0.5_dp
      if (alpha < 1.0e-10_dp) exit
   end do
   if (.not. ieee_is_finite(f_trial) .or. f_trial >= f) exit
   if (abs(trial - x) <= tol * max(1.0_dp, abs(x))) then
      x = trial
      f = f_trial
      exit
   end if
   x = trial
   f = f_trial
end do
allocate(out%estimate(1), out%gradient(1), out%hessian(1, 1))
out%estimate(1) = x
out%minimum = f
out%gradient(1) = nlm_fd_grad_scalar(fn, x)
out%hessian = 0.0_dp
if (present(hessian)) then
   if (hessian) then
      h = sqrt(epsilon(1.0_dp)) * max(1.0_dp, abs(x))
      out%hessian(1, 1) = (fn(x + h) - 2.0_dp * fn(x) + fn(x - h)) / (h * h)
   end if
end if
out%iterations = min(iter, 200)
out%code = merge(1, 4, abs(out%gradient(1)) <= 1.0e-4_dp)
end function nlm_optimize_scalar_newton_impl

function nlm_fd_grad_scalar(fn, x) result(g)
! Support nlm-style optimization for fd grad scalar.
procedure(nlm_objective_scalar) :: fn ! objective function
real(kind=dp), intent(in) :: x ! evaluation point
real(kind=dp) :: g, h
h = sqrt(epsilon(1.0_dp)) * max(1.0_dp, abs(x))
g = (fn(x + h) - fn(x - h)) / (2.0_dp * h)
end function nlm_fd_grad_scalar

subroutine nlm_fd_grad_vec(fn, x, g)
! Support nlm-style optimization for fd grad vec.
procedure(nlm_objective_vec) :: fn ! objective function
real(kind=dp), intent(in) :: x(:) ! evaluation point
real(kind=dp), intent(out) :: g(:) ! gradient values
real(kind=dp), allocatable :: xp(:), xm(:)
integer :: i
real(kind=dp) :: h
allocate(xp(size(x)), xm(size(x)))
do i = 1, size(x)
   xp = x
   xm = x
   h = sqrt(epsilon(1.0_dp)) * max(1.0_dp, abs(x(i)))
   xp(i) = xp(i) + h
   xm(i) = xm(i) - h
   g(i) = (fn(xp) - fn(xm)) / (2.0_dp * h)
end do
end subroutine nlm_fd_grad_vec

subroutine nlm_fd_hessian_vec(fn, x, hess)
! Support nlm-style optimization for fd hessian vec.
procedure(nlm_objective_vec) :: fn ! objective function
real(kind=dp), intent(in) :: x(:) ! evaluation point
real(kind=dp), intent(out) :: hess(:,:) ! Hessian matrix
real(kind=dp), allocatable :: xpp(:), xpm(:), xmp(:), xmm(:)
integer :: i, j, n
real(kind=dp) :: hi, hj
n = size(x)
allocate(xpp(n), xpm(n), xmp(n), xmm(n))
do i = 1, n
   do j = 1, n
      xpp = x
      xpm = x
      xmp = x
      xmm = x
      hi = sqrt(epsilon(1.0_dp)) * max(1.0_dp, abs(x(i)))
      hj = sqrt(epsilon(1.0_dp)) * max(1.0_dp, abs(x(j)))
      xpp(i) = xpp(i) + hi; xpp(j) = xpp(j) + hj
      xpm(i) = xpm(i) + hi; xpm(j) = xpm(j) - hj
      xmp(i) = xmp(i) - hi; xmp(j) = xmp(j) + hj
      xmm(i) = xmm(i) - hi; xmm(j) = xmm(j) - hj
      hess(i,j) = (fn(xpp) - fn(xpm) - fn(xmp) + fn(xmm)) / (4.0_dp * hi * hj)
   end do
end do
end subroutine nlm_fd_hessian_vec

subroutine print_nlm_result(fit)
! Print nlm result values in an R-like format.
type(nlm_result_t), intent(in) :: fit ! input value
write(*,"(a)")
write(*,"(a)", advance="no") "$minimum"
write(*,*)
call print_real_scalar(fit%minimum)
write(*,"(a)", advance="no") "$estimate"
write(*,*)
call print_real_vector(fit%estimate)
write(*,"(a)", advance="no") "$gradient"
write(*,*)
call print_real_vector(fit%gradient)
write(*,"(a)", advance="no") "$code"
write(*,*)
write(*,"(i0)") fit%code
write(*,"(a)", advance="no") "$iterations"
write(*,*)
write(*,"(i0)") fit%iterations
end subroutine print_nlm_result

subroutine print_table1(x, names)
! Print table1 values in an R-like format.
integer, intent(in) :: x(:) ! input vector
character(len=*), intent(in) :: names(:) ! display names
integer :: i
do i = 1, min(size(x), size(names))
   write(*,'(a,1x)', advance='no') trim(names(i))
end do
write(*,*)
do i = 1, size(x)
   write(*,'(i0,1x)', advance='no') x(i)
end do
write(*,*)
end subroutine print_table1

subroutine print_table2_int(x, row_names, col_names, digits)
! Print table2 int values in an R-like format.
integer, intent(in) :: x(:,:) ! input matrix
character(len=*), intent(in) :: row_names(:) ! display names
character(len=*), intent(in) :: col_names(:) ! display names
integer, intent(in), optional :: digits ! accepted for generic compatibility
integer :: i, j, row_w, col_w
character(len=32) :: row_fmt, col_fmt, int_fmt
row_w = 12
do i = 1, size(row_names)
   row_w = max(row_w, len_trim(row_names(i)))
end do
col_w = 12
do j = 1, size(col_names)
   col_w = max(col_w, len_trim(col_names(j)))
end do
write(row_fmt, '("(a", i0, ",1x)")') row_w
write(col_fmt, '("(a", i0, ",1x)")') col_w
write(int_fmt, '("(i", i0, ",1x)")') col_w
write(*,'(a)', advance='no') repeat(" ", row_w + 1)
do j = 1, size(x, 2)
   if (j <= size(col_names)) then
      write(*,col_fmt, advance='no') trim(col_names(j))
   else
      write(*,int_fmt, advance='no') j
   end if
end do
write(*,*)
do i = 1, size(x, 1)
   if (i <= size(row_names)) then
      write(*,row_fmt, advance='no') trim(row_names(i))
   else
      write(*,'(i0,1x)', advance='no') i
   end if
   do j = 1, size(x, 2)
      write(*,int_fmt, advance='no') x(i, j)
   end do
   write(*,*)
end do
end subroutine print_table2_int

subroutine print_table2_real(x, row_names, col_names, digits)
! Print table2 real values in an R-like format.
real(kind=dp), intent(in) :: x(:,:) ! input matrix
character(len=*), intent(in) :: row_names(:) ! display names
character(len=*), intent(in) :: col_names(:) ! display names
integer, intent(in), optional :: digits ! number of digits after decimal
integer :: i, j, row_w, col_w, digits_use
character(len=32) :: row_fmt, col_fmt, real_fmt
digits_use = 4
if (present(digits)) digits_use = max(0, min(15, digits))
row_w = 12
do i = 1, size(row_names)
   row_w = max(row_w, len_trim(row_names(i)))
end do
col_w = max(12, digits_use + 8)
do j = 1, size(col_names)
   col_w = max(col_w, len_trim(col_names(j)))
end do
write(row_fmt, '("(a", i0, ",1x)")') row_w
write(col_fmt, '("(a", i0, ",1x)")') col_w
write(real_fmt, '("(f", i0, ".", i0, ",1x)")') col_w, digits_use
write(*,'(a)', advance='no') repeat(" ", row_w + 1)
do j = 1, size(x, 2)
   if (j <= size(col_names)) then
      write(*,col_fmt, advance='no') trim(col_names(j))
   else
      write(*,'(i0,1x)', advance='no') j
   end if
end do
write(*,*)
do i = 1, size(x, 1)
   if (i <= size(row_names)) then
      write(*,row_fmt, advance='no') trim(row_names(i))
   else
      write(*,'(i0,1x)', advance='no') i
   end if
   do j = 1, size(x, 2)
      if (ieee_is_finite(x(i, j))) then
         write(*,real_fmt, advance='no') x(i, j)
      else
         write(*,col_fmt, advance='no') "NA"
      end if
   end do
   write(*,*)
end do
end subroutine print_table2_real

subroutine print_summary_vec(x)
! Print an R-like summary for a numeric vector.
real(kind=dp), intent(in) :: x(:) ! input vector
real(kind=dp), allocatable :: s(:)
s = summary_vec(x)
write(*,'(a,1x,g0)') "   Min.:", s(1)
write(*,'(a,1x,g0)') "1st Qu.:", s(2)
write(*,'(a,1x,g0)') " Median:", s(3)
write(*,'(a,1x,g0)') "   Mean:", s(4)
write(*,'(a,1x,g0)') "3rd Qu.:", s(5)
write(*,'(a,1x,g0)') "   Max.:", s(6)
end subroutine print_summary_vec

subroutine print_summary_mat(x)
! Print R-like per-column summaries for a numeric matrix.
real(kind=dp), intent(in) :: x(:,:) ! input matrix
real(kind=dp), allocatable :: s(:,:)
character(len=*), parameter :: labels(6) = [character(len=8) :: "Min.   :", "1st Qu.:", "Median :", "Mean   :", &
& "3rd Qu.:", "Max.   :"]
integer :: j
allocate(s(6, size(x, 2)))
do j = 1, size(x, 2)
   s(:, j) = summary_vec(x(:, j))
end do
write(*,'(9x)', advance='no')
do j = 1, size(x, 2)
   write(*,'(a,i0,13x)', advance='no') "V", j
end do
write(*,*)
do j = 1, 6
   write(*,'(a,1x)', advance='no') labels(j)
   write(*,'(*(es14.6,1x))') s(j, :)
end do
end subroutine print_summary_mat

pure function r_seq_int_ii(a, b) result(out)
! Return integer sequence a, a+step, ..., b with step +/-1.
integer, intent(in) :: a, b
integer, allocatable :: out(:)
integer :: i, n, step
n = abs(b - a) + 1
allocate(out(n))
step = merge(1, -1, a <= b)
do i = 1, n
   out(i) = a + (i - 1) * step
end do
end function r_seq_int_ii

pure function r_seq_int_ir(a, b) result(out)
! r_seq_int with a real upper bound (R's `:` accepts numeric operands).
integer, intent(in) :: a
real(kind=dp), intent(in) :: b
integer, allocatable :: out(:)
out = r_seq_int_ii(a, nint(b))
end function r_seq_int_ir

pure function r_seq_int_ri(a, b) result(out)
! r_seq_int with a real lower bound.
real(kind=dp), intent(in) :: a
integer, intent(in) :: b
integer, allocatable :: out(:)
out = r_seq_int_ii(nint(a), b)
end function r_seq_int_ri

pure function r_seq_int_rr(a, b) result(out)
! r_seq_int with real bounds (e.g. to:from where to, from are numeric).
real(kind=dp), intent(in) :: a, b
integer, allocatable :: out(:)
out = r_seq_int_ii(nint(a), nint(b))
end function r_seq_int_rr

pure function r_seq_len(n) result(out)
! Return integer sequence 1..n (empty for n<=0).
integer, intent(in) :: n
integer, allocatable :: out(:)
integer :: i
if (n <= 0) then
   allocate(out(0))
   return
end if
allocate(out(n))
do i = 1, n
   out(i) = i
end do
end function r_seq_len

pure function r_seq_int_by(a, b, by) result(out)
! Return integer sequence from a to b using explicit integer step by.
integer, intent(in) :: a, b, by
integer, allocatable :: out(:)
integer :: i, n
if (by == 0) then
   allocate(out(0))
   return
end if
if ((by > 0 .and. a > b) .or. (by < 0 .and. a < b)) then
   allocate(out(0))
   return
end if
n = (abs(b - a) / abs(by)) + 1
allocate(out(n))
do i = 1, n
   out(i) = a + (i - 1) * by
end do
end function r_seq_int_by

pure function r_seq_int_length(a, b, n) result(out)
! Return n integer values linearly spaced from a to b.
integer, intent(in) :: a, b, n
integer, allocatable :: out(:)
integer :: i
real(kind=dp) :: t
if (n <= 0) then
   allocate(out(0))
   return
end if
allocate(out(n))
if (n == 1) then
   out(1) = a
   return
end if
do i = 1, n
   t = real(i - 1, kind=dp) / real(n - 1, kind=dp)
   out(i) = nint((1.0_dp - t) * real(a, kind=dp) + t * real(b, kind=dp))
end do
end function r_seq_int_length

pure function r_seq_real_by(a, b, by) result(out)
! Return real sequence from a to b with real step by.
real(kind=dp), intent(in) :: a, b, by
real(kind=dp), allocatable :: out(:)
integer :: i, n
if (by == 0.0_dp) then
   allocate(out(0))
   return
end if
if ((by > 0.0_dp .and. a > b) .or. (by < 0.0_dp .and. a < b)) then
   allocate(out(0))
   return
end if
n = int(floor((b - a) / by + 1.0e-12_dp)) + 1
if (n < 0) n = 0
allocate(out(n))
do i = 1, n
   out(i) = a + real(i - 1, kind=dp) * by
end do
end function r_seq_real_by

pure function r_seq_real_length(a, b, n) result(out)
! Return n real values linearly spaced from a to b.
real(kind=dp), intent(in) :: a ! input values
real(kind=dp), intent(in) :: b ! input values
integer, intent(in) :: n
real(kind=dp), allocatable :: out(:)
integer :: i
real(kind=dp) :: t
if (n <= 0) then
   allocate(out(0))
   return
end if
allocate(out(n))
if (n == 1) then
   out(1) = a
   return
end if
do i = 1, n
   t = real(i - 1, kind=dp) / real(n - 1, kind=dp)
   out(i) = (1.0_dp - t) * a + t * b
end do
end function r_seq_real_length

pure function r_rep_real(x, times, each, len_out, times_vec) result(out)
! Repeat elements/blocks of a real vector (R-like rep subset).
real(kind=dp), intent(in) :: x(:) ! input vector
integer, intent(in), optional :: times ! integer argument
integer, intent(in), optional :: each ! integer argument
integer, intent(in), optional :: len_out ! integer argument
integer, intent(in), optional :: times_vec(:)
real(kind=dp), allocatable :: out(:)
real(kind=dp), allocatable :: y(:), z(:)
integer :: i, j, n, e, t, k, m, need, c
n = size(x)
if (n <= 0) then
   allocate(out(0))
   return
end if
if (present(each)) then
   e = each
else
   e = 1
end if
if (e < 1) e = 1
allocate(y(n * e))
k = 0
do i = 1, n
   do j = 1, e
      k = k + 1
      y(k) = x(i)
   end do
end do
if (present(len_out)) then
   need = max(0, len_out)
   allocate(out(need))
   do i = 1, need
      out(i) = y(mod(i - 1, size(y)) + 1)
   end do
   return
end if
if (present(times_vec)) then
   m = size(y)
   c = 0
   do i = 1, m
      t = times_vec(mod(i - 1, size(times_vec)) + 1)
      if (t > 0) c = c + t
   end do
   allocate(z(c))
   k = 0
   do i = 1, m
      t = times_vec(mod(i - 1, size(times_vec)) + 1)
      do j = 1, max(0, t)
         k = k + 1
         z(k) = y(i)
      end do
   end do
else
   if (present(times)) then
      t = times
   else
      t = 1
   end if
   if (t < 0) t = 0
   allocate(z(size(y) * t))
   k = 0
   do j = 1, t
      do i = 1, size(y)
         k = k + 1
         z(k) = y(i)
      end do
   end do
end if
out = z
end function r_rep_real

pure function r_rep_char(x, times, each, len_out, times_vec) result(out)
! Repeat elements/blocks of a character vector (R-like rep subset).
character(len=*), intent(in) :: x(:) ! input vector
integer, intent(in), optional :: times ! integer argument
integer, intent(in), optional :: each ! integer argument
integer, intent(in), optional :: len_out ! integer argument
integer, intent(in), optional :: times_vec(:)
character(len=:), allocatable :: out(:)
character(len=:), allocatable :: y(:), z(:)
integer :: i, j, n, e, t, k, m, need, c
integer :: item_len
n = size(x)
if (n <= 0) then
   allocate(character(len=1) :: out(0))
   return
end if
item_len = max(1, len(x(1)))
if (present(each)) then
   e = each
else
   e = 1
end if
if (e < 1) e = 1
allocate(character(len=item_len) :: y(n * e))
k = 0
do i = 1, n
   do j = 1, e
      k = k + 1
      y(k) = x(i)
   end do
end do
if (present(len_out)) then
   need = max(0, len_out)
   allocate(character(len=item_len) :: out(need))
   do i = 1, need
      out(i) = y(mod(i - 1, size(y)) + 1)
   end do
   return
end if
if (present(times_vec)) then
   m = size(y)
   c = 0
   do i = 1, m
      t = times_vec(mod(i - 1, size(times_vec)) + 1)
      if (t > 0) c = c + t
   end do
   allocate(character(len=item_len) :: z(c))
   k = 0
   do i = 1, m
      t = times_vec(mod(i - 1, size(times_vec)) + 1)
      do j = 1, max(0, t)
         k = k + 1
         z(k) = y(i)
      end do
   end do
else
   if (present(times)) then
      t = times
   else
      t = 1
   end if
   if (t < 0) t = 0
   allocate(character(len=item_len) :: z(size(y) * t))
   k = 0
   do j = 1, t
      do i = 1, size(y)
         k = k + 1
         z(k) = y(i)
      end do
   end do
end if
out = z
end function r_rep_char

pure function r_rep_int(x, times, each, len_out, times_vec) result(out)
! Repeat elements/blocks of an integer vector (R-like rep subset).
integer, intent(in) :: x(:) ! input vector
integer, intent(in), optional :: times ! integer argument
integer, intent(in), optional :: each ! integer argument
integer, intent(in), optional :: len_out ! integer argument
integer, intent(in), optional :: times_vec(:)
integer, allocatable :: out(:)
integer, allocatable :: y(:), z(:)
integer :: i, j, n, e, t, k, m, need, c
n = size(x)
if (n <= 0) then
   allocate(out(0))
   return
end if
if (present(each)) then
   e = each
else
   e = 1
end if
if (e < 1) e = 1
allocate(y(n * e))
k = 0
do i = 1, n
   do j = 1, e
      k = k + 1
      y(k) = x(i)
   end do
end do
if (present(len_out)) then
   need = max(0, len_out)
   allocate(out(need))
   do i = 1, need
      out(i) = y(mod(i - 1, size(y)) + 1)
   end do
   return
end if
if (present(times_vec)) then
   m = size(y)
   c = 0
   do i = 1, m
      t = times_vec(mod(i - 1, size(times_vec)) + 1)
      if (t > 0) c = c + t
   end do
   allocate(z(c))
   k = 0
   do i = 1, m
      t = times_vec(mod(i - 1, size(times_vec)) + 1)
      do j = 1, max(0, t)
         k = k + 1
         z(k) = y(i)
      end do
   end do
else
   if (present(times)) then
      t = times
   else
      t = 1
   end if
   if (t < 0) t = 0
   allocate(z(size(y) * t))
   k = 0
   do j = 1, t
      do i = 1, size(y)
         k = k + 1
         z(k) = y(i)
      end do
   end do
end if
out = z
end function r_rep_int

function runif1() result(u)
! Return one U(0,1) variate.
real(kind=dp) :: u
#ifdef XR2F_USE_R_RNG
u = real(xr2f_r_unif_rand(), kind=dp)
#else
call random_number(u)
#endif
end function runif1

function runif_vec(n) result(x)
! Return a length-n vector of U(0,1) variates.
integer, intent(in) :: n
real(kind=dp), allocatable :: x(:)
#ifdef XR2F_USE_R_RNG
integer :: i
#endif
allocate(x(n))
#ifdef XR2F_USE_R_RNG
do i = 1, n
   x(i) = real(xr2f_r_unif_rand(), kind=dp)
end do
#else
call random_number(x)
#endif
end function runif_vec

function rnorm1() result(x)
! Return one N(0,1) variate (Box-Muller).
real(kind=dp) :: x
real(kind=dp) :: u1, u2
#ifdef XR2F_USE_R_RNG
x = real(xr2f_r_norm_rand(), kind=dp)
#else
do
   u1 = runif1()
   u2 = runif1()
   if (u1 > tiny(1.0_dp)) exit
end do
x = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * acos(-1.0_dp) * u2)
#endif
end function rnorm1

function rnorm_vec(n) result(x)
! Return a length-n vector of N(0,1) variates.
integer, intent(in) :: n
real(kind=dp), allocatable :: x(:)
integer :: i
real(kind=dp) :: u1, u2, r, t
allocate(x(n))
#ifdef XR2F_USE_R_RNG
do i = 1, n
   x(i) = real(xr2f_r_norm_rand(), kind=dp)
end do
#else
i = 1
do while (i <= n)
   u1 = runif1()
   u2 = runif1()
   if (u1 <= tiny(1.0_dp)) cycle
   r = sqrt(-2.0_dp * log(u1))
   t = 2.0_dp * acos(-1.0_dp) * u2
   x(i) = r * cos(t)
   if (i + 1 <= n) x(i + 1) = r * sin(t)
   i = i + 2
end do
#endif
end function rnorm_vec

function rt_vec(n, df) result(x)
! Return a length-n vector of Student t variates with df degrees of freedom.
integer, intent(in) :: n
real(kind=dp), intent(in) :: df
real(kind=dp), allocatable :: x(:)
real(kind=dp), allocatable :: z(:), u(:), chi(:)
if (n <= 0) then
   allocate(x(0))
   return
end if
if (df <= 0.0_dp .or. df /= df) then
   allocate(x(n))
   x = ieee_value(0.0_dp, ieee_quiet_nan)
   return
else if (.not. ieee_is_finite(df)) then
   x = rnorm_vec(n)
   return
end if
z = rnorm_vec(n)
u = runif_vec(n)
chi = qchisq(u, df)
x = z / sqrt(chi / df)
end function rt_vec

function rexp_rng(n, rate) result(x)
! Return n exponential variates.
integer, intent(in) :: n
real(kind=dp), intent(in), optional :: rate
real(kind=dp), allocatable :: x(:)
real(kind=dp) :: rt
rt = 1.0_dp
if (present(rate)) rt = rate
allocate(x(max(0, n)))
if (rt /= rt) then
   x = ieee_value(0.0_dp, ieee_quiet_nan)
else if (.not. ieee_is_finite(rt)) then
   x = 0.0_dp
else if (rt <= 0.0_dp) then
   x = ieee_value(0.0_dp, ieee_quiet_nan)
else
   x = -log(max(tiny(1.0_dp), 1.0_dp - runif_vec(size(x)))) / rt
end if
end function rexp_rng

recursive function rgamma_one(shape, rate) result(x)
! Draw one gamma variate.
! Uses the Marsaglia-Tsang method, adapted from Alan Miller's public-domain
! rng/rgamma.f90 implementation.
real(kind=dp), intent(in) :: shape, rate
real(kind=dp) :: x
real(kind=dp) :: a, rt, d, c, u, v, z
a = shape
rt = rate
if (a == 0.0_dp) then
   x = 0.0_dp
   return
else if (a < 0.0_dp .or. a /= a .or. rt < 0.0_dp .or. rt /= rt) then
   x = ieee_value(1.0_dp, ieee_quiet_nan)
   return
else if (rt == 0.0_dp) then
   x = ieee_value(0.0_dp, ieee_positive_inf)
   return
else if (.not. ieee_is_finite(rt)) then
   x = 0.0_dp
   return
else if (.not. ieee_is_finite(a)) then
   x = ieee_value(0.0_dp, ieee_positive_inf)
   return
end if
if (a < 1.0_dp) then
   u = max(tiny(1.0_dp), runif1())
   x = rgamma_one(a + 1.0_dp, rt) * u**(1.0_dp / a)
   return
end if
d = a - 1.0_dp / 3.0_dp
c = 1.0_dp / sqrt(9.0_dp * d)
do
   do
      z = rnorm1()
      v = (1.0_dp + c * z)**3
      if (v > 0.0_dp) exit
   end do
   u = runif1()
   if (u < 1.0_dp - 0.0331_dp * z**4) exit
   if (log(max(tiny(1.0_dp), u)) < 0.5_dp * z**2 + d * (1.0_dp - v + log(v))) exit
end do
x = d * v / rt
end function rgamma_one

function rgamma_rng(n, shape, rate, scale) result(x)
! Return n gamma(shape, rate) variates.
integer, intent(in) :: n
real(kind=dp), intent(in) :: shape
real(kind=dp), intent(in), optional :: rate, scale
real(kind=dp), allocatable :: x(:)
integer :: i
real(kind=dp) :: rt
rt = 1.0_dp
if (present(rate)) rt = rate
allocate(x(max(0, n)))
if (shape == 0.0_dp) then
   x = 0.0_dp
   return
end if
if (present(scale)) then
   if (scale < 0.0_dp .or. scale /= scale) then
      x = ieee_value(0.0_dp, ieee_quiet_nan)
      return
   else if (scale == 0.0_dp) then
      x = 0.0_dp
      return
   else
      rt = 1.0_dp / scale
   end if
end if
do i = 1, size(x)
   x(i) = rgamma_one(shape, rt)
end do
end function rgamma_rng

function rbeta_rng(n, shape1, shape2) result(x)
! Return n beta(shape1, shape2) variates.
integer, intent(in) :: n
real(kind=dp), intent(in) :: shape1, shape2
real(kind=dp), allocatable :: x(:)
integer :: i
real(kind=dp) :: a, b
allocate(x(max(0, n)))
if (shape1 < 0.0_dp .or. shape2 < 0.0_dp .or. shape1 /= shape1 .or. shape2 /= shape2) then
   x = ieee_value(0.0_dp, ieee_quiet_nan)
   return
else if (.not. ieee_is_finite(shape1) .and. .not. ieee_is_finite(shape2)) then
   x = 0.5_dp
   return
else if (.not. ieee_is_finite(shape1)) then
   x = 1.0_dp
   return
else if (.not. ieee_is_finite(shape2)) then
   x = 0.0_dp
   return
else if (shape1 == 0.0_dp .and. shape2 > 0.0_dp) then
   x = 0.0_dp
   return
else if (shape2 == 0.0_dp .and. shape1 > 0.0_dp) then
   x = 1.0_dp
   return
else if (shape1 == 0.0_dp .and. shape2 == 0.0_dp) then
   x = merge(0.0_dp, 1.0_dp, runif_vec(size(x)) < 0.5_dp)
   return
end if
do i = 1, size(x)
   a = rgamma_one(shape1, 1.0_dp)
   b = rgamma_one(shape2, 1.0_dp)
   x(i) = a / max(tiny(1.0_dp), a + b)
end do
end function rbeta_rng

function rchisq_rng(n, df) result(x)
! Return n chi-square variates.
integer, intent(in) :: n
real(kind=dp), intent(in) :: df
real(kind=dp), allocatable :: x(:)
allocate(x(max(0, n)))
if (df < 0.0_dp .or. df /= df .or. .not. ieee_is_finite(df)) then
   x = ieee_value(0.0_dp, ieee_quiet_nan)
else if (df == 0.0_dp) then
   x = 0.0_dp
else
   x = rgamma(n, shape=0.5_dp * df, rate=0.5_dp)
end if
end function rchisq_rng

function rf_rng_vec(n, df1, df2) result(x)
! Return n F variates.
integer, intent(in) :: n
real(kind=dp), intent(in) :: df1, df2
real(kind=dp), allocatable :: x(:)
real(kind=dp), allocatable :: a(:), b(:)
if (df1 <= 0.0_dp .or. df2 <= 0.0_dp .or. df1 /= df1 .or. df2 /= df2) then
   allocate(x(max(0, n)))
   x = ieee_value(0.0_dp, ieee_quiet_nan)
   return
else if (.not. ieee_is_finite(df1) .and. .not. ieee_is_finite(df2)) then
   allocate(x(max(0, n)))
   x = 1.0_dp
   return
else if (.not. ieee_is_finite(df1)) then
   b = rchisq(n, df2) / df2
   x = 1.0_dp / max(tiny(1.0_dp), b)
   return
else if (.not. ieee_is_finite(df2)) then
   x = rchisq(n, df1) / df1
   return
end if
a = rchisq(n, df1) / max(tiny(1.0_dp), df1)
b = rchisq(n, df2) / max(tiny(1.0_dp), df2)
x = a / max(tiny(1.0_dp), b)
end function rf_rng_vec

function rlogis_rng(n, location, scale) result(x)
! Return n logistic variates.
integer, intent(in) :: n
real(kind=dp), intent(in), optional :: location, scale
real(kind=dp), allocatable :: x(:)
real(kind=dp), allocatable :: u(:)
real(kind=dp) :: loc, sc
loc = 0.0_dp
sc = 1.0_dp
if (present(location)) loc = location
if (present(scale)) sc = scale
if (sc /= sc .or. loc /= loc) then
   allocate(x(max(0, n)))
   x = ieee_value(0.0_dp, ieee_quiet_nan)
else if (.not. ieee_is_finite(loc)) then
   allocate(x(max(0, n)))
   if (ieee_is_finite(sc)) then
      x = loc
   else
      x = ieee_value(0.0_dp, ieee_quiet_nan)
   end if
else if (sc == 0.0_dp) then
   allocate(x(max(0, n)))
   x = loc
else
   u = runif_vec(max(0, n))
   x = loc + sc * log(max(tiny(1.0_dp), u) / max(tiny(1.0_dp), 1.0_dp - u))
end if
end function rlogis_rng

function rlnorm_rng(n, meanlog, sdlog) result(x)
! Return n lognormal variates.
integer, intent(in) :: n
real(kind=dp), intent(in), optional :: meanlog, sdlog
real(kind=dp), allocatable :: x(:)
real(kind=dp) :: mu, sig
mu = 0.0_dp
sig = 1.0_dp
if (present(meanlog)) mu = meanlog
if (present(sdlog)) sig = sdlog
if (sig < 0.0_dp .or. sig /= sig) then
   allocate(x(max(0, n)))
   x = ieee_value(0.0_dp, ieee_quiet_nan)
else if (mu /= mu) then
   allocate(x(max(0, n)))
   x = ieee_value(0.0_dp, ieee_quiet_nan)
else if (.not. ieee_is_finite(mu)) then
   allocate(x(max(0, n)))
   if (ieee_is_finite(sig)) then
      x = exp(mu)
   else
      x = ieee_value(0.0_dp, ieee_quiet_nan)
   end if
else if (sig == 0.0_dp) then
   allocate(x(max(0, n)))
   x = exp(mu)
else
   x = exp(mu + sig * rnorm_vec(max(0, n)))
end if
end function rlnorm_rng

function rweibull_rng(n, shape, scale) result(x)
! Return n Weibull variates.
integer, intent(in) :: n
real(kind=dp), intent(in) :: shape
real(kind=dp), intent(in), optional :: scale
real(kind=dp), allocatable :: x(:)
real(kind=dp) :: sc
sc = 1.0_dp
if (present(scale)) sc = scale
if (sc == 0.0_dp) then
   allocate(x(max(0, n)))
   x = 0.0_dp
else if (shape <= 0.0_dp .or. shape /= shape .or. .not. ieee_is_finite(shape) .or. &
         sc < 0.0_dp .or. sc /= sc .or. .not. ieee_is_finite(sc)) then
   allocate(x(max(0, n)))
   x = ieee_value(0.0_dp, ieee_quiet_nan)
else
   x = qweibull(runif_vec(max(0, n)), shape=shape, scale=sc)
end if
end function rweibull_rng

function rcauchy_rng(n, location, scale) result(x)
! Return n Cauchy variates.
integer, intent(in) :: n
real(kind=dp), intent(in), optional :: location, scale
real(kind=dp), allocatable :: x(:)
real(kind=dp) :: loc, sc
loc = 0.0_dp
sc = 1.0_dp
if (present(location)) loc = location
if (present(scale)) sc = scale
if (sc < 0.0_dp .or. sc /= sc) then
   allocate(x(max(0, n)))
   x = ieee_value(0.0_dp, ieee_quiet_nan)
else if (loc /= loc) then
   allocate(x(max(0, n)))
   x = ieee_value(0.0_dp, ieee_quiet_nan)
else if (.not. ieee_is_finite(loc)) then
   allocate(x(max(0, n)))
   if (ieee_is_finite(sc)) then
      x = loc
   else
      x = ieee_value(0.0_dp, ieee_quiet_nan)
   end if
else if (sc == 0.0_dp) then
   allocate(x(max(0, n)))
   x = loc
else
   x = qcauchy(runif_vec(max(0, n)), location=loc, scale=sc)
end if
end function rcauchy_rng

function rgeom_rng(n, prob) result(x)
! Return n geometric variates.
integer, intent(in) :: n
real(kind=dp), intent(in) :: prob
integer, allocatable :: x(:)
if (prob <= 0.0_dp .or. prob > 1.0_dp .or. prob /= prob) then
   allocate(x(max(0, n)))
   x = -huge(0)
else
   x = int(qgeom(runif_vec(max(0, n)), prob=prob))
end if
end function rgeom_rng

function rnbinom_rng(n, size_, prob, mu) result(x)
! Return n negative-binomial variates using the gamma-Poisson mixture.
integer, intent(in) :: n
real(kind=dp), intent(in) :: size_
real(kind=dp), intent(in), optional :: prob, mu
integer, allocatable :: x(:)
integer :: i
real(kind=dp) :: p, m, lambda
allocate(x(max(0, n)))
if (size_ <= 0.0_dp .or. size_ /= size_) then
   x = -huge(0)
   return
end if
if (present(mu)) then
   if (mu < 0.0_dp .or. mu /= mu) then
      x = -huge(0)
      return
   else if (mu == 0.0_dp) then
      x = 0
      return
   end if
else if (present(prob)) then
   if (prob <= 0.0_dp .or. prob > 1.0_dp .or. prob /= prob) then
      x = -huge(0)
      return
   else if (prob == 1.0_dp) then
      x = 0
      return
   end if
end if
p = 0.5_dp
if (present(prob)) p = max(tiny(1.0_dp), min(1.0_dp - tiny(1.0_dp), prob))
do i = 1, size(x)
   if (present(mu)) then
      m = max(0.0_dp, mu)
      lambda = rgamma_one(size_, size_ / max(tiny(1.0_dp), m))
   else
      lambda = rgamma_one(size_, p / max(tiny(1.0_dp), 1.0_dp - p))
   end if
   x(i) = rpois_one(lambda)
end do
end function rnbinom_rng

function rhyper_rng(n, m, nwhite, k) result(x)
! Return n hypergeometric variates.
integer, intent(in) :: n, m, nwhite, k
integer, allocatable :: x(:)
integer :: lower, upper
if (m < 0 .or. nwhite < 0 .or. k < 0 .or. k > m + nwhite) then
   allocate(x(max(0, n)))
   x = -huge(0)
   return
end if
lower = max(0, k - nwhite)
upper = min(k, m)
if (lower == upper) then
   allocate(x(max(0, n)))
   x = lower
   return
end if
x = int(qhyper(runif_vec(max(0, n)), m=m, n=nwhite, k=k))
end function rhyper_rng

function rwilcox_rng(n, m, n2) result(x)
! Return n Wilcoxon rank-sum variates.
integer, intent(in) :: n, m, n2
integer, allocatable :: x(:)
if (m < 0 .or. n2 < 0) then
   allocate(x(max(0, n)))
   x = -huge(0)
else if (m == 0 .or. n2 == 0) then
   allocate(x(max(0, n)))
   x = 0
else
   x = int(qwilcox(runif_vec(max(0, n)), m=m, n=n2))
end if
end function rwilcox_rng

function rsignrank_rng(n, n_obs) result(x)
! Return n Wilcoxon signed-rank variates.
integer, intent(in) :: n, n_obs
integer, allocatable :: x(:)
if (n_obs < 0) then
   allocate(x(max(0, n)))
   x = -huge(0)
else if (n_obs == 0) then
   allocate(x(max(0, n)))
   x = 0
else
   x = int(qsignrank(runif_vec(max(0, n)), n=n_obs))
end if
end function rsignrank_rng

function rmultinom_rng(n, size_, prob) result(x)
! Return a category-by-draw matrix of multinomial variates.
integer, intent(in) :: n, size_
real(kind=dp), intent(in) :: prob(:)
integer, allocatable :: x(:,:)
real(kind=dp), allocatable :: p(:), cp(:)
integer :: i, j, cat, k
real(kind=dp) :: total, u
k = size(prob)
allocate(x(k, max(0, n)))
if (k <= 0) return
p = max(0.0_dp, prob)
total = sum(p)
if (total <= 0.0_dp) p = 1.0_dp
p = p / sum(p)
allocate(cp(k))
cp = cumsum(p)
cp(k) = 1.0_dp
x = 0
if (count(p > 0.0_dp) == 1) then
   cat = maxloc(p, dim=1)
   x(cat, :) = max(0, size_)
   return
end if
do j = 1, size(x, 2)
   do i = 1, max(0, size_)
      u = runif1()
      do cat = 1, k
         if (u <= cp(cat)) exit
      end do
      x(min(cat, k), j) = x(min(cat, k), j) + 1
   end do
end do
end function rmultinom_rng

function rnorm_mat(nrow, ncol) result(x)
! Return an nrow-by-ncol matrix of N(0,1) variates.
integer, intent(in) :: nrow, ncol
real(kind=dp), allocatable :: x(:,:)
real(kind=dp), allocatable :: v(:)
v = rnorm_vec(nrow * ncol)
x = reshape(v, [nrow, ncol])
end function rnorm_mat

pure function r_filter_linear(x, filt, sides) result(out)
! Minimal stats::filter-style linear convolution for numeric vectors.
real(kind=dp), intent(in) :: x(:) ! input vector
real(kind=dp), intent(in) :: filt(:) ! input vector
integer, intent(in), optional :: sides
real(kind=dp), allocatable :: out(:)
integer :: n, nf, i, j, idx, side, before, after
real(kind=dp) :: acc, nanv
n = size(x)
nf = size(filt)
side = 2
if (present(sides)) side = sides
nanv = ieee_value(1.0_dp, ieee_quiet_nan)
allocate(out(n))
out = nanv
if (n <= 0 .or. nf <= 0) return
if (side == 1) then
   before = nf - 1
   after = 0
else
   before = nf / 2
   after = nf - before - 1
end if
do i = 1, n
   if (i - before < 1 .or. i + after > n) cycle
   acc = 0.0_dp
   do j = 1, nf
      idx = i + after - j + 1
      acc = acc + filt(j) * x(idx)
   end do
   out(i) = acc
end do
end function r_filter_linear

pure function runmed(x, k, endrule, algorithm, na_action, print_level) result(out)
! Compute R-like runmed smoothing output.
real(kind=dp), intent(in) :: x(:) ! input vector
integer, intent(in) :: k
character(len=*), intent(in), optional :: endrule ! endpoint rule: median, keep, constant
character(len=*), intent(in), optional :: algorithm ! accepted for API compatibility
character(len=*), intent(in), optional :: na_action ! accepted for API compatibility
integer, intent(in), optional :: print_level ! accepted for API compatibility
real(kind=dp), allocatable :: out(:), y_end(:)
integer :: i, j, n, kk, h, n_1, n_2
character(len=32) :: erule
n = size(x)
allocate(out(n))
if (n <= 0) return
kk = max(1, k)
if (mod(kk, 2) == 0) kk = kk + 1
if (kk > n) then
   kk = n
   if (mod(kk, 2) == 0) kk = kk - 1
   kk = max(1, kk)
end if
h = kk / 2
if (kk <= 1) then
   out = x
   return
end if
erule = "median"
if (present(endrule)) erule = trim(adjustl(endrule))
out = x
do i = h + 1, n - h
   out(i) = median(x(i - h:i + h))
end do
if (erule == "constant" .or. erule == "Constant" .or. erule == "CONSTANT") then
   out(1:h) = out(h + 1)
   out(n - h + 1:n) = out(n - h)
end if
if (erule == "median" .or. erule == "Median" .or. erule == "MEDIAN") then
   y_end = out
   n_1 = n - 1
   n_2 = n - 2
   if (h >= 2) then
      out(2) = median(y_end(1:3))
      out(n_1) = median([y_end(n), y_end(n_1), y_end(n_2)])
      if (h >= 3) then
         do i = 3, h
            j = 2 * i - 1
            out(i) = median(y_end(1:j))
            out(n - i + 1) = median(y_end(n + 1 - j:n))
         end do
      end if
   end if
   out(1) = median([y_end(1), out(2), out(2) - 2.0_dp * (out(3) - out(2))])
   out(n) = median([y_end(n), out(n_1), out(n_1) - 2.0_dp * (out(n_2) - out(n_1))])
end if
if (present(algorithm)) continue
if (present(na_action)) continue
if (present(print_level)) continue
end function runmed

pure recursive function smooth(x, kind, twiceit, endrule, do_ends) result(out)
! Compute R-like smooth smoothing output.
real(kind=dp), intent(in) :: x(:) ! input vector
character(len=*), intent(in), optional :: kind
logical, intent(in), optional :: twiceit
character(len=*), intent(in), optional :: endrule
logical, intent(in), optional :: do_ends
real(kind=dp), allocatable :: out(:), prev(:), split(:)
integer :: i, j, n, iter
logical :: use_3rs3r, repeat_medians, do_twice
character(len=:), allocatable :: k, runmed_endrule
n = size(x)
allocate(out(n))
if (n <= 0) return
k = "3RS3R"
if (present(kind)) k = kind
runmed_endrule = "median"
if (present(endrule)) then
   if (endrule == "copy") runmed_endrule = "keep"
end if
do_twice = .false.
if (present(twiceit)) do_twice = twiceit
use_3rs3r = index(k, "3RS3R") > 0
repeat_medians = index(k, "3R") > 0
out = runmed(x, 3, endrule=runmed_endrule)
if (.not. (use_3rs3r .or. repeat_medians)) then
   if (do_twice) out = out + smooth(x - out, kind=k, twiceit=.false., endrule=endrule, do_ends=do_ends)
   return
end if
if (n < 3) return
do iter = 1, max(1, n)
   prev = out
   out = runmed(out, 3, endrule=runmed_endrule)
   if (all(abs(out - prev) <= 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, max(maxval(abs(out)), maxval(abs(prev)))))) exit
end do
if (.not. use_3rs3r) then
   if (do_twice) out = out + smooth(x - out, kind=k, twiceit=.false., endrule=endrule, do_ends=do_ends)
   return
end if
split = out
i = 2
do while (i <= n - 2)
   j = i
   do
      if (j >= n) exit
      if (out(j + 1) /= out(i)) exit
      j = j + 1
   end do
   if (j == i + 1) then
      if (out(i - 1) /= out(i) .and. out(j + 1) /= out(i)) then
         if ((out(i) > out(i - 1) .and. out(i) > out(j + 1)) .or. (out(i) < out(i - 1) .and. out(i) < out(j + 1))) then
            split(i) = out(i - 1)
            split(j) = out(j + 1)
         end if
      end if
   end if
   i = j + 1
end do
out = runmed(split, 3, endrule=runmed_endrule)
if (do_twice) out = out + smooth(x - out, kind=k, twiceit=.false., endrule=endrule, do_ends=do_ends)
if (present(do_ends)) continue
end function smooth

pure function smooth_kernel_eval(x, y, x0, bandwidth, kernel) result(v)
! Compute R-like smooth kernel eval smoothing output.
real(kind=dp), intent(in) :: x(:) ! input vector
real(kind=dp), intent(in) :: y(:) ! response values
real(kind=dp), intent(in) :: x0 ! input value
real(kind=dp), intent(in) :: bandwidth ! input value
character(len=*), intent(in) :: kernel
real(kind=dp) :: v
real(kind=dp) :: bw, z, w, sw, sy, scale
integer :: i
bw = max(bandwidth, sqrt(tiny(1.0_dp)))
sw = 0.0_dp
sy = 0.0_dp
if (kernel == "box") then
   scale = 0.5_dp * bw
else
   ! R ksmooth bandwidth is the interquartile width of the normal kernel.
   scale = bw * (0.25_dp / 0.6744897501960817_dp)
end if
scale = max(scale, sqrt(tiny(1.0_dp)))
do i = 1, min(size(x), size(y))
   z = (x0 - x(i)) / scale
   if (kernel == "box") then
      if (abs(z) <= 1.0_dp) then
         w = 1.0_dp
      else
         w = 0.0_dp
      end if
   else
      if (abs(z) <= 4.0_dp) then
         w = exp(-0.5_dp * z * z)
      else
         w = 0.0_dp
      end if
   end if
   sw = sw + w
   sy = sy + w * y(i)
end do
if (sw > 0.0_dp) then
   v = sy / sw
else
   v = ieee_value(1.0_dp, ieee_quiet_nan)
end if
end function smooth_kernel_eval

pure function ksmooth(x, y, kernel, bandwidth, x_points, range_x, n_points) result(out)
! Compute R-like ksmooth smoothing output.
real(kind=dp), intent(in) :: x(:) ! input vector
real(kind=dp), intent(in) :: y(:) ! response values
character(len=*), intent(in), optional :: kernel ! input string
real(kind=dp), intent(in), optional :: bandwidth ! input value
real(kind=dp), intent(in), optional :: x_points(:) ! prediction locations
real(kind=dp), intent(in), optional :: range_x(:) ! default prediction range
integer, intent(in), optional :: n_points ! number of default prediction locations
type(smooth_xy_t) :: out
character(len=:), allocatable :: kern
real(kind=dp) :: bw, xmin, xmax
integer :: i, ngrid
kern = "normal"
if (present(kernel)) kern = kernel
bw = 0.5_dp
if (present(bandwidth)) bw = bandwidth
if (present(x_points)) then
   out%x = sort(x_points)
else
   ngrid = max(100, size(x))
   if (present(n_points)) ngrid = max(1, n_points)
   allocate(out%x(ngrid))
   if (size(x) <= 0) then
      out%x = ieee_value(0.0_dp, ieee_quiet_nan)
   else if (ngrid <= 1) then
      if (present(range_x) .and. size(range_x) >= 1) then
         out%x = range_x(1)
      else
         out%x = minval(x)
      end if
   else
      if (present(range_x) .and. size(range_x) >= 2) then
         xmin = range_x(1)
         xmax = range_x(2)
      else
         xmin = minval(x)
         xmax = maxval(x)
      end if
      do i = 1, ngrid
         out%x(i) = xmin + (xmax - xmin) * real(i - 1, kind=dp) / real(ngrid - 1, kind=dp)
      end do
   end if
end if
allocate(out%y(size(out%x)))
do i = 1, size(out%x)
   out%y(i) = smooth_kernel_eval(x, y, out%x(i), bw, kern)
end do
out%df = real(size(out%x), kind=dp)
end function ksmooth

pure function lowess(x, y, f, iter, delta) result(out)
! Compute R-like lowess smoothing output.
real(kind=dp), intent(in) :: x(:) ! input vector
real(kind=dp), intent(in) :: y(:) ! response values
real(kind=dp), intent(in), optional :: f ! callback procedure
integer, intent(in), optional :: iter
real(kind=dp), intent(in), optional :: delta
type(smooth_xy_t) :: out
real(kind=dp), allocatable :: robust(:), fitted(:), resid(:), abs_resid(:), xs(:), ys(:), xfit(:), yfit(:)
integer, allocatable :: ord(:)
real(kind=dp) :: frac, h, h1, h9, r, u, w, sw, sx, sy, sxx, sxy, den, beta, alpha, cmad, del, cut
integer :: n, i, j, ns, pass, niter, nleft, nright, last, next_i
n = min(size(x), size(y))
allocate(out%x(n), out%y(n))
if (n <= 0) return
allocate(xs(n), ys(n), ord(n))
ord = order_real(x(1:n))
do i = 1, n
   xs(i) = x(ord(i))
   ys(i) = y(ord(i))
end do
out%x = xs
frac = 2.0_dp / 3.0_dp
if (present(f)) frac = max(0.01_dp, min(1.0_dp, f))
niter = 3
if (present(iter)) niter = max(0, iter)
del = 0.01_dp * (maxval(xs) - minval(xs))
if (present(delta)) del = max(0.0_dp, delta)
ns = max(2, min(n, int(frac * real(n, kind=dp))))
allocate(robust(n), fitted(n), resid(n), abs_resid(n), xfit(n), yfit(n))
robust = 1.0_dp
fitted = ys
do pass = 0, niter
   nleft = 1
   nright = ns
   last = 0
   i = 1
   do
      if (i > n) exit
      do while (nright < n)
         if (xs(i) - xs(nleft) <= xs(nright + 1) - xs(i)) exit
         nleft = nleft + 1
         nright = nright + 1
      end do
      h = max(xs(i) - xs(nleft), xs(nright) - xs(i))
      h1 = 0.001_dp * h
      h9 = 0.999_dp * h
      sw = 0.0_dp
      sx = 0.0_dp
      sy = 0.0_dp
      sxx = 0.0_dp
      sxy = 0.0_dp
      do j = nleft, nright
         r = abs(xs(j) - xs(i))
         if (h <= sqrt(tiny(1.0_dp))) then
            w = merge(robust(j), 0.0_dp, r <= sqrt(tiny(1.0_dp)))
         else if (r <= h9) then
            if (r <= h1) then
               w = robust(j)
            else
               r = r / h
               w = (1.0_dp - r**3)**3 * robust(j)
            end if
         else
            w = 0.0_dp
         end if
         sw = sw + w
         sx = sx + w * xs(j)
         sy = sy + w * ys(j)
         sxx = sxx + w * xs(j) * xs(j)
         sxy = sxy + w * xs(j) * ys(j)
      end do
      if (sw <= sqrt(tiny(1.0_dp))) then
         fitted(i) = ys(i)
      else
         den = sw * sxx - sx * sx
         if (abs(den) <= 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(sw * sxx), abs(sx * sx))) then
            fitted(i) = sy / sw
         else
            beta = (sw * sxy - sx * sy) / den
            alpha = (sy - beta * sx) / sw
            fitted(i) = alpha + beta * xs(i)
         end if
      end if
      if (last > 0 .and. i > last + 1) then
         do j = last + 1, i - 1
            if (xs(i) == xs(last)) then
               fitted(j) = fitted(i)
            else
               fitted(j) = fitted(last) + (fitted(i) - fitted(last)) * (xs(j) - xs(last)) / (xs(i) - xs(last))
            end if
         end do
      end if
      last = i
      cut = xs(last) + del
      next_i = last + 1
      do while (next_i <= n)
         if (xs(next_i) > cut) exit
         if (xs(next_i) > xs(last)) fitted(next_i) = fitted(last)
         next_i = next_i + 1
      end do
      if (next_i > n) then
         i = n
      else
         i = max(last + 1, next_i - 1)
      end if
      if (i <= last) exit
   end do
   if (last < n) fitted(last + 1:n) = fitted(last)
   if (pass >= niter) exit
   resid = ys - fitted
   abs_resid = abs(resid)
   cmad = median(abs_resid)
   if (cmad <= 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(ys)))) exit
   do j = 1, n
      u = abs_resid(j) / (6.0_dp * cmad)
      if (u >= 1.0_dp) then
         robust(j) = 0.0_dp
      else
         robust(j) = (1.0_dp - u * u)**2
      end if
   end do
end do
out%y = fitted
out%df = real(ns, kind=dp)
end function lowess

pure function loess_fit(x, y, span, degree) result(out)
! Runtime helper for R-compatible loess fit.
real(kind=dp), intent(in) :: x(:) ! input vector
real(kind=dp), intent(in) :: y(:) ! response values
real(kind=dp), intent(in), optional :: span ! input value
integer, intent(in), optional :: degree
type(loess_fit_t) :: out
out%x = x
out%y = y
out%span = 0.75_dp
if (present(span)) out%span = span
out%degree = 2
if (present(degree)) out%degree = degree
end function loess_fit

pure function predict_loess(fit, xnew) result(yhat)
! Evaluate distribution helper predict_loess.
type(loess_fit_t), intent(in) :: fit ! input value
real(kind=dp), intent(in) :: xnew(:)
real(kind=dp), allocatable :: yhat(:)
real(kind=dp), allocatable :: dist(:), dist_sorted(:)
real(kind=dp) :: h, d, w, z, sw, sy, sx, sxx, sxy, sx3, sx4, sx2y
real(kind=dp) :: den, alpha, beta, gamma, det, b0, b1, b2
real(kind=dp) :: m00, m01, m02, m11, m12, m22
integer :: i, j, n, ns, deg
allocate(yhat(size(xnew)))
if (.not. allocated(fit%x) .or. .not. allocated(fit%y)) then
   yhat = ieee_value(1.0_dp, ieee_quiet_nan)
   return
end if
n = min(size(fit%x), size(fit%y))
if (n <= 0) then
   yhat = ieee_value(1.0_dp, ieee_quiet_nan)
   return
end if
ns = max(2, min(n, int(fit%span * real(n, kind=dp))))
deg = max(0, min(2, fit%degree))
allocate(dist(n), dist_sorted(n))
do i = 1, size(xnew)
   do j = 1, n
      dist(j) = abs(fit%x(j) - xnew(i))
   end do
   dist_sorted = sort(dist)
   h = dist_sorted(ns)
   if (h <= sqrt(tiny(1.0_dp))) h = maxval(dist_sorted)
   sw = 0.0_dp
   sy = 0.0_dp
   sx = 0.0_dp
   sxx = 0.0_dp
   sxy = 0.0_dp
   sx3 = 0.0_dp
   sx4 = 0.0_dp
   sx2y = 0.0_dp
   do j = 1, n
      z = fit%x(j) - xnew(i)
      if (h <= sqrt(tiny(1.0_dp))) then
         w = merge(1.0_dp, 0.0_dp, abs(z) <= sqrt(tiny(1.0_dp)))
      else if (abs(z) <= h) then
         d = abs(z) / h
         w = (1.0_dp - d**3)**3
      else
         w = 0.0_dp
      end if
      sw = sw + w
      sy = sy + w * fit%y(j)
      sx = sx + w * z
      sxx = sxx + w * z * z
      sxy = sxy + w * z * fit%y(j)
      if (deg >= 2) then
         sx3 = sx3 + w * z**3
         sx4 = sx4 + w * z**4
         sx2y = sx2y + w * z * z * fit%y(j)
      end if
   end do
   if (sw <= sqrt(tiny(1.0_dp))) then
      yhat(i) = fit%y(min(max(1, minloc(dist, dim=1)), n))
   else if (deg <= 0) then
      yhat(i) = sy / sw
   else if (deg == 1) then
      den = sw * sxx - sx * sx
      if (abs(den) <= 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(sw * sxx), abs(sx * sx))) then
         yhat(i) = sy / sw
      else
         beta = (sw * sxy - sx * sy) / den
         alpha = (sy - beta * sx) / sw
         yhat(i) = alpha
      end if
   else
      m00 = sw
      m01 = sx
      m02 = sxx
      m11 = sxx
      m12 = sx3
      m22 = sx4
      det = m00 * (m11 * m22 - m12 * m12) - m01 * (m01 * m22 - m12 * m02) + m02 * (m01 * m12 - m11 * m02)
      if (abs(det) <= 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(m00*m11*m22))) then
         den = sw * sxx - sx * sx
         if (abs(den) <= 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(sw * sxx), abs(sx * sx))) then
            yhat(i) = sy / sw
         else
            beta = (sw * sxy - sx * sy) / den
            alpha = (sy - beta * sx) / sw
            yhat(i) = alpha
         end if
      else
         b0 = sy
         b1 = sxy
         b2 = sx2y
         gamma = (m00 * (m11 * b2 - b1 * m12) - m01 * (m01 * b2 - b1 * m02) + b0 * (m01 * m12 - m11 * m02)) / det
         beta = (m00 * (b1 * m22 - m12 * b2) - b0 * (m01 * m22 - m12 * m02) + m02 * (m01 * b2 - b1 * m02)) / det
         alpha = (b0 * (m11 * m22 - m12 * m12) - m01 * (b1 * m22 - m12 * b2) + m02 * (b1 * m12 - m11 * b2)) / det
         yhat(i) = alpha
      end if
   end if
end do
end function predict_loess

pure function smooth_spline(x, y, df, spar) result(out)
! Compute R-like smooth spline smoothing output.
real(kind=dp), intent(in) :: x(:) ! input vector
real(kind=dp), intent(in) :: y(:) ! response values
real(kind=dp), intent(in), optional :: df, spar
type(smooth_spline_fit_t) :: out
type(loess_fit_t) :: lf
real(kind=dp) :: span_eff, df_eff, spar_eff
integer :: n
logical :: use_df, use_spar
n = min(size(x), size(y))
out%x = x(1:n)
if (n <= 0) then
   out%df = 0.0_dp
   return
end if
use_df = .false.
if (present(df)) use_df = df > 0.0_dp
use_spar = .false.
if (present(spar)) use_spar = spar >= 0.0_dp
if (.not. use_df .and. .not. use_spar) then
   out%df = real(n, kind=dp)
   out%y = y(1:n)
   return
end if
if (use_df) then
   df_eff = max(2.0_dp, min(real(n, kind=dp), df))
   out%df = df_eff
   span_eff = min(1.0_dp, max(0.08_dp, 1.78_dp / sqrt(df_eff)))
else
   spar_eff = max(0.0_dp, min(1.5_dp, spar))
   out%df = max(2.0_dp, real(n, kind=dp) * exp(-2.836_dp * spar_eff**1.32_dp))
   span_eff = min(1.0_dp, max(0.08_dp, 0.18_dp + 0.82_dp * spar_eff))
end if
lf = loess_fit(out%x, y(1:n), span=span_eff, degree=2)
out%y = predict_loess(lf, out%x)
end function smooth_spline

pure function predict_smooth_spline(fit, xnew) result(out)
! Evaluate distribution helper predict_smooth_spline.
type(smooth_spline_fit_t), intent(in) :: fit ! input value
real(kind=dp), intent(in) :: xnew(:)
type(smooth_xy_t) :: out
integer :: i, j
real(kind=dp) :: t
out%x = xnew
out%df = fit%df
allocate(out%y(size(xnew)))
if (.not. allocated(fit%x) .or. .not. allocated(fit%y)) then
   out%y = ieee_value(1.0_dp, ieee_quiet_nan)
   return
end if
if (size(fit%x) <= 0 .or. size(fit%y) <= 0) then
   out%y = ieee_value(1.0_dp, ieee_quiet_nan)
   return
end if
do i = 1, size(xnew)
   if (xnew(i) <= fit%x(1) .or. size(fit%x) == 1) then
      out%y(i) = fit%y(1)
   else if (xnew(i) >= fit%x(size(fit%x))) then
      out%y(i) = fit%y(size(fit%y))
   else
      j = 1
      do while (j < size(fit%x) - 1 .and. fit%x(j + 1) < xnew(i))
         j = j + 1
      end do
      t = (xnew(i) - fit%x(j)) / max(tiny(1.0_dp), fit%x(j + 1) - fit%x(j))
      out%y(i) = (1.0_dp - t) * fit%y(j) + t * fit%y(j + 1)
   end if
end do
end function predict_smooth_spline

function rbinom_scalar(n, size_, prob) result(x)
! Return n binomial(size, prob) variates.
integer, intent(in) :: n ! item count
integer, intent(in) :: size_ ! integer argument
real(kind=dp), intent(in) :: prob
integer, allocatable :: x(:)
integer :: i, j, s
real(kind=dp) :: u, p
allocate(x(max(0, n)))
if (size_ < 0 .or. prob < 0.0_dp .or. prob > 1.0_dp .or. prob /= prob) then
   x = -huge(0)
   return
end if
p = prob
do i = 1, size(x)
   if (size_ <= 0 .or. p <= 0.0_dp) then
      x(i) = 0
      cycle
   else if (p >= 1.0_dp) then
      x(i) = size_
      cycle
   end if
   s = 0
   do j = 1, max(0, size_)
   u = runif1()
      if (u < p) s = s + 1
   end do
   x(i) = s
end do
end function rbinom_scalar

function rbinom_vector_prob(n, size_, prob) result(x)
! Return n binomial(size, prob(i)) variates.
integer, intent(in) :: n ! item count
integer, intent(in) :: size_ ! integer argument
real(kind=dp), intent(in) :: prob(:)
integer, allocatable :: x(:)
integer :: i, j, s
real(kind=dp) :: u, p
allocate(x(max(0, n)))
if (size_ < 0) then
   x = -huge(0)
   return
end if
do i = 1, size(x)
   if (size(prob) > 0) then
      p = prob(min(i, size(prob)))
   else
      p = 0.0_dp
   end if
   if (p < 0.0_dp .or. p > 1.0_dp .or. p /= p) then
      x(i) = -huge(0)
      cycle
   end if
   if (size_ <= 0 .or. p <= 0.0_dp) then
      x(i) = 0
      cycle
   else if (p >= 1.0_dp) then
      x(i) = size_
      cycle
   end if
   s = 0
   do j = 1, max(0, size_)
   u = runif1()
      if (u < p) s = s + 1
   end do
   x(i) = s
end do
end function rbinom_vector_prob

function rpois_scalar(n, lambda) result(x)
! Return n Poisson(lambda) variates.
integer, intent(in) :: n ! item count
real(kind=dp), intent(in) :: lambda
integer, allocatable :: x(:)
integer :: i
allocate(x(max(0, n)))
do i = 1, size(x)
   x(i) = rpois_one(lambda)
end do
end function rpois_scalar

function rpois_vector(n, lambda) result(x)
! Return n Poisson(lambda(i)) variates.
integer, intent(in) :: n ! item count
real(kind=dp), intent(in) :: lambda(:)
integer, allocatable :: x(:)
integer :: i
allocate(x(max(0, n)))
do i = 1, size(x)
   if (size(lambda) > 0) then
      x(i) = rpois_one(lambda(min(i, size(lambda))))
   else
      x(i) = 0
   end if
end do
end function rpois_vector

function rpois_one(lambda) result(k)
! Draw one Poisson variate; Knuth for small lambda, normal approximation for large lambda.
real(kind=dp), intent(in) :: lambda
integer :: k
real(kind=dp) :: lprob, p, u, z
if (lambda < 0.0_dp .or. lambda /= lambda .or. .not. ieee_is_finite(lambda)) then
   k = -huge(0)
   return
else if (lambda == 0.0_dp) then
   k = 0
   return
end if
if (lambda < 40.0_dp) then
   lprob = exp(-lambda)
   k = 0
   p = 1.0_dp
   do
   u = runif1()
      p = p * u
      if (p <= lprob) exit
      k = k + 1
   end do
else
   z = rnorm1()
   k = max(0, nint(lambda + sqrt(lambda) * z))
end if
end function rpois_one

function random_choice2_prob(n, p1) result(z)
! Sample n labels in {1,2} with P(label=1)=p1.
integer, intent(in) :: n ! item count
real(kind=dp), intent(in) :: p1
integer, allocatable :: z(:)
integer :: i
real(kind=dp) :: u
allocate(z(n))
do i = 1, n
   u = runif1()
   z(i) = merge(1, 2, u < p1)
end do
end function random_choice2_prob

function randint_range(n, lo, hi) result(out)
! Sample n integers uniformly from [lo, hi].
integer, intent(in) :: n, lo, hi
integer, allocatable :: out(:)
integer :: i, span
real(kind=dp) :: u
if (hi < lo) then
   allocate(out(0))
   return
end if
span = hi - lo + 1
allocate(out(n))
do i = 1, n
   u = runif1()
   out(i) = lo + int(u * span)
   if (out(i) > hi) out(i) = hi
end do
end function randint_range

function sample_int(n, size_, replace, prob) result(out)
! R-like sample.int with optional replacement and probabilities.
integer, intent(in) :: n ! item count
integer, intent(in), optional :: size_ ! integer argument
logical, intent(in), optional :: replace ! logical flag
real(kind=dp), intent(in), optional :: prob(:)
integer, allocatable :: out(:), pool(:)
real(kind=dp), allocatable :: w(:), cdf(:)
integer :: m, i, j, pick, tmp
logical :: rep
real(kind=dp) :: u, s, acc

m = n
if (present(size_)) m = size_
rep = .false.
if (present(replace)) rep = replace

if (n < 1) then
   if (m == 0) then
      allocate(out(0))
      return
   end if
   error stop "sample_int: n must be >= 1"
end if
if (m < 0) error stop "sample_int: size must be >= 0"
if ((.not. rep) .and. (m > n)) error stop &
   & "sample_int: size > n without replacement"
if (present(prob)) then
   if (size(prob) /= n) error stop "sample_int: prob length mismatch"
   if (any(prob < 0.0_dp)) error stop "sample_int: prob must be nonnegative"
   if (sum(prob) <= 0.0_dp) error stop "sample_int: prob sum must be positive"
end if

allocate(out(m))
if (m == 0) return

if (present(prob)) then
   if (rep) then
      allocate(cdf(n))
      s = sum(prob)
      cdf(1) = prob(1) / s
      do i = 2, n
         cdf(i) = cdf(i - 1) + prob(i) / s
      end do
      cdf(n) = 1.0_dp
      do i = 1, m
   u = runif1()
         pick = 1
         do while (pick < n .and. u > cdf(pick))
            pick = pick + 1
         end do
         out(i) = pick
      end do
   else
      allocate(w(n))
      w = prob
      do i = 1, m
         s = sum(w)
         if (s <= 0.0_dp) error stop "sample_int: depleted probability mass"
   u = runif1()
         u = u * s
         acc = 0.0_dp
         pick = n
         do j = 1, n
            acc = acc + w(j)
            if (u <= acc) then
               pick = j
               exit
            end if
         end do
         out(i) = pick
         w(pick) = 0.0_dp
      end do
   end if
else
   if (rep) then
      do i = 1, m
   u = runif1()
         out(i) = 1 + int(u * real(n, kind=dp))
         if (out(i) > n) out(i) = n
      end do
   else
      allocate(pool(n))
      pool = [(j, j=1,n)]
      do i = 1, m
   u = runif1()
         pick = i + int(u * real(n - i + 1, kind=dp))
         if (pick > n) pick = n
         tmp = pool(i)
         pool(i) = pool(pick)
         pool(pick) = tmp
         out(i) = pool(i)
      end do
   end if
end if
end function sample_int

function sample_int1(n, replace, prob) result(out)
! Scalar wrapper for sample_int(..., size_=1).
integer, intent(in) :: n ! item count
logical, intent(in), optional :: replace ! logical flag
real(kind=dp), intent(in), optional :: prob(:)
integer :: out
integer, allocatable :: tmp(:)
if (present(replace) .and. present(prob)) then
   tmp = sample_int(n, size_=1, replace=replace, prob=prob)
else if (present(prob)) then
   tmp = sample_int(n, size_=1, prob=prob)
else if (present(replace)) then
    tmp = sample_int(n, size_=1, replace=replace)
else
   tmp = sample_int(n, size_=1)
end if
out = tmp(1)
end function sample_int1

function sample_value_int(x, replace) result(out)
integer, intent(in) :: x(:)
logical, intent(in), optional :: replace
integer :: out
logical :: replace_use
replace_use = .false.
if (present(replace)) replace_use = replace
out = x(sample_int1(size(x), replace=replace_use))
end function sample_value_int

pure function combn_indices(n, m) result(out)
! Return lexicographically ordered column indices for all m-combinations of n items.
integer, intent(in) :: n, m
integer, allocatable :: out(:,:)
integer(kind=int64) :: count64
integer, allocatable :: idx(:)
integer :: i, j, ncomb

if (m < 0 .or. m > n) error stop "combn: m must be between 0 and length(x)"
count64 = 1_int64
do i = 1, min(m, n - m)
   count64 = count64 * int(n - min(m, n - m) + i, kind=int64) / int(i, kind=int64)
end do
if (count64 > int(huge(ncomb), kind=int64)) error stop "combn: result has too many columns"
ncomb = int(count64)
allocate(out(m, ncomb))
if (m == 0) return
allocate(idx(m))
idx = [(i, i=1,m)]
do j = 1, ncomb
   out(:, j) = idx
   i = m
   do while (i >= 1)
      if (idx(i) /= n - m + i) exit
      i = i - 1
   end do
   if (i == 0) exit
   idx(i) = idx(i) + 1
   do while (i < m)
      idx(i + 1) = idx(i) + 1
      i = i + 1
   end do
end do
end function combn_indices

pure function combn_int(x, m) result(out)
! Return all size-m combinations of an integer vector as matrix columns.
integer, intent(in) :: x(:), m
integer, allocatable :: out(:,:), indices(:,:)
indices = combn_indices(size(x), m)
out = reshape(x(reshape(indices, [size(indices)])), shape(indices))
end function combn_int

pure function combn_real(x, m) result(out)
! Return all size-m combinations of a real vector as matrix columns.
real(kind=dp), intent(in) :: x(:)
integer, intent(in) :: m
real(kind=dp), allocatable :: out(:,:)
integer, allocatable :: indices(:,:)
indices = combn_indices(size(x), m)
out = reshape(x(reshape(indices, [size(indices)])), shape(indices))
end function combn_real

pure function combn_char(x, m) result(out)
! Return all size-m combinations of a character vector as matrix columns.
character(len=*), intent(in) :: x(:)
integer, intent(in) :: m
character(len=:), allocatable :: out(:,:)
integer, allocatable :: indices(:,:)
indices = combn_indices(size(x), m)
allocate(character(len=len(x)) :: out(size(indices, 1), size(indices, 2)))
out = reshape(x(reshape(indices, [size(indices)])), shape(indices))
end function combn_char

pure function fivenum(x) result(out)
! R-like fivenum(): Tukey five-number summary
! (minimum, lower hinge, median, upper hinge, maximum).
real(kind=dp), intent(in) :: x(:)
real(kind=dp) :: out(5)
real(kind=dp), allocatable :: xs(:)
real(kind=dp) :: n4, d(5), rn
integer :: n, k
xs = pack(x, x == x)
n = size(xs)
if (n == 0) then
   out = ieee_value(0.0_dp, ieee_quiet_nan)
   return
end if
call sort_increasing(xs)
rn = real(n, kind=dp)
n4 = floor((rn + 3.0_dp) / 2.0_dp) / 2.0_dp
d = [1.0_dp, n4, (rn + 1.0_dp) / 2.0_dp, rn + 1.0_dp - n4, rn]
do k = 1, 5
   out(k) = 0.5_dp * (xs(floor(d(k))) + xs(ceiling(d(k))))
end do
end function fivenum

pure subroutine sort_increasing(x)
! Sort a real vector in increasing order (insertion sort).
real(kind=dp), intent(inout) :: x(:) ! input vector
integer :: i, j
real(kind=dp) :: key
do i = 2, size(x)
   key = x(i)
   j = i - 1
   do
      if (j < 1) exit
      if (x(j) <= key) exit
      x(j + 1) = x(j)
      j = j - 1
   end do
   x(j + 1) = key
end do
end subroutine sort_increasing

pure logical function real_order_precedes(a, b, decreasing) result(out)
! Compare real values using R's stable ordering, with all missing values at the end.
real(kind=dp), intent(in) :: a, b
logical, intent(in) :: decreasing
integer :: a_missing, b_missing

if (a == a) then
   a_missing = 0
else
   a_missing = 1
end if
if (b == b) then
   b_missing = 0
else
   b_missing = 1
end if

if (a_missing /= b_missing) then
   out = a_missing < b_missing
else if (a_missing /= 0) then
   out = .false.
else if (decreasing) then
   out = a > b
else
   out = a < b
end if
end function real_order_precedes

pure subroutine sort_increasing_int(x)
! Sort an integer vector in increasing order (insertion sort).
integer, intent(inout) :: x(:) ! input vector
integer :: i, j, key
do i = 2, size(x)
   key = x(i)
   j = i - 1
   do
      if (j < 1) exit
      if (x(j) <= key) exit
      x(j + 1) = x(j)
      j = j - 1
   end do
   x(j + 1) = key
end do
end subroutine sort_increasing_int

pure logical function int_order_precedes(a, b, decreasing) result(out)
! Compare integer values using R's stable ordering, with integer NA at the end.
integer, intent(in) :: a, b
logical, intent(in) :: decreasing
logical :: a_missing, b_missing

a_missing = a == -huge(0)
b_missing = b == -huge(0)
if (a_missing .neqv. b_missing) then
   out = .not. a_missing
else if (a_missing) then
   out = .false.
else if (decreasing) then
   out = a > b
else
   out = a < b
end if
end function int_order_precedes

pure function sort_real(x, decreasing) result(out)
! Return a sorted copy of a real vector, omitting missing values like R sort().
real(kind=dp), intent(in) :: x(:) ! input vector
logical, intent(in), optional :: decreasing
real(kind=dp), allocatable :: out(:)
integer, allocatable :: idx(:)
logical :: dec

dec = .false.
if (present(decreasing)) dec = decreasing
out = pack(x, x == x)
idx = sort_list_real(out, dec)
out = out(idx)
end function sort_real

pure function sort_int(x, decreasing) result(out)
! Return a sorted copy of an integer vector, omitting integer NA like R sort().
integer, intent(in) :: x(:) ! input vector
logical, intent(in), optional :: decreasing
integer, allocatable :: out(:)
integer, allocatable :: idx(:)
logical :: dec

dec = .false.
if (present(decreasing)) dec = decreasing
out = pack(x, x /= -huge(0))
idx = sort_list_int(out, dec)
out = out(idx)
end function sort_int

pure function sort_char(x, decreasing) result(out)
! Return a sorted copy of a character vector.
character(len=*), intent(in) :: x(:) ! input vector
logical, intent(in), optional :: decreasing
character(len=:), allocatable :: out(:)
integer, allocatable :: idx(:)
idx = sort_list_char(x, decreasing)
allocate(character(len=len(x)) :: out(size(x)))
out = x(idx)
end function sort_char

pure function sort_list_real(x, decreasing) result(idx)
! Return 1-based indices that sort a real vector.
real(kind=dp), intent(in) :: x(:) ! input vector
logical, intent(in), optional :: decreasing
integer, allocatable :: idx(:)
integer :: i, j, t
logical :: dec
allocate(idx(size(x)))
do i = 1, size(x)
   idx(i) = i
end do
dec = .false.
if (present(decreasing)) dec = decreasing
do i = 2, size(idx)
   t = idx(i)
   j = i - 1
   do
      if (j < 1) exit
      if (.not. real_order_precedes(x(t), x(idx(j)), dec)) exit
      idx(j + 1) = idx(j)
      j = j - 1
   end do
   idx(j + 1) = t
end do
end function sort_list_real

pure function sort_list_int(x, decreasing) result(idx)
! Return 1-based indices that sort an integer vector.
integer, intent(in) :: x(:) ! input vector
logical, intent(in), optional :: decreasing
integer, allocatable :: idx(:)
integer :: i, j, t
logical :: dec
allocate(idx(size(x)))
do i = 1, size(x)
   idx(i) = i
end do
dec = .false.
if (present(decreasing)) dec = decreasing
do i = 2, size(idx)
   t = idx(i)
   j = i - 1
   do
      if (j < 1) exit
      if (.not. int_order_precedes(x(t), x(idx(j)), dec)) exit
      idx(j + 1) = idx(j)
      j = j - 1
   end do
   idx(j + 1) = t
end do
end function sort_list_int

pure function sort_list_char(x, decreasing) result(idx)
! Return 1-based indices that sort a character vector.
character(len=*), intent(in) :: x(:) ! input vector
logical, intent(in), optional :: decreasing
integer, allocatable :: idx(:)
integer :: i, j, t
logical :: dec
allocate(idx(size(x)))
do i = 1, size(x)
   idx(i) = i
end do
dec = .false.
if (present(decreasing)) dec = decreasing
do i = 2, size(idx)
   t = idx(i)
   j = i - 1
   if (dec) then
      do
         if (j < 1) exit
         if (x(idx(j)) >= x(t)) exit
         idx(j + 1) = idx(j)
         j = j - 1
      end do
   else
      do
         if (j < 1) exit
         if (x(idx(j)) <= x(t)) exit
         idx(j + 1) = idx(j)
         j = j - 1
      end do
   end if
   idx(j + 1) = t
end do
end function sort_list_char

pure function r_head_real(x, n) result(out)
! Return the first n elements of a real vector.
real(kind=dp), intent(in) :: x(:) ! input vector
integer, intent(in) :: n
real(kind=dp), allocatable :: out(:)
integer :: m
if (n < 0) then
   m = max(0, size(x) + n)
else
   m = min(n, size(x))
end if
allocate(out(m))
if (m > 0) out = x(1:m)
end function r_head_real

pure function r_head_int(x, n) result(out)
! Return the first n elements of an integer vector.
integer, intent(in) :: x(:) ! input vector
integer, intent(in) :: n
integer, allocatable :: out(:)
integer :: m
if (n < 0) then
   m = max(0, size(x) + n)
else
   m = min(n, size(x))
end if
allocate(out(m))
if (m > 0) out = x(1:m)
end function r_head_int

pure function r_head_real_mat(x, n) result(out)
! Return the first n rows of a real matrix.
real(kind=dp), intent(in) :: x(:,:) ! input matrix
integer, intent(in) :: n
real(kind=dp), allocatable :: out(:,:)
integer :: m
if (n < 0) then
   m = max(0, size(x, 1) + n)
else
   m = min(n, size(x, 1))
end if
allocate(out(m, max(0, size(x, 2))))
if (m > 0) out = x(1:m, :)
end function r_head_real_mat

pure function r_head_int_mat(x, n) result(out)
! Return the first n rows of an integer matrix.
integer, intent(in) :: x(:,:) ! input matrix
integer, intent(in) :: n
integer, allocatable :: out(:,:)
integer :: m
if (n < 0) then
   m = max(0, size(x, 1) + n)
else
   m = min(n, size(x, 1))
end if
allocate(out(m, max(0, size(x, 2))))
if (m > 0) out = x(1:m, :)
end function r_head_int_mat

pure function order_real(x) result(idx)
! Return 1-based order indices that sort a real vector increasingly.
real(kind=dp), intent(in) :: x(:)
integer, allocatable :: idx(:)
idx = sort_list_real(x)
end function order_real

pure function rank_first(x) result(out)
! Return R rank(x, ties.method="first") for a real vector.
real(kind=dp), intent(in) :: x(:)
real(kind=dp), allocatable :: out(:)
integer, allocatable :: ord(:)
integer :: i
allocate(out(size(x)))
if (size(x) <= 0) return
ord = order_real(x)
do i = 1, size(ord)
   out(ord(i)) = real(i, kind=dp)
end do
end function rank_first

pure function rank_average(x) result(out)
! Return R rank(x, ties.method="average") for a real vector.
real(kind=dp), intent(in) :: x(:)
real(kind=dp), allocatable :: out(:)
integer, allocatable :: ord(:)
integer :: first, i, last, nvalid
real(kind=dp) :: r
allocate(out(size(x)))
if (size(x) <= 0) return
ord = order_real(x)
nvalid = count(x == x)
first = 1
do while (first <= nvalid)
   last = first
   do
      if (last >= nvalid) exit
      if (.not. real_values_equal(x(ord(last + 1)), x(ord(first)))) exit
      last = last + 1
   end do
   r = 0.5_dp * real(first + last, kind=dp)
   do i = first, last
      out(ord(i)) = r
   end do
   first = last + 1
end do
do i = nvalid + 1, size(ord)
   out(ord(i)) = real(i, kind=dp)
end do
end function rank_average

pure function det_real_mat(x) result(out)
! Return the determinant of a real square matrix.
real(kind=dp), intent(in) :: x(:,:)
real(kind=dp) :: out
real(kind=dp), allocatable :: a(:)
integer :: i, j, k, n, p
real(kind=dp) :: fac, piv, t
n = size(x, 1)
out = 0.0_dp
if (n /= size(x, 2)) return
if (n == 0) then
   out = 1.0_dp
   return
end if
allocate(a(n*n))
a = reshape(x, [n*n])
out = 1.0_dp
do k = 1, n
   p = k
   piv = abs(a((k - 1)*n + k))
   do i = k + 1, n
      if (abs(a((i - 1)*n + k)) > piv) then
         p = i
         piv = abs(a((i - 1)*n + k))
      end if
   end do
   if (piv <= tiny(1.0_dp)) then
      out = 0.0_dp
      return
   end if
   if (p /= k) then
      do j = k, n
         t = a((k - 1)*n + j)
         a((k - 1)*n + j) = a((p - 1)*n + j)
         a((p - 1)*n + j) = t
      end do
      out = -out
   end if
   out = out * a((k - 1)*n + k)
   do i = k + 1, n
      fac = a((i - 1)*n + k) / a((k - 1)*n + k)
      do j = k + 1, n
         a((i - 1)*n + j) = a((i - 1)*n + j) - fac * a((k - 1)*n + j)
      end do
   end do
end do
end function det_real_mat

pure function det_real_int(x) result(out)
! Evaluate distribution helper det_real_int.
integer, intent(in) :: x(:,:)
real(kind=dp) :: out
out = det_real_mat(real(x, kind=dp))
end function det_real_int

pure function kappa_real(x) result(out)
! Return a simple 1-norm condition estimate for a square matrix.
real(kind=dp), intent(in) :: x(:,:)
real(kind=dp) :: out
real(kind=dp), allocatable :: inv_x(:,:)
integer :: n
n = size(x, 1)
out = huge(1.0_dp)
if (n <= 0 .or. size(x, 2) /= n) return
inv_x = solve_real_mat(x, real(diag(n), kind=dp))
if (all(inv_x == 0.0_dp)) return
out = maxval(sum(abs(x), dim=1)) * maxval(sum(abs(inv_x), dim=1))
end function kappa_real

pure function eigen_sym_values(x) result(vals)
! Return eigenvalues of a real symmetric matrix using Jacobi rotations.
real(kind=dp), intent(in) :: x(:,:)
real(kind=dp), allocatable :: vals(:)
real(kind=dp), allocatable :: a(:,:)
integer :: i, j, p, q, iter, n, max_iter
real(kind=dp) :: app, aqq, apq, c, off, phi, s, tau, t, tmp
n = size(x, 1)
allocate(vals(n))
vals = 0.0_dp
if (n <= 0 .or. n /= size(x, 2)) return
a = x
max_iter = max(1, 100 * n * n)
do iter = 1, max_iter
   off = 0.0_dp
   p = 1
   q = min(2, n)
   do i = 1, n - 1
      do j = i + 1, n
         if (abs(a(i, j)) > off) then
            off = abs(a(i, j))
            p = i
            q = j
         end if
      end do
   end do
   if (off <= 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(a)))) exit
   app = a(p, p)
   aqq = a(q, q)
   apq = a(p, q)
   if (abs(apq) <= tiny(1.0_dp)) cycle
   tau = (aqq - app) / (2.0_dp * apq)
   t = sign(1.0_dp, tau) / (abs(tau) + sqrt(1.0_dp + tau * tau))
   c = 1.0_dp / sqrt(1.0_dp + t * t)
   s = t * c
   do j = 1, n
      if (j /= p .and. j /= q) then
         tmp = a(j, p)
         a(j, p) = c * tmp - s * a(j, q)
         a(p, j) = a(j, p)
         a(j, q) = s * tmp + c * a(j, q)
         a(q, j) = a(j, q)
      end if
   end do
   phi = t * apq
   a(p, p) = app - phi
   a(q, q) = aqq + phi
   a(p, q) = 0.0_dp
   a(q, p) = 0.0_dp
end do
do i = 1, n
   vals(i) = a(i, i)
end do
end function eigen_sym_values

function eigen_int(x, symmetric, only_values) result(fit)
! Runtime helper for R-compatible eigen int.
integer, intent(in) :: x(:,:) ! input matrix
logical, intent(in), optional :: symmetric, only_values
type(eigen_result_t) :: fit
if (present(symmetric) .and. present(only_values)) then
   fit = eigen_real(real(x, kind=dp), symmetric=symmetric, only_values=only_values)
else if (present(symmetric)) then
   fit = eigen_real(real(x, kind=dp), symmetric=symmetric)
else if (present(only_values)) then
   fit = eigen_real(real(x, kind=dp), only_values=only_values)
else
   fit = eigen_real(real(x, kind=dp))
end if
end function eigen_int

function eigen_real(x, symmetric, only_values) result(fit)
! Small real eigen decomposition helper for R eigen() examples.
real(kind=dp), intent(in) :: x(:,:) ! input matrix
logical, intent(in), optional :: symmetric, only_values
type(eigen_result_t) :: fit
real(kind=dp), allocatable :: a(:,:), vecs(:,:), tmpv(:), qmat(:,:), rmat(:,:), v(:), w(:), bmat(:,:), coeff(:), &
& bwork(:,:)
complex(kind=dp), allocatable :: tmpcv(:)
real(kind=dp) :: off, app, aqq, apq, tau, t, c, s, phi, tmp, normv, rkk, delta, mu
real(kind=dp) :: aa, bb, cc, dd, tr, disc, root, lam
real(kind=dp) :: bound, left, right, mid, fleft, fright, fmid, xval, root_tol
complex(kind=dp) :: clam, ctmp
integer :: n, i, j, k, q, iter, max_iter, imax, m, root_count, grid, ig
logical :: do_symmetric, only_vals
n = size(x, 1)
allocate(fit%values(n), fit%vectors(n, n))
fit%values = 0.0_dp
fit%vectors = 0.0_dp
do i = 1, n
   fit%vectors(i, i) = 1.0_dp
end do
if (n <= 0 .or. n /= size(x, 2)) return
do_symmetric = maxval(abs(x - transpose(x))) <= 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(x)))
if (present(symmetric)) do_symmetric = symmetric
only_vals = .false.
if (present(only_values)) only_vals = only_values
if ((.not. do_symmetric) .and. n > 2) then
   a = x
   allocate(qmat(n, n), rmat(n, n), v(n), w(n), bmat(n, n), tmpv(n))
   max_iter = max(1, 4000 * n * n)
   do m = n, 2, -1
      do iter = 1, max_iter
         off = 1000.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(a(1:m, 1:m))))
         if (abs(a(m, m - 1)) <= off) exit
         aa = a(m - 1, m - 1)
         bb = a(m - 1, m)
         cc = a(m, m - 1)
         dd = a(m, m)
         tr = aa + dd
         disc = (aa - dd) * (aa - dd) + 4.0_dp * bb * cc
         if (disc >= 0.0_dp) then
            root = sqrt(disc)
            if (abs(0.5_dp * (tr + root) - dd) <= abs(0.5_dp * (tr - root) - dd)) then
               mu = 0.5_dp * (tr + root)
            else
               mu = 0.5_dp * (tr - root)
            end if
         else
            mu = dd
         end if
         bmat(1:m, 1:m) = a(1:m, 1:m)
         do i = 1, m
            bmat(i, i) = bmat(i, i) - mu
         end do
         qmat = 0.0_dp
         rmat = 0.0_dp
         do j = 1, m
            v = 0.0_dp
            v(1:m) = bmat(1:m, j)
            do k = 1, j - 1
               rmat(k, j) = sum(qmat(1:m, k) * v(1:m))
               v(1:m) = v(1:m) - rmat(k, j) * qmat(1:m, k)
            end do
            rkk = sqrt(sum(v(1:m) * v(1:m)))
            if (rkk <= 100.0_dp * tiny(1.0_dp)) then
               v = 0.0_dp
               v(j) = 1.0_dp
               do k = 1, j - 1
                  tmp = sum(qmat(1:m, k) * v(1:m))
                  v(1:m) = v(1:m) - tmp * qmat(1:m, k)
               end do
               rkk = sqrt(max(tiny(1.0_dp), sum(v(1:m) * v(1:m))))
            end if
            qmat(1:m, j) = v(1:m) / rkk
            rmat(j, j) = rkk
         end do
         a(1:m, 1:m) = matmul(rmat(1:m, 1:m), qmat(1:m, 1:m))
         do i = 1, m
            a(i, i) = a(i, i) + mu
         end do
      end do
      a(m, m - 1) = 0.0_dp
   end do
   do i = 1, n
      fit%values(i) = a(i, i)
   end do
   allocate(coeff(n + 1), bwork(n, n))
   coeff = 0.0_dp
   coeff(1) = 1.0_dp
   bwork = 0.0_dp
   do i = 1, n
      bwork(i, i) = 1.0_dp
   end do
   do k = 1, n
      bmat = matmul(x, bwork)
      tmp = 0.0_dp
      do i = 1, n
         tmp = tmp + bmat(i, i)
      end do
      coeff(k + 1) = -tmp / real(k, kind=dp)
      bwork = bmat
      do i = 1, n
         bwork(i, i) = bwork(i, i) + coeff(k + 1)
      end do
   end do
   bound = 1.0_dp + maxval(sum(abs(x), dim=2))
   grid = max(2000, 4000 * n)
   root_tol = 1000.0_dp * epsilon(1.0_dp) * bound
   root_count = 0
   left = -bound
   fleft = coeff(1)
   do k = 2, n + 1
      fleft = fleft * left + coeff(k)
   end do
   do ig = 1, grid
      right = -bound + 2.0_dp * bound * real(ig, kind=dp) / real(grid, kind=dp)
      fright = coeff(1)
      do k = 2, n + 1
         fright = fright * right + coeff(k)
      end do
      if (abs(fleft) <= root_tol .or. fleft * fright < 0.0_dp .or. abs(fright) <= root_tol) then
         if (abs(fleft) <= root_tol) then
            root = left
         else if (abs(fright) <= root_tol) then
            root = right
         else
            aa = left
            bb = right
            do iter = 1, 80
               mid = 0.5_dp * (aa + bb)
               fmid = coeff(1)
               do k = 2, n + 1
                  fmid = fmid * mid + coeff(k)
               end do
               if (abs(fmid) <= root_tol) exit
               if (fleft * fmid <= 0.0_dp) then
                  bb = mid
                  fright = fmid
               else
                  aa = mid
                  fleft = fmid
               end if
            end do
            root = mid
         end if
         if (root_count == 0 .or. minval(abs(fit%values(1:root_count) - root)) > 100.0_dp * root_tol) then
            root_count = root_count + 1
            if (root_count <= n) fit%values(root_count) = root
         end if
      end if
      left = right
      fleft = fright
   end do
   if (root_count > 0 .and. root_count < n) then
      do i = root_count + 1, n
         fit%values(i) = a(i, i)
      end do
   end if
   do i = 1, n - 1
      imax = i
      do j = i + 1, n
         if (abs(fit%values(j)) > abs(fit%values(imax))) imax = j
      end do
      if (imax /= i) then
         ctmp = fit%values(i)
         fit%values(i) = fit%values(imax)
         fit%values(imax) = ctmp
      end if
   end do
   do j = 1, n
      lam = real(fit%values(j), kind=dp)
      do i = 1, n
         v(i) = 1.0_dp + real(i + j - 2, kind=dp) / real(max(1, n), kind=dp)
      end do
      normv = sqrt(max(tiny(1.0_dp), sum(v * v)))
      v = v / normv
      do iter = 1, 80
         bmat = x
         delta = sqrt(epsilon(1.0_dp)) * max(1.0_dp, abs(lam), maxval(abs(x)))
         do i = 1, n
            bmat(i, i) = bmat(i, i) - lam + delta
         end do
         w = solve_real_vec(bmat, v)
         normv = sqrt(max(tiny(1.0_dp), sum(w * w)))
         if (normv <= 100.0_dp * tiny(1.0_dp)) exit
         v = w / normv
         w = matmul(x, v)
         lam = sum(v * w) / max(tiny(1.0_dp), sum(v * v))
      end do
      fit%values(j) = lam
      fit%vectors(:, j) = v
      if (real(fit%vectors(1, j), kind=dp) < 0.0_dp) fit%vectors(:, j) = -fit%vectors(:, j)
   end do
   do i = 1, n - 1
      imax = i
      do j = i + 1, n
         if (abs(fit%values(j)) > abs(fit%values(imax))) imax = j
      end do
      if (imax /= i) then
         ctmp = fit%values(i)
         fit%values(i) = fit%values(imax)
         fit%values(imax) = ctmp
         if (.not. allocated(tmpcv)) allocate(tmpcv(n))
         tmpcv = fit%vectors(:, i)
         fit%vectors(:, i) = fit%vectors(:, imax)
         fit%vectors(:, imax) = tmpcv
      end if
   end do
   if (only_vals) then
      if (allocated(fit%vectors)) deallocate(fit%vectors)
      allocate(fit%vectors(0, 0))
   end if
   return
end if
if ((.not. do_symmetric) .and. n == 2) then
   aa = x(1, 1)
   bb = x(1, 2)
   cc = x(2, 1)
   dd = x(2, 2)
   tr = aa + dd
   disc = (aa - dd) * (aa - dd) + 4.0_dp * bb * cc
   if (disc >= 0.0_dp) then
      root = sqrt(disc)
      fit%values = [(0.5_dp * (tr + root)), (0.5_dp * (tr - root))]
      do j = 1, 2
         lam = real(fit%values(j), kind=dp)
         if (abs(bb) >= abs(cc) .and. abs(bb) > tiny(1.0_dp)) then
            fit%vectors(:, j) = [bb, lam - aa]
         else if (abs(cc) > tiny(1.0_dp)) then
            fit%vectors(:, j) = [lam - dd, cc]
         else
            fit%vectors(:, j) = 0.0_dp
            fit%vectors(j, j) = 1.0_dp
         end if
         normv = sqrt(max(tiny(1.0_dp), sum(abs(fit%vectors(:, j))**2)))
         fit%vectors(:, j) = fit%vectors(:, j) / normv
      end do
      return
   end if
   ! The complex conjugate eigenpair path follows the same real/imaginary
   ! eigenvalue convention used by Alan Miller's complex_eigen module, adapted
   ! from the NSWC/EISPACK routines.
   root = sqrt(-disc)
   fit%values = [cmplx(0.5_dp * tr, 0.5_dp * root, kind=dp), cmplx(0.5_dp * tr, -0.5_dp * root, kind=dp)]
   do j = 1, 2
      clam = fit%values(j)
      if (abs(bb) >= abs(cc) .and. abs(bb) > tiny(1.0_dp)) then
         fit%vectors(:, j) = [cmplx(bb, 0.0_dp, kind=dp), clam - cmplx(aa, 0.0_dp, kind=dp)]
      else if (abs(cc) > tiny(1.0_dp)) then
         fit%vectors(:, j) = [clam - cmplx(dd, 0.0_dp, kind=dp), cmplx(cc, 0.0_dp, kind=dp)]
      else
         fit%vectors(:, j) = cmplx(0.0_dp, 0.0_dp, kind=dp)
         fit%vectors(j, j) = cmplx(1.0_dp, 0.0_dp, kind=dp)
      end if
      normv = sqrt(max(tiny(1.0_dp), sum(abs(fit%vectors(:, j))**2)))
      fit%vectors(:, j) = fit%vectors(:, j) / normv
   end do
   if (only_vals) then
      if (allocated(fit%vectors)) deallocate(fit%vectors)
      allocate(fit%vectors(0, 0))
   end if
   return
end if
a = x
allocate(vecs(n, n))
vecs = 0.0_dp
do i = 1, n
   vecs(i, i) = 1.0_dp
end do
max_iter = max(1, 100 * n * n)
do iter = 1, max_iter
   off = 0.0_dp
   i = 1
   q = min(2, n)
   do j = 1, n - 1
      do k = j + 1, n
         if (abs(a(j, k)) > off) then
            off = abs(a(j, k))
            i = j
            q = k
         end if
      end do
   end do
   if (off <= 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(a)))) exit
   app = a(i, i)
   aqq = a(q, q)
   apq = a(i, q)
   if (abs(apq) <= tiny(1.0_dp)) cycle
   tau = (aqq - app) / (2.0_dp * apq)
   t = sign(1.0_dp, tau) / (abs(tau) + sqrt(1.0_dp + tau * tau))
   c = 1.0_dp / sqrt(1.0_dp + t * t)
   s = t * c
   do k = 1, n
      if (k /= i .and. k /= q) then
         tmp = a(k, i)
         a(k, i) = c * tmp - s * a(k, q)
         a(i, k) = a(k, i)
         a(k, q) = s * tmp + c * a(k, q)
         a(q, k) = a(k, q)
      end if
      tmp = vecs(k, i)
      vecs(k, i) = c * tmp - s * vecs(k, q)
      vecs(k, q) = s * tmp + c * vecs(k, q)
   end do
   phi = t * apq
   a(i, i) = app - phi
   a(q, q) = aqq + phi
   a(i, q) = 0.0_dp
   a(q, i) = 0.0_dp
end do
do i = 1, n
   fit%values(i) = a(i, i)
end do
allocate(tmpv(n))
do i = 1, n - 1
   imax = i
   do j = i + 1, n
      if (real(fit%values(j), kind=dp) > real(fit%values(imax), kind=dp)) imax = j
   end do
   if (imax /= i) then
      ctmp = fit%values(i)
      fit%values(i) = fit%values(imax)
      fit%values(imax) = ctmp
      tmpv = vecs(:, i)
      vecs(:, i) = vecs(:, imax)
      vecs(:, imax) = tmpv
   end if
end do
do j = 1, n
   normv = sqrt(max(tiny(1.0_dp), sum(vecs(:, j)**2)))
   fit%vectors(:, j) = vecs(:, j) / normv
   if (real(fit%vectors(1, j), kind=dp) < 0.0_dp) fit%vectors(:, j) = -fit%vectors(:, j)
end do
if (only_vals) then
   if (allocated(fit%vectors)) deallocate(fit%vectors)
   allocate(fit%vectors(0, 0))
end if
end function eigen_real

subroutine print_eigen(fit)
! Print eigen values in an R-like format.
type(eigen_result_t), intent(in) :: fit ! input value
write(*,'(a)') "$values"
if (allocated(fit%values) .and. size(fit%values) > 0 .and. maxval(abs(aimag(fit%values))) <= &
   100.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(real(fit%values, kind=dp))))) then
   call print_real_vector(real(fit%values, kind=dp))
else
   call print_complex_vector(fit%values)
end if
write(*,*)
write(*,'(a)') "$vectors"
if (allocated(fit%vectors) .and. size(fit%vectors, 1) > 0 .and. size(fit%vectors, 2) > 0) then
   if (maxval(abs(aimag(fit%vectors))) <= 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(real(fit%vectors, &
   & kind=dp))))) then
      call print_matrix_rstyle(real(fit%vectors, kind=dp))
   else
      call print_matrix(fit%vectors)
   end if
else
   write(*,'(a)') "NULL"
end if
end subroutine print_eigen

function prcomp(x, center, scale_) result(fit)
! Principal components for a numeric matrix using a symmetric Jacobi eigensolver.
real(kind=dp), intent(in) :: x(:,:) ! input matrix
logical, intent(in), optional :: center, scale_
type(prcomp_fit_t) :: fit
real(kind=dp), allocatable :: z(:,:), covm(:,:), vals(:), vecs(:,:), tmpv(:)
real(kind=dp) :: off, app, aqq, apq, tau, t, c, s, phi, tmp, normv
integer :: n, p, i, j, k, q, iter, max_iter, imax
logical :: do_center, do_scale
n = size(x, 1)
p = size(x, 2)
do_center = .true.
do_scale = .false.
if (present(center)) do_center = center
if (present(scale_)) do_scale = scale_
allocate(z(n, p), fit%center(p), fit%scale(p), covm(p, p), vals(p), vecs(p, p), tmpv(p))
z = x
fit%center = 0.0_dp
fit%scale = 1.0_dp
if (do_center .and. n > 0) then
   fit%center = sum(x, dim=1) / real(n, kind=dp)
   do j = 1, p
      z(:, j) = z(:, j) - fit%center(j)
   end do
end if
if (do_scale .and. n > 1) then
   do j = 1, p
      fit%scale(j) = sqrt(max(0.0_dp, sum(z(:, j)**2) / real(n - 1, kind=dp)))
      if (fit%scale(j) > 0.0_dp) z(:, j) = z(:, j) / fit%scale(j)
   end do
end if
if (n > 1) then
   covm = matmul(transpose(z), z) / real(n - 1, kind=dp)
else
   covm = 0.0_dp
end if
vecs = 0.0_dp
do i = 1, p
   vecs(i, i) = 1.0_dp
end do
max_iter = max(1, 100 * p * p)
do iter = 1, max_iter
   off = 0.0_dp
   i = 1
   q = min(2, p)
   do j = 1, p - 1
      do k = j + 1, p
         if (abs(covm(j, k)) > off) then
            off = abs(covm(j, k))
            i = j
            q = k
         end if
      end do
   end do
   if (off <= 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(covm)))) exit
   app = covm(i, i)
   aqq = covm(q, q)
   apq = covm(i, q)
   if (abs(apq) <= tiny(1.0_dp)) cycle
   tau = (aqq - app) / (2.0_dp * apq)
   t = sign(1.0_dp, tau) / (abs(tau) + sqrt(1.0_dp + tau * tau))
   c = 1.0_dp / sqrt(1.0_dp + t * t)
   s = t * c
   do k = 1, p
      if (k /= i .and. k /= q) then
         tmp = covm(k, i)
         covm(k, i) = c * tmp - s * covm(k, q)
         covm(i, k) = covm(k, i)
         covm(k, q) = s * tmp + c * covm(k, q)
         covm(q, k) = covm(k, q)
      end if
      tmp = vecs(k, i)
      vecs(k, i) = c * tmp - s * vecs(k, q)
      vecs(k, q) = s * tmp + c * vecs(k, q)
   end do
   phi = t * apq
   covm(i, i) = app - phi
   covm(q, q) = aqq + phi
   covm(i, q) = 0.0_dp
   covm(q, i) = 0.0_dp
end do
do i = 1, p
   vals(i) = covm(i, i)
end do
do i = 1, p - 1
   imax = i
   do j = i + 1, p
      if (vals(j) > vals(imax)) imax = j
   end do
   if (imax /= i) then
      tmp = vals(i)
      vals(i) = vals(imax)
      vals(imax) = tmp
      tmpv = vecs(:, i)
      vecs(:, i) = vecs(:, imax)
      vecs(:, imax) = tmpv
   end if
end do
do j = 1, p
   normv = sqrt(max(tiny(1.0_dp), sum(vecs(:, j)**2)))
   vecs(:, j) = vecs(:, j) / normv
   if (vecs(1, j) < 0.0_dp) vecs(:, j) = -vecs(:, j)
end do
fit%sdev = sqrt(max(0.0_dp, vals))
fit%rotation = vecs
fit%x = matmul(z, vecs)
end function prcomp

subroutine print_prcomp_summary(fit)
! Print prcomp summary values in an R-like format.
type(prcomp_fit_t), intent(in) :: fit ! input value
real(kind=dp), allocatable :: var(:), prop(:), cum(:)
integer :: j, p
p = size(fit%sdev)
allocate(var(p), prop(p), cum(p))
var = fit%sdev**2
if (sum(var) > 0.0_dp) then
   prop = var / sum(var)
else
   prop = 0.0_dp
end if
cum = 0.0_dp
do j = 1, p
   if (j == 1) then
      cum(j) = prop(j)
   else
      cum(j) = cum(j - 1) + prop(j)
   end if
end do
write(*,'(a)') "Importance of components:"
write(*,'(28x)', advance='no')
do j = 1, p
   write(*,'("PC",i0,8x)', advance='no') j
end do
write(*,*)
write(*,'(a28)', advance='no') "Standard deviation"
do j = 1, p
   write(*,'(f10.4,1x)', advance='no') fit%sdev(j)
end do
write(*,*)
write(*,'(a28)', advance='no') "Proportion of Variance"
do j = 1, p
   write(*,'(f10.4,1x)', advance='no') prop(j)
end do
write(*,*)
write(*,'(a28)', advance='no') "Cumulative Proportion"
do j = 1, p
   write(*,'(f10.4,1x)', advance='no') cum(j)
end do
write(*,*)
end subroutine print_prcomp_summary

function arima_sim_scalar(ar, ma, n) result(x)
! Compute R-like time-series helper arima_sim_scalar.
real(kind=dp), intent(in) :: ar ! input value
real(kind=dp), intent(in) :: ma ! input value
integer, intent(in) :: n
real(kind=dp), allocatable :: x(:)
real(kind=dp), allocatable :: z(:)
integer :: i
allocate(x(max(0, n)), source=0.0_dp)
if (n <= 0) return
z = rnorm_vec(n + 1)
x(1) = z(2) + ma * z(1)
do i = 2, n
   x(i) = ar * x(i - 1) + z(i + 1) + ma * z(i)
end do
end function arima_sim_scalar

function arima_sim_vector(ar, ma, n) result(x)
! Compute R-like time-series helper arima_sim_vector.
real(kind=dp), intent(in) :: ar(:) ! input vector
real(kind=dp), intent(in) :: ma ! input value
integer, intent(in) :: n
real(kind=dp), allocatable :: x(:)
real(kind=dp), allocatable :: z(:)
integer :: i, j, p
allocate(x(max(0, n)), source=0.0_dp)
if (n <= 0) return
p = size(ar)
z = rnorm_vec(n + 1)
do i = 1, n
   x(i) = z(i + 1)
   if (i > 1) x(i) = x(i) + ma * z(i)
   do j = 1, min(p, i - 1)
      x(i) = x(i) + ar(j) * x(i - j)
   end do
end do
end function arima_sim_vector

function arima_sim_vector_vector(ar, ma, n) result(x)
! Compute R-like ARMA simulation with vector AR and MA coefficients.
real(kind=dp), intent(in) :: ar(:) ! autoregressive coefficients
real(kind=dp), intent(in) :: ma(:) ! moving-average coefficients
integer, intent(in) :: n
real(kind=dp), allocatable :: x(:)
real(kind=dp), allocatable :: z(:)
integer :: i, j, p, q
allocate(x(max(0, n)), source=0.0_dp)
if (n <= 0) return
p = size(ar)
q = size(ma)
z = rnorm_vec(n + q)
do i = 1, n
   x(i) = z(i + q)
   do j = 1, min(q, i - 1)
      x(i) = x(i) + ma(j) * z(i + q - j)
   end do
   do j = 1, min(p, i - 1)
      x(i) = x(i) + ar(j) * x(i - j)
   end do
end do
end function arima_sim_vector_vector

pure function arima_fit(x, order, include_mean) result(fit)
! Compute R-like time-series helper arima_fit.
real(kind=dp), intent(in) :: x(:) ! input vector
integer, intent(in) :: order(:) ! input vector
logical, intent(in), optional :: include_mean
type(arima_fit_t) :: fit
real(kind=dp) :: resid(size(x)), best_resid(size(x))
real(kind=dp), allocatable :: design(:,:), y(:), beta(:), xtx(:,:), xty(:), params(:), trial(:)
real(kind=dp) :: rss, best_rss, step, intercept, candidate, trial_rss
integer :: n, k, i, pass, p_eff, q_eff, n_eff, j, coord, sgn, max_lag
logical :: include_mean_use
n = size(x)
fit%p = merge(order(1), 0, size(order) >= 1)
fit%d = merge(order(2), 0, size(order) >= 2)
fit%q = merge(order(3), 0, size(order) >= 3)
include_mean_use = .false.
if (present(include_mean)) include_mean_use = include_mean
allocate(fit%coef(max(1, fit%p + fit%q + merge(1, 0, include_mean_use))), source=0.0_dp)
if (n <= 0) return
fit%mean = sum(x) / real(n, kind=dp)
if (fit%q == 0) then
   p_eff = max(0, fit%p)
   n_eff = max(1, n - p_eff)
   if (p_eff == 0) then
      resid = x - fit%mean
      rss = sum(resid * resid)
      if (include_mean_use) fit%coef(size(fit%coef)) = fit%mean
   else
      allocate(design(n_eff, p_eff + 1), source=1.0_dp)
      allocate(y(n_eff), source=0.0_dp)
      do i = 1, n_eff
         y(i) = x(p_eff + i)
         do j = 1, p_eff
            design(i, j + 1) = x(p_eff + i - j)
         end do
      end do
      xtx = matmul(transpose(design), design)
      xty = matmul(transpose(design), y)
      do i = 1, size(xtx, 1)
         xtx(i, i) = xtx(i, i) + 1.0e-10_dp
      end do
      beta = solve_real(xtx, xty)
      intercept = beta(1)
      fit%coef(1:p_eff) = beta(2:p_eff + 1)
      if (abs(1.0_dp - sum(fit%coef(1:p_eff))) > 100.0_dp * epsilon(1.0_dp)) then
         fit%mean = intercept / (1.0_dp - sum(fit%coef(1:p_eff)))
      else
         fit%mean = sum(x) / real(n, kind=dp)
      end if
      if (include_mean_use .and. size(fit%coef) >= p_eff + 1) fit%coef(p_eff + 1) = fit%mean
      resid = 0.0_dp
      do i = p_eff + 1, n
         resid(i) = x(i) - intercept
         do j = 1, p_eff
            resid(i) = resid(i) - fit%coef(j) * x(i - j)
         end do
      end do
      rss = sum(resid(p_eff + 1:n) * resid(p_eff + 1:n))
   end if
   fit%resid = resid
   fit%last_x = x(n)
   fit%last_resid = resid(n)
   fit%sigma2 = max(rss / real(n_eff, kind=dp), tiny(1.0_dp))
   k = fit%p + merge(1, 0, include_mean_use)
   fit%aic = real(n_eff, kind=dp) * log(fit%sigma2) + 2.0_dp * real(k, kind=dp)
   return
end if
p_eff = max(0, fit%p)
q_eff = max(0, fit%q)
max_lag = max(p_eff, q_eff)
n_eff = max(1, n - max_lag)
allocate(params(p_eff + q_eff), trial(p_eff + q_eff), source=0.0_dp)
best_rss = huge(1.0_dp)
best_resid = 0.0_dp
step = 0.5_dp
do pass = 1, 10
   do coord = 1, size(params)
      do sgn = -1, 1, 2
         trial = params
         candidate = max(-0.98_dp, min(0.98_dp, params(coord) + real(sgn, kind=dp) * step))
         trial(coord) = candidate
         resid = 0.0_dp
         trial_rss = 0.0_dp
         do i = 1, n
            resid(i) = x(i) - fit%mean
            do j = 1, min(p_eff, i - 1)
               resid(i) = resid(i) - trial(j) * (x(i - j) - fit%mean)
            end do
            do j = 1, min(q_eff, i - 1)
               resid(i) = resid(i) - trial(p_eff + j) * resid(i - j)
            end do
            if (i > max_lag) trial_rss = trial_rss + resid(i) * resid(i)
         end do
         if (trial_rss < best_rss) then
            best_rss = trial_rss
            params = trial
            best_resid = resid
         end if
      end do
   end do
   step = 0.5_dp * step
end do
if (p_eff > 0) fit%coef(1:p_eff) = params(1:p_eff)
if (q_eff > 0 .and. size(fit%coef) >= p_eff + q_eff) fit%coef(p_eff + 1:p_eff + q_eff) = params(p_eff + 1:p_eff + q_eff)
if (include_mean_use) fit%coef(size(fit%coef)) = fit%mean
fit%resid = best_resid
fit%last_x = x(n)
fit%last_resid = best_resid(n)
fit%sigma2 = max(best_rss / real(n_eff, kind=dp), tiny(1.0_dp))
k = size(fit%coef)
fit%aic = real(n_eff, kind=dp) * log(fit%sigma2) + 2.0_dp * real(k, kind=dp)
end function arima_fit

pure function arima_predict(fit, n_ahead) result(pred)
! Compute R-like time-series helper arima_predict.
type(arima_fit_t), intent(in) :: fit ! input value
integer, intent(in) :: n_ahead
real(kind=dp), allocatable :: pred(:)
integer :: i
real(kind=dp) :: phi, theta, prev, prev_resid
allocate(pred(max(0, n_ahead)), source=fit%mean)
if (n_ahead <= 0) return
phi = 0.0_dp
if (allocated(fit%coef) .and. fit%p > 0) phi = fit%coef(1)
theta = 0.0_dp
if (allocated(fit%coef) .and. fit%q > 0 .and. size(fit%coef) >= fit%p + fit%q) theta = fit%coef(fit%p + 1)
prev = fit%last_x
prev_resid = fit%last_resid
do i = 1, n_ahead
   pred(i) = fit%mean + phi * (prev - fit%mean) + theta * prev_resid
   prev = pred(i)
   prev_resid = 0.0_dp
end do
end function arima_predict

pure function arima_predict_result(fit, n_ahead) result(out)
! Compute R-like time-series helper arima_predict_result.
type(arima_fit_t), intent(in) :: fit ! input value
integer, intent(in) :: n_ahead
type(arima_predict_result_t) :: out
out%pred = arima_predict(fit, n_ahead)
allocate(out%se(max(0, n_ahead)), source=sqrt(max(fit%sigma2, 0.0_dp)))
end function arima_predict_result

function ar_fit(x, order_max, aic, method) result(fit)
! Compute R-like time-series helper ar_fit.
! The method argument is accepted for R compatibility but is currently ignored.
real(kind=dp), intent(in) :: x(:) ! input vector
integer, intent(in), optional :: order_max ! maximum AR order to consider
logical, intent(in), optional :: aic ! select order by AIC flag
character(len=*), intent(in), optional :: method ! accepted for compatibility; ignored
type(ar_fit_t) :: fit
integer :: pmax, p, n, i, j, k
real(kind=dp), allocatable :: xt(:), rhs(:), mat(:,:), coef(:)
real(kind=dp) :: mu, sse
if (present(aic)) continue
if (present(method)) continue
n = size(x)
pmax = 1
if (present(order_max)) pmax = max(0, order_max)
if (pmax <= 0 .or. n <= 1) then
   fit%order = 0
   allocate(fit%ar(0), fit%aic(1), source=0.0_dp)
   fit%var_pred = 0.0_dp
   return
end if
p = min(pmax, max(1, min(2, n - 1)))
fit%order = p
allocate(fit%ar(p), fit%aic(pmax + 1), source=0.0_dp)
mu = sum(x) / real(n, kind=dp)
allocate(xt(n), mat(p, p), rhs(p), coef(p))
xt = x - mu
mat = 0.0_dp
rhs = 0.0_dp
do i = p + 1, n
   do j = 1, p
      rhs(j) = rhs(j) + xt(i) * xt(i - j)
      mat(j, :) = mat(j, :) + xt(i - j) * xt(i - [(k, k = 1, p)])
   end do
end do
coef = solve_real(mat + diag(spread(1.0e-10_dp, dim=1, ncopies=p)), rhs)
fit%ar = coef
sse = 0.0_dp
do i = p + 1, n
   sse = sse + (xt(i) - sum(coef * xt(i - [(j, j = 1, p)])))**2
end do
fit%var_pred = sse / real(max(1, n - p), kind=dp)
do i = 1, size(fit%aic)
   fit%aic(i) = real(i - 1, kind=dp)
end do
end function ar_fit

subroutine print_arima_fit(fit)
! Print arima fit values in an R-like format.
type(arima_fit_t), intent(in) :: fit ! input value
write(*,*) "Call:"
write(*,*) "arima(x = x, order = c(", fit%p, ",", fit%d, ",", fit%q, "))"
write(*,*)
write(*,*) "Coefficients:"
if (allocated(fit%coef)) write(*,"(*(g0,1x))") fit%coef
write(*,*)
write(*,*) "sigma^2 estimated as ", fit%sigma2, ":  log likelihood omitted,  aic = ", fit%aic
end subroutine print_arima_fit

function acf_vec(x, lag_max, type, plot) result(fit)
! Compute R-like time-series helper acf_vec.
! Valid type values: "correlation" (default behavior), "covariance".
real(kind=dp), intent(in) :: x(:) ! input vector
integer, intent(in), optional :: lag_max ! maximum lag to return
character(len=*), intent(in), optional :: type ! autocorrelation output type
logical, intent(in), optional :: plot ! request plot warning flag
type(acf_fit_t) :: fit
real(kind=dp), allocatable :: xm(:,:)
allocate(xm(size(x), 1))
xm(:, 1) = x
fit = acf_mat(xm, lag_max=lag_max, type=type, plot=plot)
end function acf_vec

function acf_mat(x, lag_max, type, plot) result(fit)
! Compute R-like time-series helper acf_mat.
! Valid type values: "correlation" (default behavior), "covariance".
real(kind=dp), intent(in) :: x(:,:) ! input matrix
integer, intent(in), optional :: lag_max ! maximum lag to return
character(len=*), intent(in), optional :: type ! autocorrelation output type
logical, intent(in), optional :: plot ! request plot warning flag
type(acf_fit_t) :: fit
integer :: n, lag_n
logical :: do_cov
n = size(x, 1)
if (present(lag_max)) then
   lag_n = lag_max
else
   lag_n = min(n - 1, int(10.0_dp * log10(real(max(n, 2), kind=dp))))
end if
lag_n = max(0, min(lag_n, max(0, n - 1)))
do_cov = .false.
if (present(type)) do_cov = trim(type) == "covariance"
if (present(plot)) then
   if (plot) write(*,*) "Warning: acf plot = TRUE requested; plots are not supported."
end if
fit = acf_mat_impl(x, lag_n, do_cov)
end function acf_mat

pure function acf_mat_impl(x, lag_n, do_cov) result(fit)
! Pure worker for acf_mat after argument normalization and warning handling.
real(kind=dp), intent(in) :: x(:,:)
integer, intent(in) :: lag_n
logical, intent(in) :: do_cov
type(acf_fit_t) :: fit
integer :: n, p, h, i, j, k, cnt
real(kind=dp), allocatable :: mu(:), var0(:)
real(kind=dp) :: s
n = size(x, 1)
p = size(x, 2)
fit%n_used = n
fit%type_code = merge(2, 1, do_cov)
allocate(fit%acf(lag_n + 1, p, p), fit%lag(lag_n + 1), mu(p), var0(p))
fit%acf = ieee_value(0.0_dp, ieee_quiet_nan)
do h = 0, lag_n
   fit%lag(h + 1) = real(h, kind=dp)
end do
do j = 1, p
   s = 0.0_dp
   cnt = 0
   do i = 1, n
      if (ieee_is_finite(x(i, j))) then
         s = s + x(i, j)
         cnt = cnt + 1
      end if
   end do
   if (cnt <= 0) then
      mu(j) = ieee_value(0.0_dp, ieee_quiet_nan)
      var0(j) = ieee_value(0.0_dp, ieee_quiet_nan)
   else
      mu(j) = s / real(cnt, kind=dp)
      s = 0.0_dp
      do i = 1, n
         if (ieee_is_finite(x(i, j))) s = s + (x(i, j) - mu(j))**2
      end do
      var0(j) = s / real(n, kind=dp)
   end if
end do
do h = 0, lag_n
   do j = 1, p
      do k = 1, p
         s = 0.0_dp
         cnt = 0
         do i = 1, n - h
            if (ieee_is_finite(x(i, j)) .and. ieee_is_finite(x(i + h, k))) then
               s = s + (x(i, j) - mu(j)) * (x(i + h, k) - mu(k))
               cnt = cnt + 1
            end if
         end do
         if (cnt <= 0) cycle
         s = s / real(n, kind=dp)
         if (do_cov) then
            fit%acf(h + 1, j, k) = s
         elseif (var0(j) > 0.0_dp .and. var0(k) > 0.0_dp) then
            fit%acf(h + 1, j, k) = s / sqrt(var0(j) * var0(k))
         end if
      end do
   end do
end do
end function acf_mat_impl

function acf_values_vec(x, lag_max, type, plot) result(vals)
! Compute R-like time-series helper acf_values_vec.
! Valid type values: "correlation" (default behavior), "covariance".
real(kind=dp), intent(in) :: x(:) ! input vector
integer, intent(in), optional :: lag_max ! maximum lag to return
character(len=*), intent(in), optional :: type ! autocorrelation output type
logical, intent(in), optional :: plot ! request plot warning flag
real(kind=dp), allocatable :: vals(:)
type(acf_fit_t) :: fit
fit = acf_vec(x, lag_max=lag_max, type=type, plot=plot)
vals = reshape(fit%acf, [size(fit%acf)])
end function acf_values_vec

function acf_values_mat(x, lag_max, type, plot) result(vals)
! Compute R-like time-series helper acf_values_mat.
! Valid type values: "correlation" (default behavior), "covariance".
real(kind=dp), intent(in) :: x(:,:) ! input matrix
integer, intent(in), optional :: lag_max ! maximum lag to return
character(len=*), intent(in), optional :: type ! autocorrelation output type
logical, intent(in), optional :: plot ! request plot warning flag
real(kind=dp), allocatable :: vals(:)
type(acf_fit_t) :: fit
fit = acf_mat(x, lag_max=lag_max, type=type, plot=plot)
vals = reshape(fit%acf, [size(fit%acf)])
end function acf_values_mat

pure function ARMAacf_vecma(ar, ma, lag_max) result(vals)
! Compute R-like theoretical autocorrelations for an ARMA model.
real(kind=dp), intent(in), optional :: ar(:) ! autoregressive coefficients
real(kind=dp), intent(in), optional :: ma(:) ! moving-average coefficients
integer, intent(in), optional :: lag_max ! maximum lag to return
real(kind=dp), allocatable :: vals(:)
real(kind=dp), allocatable :: psi(:)
integer :: lag_n, h, p, q, npsi, j, i
real(kind=dp) :: denom, numer, tail_scale
p = 0
q = 0
if (present(ar)) p = size(ar)
if (present(ma)) q = size(ma)
lag_n = max(p, q + 1)
if (present(lag_max)) lag_n = max(0, lag_max)
allocate(vals(lag_n + 1), source=0.0_dp)
vals(1) = 1.0_dp
if (p == 0 .and. q == 0) return
npsi = max(lag_n + q + 1, lag_n + 10000)
allocate(psi(0:npsi), source=0.0_dp)
psi(0) = 1.0_dp
do j = 1, npsi
   if (present(ma)) then
      if (j <= q) psi(j) = psi(j) + ma(j)
   end if
   if (present(ar)) then
      do i = 1, min(p, j)
         psi(j) = psi(j) + ar(i) * psi(j - i)
      end do
   end if
end do
denom = sum(psi * psi)
if (denom <= tiny(1.0_dp)) return
do h = 1, lag_n
   if (h <= npsi) then
      numer = sum(psi(0:npsi-h) * psi(h:npsi))
      vals(h + 1) = numer / denom
   end if
end do
if (p > 0) then
   tail_scale = maxval(abs(psi(max(0, npsi - 99):npsi)))
   if (tail_scale > 1.0e-7_dp) then
      ! Near-nonstationary models can require more terms than the fixed truncation.
      ! Preserve finite output but avoid pretending to exactness beyond available precision.
      vals = max(-1.0_dp, min(1.0_dp, vals))
   end if
end if
end function ARMAacf_vecma

pure function ARMAacf_scalarma(ar, ma, lag_max) result(vals)
! Scalar-MA compatibility wrapper for calls such as ARMAacf(ar=..., ma=0.5).
real(kind=dp), intent(in), optional :: ar(:) ! autoregressive coefficients
real(kind=dp), intent(in) :: ma ! scalar moving-average coefficient
integer, intent(in), optional :: lag_max ! maximum lag to return
real(kind=dp), allocatable :: vals(:)
real(kind=dp) :: ma_vec(1)
ma_vec = [ma]
vals = ARMAacf_vecma(ar=ar, ma=ma_vec, lag_max=lag_max)
end function ARMAacf_scalarma

function ccf_vec(x, y, lag_max, type, plot) result(fit)
! Runtime helper for R-compatible ccf vec.
! Valid type values: "correlation" (default behavior), "covariance".
real(kind=dp), intent(in) :: x(:) ! input vector
real(kind=dp), intent(in) :: y(:) ! response values
integer, intent(in), optional :: lag_max ! maximum lag to return
character(len=*), intent(in), optional :: type ! cross-correlation output type
logical, intent(in), optional :: plot ! request plot warning flag
type(acf_fit_t) :: fit
integer :: n, lag_n
logical :: do_cov
n = min(size(x), size(y))
if (present(lag_max)) then
   lag_n = lag_max
else
   lag_n = min(n - 1, int(10.0_dp * log10(real(max(n, 2), kind=dp))))
end if
lag_n = max(0, min(lag_n, max(0, n - 1)))
do_cov = .false.
if (present(type)) do_cov = trim(type) == "covariance"
if (present(plot)) then
   if (plot) write(*,*) "Warning: ccf plot = TRUE requested; plots are not supported."
end if
fit = ccf_vec_impl(x, y, n, lag_n, do_cov)
end function ccf_vec

pure function ccf_vec_impl(x, y, n, lag_n, do_cov) result(fit)
! Pure worker for ccf_vec after argument normalization and warning handling.
real(kind=dp), intent(in) :: x(:)
real(kind=dp), intent(in) :: y(:)
integer, intent(in) :: n, lag_n
logical, intent(in) :: do_cov
type(acf_fit_t) :: fit
integer :: h, ii, idx, cnt
real(kind=dp) :: mux, muy, vx, vy, s
fit%n_used = n
fit%type_code = merge(2, 1, do_cov)
allocate(fit%acf(2 * lag_n + 1, 1, 1), fit%lag(2 * lag_n + 1))
fit%acf = ieee_value(0.0_dp, ieee_quiet_nan)
mux = sum(x(1:n)) / real(max(n, 1), kind=dp)
muy = sum(y(1:n)) / real(max(n, 1), kind=dp)
vx = sum((x(1:n) - mux)**2) / real(max(n, 1), kind=dp)
vy = sum((y(1:n) - muy)**2) / real(max(n, 1), kind=dp)
idx = 0
do h = -lag_n, lag_n
   idx = idx + 1
   fit%lag(idx) = real(h, kind=dp)
   s = 0.0_dp
   cnt = 0
   do ii = 1, n
      if (ii + h < 1 .or. ii + h > n) cycle
      s = s + (x(ii) - mux) * (y(ii + h) - muy)
      cnt = cnt + 1
   end do
   if (cnt <= 0) cycle
   s = s / real(n, kind=dp)
   if (do_cov) then
      fit%acf(idx, 1, 1) = s
   elseif (vx > 0.0_dp .and. vy > 0.0_dp) then
      fit%acf(idx, 1, 1) = s / sqrt(vx * vy)
   end if
end do
end function ccf_vec_impl

subroutine print_acf(fit, digits, series_name)
! Print acf values in an R-like format.
type(acf_fit_t), intent(in) :: fit ! input value
integer, intent(in), optional :: digits
character(len=*), intent(in), optional :: series_name
integer :: d, w, nlag, i, i2, j, per_line
character(len=32) :: ifmt, rfmt
d = 4
if (present(digits)) d = max(0, digits)
w = max(7, d + 4)
per_line = max(1, 80 / (w + 1))
nlag = size(fit%acf, 1)
write(ifmt, '("(i", i0, ",1x)")') w
write(rfmt, '("(f", i0, ".", i0, ",1x)")') w, d
write(*,*)
if (present(series_name)) then
   write(*,'(a)') "Autocorrelations of series '" // trim(series_name) // "', by lag"
else
   write(*,'(a)') "Autocorrelations of series, by lag"
end if
write(*,*)
do i = 1, nlag, per_line
   i2 = min(nlag, i + per_line - 1)
   do j = i, i2
      write(*, ifmt, advance='no') nint(fit%lag(j))
   end do
   write(*,*)
   do j = i, i2
      write(*, rfmt, advance='no') fit%acf(j, 1, 1)
   end do
   write(*,*)
end do
end subroutine print_acf

pure function besselJ_core(x, nu) result(out)
! Evaluate the besselJ core Bessel helper.
real(kind=dp), intent(in) :: x, nu
real(kind=dp) :: out, term, ax
integer :: k, n
ax = abs(x)
if (abs(nu - real(nint(nu), kind=dp)) <= 100.0_dp * epsilon(1.0_dp)) then
   n = abs(nint(nu))
   select case (n)
   case (0)
      out = bessel_j0(x)
   case (1)
      out = bessel_j1(x)
   case default
      out = bessel_jn(n, x)
   end select
   if (nu < 0.0_dp .and. mod(n, 2) /= 0) out = -out
   return
end if
if (ax == 0.0_dp) then
   out = merge(1.0_dp, 0.0_dp, abs(nu) <= 100.0_dp * epsilon(1.0_dp))
   return
end if
term = exp(nu * log(0.5_dp * ax) - log_gamma(nu + 1.0_dp))
out = term
do k = 1, 200
   term = -term * (0.25_dp * ax * ax) / (real(k, kind=dp) * (real(k, kind=dp) + nu))
   out = out + term
   if (abs(term) <= 1.0e-15_dp * max(1.0_dp, abs(out))) exit
end do
if (x < 0.0_dp .and. abs(nu - real(nint(nu), kind=dp)) <= 100.0_dp * epsilon(1.0_dp) .and. mod(nint(nu), &
& 2) /= 0) out = -out
end function besselJ_core

pure function besselI_core(x, nu, scaled) result(out)
! Evaluate the besselI core Bessel helper.
real(kind=dp), intent(in) :: x ! input values
real(kind=dp), intent(in) :: nu ! input value
logical, intent(in) :: scaled
real(kind=dp) :: out, term, ax, nuf
integer :: k
ax = abs(x)
nuf = nu
if (abs(nu - real(nint(nu), kind=dp)) <= 100.0_dp * epsilon(1.0_dp)) nuf = abs(nu)
if (ax == 0.0_dp) then
   out = merge(1.0_dp, 0.0_dp, abs(nuf) <= 100.0_dp * epsilon(1.0_dp))
   return
end if
term = exp(nuf * log(0.5_dp * ax) - log_gamma(nuf + 1.0_dp))
out = term
do k = 1, 300
   term = term * (0.25_dp * ax * ax) / (real(k, kind=dp) * (real(k, kind=dp) + nuf))
   out = out + term
   if (abs(term) <= 1.0e-15_dp * max(1.0_dp, abs(out))) exit
end do
if (scaled) out = out * exp(-ax)
end function besselI_core

pure function besselY_core(x, nu) result(out)
! Evaluate the besselY core Bessel helper.
real(kind=dp), intent(in) :: x, nu
real(kind=dp) :: out, eps, nuf
integer :: n
if (x <= 0.0_dp) then
   out = -huge(1.0_dp)
   return
end if
if (abs(nu - real(nint(nu), kind=dp)) <= 100.0_dp * epsilon(1.0_dp)) then
   n = abs(nint(nu))
   select case (n)
   case (0)
      out = bessel_y0(x)
   case (1)
      out = bessel_y1(x)
   case default
      out = bessel_yn(n, x)
   end select
   if (nu < 0.0_dp .and. mod(n, 2) /= 0) out = -out
   return
end if
eps = 1.0e-6_dp
nuf = nu
if (abs(sin(acos(-1.0_dp) * nuf)) < 1.0e-8_dp) nuf = nuf + eps
out = (besselJ_core(x, nuf) * cos(acos(-1.0_dp) * nuf) - besselJ_core(x, -nuf)) / sin(acos(-1.0_dp) * nuf)
end function besselY_core

pure function besselK0_core(x, scaled) result(out)
! Approximate K_0(x), optionally returning exp(x) * K_0(x).
real(kind=dp), intent(in) :: x
logical, intent(in) :: scaled
real(kind=dp) :: out, y, poly
if (x <= 0.0_dp) then
   out = huge(1.0_dp)
   return
end if
if (x <= 2.0_dp) then
   y = 0.25_dp * x * x
   poly = -0.57721566_dp + y * (0.42278420_dp + y * (0.23069756_dp + y * (0.03488590_dp + &
      y * (0.00262698_dp + y * (0.00010750_dp + y * 0.00000740_dp)))))
   out = -log(0.5_dp * x) * besselI_core(x, 0.0_dp, .false.) + poly
   if (scaled) out = out * exp(x)
else
   y = 2.0_dp / x
   poly = 1.25331414_dp + y * (-0.07832358_dp + y * (0.02189568_dp + y * (-0.01062446_dp + &
      y * (0.00587872_dp + y * (-0.00251540_dp + y * 0.00053208_dp)))))
   out = poly / sqrt(x)
   if (.not. scaled) out = out * exp(-x)
end if
end function besselK0_core

pure function besselK1_core(x, scaled) result(out)
! Approximate K_1(x), optionally returning exp(x) * K_1(x).
real(kind=dp), intent(in) :: x
logical, intent(in) :: scaled
real(kind=dp) :: out, y, poly
if (x <= 0.0_dp) then
   out = huge(1.0_dp)
   return
end if
if (x <= 2.0_dp) then
   y = 0.25_dp * x * x
   poly = 1.0_dp + y * (0.15443144_dp + y * (-0.67278579_dp + y * (-0.18156897_dp + &
      y * (-0.01919402_dp + y * (-0.00110404_dp + y * (-0.00004686_dp))))))
   out = log(0.5_dp * x) * besselI_core(x, 1.0_dp, .false.) + poly / x
   if (scaled) out = out * exp(x)
else
   y = 2.0_dp / x
   poly = 1.25331414_dp + y * (0.23498619_dp + y * (-0.03655620_dp + y * (0.01504268_dp + &
      y * (-0.00780353_dp + y * (0.00325614_dp + y * (-0.00068245_dp))))))
   out = poly / sqrt(x)
   if (.not. scaled) out = out * exp(-x)
end if
end function besselK1_core

pure function besselK_core(x, nu, scaled) result(out)
! Evaluate the besselK core Bessel helper.
real(kind=dp), intent(in) :: x ! input values
real(kind=dp), intent(in) :: nu ! input value
logical, intent(in) :: scaled
real(kind=dp) :: out, nuf, s, km1, kcur, kp1
integer :: i, n
if (x <= 0.0_dp) then
   out = huge(1.0_dp)
   return
end if
nuf = abs(nu)
if (abs(nuf - real(nint(nuf), kind=dp)) <= 100.0_dp * epsilon(1.0_dp)) then
   n = nint(nuf)
   select case (n)
   case (0)
      out = besselK0_core(x, scaled)
      return
   case (1)
      out = besselK1_core(x, scaled)
      return
   case default
      km1 = besselK0_core(x, scaled)
      kcur = besselK1_core(x, scaled)
      do i = 1, n - 1
         kp1 = km1 + (2.0_dp * real(i, kind=dp) / x) * kcur
         km1 = kcur
         kcur = kp1
      end do
      out = kcur
      return
   end select
end if
if (abs(nuf - 0.5_dp) <= 100.0_dp * epsilon(1.0_dp)) then
   out = sqrt(0.5_dp * acos(-1.0_dp) / x)
   if (.not. scaled) out = out * exp(-x)
   return
end if
if (abs(sin(acos(-1.0_dp) * nuf)) > 1.0e-7_dp .and. x <= 20.0_dp) then
   out = 0.5_dp * acos(-1.0_dp) * (besselI_core(x, -nuf, .false.) - besselI_core(x, nuf, &
   & .false.)) / sin(acos(-1.0_dp) * nuf)
else
   s = sqrt(0.5_dp * acos(-1.0_dp) / x)
   out = s * (1.0_dp + (4.0_dp * nuf * nuf - 1.0_dp) / (8.0_dp * max(x, 1.0e-6_dp)))
   if (.not. scaled) out = out * exp(-x)
end if
if (scaled .and. abs(sin(acos(-1.0_dp) * nuf)) > 1.0e-7_dp .and. x <= 20.0_dp) out = out * exp(x)
end function besselK_core

pure function besselJ_scalar_i(x, nu) result(out)
! Evaluate the besselJ scalar i Bessel helper.
real(kind=dp), intent(in) :: x ! input values
integer, intent(in) :: nu
real(kind=dp) :: out
out = besselJ_core(x, real(nu, kind=dp))
end function besselJ_scalar_i

pure function besselJ_scalar_r(x, nu) result(out)
! Evaluate the besselJ scalar r Bessel helper.
real(kind=dp), intent(in) :: x, nu
real(kind=dp) :: out
out = besselJ_core(x, nu)
end function besselJ_scalar_r

pure function besselJ_vec_i(x, nu) result(out)
! Evaluate the besselJ vec i Bessel helper.
real(kind=dp), intent(in) :: x(:) ! input vector
integer, intent(in) :: nu
real(kind=dp), allocatable :: out(:)
integer :: i
allocate(out(size(x)))
do i = 1, size(x)
   out(i) = besselJ_core(x(i), real(nu, kind=dp))
end do
end function besselJ_vec_i

pure function besselJ_vec_r(x, nu) result(out)
! Evaluate the besselJ vec r Bessel helper.
real(kind=dp), intent(in) :: x(:) ! input vector
real(kind=dp), intent(in) :: nu
real(kind=dp), allocatable :: out(:)
integer :: i
allocate(out(size(x)))
do i = 1, size(x)
   out(i) = besselJ_core(x(i), nu)
end do
end function besselJ_vec_r

pure function besselY_scalar_i(x, nu) result(out)
! Evaluate the besselY scalar i Bessel helper.
real(kind=dp), intent(in) :: x ! input values
integer, intent(in) :: nu
real(kind=dp) :: out
out = besselY_core(x, real(nu, kind=dp))
end function besselY_scalar_i

pure function besselY_scalar_r(x, nu) result(out)
! Evaluate the besselY scalar r Bessel helper.
real(kind=dp), intent(in) :: x, nu
real(kind=dp) :: out
out = besselY_core(x, nu)
end function besselY_scalar_r

pure function besselY_vec_i(x, nu) result(out)
! Evaluate the besselY vec i Bessel helper.
real(kind=dp), intent(in) :: x(:) ! input vector
integer, intent(in) :: nu
real(kind=dp), allocatable :: out(:)
integer :: i
allocate(out(size(x)))
do i = 1, size(x)
   out(i) = besselY_core(x(i), real(nu, kind=dp))
end do
end function besselY_vec_i

pure function besselY_vec_r(x, nu) result(out)
! Evaluate the besselY vec r Bessel helper.
real(kind=dp), intent(in) :: x(:) ! input vector
real(kind=dp), intent(in) :: nu
real(kind=dp), allocatable :: out(:)
integer :: i
allocate(out(size(x)))
do i = 1, size(x)
   out(i) = besselY_core(x(i), nu)
end do
end function besselY_vec_r

pure function besselI_scalar_i(x, nu, expon_scaled) result(out)
! Evaluate the besselI scalar i Bessel helper.
real(kind=dp), intent(in) :: x ! input values
integer, intent(in) :: nu ! integer argument
logical, intent(in), optional :: expon_scaled
real(kind=dp) :: out
logical :: scaled
scaled = .false.
if (present(expon_scaled)) scaled = expon_scaled
out = besselI_core(x, real(nu, kind=dp), scaled)
end function besselI_scalar_i

pure function besselI_scalar_r(x, nu, expon_scaled) result(out)
! Evaluate the besselI scalar r Bessel helper.
real(kind=dp), intent(in) :: x ! input values
real(kind=dp), intent(in) :: nu ! input value
logical, intent(in), optional :: expon_scaled
real(kind=dp) :: out
logical :: scaled
scaled = .false.
if (present(expon_scaled)) scaled = expon_scaled
out = besselI_core(x, nu, scaled)
end function besselI_scalar_r

pure function besselI_vec_i(x, nu, expon_scaled) result(out)
! Evaluate the besselI vec i Bessel helper.
real(kind=dp), intent(in) :: x(:) ! input vector
integer, intent(in) :: nu ! integer argument
logical, intent(in), optional :: expon_scaled
real(kind=dp), allocatable :: out(:)
integer :: i
logical :: scaled
allocate(out(size(x)))
scaled = .false.
if (present(expon_scaled)) scaled = expon_scaled
do i = 1, size(x)
   out(i) = besselI_core(x(i), real(nu, kind=dp), scaled)
end do
end function besselI_vec_i

pure function besselI_vec_r(x, nu, expon_scaled) result(out)
! Evaluate the besselI vec r Bessel helper.
real(kind=dp), intent(in) :: x(:) ! input vector
real(kind=dp), intent(in) :: nu ! input value
logical, intent(in), optional :: expon_scaled
real(kind=dp), allocatable :: out(:)
integer :: i
logical :: scaled
allocate(out(size(x)))
scaled = .false.
if (present(expon_scaled)) scaled = expon_scaled
do i = 1, size(x)
   out(i) = besselI_core(x(i), nu, scaled)
end do
end function besselI_vec_r

pure function besselK_scalar_i(x, nu, expon_scaled) result(out)
! Evaluate the besselK scalar i Bessel helper.
real(kind=dp), intent(in) :: x ! input values
integer, intent(in) :: nu ! integer argument
logical, intent(in), optional :: expon_scaled
real(kind=dp) :: out
logical :: scaled
scaled = .false.
if (present(expon_scaled)) scaled = expon_scaled
out = besselK_core(x, real(nu, kind=dp), scaled)
end function besselK_scalar_i

pure function besselK_scalar_r(x, nu, expon_scaled) result(out)
! Evaluate the besselK scalar r Bessel helper.
real(kind=dp), intent(in) :: x ! input values
real(kind=dp), intent(in) :: nu ! input value
logical, intent(in), optional :: expon_scaled
real(kind=dp) :: out
logical :: scaled
scaled = .false.
if (present(expon_scaled)) scaled = expon_scaled
out = besselK_core(x, nu, scaled)
end function besselK_scalar_r

pure function besselK_vec_i(x, nu, expon_scaled) result(out)
! Evaluate the besselK vec i Bessel helper.
real(kind=dp), intent(in) :: x(:) ! input vector
integer, intent(in) :: nu ! integer argument
logical, intent(in), optional :: expon_scaled
real(kind=dp), allocatable :: out(:)
integer :: i
logical :: scaled
allocate(out(size(x)))
scaled = .false.
if (present(expon_scaled)) scaled = expon_scaled
do i = 1, size(x)
   out(i) = besselK_core(x(i), real(nu, kind=dp), scaled)
end do
end function besselK_vec_i

pure function besselK_vec_r(x, nu, expon_scaled) result(out)
! Evaluate the besselK vec r Bessel helper.
real(kind=dp), intent(in) :: x(:) ! input vector
real(kind=dp), intent(in) :: nu ! input value
logical, intent(in), optional :: expon_scaled
real(kind=dp), allocatable :: out(:)
integer :: i
logical :: scaled
allocate(out(size(x)))
scaled = .false.
if (present(expon_scaled)) scaled = expon_scaled
do i = 1, size(x)
   out(i) = besselK_core(x(i), nu, scaled)
end do
end function besselK_vec_r

pure function solve_real_vec(a, b) result(x)
! Return the solution of a square linear system a %*% x = b.
real(kind=dp), intent(in) :: a(:,:) ! square coefficient matrix
real(kind=dp), intent(in) :: b(:) ! right-hand side vector
real(kind=dp), allocatable :: x(:)
real(kind=dp), allocatable :: aa(:,:), bb(:)
integer :: i, j, k, n, p
real(kind=dp) :: fac, piv, s, t
n = size(b)
allocate(x(n))
x = 0.0_dp
if (size(a, 1) /= n .or. size(a, 2) /= n) return
aa = a
bb = b
do k = 1, n
   p = k
   piv = abs(aa(k, k))
   do i = k + 1, n
      if (abs(aa(i, k)) > piv) then
         p = i
         piv = abs(aa(i, k))
      end if
   end do
   if (piv <= tiny(1.0_dp)) return
   if (p /= k) then
      do j = k, n
         t = aa(k, j)
         aa(k, j) = aa(p, j)
         aa(p, j) = t
      end do
      t = bb(k)
      bb(k) = bb(p)
      bb(p) = t
   end if
   do i = k + 1, n
      fac = aa(i, k) / aa(k, k)
      aa(i, k:n) = aa(i, k:n) - fac * aa(k, k:n)
      bb(i) = bb(i) - fac * bb(k)
   end do
end do
do i = n, 1, -1
   s = bb(i)
   if (i < n) s = s - sum(aa(i, i+1:n) * x(i+1:n))
   x(i) = s / aa(i, i)
end do
end function solve_real_vec

pure function solve_real_vec_i_r(a, b) result(x)
! Solve a real linear system with mixed numeric input kinds.
integer, intent(in) :: a(:,:) ! coefficient matrix
real(kind=dp), intent(in) :: b(:) ! right-hand side vector
real(kind=dp), allocatable :: x(:)
x = solve_real_vec(real(a, kind=dp), b)
end function solve_real_vec_i_r

pure function solve_real_vec_i_i(a, b) result(x)
! Solve a real linear system with mixed numeric input kinds.
integer, intent(in) :: a(:,:) ! coefficient matrix
integer, intent(in) :: b(:) ! right-hand side vector
real(kind=dp), allocatable :: x(:)
x = solve_real_vec(real(a, kind=dp), real(b, kind=dp))
end function solve_real_vec_i_i

pure function solve_real_mat(a, b) result(x)
! Return the solution of a square linear system a %*% x = b for matrix b.
real(kind=dp), intent(in) :: a(:,:) ! square coefficient matrix
real(kind=dp), intent(in) :: b(:,:) ! right-hand side matrix
real(kind=dp), allocatable :: x(:,:)
integer :: j, n, m
n = size(a, 1)
m = size(b, 2)
allocate(x(n, m))
x = 0.0_dp
if (size(a, 2) /= n .or. size(b, 1) /= n) return
do j = 1, m
   x(:, j) = solve_real_vec(a, b(:, j))
end do
end function solve_real_mat

pure function solve_real_mat_i_r(a, b) result(x)
! Solve a real linear system with mixed numeric input kinds.
integer, intent(in) :: a(:,:) ! coefficient matrix
real(kind=dp), intent(in) :: b(:,:) ! right-hand side matrix
real(kind=dp), allocatable :: x(:,:)
x = solve_real_mat(real(a, kind=dp), b)
end function solve_real_mat_i_r

pure function solve_real_mat_r_i(a, b) result(x)
! Solve a real linear system with mixed numeric input kinds.
real(kind=dp), intent(in) :: a(:,:) ! coefficient matrix
integer, intent(in) :: b(:,:) ! right-hand side matrix
real(kind=dp), allocatable :: x(:,:)
x = solve_real_mat(a, real(b, kind=dp))
end function solve_real_mat_r_i

pure function solve_real_mat_i_i(a, b) result(x)
! Solve a real linear system with mixed numeric input kinds.
integer, intent(in) :: a(:,:) ! coefficient matrix
integer, intent(in) :: b(:,:) ! right-hand side matrix
real(kind=dp), allocatable :: x(:,:)
x = solve_real_mat(real(a, kind=dp), real(b, kind=dp))
end function solve_real_mat_i_i

pure function mahalanobis(x, center, cov) result(out)
! Runtime helper for R-compatible mahalanobis.
real(kind=dp), intent(in) :: x(:,:) ! observations by rows
real(kind=dp), intent(in) :: center(:) ! center vector to subtract
real(kind=dp), intent(in) :: cov(:,:) ! covariance matrix
real(kind=dp), allocatable :: out(:), z(:,:), sol(:,:)
integer :: i, n, p
n = size(x, 1)
p = size(x, 2)
allocate(out(n))
out = 0.0_dp
if (size(center) /= p .or. size(cov, 1) /= p .or. size(cov, 2) /= p) return
allocate(z(n, p))
do i = 1, n
   z(i, :) = x(i, :) - center
end do
sol = solve_real_mat(cov, transpose(z))
do i = 1, n
   out(i) = sum(z(i, :) * sol(:, i))
end do
end function mahalanobis

pure function isSymmetric_real(x, tol) result(out)
! Test the R-like predicate isSymmetric_real.
! If tol is absent, use 100*epsilon(1.0_dp).
real(kind=dp), intent(in) :: x(:,:) ! matrix to test for symmetry
real(kind=dp), intent(in), optional :: tol ! absolute symmetry tolerance
logical :: out
real(kind=dp) :: eps
if (size(x, 1) /= size(x, 2)) then
   out = .false.
   return
end if
eps = 100.0_dp * epsilon(1.0_dp)
if (present(tol)) eps = tol
out = all(abs(x - transpose(x)) <= eps)
end function isSymmetric_real

pure function isSymmetric_int(x, tol) result(out)
! Test the R-like predicate isSymmetric_int.
! If tol is absent, use 100*epsilon(1.0_dp).
integer, intent(in) :: x(:,:) ! matrix to test for symmetry
real(kind=dp), intent(in), optional :: tol ! absolute symmetry tolerance
logical :: out
if (present(tol)) then
   out = isSymmetric_real(real(x, kind=dp), tol)
else
   out = isSymmetric_real(real(x, kind=dp))
end if
end function isSymmetric_int

pure function solve_real_vec_r_c(a, b) result(x)
! Solve a real linear system with mixed numeric input kinds.
real(kind=dp), intent(in) :: a(:,:) ! coefficient matrix
complex(kind=dp), intent(in) :: b(:) ! right-hand side vector
complex(kind=dp), allocatable :: x(:)
x = solve_complex_vec(cmplx(a, 0.0_dp, kind=dp), b)
end function solve_real_vec_r_c

pure function solve_real_vec_i_c(a, b) result(x)
! Solve a real linear system with mixed numeric input kinds.
integer, intent(in) :: a(:,:) ! coefficient matrix
complex(kind=dp), intent(in) :: b(:) ! right-hand side vector
complex(kind=dp), allocatable :: x(:)
x = solve_complex_vec(cmplx(real(a, kind=dp), 0.0_dp, kind=dp), b)
end function solve_real_vec_i_c

pure function solve_complex_vec(a, b) result(x)
! Return the solution of a square complex linear system a %*% x = b.
complex(kind=dp), intent(in) :: a(:,:) ! square coefficient matrix
complex(kind=dp), intent(in) :: b(:) ! right-hand side vector
complex(kind=dp), allocatable :: x(:)
complex(kind=dp), allocatable :: aa(:,:), bb(:)
integer :: i, j, k, n, p
real(kind=dp) :: piv
complex(kind=dp) :: fac, s, t
n = size(b)
allocate(x(n))
x = cmplx(0.0_dp, 0.0_dp, kind=dp)
if (size(a, 1) /= n .or. size(a, 2) /= n) return
aa = a
bb = b
do k = 1, n
   p = k
   piv = abs(aa(k, k))
   do i = k + 1, n
      if (abs(aa(i, k)) > piv) then
         p = i
         piv = abs(aa(i, k))
      end if
   end do
   if (piv <= tiny(1.0_dp)) return
   if (p /= k) then
      do j = k, n
         t = aa(k, j)
         aa(k, j) = aa(p, j)
         aa(p, j) = t
      end do
      t = bb(k)
      bb(k) = bb(p)
      bb(p) = t
   end if
   do i = k + 1, n
      fac = aa(i, k) / aa(k, k)
      aa(i, k:n) = aa(i, k:n) - fac * aa(k, k:n)
      bb(i) = bb(i) - fac * bb(k)
   end do
end do
do i = n, 1, -1
   s = bb(i)
   if (i < n) s = s - sum(aa(i, i+1:n) * x(i+1:n))
   x(i) = s / aa(i, i)
end do
end function solve_complex_vec

pure function solve_complex_mat(a, b) result(x)
! Return the solution of a square complex linear system a %*% x = b for matrix b.
complex(kind=dp), intent(in) :: a(:,:) ! square coefficient matrix
complex(kind=dp), intent(in) :: b(:,:) ! right-hand side matrix
complex(kind=dp), allocatable :: x(:,:)
integer :: j, n, m
n = size(a, 1)
m = size(b, 2)
allocate(x(n, m))
x = cmplx(0.0_dp, 0.0_dp, kind=dp)
if (size(a, 2) /= n .or. size(b, 1) /= n) return
do j = 1, m
   x(:, j) = solve_complex_vec(a, b(:, j))
end do
end function solve_complex_mat

pure function cumsum_real(x) result(out)
! Return cumulative sums of a real vector.
real(kind=dp), intent(in) :: x(:) ! input vector
real(kind=dp), allocatable :: out(:)
integer :: i
allocate(out(size(x)))
if (size(x) <= 0) return
out(1) = x(1)
if (out(1) == 0.0_dp) out(1) = 0.0_dp
do i = 2, size(x)
   out(i) = out(i - 1) + x(i)
   if (out(i) == 0.0_dp) out(i) = 0.0_dp
end do
end function cumsum_real

pure function cumsum_int(x) result(out)
! Return cumulative sums of an integer vector.
integer, intent(in) :: x(:) ! input vector
integer, allocatable :: out(:)
integer :: i
allocate(out(size(x)))
if (size(x) <= 0) return
out(1) = x(1)
do i = 2, size(x)
   if (out(i - 1) == -huge(0) .or. x(i) == -huge(0)) then
      out(i) = -huge(0)
   else
      out(i) = out(i - 1) + x(i)
   end if
end do
end function cumsum_int

pure function apply_col_cumsum(x) result(out)
! Return apply(x, 2, cumsum) for a real matrix.
real(kind=dp), intent(in) :: x(:,:) ! input matrix
real(kind=dp), allocatable :: out(:,:)
integer :: j
allocate(out(size(x, 1), size(x, 2)))
do j = 1, size(x, 2)
   out(:, j) = cumsum_real(x(:, j))
end do
end function apply_col_cumsum

pure function apply_col_sd(x) result(out)
! Return apply(x, 2, sd) for a real matrix.
real(kind=dp), intent(in) :: x(:,:) ! input matrix
real(kind=dp), allocatable :: out(:)
integer :: j
allocate(out(size(x, 2)))
do j = 1, size(x, 2)
   out(j) = sd(x(:, j))
end do
end function apply_col_sd

pure function apply_row_sd(x) result(out)
! Return apply(x, 1, sd) for a real matrix.
real(kind=dp), intent(in) :: x(:,:) ! input matrix
real(kind=dp), allocatable :: out(:)
integer :: i
allocate(out(size(x, 1)))
do i = 1, size(x, 1)
   out(i) = sd(x(i, :))
end do
end function apply_row_sd

pure function findInterval(x, vec) result(out)
! Return R-style interval counts for each x against sorted breakpoints vec.
real(kind=dp), intent(in) :: x(:) ! query values
real(kind=dp), intent(in) :: vec(:) ! sorted breakpoints
integer, allocatable :: out(:)
integer :: i, j
allocate(out(size(x)))
do i = 1, size(x)
   if (x(i) /= x(i)) then
      out(i) = -huge(0)
      cycle
   end if
   out(i) = 0
   do j = 1, size(vec)
      if (x(i) >= vec(j)) out(i) = j
   end do
end do
end function findInterval

pure function cut(x, breaks, include_lowest, labels) result(out)
! Return integer bin numbers for cut(x, breaks, labels = FALSE).
real(kind=dp), intent(in) :: x(:) ! values to classify
real(kind=dp), intent(in) :: breaks(:) ! ordered bin boundaries
logical, intent(in), optional :: include_lowest ! include first lower boundary
logical, intent(in), optional :: labels ! accepted for R call compatibility
integer, allocatable :: out(:)
logical :: inc_low
integer :: i, j, nb
inc_low = .false.
if (present(include_lowest)) inc_low = include_lowest
nb = max(0, size(breaks) - 1)
allocate(out(size(x)))
out = 0
do i = 1, size(x)
   do j = 1, nb
      if (j == 1 .and. inc_low) then
         if (x(i) >= breaks(j) .and. x(i) <= breaks(j + 1)) then
            out(i) = j
            exit
         end if
      else
         if (x(i) > breaks(j) .and. x(i) <= breaks(j + 1)) then
            out(i) = j
            exit
         end if
      end if
   end do
end do
if (present(labels)) continue
end function cut

pure function cut_n(x, breaks, include_lowest, labels) result(out)
! Return integer bin numbers for cut(x, breaks_count, labels = FALSE).
real(kind=dp), intent(in) :: x(:)
integer, intent(in) :: breaks
logical, intent(in), optional :: include_lowest
logical, intent(in), optional :: labels
integer, allocatable :: out(:)
real(kind=dp), allocatable :: edges(:)
real(kind=dp) :: dx, lo, hi
integer :: i
if (size(x) <= 0) then
   allocate(out(0))
   return
end if
if (breaks < 1) then
   allocate(out(size(x)), source=0)
   return
end if
lo = minval(x)
hi = maxval(x)
dx = hi - lo
if (dx == 0.0_dp) dx = max(abs(lo), 1.0_dp)
allocate(edges(breaks + 1))
do i = 1, breaks + 1
   edges(i) = lo + real(i - 1, kind=dp) * (hi - lo) / real(breaks, kind=dp)
end do
edges(1) = edges(1) - dx / 1000.0_dp
edges(breaks + 1) = edges(breaks + 1) + dx / 1000.0_dp
out = cut(x, edges, include_lowest=include_lowest, labels=labels)
end function cut_n

pure function outer(x, y) result(out)
! Return the default R outer(x, y) product matrix.
real(kind=dp), intent(in) :: x(:) ! row-factor vector
real(kind=dp), intent(in) :: y(:) ! column-factor vector
real(kind=dp), allocatable :: out(:,:)
integer :: i, j
allocate(out(size(x), size(y)))
do i = 1, size(x)
   do j = 1, size(y)
      out(i, j) = x(i) * y(j)
   end do
end do
end function outer

pure function outer_plus(x, y) result(out)
! R outer(x, y, "+"): all pairwise sums.
real(kind=dp), intent(in) :: x(:), y(:)
real(kind=dp), allocatable :: out(:,:)
integer :: i, j
allocate(out(size(x), size(y)))
do i = 1, size(x)
   do j = 1, size(y)
      out(i, j) = x(i) + y(j)
   end do
end do
end function outer_plus

pure function outer_minus(x, y) result(out)
! R outer(x, y, "-"): all pairwise differences.
real(kind=dp), intent(in) :: x(:), y(:)
real(kind=dp), allocatable :: out(:,:)
integer :: i, j
allocate(out(size(x), size(y)))
do i = 1, size(x)
   do j = 1, size(y)
      out(i, j) = x(i) - y(j)
   end do
end do
end function outer_minus

pure function outer_divide(x, y) result(out)
! R outer(x, y, "/"): all pairwise quotients.
real(kind=dp), intent(in) :: x(:), y(:)
real(kind=dp), allocatable :: out(:,:)
integer :: i, j
allocate(out(size(x), size(y)))
do i = 1, size(x)
   do j = 1, size(y)
      out(i, j) = x(i) / y(j)
   end do
end do
end function outer_divide

pure function outer_power(x, y) result(out)
! R outer(x, y, "^"): all pairwise powers.
real(kind=dp), intent(in) :: x(:), y(:)
real(kind=dp), allocatable :: out(:,:)
integer :: i, j
allocate(out(size(x), size(y)))
do i = 1, size(x)
   do j = 1, size(y)
      out(i, j) = x(i) ** y(j)
   end do
end do
end function outer_power

pure function cumprod_real(x) result(out)
! Return cumulative products of a real vector.
real(kind=dp), intent(in) :: x(:) ! input vector
real(kind=dp), allocatable :: out(:)
integer :: i
allocate(out(size(x)))
if (size(x) <= 0) return
out(1) = x(1)
do i = 2, size(x)
   out(i) = out(i - 1) * x(i)
end do
end function cumprod_real

pure function cumprod_int(x) result(out)
! Return cumulative products of an integer vector.
integer, intent(in) :: x(:) ! input vector
integer, allocatable :: out(:)
integer :: i
allocate(out(size(x)))
if (size(x) <= 0) return
out(1) = x(1)
do i = 2, size(x)
   if (out(i - 1) == -huge(0) .or. x(i) == -huge(0)) then
      out(i) = -huge(0)
   else
      out(i) = out(i - 1) * x(i)
   end if
end do
end function cumprod_int

pure function cummax_real(x) result(out)
! Return cumulative maxima of a real vector.
real(kind=dp), intent(in) :: x(:) ! input vector
real(kind=dp), allocatable :: out(:)
integer :: i
allocate(out(size(x)))
if (size(x) <= 0) return
out(1) = x(1)
do i = 2, size(x)
   out(i) = max(out(i - 1), x(i))
end do
end function cummax_real

pure function cummax_int(x) result(out)
! Return cumulative maxima of an integer vector.
integer, intent(in) :: x(:) ! input vector
integer, allocatable :: out(:)
integer :: i
allocate(out(size(x)))
if (size(x) <= 0) return
out(1) = x(1)
do i = 2, size(x)
   if (out(i - 1) == -huge(0) .or. x(i) == -huge(0)) then
      out(i) = -huge(0)
   else
      out(i) = max(out(i - 1), x(i))
   end if
end do
end function cummax_int

pure function diff_real(x) result(out)
! First differences of a real vector.
real(kind=dp), intent(in) :: x(:) ! input vector
real(kind=dp), allocatable :: out(:)
integer :: i, n
n = size(x)
allocate(out(max(0, n - 1)))
do i = 1, n - 1
   out(i) = diff_subtract_real(x(i + 1), x(i))
end do
end function diff_real

pure function diff_mat_real(x) result(out)
! First row differences of a real matrix, matching R diff() on matrices.
real(kind=dp), intent(in) :: x(:,:) ! input matrix
real(kind=dp), allocatable :: out(:,:)
integer :: i, j, n, p
n = size(x, 1)
p = size(x, 2)
allocate(out(max(0, n - 1), p))
do j = 1, p
   do i = 1, n - 1
      out(i, j) = diff_subtract_real(x(i + 1, j), x(i, j))
   end do
end do
end function diff_mat_real

pure elemental function diff_subtract_real(a, b) result(out)
! Subtract while preserving R's distinction between NA and ordinary NaN.
real(kind=dp), intent(in) :: a, b
real(kind=dp) :: out
if (r_is_na_payload(a)) then
   out = r_na_real()
else if (a /= a) then
   out = a
else if (r_is_na_payload(b)) then
   out = r_na_real()
else if (b /= b) then
   out = b
else
   out = a - b
end if
end function diff_subtract_real

pure function diff_int(x) result(out)
! First differences of an integer vector.
integer, intent(in) :: x(:) ! input vector
integer, allocatable :: out(:)
integer :: i, n
n = size(x)
allocate(out(max(0, n - 1)))
do i = 1, n - 1
   if (x(i + 1) == -huge(0) .or. x(i) == -huge(0)) then
      out(i) = -huge(0)
   else
      out(i) = x(i + 1) - x(i)
   end if
end do
end function diff_int

pure function diag_mat_real(a) result(out)
! Return diagonal of a real matrix.
real(kind=dp), intent(in) :: a(:,:)
real(kind=dp), allocatable :: out(:)
integer :: i, n
n = min(size(a, 1), size(a, 2))
allocate(out(n))
do i = 1, n
   out(i) = a(i, i)
end do
end function diag_mat_real

pure function diag_vec_real(v) result(out)
! Create diagonal real matrix from a real vector.
real(kind=dp), intent(in) :: v(:)
real(kind=dp), allocatable :: out(:,:)
integer :: i, n
n = size(v)
allocate(out(n, n))
out = 0.0_dp
do i = 1, n
   out(i, i) = v(i)
end do
end function diag_vec_real

pure function diag_vec_real_n(v, n) result(out)
! Create an n by n diagonal real matrix, recycling the input vector.
real(kind=dp), intent(in) :: v(:)
integer, intent(in) :: n
real(kind=dp), allocatable :: out(:,:)
integer :: i, nn
nn = max(0, n)
allocate(out(nn, nn))
out = 0.0_dp
if (size(v) <= 0) return
do i = 1, nn
   out(i, i) = v(modulo(i - 1, size(v)) + 1)
end do
end function diag_vec_real_n

pure function diag_mat_complex(a) result(out)
! Return diagonal of a complex matrix.
complex(kind=dp), intent(in) :: a(:,:)
complex(kind=dp), allocatable :: out(:)
integer :: i, n
n = min(size(a, 1), size(a, 2))
allocate(out(n))
do i = 1, n
   out(i) = a(i, i)
end do
end function diag_mat_complex

pure function diag_vec_complex(v) result(out)
! Create diagonal complex matrix from a complex vector.
complex(kind=dp), intent(in) :: v(:)
complex(kind=dp), allocatable :: out(:,:)
integer :: i, n
n = size(v)
allocate(out(n, n))
out = cmplx(0.0_dp, 0.0_dp, kind=dp)
do i = 1, n
   out(i, i) = v(i)
end do
end function diag_vec_complex

pure function diag_vec_complex_n(v, n) result(out)
! Create an n by n diagonal complex matrix, recycling the input vector.
complex(kind=dp), intent(in) :: v(:)
integer, intent(in) :: n
complex(kind=dp), allocatable :: out(:,:)
integer :: i, nn
nn = max(0, n)
allocate(out(nn, nn))
out = cmplx(0.0_dp, 0.0_dp, kind=dp)
if (size(v) <= 0) return
do i = 1, nn
   out(i, i) = v(modulo(i - 1, size(v)) + 1)
end do
end function diag_vec_complex_n

pure function diag_mat_int(a) result(out)
! Return diagonal of an integer matrix.
integer, intent(in) :: a(:,:)
integer, allocatable :: out(:)
integer :: i, n
n = min(size(a, 1), size(a, 2))
allocate(out(n))
do i = 1, n
   out(i) = a(i, i)
end do
end function diag_mat_int

pure function diag_vec_int(v) result(out)
! Create diagonal integer matrix from an integer vector.
integer, intent(in) :: v(:)
integer, allocatable :: out(:,:)
integer :: i, n
n = size(v)
allocate(out(n, n))
out = 0
do i = 1, n
   out(i, i) = v(i)
end do
end function diag_vec_int

pure function diag_vec_int_n(v, n) result(out)
! Create an n by n diagonal integer matrix, recycling the input vector.
integer, intent(in) :: v(:)
integer, intent(in) :: n
integer, allocatable :: out(:,:)
integer :: i, nn
nn = max(0, n)
allocate(out(nn, nn))
out = 0
if (size(v) <= 0) return
do i = 1, nn
   out(i, i) = v(modulo(i - 1, size(v)) + 1)
end do
end function diag_vec_int_n

pure function diag_scalar_int(n) result(out)
! Create an n by n integer identity matrix, matching R diag(n).
integer, intent(in) :: n
integer, allocatable :: out(:,:)
integer :: i
allocate(out(n, n))
out = 0
do i = 1, n
   out(i, i) = 1
end do
end function diag_scalar_int

pure function diag_scalar_real_n(x, n) result(out)
! Create an n by n real diagonal matrix with scalar value x.
real(kind=dp), intent(in) :: x ! input values
integer, intent(in) :: n
real(kind=dp), allocatable :: out(:,:)
integer :: i
allocate(out(n, n))
out = 0.0_dp
do i = 1, n
   out(i, i) = x
end do
end function diag_scalar_real_n

pure function toeplitz(x) result(out)
! Symmetric Toeplitz matrix from first column/row vector x.
real(kind=dp), intent(in) :: x(:)
real(kind=dp), allocatable :: out(:,:)
integer :: i, j, n
n = size(x)
allocate(out(n, n))
do j = 1, n
   do i = 1, n
      out(i, j) = x(abs(i - j) + 1)
   end do
end do
end function toeplitz

pure function chol_real(a) result(r)
! Upper-triangular Cholesky factor for a symmetric positive-definite matrix.
real(kind=dp), intent(in) :: a(:,:)
real(kind=dp), allocatable :: r(:,:)
real(kind=dp) :: s
integer :: j, k, n
n = size(a, 1)
allocate(r(n, n))
r = 0.0_dp
do j = 1, n
   s = a(j, j)
   if (j > 1) s = s - sum(r(1:j-1, j)**2)
   r(j, j) = sqrt(max(s, 0.0_dp))
   do k = j + 1, n
      s = a(j, k)
      if (j > 1) s = s - sum(r(1:j-1, j) * r(1:j-1, k))
      if (r(j, j) > 0.0_dp) r(j, k) = s / r(j, j)
   end do
end do
end function chol_real

pure function chol_int(a) result(r)
! Upper-triangular Cholesky factor for an integer matrix coerced to real.
integer, intent(in) :: a(:,:)
real(kind=dp), allocatable :: r(:,:)
r = chol_real(real(a, kind=dp))
end function chol_int

pure function chol2inv_real(r, size) result(out)
! Runtime helper for R-compatible chol2inv real.
real(kind=dp), intent(in) :: r(:,:) ! input matrix
integer, intent(in), optional :: size
real(kind=dp), allocatable :: out(:,:), a(:,:)
integer :: n
n = min(ubound(r, 1), ubound(r, 2))
if (present(size)) n = min(n, size)
allocate(a(n, n))
a = matmul(transpose(r(1:n, 1:n)), r(1:n, 1:n))
out = solve_real_mat(a, real(diag(n), kind=dp))
end function chol2inv_real

pure function chol2inv_int(r, size) result(out)
! Runtime helper for R-compatible chol2inv int.
integer, intent(in) :: r(:,:) ! Cholesky factor
integer, intent(in), optional :: size ! leading dimension to invert
real(kind=dp), allocatable :: out(:,:)
if (present(size)) then
   out = chol2inv_real(real(r, kind=dp), size)
else
   out = chol2inv_real(real(r, kind=dp))
end if
end function chol2inv_int

pure function forwardsolve_mat(l, b, transpose) result(x)
! Solve L x = b for lower-triangular L; transpose=.true. solves L^T x = b.
real(kind=dp), intent(in) :: l(:,:) ! lower-triangular coefficient matrix
real(kind=dp), intent(in) :: b(:,:) ! right-hand side matrix
logical, intent(in), optional :: transpose ! solve transposed system flag
real(kind=dp), allocatable :: x(:,:)
logical :: tr
integer :: i, j, n, m
real(kind=dp) :: s
n = size(l, 1)
m = size(b, 2)
allocate(x(n, m))
x = 0.0_dp
tr = .false.
if (present(transpose)) tr = transpose
if (tr) then
   do j = 1, m
      do i = n, 1, -1
         s = b(i, j)
         if (i < n) s = s - sum(l(i+1:n, i) * x(i+1:n, j))
         x(i, j) = s / l(i, i)
      end do
   end do
else
   do j = 1, m
      do i = 1, n
         s = b(i, j)
         if (i > 1) s = s - sum(l(i, 1:i-1) * x(1:i-1, j))
         x(i, j) = s / l(i, i)
      end do
   end do
end if
end function forwardsolve_mat

pure function forwardsolve_vec(l, b, transpose) result(x)
! Solve L x = b for a vector RHS.
real(kind=dp), intent(in) :: l(:,:) ! lower-triangular coefficient matrix
real(kind=dp), intent(in) :: b(:) ! right-hand side vector
logical, intent(in), optional :: transpose ! solve transposed system flag
real(kind=dp), allocatable :: x(:)
real(kind=dp), allocatable :: bm(:,:), xm(:,:)
allocate(bm(size(b), 1))
bm(:, 1) = b
if (present(transpose)) then
   xm = forwardsolve_mat(l, bm, transpose=transpose)
else
   xm = forwardsolve_mat(l, bm)
end if
x = xm(:, 1)
end function forwardsolve_vec

pure function forwardsolve_vec_i_r(l, b, transpose) result(x)
! Solve a lower-triangular system with mixed numeric input kinds.
integer, intent(in) :: l(:,:) ! lower-triangular coefficient matrix
real(kind=dp), intent(in) :: b(:) ! right-hand side vector
logical, intent(in), optional :: transpose ! solve transposed system flag
real(kind=dp), allocatable :: x(:)
if (present(transpose)) then
   x = forwardsolve_vec(real(l, kind=dp), b, transpose=transpose)
else
   x = forwardsolve_vec(real(l, kind=dp), b)
end if
end function forwardsolve_vec_i_r

pure function forwardsolve_vec_i_i(l, b, transpose) result(x)
! Solve a lower-triangular system with mixed numeric input kinds.
integer, intent(in) :: l(:,:) ! lower-triangular coefficient matrix
integer, intent(in) :: b(:) ! right-hand side vector
logical, intent(in), optional :: transpose ! solve transposed system flag
real(kind=dp), allocatable :: x(:)
if (present(transpose)) then
   x = forwardsolve_vec(real(l, kind=dp), real(b, kind=dp), transpose=transpose)
else
   x = forwardsolve_vec(real(l, kind=dp), real(b, kind=dp))
end if
end function forwardsolve_vec_i_i

pure function forwardsolve_mat_i_r(l, b, transpose) result(x)
! Solve a lower-triangular system with mixed numeric input kinds.
integer, intent(in) :: l(:,:) ! lower-triangular coefficient matrix
real(kind=dp), intent(in) :: b(:,:) ! right-hand side matrix
logical, intent(in), optional :: transpose ! solve transposed system flag
real(kind=dp), allocatable :: x(:,:)
if (present(transpose)) then
   x = forwardsolve_mat(real(l, kind=dp), b, transpose=transpose)
else
   x = forwardsolve_mat(real(l, kind=dp), b)
end if
end function forwardsolve_mat_i_r

pure function forwardsolve_mat_r_i(l, b, transpose) result(x)
! Solve a lower-triangular system with mixed numeric input kinds.
real(kind=dp), intent(in) :: l(:,:) ! lower-triangular coefficient matrix
integer, intent(in) :: b(:,:) ! right-hand side matrix
logical, intent(in), optional :: transpose ! solve transposed system flag
real(kind=dp), allocatable :: x(:,:)
if (present(transpose)) then
   x = forwardsolve_mat(l, real(b, kind=dp), transpose=transpose)
else
   x = forwardsolve_mat(l, real(b, kind=dp))
end if
end function forwardsolve_mat_r_i

pure function forwardsolve_mat_i_i(l, b, transpose) result(x)
! Solve a lower-triangular system with mixed numeric input kinds.
integer, intent(in) :: l(:,:) ! lower-triangular coefficient matrix
integer, intent(in) :: b(:,:) ! right-hand side matrix
logical, intent(in), optional :: transpose ! solve transposed system flag
real(kind=dp), allocatable :: x(:,:)
if (present(transpose)) then
   x = forwardsolve_mat(real(l, kind=dp), real(b, kind=dp), transpose=transpose)
else
   x = forwardsolve_mat(real(l, kind=dp), real(b, kind=dp))
end if
end function forwardsolve_mat_i_i

pure function backsolve_mat(r, b, transpose) result(x)
! Solve R x = b for upper-triangular R; transpose=.true. solves R^T x = b.
real(kind=dp), intent(in) :: r(:,:) ! upper-triangular coefficient matrix
real(kind=dp), intent(in) :: b(:,:) ! right-hand side matrix
logical, intent(in), optional :: transpose ! solve transposed system flag
real(kind=dp), allocatable :: x(:,:)
logical :: tr
integer :: i, j, n, m
real(kind=dp) :: s
n = size(r, 1)
m = size(b, 2)
allocate(x(n, m))
x = 0.0_dp
tr = .false.
if (present(transpose)) tr = transpose
if (tr) then
   do j = 1, m
      do i = 1, n
         s = b(i, j)
         if (i > 1) s = s - sum(r(1:i-1, i) * x(1:i-1, j))
         x(i, j) = s / r(i, i)
      end do
   end do
else
   do j = 1, m
      do i = n, 1, -1
         s = b(i, j)
         if (i < n) s = s - sum(r(i, i+1:n) * x(i+1:n, j))
         x(i, j) = s / r(i, i)
      end do
   end do
end if
end function backsolve_mat

pure function backsolve_vec(r, b, transpose) result(x)
! Solve R x = b for a vector RHS.
real(kind=dp), intent(in) :: r(:,:) ! upper-triangular coefficient matrix
real(kind=dp), intent(in) :: b(:) ! right-hand side vector
logical, intent(in), optional :: transpose ! solve transposed system flag
real(kind=dp), allocatable :: x(:)
real(kind=dp), allocatable :: bm(:,:), xm(:,:)
allocate(bm(size(b), 1))
bm(:, 1) = b
if (present(transpose)) then
   xm = backsolve_mat(r, bm, transpose=transpose)
else
   xm = backsolve_mat(r, bm)
end if
x = xm(:, 1)
end function backsolve_vec

pure function backsolve_vec_i_r(r, b, transpose) result(x)
! Solve an upper-triangular system with mixed numeric input kinds.
integer, intent(in) :: r(:,:) ! upper-triangular coefficient matrix
real(kind=dp), intent(in) :: b(:) ! right-hand side vector
logical, intent(in), optional :: transpose ! solve transposed system flag
real(kind=dp), allocatable :: x(:)
if (present(transpose)) then
   x = backsolve_vec(real(r, kind=dp), b, transpose=transpose)
else
   x = backsolve_vec(real(r, kind=dp), b)
end if
end function backsolve_vec_i_r

pure function backsolve_vec_i_i(r, b, transpose) result(x)
! Solve an upper-triangular system with mixed numeric input kinds.
integer, intent(in) :: r(:,:) ! upper-triangular coefficient matrix
integer, intent(in) :: b(:) ! right-hand side vector
logical, intent(in), optional :: transpose ! solve transposed system flag
real(kind=dp), allocatable :: x(:)
if (present(transpose)) then
   x = backsolve_vec(real(r, kind=dp), real(b, kind=dp), transpose=transpose)
else
   x = backsolve_vec(real(r, kind=dp), real(b, kind=dp))
end if
end function backsolve_vec_i_i

pure function backsolve_mat_i_r(r, b, transpose) result(x)
! Solve an upper-triangular system with mixed numeric input kinds.
integer, intent(in) :: r(:,:) ! upper-triangular coefficient matrix
real(kind=dp), intent(in) :: b(:,:) ! right-hand side matrix
logical, intent(in), optional :: transpose ! solve transposed system flag
real(kind=dp), allocatable :: x(:,:)
if (present(transpose)) then
   x = backsolve_mat(real(r, kind=dp), b, transpose=transpose)
else
   x = backsolve_mat(real(r, kind=dp), b)
end if
end function backsolve_mat_i_r

pure function backsolve_mat_r_i(r, b, transpose) result(x)
! Solve an upper-triangular system with mixed numeric input kinds.
real(kind=dp), intent(in) :: r(:,:) ! upper-triangular coefficient matrix
integer, intent(in) :: b(:,:) ! right-hand side matrix
logical, intent(in), optional :: transpose ! solve transposed system flag
real(kind=dp), allocatable :: x(:,:)
if (present(transpose)) then
   x = backsolve_mat(r, real(b, kind=dp), transpose=transpose)
else
   x = backsolve_mat(r, real(b, kind=dp))
end if
end function backsolve_mat_r_i

pure function backsolve_mat_i_i(r, b, transpose) result(x)
! Solve an upper-triangular system with mixed numeric input kinds.
integer, intent(in) :: r(:,:) ! upper-triangular coefficient matrix
integer, intent(in) :: b(:,:) ! right-hand side matrix
logical, intent(in), optional :: transpose ! solve transposed system flag
real(kind=dp), allocatable :: x(:,:)
if (present(transpose)) then
   x = backsolve_mat(real(r, kind=dp), real(b, kind=dp), transpose=transpose)
else
   x = backsolve_mat(real(r, kind=dp), real(b, kind=dp))
end if
end function backsolve_mat_i_i

pure elemental integer function nchar(s) result(out)
! Return character length (R-like nchar scalar subset).
character(len=*), intent(in) :: s ! string whose trimmed length is returned
out = len_trim(s)
end function nchar

pure elemental logical function char_ends_with(x, suffix) result(out)
! Return whether x ends with suffix, with scalar/array elemental broadcasting.
character(len=*), intent(in) :: x, suffix
integer :: nx, ns
nx = len_trim(x)
ns = len_trim(suffix)
out = nx >= ns
if (out .and. ns > 0) out = x(nx - ns + 1:nx) == suffix(1:ns)
end function char_ends_with

pure function char_join_char(x, sep) result(out)
! Implement R-like character helper char_join.
character(len=*), intent(in) :: x(:) ! strings to join
character(len=*), intent(in) :: sep ! separator between strings
character(len=:), allocatable :: out
integer :: i, total
total = 0
do i = 1, size(x)
   total = total + len_trim(x(i))
end do
if (size(x) > 1) total = total + (size(x) - 1) * len(sep)
allocate(character(len=max(0, total)) :: out)
out = ""
do i = 1, size(x)
   if (i > 1) out = out // sep
   out = out // trim(x(i))
end do
end function char_join_char

pure function char_join_int(x, sep) result(out)
! Join integer values after applying R-like character coercion.
integer, intent(in) :: x(:)
character(len=*), intent(in) :: sep
character(len=:), allocatable :: out
character(len=:), allocatable :: item
integer :: i, total
total = 0
do i = 1, size(x)
   item = int_to_string(x(i))
   total = total + len(item)
end do
if (size(x) > 1) total = total + (size(x) - 1) * len(sep)
allocate(character(len=max(0, total)) :: out)
out = ""
do i = 1, size(x)
   if (i > 1) out = out // sep
   out = out // int_to_string(x(i))
end do
end function char_join_int

pure function r_paste0_real(prefix, x) result(out)
! Prefix each numeric value with a fixed string, matching paste0(prefix, x).
character(len=*), intent(in) :: prefix ! prefix string
real(kind=dp), intent(in) :: x(:) ! values to convert and prefix
character(len=:), allocatable :: out(:)
character(len=64) :: buf
integer :: i, k
real(kind=dp) :: tol
allocate(character(len=max(1, len(prefix) + 64)) :: out(size(x)))
do i = 1, size(x)
   if (ieee_is_finite(x(i)) .and. abs(x(i)) <= real(huge(0), kind=dp)) then
      k = nint(x(i))
      tol = 1000.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(x(i)))
      if (abs(x(i) - real(k, kind=dp)) <= tol) then
         write(buf, "(i0)") k
      else
         write(buf, "(f0.10)") x(i)
         do while (index(buf, ".") > 0 .and. buf(len_trim(buf):len_trim(buf)) == "0")
            buf(len_trim(buf):len_trim(buf)) = " "
         end do
         if (buf(len_trim(buf):len_trim(buf)) == ".") buf(len_trim(buf):len_trim(buf)) = " "
         buf = adjustl(buf)
         if (len_trim(buf) >= 1 .and. buf(1:1) == ".") buf = "0" // trim(buf)
         if (len_trim(buf) >= 2 .and. buf(1:2) == "-.") buf = "-0" // trim(buf(2:))
      end if
   else
      write(buf, "(g0)") x(i)
   end if
   out(i) = prefix // trim(adjustl(buf))
end do
end function r_paste0_real

pure function r_paste0_int(prefix, x) result(out)
! Prefix each integer value with a fixed string, matching paste0(prefix, x).
character(len=*), intent(in) :: prefix ! prefix string
integer, intent(in) :: x(:) ! values to convert and prefix
character(len=:), allocatable :: out(:)
character(len=64) :: buf
integer :: i
allocate(character(len=max(1, len(prefix) + 64)) :: out(size(x)))
do i = 1, size(x)
   write(buf, "(i0)") x(i)
   out(i) = prefix // trim(adjustl(buf))
end do
end function r_paste0_int

function file_exists(path) result(out)
! Implement scalar R-like file.exists(path).
character(len=*), intent(in) :: path
logical :: out
integer :: stat
inquire(file=trim(path), exist=out)
if (out) return
if (is_windows_path_env()) then
   call execute_command_line('cmd /c if exist "' // trim(path) // '" (exit /b 0) else (exit /b 1)', &
      wait=.true., exitstat=stat)
else
   call execute_command_line('test -e "' // trim(path) // '"', wait=.true., exitstat=stat)
end if
out = stat == 0
end function file_exists

function file_create(path) result(out)
! Implement scalar R-like file.create(path).
character(len=*), intent(in) :: path
logical :: out
integer :: unit, ios
open(newunit=unit, file=trim(path), status="replace", action="write", iostat=ios)
if (ios == 0) close(unit, iostat=ios)
out = ios == 0
end function file_create

function file_remove(path) result(out)
! Implement scalar R-like file.remove(path).
character(len=*), intent(in) :: path
logical :: out
integer :: unit, ios, ios_close, stat
open(newunit=unit, file=trim(path), status="old", action="readwrite", iostat=ios)
if (ios /= 0) then
   if (is_windows_path_env()) then
      call execute_command_line('cmd /c rmdir "' // trim(path) // '" >nul 2>nul', wait=.true., exitstat=stat)
   else
      call execute_command_line('rmdir "' // trim(path) // '" >/dev/null 2>&1', wait=.true., exitstat=stat)
   end if
   out = stat == 0
   return
end if
close(unit, status="delete", iostat=ios_close)
out = ios_close == 0
end function file_remove

function file_info_scalar(path) result(out)
! Implement a compact scalar subset of R file.info(path).
character(len=*), intent(in) :: path
type(file_info_t) :: out
logical :: exists_file
integer :: stat
integer(int64) :: sz
out%path = trim(path)
inquire(file=trim(path), exist=exists_file, size=sz)
if (exists_file .and. sz >= 0_int64) then
   out%size = real(sz, kind=dp)
else if (file_exists(path)) then
   out%size = 0.0_dp
else
   out%size = ieee_value(0.0_dp, ieee_quiet_nan)
end if
if (is_windows_path_env()) then
   call execute_command_line('cmd /c if exist "' // trim(path) // '\*" (exit /b 0) else (exit /b 1)', &
      wait=.true., exitstat=stat)
else
   call execute_command_line('test -d "' // trim(path) // '"', wait=.true., exitstat=stat)
end if
out%isdir = stat == 0
out%mode = merge("directory", "file     ", out%isdir)
if (out%isdir) then
   out%exe = "no"
else if (is_windows_path_env()) then
   if (index(tolower(trim(path)), ".exe") > 0 .or. index(tolower(trim(path)), ".bat") > 0 .or. &
      & index(tolower(trim(path)), ".cmd") > 0 .or. index(tolower(trim(path)), ".com") > 0) then
      out%exe = "yes"
   else
      out%exe = "no"
   end if
else
   call execute_command_line('test -x "' // trim(path) // '"', wait=.true., exitstat=stat)
   out%exe = merge("yes", "no ", stat == 0)
end if
out%mtime = ""
out%ctime = ""
out%atime = ""
end function file_info_scalar

function file_size(path) result(out)
character(len=*), intent(in) :: path
real(kind=dp) :: out
type(file_info_t) :: info
info = file_info_scalar(path)
out = info%size
end function file_size

function file_path_value(path) result(out)
character(len=*), intent(in) :: path
character(len=:), allocatable :: out
out = trim(path)
end function file_path_value

pure function file_extension(path) result(out)
character(len=*), intent(in) :: path
character(len=:), allocatable :: out
integer :: dot_pos, i
logical :: valid_ext
dot_pos = 0
do i = len_trim(path), 1, -1
   if (path(i:i) == ".") then
      dot_pos = i
      exit
   end if
   if (path(i:i) == "/" .or. path(i:i) == achar(92)) exit
end do
valid_ext = dot_pos > 0 .and. dot_pos < len_trim(path)
if (valid_ext) then
   do i = dot_pos + 1, len_trim(path)
      if (.not. ((path(i:i) >= "A" .and. path(i:i) <= "Z") .or. &
         (path(i:i) >= "a" .and. path(i:i) <= "z") .or. &
         (path(i:i) >= "0" .and. path(i:i) <= "9"))) then
         valid_ext = .false.
         exit
      end if
   end do
end if
if (valid_ext) then
   out = path(dot_pos:len_trim(path))
else
   out = ""
end if
end function file_extension

function file_mode(path) result(out)
character(len=*), intent(in) :: path
character(len=:), allocatable :: out
type(file_info_t) :: info
info = file_info_scalar(path)
out = info%mode
end function file_mode

function file_mtime(path) result(out)
character(len=*), intent(in) :: path
character(len=:), allocatable :: out
type(file_info_t) :: info
info = file_info_scalar(path)
out = info%mtime
end function file_mtime

function file_ctime(path) result(out)
character(len=*), intent(in) :: path
character(len=:), allocatable :: out
type(file_info_t) :: info
info = file_info_scalar(path)
out = info%ctime
end function file_ctime

function file_atime(path) result(out)
character(len=*), intent(in) :: path
character(len=:), allocatable :: out
type(file_info_t) :: info
info = file_info_scalar(path)
out = info%atime
end function file_atime

function file_exe(path) result(out)
character(len=*), intent(in) :: path
character(len=:), allocatable :: out
type(file_info_t) :: info
info = file_info_scalar(path)
out = info%exe
end function file_exe

function file_info_vector(path) result(out)
! Implement a compact vector subset of R file.info(paths).
character(len=*), intent(in) :: path(:)
type(file_info_t), allocatable :: out(:)
integer :: i
allocate(out(size(path)))
do i = 1, size(path)
   out(i) = file_info_scalar(path(i))
end do
end function file_info_vector

function file_isdir(path) result(out)
! Return file.info(path)$isdir without component selection on a function result.
character(len=*), intent(in) :: path
logical :: out
type(file_info_t) :: info
info = file_info_scalar(path)
out = info%isdir
end function file_isdir

subroutine print_file_info_scalar(x)
! Print a compact R-like file.info() summary.
type(file_info_t), intent(in) :: x
integer :: name_w
name_w = max(4, len_trim(x%path))
write(*,'(a)') repeat(" ", name_w + 1) // "      size isdir mode      exe"
write(*,'(a,1x,f10.0,1x,l5,1x,a4,6x,a3)') repeat(" ", max(0, name_w - len_trim(x%path))) // trim(x%path), &
   & x%size, x%isdir, trim(x%mode), trim(x%exe)
end subroutine print_file_info_scalar

subroutine print_file_info_vector(x)
! Print a compact R-like file.info() summary for multiple files.
type(file_info_t), intent(in) :: x(:)
integer :: i, name_w
name_w = 4
do i = 1, size(x)
   name_w = max(name_w, len_trim(x(i)%path))
end do
write(*,'(a)') repeat(" ", name_w + 1) // "      size isdir mode      exe"
do i = 1, size(x)
   write(*,'(a,1x,f10.0,1x,l5,1x,a4,6x,a3)') repeat(" ", max(0, name_w - len_trim(x(i)%path))) // trim(x(i)%path), &
      & x(i)%size, x(i)%isdir, trim(x(i)%mode), trim(x(i)%exe)
end do
end subroutine print_file_info_vector

function dir_exists_scalar(path) result(out)
! Implement scalar R-like dir.exists(path).
character(len=*), intent(in) :: path
logical :: out
integer :: stat
if (is_windows_path_env()) then
   call execute_command_line('cmd /c if exist "' // trim(path) // '\*" (exit /b 0) else (exit /b 1)', &
      wait=.true., exitstat=stat)
else
   call execute_command_line('test -d "' // trim(path) // '"', wait=.true., exitstat=stat)
end if
out = stat == 0
end function dir_exists_scalar

function dir_exists_vector(path) result(out)
! Implement vectorized R-like dir.exists(path).
character(len=*), intent(in) :: path(:)
logical, allocatable :: out(:)
integer :: i
allocate(out(size(path)))
do i = 1, size(path)
   out(i) = dir_exists_scalar(path(i))
end do
end function dir_exists_vector

function dir_create(path, recursive) result(out)
! Implement scalar R-like dir.create(path).
character(len=*), intent(in) :: path
logical, intent(in), optional :: recursive
logical :: out, recur
integer :: stat
recur = .false.
if (present(recursive)) recur = recursive
if (file_exists(path)) then
   out = .false.
   return
end if
if (is_windows_path_env()) then
   if (recur) then
      call execute_command_line('cmd /c mkdir "' // trim(path) // '" >nul 2>nul', wait=.true., exitstat=stat)
   else
      call execute_command_line('cmd /c mkdir "' // trim(path) // '" >nul 2>nul', wait=.true., exitstat=stat)
   end if
else
   if (recur) then
      call execute_command_line('mkdir -p "' // trim(path) // '" >/dev/null 2>&1', wait=.true., exitstat=stat)
   else
      call execute_command_line('mkdir "' // trim(path) // '" >/dev/null 2>&1', wait=.true., exitstat=stat)
   end if
end if
out = stat == 0 .and. file_exists(path)
end function dir_create

function getwd() result(out)
! Implement scalar R-like getwd().
character(len=:), allocatable :: out
character(len=:), allocatable :: tmp
character(len=4096) :: buf
integer :: unit, ios, stat
call random_seed()
call random_number_list_files_tmp(tmp)
if (is_windows_path_env()) then
   call execute_command_line('cmd /c cd > "' // tmp // '"', wait=.true., exitstat=stat)
else
   call execute_command_line('pwd > "' // tmp // '"', wait=.true., exitstat=stat)
end if
open(newunit=unit, file=tmp, status="old", action="read", iostat=ios)
if (ios /= 0) then
   out = "."
   return
end if
read(unit, "(a)", iostat=ios) buf
close(unit, status="delete")
if (ios == 0 .and. len_trim(buf) > 0) then
   out = trim(buf)
else
   out = "."
end if
end function getwd

function tempfile(pattern) result(out)
! Return a simple temporary path with an R-like prefix.
character(len=*), intent(in), optional :: pattern
character(len=:), allocatable :: out
character(len=4096) :: root
character(len=64) :: id
character(len=1) :: sep
real(kind=dp) :: u
integer :: n, stat, got
call get_environment_variable("TEMP", root, length=got, status=stat)
if (stat /= 0 .or. got <= 0) call get_environment_variable("TMP", root, length=got, status=stat)
if (stat /= 0 .or. got <= 0) then
   root = "."
   got = 1
end if
sep = merge("\", "/", is_windows_path_env())
call random_number(u)
n = max(0, int(u * 1000000000.0_dp))
write(id, "(i0)") n
if (present(pattern)) then
   out = trim(root(1:got)) // sep // trim(pattern) // trim(id)
else
   out = trim(root(1:got)) // sep // "file" // trim(id)
end if
end function tempfile

function file_path(a, b) result(out)
! Join two path components using the platform separator.
character(len=*), intent(in) :: a, b
character(len=:), allocatable :: out
character(len=1) :: sep
sep = merge("\", "/", is_windows_path_env())
if (len_trim(a) == 0) then
   out = trim(b)
else if (a(len_trim(a):len_trim(a)) == "/" .or. a(len_trim(a):len_trim(a)) == "\") then
   out = trim(a) // trim(b)
else
   out = trim(a) // sep // trim(b)
end if
end function file_path

function list_files(path, pattern, full_names, recursive) result(out)
! Implement R-like character helper list_files.
! Defaults: path=".", pattern="", full_names=.false., recursive=.false.
character(len=*), intent(in), optional :: path ! directory to list
character(len=*), intent(in), optional :: pattern ! simple filename filter
logical, intent(in), optional :: full_names ! return full paths flag
logical, intent(in), optional :: recursive ! descend into subdirectories flag
character(len=:), allocatable :: out(:)
character(len=:), allocatable :: p, pat, tmp, cmd, line, base
character(len=4096) :: buf
logical :: fnames, recur, keep
integer :: unit, ios, stat, n, maxlen, i, slash
p = "."
if (present(path)) p = trim(path)
pat = ""
if (present(pattern)) pat = trim(pattern)
fnames = .false.
if (present(full_names)) fnames = full_names
recur = .false.
if (present(recursive)) recur = recursive
call random_seed()
call random_number_list_files_tmp(tmp)
if (is_windows_path_env()) then
   if (recur) then
      cmd = 'cmd /c dir /b /s /a-d "' // p // '" > "' // tmp // '" 2>nul'
   else
      cmd = 'cmd /c dir /b /a-d "' // p // '" > "' // tmp // '" 2>nul'
   end if
else
   if (recur) then
      cmd = 'find "' // p // '" -type f > "' // tmp // '" 2>/dev/null'
   else
      cmd = 'find "' // p // '" -maxdepth 1 -type f > "' // tmp // '" 2>/dev/null'
   end if
end if
call execute_command_line(cmd, wait=.true., exitstat=stat)
n = 0
maxlen = 1
open(newunit=unit, file=tmp, status="old", action="read", iostat=ios)
if (ios /= 0) then
   allocate(character(len=1) :: out(0))
   return
end if
do
   read(unit, "(a)", iostat=ios) buf
   if (ios /= 0) exit
   line = trim(buf)
   base = list_files_basename(line)
   keep = pat == "" .or. list_files_pattern_match(base, pat)
   if (keep) then
      n = n + 1
      if (fnames) then
         if (is_windows_path_env() .and. .not. recur) then
            maxlen = max(maxlen, len_trim(p) + 1 + len_trim(line))
         else
            maxlen = max(maxlen, len_trim(line))
         end if
      else if (recur) then
         maxlen = max(maxlen, len_trim(line))
      else
         maxlen = max(maxlen, len_trim(base))
      end if
   end if
end do
rewind(unit)
if (n <= 0) then
   allocate(character(len=1) :: out(0))
   close(unit, status="delete")
   return
end if
allocate(character(len=maxlen) :: out(n))
i = 0
do
   read(unit, "(a)", iostat=ios) buf
   if (ios /= 0) exit
   line = trim(buf)
   base = list_files_basename(line)
   keep = pat == "" .or. list_files_pattern_match(base, pat)
   if (keep) then
      i = i + 1
      out(i) = ""
      if (fnames) then
         if (is_windows_path_env() .and. .not. recur) then
            if (len_trim(p) > 0 .and. (p(len_trim(p):len_trim(p)) == "/" .or. p(len_trim(p):len_trim(p)) == "\")) then
               out(i) = trim(p) // trim(line)
            else
               out(i) = trim(p) // "/" // trim(line)
            end if
         else
            out(i) = trim(line)
         end if
      else if (recur .and. .not. is_windows_path_env()) then
         if (len_trim(p) > 0 .and. index(line, trim(p) // "/") == 1) then
            out(i) = line(len_trim(p) + 2:)
         else
            out(i) = trim(base)
         end if
      else if (recur .and. is_windows_path_env()) then
         slash = len_trim(p)
         if (slash > 0 .and. index(line, trim(p) // "\") == 1) then
            out(i) = line(slash + 2:)
         else
            out(i) = trim(base)
         end if
      else
         out(i) = trim(base)
      end if
   end if
end do
close(unit, status="delete")
end function list_files

subroutine random_number_list_files_tmp(tmp)
! Runtime helper for R-compatible random number list files tmp.
character(len=:), allocatable, intent(out) :: tmp ! generated temporary filename
real(kind=dp) :: u
integer :: k
call random_number(u)
k = max(1, int(u * 1000000000.0_dp))
tmp = "xr2f_list_files_" // int_to_string(k) // ".tmp"
end subroutine random_number_list_files_tmp

pure function int_to_string(i) result(out)
! Runtime helper for R-compatible int to string.
integer, intent(in) :: i ! value to format
character(len=:), allocatable :: out
character(len=32) :: buf
write(buf, "(i0)") i
out = trim(buf)
end function int_to_string

pure function real_to_string_f(x, digits) result(out)
! Runtime helper for fixed-format scalar sprintf real conversion.
real(kind=dp), intent(in) :: x
integer, intent(in) :: digits
character(len=:), allocatable :: out
character(len=128) :: buf
character(len=32) :: fmt
write(fmt, '("(f0.", i0, ")")') max(0, digits)
write(buf, fmt) x
out = trim(adjustl(buf))
if (len(out) >= 1 .and. out(1:1) == ".") out = "0" // out
if (len(out) >= 2 .and. out(1:2) == "-.") out = "-0" // out(2:)
end function real_to_string_f

pure function real_to_string_g_scalar(x, digits) result(out)
! Runtime helper for general-format scalar sprintf real conversion.
real(kind=dp), intent(in) :: x
integer, intent(in) :: digits
character(len=:), allocatable :: out
character(len=128) :: buf
character(len=32) :: fmt
if (digits > 0) then
   write(fmt, '("(g0.", i0, ")")') digits
else
   fmt = "(g0)"
end if
write(buf, fmt) x
out = trim(adjustl(buf))
end function real_to_string_g_scalar

pure function real_to_string_g_vector(x, digits) result(out)
real(kind=dp), intent(in) :: x(:)
integer, intent(in) :: digits
character(len=128), allocatable :: out(:)
character(len=128) :: buf
character(len=32) :: fmt
integer :: i
allocate(out(size(x)))
if (digits > 0) then
   write(fmt, '("(g0.", i0, ")")') digits
else
   fmt = "(g0)"
end if
do i = 1, size(x)
   write(buf, fmt) x(i)
   out(i) = trim(adjustl(buf))
end do
end function real_to_string_g_vector

pure function r_to_string_real(x) result(out)
real(kind=dp), intent(in) :: x
character(len=:), allocatable :: out
character(len=128) :: buf
integer :: dot
write(buf, "(g0.7)") x
out = trim(adjustl(buf))
if (index(out, "E") == 0 .and. index(out, "e") == 0) then
   dot = index(out, ".")
   if (dot > 0) then
      do while (len(out) > dot .and. out(len(out):len(out)) == "0")
         out = out(:len(out) - 1)
      end do
      if (len(out) == dot) out = out(:dot - 1)
   end if
end if
end function r_to_string_real

pure function r_substr(s, first, last) result(out)
character(len=*), intent(in) :: s
integer, intent(in) :: first, last
character(len=:), allocatable :: out
integer :: lo, hi
lo = max(1, first)
hi = min(len(s), last)
if (hi < lo) then
   out = ""
else
   out = s(lo:hi)
end if
end function r_substr

pure function r_substr_replace(s, first, last, value) result(out)
character(len=*), intent(in) :: s, value
integer, intent(in) :: first, last
character(len=:), allocatable :: out
integer :: lo, hi, nrep
out = s
lo = max(1, first)
hi = min(len(s), last)
nrep = min(max(0, hi - lo + 1), len_trim(value))
if (nrep > 0) out(lo:lo + nrep - 1) = value(1:nrep)
end function r_substr_replace

pure function ar_coef_names(nacf) result(out)
! Runtime helper for R-compatible ar coef names.
integer, intent(in) :: nacf ! number of autocorrelation coefficient names
character(len=:), allocatable :: out(:)
integer :: i, n
n = max(0, nacf)
allocate(character(len=32) :: out(n + 5))
out(1) = "order"
out(2) = "intercept"
do i = 1, n
   out(i + 2) = "phi" // int_to_string(i)
end do
out(n + 3) = "sigma2"
out(n + 4) = "aic"
out(n + 5) = "bic"
end function ar_coef_names

pure function lag_names(nlag) result(out)
! Runtime helper for R-compatible lag names.
integer, intent(in) :: nlag ! number of lag names
character(len=:), allocatable :: out(:)
integer :: i, n
n = max(0, nlag)
allocate(character(len=32) :: out(n))
do i = 1, n
   out(i) = "lag" // int_to_string(i)
end do
end function lag_names

function is_windows_path_env() result(out)
! Test the R-like predicate is_windows_path_env.
logical :: out
character(len=16) :: os
integer :: stat, n
call get_environment_variable("OS", os, length=n, status=stat)
out = stat == 0 .and. index(os(1:max(1,n)), "Windows") > 0
end function is_windows_path_env

pure function list_files_basename(s) result(out)
! Runtime helper for R-compatible list files basename.
character(len=*), intent(in) :: s ! path or filename
character(len=:), allocatable :: out
integer :: i, last
last = 0
do i = 1, len_trim(s)
   if (s(i:i) == "/" .or. s(i:i) == "\") last = i
end do
out = s(last + 1:len_trim(s))
end function list_files_basename

pure function list_files_pattern_match(name, pattern) result(out)
! Runtime helper for R-compatible list files pattern match.
! Pattern supports simple contains matching and leading/trailing "*" wildcards.
character(len=*), intent(in) :: name ! filename to test
character(len=*), intent(in) :: pattern ! simple filename pattern
logical :: out
character(len=:), allocatable :: pat, pat2, suffix, core, prefix
integer :: i, j, n
pat = trim(pattern)
allocate(character(len=max(1, len_trim(pat))) :: pat2)
pat2 = repeat(" ", len(pat2))
j = 0
do i = 1, len_trim(pat)
   if (pat(i:i) == "\") cycle
   j = j + 1
   pat2(j:j) = pat(i:i)
end do
if (j > 0) then
   pat = pat2(1:j)
else
   pat = ""
end if
n = len_trim(pat)
if (pat == "") then
   out = .true.
else if (n > 1 .and. pat(n:n) == "$") then
   core = pat(1:n-1)
   do while (len_trim(core) > 0 .and. (core(1:1) == "^" .or. core(1:1) == "*"))
      core = core(2:len_trim(core))
   end do
   j = index(core, ".*")
   if (j > 0) then
      prefix = core(1:j-1)
      suffix = core(j+2:len_trim(core))
      out = (len_trim(prefix) == 0 .or. index(name, trim(prefix)) == 1) .and. &
         & (len_trim(suffix) == 0 .or. (len_trim(name) >= len_trim(suffix) .and. &
         & name(len_trim(name)-len_trim(suffix)+1:len_trim(name)) == suffix))
   else
      suffix = core
      if (len_trim(suffix) >= 2 .and. suffix(1:2) == ".*") suffix = suffix(3:len_trim(suffix))
      out = len_trim(name) >= len_trim(suffix) .and. name(len_trim(name)-len_trim(suffix)+1:len_trim(name)) == suffix
   end if
else if (n > 2 .and. pat(1:2) == ".*") then
   out = index(name, pat(3:n)) > 0
else if (n == 2 .and. pat(1:2) == ".*") then
   out = .true.
else if (n > 1 .and. pat(1:1) == "*") then
   out = index(name, pat(2:n)) > 0
else if (pat(1:1) == "*") then
   out = .true.
else if (pat(n:n) == "*") then
   out = index(name, pat(1:n-1)) == 1
else
   out = index(name, pat) > 0
end if
end function list_files_pattern_match

pure function strsplit_fixed(s, delim) result(out)
! Implement R-like character helper strsplit_fixed.
character(len=*), intent(in) :: s ! string to split
character(len=*), intent(in) :: delim ! fixed delimiter string
character(len=:), allocatable :: out(:)
integer :: i, start, pos, n, dlen, maxlen
! R's strsplit(s, "") splits into individual characters.
if (len(delim) == 0) then
   allocate(character(len=1) :: out(len(s)))
   do i = 1, len(s)
      out(i) = s(i:i)
   end do
   return
end if
dlen = max(1, len(delim))
n = 1
start = 1
do
   pos = index(s(start:), delim)
   if (pos <= 0) exit
   n = n + 1
   start = start + pos + dlen - 1
   if (start > len(s) + 1) exit
end do
maxlen = max(1, len(s))
allocate(character(len=maxlen) :: out(n))
start = 1
do i = 1, n
   pos = index(s(start:), delim)
   if (pos <= 0) then
      out(i) = s(start:)
   else
      out(i) = s(start:start + pos - 2)
      start = start + pos + dlen - 1
   end if
end do
end function strsplit_fixed

pure function toupper(s) result(out)
! Implement R-like character helper toupper.
character(len=*), intent(in) :: s ! string to convert
character(len=len(s)) :: out
integer :: i, c
out = s
do i = 1, len(s)
   c = iachar(out(i:i))
   if (c >= iachar("a") .and. c <= iachar("z")) out(i:i) = achar(c - 32)
end do
end function toupper

pure function tolower(s) result(out)
! Implement R-like character helper tolower.
character(len=*), intent(in) :: s ! string to convert
character(len=len(s)) :: out
integer :: i, c
out = s
do i = 1, len(s)
   c = iachar(out(i:i))
   if (c >= iachar("A") .and. c <= iachar("Z")) out(i:i) = achar(c + 32)
end do
end function tolower

pure function casefold(s, upper) result(out)
! Implement R-like character helper casefold.
! If upper is absent or false, convert to lower case; if true, convert to upper case.
character(len=*), intent(in) :: s ! string to convert
logical, intent(in), optional :: upper ! choose upper-case conversion flag
character(len=len(s)) :: out
logical :: up
up = .false.
if (present(upper)) up = upper
if (up) then
   out = toupper(s)
else
   out = tolower(s)
end if
end function casefold

pure function trimws(s, which) result(out)
! Implement R-like character helper trimws.
! Valid which values: "both" (default), "left", "right".
character(len=*), intent(in) :: s ! string to trim
character(len=*), intent(in), optional :: which ! side to trim
character(len=:), allocatable :: out
character(len=16) :: w
w = "both"
if (present(which)) w = which
select case (trim(w))
case ("left")
   out = adjustl(s)
case ("right")
   out = trim(s)
case default
   out = trim(adjustl(s))
end select
end function trimws

pure function grep_value_char(pattern, x) result(out)
! Return character entries containing a fixed substring.
character(len=*), intent(in) :: pattern
character(len=*), intent(in) :: x(:)
character(len=:), allocatable :: out(:)
logical, allocatable :: keep(:)
integer :: i
allocate(keep(size(x)))
do i = 1, size(x)
   keep(i) = index(x(i), pattern) > 0
end do
out = pack(x, keep)
end function grep_value_char

function r_command_args(trailing_only, file_arg) result(out)
! Return command arguments, optionally prepending a generated --file= entry.
logical, intent(in) :: trailing_only
character(len=*), intent(in), optional :: file_arg
character(len=:), allocatable :: out(:)
integer :: i, n, stat, arg_len, out_len, prefix
character(len=4096) :: buf
n = command_argument_count()
prefix = 0
if (.not. trailing_only .and. present(file_arg)) prefix = 1
out_len = 1
if (prefix == 1) out_len = max(out_len, len(file_arg))
do i = 1, n
   call get_command_argument(i, length=arg_len, status=stat)
   if (stat == 0) out_len = max(out_len, arg_len)
end do
allocate(character(len=out_len) :: out(n + prefix))
if (prefix == 1) out(1) = file_arg
do i = 1, n
   call get_command_argument(i, buf, status=stat)
   if (stat == 0) then
      out(i + prefix) = trim(buf)
   else
      out(i + prefix) = ""
   end if
end do
end function r_command_args

pure function replace_first_fixed(s, old, new) result(out)
! Return a copy with selected first fixed entries replaced.
character(len=*), intent(in) :: s ! source string
character(len=*), intent(in) :: old ! fixed substring to replace
character(len=*), intent(in) :: new ! replacement text
character(len=:), allocatable :: out
integer :: pos
pos = index(s, old)
if (pos <= 0) then
   out = s
else
   out = s(1:pos - 1) // new // s(pos + len(old):)
end if
end function replace_first_fixed

pure function replace_all_fixed(s, old, new) result(out)
! Return a copy with selected all fixed entries replaced.
character(len=*), intent(in) :: s ! source string
character(len=*), intent(in) :: old ! fixed substring to replace
character(len=*), intent(in) :: new ! replacement text
character(len=:), allocatable :: out
character(len=:), allocatable :: rest
integer :: pos
if (len(old) == 0) then
   out = s
   return
end if
out = ""
rest = s
do
   pos = index(rest, old)
   if (pos <= 0) exit
   out = out // rest(1:pos - 1) // new
   rest = rest(pos + len(old):)
end do
out = out // rest
end function replace_all_fixed

pure function chartr(old, new, s) result(out)
! Implement R-like character helper chartr.
character(len=*), intent(in) :: old ! characters to translate from
character(len=*), intent(in) :: new ! replacement characters
character(len=*), intent(in) :: s ! source string
character(len=len(s)) :: out
integer :: i, p
out = s
do i = 1, len(s)
   p = index(old, s(i:i))
   if (p > 0 .and. p <= len(new)) out(i:i) = new(p:p)
end do
end function chartr

pure function fft(x) result(out)
! Discrete Fourier transform (forward), matching R's fft() for numeric input.
real(kind=dp), intent(in) :: x(:)
complex(kind=dp), allocatable :: out(:)
integer :: nfft, k, j
real(kind=dp) :: ang, twopi
nfft = size(x)
allocate(out(nfft))
twopi = 2.0_dp * acos(-1.0_dp)
do k = 1, nfft
   out(k) = (0.0_dp, 0.0_dp)
   do j = 1, nfft
      ang = -twopi * real((k - 1) * (j - 1), dp) / real(nfft, dp)
      out(k) = out(k) + x(j) * cmplx(cos(ang), sin(ang), kind=dp)
   end do
end do
end function fft

pure function kronecker(a, b) result(out)
! Kronecker product of two matrices (R's %x% / kronecker()).
real(kind=dp), intent(in) :: a(:,:), b(:,:)
real(kind=dp), allocatable :: out(:,:)
integer :: ma, na, mb, nb, i, j, k, l
ma = size(a, 1); na = size(a, 2)
mb = size(b, 1); nb = size(b, 2)
allocate(out(ma * mb, na * nb))
do j = 1, na
   do i = 1, ma
      do l = 1, nb
         do k = 1, mb
            out((i - 1) * mb + k, (j - 1) * nb + l) = a(i, j) * b(k, l)
         end do
      end do
   end do
end do
end function kronecker

pure function nextn(n) result(out)
! Smallest integer >= n that is a product of factors 2, 3 and 5 (R's nextn default).
integer, intent(in) :: n
integer :: out, m
out = max(n, 1)
do
   m = out
   do while (mod(m, 2) == 0)
      m = m / 2
   end do
   do while (mod(m, 3) == 0)
      m = m / 3
   end do
   do while (mod(m, 5) == 0)
      m = m / 5
   end do
   if (m == 1) exit
   out = out + 1
end do
end function nextn

pure function urldecode(s) result(out)
! Decode percent-encoded (%XX) escapes as R's utils::URLdecode does.
character(len=*), intent(in) :: s
character(len=:), allocatable :: out
character(len=len(s)) :: buf
integer :: i, n, code, ios
n = 0
i = 1
do while (i <= len(s))
   if (s(i:i) == "%" .and. i + 2 <= len(s)) then
      read(s(i+1:i+2), "(z2)", iostat=ios) code
      if (ios == 0) then
         n = n + 1
         buf(n:n) = achar(code)
         i = i + 3
         cycle
      end if
   end if
   n = n + 1
   buf(n:n) = s(i:i)
   i = i + 1
end do
out = buf(1:n)
end function urldecode

pure function lower_tri(x, diag) result(out)
! Return R-like matrix index helper lower_tri.
! If diag is absent or false, exclude diagonal positions.
real(kind=dp), intent(in) :: x(:,:) ! matrix whose shape defines the mask
logical, intent(in), optional :: diag ! include diagonal positions flag
logical, allocatable :: out(:,:)
integer :: i, j
logical :: include_diag
include_diag = .false.
if (present(diag)) include_diag = diag
allocate(out(size(x,1), size(x,2)))
do j = 1, size(x,2)
   do i = 1, size(x,1)
      if (include_diag) then
         out(i,j) = i >= j
      else
         out(i,j) = i > j
      end if
   end do
end do
end function lower_tri

pure function upper_tri(x, diag) result(out)
! Return R-like matrix index helper upper_tri.
! If diag is absent or false, exclude diagonal positions.
real(kind=dp), intent(in) :: x(:,:) ! matrix whose shape defines the mask
logical, intent(in), optional :: diag ! include diagonal positions flag
logical, allocatable :: out(:,:)
integer :: i, j
logical :: include_diag
include_diag = .false.
if (present(diag)) include_diag = diag
allocate(out(size(x,1), size(x,2)))
do j = 1, size(x,2)
   do i = 1, size(x,1)
      if (include_diag) then
         out(i,j) = i <= j
      else
         out(i,j) = i < j
      end if
   end do
end do
end function upper_tri

pure function row_index_mat(x) result(out)
! Return R-like matrix index helper row_index_mat.
real(kind=dp), intent(in) :: x(:,:) ! matrix whose shape defines row indices
integer, allocatable :: out(:,:)
integer :: i
allocate(out(size(x,1), size(x,2)))
do i = 1, size(x,1)
   out(i, :) = i
end do
end function row_index_mat

pure function col_index_mat(x) result(out)
! Return R-like matrix index helper col_index_mat.
real(kind=dp), intent(in) :: x(:,:) ! matrix whose shape defines column indices
integer, allocatable :: out(:,:)
integer :: j
allocate(out(size(x,1), size(x,2)))
do j = 1, size(x,2)
   out(:, j) = j
end do
end function col_index_mat

pure subroutine matrix_set_grow_real(x, row, col, value)
! Assign a matrix element, growing matrix-backed R objects as needed.
real(kind=dp), allocatable, intent(inout) :: x(:,:)
integer, intent(in) :: row, col
real(kind=dp), intent(in) :: value
real(kind=dp), allocatable :: grown(:,:)
integer :: nr, nc
nr = max(row, size(x, 1))
nc = max(col, size(x, 2))
if (nr /= size(x, 1) .or. nc /= size(x, 2)) then
   allocate(grown(nr, nc))
   grown = r_na_real()
   grown(1:size(x, 1), 1:size(x, 2)) = x
   call move_alloc(grown, x)
end if
x(row, col) = value
end subroutine matrix_set_grow_real

pure elemental logical function is_na_real_scalar(x) result(out)
! True when real scalar is NA/NaN in this subset.
real(kind=dp), intent(in) :: x ! value to test
out = (x /= x)
end function is_na_real_scalar

pure function is_na_real_vec(x) result(out)
! Elementwise NA test for a real vector.
real(kind=dp), intent(in) :: x(:) ! values to test
logical, allocatable :: out(:)
allocate(out(size(x)))
out = (x /= x)
end function is_na_real_vec

pure function is_na_real_mat(x) result(out)
! Elementwise NA test for a real matrix.
real(kind=dp), intent(in) :: x(:,:)
logical, allocatable :: out(:,:)
allocate(out(size(x, 1), size(x, 2)))
out = (x /= x)
end function is_na_real_mat

pure elemental logical function is_na_int_scalar(x) result(out)
! True when integer scalar uses NA sentinel.
integer, intent(in) :: x ! value to test
out = (x == -huge(0))
end function is_na_int_scalar

pure function is_na_int_vec(x) result(out)
! Elementwise NA test for integer vector.
integer, intent(in) :: x(:) ! values to test
logical, allocatable :: out(:)
allocate(out(size(x)))
out = (x == -huge(0))
end function is_na_int_vec

pure function is_na_int_mat(x) result(out)
! Elementwise NA test for an integer matrix.
integer, intent(in) :: x(:,:)
logical, allocatable :: out(:,:)
allocate(out(size(x, 1), size(x, 2)))
out = (x == -huge(0))
end function is_na_int_mat

pure elemental logical function is_na_logical_scalar(x) result(out)
! Logical values have no NA sentinel in this subset.
logical, intent(in) :: x ! value to test
out = .false.
end function is_na_logical_scalar

pure function is_na_logical_vec(x) result(out)
! Elementwise NA test for a logical vector.
logical, intent(in) :: x(:) ! values to test
logical, allocatable :: out(:)
allocate(out(size(x)))
out = .false.
end function is_na_logical_vec

pure function is_na_logical_mat(x) result(out)
! Elementwise NA test for a logical matrix.
logical, intent(in) :: x(:,:)
logical, allocatable :: out(:,:)
allocate(out(size(x, 1), size(x, 2)))
out = .false.
end function is_na_logical_mat

pure elemental logical function is_na_complex_scalar(x) result(out)
! True when either complex component is NA/NaN.
complex(kind=dp), intent(in) :: x ! value to test
out = (real(x, kind=dp) /= real(x, kind=dp)) .or. (aimag(x) /= aimag(x))
end function is_na_complex_scalar

pure function is_na_complex_vec(x) result(out)
! Elementwise NA test for a complex vector.
complex(kind=dp), intent(in) :: x(:) ! values to test
logical, allocatable :: out(:)
allocate(out(size(x)))
out = (real(x, kind=dp) /= real(x, kind=dp)) .or. (aimag(x) /= aimag(x))
end function is_na_complex_vec

pure function is_na_complex_mat(x) result(out)
! Elementwise NA test for a complex matrix.
complex(kind=dp), intent(in) :: x(:,:)
logical, allocatable :: out(:,:)
allocate(out(size(x, 1), size(x, 2)))
out = (real(x, kind=dp) /= real(x, kind=dp)) .or. (aimag(x) /= aimag(x))
end function is_na_complex_mat

pure elemental logical function is_na_char_scalar(x) result(out)
! True when character scalar uses NA sentinel in this subset.
character(len=*), intent(in) :: x ! string to test
out = (x == "")
end function is_na_char_scalar

pure function is_na_char_vec(x) result(out)
! Elementwise NA test for character vector.
character(len=*), intent(in) :: x(:) ! strings to test
logical, allocatable :: out(:)
integer :: i
allocate(out(size(x)))
do i = 1, size(x)
   out(i) = (x(i) == "")
end do
end function is_na_char_vec

pure function is_na_char_mat(x) result(out)
! Elementwise NA test for a character matrix.
character(len=*), intent(in) :: x(:,:)
logical, allocatable :: out(:,:)
allocate(out(size(x, 1), size(x, 2)))
out = (x == "")
end function is_na_char_mat

pure function which_logical(x) result(out)
! Runtime helper for R-compatible which logical.
logical, intent(in) :: x(:) ! selection mask
integer, allocatable :: out(:)
integer :: i, n
n = count(x)
allocate(out(n))
n = 0
do i = 1, size(x)
   if (x(i)) then
      n = n + 1
      out(n) = i
   end if
end do
end function which_logical

pure function which_logical_mat(x) result(out)
! Runtime helper for R-compatible which on matrix masks.
logical, intent(in) :: x(:,:) ! selection mask matrix
integer, allocatable :: out(:)
integer :: i, j, k, n
n = count(x)
allocate(out(n))
k = 0
do j = 1, size(x, 2)
   do i = 1, size(x, 1)
      if (x(i, j)) then
         k = k + 1
         out(k) = i + (j - 1) * size(x, 1)
      end if
   end do
end do
end function which_logical_mat

pure function which_first(x) result(out)
! Return the first R-style which() index, or 0 when no element is true.
logical, intent(in) :: x(:) ! selection mask
integer :: out
integer :: i
out = 0
do i = 1, size(x)
   if (x(i)) then
      out = i
      return
   end if
end do
end function which_first

pure function which_last(x) result(out)
! Return the last R-style which() index, or 0 when no element is true.
logical, intent(in) :: x(:) ! selection mask
integer :: out
integer :: i
out = 0
do i = size(x), 1, -1
   if (x(i)) then
      out = i
      return
   end if
end do
end function which_last

pure function which_arr_ind(x) result(out)
! Return R-like matrix index helper which_arr_ind.
logical, intent(in) :: x(:,:) ! selection mask matrix
integer, allocatable :: out(:,:)
integer :: i, j, k, n
n = count(x)
allocate(out(n, 2))
k = 0
do j = 1, size(x, 2)
   do i = 1, size(x, 1)
      if (x(i, j)) then
         k = k + 1
         out(k, 1) = i
         out(k, 2) = j
      end if
   end do
end do
end function which_arr_ind

pure function replace_real_idx_scalar(x, idx, values) result(out)
! Return a copy with selected real idx scalar entries replaced.
real(kind=dp), intent(in) :: x(:) ! source vector
real(kind=dp), intent(in) :: values ! replacement value
integer, intent(in) :: idx(:) ! one-based replacement positions
real(kind=dp), allocatable :: out(:)
integer :: i
out = x
do i = 1, size(idx)
   if (idx(i) >= 1 .and. idx(i) <= size(out)) out(idx(i)) = values
end do
end function replace_real_idx_scalar

pure function replace_real_idx_vec(x, idx, values) result(out)
! Return a copy with selected real idx vec entries replaced.
real(kind=dp), intent(in) :: x(:) ! source vector
real(kind=dp), intent(in) :: values(:) ! replacement values, recycled over positions
integer, intent(in) :: idx(:) ! one-based replacement positions
real(kind=dp), allocatable :: out(:)
integer :: i
out = x
if (size(values) <= 0) return
do i = 1, size(idx)
   if (idx(i) >= 1 .and. idx(i) <= size(out)) out(idx(i)) = values(1 + mod(i - 1, size(values)))
end do
end function replace_real_idx_vec

pure function replace_real_mask_scalar(x, mask, values) result(out)
! Return a copy with selected real mask scalar entries replaced.
real(kind=dp), intent(in) :: x(:) ! source vector
real(kind=dp), intent(in) :: values ! replacement value
logical, intent(in) :: mask(:) ! replacement mask
real(kind=dp), allocatable :: out(:)
out = x
where (mask(1:min(size(mask), size(out)))) out(1:min(size(mask), size(out))) = values
end function replace_real_mask_scalar

pure function replace_real_mask_vec(x, mask, values) result(out)
! Return a copy with selected real mask vec entries replaced.
real(kind=dp), intent(in) :: x(:) ! source vector
real(kind=dp), intent(in) :: values(:) ! replacement values, recycled over selected entries
logical, intent(in) :: mask(:) ! replacement mask
real(kind=dp), allocatable :: out(:)
integer :: i, j
out = x
if (size(values) <= 0) return
j = 0
do i = 1, min(size(mask), size(out))
   if (mask(i)) then
      j = j + 1
      out(i) = values(1 + mod(j - 1, size(values)))
   end if
end do
end function replace_real_mask_vec

pure function replace_int_idx_scalar(x, idx, values) result(out)
! Return a copy with selected int idx scalar entries replaced.
integer, intent(in) :: x(:) ! source vector
integer, intent(in) :: values ! replacement value
integer, intent(in) :: idx(:) ! one-based replacement positions
integer, allocatable :: out(:)
integer :: i
out = x
do i = 1, size(idx)
   if (idx(i) >= 1 .and. idx(i) <= size(out)) out(idx(i)) = values
end do
end function replace_int_idx_scalar

pure function replace_int_idx_vec(x, idx, values) result(out)
! Return a copy with selected int idx vec entries replaced.
integer, intent(in) :: x(:) ! source vector
integer, intent(in) :: values(:) ! replacement values, recycled over positions
integer, intent(in) :: idx(:) ! one-based replacement positions
integer, allocatable :: out(:)
integer :: i
out = x
if (size(values) <= 0) return
do i = 1, size(idx)
   if (idx(i) >= 1 .and. idx(i) <= size(out)) out(idx(i)) = values(1 + mod(i - 1, size(values)))
end do
end function replace_int_idx_vec

pure function replace_int_idx_scalar_real(x, idx, values) result(out)
! Return a copy with selected int idx scalar real entries replaced.
integer, intent(in) :: x(:) ! source vector
integer, intent(in) :: idx(:) ! one-based replacement positions
real(kind=dp), intent(in) :: values ! replacement value
real(kind=dp), allocatable :: out(:)
integer :: i
out = real(x, kind=dp)
do i = 1, size(idx)
   if (idx(i) >= 1 .and. idx(i) <= size(out)) out(idx(i)) = values
end do
end function replace_int_idx_scalar_real

pure function replace_int_idx_vec_real(x, idx, values) result(out)
! Return a copy with selected int idx vec real entries replaced.
integer, intent(in) :: x(:) ! source vector
integer, intent(in) :: idx(:) ! one-based replacement positions
real(kind=dp), intent(in) :: values(:) ! replacement values, recycled over positions
real(kind=dp), allocatable :: out(:)
integer :: i
out = real(x, kind=dp)
if (size(values) <= 0) return
do i = 1, size(idx)
   if (idx(i) >= 1 .and. idx(i) <= size(out)) out(idx(i)) = values(1 + mod(i - 1, size(values)))
end do
end function replace_int_idx_vec_real

pure function replace_int_mask_scalar(x, mask, values) result(out)
! Return a copy with selected int mask scalar entries replaced.
integer, intent(in) :: x(:) ! source vector
integer, intent(in) :: values ! replacement value
logical, intent(in) :: mask(:) ! replacement mask
integer, allocatable :: out(:)
out = x
where (mask(1:min(size(mask), size(out)))) out(1:min(size(mask), size(out))) = values
end function replace_int_mask_scalar

pure function replace_int_mask_vec(x, mask, values) result(out)
! Return a copy with selected int mask vec entries replaced.
integer, intent(in) :: x(:) ! source vector
integer, intent(in) :: values(:) ! replacement values, recycled over selected entries
logical, intent(in) :: mask(:) ! replacement mask
integer, allocatable :: out(:)
integer :: i, j
out = x
if (size(values) <= 0) return
j = 0
do i = 1, min(size(mask), size(out))
   if (mask(i)) then
      j = j + 1
      out(i) = values(1 + mod(j - 1, size(values)))
   end if
end do
end function replace_int_mask_vec

pure function replace_int_mask_scalar_real(x, mask, values) result(out)
! Return a copy with selected int mask scalar real entries replaced.
integer, intent(in) :: x(:) ! source vector
logical, intent(in) :: mask(:) ! replacement mask
real(kind=dp), intent(in) :: values ! replacement value
real(kind=dp), allocatable :: out(:)
out = real(x, kind=dp)
where (mask(1:min(size(mask), size(out)))) out(1:min(size(mask), size(out))) = values
end function replace_int_mask_scalar_real

pure function replace_int_mask_vec_real(x, mask, values) result(out)
! Return a copy with selected int mask vec real entries replaced.
integer, intent(in) :: x(:) ! source vector
logical, intent(in) :: mask(:) ! replacement mask
real(kind=dp), intent(in) :: values(:) ! replacement values, recycled over selected entries
real(kind=dp), allocatable :: out(:)
integer :: i, j
out = real(x, kind=dp)
if (size(values) <= 0) return
j = 0
do i = 1, min(size(mask), size(out))
   if (mask(i)) then
      j = j + 1
      out(i) = values(1 + mod(j - 1, size(values)))
   end if
end do
end function replace_int_mask_vec_real

pure function replace_char_mask_scalar(x, mask, values) result(out)
! Return a copy with selected char mask scalar entries replaced.
character(len=*), intent(in) :: x(:) ! source strings
logical, intent(in) :: mask(:) ! replacement mask
character(len=*), intent(in) :: values ! replacement string
character(len=:), allocatable :: out(:)
integer :: i, lch
lch = max(1, len(x), len(values))
allocate(character(len=lch) :: out(size(x)))
do i = 1, size(x)
   out(i) = x(i)
end do
do i = 1, min(size(mask), size(out))
   if (mask(i)) out(i) = values
end do
end function replace_char_mask_scalar

pure function replace_char_mask_vec(x, mask, values, value_len) result(out)
! Return a copy with selected char mask vec entries replaced.
character(len=*), intent(in) :: x(:) ! source strings
logical, intent(in) :: mask(:) ! replacement mask
character(len=*), intent(in) :: values(:) ! replacement strings, recycled over selected entries
integer, intent(in), optional :: value_len ! minimum output string length
character(len=:), allocatable :: out(:)
integer :: i, j, lch
lch = max(1, len(x), len(values))
if (present(value_len)) lch = max(lch, value_len)
allocate(character(len=lch) :: out(size(x)))
do i = 1, size(x)
   out(i) = x(i)
end do
if (size(values) <= 0) return
j = 0
do i = 1, min(size(mask), size(out))
   if (mask(i)) then
      j = j + 1
      out(i) = values(1 + mod(j - 1, size(values)))
   end if
end do
end function replace_char_mask_vec

pure function replace_logical_mask_scalar(x, mask, values) result(out)
! Return a copy with selected logical mask scalar entries replaced.
logical, intent(in) :: x(:) ! source values
logical, intent(in) :: mask(:) ! replacement mask
logical, intent(in) :: values ! replacement value
logical, allocatable :: out(:)
out = x
where (mask(1:min(size(mask), size(out)))) out(1:min(size(mask), size(out))) = values
end function replace_logical_mask_scalar

pure function rle_real(x) result(out)
! Compute run-length encoding for real input.
real(kind=dp), intent(in) :: x(:) ! values to encode
type(rle_real_t) :: out
integer :: i, j, nr
nr = 0
if (size(x) > 0) then
   nr = 1
   do i = 2, size(x)
      if (x(i) /= x(i - 1)) nr = nr + 1
   end do
end if
allocate(out%lengths(nr), out%values(nr))
if (nr == 0) return
j = 1
out%lengths(1) = 1
out%values(1) = x(1)
do i = 2, size(x)
   if (x(i) == x(i - 1)) then
      out%lengths(j) = out%lengths(j) + 1
   else
      j = j + 1
      out%lengths(j) = 1
      out%values(j) = x(i)
   end if
end do
end function rle_real

pure elemental logical function rle_int_values_equal(a, b) result(out)
! Integer NA never continues an rle() run, including another integer NA.
integer, intent(in) :: a, b
out = a == b .and. a /= -huge(0)
end function rle_int_values_equal

pure function rle_int(x) result(out)
! Compute run-length encoding for int input.
integer, intent(in) :: x(:) ! values to encode
type(rle_int_t) :: out
integer :: i, j, nr
nr = 0
if (size(x) > 0) then
   nr = 1
   do i = 2, size(x)
      if (.not. rle_int_values_equal(x(i), x(i - 1))) nr = nr + 1
   end do
end if
allocate(out%lengths(nr), out%values(nr))
if (nr == 0) return
j = 1
out%lengths(1) = 1
out%values(1) = x(1)
do i = 2, size(x)
   if (rle_int_values_equal(x(i), x(i - 1))) then
      out%lengths(j) = out%lengths(j) + 1
   else
      j = j + 1
      out%lengths(j) = 1
      out%values(j) = x(i)
   end if
end do
end function rle_int

pure function rle_char(x) result(out)
! Compute run-length encoding for char input.
character(len=*), intent(in) :: x(:) ! strings to encode
type(rle_char_t) :: out
integer :: i, j, nr
nr = 0
if (size(x) > 0) then
   nr = 1
   do i = 2, size(x)
      if (x(i) /= x(i - 1)) nr = nr + 1
   end do
end if
allocate(out%lengths(nr))
allocate(character(len=max(1, len(x))) :: out%values(nr))
if (nr == 0) return
j = 1
out%lengths(1) = 1
out%values(1) = x(1)
do i = 2, size(x)
   if (x(i) == x(i - 1)) then
      out%lengths(j) = out%lengths(j) + 1
   else
      j = j + 1
      out%lengths(j) = 1
      out%values(j) = x(i)
   end if
end do
end function rle_char

pure function rle_logical(x) result(out)
! Compute run-length encoding for logical input.
logical, intent(in) :: x(:) ! values to encode
type(rle_logical_t) :: out
integer :: i, j, nr
nr = 0
if (size(x) > 0) then
   nr = 1
   do i = 2, size(x)
      if (x(i) .neqv. x(i - 1)) nr = nr + 1
   end do
end if
allocate(out%lengths(nr), out%values(nr))
if (nr == 0) return
j = 1
out%lengths(1) = 1
out%values(1) = x(1)
do i = 2, size(x)
   if (x(i) .eqv. x(i - 1)) then
      out%lengths(j) = out%lengths(j) + 1
   else
      j = j + 1
      out%lengths(j) = 1
      out%values(j) = x(i)
   end if
end do
end function rle_logical

pure function inverse_rle_real(fit) result(out)
! Expand run-length encoded real values.
type(rle_real_t), intent(in) :: fit ! run-length encoding to expand
real(kind=dp), allocatable :: out(:)
integer :: i, p, n
n = sum(fit%lengths)
allocate(out(n))
p = 1
do i = 1, size(fit%lengths)
   if (fit%lengths(i) > 0) out(p:p + fit%lengths(i) - 1) = fit%values(i)
   p = p + fit%lengths(i)
end do
end function inverse_rle_real

pure function inverse_rle_int(fit) result(out)
! Expand run-length encoded int values.
type(rle_int_t), intent(in) :: fit ! run-length encoding to expand
integer, allocatable :: out(:)
integer :: i, p, n
n = sum(fit%lengths)
allocate(out(n))
p = 1
do i = 1, size(fit%lengths)
   if (fit%lengths(i) > 0) out(p:p + fit%lengths(i) - 1) = fit%values(i)
   p = p + fit%lengths(i)
end do
end function inverse_rle_int

pure function inverse_rle_char(fit) result(out)
! Expand run-length encoded char values.
type(rle_char_t), intent(in) :: fit ! run-length encoding to expand
character(len=:), allocatable :: out(:)
integer :: i, p, n
n = sum(fit%lengths)
allocate(character(len=max(1, len(fit%values))) :: out(n))
p = 1
do i = 1, size(fit%lengths)
   if (fit%lengths(i) > 0) out(p:p + fit%lengths(i) - 1) = fit%values(i)
   p = p + fit%lengths(i)
end do
end function inverse_rle_char

pure function inverse_rle_logical(fit) result(out)
! Expand run-length encoded logical values.
type(rle_logical_t), intent(in) :: fit ! run-length encoding to expand
logical, allocatable :: out(:)
integer :: i, p, n
n = sum(fit%lengths)
allocate(out(n))
p = 1
do i = 1, size(fit%lengths)
   if (fit%lengths(i) > 0) out(p:p + fit%lengths(i) - 1) = fit%values(i)
   p = p + fit%lengths(i)
end do
end function inverse_rle_logical

subroutine print_rle_real(fit)
! Print rle real values in an R-like format.
type(rle_real_t), intent(in) :: fit ! run-length encoding to print
write(*,'(a)') "Run Length Encoding"
write(*,'(a)', advance='no') "  lengths: "
write(*,"(*(1x,i0))") fit%lengths
write(*,'(a)', advance='no') "  values : "
call print_real_vector(fit%values)
end subroutine print_rle_real

subroutine print_rle_int(fit)
! Print rle int values in an R-like format.
type(rle_int_t), intent(in) :: fit ! run-length encoding to print
write(*,'(a)') "Run Length Encoding"
write(*,'(a)', advance='no') "  lengths: "
write(*,"(*(1x,i0))") fit%lengths
write(*,'(a)', advance='no') "  values : "
write(*,"(*(1x,i0))") fit%values
end subroutine print_rle_int

subroutine print_rle_char(fit)
! Print rle char values in an R-like format.
type(rle_char_t), intent(in) :: fit ! run-length encoding to print
write(*,'(a)') "Run Length Encoding"
write(*,'(a)', advance='no') "  lengths: "
write(*,"(*(1x,i0))") fit%lengths
write(*,'(a)', advance='no') "  values : "
call print_char_vector(fit%values)
end subroutine print_rle_char

subroutine print_rle_logical(fit)
! Print rle logical values in an R-like format.
type(rle_logical_t), intent(in) :: fit ! run-length encoding to print
write(*,'(a)') "Run Length Encoding"
write(*,'(a)', advance='no') "  lengths: "
write(*,"(*(1x,i0))") fit%lengths
write(*,'(a)', advance='no') "  values : "
write(*,"(*(g0,1x))") fit%values
end subroutine print_rle_logical

pure function r_typeof_real_scalar(x) result(out)
! Return R-like type label for real scalar.
real(kind=dp), intent(in) :: x
character(len=:), allocatable :: out
out = "double"
if (storage_size(x) < 0) out = ""
end function r_typeof_real_scalar

pure function r_typeof_real_vec(x) result(out)
! Return R-like type label for real vector.
real(kind=dp), intent(in) :: x(:)
character(len=:), allocatable :: out
out = "double"
if (size(x) < 0) out = ""
end function r_typeof_real_vec

pure function r_typeof_real_mat(x) result(out)
! Return R-like type label for real matrix.
real(kind=dp), intent(in) :: x(:,:)
character(len=:), allocatable :: out
out = "double"
if (size(x) < 0) out = ""
end function r_typeof_real_mat

pure function r_typeof_complex_scalar(x) result(out)
! Return R-like type label for complex scalar.
complex(kind=dp), intent(in) :: x
character(len=:), allocatable :: out
out = "complex"
if (storage_size(x) < 0) out = ""
end function r_typeof_complex_scalar

pure function r_typeof_complex_vec(x) result(out)
! Return R-like type label for complex vector.
complex(kind=dp), intent(in) :: x(:)
character(len=:), allocatable :: out
out = "complex"
if (size(x) < 0) out = ""
end function r_typeof_complex_vec

pure function r_typeof_complex_mat(x) result(out)
! Return R-like type label for complex matrix.
complex(kind=dp), intent(in) :: x(:,:)
character(len=:), allocatable :: out
out = "complex"
if (size(x) < 0) out = ""
end function r_typeof_complex_mat

pure function r_typeof_int_scalar(x) result(out)
! Return R-like type label for integer scalar.
integer, intent(in) :: x
character(len=:), allocatable :: out
out = "integer"
if (storage_size(x) < 0) out = ""
end function r_typeof_int_scalar

pure function r_typeof_int_vec(x) result(out)
! Return R-like type label for integer vector.
integer, intent(in) :: x(:)
character(len=:), allocatable :: out
out = "integer"
if (size(x) < 0) out = ""
end function r_typeof_int_vec

pure function r_typeof_int_mat(x) result(out)
! Return R-like type label for integer matrix.
integer, intent(in) :: x(:,:)
character(len=:), allocatable :: out
out = "integer"
if (size(x) < 0) out = ""
end function r_typeof_int_mat

pure function r_typeof_char_scalar(x) result(out)
! Return R-like type label for character scalar.
character(len=*), intent(in) :: x
character(len=:), allocatable :: out
out = "character"
if (len(x) < 0) out = ""
end function r_typeof_char_scalar

pure function r_typeof_char_vec(x) result(out)
! Return R-like type label for character vector.
character(len=*), intent(in) :: x(:)
character(len=:), allocatable :: out
out = "character"
if (size(x) < 0) out = ""
end function r_typeof_char_vec

pure function r_typeof_char_mat(x) result(out)
! Return R-like type label for character matrix.
character(len=*), intent(in) :: x(:,:)
character(len=:), allocatable :: out
out = "character"
if (size(x) < 0) out = ""
end function r_typeof_char_mat

pure function r_typeof_logical_scalar(x) result(out)
! Return R-like type label for logical scalar.
logical, intent(in) :: x
character(len=:), allocatable :: out
out = "logical"
if (merge(0, 0, x) < 0) out = ""
end function r_typeof_logical_scalar

pure function r_typeof_logical_vec(x) result(out)
! Return R-like type label for logical vector.
logical, intent(in) :: x(:)
character(len=:), allocatable :: out
out = "logical"
if (size(x) < 0) out = ""
end function r_typeof_logical_vec

pure function r_typeof_logical_mat(x) result(out)
! Return R-like type label for logical matrix.
logical, intent(in) :: x(:,:)
character(len=:), allocatable :: out
out = "logical"
if (size(x) < 0) out = ""
end function r_typeof_logical_mat

pure function hist_nbreaks(x, breaks, plot) result(out)
! Compute a small R-like histogram object for numeric vector input.
real(kind=dp), intent(in) :: x(:)
integer, intent(in), optional :: breaks
logical, intent(in), optional :: plot
type(hist_result_t) :: out
integer :: nb, i, j
real(kind=dp) :: xmin, xmax, width, total, raw_width, base_width, step_width
logical :: ignored_plot
ignored_plot = .false.
if (present(plot)) ignored_plot = plot
nb = 10
if (present(breaks)) nb = max(1, breaks)
allocate(out%breaks(nb + 1), out%mids(nb), out%density(nb), out%counts(nb))
out%counts = 0
out%density = 0.0_dp
if (size(x) <= 0) then
   out%breaks = 0.0_dp
   out%mids = 0.0_dp
   return
end if
xmin = minval(x)
xmax = maxval(x)
if (xmax <= xmin) xmax = xmin + 1.0_dp
raw_width = (xmax - xmin) / real(nb, kind=dp)
base_width = 10.0_dp ** floor(log10(raw_width))
step_width = base_width
if (step_width < raw_width) step_width = 2.0_dp * base_width
if (step_width < raw_width) step_width = 5.0_dp * base_width
if (step_width < raw_width) step_width = 10.0_dp * base_width
xmin = floor(xmin / step_width) * step_width
xmax = ceiling(xmax / step_width) * step_width
nb = max(1, nint((xmax - xmin) / step_width))
deallocate(out%breaks, out%mids, out%density, out%counts)
allocate(out%breaks(nb + 1), out%mids(nb), out%density(nb), out%counts(nb))
out%counts = 0
out%density = 0.0_dp
width = step_width
do i = 1, nb + 1
   out%breaks(i) = xmin + real(i - 1, kind=dp) * width
end do
out%breaks(nb + 1) = xmax
do i = 1, nb
   out%mids(i) = 0.5_dp * (out%breaks(i) + out%breaks(i + 1))
end do
do i = 1, size(x)
   if (x(i) < out%breaks(1) .or. x(i) > out%breaks(nb + 1)) cycle
   j = ceiling((x(i) - out%breaks(1)) / width)
   if (j < 1) j = 1
   if (j > nb) j = nb
   out%counts(j) = out%counts(j) + 1
end do
total = real(sum(out%counts), kind=dp) * width
if (total > 0.0_dp) out%density = real(out%counts, kind=dp) / total
if (ignored_plot .and. size(out%counts) < 0) out%counts = out%counts
end function hist_nbreaks

pure function hist_breaks_real(x, breaks, plot) result(out)
! Compute a small R-like histogram object using explicit numeric breaks.
real(kind=dp), intent(in) :: x(:), breaks(:)
logical, intent(in), optional :: plot
type(hist_result_t) :: out
integer :: nb, i, j
real(kind=dp) :: total, width
logical :: ignored_plot
ignored_plot = .false.
if (present(plot)) ignored_plot = plot
nb = max(0, size(breaks) - 1)
allocate(out%breaks(size(breaks)), out%mids(nb), out%density(nb), out%counts(nb))
out%breaks = breaks
out%counts = 0
out%density = 0.0_dp
if (nb <= 0) return
do i = 1, nb
   out%mids(i) = 0.5_dp * (out%breaks(i) + out%breaks(i + 1))
end do
do i = 1, size(x)
   if (x(i) < out%breaks(1) .or. x(i) > out%breaks(nb + 1)) cycle
   j = 0
   if (x(i) == out%breaks(1)) then
      j = 1
   else
      do while (j < nb)
         j = j + 1
         if (x(i) > out%breaks(j) .and. x(i) <= out%breaks(j + 1)) exit
      end do
   end if
   if (j >= 1 .and. j <= nb) out%counts(j) = out%counts(j) + 1
end do
total = real(sum(out%counts), kind=dp)
if (total > 0.0_dp) then
   do i = 1, nb
      width = out%breaks(i + 1) - out%breaks(i)
      if (width > 0.0_dp) out%density(i) = real(out%counts(i), kind=dp) / (total * width)
   end do
end if
if (ignored_plot .and. size(out%counts) < 0) out%counts = out%counts
end function hist_breaks_real

subroutine print_hist(h)
! Print a compact summary for an R-like hist result.
type(hist_result_t), intent(in) :: h
write(*,'(a)') "$breaks"
call print_real_vector(h%breaks)
write(*,'(a)') "$counts"
call print_integer_vector(h%counts)
write(*,'(a)') "$density"
call print_real_vector(h%density)
write(*,'(a)') "$mids"
call print_real_vector(h%mids)
end subroutine print_hist

pure function quantile(x, probs, names, type, na_rm) result(out)
! Compute Type-7 quantiles for a numeric vector.
! names and type are accepted for API compatibility; this subset always uses Type 7.
real(kind=dp), intent(in) :: x(:) ! sample values
real(kind=dp), intent(in) :: probs(:) ! probabilities, clamped to [0,1]
logical, intent(in), optional :: names ! accepted for compatibility; ignored
integer, intent(in), optional :: type ! accepted for compatibility; ignored
logical, intent(in), optional :: na_rm ! remove missing values before calculation
real(kind=dp), allocatable :: out(:), xs(:)
integer :: n, i, j
real(kind=dp) :: p, h, g
logical :: remove_missing
allocate(out(size(probs)))
remove_missing = .false.
if (present(na_rm)) remove_missing = na_rm
if (.not. remove_missing .and. any(x /= x)) then
   out = ieee_value(0.0_dp, ieee_quiet_nan)
   return
end if
if (remove_missing) then
   xs = pack(x, x == x)
else
   xs = x
end if
n = size(xs)
if (n <= 0) then
   out = ieee_value(0.0_dp, ieee_quiet_nan)
   return
end if
call sort_increasing(xs)
do i = 1, size(probs)
   if (probs(i) /= probs(i)) then
      out(i) = ieee_value(0.0_dp, ieee_quiet_nan)
      cycle
   end if
   p = min(1.0_dp, max(0.0_dp, probs(i)))
   h = (n - 1) * p + 1.0_dp
   j = int(floor(h))
   g = h - j
   if (j < 1) then
      out(i) = xs(1)
   else if (j >= n) then
      out(i) = xs(n)
   else if (g == 0.0_dp) then
      out(i) = xs(j)
   else
      out(i) = (1.0_dp - g) * xs(j) + g * xs(j + 1)
   end if
end do
! names/type accepted for API compatibility in this subset.
if (present(names)) continue
if (present(type)) continue
end function quantile

pure function median(x, na_rm) result(out)
! Compute the median of a numeric vector.
real(kind=dp), intent(in) :: x(:) ! sample values
logical, intent(in), optional :: na_rm ! remove missing values before calculation
real(kind=dp) :: out
real(kind=dp), allocatable :: xs(:)
integer :: n, mid
logical :: remove_missing
remove_missing = .false.
if (present(na_rm)) remove_missing = na_rm
if (.not. remove_missing .and. any(x /= x)) then
   out = ieee_value(0.0_dp, ieee_quiet_nan)
   return
end if
if (remove_missing) then
   xs = pack(x, x == x)
else
   xs = x
end if
n = size(xs)
if (n <= 0) then
   out = ieee_value(0.0_dp, ieee_quiet_nan)
   return
end if
call sort_increasing(xs)
mid = (n + 1) / 2
if (mod(n, 2) == 1) then
   out = xs(mid)
else
   out = 0.5_dp * (xs(mid) + xs(mid + 1))
end if
end function median

pure function summary_vec(x) result(out)
! Return R-like numeric summary: Min, 1st Qu., Median, Mean, 3rd Qu., Max.
real(kind=dp), intent(in) :: x(:) ! sample values
real(kind=dp), allocatable :: out(:), qs(:), clean(:)
integer :: n
clean = pack(x, x == x)
n = size(clean)
allocate(out(6))
if (n <= 0) then
   out = ieee_value(0.0_dp, ieee_quiet_nan)
   return
end if
qs = quantile(clean, [0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp], .false., 7)
out = [qs(1), qs(2), qs(3), sum(clean) / real(n, kind=dp), qs(4), qs(5)]
end function summary_vec

pure function summary_mat(x) result(out)
! Return R-like numeric summary over all elements of a numeric matrix.
real(kind=dp), intent(in) :: x(:,:) ! sample matrix
real(kind=dp), allocatable :: out(:)
out = summary_vec(reshape(x, [size(x)]))
end function summary_mat

pure function dnorm_vec(x, mean, sd, log_) result(out)
! Evaluate normal density (or log-density) elementwise.
! Defaults: mean=0, sd=1, log_=.false.; nonpositive sd returns NaN values.
real(kind=dp), intent(in) :: x(:) ! quantiles
real(kind=dp), intent(in), optional :: mean ! distribution mean
real(kind=dp), intent(in), optional :: sd ! distribution standard deviation
logical, intent(in), optional :: log_ ! return log-density flag
real(kind=dp), allocatable :: out(:), z(:)
real(kind=dp) :: mu, sig
logical :: l
l = .false.
if (present(log_)) l = log_
mu = 0.0_dp
sig = 1.0_dp
if (present(mean)) mu = mean
if (present(sd)) sig = sd
if (sig <= 0.0_dp) then
   allocate(out(size(x)))
   out = ieee_value(0.0_dp, ieee_quiet_nan)
   return
end if
z = (x - mu) / sig
out = -0.5_dp * log(2.0_dp * acos(-1.0_dp)) - log(sig) - 0.5_dp * z**2
if (.not. l) out = exp(out)
end function dnorm_vec

pure function dnorm_scalar(x, mean, sd, log_) result(out)
! Evaluate normal density (or log-density) for one value.
! Defaults: mean=0, sd=1, log_=.false.; nonpositive sd returns NaN.
real(kind=dp), intent(in) :: x ! quantile
real(kind=dp), intent(in), optional :: mean ! distribution mean
real(kind=dp), intent(in), optional :: sd ! distribution standard deviation
logical, intent(in), optional :: log_ ! return log-density flag
real(kind=dp) :: out
real(kind=dp) :: mu, sig
real(kind=dp), allocatable :: tmp(:)
mu = 0.0_dp
sig = 1.0_dp
if (present(mean)) mu = mean
if (present(sd)) sig = sd
tmp = dnorm_vec([x], mu, sig, log_)
out = tmp(1)
end function dnorm_scalar

pure function dnorm_int_scalar(x, mean, sd, log_) result(out)
! Evaluate normal density (or log-density) for one integer value.
! Defaults: mean=0, sd=1, log_=.false.; nonpositive sd returns NaN.
integer, intent(in) :: x ! quantile
real(kind=dp), intent(in), optional :: mean ! distribution mean
real(kind=dp), intent(in), optional :: sd ! distribution standard deviation
logical, intent(in), optional :: log_ ! return log-density flag
real(kind=dp) :: out
out = dnorm_scalar(real(x, kind=dp), mean, sd, log_)
end function dnorm_int_scalar

pure function tail(x, n) result(out)
! Return the last n elements of a vector.
real(kind=dp), intent(in) :: x(:) ! source vector
integer, intent(in) :: n ! requested tail length, clamped to vector length
real(kind=dp), allocatable :: out(:)
integer :: m, n0
n0 = size(x)
if (n < 0) then
   m = max(0, n0 + n)
else
   m = min(n, n0)
end if
allocate(out(m))
if (m > 0) out = x(n0 - m + 1:n0)
end function tail

pure function cbind2(a, b) result(out)
! Bind two vectors as columns of a 2D array.
real(kind=dp), intent(in) :: a(:) ! first column source
real(kind=dp), intent(in) :: b(:) ! second column source
real(kind=dp), allocatable :: out(:,:)
integer :: i, j, n, nc
n = max(size(a), size(b))
if (n == 0) then
   nc = 2
else
   nc = merge(1, 0, size(a) > 0) + merge(1, 0, size(b) > 0)
end if
allocate(out(n, nc))
j = 0
if (size(a) > 0) then
   j = j + 1
   do i = 1, n
      out(i, j) = a(modulo(i - 1, size(a)) + 1)
   end do
end if
if (size(b) > 0) then
   j = j + 1
   do i = 1, n
      out(i, j) = b(modulo(i - 1, size(b)) + 1)
   end do
end if
end function cbind2

pure function cbind(r1, r2, r3) result(out)
! Bind two or three vectors as columns of a 2D array.
real(kind=dp), intent(in) :: r1(:) ! first column source
real(kind=dp), intent(in) :: r2(:) ! second column source
real(kind=dp), intent(in), optional :: r3(:) ! optional third column source
real(kind=dp), allocatable :: out(:,:)
integer :: i, j, n, nc
if (.not. present(r3)) then
   out = cbind2(r1, r2)
   return
end if
n = max(size(r1), max(size(r2), size(r3)))
if (n == 0) then
   nc = 3
else
   nc = merge(1, 0, size(r1) > 0) + merge(1, 0, size(r2) > 0) + merge(1, 0, size(r3) > 0)
end if
allocate(out(n, nc))
j = 0
if (size(r1) > 0) then
   j = j + 1
   do i = 1, n
      out(i, j) = r1(modulo(i - 1, size(r1)) + 1)
   end do
end if
if (size(r2) > 0) then
   j = j + 1
   do i = 1, n
      out(i, j) = r2(modulo(i - 1, size(r2)) + 1)
   end do
end if
if (size(r3) > 0) then
   j = j + 1
   do i = 1, n
      out(i, j) = r3(modulo(i - 1, size(r3)) + 1)
   end do
end if
end function cbind

pure function rbind_vec(a, b) result(out)
! Bind two vectors as rows of a 2D array.
real(kind=dp), intent(in) :: a(:) ! first row source
real(kind=dp), intent(in) :: b(:) ! second row source
real(kind=dp), allocatable :: out(:,:)
integer :: i, j, n, nr
n = max(size(a), size(b))
if (n == 0) then
   nr = 2
else
   nr = merge(1, 0, size(a) > 0) + merge(1, 0, size(b) > 0)
end if
allocate(out(nr, n))
j = 0
if (size(a) > 0) then
   j = j + 1
   do i = 1, n
      out(j, i) = a(modulo(i - 1, size(a)) + 1)
   end do
end if
if (size(b) > 0) then
   j = j + 1
   do i = 1, n
      out(j, i) = b(modulo(i - 1, size(b)) + 1)
   end do
end if
end function rbind_vec

pure function rbind_mat(a, b) result(out)
! Bind two matrices by concatenating rows.
real(kind=dp), intent(in) :: a(:,:) ! top matrix
real(kind=dp), intent(in) :: b(:,:) ! bottom matrix
real(kind=dp), allocatable :: out(:,:)
integer :: n1, n2, ncol
n1 = size(a, 1)
n2 = size(b, 1)
ncol = min(size(a, 2), size(b, 2))
allocate(out(n1 + n2, ncol))
if (ncol == 0 .or. (n1 == 0 .and. n2 == 0)) return
if (n1 > 0) out(1:n1, 1:ncol) = a(:, 1:ncol)
if (n2 > 0) out(n1 + 1:n1 + n2, 1:ncol) = b(:, 1:ncol)
end function rbind_mat

pure function rbind_vec_mat(a, b) result(out)
! Bind a vector as a row and a matrix below it.
real(kind=dp), intent(in) :: a(:) ! top row source
real(kind=dp), intent(in) :: b(:,:) ! rows appended below
real(kind=dp), allocatable :: out(:,:)
integer :: i, ncol, nrow
ncol = size(b, 2)
nrow = size(b, 1)
if (size(a) <= 0) then
   out = b
   return
end if
allocate(out(nrow + 1, ncol))
if (ncol == 0) return
do i = 1, ncol
   out(1, i) = a(modulo(i - 1, size(a)) + 1)
end do
if (nrow > 0) out(2:, 1:ncol) = b(:, 1:ncol)
end function rbind_vec_mat

pure function rbind_mat_vec(a, b) result(out)
! Bind a matrix above a vector row.
real(kind=dp), intent(in) :: a(:,:) ! rows kept above
real(kind=dp), intent(in) :: b(:) ! bottom row source
real(kind=dp), allocatable :: out(:,:)
integer :: i, ncol, nrow
ncol = size(a, 2)
nrow = size(a, 1)
if (size(b) <= 0) then
   out = a
   return
end if
allocate(out(nrow + 1, ncol))
if (ncol == 0) return
if (nrow > 0) out(1:nrow, 1:ncol) = a(:, 1:ncol)
do i = 1, ncol
   out(nrow + 1, i) = b(modulo(i - 1, size(b)) + 1)
end do
end function rbind_mat_vec

pure function rbind_int_vec(a, b) result(out)
! Bind two integer vectors as rows of a 2D array.
integer, intent(in) :: a(:)
integer, intent(in) :: b(:)
integer, allocatable :: out(:,:)
integer :: i, j, n, nr
n = max(size(a), size(b))
if (n == 0) then
   nr = 2
else
   nr = merge(1, 0, size(a) > 0) + merge(1, 0, size(b) > 0)
end if
allocate(out(nr, n))
j = 0
if (size(a) > 0) then
   j = j + 1
   do i = 1, n
      out(j, i) = a(modulo(i - 1, size(a)) + 1)
   end do
end if
if (size(b) > 0) then
   j = j + 1
   do i = 1, n
      out(j, i) = b(modulo(i - 1, size(b)) + 1)
   end do
end if
end function rbind_int_vec

pure function rbind_int_mat(a, b) result(out)
! Bind two integer matrices by concatenating rows.
integer, intent(in) :: a(:,:)
integer, intent(in) :: b(:,:)
integer, allocatable :: out(:,:)
integer :: n1, n2, ncol
n1 = size(a, 1)
n2 = size(b, 1)
ncol = min(size(a, 2), size(b, 2))
allocate(out(n1 + n2, ncol))
if (ncol == 0 .or. (n1 == 0 .and. n2 == 0)) return
if (n1 > 0) out(1:n1, 1:ncol) = a(:, 1:ncol)
if (n2 > 0) out(n1 + 1:n1 + n2, 1:ncol) = b(:, 1:ncol)
end function rbind_int_mat

pure function rbind_int_vec_mat(a, b) result(out)
! Bind an integer vector as a row and an integer matrix below it.
integer, intent(in) :: a(:)
integer, intent(in) :: b(:,:)
integer, allocatable :: out(:,:)
integer :: i, ncol, nrow
ncol = size(b, 2)
nrow = size(b, 1)
if (size(a) <= 0) then
   out = b
   return
end if
allocate(out(nrow + 1, ncol))
if (ncol == 0) return
do i = 1, ncol
   out(1, i) = a(modulo(i - 1, size(a)) + 1)
end do
if (nrow > 0) out(2:, 1:ncol) = b(:, 1:ncol)
end function rbind_int_vec_mat

pure function rbind_int_mat_vec(a, b) result(out)
! Bind an integer matrix above an integer vector row.
integer, intent(in) :: a(:,:)
integer, intent(in) :: b(:)
integer, allocatable :: out(:,:)
integer :: i, ncol, nrow
ncol = size(a, 2)
nrow = size(a, 1)
if (size(b) <= 0) then
   out = a
   return
end if
allocate(out(nrow + 1, ncol))
if (ncol == 0) return
if (nrow > 0) out(1:nrow, 1:ncol) = a(:, 1:ncol)
do i = 1, ncol
   out(nrow + 1, i) = b(modulo(i - 1, size(b)) + 1)
end do
end function rbind_int_mat_vec

pure function matrix_real(x, nrow, ncol) result(out)
! Build matrix with R-like recycling in column-major order.
real(kind=dp), intent(in) :: x(:) ! values recycled into the matrix
integer, intent(in) :: nrow ! requested row count
integer, intent(in), optional :: ncol ! requested column count; inferred from x and nrow if absent
real(kind=dp), allocatable :: out(:,:)
real(kind=dp), allocatable :: buf(:)
integer :: i, need_n, nx, nc
nx = size(x)
if (present(ncol)) then
   nc = ncol
else if (nrow > 0) then
   nc = (nx + nrow - 1) / nrow
else
   nc = 0
end if
if (nrow < 0 .or. nc < 0) then
   allocate(out(0, 0))
   return
end if
need_n = nrow * nc
if (need_n <= 0) then
   allocate(out(nrow, nc))
   return
end if
if (nx <= 0) then
   allocate(out(nrow, nc))
   out = r_na_real()
   return
end if
allocate(buf(need_n))
do i = 1, need_n
   buf(i) = x(modulo(i - 1, nx) + 1)
end do
out = reshape(buf, [nrow, nc])
end function matrix_real

pure function matrix_int(x, nrow, ncol) result(out)
! Integer variant of matrix() with R-like recycling.
integer, intent(in) :: x(:) ! values recycled into the matrix
integer, intent(in) :: nrow ! requested row count
integer, intent(in), optional :: ncol ! requested column count; inferred from x and nrow if absent
integer, allocatable :: out(:,:)
integer, allocatable :: buf(:)
integer :: i, need_n, nx, nc
nx = size(x)
if (present(ncol)) then
   nc = ncol
else if (nrow > 0) then
   nc = (nx + nrow - 1) / nrow
else
   nc = 0
end if
if (nrow < 0 .or. nc < 0) then
   allocate(out(0, 0))
   return
end if
need_n = nrow * nc
if (need_n <= 0) then
   allocate(out(nrow, nc))
   return
end if
if (nx <= 0) then
   allocate(out(nrow, nc))
   out = -huge(0)
   return
end if
allocate(buf(need_n))
do i = 1, need_n
   buf(i) = x(modulo(i - 1, nx) + 1)
end do
out = reshape(buf, [nrow, nc])
end function matrix_int

pure function r_matmul_vv_real(a, b) result(out)
! Matrix-product helper for real vectors: dot product.
real(kind=dp), intent(in) :: a(:) ! left operand vector
real(kind=dp), intent(in) :: b(:) ! right operand vector
real(kind=dp) :: out
out = dot_product(a, b)
end function r_matmul_vv_real

pure function r_matmul_vv_int(a, b) result(out)
! Matrix-product helper for integer vectors: dot product (real result).
integer, intent(in) :: a(:) ! left operand vector
integer, intent(in) :: b(:) ! right operand vector
real(kind=dp) :: out
out = dot_product(real(a, kind=dp), real(b, kind=dp))
end function r_matmul_vv_int

pure function r_matmul_vv_real_int(a, b) result(out)
! Matrix-product helper for mixed real/int vectors: dot product.
real(kind=dp), intent(in) :: a(:) ! left operand vector
integer, intent(in) :: b(:) ! right operand vector
real(kind=dp) :: out
out = dot_product(a, real(b, kind=dp))
end function r_matmul_vv_real_int

pure function r_matmul_vv_int_real(a, b) result(out)
! Matrix-product helper for mixed int/real vectors: dot product.
integer, intent(in) :: a(:) ! left operand vector
real(kind=dp), intent(in) :: b(:) ! right operand vector
real(kind=dp) :: out
out = dot_product(real(a, kind=dp), b)
end function r_matmul_vv_int_real

pure function r_matmul_mv_real(a, b) result(out)
! Matrix-product helper for real matrix-vector multiplication.
real(kind=dp), intent(in) :: a(:,:) ! left operand matrix
real(kind=dp), intent(in) :: b(:) ! right operand vector
real(kind=dp), allocatable :: out(:)
allocate(out(size(a, 1)))
out = matmul(a, b)
end function r_matmul_mv_real

pure function r_matmul_mv_int(a, b) result(out)
! Matrix-product helper for integer matrix-vector multiplication.
integer, intent(in) :: a(:,:) ! left operand matrix
integer, intent(in) :: b(:) ! right operand vector
real(kind=dp), allocatable :: out(:)
allocate(out(size(a, 1)))
out = matmul(real(a, kind=dp), real(b, kind=dp))
end function r_matmul_mv_int

pure function r_matmul_mv_real_int(a, b) result(out)
! Matrix-product helper for mixed real/int matrix-vector multiplication.
real(kind=dp), intent(in) :: a(:,:) ! left operand matrix
integer, intent(in) :: b(:) ! right operand vector
real(kind=dp), allocatable :: out(:)
allocate(out(size(a, 1)))
out = matmul(a, real(b, kind=dp))
end function r_matmul_mv_real_int

pure function r_matmul_mv_int_real(a, b) result(out)
! Matrix-product helper for mixed int/real matrix-vector multiplication.
integer, intent(in) :: a(:,:) ! left operand matrix
real(kind=dp), intent(in) :: b(:) ! right operand vector
real(kind=dp), allocatable :: out(:)
allocate(out(size(a, 1)))
out = matmul(real(a, kind=dp), b)
end function r_matmul_mv_int_real

pure function r_matmul_mv_complex(a, b) result(out)
! Matrix-product helper for complex matrix-vector multiplication.
complex(kind=dp), intent(in) :: a(:,:) ! left operand matrix
complex(kind=dp), intent(in) :: b(:) ! right operand vector
complex(kind=dp), allocatable :: out(:)
allocate(out(size(a, 1)))
out = matmul(a, b)
end function r_matmul_mv_complex

pure function r_matmul_mv_real_complex(a, b) result(out)
! Matrix-product helper for real matrix and complex vector multiplication.
real(kind=dp), intent(in) :: a(:,:) ! left operand matrix
complex(kind=dp), intent(in) :: b(:) ! right operand vector
complex(kind=dp), allocatable :: out(:)
allocate(out(size(a, 1)))
out = matmul(cmplx(a, 0.0_dp, kind=dp), b)
end function r_matmul_mv_real_complex

pure function r_matmul_mv_int_complex(a, b) result(out)
! Matrix-product helper for integer matrix and complex vector multiplication.
integer, intent(in) :: a(:,:) ! left operand matrix
complex(kind=dp), intent(in) :: b(:) ! right operand vector
complex(kind=dp), allocatable :: out(:)
allocate(out(size(a, 1)))
out = matmul(cmplx(real(a, kind=dp), 0.0_dp, kind=dp), b)
end function r_matmul_mv_int_complex

pure function r_matmul_vm_real(a, b) result(out)
! Matrix-product helper for real vector-matrix multiplication.
real(kind=dp), intent(in) :: a(:) ! left operand vector
real(kind=dp), intent(in) :: b(:,:) ! right operand matrix
real(kind=dp), allocatable :: out(:)
allocate(out(size(b, 2)))
out = matmul(a, b)
end function r_matmul_vm_real

pure function r_matmul_vm_int(a, b) result(out)
! Matrix-product helper for integer vector-matrix multiplication.
integer, intent(in) :: a(:) ! left operand vector
integer, intent(in) :: b(:,:) ! right operand matrix
real(kind=dp), allocatable :: out(:)
allocate(out(size(b, 2)))
out = matmul(real(a, kind=dp), real(b, kind=dp))
end function r_matmul_vm_int

pure function r_matmul_vm_real_int(a, b) result(out)
! Matrix-product helper for mixed real/int vector-matrix multiplication.
real(kind=dp), intent(in) :: a(:) ! left operand vector
integer, intent(in) :: b(:,:) ! right operand matrix
real(kind=dp), allocatable :: out(:)
allocate(out(size(b, 2)))
out = matmul(a, real(b, kind=dp))
end function r_matmul_vm_real_int

pure function r_matmul_vm_int_real(a, b) result(out)
! Matrix-product helper for mixed int/real vector-matrix multiplication.
integer, intent(in) :: a(:) ! left operand vector
real(kind=dp), intent(in) :: b(:,:) ! right operand matrix
real(kind=dp), allocatable :: out(:)
allocate(out(size(b, 2)))
out = matmul(real(a, kind=dp), b)
end function r_matmul_vm_int_real

pure function r_matmul_mm_real(a, b) result(out)
! Matrix-product helper for real matrix-matrix multiplication.
real(kind=dp), intent(in) :: a(:,:) ! left operand matrix
real(kind=dp), intent(in) :: b(:,:) ! right operand matrix
real(kind=dp), allocatable :: out(:,:)
allocate(out(size(a, 1), size(b, 2)))
out = matmul(a, b)
end function r_matmul_mm_real

pure function r_matmul_mm_int(a, b) result(out)
! Matrix-product helper for integer matrix-matrix multiplication.
integer, intent(in) :: a(:,:) ! left operand matrix
integer, intent(in) :: b(:,:) ! right operand matrix
real(kind=dp), allocatable :: out(:,:)
allocate(out(size(a, 1), size(b, 2)))
out = matmul(real(a, kind=dp), real(b, kind=dp))
end function r_matmul_mm_int

pure function r_matmul_mm_real_int(a, b) result(out)
! Matrix-product helper for mixed real/int matrix-matrix multiplication.
real(kind=dp), intent(in) :: a(:,:) ! left operand matrix
integer, intent(in) :: b(:,:) ! right operand matrix
real(kind=dp), allocatable :: out(:,:)
allocate(out(size(a, 1), size(b, 2)))
out = matmul(a, real(b, kind=dp))
end function r_matmul_mm_real_int

pure function r_matmul_mm_int_real(a, b) result(out)
! Matrix-product helper for mixed int/real matrix-matrix multiplication.
integer, intent(in) :: a(:,:) ! left operand matrix
real(kind=dp), intent(in) :: b(:,:) ! right operand matrix
real(kind=dp), allocatable :: out(:,:)
allocate(out(size(a, 1), size(b, 2)))
out = matmul(real(a, kind=dp), b)
end function r_matmul_mm_int_real

pure function r_matmul_mm_complex(a, b) result(out)
! Matrix-product helper for complex matrix-matrix multiplication.
complex(kind=dp), intent(in) :: a(:,:) ! left operand matrix
complex(kind=dp), intent(in) :: b(:,:) ! right operand matrix
complex(kind=dp), allocatable :: out(:,:)
allocate(out(size(a, 1), size(b, 2)))
out = matmul(a, b)
end function r_matmul_mm_complex

pure function r_matmul_mm_complex_real(a, b) result(out)
! Matrix-product helper for complex matrix and real matrix multiplication.
complex(kind=dp), intent(in) :: a(:,:) ! left operand matrix
real(kind=dp), intent(in) :: b(:,:) ! right operand matrix
complex(kind=dp), allocatable :: out(:,:)
allocate(out(size(a, 1), size(b, 2)))
out = matmul(a, cmplx(b, 0.0_dp, kind=dp))
end function r_matmul_mm_complex_real

pure function r_matmul_mm_real_complex(a, b) result(out)
! Matrix-product helper for real matrix and complex matrix multiplication.
real(kind=dp), intent(in) :: a(:,:) ! left operand matrix
complex(kind=dp), intent(in) :: b(:,:) ! right operand matrix
complex(kind=dp), allocatable :: out(:,:)
allocate(out(size(a, 1), size(b, 2)))
out = matmul(cmplx(a, 0.0_dp, kind=dp), b)
end function r_matmul_mm_real_complex

function r_add_vv(a, b) result(out)
! Recycle and add two vectors (R-style recycling).
real(kind=dp), intent(in) :: a(:) ! left operand vector
real(kind=dp), intent(in) :: b(:) ! right operand vector
real(kind=dp), allocatable :: out(:)
integer :: i, n, na, nb
na = size(a)
nb = size(b)
if (na <= 0 .or. nb <= 0) then
   allocate(out(0))
   return
end if
n = max(na, nb)
allocate(out(n))
call maybe_warn_recycle("x + y", na, nb)
do i = 1, n
   out(i) = a(modulo(i - 1, na) + 1) + b(modulo(i - 1, nb) + 1)
end do
end function r_add_vv

pure function r_add_vs(a, b) result(out)
! Add scalar to vector.
real(kind=dp), intent(in) :: a(:) ! left operand vector
real(kind=dp), intent(in) :: b ! right operand scalar
real(kind=dp), allocatable :: out(:)
integer :: n
n = size(a)
allocate(out(n))
if (n > 0) out = a + b
end function r_add_vs

pure function r_add_sv(a, b) result(out)
! Add vector to scalar.
real(kind=dp), intent(in) :: a ! left operand scalar
real(kind=dp), intent(in) :: b(:) ! right operand vector
real(kind=dp), allocatable :: out(:)
integer :: n
n = size(b)
allocate(out(n))
if (n > 0) out = a + b
end function r_add_sv

function r_sub_vv(a, b) result(out)
! Recycle and subtract two vectors (a - b).
real(kind=dp), intent(in) :: a(:) ! left operand vector
real(kind=dp), intent(in) :: b(:) ! right operand vector
real(kind=dp), allocatable :: out(:)
integer :: i, n, na, nb
na = size(a)
nb = size(b)
if (na <= 0 .or. nb <= 0) then
   allocate(out(0))
   return
end if
n = max(na, nb)
allocate(out(n))
call maybe_warn_recycle("x - y", na, nb)
do i = 1, n
   out(i) = a(modulo(i - 1, na) + 1) - b(modulo(i - 1, nb) + 1)
end do
end function r_sub_vv

pure function r_sub_vs(a, b) result(out)
! Subtract scalar from vector.
real(kind=dp), intent(in) :: a(:) ! left operand vector
real(kind=dp), intent(in) :: b ! right operand scalar
real(kind=dp), allocatable :: out(:)
integer :: n
n = size(a)
allocate(out(n))
if (n > 0) out = a - b
end function r_sub_vs

pure function r_sub_sv(a, b) result(out)
! Subtract vector from scalar.
real(kind=dp), intent(in) :: a ! left operand scalar
real(kind=dp), intent(in) :: b(:) ! right operand vector
real(kind=dp), allocatable :: out(:)
integer :: n
n = size(b)
allocate(out(n))
if (n > 0) out = a - b
end function r_sub_sv

function r_mul_vv(a, b) result(out)
! Recycle and multiply two vectors.
real(kind=dp), intent(in) :: a(:) ! left operand vector
real(kind=dp), intent(in) :: b(:) ! right operand vector
real(kind=dp), allocatable :: out(:)
integer :: i, n, na, nb
na = size(a)
nb = size(b)
if (na <= 0 .or. nb <= 0) then
   allocate(out(0))
   return
end if
n = max(na, nb)
allocate(out(n))
call maybe_warn_recycle("x * y", na, nb)
do i = 1, n
   out(i) = a(modulo(i - 1, na) + 1) * b(modulo(i - 1, nb) + 1)
end do
end function r_mul_vv

pure function r_mul_vs(a, b) result(out)
! Multiply vector by scalar.
real(kind=dp), intent(in) :: a(:) ! left operand vector
real(kind=dp), intent(in) :: b ! right operand scalar
real(kind=dp), allocatable :: out(:)
integer :: n
n = size(a)
allocate(out(n))
if (n > 0) out = a * b
end function r_mul_vs

pure function r_mul_sv(a, b) result(out)
! Multiply scalar by vector.
real(kind=dp), intent(in) :: a ! left operand scalar
real(kind=dp), intent(in) :: b(:) ! right operand vector
real(kind=dp), allocatable :: out(:)
integer :: n
n = size(b)
allocate(out(n))
if (n > 0) out = a * b
end function r_mul_sv

function r_div_vv(a, b) result(out)
! Recycle and divide two vectors (a / b).
real(kind=dp), intent(in) :: a(:) ! left operand vector
real(kind=dp), intent(in) :: b(:) ! right operand vector
real(kind=dp), allocatable :: out(:)
integer :: i, n, na, nb
na = size(a)
nb = size(b)
if (na <= 0 .or. nb <= 0) then
   allocate(out(0))
   return
end if
n = max(na, nb)
allocate(out(n))
call maybe_warn_recycle("x / y", na, nb)
do i = 1, n
   out(i) = a(modulo(i - 1, na) + 1) / b(modulo(i - 1, nb) + 1)
end do
end function r_div_vv

pure function r_div_vs(a, b) result(out)
! Divide vector by scalar.
real(kind=dp), intent(in) :: a(:) ! left operand vector
real(kind=dp), intent(in) :: b ! right operand scalar
real(kind=dp), allocatable :: out(:)
integer :: n
n = size(a)
allocate(out(n))
if (n > 0) out = a / b
end function r_div_vs

pure function r_div_sv(a, b) result(out)
! Divide scalar by vector.
real(kind=dp), intent(in) :: a ! left operand scalar
real(kind=dp), intent(in) :: b(:) ! right operand vector
real(kind=dp), allocatable :: out(:)
integer :: n
n = size(b)
allocate(out(n))
if (n > 0) out = a / b
end function r_div_sv

function r_add_mv(a, b) result(out)
! Recycle a vector across a matrix in R column-major order and add.
real(kind=dp), intent(in) :: a(:,:) ! left operand matrix
real(kind=dp), intent(in) :: b(:) ! right operand vector
real(kind=dp), allocatable :: out(:,:)
integer :: i, na, nb
na = size(a)
nb = size(b)
if (na <= 0) then
   allocate(out(size(a, 1), size(a, 2)))
   return
end if
if (nb <= 0) then
   allocate(out(0, 0))
   return
end if
allocate(out(size(a, 1), size(a, 2)))
call maybe_warn_recycle("x + y", na, nb)
do i = 1, na
   out(modulo(i - 1, size(a, 1)) + 1, ((i - 1) / size(a, 1)) + 1) = &
      a(modulo(i - 1, size(a, 1)) + 1, ((i - 1) / size(a, 1)) + 1) + b(modulo(i - 1, nb) + 1)
end do
end function r_add_mv

function r_add_vm(a, b) result(out)
! Recycle a vector across a matrix in R column-major order and add.
real(kind=dp), intent(in) :: a(:) ! left operand vector
real(kind=dp), intent(in) :: b(:,:) ! right operand matrix
real(kind=dp), allocatable :: out(:,:)
out = r_add_mv(b, a)
end function r_add_vm

function r_sub_mv(a, b) result(out)
! Recycle a vector across a matrix in R column-major order and subtract.
real(kind=dp), intent(in) :: a(:,:) ! left operand matrix
real(kind=dp), intent(in) :: b(:) ! right operand vector
real(kind=dp), allocatable :: out(:,:)
integer :: i, na, nb
na = size(a)
nb = size(b)
if (na <= 0) then
   allocate(out(size(a, 1), size(a, 2)))
   return
end if
if (nb <= 0) then
   allocate(out(0, 0))
   return
end if
allocate(out(size(a, 1), size(a, 2)))
call maybe_warn_recycle("x - y", na, nb)
do i = 1, na
   out(modulo(i - 1, size(a, 1)) + 1, ((i - 1) / size(a, 1)) + 1) = &
      a(modulo(i - 1, size(a, 1)) + 1, ((i - 1) / size(a, 1)) + 1) - b(modulo(i - 1, nb) + 1)
end do
end function r_sub_mv

function r_sub_vm(a, b) result(out)
! Recycle a vector across a matrix in R column-major order and subtract from it.
real(kind=dp), intent(in) :: a(:) ! left operand vector
real(kind=dp), intent(in) :: b(:,:) ! right operand matrix
real(kind=dp), allocatable :: out(:,:)
integer :: i, nb, na
na = size(a)
nb = size(b)
if (nb <= 0) then
   allocate(out(size(b, 1), size(b, 2)))
   return
end if
if (na <= 0) then
   allocate(out(0, 0))
   return
end if
allocate(out(size(b, 1), size(b, 2)))
call maybe_warn_recycle("x - y", na, nb)
do i = 1, nb
   out(modulo(i - 1, size(b, 1)) + 1, ((i - 1) / size(b, 1)) + 1) = &
      a(modulo(i - 1, na) + 1) - b(modulo(i - 1, size(b, 1)) + 1, ((i - 1) / size(b, 1)) + 1)
end do
end function r_sub_vm

function r_mul_mv(a, b) result(out)
! Recycle a vector across a matrix in R column-major order and multiply.
real(kind=dp), intent(in) :: a(:,:) ! left operand matrix
real(kind=dp), intent(in) :: b(:) ! right operand vector
real(kind=dp), allocatable :: out(:,:)
integer :: i, na, nb
na = size(a)
nb = size(b)
if (na <= 0) then
   allocate(out(size(a, 1), size(a, 2)))
   return
end if
if (nb <= 0) then
   allocate(out(0, 0))
   return
end if
allocate(out(size(a, 1), size(a, 2)))
call maybe_warn_recycle("x * y", na, nb)
do i = 1, na
   out(modulo(i - 1, size(a, 1)) + 1, ((i - 1) / size(a, 1)) + 1) = &
      a(modulo(i - 1, size(a, 1)) + 1, ((i - 1) / size(a, 1)) + 1) * b(modulo(i - 1, nb) + 1)
end do
end function r_mul_mv

function r_mul_vm(a, b) result(out)
! Recycle a vector across a matrix in R column-major order and multiply.
real(kind=dp), intent(in) :: a(:) ! left operand vector
real(kind=dp), intent(in) :: b(:,:) ! right operand matrix
real(kind=dp), allocatable :: out(:,:)
out = r_mul_mv(b, a)
end function r_mul_vm

function r_div_mv(a, b) result(out)
! Recycle a vector across a matrix in R column-major order and divide.
real(kind=dp), intent(in) :: a(:,:) ! left operand matrix
real(kind=dp), intent(in) :: b(:) ! right operand vector
real(kind=dp), allocatable :: out(:,:)
integer :: i, na, nb
na = size(a)
nb = size(b)
if (na <= 0) then
   allocate(out(size(a, 1), size(a, 2)))
   return
end if
if (nb <= 0) then
   allocate(out(0, 0))
   return
end if
allocate(out(size(a, 1), size(a, 2)))
call maybe_warn_recycle("x / y", na, nb)
do i = 1, na
   out(modulo(i - 1, size(a, 1)) + 1, ((i - 1) / size(a, 1)) + 1) = &
      a(modulo(i - 1, size(a, 1)) + 1, ((i - 1) / size(a, 1)) + 1) / b(modulo(i - 1, nb) + 1)
end do
end function r_div_mv

function r_div_vm(a, b) result(out)
! Recycle a vector across a matrix in R column-major order and divide by it.
real(kind=dp), intent(in) :: a(:) ! left operand vector
real(kind=dp), intent(in) :: b(:,:) ! right operand matrix
real(kind=dp), allocatable :: out(:,:)
integer :: i, nb, na
na = size(a)
nb = size(b)
if (nb <= 0) then
   allocate(out(size(b, 1), size(b, 2)))
   return
end if
if (na <= 0) then
   allocate(out(0, 0))
   return
end if
allocate(out(size(b, 1), size(b, 2)))
call maybe_warn_recycle("x / y", na, nb)
do i = 1, nb
   out(modulo(i - 1, size(b, 1)) + 1, ((i - 1) / size(b, 1)) + 1) = &
      a(modulo(i - 1, na) + 1) / b(modulo(i - 1, size(b, 1)) + 1, ((i - 1) / size(b, 1)) + 1)
end do
end function r_div_vm

pure function r_array_real(x, dim) result(out)
! Build 2D real array with R-like recycling.
real(kind=dp), intent(in) :: x(:) ! values recycled into the array
integer, intent(in) :: dim(:) ! requested dimensions; first two entries are used
real(kind=dp), allocatable :: out(:,:)
real(kind=dp), allocatable :: buf(:)
integer :: i, need_n, nx, nrow, ncol
if (size(dim) < 2) then
   allocate(out(0, 0))
   return
end if
nrow = dim(1)
ncol = dim(2)
nx = size(x)
if (nrow < 0 .or. ncol < 0) then
   allocate(out(0, 0))
   return
end if
need_n = nrow * ncol
if (need_n <= 0) then
   allocate(out(nrow, ncol))
   return
end if
if (nx <= 0) then
   allocate(out(nrow, ncol))
   out = r_na_real()
   return
end if
allocate(buf(need_n))
do i = 1, need_n
   buf(i) = x(modulo(i - 1, nx) + 1)
end do
out = reshape(buf, [nrow, ncol])
end function r_array_real

pure function r_array_int(x, dim) result(out)
! Build 2D integer array with R-like recycling.
integer, intent(in) :: x(:) ! values recycled into the array
integer, intent(in) :: dim(:) ! requested dimensions; first two entries are used
integer, allocatable :: out(:,:)
integer, allocatable :: buf(:)
integer :: i, need_n, nx, nrow, ncol
if (size(dim) < 2) then
   allocate(out(0, 0))
   return
end if
nrow = dim(1)
ncol = dim(2)
nx = size(x)
if (nrow < 0 .or. ncol < 0) then
   allocate(out(0, 0))
   return
end if
need_n = nrow * ncol
if (need_n <= 0) then
   allocate(out(nrow, ncol))
   return
end if
if (nx <= 0) then
   allocate(out(nrow, ncol))
   out = -huge(0)
   return
end if
allocate(buf(need_n))
do i = 1, need_n
   buf(i) = x(modulo(i - 1, nx) + 1)
end do
out = reshape(buf, [nrow, ncol])
end function r_array_int

pure function r_array_char(x, dim) result(out)
! Build 2D character array with R-like recycling.
character(len=*), intent(in) :: x(:) ! strings recycled into the array
integer, intent(in) :: dim(:) ! requested dimensions; first two entries are used
character(len=:), allocatable :: out(:,:)
character(len=:), allocatable :: buf(:)
integer :: i, need_n, nx, nrow, ncol, lch
if (size(dim) < 2) then
   allocate(character(len=1) :: out(0, 0))
   return
end if
nrow = dim(1)
ncol = dim(2)
nx = size(x)
lch = max(1, len(x))
if (nrow < 0 .or. ncol < 0) then
   allocate(character(len=lch) :: out(0, 0))
   return
end if
need_n = nrow * ncol
if (need_n <= 0) then
   allocate(character(len=lch) :: out(nrow, ncol))
   return
end if
if (nx <= 0) then
   allocate(character(len=lch) :: out(nrow, ncol))
   out = ""
   return
end if
allocate(character(len=lch) :: buf(need_n))
do i = 1, need_n
   buf(i) = x(modulo(i - 1, nx) + 1)
end do
allocate(character(len=lch) :: out(nrow, ncol))
out = reshape(buf, [nrow, ncol])
end function r_array_char

pure function numeric(n) result(out)
! Allocate a length-n real vector initialized to zero.
integer, intent(in) :: n ! requested vector length
real(kind=dp), allocatable :: out(:)
allocate(out(max(0, n)))
if (n > 0) out = 0.0_dp
end function numeric

pure elemental function pmax(a, b) result(out)
! Elementwise maximum of two real scalars.
real(kind=dp), intent(in) :: a ! first value
real(kind=dp), intent(in) :: b ! second value
real(kind=dp) :: out
if (b /= b) then
   out = b
else if (a /= a) then
   out = a
else if (b > a) then
   out = b
else
   out = a
end if
end function pmax

pure elemental function r_round(x, digits) result(out)
! Round to the requested number of decimal digits using R's ties-to-even rule.
real(kind=dp), intent(in) :: x
integer, intent(in) :: digits
real(kind=dp) :: out
real(kind=dp) :: base, frac, rounded, scale, scaled
if (.not. ieee_is_finite(x)) then
   out = x
   return
end if
scale = 10.0_dp ** digits
if (.not. ieee_is_finite(scale)) then
   out = x
   return
else if (scale == 0.0_dp) then
   out = 0.0_dp
   return
end if
scaled = abs(x * scale)
if (.not. ieee_is_finite(scaled)) then
   out = x
   return
end if
base = aint(scaled)
frac = scaled - base
if (frac < 0.5_dp) then
   rounded = base
else if (frac > 0.5_dp) then
   rounded = base + 1.0_dp
else if (modulo(base, 2.0_dp) == 0.0_dp) then
   rounded = base
else
   rounded = base + 1.0_dp
end if
out = sign(rounded / scale, x)
end function r_round

pure function sd_vec(x) result(out)
! Sample standard deviation (n-1 denominator).
real(kind=dp), intent(in) :: x(:) ! sample values
real(kind=dp) :: out, m
integer :: n
n = size(x)
if (n <= 1) then
   out = ieee_value(0.0_dp, ieee_quiet_nan)
   return
end if
m = sum(x) / real(n, kind=dp)
out = sqrt(sum((x - m)**2) / real(n - 1, kind=dp))
end function sd_vec

pure function sd_mat(x) result(out)
! R sd(matrix) treats the matrix as a vector.
real(kind=dp), intent(in) :: x(:,:) ! sample matrix
real(kind=dp) :: out
out = sd_vec(pack(x, .true.))
end function sd_mat

pure function var_vec(x) result(out)
! Sample variance (n-1 denominator).
real(kind=dp), intent(in) :: x(:) ! sample values
real(kind=dp) :: out, m
integer :: n
n = size(x)
if (n <= 1) then
   out = ieee_value(0.0_dp, ieee_quiet_nan)
   return
end if
m = sum(x) / real(n, kind=dp)
out = sum((x - m)**2) / real(n - 1, kind=dp)
end function var_vec

pure function var_mat(x) result(out)
! Column covariance matrix, matching R var(matrix).
real(kind=dp), intent(in) :: x(:,:) ! observations by rows, variables by columns
real(kind=dp), allocatable :: out(:,:)
out = cov_mat(x)
end function var_mat

pure function r_sd(x) result(out)
! Alias for sd(x), used by transpiled code to avoid local-name collisions.
real(kind=dp), intent(in) :: x(:) ! sample values
real(kind=dp) :: out
out = sd(x)
end function r_sd

pure elemental function r_gamma(x) result(out)
! Elementwise gamma function.
real(kind=dp), intent(in) :: x ! function argument
real(kind=dp) :: out
out = gamma(x)
end function r_gamma

pure elemental function r_lgamma(x) result(out)
! Elementwise log gamma function.
real(kind=dp), intent(in) :: x ! function argument
real(kind=dp) :: out
out = log_gamma(x)
end function r_lgamma

pure elemental function r_beta(a, b) result(out)
! Elementwise beta function.
real(kind=dp), intent(in) :: a ! first shape parameter
real(kind=dp), intent(in) :: b ! second shape parameter
real(kind=dp) :: out
out = exp(log_gamma(a) + log_gamma(b) - log_gamma(a + b))
end function r_beta

pure elemental function r_lbeta(a, b) result(out)
! Elementwise log beta function.
real(kind=dp), intent(in) :: a ! first shape parameter
real(kind=dp), intent(in) :: b ! second shape parameter
real(kind=dp) :: out
out = log_gamma(a) + log_gamma(b) - log_gamma(a + b)
end function r_lbeta

pure elemental function r_choose(n, k) result(out)
! Elementwise binomial coefficient function.
real(kind=dp), intent(in) :: n ! number of items
real(kind=dp), intent(in) :: k ! number chosen
real(kind=dp) :: out
integer :: i, kk
real(kind=dp) :: fraction
if (k < 0.0_dp) then
   out = 0.0_dp
else
   kk = floor(k)
   fraction = k - real(kk, kind=dp)
   if (fraction > 0.5_dp .or. (fraction == 0.5_dp .and. mod(kk, 2) /= 0)) kk = kk + 1
   out = 1.0_dp
   do i = 1, kk
      out = out * (n - real(i - 1, kind=dp)) / real(i, kind=dp)
   end do
end if
end function r_choose

pure elemental function r_lchoose(n, k) result(out)
! Elementwise log binomial coefficient.
real(kind=dp), intent(in) :: n ! number of items
real(kind=dp), intent(in) :: k ! number chosen
real(kind=dp) :: out
integer :: i, kk
real(kind=dp) :: factor, fraction
if (k < 0.0_dp) then
   out = -ieee_value(0.0_dp, ieee_positive_inf)
else
   kk = floor(k)
   fraction = k - real(kk, kind=dp)
   if (fraction > 0.5_dp .or. (fraction == 0.5_dp .and. mod(kk, 2) /= 0)) kk = kk + 1
   out = 0.0_dp
   do i = 1, kk
      factor = abs(n - real(i - 1, kind=dp)) / real(i, kind=dp)
      if (factor == 0.0_dp) then
         out = -ieee_value(0.0_dp, ieee_positive_inf)
         return
      end if
      out = out + log(factor)
   end do
end if
end function r_lchoose

pure elemental function r_factorial(x) result(out)
! Elementwise factorial using gamma(x + 1) for compatibility.
real(kind=dp), intent(in) :: x ! function argument
real(kind=dp) :: out
if (x < 0.0_dp .and. x == floor(x)) then
   out = ieee_value(0.0_dp, ieee_quiet_nan)
else
   out = r_gamma(x + 1.0_dp)
end if
end function r_factorial

pure elemental function r_lfactorial(x) result(out)
! Elementwise log factorial.
real(kind=dp), intent(in) :: x ! function argument
real(kind=dp) :: out
out = r_lgamma(x + 1.0_dp)
end function r_lfactorial

pure elemental integer function r_digit_power_sum(x, power) result(out)
! Sum the decimal digits of x raised to a fixed positive power.
integer, intent(in) :: x, power
integer :: work, digit
out = 0
work = abs(x)
do
   digit = mod(work, 10)
   out = out + digit**power
   work = work / 10
   if (work == 0) exit
end do
end function r_digit_power_sum

recursive pure elemental function r_digamma(x) result(out)
! Elementwise digamma using recurrence, reflection, and an asymptotic series.
real(kind=dp), intent(in) :: x ! function argument
real(kind=dp) :: out, inv, inv2, y
real(kind=dp), parameter :: pi = acos(-1.0_dp)
if (.not. ieee_is_finite(x)) then
   if (x > 0.0_dp) then
      out = x
   else
      out = ieee_value(0.0_dp, ieee_quiet_nan)
   end if
   return
end if
if (x <= 0.0_dp) then
   if (x == anint(x)) then
      out = ieee_value(0.0_dp, ieee_quiet_nan)
   else
      out = r_digamma(1.0_dp - x) - pi / tan(pi * x)
   end if
   return
end if
out = 0.0_dp
y = x
do while (y < 8.0_dp)
   out = out - 1.0_dp / y
   y = y + 1.0_dp
end do
inv = 1.0_dp / y
inv2 = inv * inv
out = out + log(y) - 0.5_dp * inv - inv2 * (1.0_dp / 12.0_dp - &
   inv2 * (1.0_dp / 120.0_dp - inv2 * (1.0_dp / 252.0_dp - &
   inv2 * (1.0_dp / 240.0_dp - inv2 * (1.0_dp / 132.0_dp)))))
end function r_digamma

recursive pure elemental function r_trigamma(x) result(out)
! Elementwise trigamma using recurrence, reflection, and an asymptotic series.
real(kind=dp), intent(in) :: x ! function argument
real(kind=dp) :: out, inv, inv2, y
real(kind=dp), parameter :: pi = acos(-1.0_dp)
if (.not. ieee_is_finite(x)) then
   if (x > 0.0_dp) then
      out = 0.0_dp
   else
      out = ieee_value(0.0_dp, ieee_quiet_nan)
   end if
   return
end if
if (x <= 0.0_dp) then
   if (x == anint(x)) then
      out = ieee_value(0.0_dp, ieee_quiet_nan)
   else
      out = (pi / sin(pi * x))**2 - r_trigamma(1.0_dp - x)
   end if
   return
end if
out = 0.0_dp
y = x
do while (y < 8.0_dp)
   out = out + 1.0_dp / (y * y)
   y = y + 1.0_dp
end do
inv = 1.0_dp / y
inv2 = inv * inv
out = out + inv + 0.5_dp * inv2 + inv * inv2 * (1.0_dp / 6.0_dp - &
   inv2 * (1.0_dp / 30.0_dp - inv2 * (1.0_dp / 42.0_dp - &
   inv2 * (1.0_dp / 30.0_dp - inv2 * (5.0_dp / 66.0_dp)))))
end function r_trigamma

pure elemental function r_psigamma(x, deriv) result(out)
! Elementwise poly-gamma approximation.
! deriv=0 returns digamma; positive orders use a convergent Hurwitz-zeta series.
real(kind=dp), intent(in) :: x ! function argument
integer, intent(in) :: deriv ! derivative order
real(kind=dp) :: out
integer :: i, p
real(kind=dp) :: factorial, series, sign, y
if (deriv < 0) then
   out = ieee_value(0.0_dp, ieee_quiet_nan)
elseif (deriv == 0) then
   out = r_digamma(x)
elseif (deriv == 1) then
   out = r_trigamma(x)
elseif (x <= 0.0_dp) then
   out = ieee_value(0.0_dp, ieee_quiet_nan)
else
   p = deriv + 1
   factorial = 1.0_dp
   do i = 2, deriv
      factorial = factorial * real(i, kind=dp)
   end do
   series = 0.0_dp
   do i = 0, 31
      series = series + 1.0_dp / (x + real(i, kind=dp))**p
   end do
   y = x + 32.0_dp
   series = series + y**(1 - p) / real(p - 1, kind=dp) + 0.5_dp * y**(-p) + &
      real(p, kind=dp) * y**(-p - 1) / 12.0_dp - &
      real(p * (p + 1) * (p + 2), kind=dp) * y**(-p - 3) / 720.0_dp
   sign = merge(1.0_dp, -1.0_dp, mod(deriv, 2) == 1)
   out = sign * factorial * series
end if
end function r_psigamma

pure function colMeans(x) result(out)
! Column means of a numeric matrix.
real(kind=dp), intent(in) :: x(:,:) ! observations by rows, variables by columns
real(kind=dp), allocatable :: out(:)
integer :: n
n = size(x, 1)
if (size(x, 2) <= 0) then
   allocate(out(0))
   return
end if
if (n <= 0) then
   allocate(out(size(x, 2)))
   out = ieee_value(0.0_dp, ieee_quiet_nan)
   return
end if
out = sum(x, dim=1) / real(n, kind=dp)
end function colMeans

pure function cov_vec(x, y) result(out)
! Sample covariance of two vectors (n-1 denominator).
real(kind=dp), intent(in) :: x(:) ! first sample values
real(kind=dp), intent(in) :: y(:) ! second sample values
real(kind=dp) :: out, mx, my
integer :: n
n = min(size(x), size(y))
if (n <= 1) then
   out = ieee_value(0.0_dp, ieee_quiet_nan)
   return
end if
mx = sum(x(1:n)) / real(n, kind=dp)
my = sum(y(1:n)) / real(n, kind=dp)
out = sum((x(1:n) - mx) * (y(1:n) - my)) / real(n - 1, kind=dp)
end function cov_vec

pure function cov_mat(x) result(out)
! Sample covariance matrix of columns of x (n-1 denominator).
real(kind=dp), intent(in) :: x(:,:) ! observations by rows, variables by columns
real(kind=dp), allocatable :: out(:,:)
real(kind=dp), allocatable :: xc(:,:), mu(:)
integer :: n, p
n = size(x, 1)
p = size(x, 2)
if (p <= 0) then
   allocate(out(0,0))
   return
end if
if (n <= 1) then
   allocate(out(p, p))
   out = ieee_value(0.0_dp, ieee_quiet_nan)
   return
end if
mu = sum(x, dim=1) / real(n, kind=dp)
xc = x - spread(mu, dim=1, ncopies=n)
out = matmul(transpose(xc), xc) / real(n - 1, kind=dp)
end function cov_mat

pure function cov2cor(x) result(out)
! Convert a covariance matrix to a correlation matrix.
real(kind=dp), intent(in) :: x(:,:) ! covariance matrix
real(kind=dp), allocatable :: out(:,:)
real(kind=dp), allocatable :: d(:)
integer :: i, j, n
n = min(size(x, 1), size(x, 2))
allocate(out(size(x, 1), size(x, 2)))
if (size(x, 1) == 0 .or. size(x, 2) == 0) return
allocate(d(n))
do i = 1, n
   if (x(i, i) > 0.0_dp) then
      d(i) = sqrt(x(i, i))
   else
      d(i) = ieee_value(0.0_dp, ieee_quiet_nan)
   end if
end do
do j = 1, size(x, 2)
   do i = 1, size(x, 1)
      if (i <= n .and. j <= n .and. d(i) > 0.0_dp .and. d(j) > 0.0_dp) then
         out(i, j) = x(i, j) / (d(i) * d(j))
      else
         out(i, j) = ieee_value(0.0_dp, ieee_quiet_nan)
      end if
   end do
end do
end function cov2cor

pure function scale_vec(x, center, scale) result(out)
! R scale() on a vector returns an n-by-1 matrix.
! Defaults: center=.true., scale=.true.
real(kind=dp), intent(in) :: x(:) ! values to scale as one column
logical, intent(in), optional :: center ! subtract column mean flag
logical, intent(in), optional :: scale ! divide by column scale flag
real(kind=dp), allocatable :: out(:,:)
real(kind=dp), allocatable :: xmat(:,:)
integer :: n
n = size(x)
allocate(xmat(n, 1))
if (n > 0) xmat(:, 1) = x
out = scale_mat(xmat, center=center, scale=scale)
end function scale_vec

pure function scale_mat(x, center, scale) result(out)
! Column-center and optionally column-scale a matrix, matching common R scale() use.
! Defaults: center=.true., scale=.true.
real(kind=dp), intent(in) :: x(:,:) ! observations by rows, variables by columns
logical, intent(in), optional :: center ! subtract column means flag
logical, intent(in), optional :: scale ! divide by column scales flag
real(kind=dp), allocatable :: out(:,:)
real(kind=dp), allocatable :: mu(:), sig(:)
integer :: n, p, j
logical :: do_center, do_scale
n = size(x, 1)
p = size(x, 2)
do_center = .true.
do_scale = .true.
if (present(center)) do_center = center
if (present(scale)) do_scale = scale
allocate(out(n, p))
out = x
if (p <= 0 .or. n <= 0) return
if (do_center) then
   allocate(mu(p))
   mu = sum(x, dim=1) / real(n, kind=dp)
   out = out - spread(mu, dim=1, ncopies=n)
end if
if (do_scale) then
   allocate(sig(p))
   do j = 1, p
      if (n > 1) then
         sig(j) = sqrt(sum(out(:, j)**2) / real(n - 1, kind=dp))
      else
         sig(j) = 0.0_dp
      end if
      if (sig(j) > 0.0_dp) then
         out(:, j) = out(:, j) / sig(j)
      else
         out(:, j) = ieee_value(0.0_dp, ieee_quiet_nan)
      end if
   end do
end if
end function scale_mat

pure function all_equal_real_scalar(a, b, tolerance) result(out)
! Compare real scalar inputs with R all.equal-style semantics.
! Default tolerance is sqrt(epsilon(1.0_dp)).
real(kind=dp), intent(in) :: a ! first value
real(kind=dp), intent(in) :: b ! second value
real(kind=dp), intent(in), optional :: tolerance ! relative comparison tolerance
logical :: out
real(kind=dp) :: tol, scale_ab
tol = sqrt(epsilon(1.0_dp))
if (present(tolerance)) tol = tolerance
if (a /= a .or. b /= b) then
   out = (a /= a) .and. (b /= b)
else if (.not. ieee_is_finite(a) .or. .not. ieee_is_finite(b)) then
   out = a == b
else
   scale_ab = max(1.0_dp, abs(a), abs(b))
   out = abs(a - b) <= tol * scale_ab
end if
end function all_equal_real_scalar

pure function all_equal_real_vec(a, b, tolerance) result(out)
! Compare real vec inputs with R all.equal-style semantics.
! Default tolerance is sqrt(epsilon(1.0_dp)).
real(kind=dp), intent(in) :: a(:) ! first values
real(kind=dp), intent(in) :: b(:) ! second values
real(kind=dp), intent(in), optional :: tolerance ! relative comparison tolerance
logical :: out
real(kind=dp) :: tol
integer :: i
tol = sqrt(epsilon(1.0_dp))
if (present(tolerance)) tol = tolerance
if (size(a) /= size(b)) then
   out = .false.
else
   out = .true.
   do i = 1, size(a)
      if (.not. all_equal_real_scalar(a(i), b(i), tolerance=tol)) then
         out = .false.
         exit
      end if
   end do
end if
end function all_equal_real_vec

pure function all_equal_real_mat(a, b, tolerance) result(out)
! Compare real mat inputs with R all.equal-style semantics.
! Default tolerance is sqrt(epsilon(1.0_dp)).
real(kind=dp), intent(in) :: a(:,:) ! first matrix
real(kind=dp), intent(in) :: b(:,:) ! second matrix
real(kind=dp), intent(in), optional :: tolerance ! relative comparison tolerance
logical :: out
real(kind=dp) :: tol
integer :: i, j
tol = sqrt(epsilon(1.0_dp))
if (present(tolerance)) tol = tolerance
if (any(shape(a) /= shape(b))) then
   out = .false.
else
   out = .true.
   do j = 1, size(a, 2)
      do i = 1, size(a, 1)
         if (.not. all_equal_real_scalar(a(i, j), b(i, j), tolerance=tol)) then
            out = .false.
            return
         end if
      end do
   end do
end if
end function all_equal_real_mat

pure function all_equal_int_scalar(a, b, tolerance) result(out)
! Compare int scalar inputs with R all.equal-style semantics.
integer, intent(in) :: a ! first value
integer, intent(in) :: b ! second value
real(kind=dp), intent(in), optional :: tolerance ! relative comparison tolerance
logical :: out
out = all_equal_real_scalar(real(a, kind=dp), real(b, kind=dp), tolerance=tolerance)
end function all_equal_int_scalar

pure function all_equal_int_vec(a, b, tolerance) result(out)
! Compare int vec inputs with R all.equal-style semantics.
integer, intent(in) :: a(:) ! first values
integer, intent(in) :: b(:) ! second values
real(kind=dp), intent(in), optional :: tolerance ! relative comparison tolerance
logical :: out
out = all_equal_real_vec(real(a, kind=dp), real(b, kind=dp), tolerance=tolerance)
end function all_equal_int_vec

pure function all_equal_int_mat(a, b, tolerance) result(out)
! Compare int mat inputs with R all.equal-style semantics.
integer, intent(in) :: a(:,:) ! first matrix
integer, intent(in) :: b(:,:) ! second matrix
real(kind=dp), intent(in), optional :: tolerance ! relative comparison tolerance
logical :: out
out = all_equal_real_mat(real(a, kind=dp), real(b, kind=dp), tolerance=tolerance)
end function all_equal_int_mat

pure function all_equal_complex_scalar(a, b, tolerance) result(out)
! Compare complex scalars componentwise with R all.equal-style tolerance.
complex(kind=dp), intent(in) :: a, b
real(kind=dp), intent(in), optional :: tolerance
logical :: out
out = all_equal_real_scalar(real(a, kind=dp), real(b, kind=dp), tolerance=tolerance) .and. &
   all_equal_real_scalar(aimag(a), aimag(b), tolerance=tolerance)
end function all_equal_complex_scalar

pure function all_equal_complex_vec(a, b, tolerance) result(out)
! Compare complex vectors componentwise with R all.equal-style tolerance.
complex(kind=dp), intent(in) :: a(:), b(:)
real(kind=dp), intent(in), optional :: tolerance
logical :: out
integer :: i
if (size(a) /= size(b)) then
   out = .false.
   return
end if
out = .true.
do i = 1, size(a)
   if (.not. all_equal_complex_scalar(a(i), b(i), tolerance=tolerance)) then
      out = .false.
      exit
   end if
end do
end function all_equal_complex_vec

pure function all_equal_complex_mat(a, b, tolerance) result(out)
! Compare complex matrices componentwise with R all.equal-style tolerance.
complex(kind=dp), intent(in) :: a(:,:), b(:,:)
real(kind=dp), intent(in), optional :: tolerance
logical :: out
integer :: i, j
if (any(shape(a) /= shape(b))) then
   out = .false.
   return
end if
out = .true.
do j = 1, size(a, 2)
   do i = 1, size(a, 1)
      if (.not. all_equal_complex_scalar(a(i, j), b(i, j), tolerance=tolerance)) then
         out = .false.
         return
      end if
   end do
end do
end function all_equal_complex_mat

pure function all_equal_logical_scalar(a, b, tolerance) result(out)
! Compare logical scalars; tolerance is accepted for API compatibility.
logical, intent(in) :: a, b
real(kind=dp), intent(in), optional :: tolerance
logical :: out
if (present(tolerance)) continue
out = a .eqv. b
end function all_equal_logical_scalar

pure function all_equal_logical_vec(a, b, tolerance) result(out)
! Compare logical vec inputs with R all.equal-style semantics.
! tolerance is accepted for API compatibility and ignored.
logical, intent(in) :: a(:) ! first values
logical, intent(in) :: b(:) ! second values
real(kind=dp), intent(in), optional :: tolerance ! accepted for compatibility; ignored
logical :: out
if (present(tolerance)) continue
out = size(a) == size(b)
if (out) out = all(a .eqv. b)
end function all_equal_logical_vec

pure function all_equal_logical_mat(a, b, tolerance) result(out)
! Compare logical mat inputs with R all.equal-style semantics.
! tolerance is accepted for API compatibility and ignored.
logical, intent(in) :: a(:,:) ! first matrix
logical, intent(in) :: b(:,:) ! second matrix
real(kind=dp), intent(in), optional :: tolerance ! accepted for compatibility; ignored
logical :: out
if (present(tolerance)) continue
out = all(shape(a) == shape(b))
if (out) out = all(a .eqv. b)
end function all_equal_logical_mat

pure function all_equal_char_scalar(a, b, tolerance) result(out)
character(len=*), intent(in) :: a, b
real(kind=dp), intent(in), optional :: tolerance
logical :: out
if (present(tolerance)) continue
out = a == b
end function all_equal_char_scalar

pure function all_equal_char_vec(a, b, tolerance) result(out)
character(len=*), intent(in) :: a(:), b(:)
real(kind=dp), intent(in), optional :: tolerance
logical :: out
if (present(tolerance)) continue
out = size(a) == size(b)
if (out) out = all(a == b)
end function all_equal_char_vec

pure function all_equal_char_mat(a, b, tolerance) result(out)
character(len=*), intent(in) :: a(:,:), b(:,:)
real(kind=dp), intent(in), optional :: tolerance
logical :: out
if (present(tolerance)) continue
out = all(shape(a) == shape(b))
if (out) out = all(a == b)
end function all_equal_char_mat

elemental function r_log_scalar(x) result(out)
! Apply R-like log handling to scalar input.
real(kind=dp), intent(in) :: x ! value whose natural logarithm is requested
real(kind=dp) :: out
if (r_is_na_payload(x)) then
   out = r_na_real()
else if (x /= x) then
   out = x
else if (x < 0.0_dp) then
   out = ieee_value(0.0_dp, ieee_quiet_nan)
else
   out = log(x)
end if
end function r_log_scalar

pure function r_log_vec(x) result(out)
! Apply R-like log handling to vec input.
real(kind=dp), intent(in) :: x(:) ! values whose natural logarithms are requested
real(kind=dp), allocatable :: out(:)
integer :: i
allocate(out(size(x)))
do i = 1, size(x)
   out(i) = r_log_scalar(x(i))
end do
end function r_log_vec

pure function r_log_mat(x) result(out)
! Apply R-like log handling to mat input.
real(kind=dp), intent(in) :: x(:,:) ! matrix whose natural logarithms are requested
real(kind=dp), allocatable :: out(:,:)
integer :: i, j
allocate(out(size(x, 1), size(x, 2)))
do j = 1, size(x, 2)
   do i = 1, size(x, 1)
      out(i, j) = r_log_scalar(x(i, j))
   end do
end do
end function r_log_mat

pure function cor_vec(x, y) result(out)
! Sample correlation of two vectors.
real(kind=dp), intent(in) :: x(:) ! first sample values
real(kind=dp), intent(in) :: y(:) ! second sample values
real(kind=dp) :: out, sdx, sdy
integer :: n
n = min(size(x), size(y))
if (n <= 1) then
   out = ieee_value(0.0_dp, ieee_quiet_nan)
   return
end if
sdx = sd(x(1:n))
sdy = sd(y(1:n))
if (.not. ieee_is_finite(sdx) .or. .not. ieee_is_finite(sdy) .or. sdx <= 0.0_dp .or. sdy <= 0.0_dp) then
   out = ieee_value(0.0_dp, ieee_quiet_nan)
   return
end if
out = cov_vec(x(1:n), y(1:n)) / (sdx * sdy)
end function cor_vec

pure function cor_mat(x) result(out)
! Sample correlation matrix of columns of x.
real(kind=dp), intent(in) :: x(:,:) ! observations by rows, variables by columns
real(kind=dp), allocatable :: out(:,:), s(:)
real(kind=dp), allocatable :: c(:,:)
integer :: i, j, p
c = cov_mat(x)
p = size(c, 1)
if (p <= 0) then
   allocate(out(0,0))
   return
end if
allocate(out(p, p), s(p))
do i = 1, p
   s(i) = sqrt(c(i, i))
end do
do i = 1, p
   do j = 1, p
      if (.not. ieee_is_finite(s(i)) .or. .not. ieee_is_finite(s(j)) .or. s(i) <= 0.0_dp .or. s(j) <= 0.0_dp) then
         out(i, j) = ieee_value(0.0_dp, ieee_quiet_nan)
      else
         out(i, j) = c(i, j) / (s(i) * s(j))
      end if
   end do
end do
end function cor_mat

pure function cor_mat_pair(x, y) result(out)
! Sample column-pair correlation matrix between x and y.
real(kind=dp), intent(in) :: x(:,:) ! first observation matrix
real(kind=dp), intent(in) :: y(:,:) ! second observation matrix
real(kind=dp), allocatable :: out(:,:)
integer :: i, j, n, px, py
n = min(size(x, 1), size(y, 1))
px = size(x, 2)
py = size(y, 2)
allocate(out(px, py))
if (n <= 1) then
   out = ieee_value(0.0_dp, ieee_quiet_nan)
   return
end if
do i = 1, px
   do j = 1, py
      out(i, j) = cor_vec(x(1:n, i), y(1:n, j))
   end do
end do
end function cor_mat_pair

pure integer function count_ws_tokens(line) result(n_tok)
! Count whitespace-separated tokens in one text line.
character(len=*), intent(in) :: line ! text line to scan
integer :: i, n
logical :: in_tok
n = len_trim(line)
n_tok = 0
in_tok = .false.
do i = 1, n
   if (line(i:i) /= " " .and. line(i:i) /= char(9)) then
      if (.not. in_tok) then
         n_tok = n_tok + 1
         in_tok = .true.
      end if
   else
      in_tok = .false.
   end if
end do
end function count_ws_tokens

subroutine read_real_vector(file_path, x)
! Read whitespace-delimited real values into a vector.
character(len=*), intent(in) :: file_path ! input file path
real(kind=dp), allocatable, intent(out) :: x(:) ! values read from file
integer :: fp, ios, n, cap, new_cap, n_tok
character(len=4096) :: line
real(kind=dp), allocatable :: vals(:)
n = 0
cap = 0
open(newunit=fp, file=file_path, status="old", action="read")
do
   read(fp, "(a)", iostat=ios) line
   if (ios /= 0) exit
   n_tok = count_ws_tokens(line)
   if (n_tok <= 0) cycle
   allocate(vals(n_tok))
   read(line, *, iostat=ios) vals
   if (ios /= 0) then
      deallocate(vals)
      cycle
   end if
   if (n + n_tok > cap) then
      new_cap = max(n + n_tok, merge(1024, 2 * cap, cap == 0))
      block
         real(kind=dp), allocatable :: tmp(:)
         allocate(tmp(new_cap))
         if (allocated(x) .and. n > 0) tmp(1:n) = x(1:n)
         call move_alloc(tmp, x)
      end block
      cap = new_cap
   end if
   x(n + 1:n + n_tok) = vals
   n = n + n_tok
   deallocate(vals)
end do
close(fp)
if (n == 0) then
   allocate(x(0))
else if (n < size(x)) then
   x = x(1:n)
end if
end subroutine read_real_vector

function scan_real(file_path) result(x)
! Expression-form numeric scan(file_path) helper.
character(len=*), intent(in) :: file_path
real(kind=dp), allocatable :: x(:)
call read_real_vector(file_path, x)
end function scan_real

subroutine read_table_real_matrix(file_path, x, header)
! Read a whitespace-delimited numeric table into a matrix.
! If header is present and true, skip the first nonblank line.
character(len=*), intent(in) :: file_path ! input file path
real(kind=dp), allocatable, intent(out) :: x(:,:) ! numeric table values
logical, intent(in), optional :: header ! header row flag
integer :: fp, ios, nrow, ncol, i
character(len=4096) :: line
logical :: has_header, skipped_header
nrow = 0
ncol = 0
has_header = .false.
if (present(header)) has_header = header
skipped_header = .false.
open(newunit=fp, file=file_path, status="old", action="read")
do
   read(fp, "(A)", iostat=ios) line
   if (ios /= 0) exit
   if (len_trim(line) == 0) cycle
   if (ncol == 0) ncol = count_ws_tokens(line)
   if (has_header .and. .not. skipped_header) then
      skipped_header = .true.
      cycle
   end if
   nrow = nrow + 1
end do
if (nrow <= 0 .or. ncol <= 0) then
   allocate(x(0,0))
   close(fp)
   return
end if
allocate(x(nrow, ncol))
rewind(fp)
i = 0
skipped_header = .false.
do
   read(fp, "(A)", iostat=ios) line
   if (ios /= 0) exit
   if (len_trim(line) == 0) cycle
   if (has_header .and. .not. skipped_header) then
      skipped_header = .true.
      cycle
   end if
   i = i + 1
   read(line, *) x(i, 1:ncol)
end do
close(fp)
end subroutine read_table_real_matrix

subroutine read_csv_real_matrix(file_path, x, max_rows, max_cols)
! Read a comma-delimited numeric CSV with one header row into a matrix.
character(len=*), intent(in) :: file_path ! input CSV file path
real(kind=dp), allocatable, intent(out) :: x(:,:) ! numeric data rows
integer, intent(in), optional :: max_rows ! maximum number of data rows to read
integer, intent(in), optional :: max_cols ! maximum number of leading columns to read
character(len=4096) :: line
character(len=256) :: token
integer :: fp, ios, nrow, ncol, i, j, k, start, stop, yy, mm, dd
integer :: row_limit, col_limit

row_limit = huge(1)
col_limit = huge(1)
if (present(max_rows)) then
   if (max_rows < 0) error stop "read_csv_real_matrix: max_rows must be nonnegative"
   row_limit = max_rows
end if
if (present(max_cols)) then
   if (max_cols < 0) error stop "read_csv_real_matrix: max_cols must be nonnegative"
   col_limit = max_cols
end if
nrow = 0
ncol = 0
open(newunit=fp, file=file_path, status="old", action="read")
read(fp, "(A)", iostat=ios) line
if (ios /= 0) then
   allocate(x(0,0))
   close(fp)
   return
end if
ncol = 1
do k = 1, len_trim(line)
   if (line(k:k) == ",") ncol = ncol + 1
end do
ncol = min(ncol, col_limit)
do
   if (nrow >= row_limit) exit
   read(fp, "(A)", iostat=ios) line
   if (ios /= 0) exit
   if (len_trim(line) == 0) cycle
   nrow = nrow + 1
end do
if (nrow <= 0 .or. ncol <= 0) then
   allocate(x(nrow, ncol))
   close(fp)
   return
end if
allocate(x(nrow, ncol))
rewind(fp)
read(fp, "(A)", iostat=ios) line
i = 0
do
   read(fp, "(A)", iostat=ios) line
   if (ios /= 0) exit
   if (len_trim(line) == 0) cycle
   if (i >= nrow) exit
   i = i + 1
   start = 1
   do j = 1, ncol
      stop = len_trim(line) + 1
      do k = start, len_trim(line)
         if (line(k:k) == ",") then
            stop = k
            exit
         end if
      end do
      token = adjustl(line(start:stop - 1))
      read(token, *, iostat=ios) x(i, j)
      if (ios /= 0) then
         if (len_trim(token) >= 10 .and. token(5:5) == "-" .and. token(8:8) == "-") then
            read(token(1:4), *, iostat=ios) yy
            if (ios == 0) read(token(6:7), *, iostat=ios) mm
            if (ios == 0) read(token(9:10), *, iostat=ios) dd
            if (ios == 0) then
               x(i, j) = real(yy * 10000 + mm * 100 + dd, kind=dp)
            else
               x(i, j) = ieee_value(0.0_dp, ieee_quiet_nan)
            end if
         else
            x(i, j) = ieee_value(0.0_dp, ieee_quiet_nan)
         end if
      end if
      start = stop + 1
   end do
end do
close(fp)
end subroutine read_csv_real_matrix

function read_csv_text_column(file_path, column, max_rows) result(values)
! Read one CSV column as text, excluding the header row.
character(len=*), intent(in) :: file_path
integer, intent(in) :: column
integer, intent(in), optional :: max_rows
character(len=:), allocatable :: values(:)
character(len=4096) :: line
integer :: fp, ios, nrow, i, row_limit

if (column < 1) error stop "read_csv_text_column: column must be positive"
row_limit = huge(1)
if (present(max_rows)) then
   if (max_rows < 0) error stop "read_csv_text_column: max_rows must be nonnegative"
   row_limit = max_rows
end if
open(newunit=fp, file=file_path, status="old", action="read")
read(fp, "(A)", iostat=ios) line
nrow = 0
do
   if (nrow >= row_limit) exit
   read(fp, "(A)", iostat=ios) line
   if (ios /= 0) exit
   if (len_trim(line) == 0) cycle
   nrow = nrow + 1
end do
allocate(character(len=4096) :: values(nrow))
rewind(fp)
read(fp, "(A)", iostat=ios) line
i = 0
do while (i < nrow)
   read(fp, "(A)", iostat=ios) line
   if (ios /= 0) exit
   if (len_trim(line) == 0) cycle
   i = i + 1
   values(i) = csv_text_field(line, column)
end do
close(fp)
end function read_csv_text_column

function csv_text_field(line, column) result(value)
! Extract one field from a CSV record, respecting simple quoted fields.
character(len=*), intent(in) :: line
integer, intent(in) :: column
character(len=:), allocatable :: value
character(len=1) :: ch, quote
integer :: current, i, start, line_len, value_len
logical :: in_quotes

value = ""
current = 1
start = 1
line_len = len_trim(line)
in_quotes = .false.
quote = " "
do i = 1, line_len + 1
   if (i <= line_len) then
      ch = line(i:i)
   else
      ch = ","
   end if
   if (in_quotes) then
      if (ch == quote) in_quotes = .false.
   else if (ch == achar(34) .or. ch == achar(39)) then
      in_quotes = .true.
      quote = ch
   else if (ch == ",") then
      if (current == column) then
         if (i > start) value = trim(adjustl(line(start:i - 1)))
         value_len = len(value)
         if (value_len >= 2) then
            if ((value(1:1) == achar(34) .and. value(value_len:value_len) == achar(34)) .or. &
                (value(1:1) == achar(39) .and. value(value_len:value_len) == achar(39))) then
               if (value_len > 2) then
                  value = value(2:value_len - 1)
               else
                  value = ""
               end if
            end if
         end if
         return
      end if
      current = current + 1
      start = i + 1
   end if
end do
error stop "csv_text_field: column not found"
end function csv_text_field

function read_csv_header_names(file_path) result(names)
! Read the first CSV row as a character vector of header names.
character(len=*), intent(in) :: file_path ! input CSV file path
character(len=:), allocatable :: names(:)
character(len=4096) :: line
integer :: fp, ios, ncol, i, start, pos
open(newunit=fp, file=file_path, status="old", action="read")
read(fp, "(A)", iostat=ios) line
close(fp)
if (ios /= 0 .or. len_trim(line) == 0) then
   allocate(character(len=1) :: names(0))
   return
end if
ncol = 1
do i = 1, len_trim(line)
   if (line(i:i) == ",") ncol = ncol + 1
end do
allocate(character(len=128) :: names(ncol))
start = 1
pos = 1
do i = 1, len_trim(line) + 1
   if (i > len_trim(line) .or. line(i:i) == ",") then
      if (i > start) then
         names(pos) = adjustl(line(start:i-1))
         names(pos) = trim(names(pos))
         if (len_trim(names(pos)) >= 2) then
            if ((names(pos)(1:1) == '"' .and. names(pos)(len_trim(names(pos)):len_trim(names(pos))) == '"') .or. &
                (names(pos)(1:1) == "'" .and. names(pos)(len_trim(names(pos)):len_trim(names(pos))) == "'")) then
               names(pos) = names(pos)(2:len_trim(names(pos))-1)
            end if
         end if
      else
         names(pos) = ""
      end if
      pos = pos + 1
      start = i + 1
   end if
end do
end function read_csv_header_names

subroutine write_table_real_matrix(file_path, x, names)
! Write a numeric matrix as a whitespace-delimited table.
character(len=*), intent(in) :: file_path ! output file path
real(kind=dp), intent(in) :: x(:,:) ! matrix to write
character(len=*), intent(in), optional :: names(:) ! optional column names
integer :: fp, i, j
open(newunit=fp, file=file_path, status="replace", action="write")
if (present(names)) then
   do j = 1, size(x, 2)
      if (j > 1) write(fp, '(1x)', advance='no')
      if (j <= size(names)) then
         write(fp, '(a)', advance='no') trim(names(j))
      else
         write(fp, '("V",i0)', advance='no') j
      end if
   end do
   write(fp,*)
end if
do i = 1, size(x, 1)
   write(fp, *) x(i, 1:size(x, 2))
end do
close(fp)
end subroutine write_table_real_matrix

subroutine write_table_real_vector(file_path, x, name)
! Write a numeric vector as a one-column whitespace-delimited table.
character(len=*), intent(in) :: file_path ! output file path
real(kind=dp), intent(in) :: x(:) ! vector to write
character(len=*), intent(in), optional :: name ! optional column name
integer :: fp, i
open(newunit=fp, file=file_path, status="replace", action="write")
if (present(name)) write(fp, '(a)') trim(name)
do i = 1, size(x)
   write(fp, *) x(i)
end do
close(fp)
end subroutine write_table_real_vector

subroutine print_matrix_real(x, int_like)
! Print a real matrix row-by-row; use integer format when all values are integer-like.
real(kind=dp), intent(in) :: x(:,:) ! matrix to print
logical, intent(in), optional :: int_like ! force integer-like formatting
logical :: use_int_like, all_int
integer :: i, j
integer(kind=int64) :: k
real(kind=dp) :: r, tol
use_int_like = print_int_like_default
if (present(int_like)) use_int_like = int_like
all_int = .false.
if (use_int_like) then
   all_int = .true.
   do i = 1, size(x, 1)
      do j = 1, size(x, 2)
         r = x(i, j)
         if (.not. ieee_is_finite(r)) then
            all_int = .false.
            exit
         end if
         if (abs(r) > real(huge(0_int64), kind=dp)) then
            all_int = .false.
            exit
         end if
         k = nint(r, kind=int64)
         tol = print_int_like_tol * max(1.0_dp, abs(r))
         if (abs(r - real(k, kind=dp)) > tol) then
            all_int = .false.
            exit
         end if
      end do
      if (.not. all_int) exit
   end do
end if
if (all_int) then
   do i = 1, size(x, 1)
      do j = 1, size(x, 2)
         k = nint(x(i, j), kind=int64)
         write(*,"(i0)", advance="no") k
         if (j < size(x, 2)) write(*,"(a)", advance="no") " "
      end do
      write(*,*)
   end do
else
   do i = 1, size(x, 1)
      write(*,"(*(g0,1x))") x(i, :)
   end do
end if
end subroutine print_matrix_real

subroutine print_matrix_rstyle_real(x, digits)
! Print a real matrix with R-like column and row labels.
real(kind=dp), intent(in) :: x(:,:) ! matrix to print
integer, intent(in), optional :: digits ! digits after decimal point for display
integer :: i, j
integer(kind=int64) :: k
logical :: all_int
real(kind=dp) :: r, tol
character(len=12) :: col_label
character(len=32) :: fmt
all_int = .true.
if (present(digits)) then
   all_int = .false.
   write(fmt, '("(f12.",i0,",1x)")') max(0, digits)
else
   do i = 1, size(x, 1)
      do j = 1, size(x, 2)
         r = x(i, j)
         if (.not. ieee_is_finite(r)) then
            all_int = .false.
            exit
         end if
         if (abs(r) > real(huge(0_int64), kind=dp)) then
            all_int = .false.
            exit
         end if
         k = nint(r, kind=int64)
         tol = print_int_like_tol * max(1.0_dp, abs(r))
         if (abs(r - real(k, kind=dp)) > tol) then
            all_int = .false.
            exit
         end if
      end do
      if (.not. all_int) exit
   end do
end if
write(*,'(5x)', advance='no')
do i = 1, size(x, 2)
   write(col_label, '("[,",i0,"]")') i
   write(*,'(a12,1x)', advance='no') adjustr(col_label)
end do
write(*,*)
do i = 1, size(x, 1)
   write(*,'("[",i0,",]",1x)', advance='no') i
   do j = 1, size(x, 2)
      if (all_int) then
         write(*,'(i12,1x)', advance='no') nint(x(i, j), kind=int64)
      else if (.not. ieee_is_finite(x(i, j))) then
         write(*,'(a12,1x)', advance='no') "NA"
      else if (present(digits)) then
         write(*,fmt, advance='no') x(i, j)
      else if (x(i, j) == 0.0_dp .or. (abs(x(i, j)) >= 1.0e-4_dp .and. abs(x(i, j)) < 1.0e6_dp)) then
         write(*,'(f12.4,1x)', advance='no') x(i, j)
      else
         write(*,'(es12.5,1x)', advance='no') x(i, j)
      end if
   end do
   write(*,*)
end do
end subroutine print_matrix_rstyle_real

subroutine print_matrix_rstyle_named_real(x, names, int_cols, row_names, digits)
! Print a real matrix with R-like row labels and provided column names.
real(kind=dp), intent(in) :: x(:,:) ! matrix to print
character(len=*), intent(in) :: names(:) ! column names
logical, intent(in), optional :: int_cols(:) ! columns to print as integers
character(len=*), intent(in), optional :: row_names(:) ! row names
integer, intent(in), optional :: digits ! digits after decimal point for display
integer :: i, j
logical :: as_int_col
character(len=32) :: fmt
if (present(digits)) write(fmt, '("(f12.",i0,",1x)")') max(0, digits)
write(*,'(5x)', advance='no')
do i = 1, size(x, 2)
   if (i <= size(names)) then
      write(*,'(a12,1x)', advance='no') trim(names(i))
   else
      write(*,'("[,",i0,"]",8x)', advance='no') i
   end if
end do
write(*,*)
do i = 1, size(x, 1)
   if (present(row_names) .and. i <= size(row_names)) then
      write(*,'(a6,1x)', advance='no') trim(row_names(i))
   else
      write(*,'("[",i0,",]",1x)', advance='no') i
   end if
   do j = 1, size(x, 2)
      as_int_col = .false.
      if (present(int_cols)) then
         if (j <= size(int_cols)) as_int_col = int_cols(j)
      end if
      if (as_int_col .and. ieee_is_finite(x(i, j))) then
         write(*,'(i12,1x)', advance='no') nint(x(i, j))
      else if (.not. ieee_is_finite(x(i, j))) then
         write(*,'(a12,1x)', advance='no') "NA"
      else if (present(digits)) then
         write(*,fmt, advance='no') x(i, j)
      else
         write(*,'(f12.4,1x)', advance='no') x(i, j)
      end if
   end do
   write(*,*)
end do
end subroutine print_matrix_rstyle_named_real

subroutine print_matrix_rstyle_named_int(x, names)
! Print an integer matrix with R-like row labels and provided column names.
integer, intent(in) :: x(:,:) ! matrix to print
character(len=*), intent(in) :: names(:) ! column names
integer :: i, j
write(*,'(5x)', advance='no')
do i = 1, size(x, 2)
   if (i <= size(names)) then
      write(*,'(a12,1x)', advance='no') trim(names(i))
   else
      write(*,'("[,",i0,"]",8x)', advance='no') i
   end if
end do
write(*,*)
do i = 1, size(x, 1)
   write(*,'("[",i0,",]",1x)', advance='no') i
   do j = 1, size(x, 2)
      write(*,'(i12,1x)', advance='no') x(i, j)
   end do
   write(*,*)
end do
end subroutine print_matrix_rstyle_named_int

subroutine print_matrix_int(x)
! Print an integer matrix row-by-row.
integer, intent(in) :: x(:,:) ! matrix to print
integer :: i
do i = 1, size(x, 1)
   write(*,"(*(i0,1x))") x(i, :)
end do
end subroutine print_matrix_int

subroutine print_matrix_complex(x)
! Print a complex matrix row-by-row.
complex(kind=dp), intent(in) :: x(:,:) ! matrix to print
integer :: i
do i = 1, size(x, 1)
   call print_complex_vector(x(i, :))
end do
end subroutine print_matrix_complex

subroutine print_matrix_logical(x)
! Print a logical matrix row-by-row.
logical, intent(in) :: x(:,:) ! matrix to print
integer :: i, j
do i = 1, size(x, 1)
   do j = 1, size(x, 2)
      if (x(i, j)) then
         write(*,'(a)', advance='no') 'T'
      else
         write(*,'(a)', advance='no') 'F'
      end if
      if (j < size(x, 2)) write(*,'(a)', advance='no') ' '
   end do
   write(*,*)
end do
end subroutine print_matrix_logical

subroutine print_matrix_rstyle_int(x)
! Print an integer matrix with R-like column and row labels.
integer, intent(in) :: x(:,:) ! matrix to print
integer :: i
character(len=12) :: col_label
write(*,'(5x)', advance='no')
do i = 1, size(x, 2)
   write(col_label, '("[,",i0,"]")') i
   write(*,'(a12,1x)', advance='no') adjustr(col_label)
end do
write(*,*)
do i = 1, size(x, 1)
   write(*,'("[",i0,",]",1x)', advance='no') i
   write(*,'(*(i12,1x))') x(i, :)
end do
end subroutine print_matrix_rstyle_int

subroutine print_matrix_rstyle_logical(x)
! Print a logical matrix with R-like column and row labels.
logical, intent(in) :: x(:,:) ! matrix to print
call print_matrix_logical(x)
end subroutine print_matrix_rstyle_logical

pure function lm_predict_general(fit, xpred) result(yhat)
! Predict responses for a fitted linear model.
type(lm_fit_t), intent(in) :: fit ! fitted linear model
real(kind=dp), intent(in) :: xpred(:,:) ! predictor rows without intercept column
real(kind=dp), allocatable :: yhat(:)
integer :: p
p = size(xpred, 2)
if (fit%has_intercept) then
   if (size(fit%coef) /= p + 1) error stop "error: predictor count mismatch"
else
   if (size(fit%coef) /= p) error stop "error: predictor count mismatch"
end if
allocate(yhat(size(xpred, 1)))
if (fit%has_intercept) then
   yhat = fit%coef(1) + matmul(xpred, fit%coef(2:p+1))
else
   yhat = matmul(xpred, fit%coef)
end if
end function lm_predict_general

pure function lm_predict_interval(fit, xpred) result(out)
! Prediction intervals for a fitted linear model.
type(lm_fit_t), intent(in) :: fit ! fitted linear model
real(kind=dp), intent(in) :: xpred(:,:) ! predictor rows without intercept column
real(kind=dp), allocatable :: out(:,:)
real(kind=dp), allocatable :: x0(:), yhat(:)
real(kind=dp) :: crit, se_pred
integer :: i, k, p
p = size(xpred, 2)
yhat = lm_predict_general(fit, xpred)
k = p + merge(1, 0, fit%has_intercept)
allocate(out(size(xpred, 1), 3), x0(k))
crit = t_crit_975(fit%df)
do i = 1, size(xpred, 1)
   if (fit%has_intercept) then
      x0(1) = 1.0_dp
      if (p > 0) x0(2:p+1) = xpred(i, :)
   else if (p > 0) then
      x0 = xpred(i, :)
   end if
   se_pred = fit%sigma
   if (allocated(fit%cov_unscaled) .and. size(fit%cov_unscaled, 1) >= k .and. size(fit%cov_unscaled, 2) >= k) then
      se_pred = fit%sigma * sqrt(max(0.0_dp, 1.0_dp + dot_product(x0, matmul(fit%cov_unscaled, x0))))
   end if
   out(i, 1) = yhat(i)
   out(i, 2) = yhat(i) - crit * se_pred
   out(i, 3) = yhat(i) + crit * se_pred
end do
end function lm_predict_interval

subroutine print_lm_prediction_interval(fit, xpred)
! Print lm prediction interval values in an R-like format.
type(lm_fit_t), intent(in) :: fit ! fitted linear model
real(kind=dp), intent(in) :: xpred(:,:) ! predictor rows without intercept column
real(kind=dp), allocatable :: pred(:,:)
integer :: i
pred = lm_predict_interval(fit, xpred)
write(*,'(5x,a9,1x,a9,1x,a9)') "fit", "lwr", "upr"
do i = 1, size(pred, 1)
   write(*,'(i0,1x,f9.6,1x,f9.6,1x,f9.5)') i, pred(i, 1), pred(i, 2), pred(i, 3)
end do
end subroutine print_lm_prediction_interval

pure function lm_confint(fit, level) result(out)
! Coefficient confidence intervals for a fitted linear model.
! level is accepted for API compatibility; this subset uses 95% intervals.
type(lm_fit_t), intent(in) :: fit ! fitted linear model
real(kind=dp), intent(in), optional :: level ! accepted for compatibility; ignored
real(kind=dp), allocatable :: out(:,:)
real(kind=dp) :: crit, se
integer :: j, p
if (present(level)) continue
p = size(fit%coef)
allocate(out(p, 2))
if (p <= 0) return
crit = t_crit_975(fit%df)
do j = 1, p
   se = fit%sigma
   if (allocated(fit%cov_unscaled) .and. size(fit%cov_unscaled, 1) >= j .and. size(fit%cov_unscaled, 2) >= j) then
      se = fit%sigma * sqrt(max(0.0_dp, fit%cov_unscaled(j, j)))
   end if
   out(j, 1) = fit%coef(j) - crit * se
   out(j, 2) = fit%coef(j) + crit * se
end do
end function lm_confint

pure function t_crit_975(df) result(out)
! Two-sided 95% Student-t critical value for common residual degrees of freedom.
integer, intent(in) :: df ! residual degrees of freedom
real(kind=dp) :: out
select case (df)
case (1); out = 12.7062047364321_dp
case (2); out = 4.30265272974946_dp
case (3); out = 3.18244630528426_dp
case (4); out = 2.77644510519779_dp
case (5); out = 2.57058183563631_dp
case (6); out = 2.44691184879168_dp
case (7); out = 2.36462425159278_dp
case (8); out = 2.30600413503337_dp
case (9); out = 2.26215716274099_dp
case (10); out = 2.22813885196494_dp
case (11); out = 2.20098516008295_dp
case (12); out = 2.17881282966342_dp
case (13); out = 2.16036865646101_dp
case (14); out = 2.14478668791693_dp
case (15); out = 2.13144954555978_dp
case (16); out = 2.11990529922125_dp
case (17); out = 2.10981557783318_dp
case (18); out = 2.10092204024096_dp
case (19); out = 2.09302405440831_dp
case (20); out = 2.08596344726586_dp
case (21); out = 2.07961384472766_dp
case (22); out = 2.07387306790401_dp
case (23); out = 2.06865761041904_dp
case (24); out = 2.06389856162802_dp
case (25); out = 2.05953855275329_dp
case (26); out = 2.05552943864287_dp
case (27); out = 2.05183051648028_dp
case (28); out = 2.04840714179524_dp
case (29); out = 2.04522964213270_dp
case (30); out = 2.04227245630124_dp
case default
   if (df <= 0) then
      out = ieee_value(0.0_dp, ieee_quiet_nan)
   else if (df <= 40) then
      out = 2.02107539030627_dp
   else if (df <= 60) then
      out = 2.00029782105826_dp
   else if (df <= 120) then
      out = 1.97993040505278_dp
   else
      out = 1.959963984540054_dp
   end if
end select
end function t_crit_975

pure function lm_cooks_distance(fit) result(out)
! Cook's distance for each observation in a fitted linear model.
type(lm_fit_t), intent(in) :: fit ! fitted linear model
real(kind=dp), allocatable :: out(:)
real(kind=dp), allocatable :: x0(:)
real(kind=dp) :: h, denom
integer :: i, p
if (.not. allocated(fit%resid)) then
   allocate(out(0))
   return
end if
p = size(fit%coef)
allocate(out(size(fit%resid)))
out = 0.0_dp
if (.not. allocated(fit%xpred) .or. .not. allocated(fit%cov_unscaled) .or. p <= 0 .or. fit%sigma <= 0.0_dp) return
allocate(x0(p))
do i = 1, size(fit%resid)
   x0(1) = 1.0_dp
   if (p > 1) x0(2:p) = fit%xpred(i, 1:p-1)
   h = dot_product(x0, matmul(fit%cov_unscaled, x0))
   denom = real(p, kind=dp) * fit%sigma * fit%sigma * max(tiny(1.0_dp), (1.0_dp - h) * (1.0_dp - h))
   out(i) = fit%resid(i) * fit%resid(i) * h / denom
end do
end function lm_cooks_distance

function qr(x, tol, lapack) result(fit)
! Compact QR decomposition object for base R-style qr(x) uses.
! Defaults: tol=1.0e-7; lapack is accepted for compatibility but ignored.
real(kind=dp), intent(in) :: x(:,:) ! matrix to decompose
real(kind=dp), intent(in), optional :: tol ! rank tolerance
logical, intent(in), optional :: lapack ! accepted for compatibility; ignored
type(qr_fit_t) :: fit
integer :: i, j, m, n, k, rk
real(kind=dp) :: nrm, eps, max_diag, t
real(kind=dp), allocatable :: v(:)
logical :: lapack_dummy
m = size(x, 1)
n = size(x, 2)
k = min(m, n)
eps = 1.0e-7_dp
if (present(tol)) eps = tol
lapack_dummy = .false.
if (present(lapack)) lapack_dummy = lapack
allocate(fit%qr(m, n), fit%q(m, k), fit%r(k, n), fit%qraux(k), fit%pivot(n), v(m))
fit%qr = x
fit%q = 0.0_dp
fit%r = 0.0_dp
fit%qraux = 0.0_dp
fit%pivot = [(i, i = 1, n)]
rk = 0
max_diag = 0.0_dp
do j = 1, k
   nrm = sqrt(max(0.0_dp, dot_product(fit%qr(j:m, j), fit%qr(j:m, j))))
   if (nrm /= 0.0_dp) then
      if (fit%qr(j, j) /= 0.0_dp) nrm = sign(nrm, fit%qr(j, j))
      fit%qr(j:m, j) = fit%qr(j:m, j) / nrm
      fit%qr(j, j) = 1.0_dp + fit%qr(j, j)
      fit%qraux(j) = fit%qr(j, j)
      do i = j + 1, n
         t = -dot_product(fit%qr(j:m, j), fit%qr(j:m, i)) / fit%qr(j, j)
         fit%qr(j:m, i) = fit%qr(j:m, i) + t * fit%qr(j:m, j)
      end do
      fit%qr(j, j) = -nrm
   end if
   max_diag = max(max_diag, abs(fit%qr(j, j)))
   if (abs(fit%qr(j, j)) > eps * max(1.0_dp, max_diag)) rk = rk + 1
end do
fit%rank = rk
fit%r = 0.0_dp
do i = 1, k
   fit%r(i, i:n) = fit%qr(i, i:n)
end do
end function qr

pure function qr_Q(fit, complete) result(q)
! Reconstruct the Q matrix from compact qr() Householder storage.
! If complete is absent or false, return the thin Q matrix.
type(qr_fit_t), intent(in) :: fit ! QR decomposition result
logical, intent(in), optional :: complete ! request full Q matrix flag
real(kind=dp), allocatable :: q(:,:)
integer :: i, j, col, m, n, nq
real(kind=dp) :: t
real(kind=dp), allocatable :: v(:)
logical :: full
m = size(fit%qr, 1)
n = size(fit%qr, 2)
full = .false.
if (present(complete)) full = complete
nq = merge(m, min(m, n), full)
allocate(q(m, nq))
q = 0.0_dp
do i = 1, min(m, nq)
   q(i, i) = 1.0_dp
end do
do j = min(m, n), 1, -1
   if (j > size(fit%qraux) .or. fit%qraux(j) == 0.0_dp) cycle
   allocate(v(m - j + 1))
   v = fit%qr(j:m, j)
   v(1) = fit%qraux(j)
   do col = 1, nq
      t = -dot_product(v, q(j:m, col)) / v(1)
      q(j:m, col) = q(j:m, col) + t * v
   end do
   deallocate(v)
end do
end function qr_Q

pure function qr_R(fit, complete) result(r)
! Extract the upper-triangular R matrix from compact qr() storage.
! If complete is absent or false, return the thin R matrix.
type(qr_fit_t), intent(in) :: fit ! QR decomposition result
logical, intent(in), optional :: complete ! request full R matrix flag
real(kind=dp), allocatable :: r(:,:)
integer :: i, m, n, nr
logical :: full
m = size(fit%qr, 1)
n = size(fit%qr, 2)
full = .false.
if (present(complete)) full = complete
nr = merge(m, min(m, n), full)
allocate(r(nr, n))
r = 0.0_dp
do i = 1, min(nr, n)
   r(i, i:n) = fit%qr(i, i:n)
end do
end function qr_R

pure function qr_rank(fit) result(out)
! Return QR decomposition rank values.
type(qr_fit_t), intent(in) :: fit ! QR decomposition result
integer :: out
out = fit%rank
end function qr_rank

pure function qr_pivot(fit) result(out)
! Return QR decomposition pivot values.
type(qr_fit_t), intent(in) :: fit ! QR decomposition result
integer, allocatable :: out(:)
allocate(out(size(fit%pivot)))
out = fit%pivot
end function qr_pivot

pure function qr_coef_vec(fit, y) result(coef)
! Return QR decomposition coef vec values.
type(qr_fit_t), intent(in) :: fit ! QR decomposition result
real(kind=dp), intent(in) :: y(:) ! response vector
real(kind=dp), allocatable :: coef(:)
real(kind=dp), allocatable :: q(:,:), r(:,:), qty(:)
integer :: k
k = max(0, fit%rank)
if (k <= 0) then
   allocate(coef(0))
   return
end if
q = qr_Q(fit)
r = qr_R(fit)
qty = matmul(transpose(q(:, 1:k)), y)
coef = backsolve_vec(r(1:k, 1:k), qty(1:k))
end function qr_coef_vec

pure function qr_coef_mat(fit, y) result(coef)
! Return QR decomposition coef mat values.
type(qr_fit_t), intent(in) :: fit ! QR decomposition result
real(kind=dp), intent(in) :: y(:,:) ! response matrix
real(kind=dp), allocatable :: coef(:,:)
real(kind=dp), allocatable :: q(:,:), r(:,:), qty(:,:)
integer :: k
k = max(0, fit%rank)
if (k <= 0) then
   allocate(coef(0, size(y, 2)))
   return
end if
q = qr_Q(fit)
r = qr_R(fit)
qty = matmul(transpose(q(:, 1:k)), y)
coef = backsolve_mat(r(1:k, 1:k), qty(1:k, :))
end function qr_coef_mat

pure function qr_fitted_vec(fit, y) result(out)
! Return QR decomposition fitted vec values.
type(qr_fit_t), intent(in) :: fit ! QR decomposition result
real(kind=dp), intent(in) :: y(:) ! response vector
real(kind=dp), allocatable :: out(:)
real(kind=dp), allocatable :: q(:,:), qty(:)
integer :: k
k = max(0, fit%rank)
allocate(out(size(y)))
if (k <= 0) then
   out = 0.0_dp
   return
end if
q = qr_Q(fit)
qty = matmul(transpose(q(:, 1:k)), y)
out = matmul(q(:, 1:k), qty(1:k))
end function qr_fitted_vec

pure function qr_fitted_mat(fit, y) result(out)
! Return QR decomposition fitted mat values.
type(qr_fit_t), intent(in) :: fit ! QR decomposition result
real(kind=dp), intent(in) :: y(:,:) ! response matrix
real(kind=dp), allocatable :: out(:,:)
real(kind=dp), allocatable :: q(:,:), qty(:,:)
integer :: k
k = max(0, fit%rank)
allocate(out(size(y, 1), size(y, 2)))
if (k <= 0) then
   out = 0.0_dp
   return
end if
q = qr_Q(fit)
qty = matmul(transpose(q(:, 1:k)), y)
out = matmul(q(:, 1:k), qty(1:k, :))
end function qr_fitted_mat

pure function qr_resid_vec(fit, y) result(out)
! Return QR decomposition resid vec values.
type(qr_fit_t), intent(in) :: fit ! QR decomposition result
real(kind=dp), intent(in) :: y(:) ! response vector
real(kind=dp), allocatable :: out(:)
out = y - qr_fitted(fit, y)
end function qr_resid_vec

pure function qr_resid_mat(fit, y) result(out)
! Return QR decomposition resid mat values.
type(qr_fit_t), intent(in) :: fit ! QR decomposition result
real(kind=dp), intent(in) :: y(:,:) ! response matrix
real(kind=dp), allocatable :: out(:,:)
out = y - qr_fitted(fit, y)
end function qr_resid_mat

pure function qr_qty_vec(fit, y) result(out)
! Return QR decomposition qty vec values.
type(qr_fit_t), intent(in) :: fit ! QR decomposition result
real(kind=dp), intent(in) :: y(:) ! vector multiplied by transpose(Q)
real(kind=dp), allocatable :: out(:)
real(kind=dp), allocatable :: q(:,:)
q = qr_Q(fit, complete=.true.)
out = matmul(transpose(q), y)
end function qr_qty_vec

pure function qr_qty_mat(fit, y) result(out)
! Return QR decomposition qty mat values.
type(qr_fit_t), intent(in) :: fit ! QR decomposition result
real(kind=dp), intent(in) :: y(:,:) ! matrix multiplied by transpose(Q)
real(kind=dp), allocatable :: out(:,:)
real(kind=dp), allocatable :: q(:,:)
q = qr_Q(fit, complete=.true.)
out = matmul(transpose(q), y)
end function qr_qty_mat

pure function qr_qy_vec(fit, y) result(out)
! Return QR decomposition qy vec values.
type(qr_fit_t), intent(in) :: fit ! QR decomposition result
real(kind=dp), intent(in) :: y(:) ! vector multiplied by Q
real(kind=dp), allocatable :: out(:)
real(kind=dp), allocatable :: q(:,:)
q = qr_Q(fit, complete=.true.)
out = matmul(q, y)
end function qr_qy_vec

pure function qr_qy_mat(fit, y) result(out)
! Return QR decomposition qy mat values.
type(qr_fit_t), intent(in) :: fit ! QR decomposition result
real(kind=dp), intent(in) :: y(:,:) ! matrix multiplied by Q
real(kind=dp), allocatable :: out(:,:)
real(kind=dp), allocatable :: q(:,:)
q = qr_Q(fit, complete=.true.)
out = matmul(q, y)
end function qr_qy_mat

subroutine print_qr(fit)
! Print qr values in an R-like format.
type(qr_fit_t), intent(in) :: fit ! QR decomposition result
write(*,'(a)') "$qr"
call print_matrix_rstyle(fit%qr)
write(*,'(a)') ""
write(*,'(a)') "$rank"
call print_real_scalar(real(fit%rank, kind=dp))
write(*,'(a)') ""
write(*,'(a)') "$qraux"
call print_real_vector(fit%qraux)
write(*,'(a)') ""
write(*,'(a)') "$pivot"
call print_real_vector(real(fit%pivot, kind=dp))
end subroutine print_qr

subroutine print_lm_cooks_top(fit, n)
! Print the largest Cook's distances with original observation numbers.
type(lm_fit_t), intent(in) :: fit ! fitted linear model
integer, intent(in) :: n ! number of largest distances to print
real(kind=dp), allocatable :: cook(:)
integer, allocatable :: ord(:)
integer :: i, j, m, chunk
cook = lm_cooks_distance(fit)
ord = order_real(cook)
m = max(0, min(n, size(cook)))
if (m <= 0) then
   write(*,*)
   return
end if
i = 1
do while (i <= m)
   chunk = min(7, m - i + 1)
   do j = 0, chunk - 1
      write(*,'(i10,1x)', advance='no') ord(size(ord) - (i + j) + 1)
   end do
   write(*,*)
   do j = 0, chunk - 1
      write(*,'(f10.8,1x)', advance='no') cook(ord(size(ord) - (i + j) + 1))
   end do
   write(*,*)
   i = i + chunk
end do
end subroutine print_lm_cooks_top

pure subroutine solve_linear(a, b, x, ok)
! Solve Ax=b by Gaussian elimination with partial pivoting.
real(kind=dp), intent(inout) :: a(:,:) ! coefficient matrix, overwritten during factorization
real(kind=dp), intent(inout) :: b(:) ! right-hand side vector, overwritten during factorization
real(kind=dp), intent(out) :: x(:) ! solution vector
logical, intent(out) :: ok ! success flag
integer :: i, j, k, p, n
real(kind=dp) :: piv, fac, t
ok = .true.
n = size(b)
if (size(a,1) /= n .or. size(a,2) /= n .or. size(x) /= n) then
   ok = .false.
   x = 0.0_dp
   return
end if
do k = 1, n
   p = k
   piv = abs(a(k,k))
   do i = k + 1, n
      if (abs(a(i,k)) > piv) then
         p = i
         piv = abs(a(i,k))
      end if
   end do
   if (piv <= tiny(1.0_dp)) then
      ok = .false.
      x = 0.0_dp
      return
   end if
   if (p /= k) then
      do j = k, n
         t = a(k,j)
         a(k,j) = a(p,j)
         a(p,j) = t
      end do
      t = b(k)
      b(k) = b(p)
      b(p) = t
   end if
   do i = k + 1, n
      fac = a(i,k) / a(k,k)
      a(i,k:n) = a(i,k:n) - fac * a(k,k:n)
      b(i) = b(i) - fac * b(k)
   end do
end do
x(n) = b(n) / a(n,n)
do i = n - 1, 1, -1
   if (i < n) then
      x(i) = (b(i) - sum(a(i,i+1:n) * x(i+1:n))) / a(i,i)
   else
      x(i) = b(i) / a(i,i)
   end if
end do
end subroutine solve_linear

function lm_fit_general(y, xpred, intercept) result(fit)
! Fit linear regression with optional intercept and arbitrary predictors.
real(kind=dp), intent(in) :: y(:) ! response values
real(kind=dp), intent(in) :: xpred(:,:) ! predictor rows without intercept column
logical, intent(in), optional :: intercept ! include an intercept column; default true
type(lm_fit_t) :: fit
integer :: i, j, n, p, k, dof
real(kind=dp), allocatable :: a(:,:), a_xtx(:,:), a_work(:,:), b(:), b_work(:), beta(:), sol(:)
real(kind=dp) :: ybar, sse, sst
logical :: ok
if (size(y) /= size(xpred,1)) then
   error stop "error: need size(y) == size(xpred,1)"
end if
n = size(y)
p = size(xpred,2)
fit%has_intercept = .true.
if (present(intercept)) fit%has_intercept = intercept
k = p + merge(1, 0, fit%has_intercept)
if (n < k) error stop "error: need n >= number of parameters"
allocate(a(k,k), b(k), beta(k))
a = 0.0_dp
b = 0.0_dp
if (fit%has_intercept) then
   a(1,1) = n
   b(1) = sum(y)
   do j = 1, p
      a(1,j+1) = sum(xpred(:,j))
      a(j+1,1) = a(1,j+1)
      b(j+1) = sum(xpred(:,j) * y)
   end do
   do i = 1, p
      do j = i, p
         a(i+1,j+1) = sum(xpred(:,i) * xpred(:,j))
         a(j+1,i+1) = a(i+1,j+1)
      end do
   end do
else
   do j = 1, p
      b(j) = sum(xpred(:,j) * y)
   end do
   do i = 1, p
      do j = i, p
         a(i,j) = sum(xpred(:,i) * xpred(:,j))
         a(j,i) = a(i,j)
      end do
   end do
end if
a_xtx = a
call solve_linear(a, b, beta, ok)
if (.not. ok) error stop "error: singular normal equations"
allocate(fit%cov_unscaled(k, k), a_work(k, k), b_work(k), sol(k))
do j = 1, k
   a_work = a_xtx
   b_work = 0.0_dp
   b_work(j) = 1.0_dp
   call solve_linear(a_work, b_work, sol, ok)
   if (.not. ok) error stop "error: singular normal equations"
   fit%cov_unscaled(:, j) = sol
end do
fit%coef = beta
fit%y = y
fit%xpred = xpred
if (fit%has_intercept) then
   fit%fitted = beta(1) + matmul(xpred, beta(2:k))
else
   fit%fitted = matmul(xpred, beta)
end if
fit%resid = y - fit%fitted
sse = sum(fit%resid**2)
ybar = sum(y) / n
sst = sum((y - ybar)**2)
if (sst > 0.0_dp) then
   fit%r_squared = 1.0_dp - sse / sst
else
   fit%r_squared = 0.0_dp
end if
dof = max(1, n - k)
fit%df = n - k
fit%sigma = sqrt(sse / dof)
fit%adj_r_squared = 1.0_dp - (1.0_dp - fit%r_squared) * (n - 1) / dof
end function lm_fit_general

function step_lm(lower, upper, k) result(best_fit)
! Stepwise add/drop search over upper%xpred columns using an AIC-like score.
! Default k=2 gives AIC-style scoring; k=log(n) gives BIC-style scoring.
type(lm_fit_t), intent(in) :: lower ! lower-scope fitted model
type(lm_fit_t), intent(in) :: upper ! upper-scope fitted model
real(kind=dp), intent(in), optional :: k ! penalty per fitted parameter
type(lm_fit_t) :: best_fit, cand_fit
logical, allocatable :: selected(:), cand_selected(:)
real(kind=dp), allocatable :: xsel(:,:), xcand(:,:)
real(kind=dp) :: kval, best_score, cand_score
integer :: p, j, iter
logical :: improved
kval = 2.0_dp
if (present(k)) kval = k
p = size(upper%xpred, 2)
allocate(selected(p), cand_selected(p))
selected = .false.
if (allocated(lower%xpred)) then
   do j = 1, min(size(lower%xpred, 2), p)
      selected(j) = .true.
   end do
end if
call build_lm_design(upper%xpred, selected, xsel)
best_fit = lm_fit_general(upper%y, xsel)
best_score = lm_aic_score(best_fit, kval)
do iter = 1, max(1, 2 * p + 2)
   improved = .false.
   cand_selected = selected
   do j = 1, p
      if (.not. selected(j)) then
         cand_selected = selected
         cand_selected(j) = .true.
         call build_lm_design(upper%xpred, cand_selected, xcand)
         cand_fit = lm_fit_general(upper%y, xcand)
         cand_score = lm_aic_score(cand_fit, kval)
         if (cand_score < best_score - 1.0e-8_dp) then
            best_score = cand_score
            best_fit = cand_fit
            selected = cand_selected
            improved = .true.
         end if
      end if
   end do
   do j = 1, p
      if (selected(j) .and. count(selected) > 0) then
         cand_selected = selected
         cand_selected(j) = .false.
         call build_lm_design(upper%xpred, cand_selected, xcand)
         cand_fit = lm_fit_general(upper%y, xcand)
         cand_score = lm_aic_score(cand_fit, kval)
         if (cand_score < best_score - 1.0e-8_dp) then
            best_score = cand_score
            best_fit = cand_fit
            selected = cand_selected
            improved = .true.
         end if
      end if
   end do
   if (.not. improved) exit
end do
end function step_lm

pure function lm_aic_score(fit, k) result(score)
! Support linear-model helper lm_aic_score.
type(lm_fit_t), intent(in) :: fit ! fitted linear model
real(kind=dp), intent(in) :: k ! penalty per fitted parameter
real(kind=dp) :: score, rss, n
n = real(size(fit%y), kind=dp)
rss = max(tiny(1.0_dp), sum(fit%resid**2))
score = n * log(rss / max(1.0_dp, n)) + k * real(size(fit%coef), kind=dp)
end function lm_aic_score

pure subroutine build_lm_design(x, selected, out)
! Support linear-model helper build_lm_design.
real(kind=dp), intent(in) :: x(:,:) ! full predictor matrix
logical, intent(in) :: selected(:) ! selected predictor columns
real(kind=dp), allocatable, intent(out) :: out(:,:) ! reduced predictor matrix
integer :: j, jj, n, p
n = size(x, 1)
p = count(selected)
allocate(out(n, p))
jj = 0
do j = 1, min(size(x, 2), size(selected))
   if (selected(j)) then
      jj = jj + 1
      out(:, jj) = x(:, j)
   end if
end do
end subroutine build_lm_design

function lm_coef(y, xpred, intercept) result(coef)
! Fit linear model and return only coefficient vector.
real(kind=dp), intent(in) :: y(:) ! response values
real(kind=dp), intent(in) :: xpred(:,:) ! predictor rows without intercept column
logical, intent(in), optional :: intercept ! include an intercept column; default true
real(kind=dp), allocatable :: coef(:)
type(lm_fit_t) :: fit
logical :: use_intercept
use_intercept = .true.
if (present(intercept)) use_intercept = intercept
fit = lm_fit_general(y, xpred, intercept=use_intercept)
coef = fit%coef
end function lm_coef

function lm_r_squared_general(y, xpred, intercept) result(out)
! Fit linear model and return only R-squared.
real(kind=dp), intent(in) :: y(:) ! response values
real(kind=dp), intent(in) :: xpred(:,:) ! predictor rows without intercept column
logical, intent(in), optional :: intercept ! include an intercept column; default true
real(kind=dp) :: out
type(lm_fit_t) :: fit
logical :: use_intercept
use_intercept = .true.
if (present(intercept)) use_intercept = intercept
fit = lm_fit_general(y, xpred, intercept=use_intercept)
out = fit%r_squared
end function lm_r_squared_general

function glm_binomial_fit_int(y, xpred) result(fit)
! Fit a binomial-logit GLM for integer 0/1 responses.
integer, intent(in) :: y(:) ! binary response values
real(kind=dp), intent(in) :: xpred(:,:) ! predictor rows without intercept column
type(glm_fit_t) :: fit
fit = glm_binomial_fit_real(real(y, kind=dp), xpred)
end function glm_binomial_fit_int

function glm_binomial_fit_real(y, xpred) result(fit)
! Fit a binomial-logit GLM by Newton/IRLS normal equations.
real(kind=dp), intent(in) :: y(:) ! binary response values
real(kind=dp), intent(in) :: xpred(:,:) ! predictor rows without intercept column
type(glm_fit_t) :: fit
integer :: i, j, l, n, p, k, iter
real(kind=dp), allocatable :: beta(:), grad(:), hess(:,:), hess0(:,:), step(:), eta(:), mu(:), w(:), rhs(:)
real(kind=dp) :: e, max_step
logical :: ok
if (size(y) /= size(xpred, 1)) error stop "error: glm response/predictor size mismatch"
n = size(y)
p = size(xpred, 2)
k = p + 1
allocate(beta(k), grad(k), hess(k,k), hess0(k,k), step(k), eta(n), mu(n), w(n), rhs(k))
beta = 0.0_dp
do iter = 1, 50
   eta = beta(1)
   if (p > 0) eta = eta + matmul(xpred, beta(2:k))
   eta = max(-35.0_dp, min(35.0_dp, eta))
   mu = 1.0_dp / (1.0_dp + exp(-eta))
   w = max(1.0e-12_dp, mu * (1.0_dp - mu))
   grad = 0.0_dp
   hess = 0.0_dp
   do i = 1, n
      e = y(i) - mu(i)
      grad(1) = grad(1) + e
      hess(1,1) = hess(1,1) + w(i)
      do j = 1, p
         grad(j+1) = grad(j+1) + xpred(i,j) * e
         hess(1,j+1) = hess(1,j+1) + w(i) * xpred(i,j)
         hess(j+1,1) = hess(1,j+1)
         do l = 1, p
            hess(j+1,l+1) = hess(j+1,l+1) + w(i) * xpred(i,j) * xpred(i,l)
         end do
      end do
   end do
   hess0 = hess
   rhs = grad
   call solve_linear(hess, rhs, step, ok)
   if (.not. ok) exit
   beta = beta + step
   max_step = maxval(abs(step))
   if (max_step < 1.0e-10_dp) then
      fit%convergence = 0
      exit
   end if
end do
fit%iter = iter
fit%family = 1
fit%coef = beta
fit%y = y
fit%xpred = xpred
allocate(fit%offset(n))
fit%offset = 0.0_dp
fit%fitted = glm_predict_response(fit, xpred)
fit%resid = y - fit%fitted
fit%df = n - k
allocate(fit%se(k), fit%z_value(k), fit%p_value(k))
fit%se = ieee_value(0.0_dp, ieee_quiet_nan)
fit%z_value = ieee_value(0.0_dp, ieee_quiet_nan)
fit%p_value = ieee_value(0.0_dp, ieee_quiet_nan)
do j = 1, k
   hess = hess0
   rhs = 0.0_dp
   rhs(j) = 1.0_dp
   call solve_linear(hess, rhs, step, ok)
   if (ok .and. step(j) >= 0.0_dp) then
      fit%se(j) = sqrt(step(j))
      if (fit%se(j) > 0.0_dp) then
         fit%z_value(j) = fit%coef(j) / fit%se(j)
         fit%p_value(j) = 2.0_dp * max(0.0_dp, min(1.0_dp, 1.0_dp - normal_cdf(abs(fit%z_value(j)))))
      end if
   end if
end do
end function glm_binomial_fit_real

function glm_poisson_fit_int(y, xpred, offset) result(fit)
! Fit a Poisson-log GLM for integer count responses.
integer, intent(in) :: y(:) ! count response values
real(kind=dp), intent(in) :: xpred(:,:) ! predictor rows without intercept column
real(kind=dp), intent(in), optional :: offset(:) ! additive linear-predictor offset
type(glm_fit_t) :: fit
fit = glm_poisson_fit_real(real(y, kind=dp), xpred, offset)
end function glm_poisson_fit_int

function glm_poisson_fit_real(y, xpred, offset) result(fit)
! Fit a Poisson-log GLM by Newton/IRLS normal equations.
real(kind=dp), intent(in) :: y(:) ! count response values
real(kind=dp), intent(in) :: xpred(:,:) ! predictor rows without intercept column
real(kind=dp), intent(in), optional :: offset(:) ! additive linear-predictor offset
type(glm_fit_t) :: fit
integer :: i, j, l, n, p, k, iter
real(kind=dp), allocatable :: beta(:), grad(:), hess(:,:), hess0(:,:), step(:), eta(:), mu(:), off(:), rhs(:)
real(kind=dp) :: e, max_step, base_rate
logical :: ok
if (size(y) /= size(xpred, 1)) error stop "error: glm response/predictor size mismatch"
n = size(y)
p = size(xpred, 2)
k = p + 1
allocate(beta(k), grad(k), hess(k,k), hess0(k,k), step(k), eta(n), mu(n), off(n), rhs(k))
off = 0.0_dp
if (present(offset)) then
   if (size(offset) /= n) error stop "error: glm offset size mismatch"
   off = offset
end if
base_rate = max(tiny(1.0_dp), (sum(y) + 0.1_dp) / real(max(1, n), kind=dp))
beta = 0.0_dp
beta(1) = log(base_rate)
do iter = 1, 80
   eta = off + beta(1)
   if (p > 0) eta = eta + matmul(xpred, beta(2:k))
   eta = max(-35.0_dp, min(35.0_dp, eta))
   mu = max(tiny(1.0_dp), exp(eta))
   grad = 0.0_dp
   hess = 0.0_dp
   do i = 1, n
      e = y(i) - mu(i)
      grad(1) = grad(1) + e
      hess(1,1) = hess(1,1) + mu(i)
      do j = 1, p
         grad(j+1) = grad(j+1) + xpred(i,j) * e
         hess(1,j+1) = hess(1,j+1) + mu(i) * xpred(i,j)
         hess(j+1,1) = hess(1,j+1)
         do l = 1, p
            hess(j+1,l+1) = hess(j+1,l+1) + mu(i) * xpred(i,j) * xpred(i,l)
         end do
      end do
   end do
   hess0 = hess
   rhs = grad
   call solve_linear(hess, rhs, step, ok)
   if (.not. ok) exit
   beta = beta + step
   max_step = maxval(abs(step))
   if (max_step < 1.0e-10_dp) then
      fit%convergence = 0
      exit
   end if
end do
fit%family = 2
fit%iter = iter
fit%coef = beta
fit%y = y
fit%xpred = xpred
fit%offset = off
fit%fitted = glm_predict_response(fit, xpred)
fit%resid = y - fit%fitted
fit%df = n - k
allocate(fit%se(k), fit%z_value(k), fit%p_value(k))
fit%se = ieee_value(0.0_dp, ieee_quiet_nan)
fit%z_value = ieee_value(0.0_dp, ieee_quiet_nan)
fit%p_value = ieee_value(0.0_dp, ieee_quiet_nan)
do j = 1, k
   hess = hess0
   rhs = 0.0_dp
   rhs(j) = 1.0_dp
   call solve_linear(hess, rhs, step, ok)
   if (ok .and. step(j) >= 0.0_dp) then
      fit%se(j) = sqrt(step(j))
      if (fit%se(j) > 0.0_dp) then
         fit%z_value(j) = fit%coef(j) / fit%se(j)
         fit%p_value(j) = 2.0_dp * max(0.0_dp, min(1.0_dp, 1.0_dp - normal_cdf(abs(fit%z_value(j)))))
      end if
   end if
end do
end function glm_poisson_fit_real

pure function glm_predict_response(fit, xpred) result(out)
! Predicted responses for a GLM.
type(glm_fit_t), intent(in) :: fit ! fitted GLM
real(kind=dp), intent(in) :: xpred(:,:) ! predictor rows without intercept column
real(kind=dp), allocatable :: out(:)
real(kind=dp), allocatable :: eta(:)
integer :: p
p = size(xpred, 2)
allocate(out(size(xpred, 1)), eta(size(xpred, 1)))
eta = fit%coef(1)
if (p > 0) eta = eta + matmul(xpred, fit%coef(2:p+1))
eta = max(-35.0_dp, min(35.0_dp, eta))
if (fit%family == 2) then
   if (allocated(fit%offset) .and. size(fit%offset) == size(eta)) eta = eta + fit%offset
   out = exp(eta)
else
   out = 1.0_dp / (1.0_dp + exp(-eta))
end if
end function glm_predict_response

pure function glm_pearson_resid(fit) result(out)
! Pearson residuals for supported GLM families.
type(glm_fit_t), intent(in) :: fit ! fitted GLM
real(kind=dp), allocatable :: out(:)
real(kind=dp), allocatable :: v(:)
allocate(out(size(fit%fitted)), v(size(fit%fitted)))
if (fit%family == 2) then
   v = max(tiny(1.0_dp), fit%fitted)
else
   v = max(tiny(1.0_dp), fit%fitted * (1.0_dp - fit%fitted))
end if
out = (fit%y - fit%fitted) / sqrt(v)
end function glm_pearson_resid

subroutine print_glm_summary(fit, term_names)
! Print a compact coefficient table for a binomial-logit GLM.
type(glm_fit_t), intent(in) :: fit ! fitted GLM
character(len=*), intent(in), optional :: term_names(:) ! predictor names excluding intercept
integer :: j, p
character(len=32) :: lbl
write(*,'(a)') "glm summary:"
write(*,'(a)') "Coefficients:"
write(*,'(a12,1x,a10,1x,a10,1x,a8,1x,a10)') "", "Estimate", "Std. Error", "z value", "Pr(>|z|)"
p = size(fit%coef)
do j = 1, p
   if (j == 1) then
      lbl = "(Intercept)"
   else if (present(term_names) .and. size(term_names) >= j - 1) then
      lbl = trim(term_names(j - 1))
   else
      write(lbl,'(a,i0)') "x", j - 1
   end if
   write(*,'(a12,1x,f10.4,1x,f10.4,1x,f8.3,1x,es10.3)') &
      & trim(lbl), fit%coef(j), fit%se(j), fit%z_value(j), fit%p_value(j)
end do
write(*,'(a,i0)') "convergence: ", fit%convergence
end subroutine print_glm_summary

pure function polyroot(coef) result(root_mod)
! Return moduli of roots for coefficients in R polyroot order.
real(kind=dp), intent(in) :: coef(:) ! polynomial coefficients in increasing order
real(kind=dp), allocatable :: root_mod(:)
integer :: n, i, j, iter
complex(kind=dp), allocatable :: z(:), znew(:)
complex(kind=dp) :: pz, denom
real(kind=dp) :: theta, max_delta
real(kind=dp), parameter :: twopi = 6.2831853071795864769252867665590058_dp
n = size(coef) - 1
allocate(root_mod(max(0, n)))
if (n <= 0) return
if (abs(coef(n + 1)) <= tiny(1.0_dp)) then
   root_mod = huge(1.0_dp)
   return
end if
allocate(z(n), znew(n))
do i = 1, n
   theta = twopi * real(i - 1, kind=dp) / real(n, kind=dp)
   z(i) = cmplx(cos(theta), sin(theta), kind=dp)
end do
do iter = 1, 200
   max_delta = 0.0_dp
   znew = z
   do i = 1, n
      pz = cmplx(coef(n + 1), 0.0_dp, kind=dp)
      do j = n, 1, -1
         pz = pz * z(i) + cmplx(coef(j), 0.0_dp, kind=dp)
      end do
      denom = cmplx(1.0_dp, 0.0_dp, kind=dp)
      do j = 1, n
         if (j /= i) denom = denom * (z(i) - z(j))
      end do
      if (abs(denom) > tiny(1.0_dp)) znew(i) = z(i) - (pz / cmplx(coef(n + 1), 0.0_dp, kind=dp)) / denom
      max_delta = max(max_delta, abs(znew(i) - z(i)))
   end do
   z = znew
   if (max_delta <= 100.0_dp * epsilon(1.0_dp)) exit
end do
root_mod = abs(z)
end function polyroot

subroutine print_lm_summary(fit, term_names)
! Print a compact summary of fitted linear model diagnostics.
type(lm_fit_t), intent(in) :: fit ! fitted linear model
character(len=*), intent(in), optional :: term_names(:) ! predictor names excluding intercept
integer :: j, p
real(kind=dp) :: se, tval, pval
real(kind=dp), allocatable :: rq(:)
character(len=32) :: lbl
write(*,'(a)') "lm summary:"
if (allocated(fit%resid) .and. size(fit%resid) > 0) then
   rq = quantile(fit%resid, [0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp], .false., 7)
   write(*,'(a)') "Residuals:"
   write(*,'(a9,1x,a9,1x,a9,1x,a9,1x,a9)') "Min", "1Q", "Median", "3Q", "Max"
   write(*,'(f9.4,1x,f9.4,1x,f9.4,1x,f9.4,1x,f9.4)') rq
end if
write(*,'(a)') "Coefficients:"
write(*,'(a12,1x,a10,1x,a10,1x,a8,1x,a10)') "", "Estimate", "Std. Error", "t value", "Pr(>|t|)"
p = size(fit%coef)
do j = 1, p
   if (fit%has_intercept .and. j == 1) then
      lbl = "(Intercept)"
   else if (present(term_names) .and. size(term_names) >= j - merge(1, 0, fit%has_intercept)) then
      lbl = trim(term_names(j - merge(1, 0, fit%has_intercept)))
   else if (fit%has_intercept) then
      write(lbl,'(a,i0)') "x", j - 1
   else
      write(lbl,'(a,i0)') "x", j
   end if
   se = fit%sigma
   if (allocated(fit%cov_unscaled) .and. size(fit%cov_unscaled, 1) >= j .and. size(fit%cov_unscaled, 2) >= j) then
      se = fit%sigma * sqrt(max(0.0_dp, fit%cov_unscaled(j, j)))
   end if
   if (se > 0.0_dp) then
      tval = fit%coef(j) / se
      pval = student_t_two_sided_pvalue(tval, fit%df)
   else
      tval = ieee_value(0.0_dp, ieee_quiet_nan)
      pval = ieee_value(0.0_dp, ieee_quiet_nan)
   end if
   write(*,'(a12,1x,f10.4,1x,f10.4,1x,f8.3,1x,es10.3)') trim(lbl), fit%coef(j), se, tval, pval
end do
write(*,*)
write(*,'(a,g0)') "sigma: ", fit%sigma
write(*,'(a,g0)') "r.squared: ", fit%r_squared
write(*,'(a,g0)') "adj.r.squared: ", fit%adj_r_squared
end subroutine print_lm_summary

pure function student_t_two_sided_pvalue(t, df) result(p)
! Evaluate numerical special-function helper student_t_two_sided_pvalue.
real(kind=dp), intent(in) :: t ! t statistic
integer, intent(in) :: df ! degrees of freedom
real(kind=dp) :: p, x, a, b
if (df <= 0) then
   p = ieee_value(0.0_dp, ieee_quiet_nan)
   return
end if
x = real(df, kind=dp) / (real(df, kind=dp) + t * t)
a = 0.5_dp * real(df, kind=dp)
b = 0.5_dp
p = regularized_beta(x, a, b)
p = max(0.0_dp, min(1.0_dp, p))
end function student_t_two_sided_pvalue

pure function regularized_beta(x, a, b) result(bt)
! Evaluate numerical special-function helper regularized_beta.
real(kind=dp), intent(in) :: x ! evaluation point in [0,1]
real(kind=dp), intent(in) :: a ! first shape parameter
real(kind=dp), intent(in) :: b ! second shape parameter
real(kind=dp) :: bt, front
if (x <= 0.0_dp) then
   bt = 0.0_dp
   return
end if
if (x >= 1.0_dp) then
   bt = 1.0_dp
   return
end if
front = exp(log_gamma(a + b) - log_gamma(a) - log_gamma(b) + a * log(x) + b * log(1.0_dp - x))
if (x < (a + 1.0_dp) / (a + b + 2.0_dp)) then
   bt = front * beta_cont_frac(a, b, x) / a
else
   bt = 1.0_dp - front * beta_cont_frac(b, a, 1.0_dp - x) / b
end if
end function regularized_beta

pure function beta_cont_frac(a, b, x) result(h)
! Evaluate numerical special-function helper beta_cont_frac.
real(kind=dp), intent(in) :: a ! first shape parameter
real(kind=dp), intent(in) :: b ! second shape parameter
real(kind=dp), intent(in) :: x ! evaluation point
real(kind=dp) :: h, aa, c, d, del, qab, qam, qap
integer :: m, m2
real(kind=dp), parameter :: fpmin = tiny(1.0_dp) / epsilon(1.0_dp)
qab = a + b
qap = a + 1.0_dp
qam = a - 1.0_dp
c = 1.0_dp
d = 1.0_dp - qab * x / qap
if (abs(d) < fpmin) d = fpmin
d = 1.0_dp / d
h = d
do m = 1, 200
   m2 = 2 * m
   aa = real(m, kind=dp) * (b - real(m, kind=dp)) * x / ((qam + real(m2, kind=dp)) * (a + real(m2, kind=dp)))
   d = 1.0_dp + aa * d
   if (abs(d) < fpmin) d = fpmin
   c = 1.0_dp + aa / c
   if (abs(c) < fpmin) c = fpmin
   d = 1.0_dp / d
   h = h * d * c
   aa = -(a + real(m, kind=dp)) * (qab + real(m, kind=dp)) * x / ((a + real(m2, kind=dp)) * (qap + real(m2, kind=dp)))
   d = 1.0_dp + aa * d
   if (abs(d) < fpmin) d = fpmin
   c = 1.0_dp + aa / c
   if (abs(c) < fpmin) c = fpmin
   d = 1.0_dp / d
   del = d * c
   h = h * del
   if (abs(del - 1.0_dp) <= 3.0e-14_dp) exit
end do
end function beta_cont_frac

subroutine print_lm_coef_rstyle(fit, term_names)
! Print coefficients with R-like aligned header/value rows.
type(lm_fit_t), intent(in) :: fit ! fitted model or test result
character(len=*), intent(in), optional :: term_names(:) ! predictor names excluding intercept
integer :: j, p
character(len=32) :: lbl
p = max(0, size(fit%coef) - 1)
write(*,'(a14)', advance='no') "(Intercept)"
do j = 1, p
   if (present(term_names) .and. size(term_names) >= j) then
      write(*,'(a14)', advance='no') trim(term_names(j))
   else
      write(lbl,'(a,i0)') "x", j
      write(*,'(a14)', advance='no') trim(lbl)
   end if
end do
write(*,*)
if (size(fit%coef) > 0) then
   write(*,'(*(f14.7))') fit%coef
else
   write(*,*)
end if
end subroutine print_lm_coef_rstyle

subroutine print_lm_confint(fit, term_names)
! Print coefficient confidence intervals with R-like row and column labels.
type(lm_fit_t), intent(in) :: fit ! fitted model or test result
character(len=*), intent(in), optional :: term_names(:) ! predictor names excluding intercept
real(kind=dp), allocatable :: ci(:,:)
integer :: j, p
character(len=32) :: lbl
ci = lm_confint(fit)
p = size(ci, 1)
write(*,'(13x,a12,1x,a12)') "2.5 %", "97.5 %"
do j = 1, p
   if (j == 1) then
      lbl = "(Intercept)"
   else if (present(term_names) .and. size(term_names) >= j - 1) then
      lbl = trim(term_names(j - 1))
   else
      write(lbl,'(a,i0)') "x", j - 1
   end if
   write(*,'(a12,1x,f12.7,1x,f12.7)') trim(lbl), ci(j, 1), ci(j, 2)
end do
end subroutine print_lm_confint

subroutine print_lm_anova(fit, term_names, term_df)
! Print a sequential analysis-of-variance table for an lm fit.
type(lm_fit_t), intent(in) :: fit ! fitted model or test result
character(len=*), intent(in), optional :: term_names(:) ! predictor names excluding intercept
integer, intent(in), optional :: term_df(:) ! input vector
type(lm_fit_t) :: fit_j
real(kind=dp) :: rss, rss_prev, ss_term, ms_term, ms_resid, fval, pval, ybar
integer :: j, p, nterms, df_term, first_col, last_col
character(len=32) :: lbl
rss = sum(fit%resid**2)
if (fit%df > 0) then
   ms_resid = rss / real(fit%df, kind=dp)
else
   ms_resid = ieee_value(0.0_dp, ieee_quiet_nan)
end if
write(*,'(a)') "Analysis of Variance Table"
write(*,'(a12,1x,a6,1x,a12,1x,a12,1x,a10,1x,a10)') "", "Df", "Sum Sq", "Mean Sq", "F value", "Pr(>F)"
if (allocated(fit%y) .and. allocated(fit%xpred)) then
   p = size(fit%xpred, 2)
   if (present(term_df)) then
      nterms = size(term_df)
   else
      nterms = p
   end if
   ybar = sum(fit%y) / real(size(fit%y), kind=dp)
   rss_prev = sum((fit%y - ybar)**2)
   first_col = 1
   do j = 1, nterms
      if (present(term_df)) then
         df_term = max(1, term_df(j))
      else
         df_term = 1
      end if
      last_col = min(p, first_col + df_term - 1)
      if (last_col < first_col) exit
      fit_j = lm_fit_general(fit%y, fit%xpred(:, 1:last_col))
      ss_term = max(0.0_dp, rss_prev - sum(fit_j%resid**2))
      ms_term = ss_term / real(max(1, df_term), kind=dp)
      if (ms_resid > 0.0_dp) then
         fval = ms_term / ms_resid
         pval = f_upper_tail_approx(fval, real(df_term, kind=dp), real(fit%df, kind=dp))
      else
         fval = ieee_value(0.0_dp, ieee_quiet_nan)
         pval = ieee_value(0.0_dp, ieee_quiet_nan)
      end if
      if (present(term_names) .and. size(term_names) >= j) then
         lbl = trim(term_names(j))
      else
         write(lbl,'(a,i0)') "x", j
      end if
      write(*,'(a12,1x,i6,1x,f12.4,1x,f12.4,1x,f10.2,1x,es10.3)') &
         & trim(lbl), df_term, ss_term, ms_term, fval, pval
      rss_prev = sum(fit_j%resid**2)
      first_col = last_col + 1
   end do
end if
write(*,'(a12,1x,i6,1x,f12.4,1x,f12.4)') "Residuals", fit%df, rss, ms_resid
end subroutine print_lm_anova

pure function f_upper_tail_approx(f, df1, df2) result(p)
! Runtime helper for R-compatible f upper tail approx.
real(kind=dp), intent(in) :: f, df1, df2
real(kind=dp) :: p, x
if (f <= 0.0_dp) then
   p = 1.0_dp
   return
end if
if (df1 <= 0.0_dp .or. df2 <= 0.0_dp) then
   p = ieee_value(0.0_dp, ieee_quiet_nan)
   return
end if
x = df2 / (df2 + df1 * f)
p = max(0.0_dp, min(1.0_dp, regularized_beta(x, 0.5_dp * df2, 0.5_dp * df1)))
end function f_upper_tail_approx

pure elemental function normal_cdf(x) result(p)
! Evaluate numerical special-function helper normal_cdf.
real(kind=dp), intent(in) :: x
real(kind=dp) :: p
p = 0.5_dp * (1.0_dp + erf(x / sqrt(2.0_dp)))
end function normal_cdf

pure function chisq_upper_tail_approx(x, df) result(p)
! Runtime helper for R-compatible chisq upper tail approx.
real(kind=dp), intent(in) :: x, df
real(kind=dp) :: p, z
if (df <= 0.0_dp) then
   p = 1.0_dp
   return
end if
if (x <= 0.0_dp) then
   p = 1.0_dp
   return
end if
if (abs(df - 1.0_dp) <= epsilon(1.0_dp)) then
   p = erfc(sqrt(0.5_dp * x))
   return
end if
z = ((x / df)**(1.0_dp / 3.0_dp) - (1.0_dp - 2.0_dp / (9.0_dp * df))) / &
   & sqrt(2.0_dp / (9.0_dp * df))
p = max(0.0_dp, min(1.0_dp, 1.0_dp - normal_cdf(z)))
end function chisq_upper_tail_approx

pure function pchisq_scalar(x, df, ncp) result(p)
! Lower-tail chi-square probability.
real(kind=dp), intent(in) :: x, df
real(kind=dp), intent(in), optional :: ncp
real(kind=dp) :: p
real(kind=dp) :: half_ncp, weight, sumw, term
integer :: j
if (x /= x .or. df <= 0.0_dp) then
   p = ieee_value(0.0_dp, ieee_quiet_nan)
   return
end if
if (present(ncp)) then
   if (ncp /= ncp .or. ncp < 0.0_dp .or. .not. ieee_is_finite(ncp)) then
      p = ieee_value(0.0_dp, ieee_quiet_nan)
      return
   end if
   if (ncp > 0.0_dp) then
      if (x <= 0.0_dp) then
         p = 0.0_dp
         return
      else if (.not. ieee_is_finite(x)) then
         p = 1.0_dp
         return
      end if
      half_ncp = 0.5_dp * ncp
      weight = exp(-half_ncp)
      p = weight * (1.0_dp - chisq_upper_tail_approx(x, df))
      sumw = weight
      do j = 1, 10000
         weight = weight * half_ncp / real(j, kind=dp)
         term = weight * (1.0_dp - chisq_upper_tail_approx(x, df + 2.0_dp * real(j, kind=dp)))
         p = p + term
         sumw = sumw + weight
         if (abs(term) <= 1.0e-13_dp * max(1.0_dp, abs(p)) .and. abs(1.0_dp - sumw) <= 1.0e-12_dp) exit
      end do
      p = max(0.0_dp, min(1.0_dp, p))
      return
   end if
end if
p = 1.0_dp - chisq_upper_tail_approx(x, df)
p = max(0.0_dp, min(1.0_dp, p))
end function pchisq_scalar

pure function qnorm_vec(pv, mean, sd, lower_tail) result(out)
! Approximate normal quantiles by bisection.
! Defaults: mean=0, sd=1, lower_tail=.true.
real(kind=dp), intent(in) :: pv(:) ! probabilities, clamped to [0,1]
real(kind=dp), intent(in), optional :: mean ! distribution mean
real(kind=dp), intent(in), optional :: sd ! distribution standard deviation
logical, intent(in), optional :: lower_tail ! lower-tail probability flag
real(kind=dp), allocatable :: out(:)
real(kind=dp) :: mu, sig, prob, lo, hi, mid
logical :: lower
integer :: i, iter
mu = 0.0_dp
sig = 1.0_dp
lower = .true.
if (present(mean)) mu = mean
if (present(sd)) sig = sd
if (present(lower_tail)) lower = lower_tail
allocate(out(size(pv)))
do i = 1, size(pv)
   if (pv(i) /= pv(i) .or. pv(i) < 0.0_dp .or. pv(i) > 1.0_dp .or. sig < 0.0_dp) then
      out(i) = ieee_value(0.0_dp, ieee_quiet_nan)
      cycle
   end if
   prob = pv(i)
   if (.not. lower) prob = 1.0_dp - prob
   if (prob <= 0.0_dp) then
      out(i) = -ieee_value(0.0_dp, ieee_positive_inf)
   else if (prob >= 1.0_dp) then
      out(i) = ieee_value(0.0_dp, ieee_positive_inf)
   else
      lo = -8.0_dp
      hi = 8.0_dp
      do iter = 1, 80
         mid = 0.5_dp * (lo + hi)
         if (normal_cdf(mid) < prob) then
            lo = mid
         else
            hi = mid
         end if
      end do
      out(i) = mu + sig * 0.5_dp * (lo + hi)
   end if
end do
end function qnorm_vec

pure function ppois_vec(q, lambda, lower_tail) result(out)
! Lower-tail Poisson CDF for vector q.
! Defaults: lambda=1, lower_tail=.true.
real(kind=dp), intent(in) :: q(:) ! quantiles
real(kind=dp), intent(in), optional :: lambda ! rate parameter
logical, intent(in), optional :: lower_tail ! lower-tail probability flag
real(kind=dp), allocatable :: out(:)
real(kind=dp) :: lam, term, cdf
logical :: lower
integer :: i, k, kk
lam = 1.0_dp
lower = .true.
if (present(lambda)) lam = lambda
if (present(lower_tail)) lower = lower_tail
allocate(out(size(q)))
do i = 1, size(q)
   if (lam < 0.0_dp) then
      out(i) = ieee_value(0.0_dp, ieee_quiet_nan)
   else if (q(i) /= q(i)) then
      out(i) = q(i)
   else if (q(i) < 0.0_dp) then
      out(i) = merge(0.0_dp, 1.0_dp, lower)
   else if (.not. ieee_is_finite(q(i))) then
      out(i) = merge(1.0_dp, 0.0_dp, lower)
   else
      kk = int(floor(q(i)))
      term = exp(-lam)
      cdf = term
      do k = 1, kk
         term = term * lam / real(k, kind=dp)
         cdf = cdf + term
      end do
      cdf = max(0.0_dp, min(1.0_dp, cdf))
      out(i) = merge(cdf, 1.0_dp - cdf, lower)
   end if
end do
end function ppois_vec

pure function qpois_vec(pv, lambda, lower_tail) result(out)
! Poisson quantiles as real values for matrix/vector composition.
! Defaults: lambda=1, lower_tail=.true.
real(kind=dp), intent(in) :: pv(:) ! probabilities, clamped to [0,1]
real(kind=dp), intent(in), optional :: lambda ! rate parameter
logical, intent(in), optional :: lower_tail ! lower-tail probability flag
real(kind=dp), allocatable :: out(:)
real(kind=dp) :: lam, prob, term, cdf
logical :: lower
integer :: i, k
lam = 1.0_dp
lower = .true.
if (present(lambda)) lam = lambda
if (present(lower_tail)) lower = lower_tail
allocate(out(size(pv)))
do i = 1, size(pv)
   if (lam < 0.0_dp .or. pv(i) /= pv(i) .or. pv(i) < 0.0_dp .or. pv(i) > 1.0_dp) then
      out(i) = ieee_value(0.0_dp, ieee_quiet_nan)
      cycle
   end if
   prob = pv(i)
   if (.not. lower) prob = 1.0_dp - prob
   if (prob >= 1.0_dp) then
      out(i) = ieee_value(0.0_dp, ieee_positive_inf)
   else if (prob <= 0.0_dp) then
      out(i) = 0.0_dp
   else
      term = exp(-lam)
      cdf = term
      k = 0
      do while (cdf < prob .and. k < 100000)
         k = k + 1
         term = term * lam / real(k, kind=dp)
         cdf = cdf + term
      end do
      out(i) = real(k, kind=dp)
   end if
end do
end function qpois_vec

pure elemental function r_choose_real(n, k) result(out)
! Runtime helper for R-compatible r choose real.
real(kind=dp), intent(in) :: n, k
real(kind=dp) :: out
if (k < 0.0_dp .or. k > n) then
   out = 0.0_dp
else
   out = exp(log_gamma(n + 1.0_dp) - log_gamma(k + 1.0_dp) - log_gamma(n - k + 1.0_dp))
end if
end function r_choose_real

pure elemental function gamma_p(a, x) result(gp)
! Runtime helper for R-compatible gamma p.
real(kind=dp), intent(in) :: a, x
real(kind=dp) :: gp, gln, ap, del, sumv, b, c, d, h, an
integer :: n
if (a <= 0.0_dp .or. x < 0.0_dp) then
   gp = ieee_value(0.0_dp, ieee_quiet_nan)
   return
end if
if (x == 0.0_dp) then
   gp = 0.0_dp
   return
end if
gln = log_gamma(a)
if (x < a + 1.0_dp) then
   ap = a
   sumv = 1.0_dp / a
   del = sumv
   do n = 1, 200
      ap = ap + 1.0_dp
      del = del * x / ap
      sumv = sumv + del
      if (abs(del) < abs(sumv) * 3.0e-14_dp) exit
   end do
   gp = sumv * exp(-x + a * log(x) - gln)
else
   b = x + 1.0_dp - a
   c = 1.0e30_dp
   d = 1.0_dp / b
   h = d
   do n = 1, 200
      an = -real(n, kind=dp) * (real(n, kind=dp) - a)
      b = b + 2.0_dp
      d = an * d + b
      if (abs(d) < 1.0e-30_dp) d = 1.0e-30_dp
      c = b + an / c
      if (abs(c) < 1.0e-30_dp) c = 1.0e-30_dp
      d = 1.0_dp / d
      del = d * c
      h = h * del
      if (abs(del - 1.0_dp) < 3.0e-14_dp) exit
   end do
   gp = 1.0_dp - exp(-x + a * log(x) - gln) * h
end if
gp = max(0.0_dp, min(1.0_dp, gp))
end function gamma_p

pure subroutine maybe_log_density(out, log_)
real(kind=dp), intent(inout) :: out(:)
logical, intent(in), optional :: log_
logical :: l
l = .false.
if (present(log_)) l = log_
if (l) out = log(out)
end subroutine maybe_log_density

pure function dunif_vec(x, min, max, log_) result(out)
! Evaluate distribution helper dunif.
real(kind=dp), intent(in) :: x(:) ! input vector
real(kind=dp), intent(in), optional :: min, max
logical, intent(in), optional :: log_
real(kind=dp), allocatable :: out(:)
real(kind=dp) :: lo, hi
integer :: i
lo = 0.0_dp; hi = 1.0_dp
if (present(min)) lo = min
if (present(max)) hi = max
allocate(out(size(x)))
if (hi <= lo) then
   out = ieee_value(0.0_dp, ieee_quiet_nan)
else
   do i = 1, size(x)
      if (x(i) /= x(i)) then
         out(i) = ieee_value(0.0_dp, ieee_quiet_nan)
      else if (x(i) >= lo .and. x(i) <= hi) then
         out(i) = 1.0_dp / (hi - lo)
      else
         out(i) = 0.0_dp
      end if
   end do
end if
call maybe_log_density(out, log_)
end function dunif_vec

pure function punif_vec(x, min, max) result(out)
! Evaluate distribution helper punif.
real(kind=dp), intent(in) :: x(:) ! input vector
real(kind=dp), intent(in), optional :: min, max
real(kind=dp), allocatable :: out(:)
real(kind=dp) :: lo, hi
lo = 0.0_dp; hi = 1.0_dp
if (present(min)) lo = min
if (present(max)) hi = max
allocate(out(size(x)))
if (hi <= lo) then
   out = ieee_value(0.0_dp, ieee_quiet_nan)
else
   out = (x - lo) / (hi - lo)
   where (out < 0.0_dp) out = 0.0_dp
   where (out > 1.0_dp) out = 1.0_dp
end if
end function punif_vec

pure function qunif_vec(p, min, max) result(out)
! Evaluate distribution helper qunif.
real(kind=dp), intent(in) :: p(:) ! probability value
real(kind=dp), intent(in), optional :: min, max
real(kind=dp), allocatable :: out(:)
real(kind=dp) :: lo, hi
integer :: i
lo = 0.0_dp; hi = 1.0_dp
if (present(min)) lo = min
if (present(max)) hi = max
allocate(out(size(p)))
do i = 1, size(p)
   if (p(i) /= p(i) .or. p(i) < 0.0_dp .or. p(i) > 1.0_dp .or. hi < lo) then
      out(i) = ieee_value(0.0_dp, ieee_quiet_nan)
   else
      out(i) = lo + p(i) * (hi - lo)
   end if
end do
end function qunif_vec

pure function dexp_vec(x, rate, log_) result(out)
! Evaluate distribution helper dexp.
real(kind=dp), intent(in) :: x(:) ! input vector
real(kind=dp), intent(in), optional :: rate ! rate parameter
logical, intent(in), optional :: log_
real(kind=dp), allocatable :: out(:)
real(kind=dp) :: r
logical :: l
integer :: i
r = 1.0_dp; if (present(rate)) r = rate
l = .false.; if (present(log_)) l = log_
allocate(out(size(x)))
do i = 1, size(x)
   if (x(i) /= x(i) .or. r < 0.0_dp) then
      out(i) = ieee_value(0.0_dp, ieee_quiet_nan)
   else if (x(i) < 0.0_dp) then
      out(i) = merge(-ieee_value(0.0_dp, ieee_positive_inf), 0.0_dp, l)
   else if (l) then
      out(i) = log(r) - r * x(i)
   else
      out(i) = r * exp(-r * x(i))
   end if
end do
end function dexp_vec

pure function pexp_vec(x, rate) result(out)
! Evaluate distribution helper pexp.
real(kind=dp), intent(in) :: x(:) ! input vector
real(kind=dp), intent(in), optional :: rate
real(kind=dp), allocatable :: out(:)
real(kind=dp) :: r
integer :: i
r = 1.0_dp; if (present(rate)) r = rate
allocate(out(size(x)))
do i = 1, size(x)
   if (x(i) /= x(i) .or. r < 0.0_dp) then
      out(i) = ieee_value(0.0_dp, ieee_quiet_nan)
   else if (x(i) <= 0.0_dp) then
      out(i) = 0.0_dp
   else
      out(i) = 1.0_dp - exp(-r * x(i))
   end if
end do
end function pexp_vec

pure function qexp_vec(p, rate) result(out)
! Evaluate distribution helper qexp.
real(kind=dp), intent(in) :: p(:) ! probability value
real(kind=dp), intent(in), optional :: rate
real(kind=dp), allocatable :: out(:)
real(kind=dp) :: r
integer :: i
r = 1.0_dp; if (present(rate)) r = rate
allocate(out(size(p)))
do i = 1, size(p)
   if (p(i) /= p(i) .or. p(i) < 0.0_dp .or. p(i) > 1.0_dp .or. r <= 0.0_dp) then
      out(i) = ieee_value(0.0_dp, ieee_quiet_nan)
   else if (p(i) >= 1.0_dp) then
      out(i) = ieee_value(0.0_dp, ieee_positive_inf)
   else
      out(i) = -log(1.0_dp - p(i)) / r
   end if
end do
end function qexp_vec

pure function dgamma_vec(x, shape, rate, log_) result(out)
! Evaluate distribution helper dgamma.
real(kind=dp), intent(in) :: x(:) ! input vector
real(kind=dp), intent(in) :: shape ! shape parameter
real(kind=dp), intent(in), optional :: rate
logical, intent(in), optional :: log_
real(kind=dp), allocatable :: out(:)
real(kind=dp) :: r, logc
integer :: i
r = 1.0_dp; if (present(rate)) r = rate
allocate(out(size(x)))
if (shape <= 0.0_dp .or. r <= 0.0_dp) then
   out = ieee_value(0.0_dp, ieee_quiet_nan)
else
   logc = shape * log(r) - log_gamma(shape)
   do i = 1, size(x)
      if (x(i) /= x(i)) then
         out(i) = ieee_value(0.0_dp, ieee_quiet_nan)
      else if (x(i) < 0.0_dp .or. .not. ieee_is_finite(x(i))) then
         out(i) = 0.0_dp
      else if (x(i) == 0.0_dp .and. shape < 1.0_dp) then
         out(i) = ieee_value(0.0_dp, ieee_positive_inf)
      else if (x(i) == 0.0_dp .and. shape == 1.0_dp) then
         out(i) = r
      else if (x(i) == 0.0_dp) then
         out(i) = 0.0_dp
      else
         out(i) = exp(logc + (shape - 1.0_dp) * log(x(i)) - r * x(i))
      end if
   end do
end if
call maybe_log_density(out, log_)
end function dgamma_vec

pure function pgamma_vec(x, shape, rate) result(out)
! Evaluate distribution helper pgamma.
real(kind=dp), intent(in) :: x(:) ! input vector
real(kind=dp), intent(in) :: shape ! shape parameter
real(kind=dp), intent(in), optional :: rate
real(kind=dp), allocatable :: out(:)
real(kind=dp) :: r
integer :: i
r = 1.0_dp; if (present(rate)) r = rate
allocate(out(size(x)))
do i = 1, size(x)
   if (x(i) /= x(i) .or. shape <= 0.0_dp .or. r <= 0.0_dp) then
      out(i) = ieee_value(0.0_dp, ieee_quiet_nan)
   else if (x(i) <= 0.0_dp) then
      out(i) = 0.0_dp
   else if (.not. ieee_is_finite(x(i))) then
      out(i) = 1.0_dp
   else
      out(i) = gamma_p(shape, r * x(i))
   end if
end do
end function pgamma_vec

pure function qgamma_vec(p, shape, rate) result(out)
! Evaluate distribution helper qgamma.
real(kind=dp), intent(in) :: p(:) ! probability value
real(kind=dp), intent(in) :: shape ! shape parameter
real(kind=dp), intent(in), optional :: rate
real(kind=dp), allocatable :: out(:)
real(kind=dp) :: r, lo, hi, mid
integer :: i, it
r = 1.0_dp; if (present(rate)) r = rate
allocate(out(size(p)))
do i = 1, size(p)
   if (p(i) /= p(i) .or. p(i) < 0.0_dp .or. p(i) > 1.0_dp .or. &
      shape <= 0.0_dp .or. r <= 0.0_dp) then
      out(i) = ieee_value(0.0_dp, ieee_quiet_nan)
      cycle
   else if (p(i) == 0.0_dp) then
      out(i) = 0.0_dp
      cycle
   else if (p(i) == 1.0_dp) then
      out(i) = ieee_value(0.0_dp, ieee_positive_inf)
      cycle
   end if
   lo = 0.0_dp; hi = max(1.0_dp, shape / r * 10.0_dp)
   do while (gamma_p(shape, r * hi) < p(i))
      hi = hi * 2.0_dp
      if (hi > 1.0e12_dp) exit
   end do
   do it = 1, 70
      mid = 0.5_dp * (lo + hi)
      if (gamma_p(shape, r * mid) < p(i)) lo = mid
      if (gamma_p(shape, r * mid) >= p(i)) hi = mid
   end do
   out(i) = 0.5_dp * (lo + hi)
end do
end function qgamma_vec

pure function dbeta_vec(x, shape1, shape2, log_) result(out)
! Evaluate distribution helper dbeta.
real(kind=dp), intent(in) :: x(:) ! quantiles or observed values
real(kind=dp), intent(in) :: shape1 ! first shape parameter
real(kind=dp), intent(in) :: shape2 ! second shape parameter
logical, intent(in), optional :: log_
real(kind=dp), allocatable :: out(:)
real(kind=dp) :: logc
integer :: i
allocate(out(size(x)))
if (shape1 <= 0.0_dp .or. shape2 <= 0.0_dp) then
   out = ieee_value(0.0_dp, ieee_quiet_nan)
else
   logc = log_gamma(shape1 + shape2) - log_gamma(shape1) - log_gamma(shape2)
   do i = 1, size(x)
      if (x(i) /= x(i)) then
         out(i) = ieee_value(0.0_dp, ieee_quiet_nan)
      else if (x(i) < 0.0_dp .or. x(i) > 1.0_dp) then
         out(i) = 0.0_dp
      else if (x(i) == 0.0_dp .and. shape1 < 1.0_dp) then
         out(i) = ieee_value(0.0_dp, ieee_positive_inf)
      else if (x(i) == 0.0_dp .and. shape1 == 1.0_dp) then
         out(i) = shape2
      else if (x(i) == 0.0_dp) then
         out(i) = 0.0_dp
      else if (x(i) == 1.0_dp .and. shape2 < 1.0_dp) then
         out(i) = ieee_value(0.0_dp, ieee_positive_inf)
      else if (x(i) == 1.0_dp .and. shape2 == 1.0_dp) then
         out(i) = shape1
      else if (x(i) == 1.0_dp) then
         out(i) = 0.0_dp
      else
         out(i) = exp(logc + (shape1 - 1.0_dp) * log(x(i)) + &
            (shape2 - 1.0_dp) * log(1.0_dp - x(i)))
      end if
   end do
end if
call maybe_log_density(out, log_)
end function dbeta_vec

pure function pbeta_vec(x, shape1, shape2) result(out)
! Evaluate distribution helper pbeta.
real(kind=dp), intent(in) :: x(:) ! quantiles or observed values
real(kind=dp), intent(in) :: shape1 ! first shape parameter
real(kind=dp), intent(in) :: shape2 ! second shape parameter
real(kind=dp), allocatable :: out(:)
integer :: i
allocate(out(size(x)))
do i = 1, size(x)
   if (x(i) /= x(i) .or. shape1 <= 0.0_dp .or. shape2 <= 0.0_dp) then
      out(i) = ieee_value(0.0_dp, ieee_quiet_nan)
   else if (x(i) <= 0.0_dp) then
      out(i) = 0.0_dp
   else if (x(i) >= 1.0_dp) then
      out(i) = 1.0_dp
   else
      out(i) = regularized_beta(x(i), shape1, shape2)
   end if
end do
end function pbeta_vec

pure function qbeta_vec(p, shape1, shape2) result(out)
! Evaluate distribution helper qbeta.
real(kind=dp), intent(in) :: p(:) ! probabilities, clamped to [0,1]
real(kind=dp), intent(in) :: shape1 ! first shape parameter
real(kind=dp), intent(in) :: shape2 ! second shape parameter
real(kind=dp), allocatable :: out(:)
real(kind=dp) :: lo, hi, mid
integer :: i, it
allocate(out(size(p)))
do i = 1, size(p)
   if (p(i) /= p(i) .or. p(i) < 0.0_dp .or. p(i) > 1.0_dp .or. &
      shape1 <= 0.0_dp .or. shape2 <= 0.0_dp) then
      out(i) = ieee_value(0.0_dp, ieee_quiet_nan)
      cycle
   else if (p(i) == 0.0_dp) then
      out(i) = 0.0_dp
      cycle
   else if (p(i) == 1.0_dp) then
      out(i) = 1.0_dp
      cycle
   end if
   lo = 0.0_dp; hi = 1.0_dp
   do it = 1, 70
      mid = 0.5_dp * (lo + hi)
      if (regularized_beta(mid, shape1, shape2) < p(i)) lo = mid
      if (regularized_beta(mid, shape1, shape2) >= p(i)) hi = mid
   end do
   out(i) = 0.5_dp * (lo + hi)
end do
end function qbeta_vec

pure function dchisq_vec(x, df, log_) result(out)
! Evaluate distribution helper dchisq.
real(kind=dp), intent(in) :: x(:) ! quantiles or observed values
real(kind=dp), intent(in) :: df ! degrees of freedom
logical, intent(in), optional :: log_
real(kind=dp), allocatable :: out(:)
out = dgamma(x, 0.5_dp * df, rate=0.5_dp, log_=log_)
end function dchisq_vec

pure function pchisq_vec(x, df, ncp) result(out)
! Evaluate distribution helper pchisq_vec.
real(kind=dp), intent(in) :: x(:) ! quantiles or observed values
real(kind=dp), intent(in) :: df ! degrees of freedom
real(kind=dp), intent(in), optional :: ncp
real(kind=dp), allocatable :: out(:)
integer :: i
allocate(out(size(x)))
do i = 1, size(x)
   if (present(ncp)) then
      out(i) = pchisq_scalar(x(i), df, ncp=ncp)
   else
      out(i) = pchisq_scalar(x(i), df)
   end if
end do
end function pchisq_vec

pure function qchisq_vec(p, df) result(out)
! Evaluate distribution helper qchisq_vec.
real(kind=dp), intent(in) :: p(:) ! probabilities, clamped to [0,1]
real(kind=dp), intent(in) :: df ! degrees of freedom
real(kind=dp), allocatable :: out(:)
out = qgamma(p, 0.5_dp * df, rate=0.5_dp)
end function qchisq_vec

pure function qchisq_scalar(p, df) result(out)
! Evaluate distribution helper qchisq_scalar.
real(kind=dp), intent(in) :: p ! probabilities, clamped to [0,1]
real(kind=dp), intent(in) :: df ! degrees of freedom
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = qchisq_vec([p], df)
out = tmp(1)
end function qchisq_scalar

pure function qchisq_scalar_i(p, df) result(out)
! Evaluate distribution helper qchisq_scalar_i.
real(kind=dp), intent(in) :: p ! probability value
integer, intent(in) :: df ! degrees of freedom
real(kind=dp) :: out
out = qchisq_scalar(p, real(df, kind=dp))
end function qchisq_scalar_i

pure function dt_vec(x, df, log_) result(out)
! Evaluate distribution helper dt.
real(kind=dp), intent(in) :: x(:) ! quantiles or observed values
real(kind=dp), intent(in) :: df ! degrees of freedom
logical, intent(in), optional :: log_
real(kind=dp), allocatable :: out(:)
real(kind=dp) :: logc
logical :: l
allocate(out(size(x)))
l = .false.; if (present(log_)) l = log_
if (df <= 0.0_dp) then
   out = ieee_value(0.0_dp, ieee_quiet_nan)
   return
end if
logc = log_gamma(0.5_dp * (df + 1.0_dp)) - log_gamma(0.5_dp * df) - 0.5_dp * log(df * acos(-1.0_dp))
out = exp(logc - 0.5_dp * (df + 1.0_dp) * log(1.0_dp + x * x / df))
if (l) out = log(max(tiny(1.0_dp), out))
end function dt_vec

pure function pt_vec(x, df) result(out)
! Evaluate distribution helper pt.
real(kind=dp), intent(in) :: x(:) ! quantiles or observed values
real(kind=dp), intent(in) :: df ! degrees of freedom
real(kind=dp), allocatable :: out(:)
real(kind=dp) :: z
integer :: i
allocate(out(size(x)))
do i = 1, size(x)
   if (x(i) /= x(i) .or. df <= 0.0_dp) then
      out(i) = ieee_value(0.0_dp, ieee_quiet_nan)
      cycle
   else if (.not. ieee_is_finite(x(i))) then
      out(i) = merge(1.0_dp, 0.0_dp, x(i) > 0.0_dp)
      cycle
   end if
   z = df / (df + x(i) * x(i))
   if (x(i) >= 0.0_dp) then
      out(i) = 1.0_dp - 0.5_dp * regularized_beta(z, 0.5_dp * df, 0.5_dp)
   else
      out(i) = 0.5_dp * regularized_beta(z, 0.5_dp * df, 0.5_dp)
   end if
end do
end function pt_vec

pure function qt_vec(p, df) result(out)
! Evaluate distribution helper qt.
real(kind=dp), intent(in) :: p(:) ! probabilities, clamped to [0,1]
real(kind=dp), intent(in) :: df ! degrees of freedom
real(kind=dp), allocatable :: out(:)
real(kind=dp) :: lo, hi, mid, pmid(1)
integer :: i, it
allocate(out(size(p)))
do i = 1, size(p)
   if (p(i) /= p(i) .or. p(i) < 0.0_dp .or. p(i) > 1.0_dp .or. df <= 0.0_dp) then
      out(i) = ieee_value(0.0_dp, ieee_quiet_nan)
      cycle
   else if (p(i) == 0.0_dp) then
      out(i) = -ieee_value(0.0_dp, ieee_positive_inf)
      cycle
   else if (p(i) == 1.0_dp) then
      out(i) = ieee_value(0.0_dp, ieee_positive_inf)
      cycle
   else if (p(i) == 0.5_dp) then
      out(i) = 0.0_dp
      cycle
   end if
   lo = -32.0_dp; hi = 32.0_dp
   do it = 1, 80
      mid = 0.5_dp * (lo + hi)
      pmid = pt([mid], df)
      if (pmid(1) < p(i)) lo = mid
      if (pmid(1) >= p(i)) hi = mid
   end do
   out(i) = 0.5_dp * (lo + hi)
end do
end function qt_vec

pure function df_vec(x, df1, df2, log_) result(out)
! Evaluate distribution helper df.
real(kind=dp), intent(in) :: x(:) ! quantiles or observed values
real(kind=dp), intent(in) :: df1 ! numerator degrees of freedom
real(kind=dp), intent(in) :: df2 ! denominator degrees of freedom
logical, intent(in), optional :: log_
real(kind=dp), allocatable :: out(:)
real(kind=dp) :: a, b, logc
integer :: i
allocate(out(size(x)))
if (df1 <= 0.0_dp .or. df2 <= 0.0_dp) then
   out = ieee_value(0.0_dp, ieee_quiet_nan)
   return
end if
a = 0.5_dp * df1; b = 0.5_dp * df2
logc = a * log(df1 / df2) - (log_gamma(a) + log_gamma(b) - log_gamma(a + b))
do i = 1, size(x)
   if (x(i) /= x(i)) then
      out(i) = ieee_value(0.0_dp, ieee_quiet_nan)
   else if (x(i) < 0.0_dp .or. .not. ieee_is_finite(x(i))) then
      out(i) = 0.0_dp
   else if (x(i) == 0.0_dp .and. a < 1.0_dp) then
      out(i) = ieee_value(0.0_dp, ieee_positive_inf)
   else if (x(i) == 0.0_dp .and. a == 1.0_dp) then
      out(i) = 1.0_dp
   else if (x(i) == 0.0_dp) then
      out(i) = 0.0_dp
   else
      out(i) = exp(logc + (a - 1.0_dp) * log(x(i)) - &
         (a + b) * log(1.0_dp + (df1 / df2) * x(i)))
   end if
end do
call maybe_log_density(out, log_)
end function df_vec

pure function pf_vec(x, df1, df2) result(out)
! Evaluate distribution helper pf.
real(kind=dp), intent(in) :: x(:) ! quantiles or observed values
real(kind=dp), intent(in) :: df1 ! numerator degrees of freedom
real(kind=dp), intent(in) :: df2 ! denominator degrees of freedom
real(kind=dp), allocatable :: out(:)
real(kind=dp) :: z
integer :: i
allocate(out(size(x)))
do i = 1, size(x)
   if (x(i) /= x(i) .or. df1 <= 0.0_dp .or. df2 <= 0.0_dp) then
      out(i) = ieee_value(0.0_dp, ieee_quiet_nan)
   else if (x(i) <= 0.0_dp) then
      out(i) = 0.0_dp
   else if (.not. ieee_is_finite(x(i))) then
      out(i) = 1.0_dp
   else
      z = (df1 * x(i)) / (df1 * x(i) + df2)
      out(i) = regularized_beta(z, 0.5_dp * df1, 0.5_dp * df2)
   end if
end do
end function pf_vec

pure function qf_vec(p, df1, df2) result(out)
! Evaluate distribution helper qf.
real(kind=dp), intent(in) :: p(:) ! probabilities, clamped to [0,1]
real(kind=dp), intent(in) :: df1 ! numerator degrees of freedom
real(kind=dp), intent(in) :: df2 ! denominator degrees of freedom
real(kind=dp), allocatable :: out(:)
real(kind=dp) :: lo, hi, mid, tmp(1)
integer :: i, it
allocate(out(size(p)))
do i = 1, size(p)
   if (p(i) /= p(i) .or. p(i) < 0.0_dp .or. p(i) > 1.0_dp .or. &
      df1 <= 0.0_dp .or. df2 <= 0.0_dp) then
      out(i) = ieee_value(0.0_dp, ieee_quiet_nan)
      cycle
   else if (p(i) == 0.0_dp) then
      out(i) = 0.0_dp
      cycle
   else if (p(i) == 1.0_dp) then
      out(i) = ieee_value(0.0_dp, ieee_positive_inf)
      cycle
   end if
   lo = 0.0_dp; hi = 1.0_dp
   do
      tmp = pf([hi], df1, df2)
      if (tmp(1) >= p(i) .or. hi > 1.0e12_dp) exit
      hi = hi * 2.0_dp
   end do
   do it = 1, 70
      mid = 0.5_dp * (lo + hi)
      tmp = pf([mid], df1, df2)
      if (tmp(1) < p(i)) lo = mid
      if (tmp(1) >= p(i)) hi = mid
   end do
   out(i) = 0.5_dp * (lo + hi)
end do
end function qf_vec

pure function dlogis_vec(x, location, scale, log_) result(out)
! Evaluate distribution helper dlogis.
real(kind=dp), intent(in) :: x(:) ! input vector
real(kind=dp), intent(in), optional :: location ! location parameter
real(kind=dp), intent(in), optional :: scale ! scale parameter
logical, intent(in), optional :: log_
real(kind=dp), allocatable :: out(:)
real(kind=dp) :: loc, sc, z, ez
integer :: i
loc = 0.0_dp; sc = 1.0_dp
if (present(location)) loc = location
if (present(scale)) sc = scale
allocate(out(size(x)))
if (sc <= 0.0_dp) then
   out = ieee_value(0.0_dp, ieee_quiet_nan)
else
   do i = 1, size(x)
      z = abs((x(i) - loc) / sc)
      ez = exp(-z)
      out(i) = ez / (sc * (1.0_dp + ez)**2)
   end do
end if
call maybe_log_density(out, log_)
end function dlogis_vec

pure function plogis_vec(x, location, scale) result(out)
! Evaluate distribution helper plogis.
real(kind=dp), intent(in) :: x(:) ! input vector
real(kind=dp), intent(in), optional :: location ! location parameter
real(kind=dp), intent(in), optional :: scale ! scale parameter
real(kind=dp), allocatable :: out(:)
real(kind=dp) :: loc, sc, z, ez
integer :: i
loc = 0.0_dp; sc = 1.0_dp
if (present(location)) loc = location
if (present(scale)) sc = scale
allocate(out(size(x)))
if (sc <= 0.0_dp) then
   out = ieee_value(0.0_dp, ieee_quiet_nan)
else
   do i = 1, size(x)
      z = (x(i) - loc) / sc
      if (z >= 0.0_dp) then
         out(i) = 1.0_dp / (1.0_dp + exp(-z))
      else
         ez = exp(z)
         out(i) = ez / (1.0_dp + ez)
      end if
   end do
end if
end function plogis_vec

pure function qlogis_vec(p, location, scale) result(out)
! Evaluate distribution helper qlogis.
real(kind=dp), intent(in) :: p(:) ! probability value
real(kind=dp), intent(in), optional :: location ! location parameter
real(kind=dp), intent(in), optional :: scale ! scale parameter
real(kind=dp), allocatable :: out(:)
real(kind=dp) :: loc, sc
integer :: i
loc = 0.0_dp; sc = 1.0_dp
if (present(location)) loc = location
if (present(scale)) sc = scale
allocate(out(size(p)))
do i = 1, size(p)
   if (p(i) /= p(i) .or. p(i) < 0.0_dp .or. p(i) > 1.0_dp .or. sc <= 0.0_dp) then
      out(i) = ieee_value(0.0_dp, ieee_quiet_nan)
   else if (p(i) <= 0.0_dp) then
      out(i) = -ieee_value(0.0_dp, ieee_positive_inf)
   else if (p(i) >= 1.0_dp) then
      out(i) = ieee_value(0.0_dp, ieee_positive_inf)
   else
      out(i) = loc + sc * log(p(i) / (1.0_dp - p(i)))
   end if
end do
end function qlogis_vec

pure function dlnorm_vec(x, meanlog, sdlog, log_) result(out)
! Evaluate distribution helper dlnorm.
real(kind=dp), intent(in) :: x(:) ! input vector
real(kind=dp), intent(in), optional :: meanlog ! mean on log scale
real(kind=dp), intent(in), optional :: sdlog ! standard deviation on log scale
logical, intent(in), optional :: log_
real(kind=dp), allocatable :: out(:)
real(kind=dp) :: mu, sig
integer :: i
mu = 0.0_dp; sig = 1.0_dp
if (present(meanlog)) mu = meanlog
if (present(sdlog)) sig = sdlog
allocate(out(size(x)))
if (sig <= 0.0_dp) then
   out = ieee_value(0.0_dp, ieee_quiet_nan)
else
   do i = 1, size(x)
      if (x(i) /= x(i)) then
         out(i) = ieee_value(0.0_dp, ieee_quiet_nan)
      else if (x(i) <= 0.0_dp .or. .not. ieee_is_finite(x(i))) then
         out(i) = 0.0_dp
      else
         out(i) = exp(-log(x(i) * sig * sqrt(2.0_dp * acos(-1.0_dp))) - &
            0.5_dp * ((log(x(i)) - mu) / sig)**2)
      end if
   end do
end if
call maybe_log_density(out, log_)
end function dlnorm_vec

pure function plnorm_vec(x, meanlog, sdlog) result(out)
! Evaluate distribution helper plnorm.
real(kind=dp), intent(in) :: x(:) ! input vector
real(kind=dp), intent(in), optional :: meanlog ! mean on log scale
real(kind=dp), intent(in), optional :: sdlog ! standard deviation on log scale
real(kind=dp), allocatable :: out(:)
real(kind=dp) :: mu, sig
integer :: i
mu = 0.0_dp; sig = 1.0_dp
if (present(meanlog)) mu = meanlog
if (present(sdlog)) sig = sdlog
allocate(out(size(x)))
if (sig <= 0.0_dp) then
   out = ieee_value(0.0_dp, ieee_quiet_nan)
else
   do i = 1, size(x)
      if (x(i) /= x(i)) then
         out(i) = ieee_value(0.0_dp, ieee_quiet_nan)
      else if (x(i) <= 0.0_dp) then
         out(i) = 0.0_dp
      else
         out(i) = normal_cdf((log(x(i)) - mu) / sig)
      end if
   end do
end if
end function plnorm_vec

pure function qlnorm_vec(p, meanlog, sdlog) result(out)
! Evaluate distribution helper qlnorm.
real(kind=dp), intent(in) :: p(:) ! probability value
real(kind=dp), intent(in), optional :: meanlog ! mean on log scale
real(kind=dp), intent(in), optional :: sdlog ! standard deviation on log scale
real(kind=dp), allocatable :: out(:)
real(kind=dp) :: mu, sig
mu = 0.0_dp; sig = 1.0_dp
if (present(meanlog)) mu = meanlog
if (present(sdlog)) sig = sdlog
out = exp(qnorm(p, mean=mu, sd=sig))
end function qlnorm_vec

pure function dweibull_vec(x, shape, scale, log_) result(out)
! Evaluate distribution helper dweibull.
real(kind=dp), intent(in) :: x(:) ! input vector
real(kind=dp), intent(in) :: shape ! shape parameter
real(kind=dp), intent(in), optional :: scale ! scale parameter
logical, intent(in), optional :: log_
real(kind=dp), allocatable :: out(:)
real(kind=dp) :: sc
integer :: i
sc = 1.0_dp; if (present(scale)) sc = scale
allocate(out(size(x)))
if (shape <= 0.0_dp .or. sc <= 0.0_dp) then
   out = ieee_value(0.0_dp, ieee_quiet_nan)
else
   do i = 1, size(x)
      if (x(i) /= x(i)) then
         out(i) = ieee_value(0.0_dp, ieee_quiet_nan)
      else if (x(i) < 0.0_dp .or. (x(i) > 0.0_dp .and. .not. ieee_is_finite(x(i)))) then
         out(i) = 0.0_dp
      else if (x(i) == 0.0_dp .and. shape < 1.0_dp) then
         out(i) = ieee_value(0.0_dp, ieee_positive_inf)
      else if (x(i) == 0.0_dp .and. shape == 1.0_dp) then
         out(i) = 1.0_dp / sc
      else if (x(i) == 0.0_dp) then
         out(i) = 0.0_dp
      else
         out(i) = (shape / sc) * (x(i) / sc)**(shape - 1.0_dp) * exp(-(x(i) / sc)**shape)
      end if
   end do
end if
call maybe_log_density(out, log_)
end function dweibull_vec

pure function pweibull_vec(x, shape, scale) result(out)
! Evaluate distribution helper pweibull.
real(kind=dp), intent(in) :: x(:) ! input vector
real(kind=dp), intent(in) :: shape ! shape parameter
real(kind=dp), intent(in), optional :: scale ! scale parameter
real(kind=dp), allocatable :: out(:)
real(kind=dp) :: sc
integer :: i
sc = 1.0_dp; if (present(scale)) sc = scale
allocate(out(size(x)))
if (shape <= 0.0_dp .or. sc <= 0.0_dp) then
   out = ieee_value(0.0_dp, ieee_quiet_nan)
else
   do i = 1, size(x)
      if (x(i) /= x(i)) then
         out(i) = ieee_value(0.0_dp, ieee_quiet_nan)
      else if (x(i) <= 0.0_dp) then
         out(i) = 0.0_dp
      else
         out(i) = 1.0_dp - exp(-(x(i) / sc)**shape)
      end if
   end do
end if
end function pweibull_vec

pure function qweibull_vec(p, shape, scale) result(out)
! Evaluate distribution helper qweibull.
real(kind=dp), intent(in) :: p(:) ! probability value
real(kind=dp), intent(in) :: shape ! shape parameter
real(kind=dp), intent(in), optional :: scale ! scale parameter
real(kind=dp), allocatable :: out(:)
real(kind=dp) :: sc
integer :: i
sc = 1.0_dp; if (present(scale)) sc = scale
allocate(out(size(p)))
do i = 1, size(p)
   if (p(i) /= p(i) .or. p(i) < 0.0_dp .or. p(i) > 1.0_dp .or. &
      shape <= 0.0_dp .or. sc <= 0.0_dp) then
      out(i) = ieee_value(0.0_dp, ieee_quiet_nan)
   else if (p(i) >= 1.0_dp) then
      out(i) = ieee_value(0.0_dp, ieee_positive_inf)
   else
      out(i) = sc * (-log(1.0_dp - p(i)))**(1.0_dp / shape)
   end if
end do
end function qweibull_vec

pure function dcauchy_vec(x, location, scale, log_) result(out)
! Evaluate distribution helper dcauchy.
real(kind=dp), intent(in) :: x(:) ! input vector
real(kind=dp), intent(in), optional :: location ! location parameter
real(kind=dp), intent(in), optional :: scale ! scale parameter
logical, intent(in), optional :: log_
real(kind=dp), allocatable :: out(:)
real(kind=dp) :: loc, sc, z(size(x))
loc = 0.0_dp; sc = 1.0_dp
if (present(location)) loc = location
if (present(scale)) sc = scale
allocate(out(size(x)))
z = (x - loc) / sc
if (sc <= 0.0_dp) then
   out = ieee_value(0.0_dp, ieee_quiet_nan)
else
   out = 1.0_dp / (acos(-1.0_dp) * sc * (1.0_dp + z * z))
end if
call maybe_log_density(out, log_)
end function dcauchy_vec

pure function pcauchy_vec(x, location, scale) result(out)
! Evaluate distribution helper pcauchy.
real(kind=dp), intent(in) :: x(:) ! input vector
real(kind=dp), intent(in), optional :: location ! location parameter
real(kind=dp), intent(in), optional :: scale ! scale parameter
real(kind=dp), allocatable :: out(:)
real(kind=dp) :: loc, sc
loc = 0.0_dp; sc = 1.0_dp
if (present(location)) loc = location
if (present(scale)) sc = scale
allocate(out(size(x)))
if (sc <= 0.0_dp) then
   out = ieee_value(0.0_dp, ieee_quiet_nan)
else
   out = 0.5_dp + atan((x - loc) / sc) / acos(-1.0_dp)
end if
end function pcauchy_vec

pure function qcauchy_vec(p, location, scale) result(out)
! Evaluate distribution helper qcauchy.
real(kind=dp), intent(in) :: p(:) ! probability value
real(kind=dp), intent(in), optional :: location ! location parameter
real(kind=dp), intent(in), optional :: scale ! scale parameter
real(kind=dp), allocatable :: out(:)
real(kind=dp) :: loc, sc
integer :: i
loc = 0.0_dp; sc = 1.0_dp
if (present(location)) loc = location
if (present(scale)) sc = scale
allocate(out(size(p)))
do i = 1, size(p)
   if (p(i) /= p(i) .or. p(i) < 0.0_dp .or. p(i) > 1.0_dp .or. sc <= 0.0_dp) then
      out(i) = ieee_value(0.0_dp, ieee_quiet_nan)
   else if (p(i) <= 0.0_dp) then
      out(i) = -ieee_value(0.0_dp, ieee_positive_inf)
   else if (p(i) >= 1.0_dp) then
      out(i) = ieee_value(0.0_dp, ieee_positive_inf)
   else
      out(i) = loc + sc * tan(acos(-1.0_dp) * (p(i) - 0.5_dp))
   end if
end do
end function qcauchy_vec

pure function dbinom_vec(x, nsize, prob, log_) result(out)
! Evaluate distribution helper dbinom.
real(kind=dp), intent(in) :: x(:) ! input vector
integer, intent(in) :: nsize ! input value
real(kind=dp), intent(in) :: prob ! success probability
logical, intent(in), optional :: log_
real(kind=dp), allocatable :: out(:)
integer :: i, k
allocate(out(size(x)))
do i = 1, size(x)
   if (.not. ieee_is_finite(x(i))) then
      out(i) = merge(x(i), 0.0_dp, x(i) /= x(i))
      cycle
   end if
   k = int(x(i))
   if (x(i) /= real(k, kind=dp) .or. k < 0 .or. k > nsize .or. &
      prob < 0.0_dp .or. prob > 1.0_dp) then
      out(i) = 0.0_dp
   else
      out(i) = r_choose_real(real(nsize, kind=dp), real(k, kind=dp)) * prob**k * (1.0_dp - prob)**(nsize - k)
   end if
end do
call maybe_log_density(out, log_)
end function dbinom_vec

pure function pbinom_vec(x, nsize, prob) result(out)
! Evaluate distribution helper pbinom.
real(kind=dp), intent(in) :: x(:) ! input vector
integer, intent(in) :: nsize ! input value
real(kind=dp), intent(in) :: prob ! success probability
real(kind=dp), allocatable :: out(:)
integer :: i, k
allocate(out(size(x)))
do i = 1, size(x)
   if (x(i) /= x(i)) then
      out(i) = x(i)
   else if (x(i) < 0.0_dp) then
      out(i) = 0.0_dp
   else if (.not. ieee_is_finite(x(i))) then
      out(i) = 1.0_dp
   else
      out(i) = 0.0_dp
      do k = 0, min(nsize, int(floor(x(i))))
         out(i) = out(i) + r_choose_real(real(nsize, kind=dp), real(k, kind=dp)) * &
            prob**k * (1.0_dp - prob)**(nsize - k)
      end do
   end if
end do
end function pbinom_vec

pure function qbinom_vec(p, nsize, prob) result(out)
! Evaluate distribution helper qbinom.
real(kind=dp), intent(in) :: p(:) ! probability value
integer, intent(in) :: nsize ! input value
real(kind=dp), intent(in) :: prob ! success probability
real(kind=dp), allocatable :: out(:)
integer :: i, k
real(kind=dp) :: cdf
allocate(out(size(p)))
do i = 1, size(p)
   if (p(i) /= p(i) .or. p(i) < 0.0_dp .or. p(i) > 1.0_dp) then
      out(i) = ieee_value(0.0_dp, ieee_quiet_nan)
   else if (p(i) <= 0.0_dp) then
      out(i) = 0.0_dp
   else if (p(i) >= 1.0_dp) then
      out(i) = real(nsize, kind=dp)
   else
      cdf = 0.0_dp
      do k = 0, nsize
         cdf = cdf + r_choose_real(real(nsize, kind=dp), real(k, kind=dp)) * &
            prob**k * (1.0_dp - prob)**(nsize - k)
         if (cdf >= p(i)) exit
      end do
      out(i) = real(k, kind=dp)
   end if
end do
end function qbinom_vec

pure function dpois_vec(x, lambda, log_) result(out)
! Evaluate distribution helper dpois.
real(kind=dp), intent(in) :: x(:) ! input vector
real(kind=dp), intent(in), optional :: lambda ! rate parameter
logical, intent(in), optional :: log_
real(kind=dp), allocatable :: out(:)
real(kind=dp) :: lam
integer :: i, k
lam = 1.0_dp; if (present(lambda)) lam = lambda
allocate(out(size(x)))
do i = 1, size(x)
   if (.not. ieee_is_finite(x(i))) then
      out(i) = merge(x(i), 0.0_dp, x(i) /= x(i))
      cycle
   end if
   k = int(x(i))
   if (x(i) /= real(k, kind=dp) .or. k < 0) then
      out(i) = 0.0_dp
   else
      out(i) = exp(-lam + real(k, kind=dp) * log(lam) - log_gamma(real(k + 1, kind=dp)))
   end if
end do
call maybe_log_density(out, log_)
end function dpois_vec

pure function dgeom_vec(x, prob, log_) result(out)
! Evaluate distribution helper dgeom.
real(kind=dp), intent(in) :: x(:) ! quantiles or observed values
real(kind=dp), intent(in) :: prob ! success probability
logical, intent(in), optional :: log_
real(kind=dp), allocatable :: out(:)
integer :: i, k
allocate(out(size(x)))
do i = 1, size(x)
   if (.not. ieee_is_finite(x(i))) then
      out(i) = merge(x(i), 0.0_dp, x(i) /= x(i))
      cycle
   end if
   k = int(x(i))
   if (x(i) /= real(k, kind=dp) .or. k < 0) then
      out(i) = 0.0_dp
   else
      out(i) = prob * (1.0_dp - prob)**k
   end if
end do
call maybe_log_density(out, log_)
end function dgeom_vec

pure function pgeom_vec(x, prob) result(out)
! Evaluate distribution helper pgeom.
real(kind=dp), intent(in) :: x(:) ! quantiles or observed values
real(kind=dp), intent(in) :: prob ! success probability
real(kind=dp), allocatable :: out(:)
integer :: i
allocate(out(size(x)))
do i = 1, size(x)
   if (x(i) /= x(i)) then
      out(i) = x(i)
   else if (x(i) < 0.0_dp) then
      out(i) = 0.0_dp
   else if (.not. ieee_is_finite(x(i))) then
      out(i) = 1.0_dp
   else
      out(i) = 1.0_dp - (1.0_dp - prob)**(int(floor(x(i))) + 1)
   end if
end do
end function pgeom_vec

pure function qgeom_vec(p, prob) result(out)
! Evaluate distribution helper qgeom.
real(kind=dp), intent(in) :: p(:) ! probabilities, clamped to [0,1]
real(kind=dp), intent(in) :: prob ! success probability
real(kind=dp), allocatable :: out(:)
integer :: i
allocate(out(size(p)))
do i = 1, size(p)
   if (p(i) /= p(i) .or. p(i) < 0.0_dp .or. p(i) > 1.0_dp) then
      out(i) = ieee_value(0.0_dp, ieee_quiet_nan)
   else if (p(i) <= 0.0_dp) then
      out(i) = 0.0_dp
   else if (p(i) >= 1.0_dp) then
      out(i) = ieee_value(0.0_dp, ieee_positive_inf)
   else
      out(i) = max(0.0_dp, real(ceiling(log(1.0_dp - p(i)) / &
         log(1.0_dp - prob) - 1.0_dp), kind=dp))
   end if
end do
end function qgeom_vec

pure function dnbinom_vec(x, nsize, prob, log_) result(out)
! Evaluate distribution helper dnbinom.
real(kind=dp), intent(in) :: x(:) ! quantiles or observed values
real(kind=dp), intent(in) :: prob ! success probability
integer, intent(in) :: nsize ! input value
logical, intent(in), optional :: log_
real(kind=dp), allocatable :: out(:)
integer :: i, k
allocate(out(size(x)))
do i = 1, size(x)
   if (.not. ieee_is_finite(x(i))) then
      out(i) = merge(x(i), 0.0_dp, x(i) /= x(i))
      cycle
   end if
   k = int(x(i))
   if (x(i) /= real(k, kind=dp) .or. k < 0) then
      out(i) = 0.0_dp
   else
      out(i) = r_choose_real(real(k + nsize - 1, kind=dp), real(k, kind=dp)) * &
         prob**nsize * (1.0_dp - prob)**k
   end if
end do
call maybe_log_density(out, log_)
end function dnbinom_vec

pure function pnbinom_vec(x, nsize, prob) result(out)
! Evaluate distribution helper pnbinom.
real(kind=dp), intent(in) :: x(:) ! quantiles or observed values
real(kind=dp), intent(in) :: prob ! success probability
integer, intent(in) :: nsize ! input value
real(kind=dp), allocatable :: out(:)
integer :: i, k
allocate(out(size(x)))
do i = 1, size(x)
   if (x(i) /= x(i)) then
      out(i) = x(i)
   else if (x(i) < 0.0_dp) then
      out(i) = 0.0_dp
   else if (.not. ieee_is_finite(x(i))) then
      out(i) = 1.0_dp
   else
      out(i) = 0.0_dp
      do k = 0, int(floor(x(i)))
         out(i) = out(i) + r_choose_real(real(k + nsize - 1, kind=dp), &
            real(k, kind=dp)) * prob**nsize * (1.0_dp - prob)**k
      end do
   end if
end do
end function pnbinom_vec

pure function qnbinom_vec(p, nsize, prob) result(out)
! Evaluate distribution helper qnbinom.
real(kind=dp), intent(in) :: p(:) ! probabilities, clamped to [0,1]
real(kind=dp), intent(in) :: prob ! success probability
integer, intent(in) :: nsize ! input value
real(kind=dp), allocatable :: out(:)
integer :: i, k
real(kind=dp) :: cdf
allocate(out(size(p)))
do i = 1, size(p)
   if (p(i) /= p(i) .or. p(i) < 0.0_dp .or. p(i) > 1.0_dp) then
      out(i) = ieee_value(0.0_dp, ieee_quiet_nan)
   else if (p(i) <= 0.0_dp) then
      out(i) = 0.0_dp
   else if (p(i) >= 1.0_dp) then
      out(i) = ieee_value(0.0_dp, ieee_positive_inf)
   else
      cdf = 0.0_dp; k = 0
      do while (cdf < p(i) .and. k < 100000)
         cdf = cdf + r_choose_real(real(k + nsize - 1, kind=dp), real(k, kind=dp)) * &
            prob**nsize * (1.0_dp - prob)**k
         if (cdf >= p(i)) exit
         k = k + 1
      end do
      out(i) = real(k, kind=dp)
   end if
end do
end function qnbinom_vec

pure function dhyper_vec(x, m, n, k, log_) result(out)
! Evaluate distribution helper dhyper.
real(kind=dp), intent(in) :: x(:) ! input vector
integer, intent(in) :: m ! white ball count
integer, intent(in) :: n ! black ball count or sample size
integer, intent(in) :: k ! draw count or group count
logical, intent(in), optional :: log_
real(kind=dp), allocatable :: out(:)
integer :: i, xx
real(kind=dp) :: den
allocate(out(size(x)))
den = r_choose_real(real(m + n, kind=dp), real(k, kind=dp))
do i = 1, size(x)
   if (.not. ieee_is_finite(x(i))) then
      out(i) = merge(x(i), 0.0_dp, x(i) /= x(i))
      cycle
   end if
   xx = int(x(i))
   if (x(i) /= real(xx, kind=dp)) then
      out(i) = 0.0_dp
   else
      out(i) = r_choose_real(real(m, kind=dp), real(xx, kind=dp)) * &
         & r_choose_real(real(n, kind=dp), real(k - xx, kind=dp)) / den
   end if
end do
call maybe_log_density(out, log_)
end function dhyper_vec

pure function phyper_vec(x, m, n, k) result(out)
! Evaluate distribution helper phyper.
real(kind=dp), intent(in) :: x(:) ! input vector
integer, intent(in) :: m ! white ball count
integer, intent(in) :: n ! black ball count or sample size
integer, intent(in) :: k ! draw count or group count
real(kind=dp), allocatable :: out(:)
integer :: i, xx
allocate(out(size(x)))
do i = 1, size(x)
   if (x(i) /= x(i)) then
      out(i) = x(i)
   else if (x(i) < 0.0_dp) then
      out(i) = 0.0_dp
   else if (.not. ieee_is_finite(x(i))) then
      out(i) = 1.0_dp
   else
      out(i) = 0.0_dp
      do xx = 0, int(floor(x(i)))
         out(i) = out(i) + r_choose_real(real(m, kind=dp), real(xx, kind=dp)) * &
            & r_choose_real(real(n, kind=dp), real(k - xx, kind=dp)) / &
            & r_choose_real(real(m + n, kind=dp), real(k, kind=dp))
      end do
   end if
end do
end function phyper_vec

pure function qhyper_vec(p, m, n, k) result(out)
! Evaluate distribution helper qhyper.
real(kind=dp), intent(in) :: p(:) ! probability value
integer, intent(in) :: m ! white ball count
integer, intent(in) :: n ! black ball count or sample size
integer, intent(in) :: k ! draw count or group count
real(kind=dp), allocatable :: out(:)
integer :: i, xx
real(kind=dp) :: cdf
allocate(out(size(p)))
do i = 1, size(p)
   if (p(i) /= p(i) .or. p(i) < 0.0_dp .or. p(i) > 1.0_dp) then
      out(i) = ieee_value(0.0_dp, ieee_quiet_nan)
   else if (p(i) <= 0.0_dp) then
      out(i) = real(max(0, k - n), kind=dp)
   else if (p(i) >= 1.0_dp) then
      out(i) = real(min(k, m), kind=dp)
   else
      cdf = 0.0_dp
      do xx = max(0, k - n), min(k, m)
         cdf = cdf + r_choose_real(real(m, kind=dp), real(xx, kind=dp)) * &
            & r_choose_real(real(n, kind=dp), real(k - xx, kind=dp)) / &
            & r_choose_real(real(m + n, kind=dp), real(k, kind=dp))
         if (cdf >= p(i)) exit
      end do
      out(i) = real(xx, kind=dp)
   end if
end do
end function qhyper_vec

pure function wilcox_counts(m, n) result(counts)
! Runtime helper for R-compatible wilcox counts.
integer, intent(in) :: m ! white ball count
integer, intent(in) :: n ! black ball count or sample size
real(kind=dp), allocatable :: counts(:)
integer :: maxw, offset, maxrank, r, j, s, w
real(kind=dp), allocatable :: dpw(:,:)
maxw = m * n
offset = m * (m + 1) / 2
maxrank = offset + maxw
allocate(dpw(0:m, 0:maxrank))
dpw = 0.0_dp
dpw(0, 0) = 1.0_dp
do r = 1, m + n
   do j = min(r, m), 1, -1
      do s = maxrank, r, -1
         dpw(j, s) = dpw(j, s) + dpw(j - 1, s - r)
      end do
   end do
end do
allocate(counts(maxw + 1))
counts = 0.0_dp
do w = 0, maxw
   counts(w + 1) = dpw(m, w + offset)
end do
end function wilcox_counts

pure function dwilcox_vec(x, m, n, log_) result(out)
! Evaluate distribution helper dwilcox.
real(kind=dp), intent(in) :: x(:) ! input vector
integer, intent(in) :: m ! white ball count
integer, intent(in) :: n ! black ball count or sample size
logical, intent(in), optional :: log_
real(kind=dp), allocatable :: out(:), counts(:)
integer :: i, xx, offset
real(kind=dp) :: den
allocate(out(size(x)))
if (m <= 0 .or. n <= 0) then
   out = ieee_value(0.0_dp, ieee_quiet_nan)
   return
end if
counts = wilcox_counts(m, n)
offset = m * (m + 1) / 2
den = r_choose_real(real(m + n, kind=dp), real(m, kind=dp))
do i = 1, size(x)
   if (.not. ieee_is_finite(x(i))) then
      out(i) = merge(x(i), 0.0_dp, x(i) /= x(i))
      cycle
   end if
   xx = int(x(i))
   if (x(i) == real(xx, kind=dp) .and. xx >= 0 .and. xx <= m * n) then
      out(i) = counts(xx + 1) / den
   else
      out(i) = 0.0_dp
   end if
end do
call maybe_log_density(out, log_)
end function dwilcox_vec

pure function pwilcox_vec(x, m, n) result(out)
! Evaluate distribution helper pwilcox.
real(kind=dp), intent(in) :: x(:) ! input vector
integer, intent(in) :: m ! white ball count
integer, intent(in) :: n ! black ball count or sample size
real(kind=dp), allocatable :: out(:), d(:)
integer :: i, xx
allocate(out(size(x)))
if (m <= 0 .or. n <= 0) then
   out = ieee_value(0.0_dp, ieee_quiet_nan)
   return
end if
d = dwilcox([(real(xx, kind=dp), xx = 0, m * n)], m, n)
do i = 1, size(x)
   if (x(i) /= x(i)) then
      out(i) = x(i)
   else if (x(i) < 0.0_dp) then
      out(i) = 0.0_dp
   else if (x(i) >= real(m * n, kind=dp)) then
      out(i) = 1.0_dp
   else
      xx = int(floor(x(i)))
      out(i) = sum(d(1:xx + 1))
   end if
end do
end function pwilcox_vec

pure function qwilcox_vec(p, m, n) result(out)
! Evaluate distribution helper qwilcox.
real(kind=dp), intent(in) :: p(:) ! probability value
integer, intent(in) :: m ! white ball count
integer, intent(in) :: n ! black ball count or sample size
real(kind=dp), allocatable :: out(:), d(:)
integer :: i, xx
real(kind=dp) :: cdf
allocate(out(size(p)))
if (m <= 0 .or. n <= 0) then
   out = ieee_value(0.0_dp, ieee_quiet_nan)
   return
end if
d = dwilcox([(real(xx, kind=dp), xx = 0, m * n)], m, n)
do i = 1, size(p)
   if (p(i) /= p(i) .or. p(i) < 0.0_dp .or. p(i) > 1.0_dp) then
      out(i) = ieee_value(0.0_dp, ieee_quiet_nan)
   else if (p(i) <= 0.0_dp) then
      out(i) = 0.0_dp
   else if (p(i) >= 1.0_dp) then
      out(i) = real(m * n, kind=dp)
   else
      cdf = 0.0_dp
      do xx = 0, m * n
         cdf = cdf + d(xx + 1)
         if (cdf >= p(i)) exit
      end do
      out(i) = real(xx, kind=dp)
   end if
end do
end function qwilcox_vec

pure function signrank_counts(n) result(counts)
! Runtime helper for R-compatible signrank counts.
integer, intent(in) :: n ! black ball count or sample size
real(kind=dp), allocatable :: counts(:)
integer :: maxs, r, s
maxs = n * (n + 1) / 2
allocate(counts(maxs + 1))
counts = 0.0_dp
counts(1) = 1.0_dp
do r = 1, n
   do s = maxs, r, -1
      counts(s + 1) = counts(s + 1) + counts(s - r + 1)
   end do
end do
end function signrank_counts

pure function dsignrank_vec(x, n, log_) result(out)
! Evaluate distribution helper dsignrank.
real(kind=dp), intent(in) :: x(:) ! input vector
integer, intent(in) :: n ! black ball count or sample size
logical, intent(in), optional :: log_
real(kind=dp), allocatable :: out(:), counts(:)
integer :: i, xx
real(kind=dp) :: den
allocate(out(size(x)))
if (n <= 0) then
   out = ieee_value(0.0_dp, ieee_quiet_nan)
   return
end if
counts = signrank_counts(n)
den = 2.0_dp**n
do i = 1, size(x)
   if (.not. ieee_is_finite(x(i))) then
      out(i) = merge(x(i), 0.0_dp, x(i) /= x(i))
      cycle
   end if
   xx = int(x(i))
   if (x(i) == real(xx, kind=dp) .and. xx >= 0 .and. xx <= n * (n + 1) / 2) then
      out(i) = counts(xx + 1) / den
   else
      out(i) = 0.0_dp
   end if
end do
call maybe_log_density(out, log_)
end function dsignrank_vec

pure function psignrank_vec(x, n) result(out)
! Evaluate distribution helper psignrank.
real(kind=dp), intent(in) :: x(:) ! input vector
integer, intent(in) :: n ! black ball count or sample size
real(kind=dp), allocatable :: out(:), d(:)
integer :: i, xx, maxs
maxs = n * (n + 1) / 2
allocate(out(size(x)))
if (n <= 0) then
   out = ieee_value(0.0_dp, ieee_quiet_nan)
   return
end if
d = dsignrank([(real(xx, kind=dp), xx = 0, maxs)], n)
do i = 1, size(x)
   if (x(i) /= x(i)) then
      out(i) = x(i)
   else if (x(i) < 0.0_dp) then
      out(i) = 0.0_dp
   else if (x(i) >= real(maxs, kind=dp)) then
      out(i) = 1.0_dp
   else
      xx = int(floor(x(i)))
      out(i) = sum(d(1:xx + 1))
   end if
end do
end function psignrank_vec

pure function qsignrank_vec(p, n) result(out)
! Evaluate distribution helper qsignrank.
real(kind=dp), intent(in) :: p(:) ! probability value
integer, intent(in) :: n ! black ball count or sample size
real(kind=dp), allocatable :: out(:), d(:)
integer :: i, xx, maxs
real(kind=dp) :: cdf
maxs = n * (n + 1) / 2
allocate(out(size(p)))
if (n <= 0) then
   out = ieee_value(0.0_dp, ieee_quiet_nan)
   return
end if
d = dsignrank([(real(xx, kind=dp), xx = 0, maxs)], n)
do i = 1, size(p)
   if (p(i) /= p(i) .or. p(i) < 0.0_dp .or. p(i) > 1.0_dp) then
      out(i) = ieee_value(0.0_dp, ieee_quiet_nan)
   else if (p(i) <= 0.0_dp) then
      out(i) = 0.0_dp
   else if (p(i) >= 1.0_dp) then
      out(i) = real(maxs, kind=dp)
   else
      cdf = 0.0_dp
      do xx = 0, maxs
         cdf = cdf + d(xx + 1)
         if (cdf >= p(i)) exit
      end do
      out(i) = real(xx, kind=dp)
   end if
end do
end function qsignrank_vec

pure function chisq_test_real_vec(x, p) result(out)
! Runtime helper for R-compatible chisq test real vec.
real(kind=dp), intent(in) :: x(:) ! input vector
real(kind=dp), intent(in), optional :: p(:) ! probabilities, clamped to [0,1]
type(chisq_test_result_t) :: out
real(kind=dp), allocatable :: expected(:), prob(:)
real(kind=dp) :: total, psum
integer :: n
n = size(x)
if (n <= 1) return
total = sum(x)
if (total <= 0.0_dp) return
if (present(p)) then
   if (size(p) /= n) return
   psum = sum(p)
   if (psum <= 0.0_dp) return
   prob = p / psum
else
   allocate(prob(n))
   prob = 1.0_dp / real(n, kind=dp)
end if
expected = total * prob
if (any(expected <= 0.0_dp)) return
out%statistic = sum((x - expected)**2 / expected)
out%parameter = n - 1
out%p_value = chisq_upper_tail_approx(out%statistic, real(out%parameter, kind=dp))
out%method = 1
end function chisq_test_real_vec

pure function chisq_test_int_vec(x, p) result(out)
! Runtime helper for R-compatible chisq test int vec.
integer, intent(in) :: x(:) ! input vector
real(kind=dp), intent(in), optional :: p(:) ! probabilities, clamped to [0,1]
type(chisq_test_result_t) :: out
if (present(p)) then
   out = chisq_test_real_vec(real(x, kind=dp), p)
else
   out = chisq_test_real_vec(real(x, kind=dp))
end if
end function chisq_test_int_vec

pure function chisq_test_int_mat(x) result(out)
! Runtime helper for R-compatible chisq test int mat.
integer, intent(in) :: x(:,:) ! quantiles or observed values
type(chisq_test_result_t) :: out
out = chisq_test_real_mat(real(x, kind=dp))
end function chisq_test_int_mat

pure function chisq_test_real_mat(x) result(out)
! Runtime helper for R-compatible chisq test real mat.
real(kind=dp), intent(in) :: x(:,:) ! quantiles or observed values
type(chisq_test_result_t) :: out
real(kind=dp), allocatable :: xr(:,:), row_tot(:), col_tot(:), expected(:,:)
real(kind=dp) :: total
integer :: nr, nc, i, j
nr = size(x, 1)
nc = size(x, 2)
if (nr <= 1 .or. nc <= 1) return
xr = x
total = sum(xr)
if (total <= 0.0_dp) return
allocate(row_tot(nr), col_tot(nc), expected(nr, nc))
do i = 1, nr
   row_tot(i) = sum(xr(i, :))
end do
do j = 1, nc
   col_tot(j) = sum(xr(:, j))
end do
do i = 1, nr
   do j = 1, nc
      expected(i, j) = row_tot(i) * col_tot(j) / total
   end do
end do
if (any(expected <= 0.0_dp)) return
out%statistic = sum((xr - expected)**2 / expected)
out%parameter = (nr - 1) * (nc - 1)
out%p_value = chisq_upper_tail_approx(out%statistic, real(out%parameter, kind=dp))
out%method = 2
end function chisq_test_real_mat

pure function prop_test_int_scalar(x, n, p, correct) result(out)
! Evaluate distribution helper prop_test_int_scalar.
integer, intent(in) :: x ! input values
integer, intent(in) :: n ! count or sample size
real(kind=dp), intent(in), optional :: p ! probability value
logical, intent(in), optional :: correct ! continuity-correction flag
type(prop_test_result_t) :: out
out = prop_test_real_scalar(real(x, kind=dp), real(n, kind=dp), p, correct)
end function prop_test_int_scalar

pure function prop_test_real_scalar(x, n, p, correct) result(out)
! Evaluate distribution helper prop_test_real_scalar.
real(kind=dp), intent(in) :: x ! input values
real(kind=dp), intent(in) :: n ! count or sample size
real(kind=dp), intent(in), optional :: p ! probability value
logical, intent(in), optional :: correct ! continuity-correction flag
type(prop_test_result_t) :: out
real(kind=dp) :: p0, diff, denom
logical :: use_correct
p0 = 0.5_dp
if (present(p)) p0 = p
if (n <= 0.0_dp .or. p0 <= 0.0_dp .or. p0 >= 1.0_dp) return
use_correct = .false.
if (present(correct)) use_correct = correct
out%estimate = x / n
out%null_value = p0
out%parameter = 1
diff = abs(x - n * p0)
if (use_correct) diff = max(0.0_dp, diff - 0.5_dp)
denom = n * p0 * (1.0_dp - p0)
if (denom > 0.0_dp) out%statistic = diff * diff / denom
out%p_value = chisq_upper_tail_approx(out%statistic, real(out%parameter, kind=dp))
out%method = 1
end function prop_test_real_scalar

pure function prop_test_int_vec(x, n, p, correct) result(out)
! Evaluate distribution helper prop_test_int_vec.
integer, intent(in) :: x(:) ! input vector
integer, intent(in) :: n(:) ! count or sample size
real(kind=dp), intent(in), optional :: p ! probability value
logical, intent(in), optional :: correct ! continuity-correction flag
type(prop_test_result_t) :: out
out = prop_test_real_vec(real(x, kind=dp), real(n, kind=dp), p, correct)
end function prop_test_int_vec

pure function prop_test_real_vec(x, n, p, correct) result(out)
! Evaluate distribution helper prop_test_real_vec.
real(kind=dp), intent(in) :: x(:) ! input vector
real(kind=dp), intent(in) :: n(:) ! count or sample size
real(kind=dp), intent(in), optional :: p ! probability value
logical, intent(in), optional :: correct ! continuity-correction flag
type(prop_test_result_t) :: out
real(kind=dp) :: p0, pooled, denom, diff
integer :: k, i
logical :: use_correct
k = size(x)
if (k /= size(n) .or. k < 1) return
if (any(n <= 0.0_dp)) return
use_correct = .false.
if (present(correct)) use_correct = correct
if (k == 1) then
   if (present(p)) then
      out = prop_test_real_scalar(x(1), n(1), p, correct)
   else
      out = prop_test_real_scalar(x(1), n(1), correct=correct)
   end if
   return
end if
out%parameter = k - 1
out%estimate = x(1) / n(1)
if (k >= 2) out%estimate2 = x(2) / n(2)
pooled = sum(x) / sum(n)
out%null_value = pooled
if (pooled <= 0.0_dp .or. pooled >= 1.0_dp) return
do i = 1, k
   p0 = x(i) / n(i)
   diff = abs(p0 - pooled)
   if (use_correct .and. k == 2) diff = max(0.0_dp, diff - 0.5_dp / n(i))
   denom = pooled * (1.0_dp - pooled) / n(i)
   if (denom > 0.0_dp) out%statistic = out%statistic + diff * diff / denom
end do
out%p_value = chisq_upper_tail_approx(out%statistic, real(out%parameter, kind=dp))
out%method = 2
end function prop_test_real_vec

pure function cor_test_int_vec(x, y, method) result(out)
! Runtime helper for R-compatible cor test int vec.
! method values starting with "s" or "S" use Spearman ranks; otherwise Pearson.
integer, intent(in) :: x(:) ! first sample values
integer, intent(in) :: y(:) ! second sample values
character(len=*), intent(in), optional :: method ! correlation method selector
type(cor_test_result_t) :: out
if (present(method)) then
   out = cor_test_real_vec(real(x, kind=dp), real(y, kind=dp), method)
else
   out = cor_test_real_vec(real(x, kind=dp), real(y, kind=dp))
end if
end function cor_test_int_vec

pure function cor_test_real_vec(x, y, method) result(out)
! Runtime helper for R-compatible cor test real vec.
! method values starting with "s" or "S" use Spearman ranks; otherwise Pearson.
real(kind=dp), intent(in) :: x(:) ! first sample values
real(kind=dp), intent(in) :: y(:) ! second sample values
character(len=*), intent(in), optional :: method ! correlation method selector
type(cor_test_result_t) :: out
real(kind=dp), allocatable :: rx(:), ry(:)
real(kind=dp) :: r, denom
integer :: n
n = min(size(x), size(y))
if (n <= 2) return
if (present(method)) then
   if (len_trim(method) >= 1 .and. (method(1:1) == "s" .or. method(1:1) == "S")) then
      rx = rank_average(x(1:n))
      ry = rank_average(y(1:n))
      r = cor_vec(rx, ry)
      out%method = 2
   else
      r = cor_vec(x(1:n), y(1:n))
      out%method = 1
   end if
else
   r = cor_vec(x(1:n), y(1:n))
   out%method = 1
end if
out%estimate = r
out%parameter = n - 2
denom = max(0.0_dp, 1.0_dp - r * r)
if (denom > 0.0_dp) out%statistic = r * sqrt(real(out%parameter, kind=dp) / denom)
out%p_value = 2.0_dp * max(0.0_dp, min(1.0_dp, 1.0_dp - normal_cdf(abs(out%statistic))))
end function cor_test_real_vec

pure function fisher_log_choose(n, k) result(out)
! Runtime helper for R-compatible fisher log choose.
integer, intent(in) :: n ! black ball count or sample size
integer, intent(in) :: k ! draw count or group count
real(kind=dp) :: out
if (k < 0 .or. k > n) then
   out = -huge(1.0_dp)
else
   out = log_gamma(real(n + 1, kind=dp)) - log_gamma(real(k + 1, kind=dp)) - &
      & log_gamma(real(n - k + 1, kind=dp))
end if
end function fisher_log_choose

pure function fisher_hyper_prob(a, r1, r2, c1, n) result(out)
! Runtime helper for R-compatible fisher hyper prob.
integer, intent(in) :: a ! input value
integer, intent(in) :: r1 ! input value
integer, intent(in) :: r2 ! input value
integer, intent(in) :: c1 ! input value
integer, intent(in) :: n ! black ball count or sample size
real(kind=dp) :: out, lp
lp = fisher_log_choose(r1, a) + fisher_log_choose(r2, c1 - a) - fisher_log_choose(n, c1)
out = exp(lp)
end function fisher_hyper_prob

pure function fisher_noncentral_mean(log_odds, r1, r2, c1) result(out)
! Runtime helper for R-compatible fisher noncentral mean.
real(kind=dp), intent(in) :: log_odds ! input value
integer, intent(in) :: r1 ! input value
integer, intent(in) :: r2 ! input value
integer, intent(in) :: c1 ! input value
real(kind=dp) :: out, lp, max_lp, w, wsum
integer :: lo, hi, aa
lo = max(0, c1 - r2)
hi = min(r1, c1)
max_lp = -huge(1.0_dp)
do aa = lo, hi
   lp = fisher_log_choose(r1, aa) + fisher_log_choose(r2, c1 - aa) + real(aa, kind=dp) * log_odds
   max_lp = max(max_lp, lp)
end do
out = 0.0_dp
wsum = 0.0_dp
do aa = lo, hi
   lp = fisher_log_choose(r1, aa) + fisher_log_choose(r2, c1 - aa) + real(aa, kind=dp) * log_odds
   w = exp(lp - max_lp)
   out = out + real(aa, kind=dp) * w
   wsum = wsum + w
end do
if (wsum > 0.0_dp) out = out / wsum
end function fisher_noncentral_mean

pure function fisher_odds_mle(a, r1, r2, c1) result(out)
! Runtime helper for R-compatible fisher odds mle.
integer, intent(in) :: a ! input value
integer, intent(in) :: r1 ! input value
integer, intent(in) :: r2 ! input value
integer, intent(in) :: c1 ! input value
real(kind=dp) :: out, lo_eta, hi_eta, mid_eta, target
integer :: lo, hi, iter
lo = max(0, c1 - r2)
hi = min(r1, c1)
if (a <= lo) then
   out = 0.0_dp
   return
end if
if (a >= hi) then
   out = huge(1.0_dp)
   return
end if
target = real(a, kind=dp)
lo_eta = -60.0_dp
hi_eta = 60.0_dp
do iter = 1, 120
   mid_eta = 0.5_dp * (lo_eta + hi_eta)
   if (fisher_noncentral_mean(mid_eta, r1, r2, c1) < target) then
      lo_eta = mid_eta
   else
      hi_eta = mid_eta
   end if
end do
out = exp(0.5_dp * (lo_eta + hi_eta))
end function fisher_odds_mle

pure function fisher_test_int_mat(x) result(out)
! Runtime helper for R-compatible fisher test int mat.
integer, intent(in) :: x(:,:) ! quantiles or observed values
type(fisher_test_result_t) :: out
out = fisher_test_real_mat(real(x, kind=dp))
end function fisher_test_int_mat

pure function fisher_test_real_mat(x) result(out)
! Runtime helper for R-compatible fisher test real mat.
real(kind=dp), intent(in) :: x(:,:) ! quantiles or observed values
type(fisher_test_result_t) :: out
integer :: a, b, c, d, r1, r2, c1, n, lo, hi, aa
real(kind=dp) :: p_obs, p_aa
if (size(x, 1) /= 2 .or. size(x, 2) /= 2) return
if (any(x < 0.0_dp)) return
a = nint(x(1, 1))
b = nint(x(1, 2))
c = nint(x(2, 1))
d = nint(x(2, 2))
if (abs(x(1, 1) - real(a, kind=dp)) > print_int_like_tol) return
if (abs(x(1, 2) - real(b, kind=dp)) > print_int_like_tol) return
if (abs(x(2, 1) - real(c, kind=dp)) > print_int_like_tol) return
if (abs(x(2, 2) - real(d, kind=dp)) > print_int_like_tol) return
r1 = a + b
r2 = c + d
c1 = a + c
n = r1 + r2
if (n <= 0) return
lo = max(0, c1 - r2)
hi = min(r1, c1)
p_obs = fisher_hyper_prob(a, r1, r2, c1, n)
out%p_value = 0.0_dp
do aa = lo, hi
   p_aa = fisher_hyper_prob(aa, r1, r2, c1, n)
   if (p_aa <= p_obs * (1.0_dp + 1000.0_dp * epsilon(1.0_dp))) then
      out%p_value = out%p_value + p_aa
   end if
end do
out%p_value = max(0.0_dp, min(1.0_dp, out%p_value))
out%estimate = fisher_odds_mle(a, r1, r2, c1)
out%method = 1
end function fisher_test_real_mat

pure function wilcox_test_two_sample(x, y, paired) result(out)
! Approximate Wilcoxon rank-sum or paired signed-rank test.
real(kind=dp), intent(in) :: x(:) ! input vector
real(kind=dp), intent(in) :: y(:) ! response values
logical, intent(in), optional :: paired ! paired-sample flag
type(wilcox_test_result_t) :: out
real(kind=dp), allocatable :: xy(:), r(:), d(:), ad(:), rd(:)
real(kind=dp) :: w, mu, sig, z
integer :: nx, ny, n, i
logical :: is_paired
is_paired = .false.
if (present(paired)) is_paired = paired
if (is_paired) then
   n = min(size(x), size(y))
   if (n <= 0) return
   allocate(d(n))
   d = x(1:n) - y(1:n)
   n = count(abs(d) > tiny(1.0_dp))
   if (n <= 0) return
   allocate(ad(n), rd(n))
   ad = pack(abs(d), abs(d) > tiny(1.0_dp))
   rd = rank_average(ad)
   w = 0.0_dp
   i = 0
   do nx = 1, min(size(x), size(y))
      if (abs(x(nx) - y(nx)) > tiny(1.0_dp)) then
         i = i + 1
         if (x(nx) - y(nx) > 0.0_dp) w = w + rd(i)
      end if
   end do
   mu = real(n * (n + 1), kind=dp) / 4.0_dp
   sig = sqrt(real(n * (n + 1) * (2 * n + 1), kind=dp) / 24.0_dp)
   out%statistic = w
   out%method = 2
else
   nx = size(x)
   ny = size(y)
   if (nx <= 0 .or. ny <= 0) return
   allocate(xy(nx + ny))
   xy(1:nx) = x
   xy(nx+1:nx+ny) = y
   r = rank_average(xy)
   w = sum(r(1:nx)) - real(nx * (nx + 1), kind=dp) / 2.0_dp
   mu = real(nx * ny, kind=dp) / 2.0_dp
   sig = sqrt(real(nx * ny * (nx + ny + 1), kind=dp) / 12.0_dp)
   out%statistic = w
   out%method = 1
end if
if (sig > 0.0_dp) then
   z = (abs(w - mu) - 0.5_dp) / sig
   out%p_value = 2.0_dp * max(0.0_dp, min(1.0_dp, 1.0_dp - normal_cdf(abs(z))))
end if
end function wilcox_test_two_sample

pure function kruskal_test(x, g) result(out)
! Approximate Kruskal-Wallis rank-sum test for integer-coded groups.
real(kind=dp), intent(in) :: x(:) ! input vector
integer, intent(in) :: g(:) ! input vector
type(kruskal_test_result_t) :: out
real(kind=dp), allocatable :: r(:), sum_r(:)
integer, allocatable :: n_g(:)
integer :: n, k, i, gi
real(kind=dp) :: h
n = min(size(x), size(g))
if (n <= 1) return
k = maxval(g(1:n))
if (k <= 1) return
allocate(sum_r(k), n_g(k), r(n))
sum_r = 0.0_dp
n_g = 0
r = rank_average(x(1:n))
do i = 1, n
   gi = g(i)
   if (gi >= 1 .and. gi <= k) then
      sum_r(gi) = sum_r(gi) + r(i)
      n_g(gi) = n_g(gi) + 1
   end if
end do
h = 0.0_dp
do gi = 1, k
   if (n_g(gi) > 0) h = h + sum_r(gi) * sum_r(gi) / real(n_g(gi), kind=dp)
end do
h = 12.0_dp * h / real(n * (n + 1), kind=dp) - 3.0_dp * real(n + 1, kind=dp)
out%statistic = max(0.0_dp, h)
out%parameter = k - 1
out%p_value = chisq_upper_tail_approx(out%statistic, real(out%parameter, kind=dp))
end function kruskal_test

pure function t_test_one(x, mu) result(out)
! Runtime helper for R-compatible t test one.
real(kind=dp), intent(in) :: x(:) ! input vector
real(kind=dp), intent(in), optional :: mu ! mean parameter
type(t_test_result_t) :: out
real(kind=dp) :: xbar, s, mu0
integer :: n
n = size(x)
mu0 = 0.0_dp
if (present(mu)) mu0 = mu
if (n <= 1) return
xbar = sum(x) / real(n, kind=dp)
s = sd(x)
out%stderr = s / sqrt(real(n, kind=dp))
out%estimate = xbar
out%null_value = mu0
out%parameter = real(n - 1, kind=dp)
if (out%stderr > 0.0_dp) out%statistic = (xbar - mu0) / out%stderr
out%p_value = 2.0_dp * max(0.0_dp, min(1.0_dp, 1.0_dp - normal_cdf(abs(out%statistic))))
out%conf_low = xbar - 1.96_dp * out%stderr
out%conf_high = xbar + 1.96_dp * out%stderr
out%method = 1
end function t_test_one

pure function t_test_two(x, y, paired, var_equal) result(out)
! Runtime helper for R-compatible t test two.
real(kind=dp), intent(in) :: x(:) ! input vector
real(kind=dp), intent(in) :: y(:) ! response values
logical, intent(in), optional :: paired ! paired-sample flag
logical, intent(in), optional :: var_equal ! equal-variance flag
type(t_test_result_t) :: out
logical :: is_paired, is_equal
real(kind=dp) :: nx, ny, mx, my, vx, vy, sp2, denom
real(kind=dp), allocatable :: d(:)
is_paired = .false.
if (present(paired)) is_paired = paired
is_equal = .false.
if (present(var_equal)) is_equal = var_equal
if (is_paired) then
   if (size(x) /= size(y)) return
   d = x - y
   out = t_test_one(d, 0.0_dp)
   out%method = 4
   return
end if
if (size(x) <= 1 .or. size(y) <= 1) return
nx = real(size(x), kind=dp)
ny = real(size(y), kind=dp)
mx = sum(x) / nx
my = sum(y) / ny
vx = var(x)
vy = var(y)
out%estimate = mx
out%estimate2 = my
out%null_value = 0.0_dp
if (is_equal) then
   sp2 = ((nx - 1.0_dp) * vx + (ny - 1.0_dp) * vy) / (nx + ny - 2.0_dp)
   out%stderr = sqrt(sp2 * (1.0_dp / nx + 1.0_dp / ny))
   out%parameter = nx + ny - 2.0_dp
   out%method = 3
else
   out%stderr = sqrt(vx / nx + vy / ny)
   denom = (vx / nx)**2 / max(1.0_dp, nx - 1.0_dp) + (vy / ny)**2 / max(1.0_dp, ny - 1.0_dp)
   if (denom > 0.0_dp) out%parameter = (vx / nx + vy / ny)**2 / denom
   out%method = 2
end if
if (out%stderr > 0.0_dp) out%statistic = ((mx - my) - out%null_value) / out%stderr
out%p_value = 2.0_dp * max(0.0_dp, min(1.0_dp, 1.0_dp - normal_cdf(abs(out%statistic))))
out%conf_low = (mx - my) - 1.96_dp * out%stderr
out%conf_high = (mx - my) + 1.96_dp * out%stderr
end function t_test_two

pure function t_test_p_value_one(x, mu) result(p)
! Runtime helper for R-compatible t test p value one.
real(kind=dp), intent(in) :: x(:) ! input vector
real(kind=dp), intent(in), optional :: mu ! mean parameter
real(kind=dp) :: p
type(t_test_result_t) :: fit
fit = t_test_one(x, mu)
p = fit%p_value
end function t_test_p_value_one

pure function t_test_p_value_two(x, y, paired, var_equal) result(p)
! Runtime helper for R-compatible t test p value two.
real(kind=dp), intent(in) :: x(:) ! input vector
real(kind=dp), intent(in) :: y(:) ! response values
logical, intent(in), optional :: paired ! paired-sample flag
logical, intent(in), optional :: var_equal ! equal-variance flag
real(kind=dp) :: p
type(t_test_result_t) :: fit
fit = t_test_two(x, y, paired, var_equal)
p = fit%p_value
end function t_test_p_value_two

subroutine print_t_test(fit)
! Print t test values in an R-like format.
type(t_test_result_t), intent(in) :: fit ! fitted model or test result
select case (fit%method)
case (1)
   write(*,'(a)') "One Sample t-test"
case (2)
   write(*,'(a)') "Welch Two Sample t-test"
case (3)
   write(*,'(a)') "Two Sample t-test"
case (4)
   write(*,'(a)') "Paired t-test"
case default
   write(*,'(a)') "t-test"
end select
write(*,'(a,g0,a,g0,a,g0)') "t = ", fit%statistic, ", df = ", fit%parameter, ", p-value = ", fit%p_value
write(*,'(a,g0,a,g0)') "95 percent confidence interval: ", fit%conf_low, " ", fit%conf_high
if (fit%method == 2 .or. fit%method == 3) then
   write(*,'(a,g0,a,g0)') "sample estimates: ", fit%estimate, " ", fit%estimate2
else
   write(*,'(a,g0)') "sample estimate: ", fit%estimate
end if
end subroutine print_t_test

subroutine print_chisq_test(fit)
! Print chisq test values in an R-like format.
type(chisq_test_result_t), intent(in) :: fit ! fitted model or test result
select case (fit%method)
case (1)
   write(*,'(a)') "Chi-squared goodness-of-fit test"
case (2)
   write(*,'(a)') "Pearson's Chi-squared test"
case default
   write(*,'(a)') "Chi-squared test"
end select
write(*,'(a,g0,a,i0,a,g0)') "X-squared = ", fit%statistic, ", df = ", fit%parameter, &
   & ", p-value = ", fit%p_value
end subroutine print_chisq_test

subroutine print_prop_test(fit)
! Print prop test values in an R-like format.
type(prop_test_result_t), intent(in) :: fit ! fitted model or test result
select case (fit%method)
case (1)
   write(*,'(a)') "1-sample proportions test"
case (2)
   write(*,'(a)') "2-sample test for equality of proportions"
case default
   write(*,'(a)') "proportions test"
end select
write(*,'(a,g0,a,i0,a,g0)') "X-squared = ", fit%statistic, ", df = ", fit%parameter, &
   & ", p-value = ", fit%p_value
if (fit%method == 2) then
   write(*,'(a,g0,a,g0)') "sample estimates: ", fit%estimate, " ", fit%estimate2
else
   write(*,'(a,g0)') "sample estimate: ", fit%estimate
   write(*,'(a,g0)') "null value: ", fit%null_value
end if
end subroutine print_prop_test

subroutine print_cor_test(fit)
! Print cor test values in an R-like format.
type(cor_test_result_t), intent(in) :: fit ! fitted model or test result
select case (fit%method)
case (2)
   write(*,'(a)') "Spearman's rank correlation rho"
   write(*,'(a,g0,a,i0,a,g0)') "t = ", fit%statistic, ", df = ", fit%parameter, &
      & ", p-value = ", fit%p_value
   write(*,'(a,g0)') "sample estimate rho: ", fit%estimate
case default
   write(*,'(a)') "Pearson's product-moment correlation"
   write(*,'(a,g0,a,i0,a,g0)') "t = ", fit%statistic, ", df = ", fit%parameter, &
      & ", p-value = ", fit%p_value
   write(*,'(a,g0)') "sample estimate cor: ", fit%estimate
end select
end subroutine print_cor_test

subroutine print_fisher_test(fit)
! Print fisher test values in an R-like format.
type(fisher_test_result_t), intent(in) :: fit ! fitted model or test result
write(*,'(a)') "Fisher's Exact Test for Count Data"
write(*,'(a,g0)') "p-value = ", fit%p_value
write(*,'(a,g0)') "sample estimate odds ratio: ", fit%estimate
end subroutine print_fisher_test

subroutine print_wilcox_test(fit)
! Print wilcox test values in an R-like format.
type(wilcox_test_result_t), intent(in) :: fit ! fitted model or test result
select case (fit%method)
case (2)
   write(*,'(a)') "Wilcoxon signed rank test with continuity correction"
   if (abs(fit%statistic - real(nint(fit%statistic), kind=dp)) <= print_int_like_tol) then
      write(*,'(a,i0,a,g0)') "V = ", nint(fit%statistic), ", p-value = ", fit%p_value
   else
      write(*,'(a,g0,a,g0)') "V = ", fit%statistic, ", p-value = ", fit%p_value
   end if
case default
   write(*,'(a)') "Wilcoxon rank sum test with continuity correction"
   if (abs(fit%statistic - real(nint(fit%statistic), kind=dp)) <= print_int_like_tol) then
      write(*,'(a,i0,a,g0)') "W = ", nint(fit%statistic), ", p-value = ", fit%p_value
   else
      write(*,'(a,g0,a,g0)') "W = ", fit%statistic, ", p-value = ", fit%p_value
   end if
end select
end subroutine print_wilcox_test

subroutine print_kruskal_test(fit)
! Print kruskal test values in an R-like format.
type(kruskal_test_result_t), intent(in) :: fit ! fitted model or test result
write(*,'(a)') "Kruskal-Wallis rank sum test"
write(*,'(a,g0,a,i0,a,g0)') "Kruskal-Wallis chi-squared = ", fit%statistic, &
   & ", df = ", fit%parameter, ", p-value = ", fit%p_value
end subroutine print_kruskal_test

pure function ks_test(x, mean, sd) result(out)
! Runtime helper for R-compatible ks test.
real(kind=dp), intent(in) :: x(:) ! input vector
real(kind=dp), intent(in), optional :: mean ! input value
real(kind=dp), intent(in), optional :: sd ! input value
type(ks_test_result_t) :: out
real(kind=dp), allocatable :: xs(:)
real(kind=dp) :: mu, sig, fi, dplus, dminus, z
integer :: i, n
mu = 0.0_dp
sig = 1.0_dp
if (present(mean)) mu = mean
if (present(sd)) sig = sd
n = size(x)
out%n = n
if (n <= 0 .or. sig <= 0.0_dp) return
xs = sort(x)
do i = 1, n
   fi = normal_cdf((xs(i) - mu) / sig)
   dplus = real(i, kind=dp) / real(n, kind=dp) - fi
   dminus = fi - real(i - 1, kind=dp) / real(n, kind=dp)
   out%statistic = max(out%statistic, max(dplus, dminus))
end do
z = (sqrt(real(n, kind=dp)) + 0.12_dp + 0.11_dp / sqrt(real(n, kind=dp))) * out%statistic
out%p_value = max(0.0_dp, min(1.0_dp, 2.0_dp * exp(-2.0_dp * z * z)))
end function ks_test

subroutine print_ks_test(fit)
! Print ks test values in an R-like format.
type(ks_test_result_t), intent(in) :: fit ! fitted model or test result
write(*,'(a)') "One-sample Kolmogorov-Smirnov test"
write(*,'(a,g0,a,g0)') "D = ", fit%statistic, ", p-value = ", fit%p_value
end subroutine print_ks_test

subroutine print_factanal(x, factors)
! Print factanal values in an R-like format.
real(kind=dp), intent(in) :: x(:,:) ! input matrix
integer, intent(in) :: factors ! input value
type(eigen_result_t) :: eg
real(kind=dp), allocatable :: loadings(:,:)
integer :: j, nf
nf = max(1, min(factors, size(x, 2)))
eg = eigen(cor(x))
allocate(loadings(size(x, 2), nf))
do j = 1, nf
   loadings(:, j) = real(eg%vectors(:, j), kind=dp) * sqrt(max(0.0_dp, real(eg%values(j), kind=dp)))
end do
write(*,'(a)') "Factor Analysis (principal-factor approximation)"
write(*,'(a)') "Loadings:"
call print_matrix_rstyle(loadings)
write(*,'(a)') "Uniquenesses:"
call print_real_vector(max(0.0_dp, 1.0_dp - sum(loadings**2, dim=2)))
end subroutine print_factanal

pure function decompose(x, type, frequency) result(out)
! Evaluate distribution helper decompose.
! type is accepted for API compatibility and ignored; frequency defaults to 1.
real(kind=dp), intent(in) :: x(:) ! time-series values
character(len=*), intent(in), optional :: type ! accepted for compatibility; ignored
integer, intent(in), optional :: frequency ! seasonal frequency
type(decompose_result_t) :: out
integer :: n, f, i, j, lo, hi, cnt
real(kind=dp), allocatable :: detr(:)
if (present(type)) continue
n = size(x)
f = 1
if (present(frequency)) f = max(1, frequency)
allocate(out%trend(n), out%seasonal(n), out%random(n), out%figure(f), detr(n))
out%trend = 0.0_dp
do i = 1, n
   lo = max(1, i - f/2)
   hi = min(n, i + f/2)
   out%trend(i) = sum(x(lo:hi)) / real(hi - lo + 1, kind=dp)
end do
detr = x - out%trend
out%figure = 0.0_dp
do j = 1, f
   cnt = 0
   do i = j, n, f
      out%figure(j) = out%figure(j) + detr(i)
      cnt = cnt + 1
   end do
   if (cnt > 0) out%figure(j) = out%figure(j) / real(cnt, kind=dp)
end do
out%figure = out%figure - sum(out%figure) / real(f, kind=dp)
do i = 1, n
   j = mod(i - 1, f) + 1
   out%seasonal(i) = out%figure(j)
end do
out%random = x - out%trend - out%seasonal
end function decompose

pure function ecdf_eval(x, q) result(out)
! Runtime helper for R-compatible ecdf eval.
real(kind=dp), intent(in) :: x(:) ! quantiles or observed values
real(kind=dp), intent(in) :: q(:) ! quantiles
real(kind=dp), allocatable :: out(:)
integer :: i, j, n_le, n_valid
allocate(out(size(q)))
n_valid = 0
do j = 1, size(x)
   if (x(j) == x(j)) n_valid = n_valid + 1
end do
do i = 1, size(q)
   if (r_is_na_payload(q(i))) then
      out(i) = r_na_real()
   else if (q(i) /= q(i)) then
      out(i) = q(i)
   else if (n_valid <= 0) then
      out(i) = ieee_value(0.0_dp, ieee_quiet_nan)
   else
      n_le = 0
      do j = 1, size(x)
         if (x(j) == x(j)) then
            if (x(j) <= q(i)) n_le = n_le + 1
         end if
      end do
      out(i) = real(n_le, kind=dp) / real(n_valid, kind=dp)
   end if
end do
end function ecdf_eval

pure function qnorm_scalar(pv, mean, sd, lower_tail) result(out)
real(kind=dp), intent(in) :: pv ! probabilities, clamped to [0,1]
real(kind=dp), intent(in), optional :: mean ! distribution mean
real(kind=dp), intent(in), optional :: sd ! distribution standard deviation
logical, intent(in), optional :: lower_tail ! lower-tail probability flag
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = qnorm_vec([pv], mean=mean, sd=sd, lower_tail=lower_tail)
out = tmp(1)
end function qnorm_scalar

pure function ppois_scalar(q, lambda, lower_tail) result(out)
real(kind=dp), intent(in) :: q ! quantiles
real(kind=dp), intent(in), optional :: lambda ! rate parameter
logical, intent(in), optional :: lower_tail ! lower-tail probability flag
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = ppois_vec([q], lambda=lambda, lower_tail=lower_tail)
out = tmp(1)
end function ppois_scalar

pure function qpois_scalar(pv, lambda, lower_tail) result(out)
real(kind=dp), intent(in) :: pv ! probabilities, clamped to [0,1]
real(kind=dp), intent(in), optional :: lambda ! rate parameter
logical, intent(in), optional :: lower_tail ! lower-tail probability flag
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = qpois_vec([pv], lambda=lambda, lower_tail=lower_tail)
out = tmp(1)
end function qpois_scalar

pure function dunif_scalar(x, min, max, log_) result(out)
real(kind=dp), intent(in) :: x ! input vector
real(kind=dp), intent(in), optional :: min, max
logical, intent(in), optional :: log_
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = dunif_vec([x], min=min, max=max, log_=log_)
out = tmp(1)
end function dunif_scalar

pure function punif_scalar(x, min, max) result(out)
real(kind=dp), intent(in) :: x ! input vector
real(kind=dp), intent(in), optional :: min, max
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = punif_vec([x], min=min, max=max)
out = tmp(1)
end function punif_scalar

pure function qunif_scalar(p, min, max) result(out)
real(kind=dp), intent(in) :: p ! probability value
real(kind=dp), intent(in), optional :: min, max
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = qunif_vec([p], min=min, max=max)
out = tmp(1)
end function qunif_scalar

pure function dexp_scalar(x, rate, log_) result(out)
real(kind=dp), intent(in) :: x ! input vector
real(kind=dp), intent(in), optional :: rate ! rate parameter
logical, intent(in), optional :: log_
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = dexp_vec([x], rate=rate, log_=log_)
out = tmp(1)
end function dexp_scalar

pure function pexp_scalar(x, rate) result(out)
real(kind=dp), intent(in) :: x ! input vector
real(kind=dp), intent(in), optional :: rate
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = pexp_vec([x], rate=rate)
out = tmp(1)
end function pexp_scalar

pure function qexp_scalar(p, rate) result(out)
real(kind=dp), intent(in) :: p ! probability value
real(kind=dp), intent(in), optional :: rate
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = qexp_vec([p], rate=rate)
out = tmp(1)
end function qexp_scalar

pure function dgamma_scalar(x, shape, rate, log_) result(out)
real(kind=dp), intent(in) :: x ! input vector
real(kind=dp), intent(in) :: shape ! shape parameter
real(kind=dp), intent(in), optional :: rate
logical, intent(in), optional :: log_
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = dgamma_vec([x], shape=shape, rate=rate, log_=log_)
out = tmp(1)
end function dgamma_scalar

pure function pgamma_scalar(x, shape, rate) result(out)
real(kind=dp), intent(in) :: x ! input vector
real(kind=dp), intent(in) :: shape ! shape parameter
real(kind=dp), intent(in), optional :: rate
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = pgamma_vec([x], shape=shape, rate=rate)
out = tmp(1)
end function pgamma_scalar

pure function qgamma_scalar(p, shape, rate) result(out)
real(kind=dp), intent(in) :: p ! probability value
real(kind=dp), intent(in) :: shape ! shape parameter
real(kind=dp), intent(in), optional :: rate
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = qgamma_vec([p], shape=shape, rate=rate)
out = tmp(1)
end function qgamma_scalar

pure function dbeta_scalar(x, shape1, shape2, log_) result(out)
real(kind=dp), intent(in) :: x ! quantiles or observed values
real(kind=dp), intent(in) :: shape1 ! first shape parameter
real(kind=dp), intent(in) :: shape2 ! second shape parameter
logical, intent(in), optional :: log_
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = dbeta_vec([x], shape1=shape1, shape2=shape2, log_=log_)
out = tmp(1)
end function dbeta_scalar

pure function pbeta_scalar(x, shape1, shape2) result(out)
real(kind=dp), intent(in) :: x ! quantiles or observed values
real(kind=dp), intent(in) :: shape1 ! first shape parameter
real(kind=dp), intent(in) :: shape2 ! second shape parameter
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = pbeta_vec([x], shape1=shape1, shape2=shape2)
out = tmp(1)
end function pbeta_scalar

pure function qbeta_scalar(p, shape1, shape2) result(out)
real(kind=dp), intent(in) :: p ! probabilities, clamped to [0,1]
real(kind=dp), intent(in) :: shape1 ! first shape parameter
real(kind=dp), intent(in) :: shape2 ! second shape parameter
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = qbeta_vec([p], shape1=shape1, shape2=shape2)
out = tmp(1)
end function qbeta_scalar

pure function dchisq_scalar(x, df, log_) result(out)
real(kind=dp), intent(in) :: x ! quantiles or observed values
real(kind=dp), intent(in) :: df ! degrees of freedom
logical, intent(in), optional :: log_
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = dchisq_vec([x], df=df, log_=log_)
out = tmp(1)
end function dchisq_scalar

pure function dt_scalar(x, df, log_) result(out)
real(kind=dp), intent(in) :: x ! quantiles or observed values
real(kind=dp), intent(in) :: df ! degrees of freedom
logical, intent(in), optional :: log_
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = dt_vec([x], df=df, log_=log_)
out = tmp(1)
end function dt_scalar

pure function pt_scalar(x, df) result(out)
real(kind=dp), intent(in) :: x ! quantiles or observed values
real(kind=dp), intent(in) :: df ! degrees of freedom
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = pt_vec([x], df=df)
out = tmp(1)
end function pt_scalar

pure function qt_scalar(p, df) result(out)
real(kind=dp), intent(in) :: p ! probabilities, clamped to [0,1]
real(kind=dp), intent(in) :: df ! degrees of freedom
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = qt_vec([p], df=df)
out = tmp(1)
end function qt_scalar

pure function df_scalar(x, df1, df2, log_) result(out)
real(kind=dp), intent(in) :: x ! quantiles or observed values
real(kind=dp), intent(in) :: df1 ! numerator degrees of freedom
real(kind=dp), intent(in) :: df2 ! denominator degrees of freedom
logical, intent(in), optional :: log_
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = df_vec([x], df1=df1, df2=df2, log_=log_)
out = tmp(1)
end function df_scalar

pure function pf_scalar(x, df1, df2) result(out)
real(kind=dp), intent(in) :: x ! quantiles or observed values
real(kind=dp), intent(in) :: df1 ! numerator degrees of freedom
real(kind=dp), intent(in) :: df2 ! denominator degrees of freedom
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = pf_vec([x], df1=df1, df2=df2)
out = tmp(1)
end function pf_scalar

pure function qf_scalar(p, df1, df2) result(out)
real(kind=dp), intent(in) :: p ! probabilities, clamped to [0,1]
real(kind=dp), intent(in) :: df1 ! numerator degrees of freedom
real(kind=dp), intent(in) :: df2 ! denominator degrees of freedom
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = qf_vec([p], df1=df1, df2=df2)
out = tmp(1)
end function qf_scalar

pure function dlogis_scalar(x, location, scale, log_) result(out)
real(kind=dp), intent(in) :: x ! input vector
real(kind=dp), intent(in), optional :: location ! location parameter
real(kind=dp), intent(in), optional :: scale ! scale parameter
logical, intent(in), optional :: log_
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = dlogis_vec([x], location=location, scale=scale, log_=log_)
out = tmp(1)
end function dlogis_scalar

pure function plogis_scalar(x, location, scale) result(out)
real(kind=dp), intent(in) :: x ! input vector
real(kind=dp), intent(in), optional :: location ! location parameter
real(kind=dp), intent(in), optional :: scale ! scale parameter
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = plogis_vec([x], location=location, scale=scale)
out = tmp(1)
end function plogis_scalar

pure function qlogis_scalar(p, location, scale) result(out)
real(kind=dp), intent(in) :: p ! probability value
real(kind=dp), intent(in), optional :: location ! location parameter
real(kind=dp), intent(in), optional :: scale ! scale parameter
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = qlogis_vec([p], location=location, scale=scale)
out = tmp(1)
end function qlogis_scalar

pure function dlnorm_scalar(x, meanlog, sdlog, log_) result(out)
real(kind=dp), intent(in) :: x ! input vector
real(kind=dp), intent(in), optional :: meanlog ! mean on log scale
real(kind=dp), intent(in), optional :: sdlog ! standard deviation on log scale
logical, intent(in), optional :: log_
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = dlnorm_vec([x], meanlog=meanlog, sdlog=sdlog, log_=log_)
out = tmp(1)
end function dlnorm_scalar

pure function plnorm_scalar(x, meanlog, sdlog) result(out)
real(kind=dp), intent(in) :: x ! input vector
real(kind=dp), intent(in), optional :: meanlog ! mean on log scale
real(kind=dp), intent(in), optional :: sdlog ! standard deviation on log scale
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = plnorm_vec([x], meanlog=meanlog, sdlog=sdlog)
out = tmp(1)
end function plnorm_scalar

pure function qlnorm_scalar(p, meanlog, sdlog) result(out)
real(kind=dp), intent(in) :: p ! probability value
real(kind=dp), intent(in), optional :: meanlog ! mean on log scale
real(kind=dp), intent(in), optional :: sdlog ! standard deviation on log scale
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = qlnorm_vec([p], meanlog=meanlog, sdlog=sdlog)
out = tmp(1)
end function qlnorm_scalar

pure function dweibull_scalar(x, shape, scale, log_) result(out)
real(kind=dp), intent(in) :: x ! input vector
real(kind=dp), intent(in) :: shape ! shape parameter
real(kind=dp), intent(in), optional :: scale ! scale parameter
logical, intent(in), optional :: log_
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = dweibull_vec([x], shape=shape, scale=scale, log_=log_)
out = tmp(1)
end function dweibull_scalar

pure function pweibull_scalar(x, shape, scale) result(out)
real(kind=dp), intent(in) :: x ! input vector
real(kind=dp), intent(in) :: shape ! shape parameter
real(kind=dp), intent(in), optional :: scale ! scale parameter
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = pweibull_vec([x], shape=shape, scale=scale)
out = tmp(1)
end function pweibull_scalar

pure function qweibull_scalar(p, shape, scale) result(out)
real(kind=dp), intent(in) :: p ! probability value
real(kind=dp), intent(in) :: shape ! shape parameter
real(kind=dp), intent(in), optional :: scale ! scale parameter
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = qweibull_vec([p], shape=shape, scale=scale)
out = tmp(1)
end function qweibull_scalar

pure function dcauchy_scalar(x, location, scale, log_) result(out)
real(kind=dp), intent(in) :: x ! input vector
real(kind=dp), intent(in), optional :: location ! location parameter
real(kind=dp), intent(in), optional :: scale ! scale parameter
logical, intent(in), optional :: log_
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = dcauchy_vec([x], location=location, scale=scale, log_=log_)
out = tmp(1)
end function dcauchy_scalar

pure function pcauchy_scalar(x, location, scale) result(out)
real(kind=dp), intent(in) :: x ! input vector
real(kind=dp), intent(in), optional :: location ! location parameter
real(kind=dp), intent(in), optional :: scale ! scale parameter
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = pcauchy_vec([x], location=location, scale=scale)
out = tmp(1)
end function pcauchy_scalar

pure function qcauchy_scalar(p, location, scale) result(out)
real(kind=dp), intent(in) :: p ! probability value
real(kind=dp), intent(in), optional :: location ! location parameter
real(kind=dp), intent(in), optional :: scale ! scale parameter
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = qcauchy_vec([p], location=location, scale=scale)
out = tmp(1)
end function qcauchy_scalar

pure function dbinom_scalar(x, nsize, prob, log_) result(out)
real(kind=dp), intent(in) :: x ! input vector
integer, intent(in) :: nsize ! input value
real(kind=dp), intent(in) :: prob ! success probability
logical, intent(in), optional :: log_
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = dbinom_vec([x], nsize=nsize, prob=prob, log_=log_)
out = tmp(1)
end function dbinom_scalar

pure function pbinom_scalar(x, nsize, prob) result(out)
real(kind=dp), intent(in) :: x ! input vector
integer, intent(in) :: nsize ! input value
real(kind=dp), intent(in) :: prob ! success probability
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = pbinom_vec([x], nsize=nsize, prob=prob)
out = tmp(1)
end function pbinom_scalar

pure function qbinom_scalar(p, nsize, prob) result(out)
real(kind=dp), intent(in) :: p ! probability value
integer, intent(in) :: nsize ! input value
real(kind=dp), intent(in) :: prob ! success probability
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = qbinom_vec([p], nsize=nsize, prob=prob)
out = tmp(1)
end function qbinom_scalar

pure function dpois_scalar(x, lambda, log_) result(out)
real(kind=dp), intent(in) :: x ! input vector
real(kind=dp), intent(in), optional :: lambda ! rate parameter
logical, intent(in), optional :: log_
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = dpois_vec([x], lambda=lambda, log_=log_)
out = tmp(1)
end function dpois_scalar

pure function dgeom_scalar(x, prob, log_) result(out)
real(kind=dp), intent(in) :: x ! quantiles or observed values
real(kind=dp), intent(in) :: prob ! success probability
logical, intent(in), optional :: log_
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = dgeom_vec([x], prob=prob, log_=log_)
out = tmp(1)
end function dgeom_scalar

pure function pgeom_scalar(x, prob) result(out)
real(kind=dp), intent(in) :: x ! quantiles or observed values
real(kind=dp), intent(in) :: prob ! success probability
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = pgeom_vec([x], prob=prob)
out = tmp(1)
end function pgeom_scalar

pure function qgeom_scalar(p, prob) result(out)
real(kind=dp), intent(in) :: p ! probabilities, clamped to [0,1]
real(kind=dp), intent(in) :: prob ! success probability
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = qgeom_vec([p], prob=prob)
out = tmp(1)
end function qgeom_scalar

pure function dnbinom_scalar(x, nsize, prob, log_) result(out)
real(kind=dp), intent(in) :: x ! quantiles or observed values
real(kind=dp), intent(in) :: prob ! success probability
integer, intent(in) :: nsize ! input value
logical, intent(in), optional :: log_
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = dnbinom_vec([x], nsize=nsize, prob=prob, log_=log_)
out = tmp(1)
end function dnbinom_scalar

pure function pnbinom_scalar(x, nsize, prob) result(out)
real(kind=dp), intent(in) :: x ! quantiles or observed values
real(kind=dp), intent(in) :: prob ! success probability
integer, intent(in) :: nsize ! input value
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = pnbinom_vec([x], nsize=nsize, prob=prob)
out = tmp(1)
end function pnbinom_scalar

pure function qnbinom_scalar(p, nsize, prob) result(out)
real(kind=dp), intent(in) :: p ! probabilities, clamped to [0,1]
real(kind=dp), intent(in) :: prob ! success probability
integer, intent(in) :: nsize ! input value
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = qnbinom_vec([p], nsize=nsize, prob=prob)
out = tmp(1)
end function qnbinom_scalar

pure function dhyper_scalar(x, m, n, k, log_) result(out)
real(kind=dp), intent(in) :: x ! input vector
integer, intent(in) :: m ! white ball count
integer, intent(in) :: n ! black ball count or sample size
integer, intent(in) :: k ! draw count or group count
logical, intent(in), optional :: log_
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = dhyper_vec([x], m=m, n=n, k=k, log_=log_)
out = tmp(1)
end function dhyper_scalar

pure function phyper_scalar(x, m, n, k) result(out)
real(kind=dp), intent(in) :: x ! input vector
integer, intent(in) :: m ! white ball count
integer, intent(in) :: n ! black ball count or sample size
integer, intent(in) :: k ! draw count or group count
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = phyper_vec([x], m=m, n=n, k=k)
out = tmp(1)
end function phyper_scalar

pure function qhyper_scalar(p, m, n, k) result(out)
real(kind=dp), intent(in) :: p ! probability value
integer, intent(in) :: m ! white ball count
integer, intent(in) :: n ! black ball count or sample size
integer, intent(in) :: k ! draw count or group count
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = qhyper_vec([p], m=m, n=n, k=k)
out = tmp(1)
end function qhyper_scalar

pure function dwilcox_scalar(x, m, n, log_) result(out)
real(kind=dp), intent(in) :: x ! input vector
integer, intent(in) :: m ! white ball count
integer, intent(in) :: n ! black ball count or sample size
logical, intent(in), optional :: log_
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = dwilcox_vec([x], m=m, n=n, log_=log_)
out = tmp(1)
end function dwilcox_scalar

pure function pwilcox_scalar(x, m, n) result(out)
real(kind=dp), intent(in) :: x ! input vector
integer, intent(in) :: m ! white ball count
integer, intent(in) :: n ! black ball count or sample size
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = pwilcox_vec([x], m=m, n=n)
out = tmp(1)
end function pwilcox_scalar

pure function qwilcox_scalar(p, m, n) result(out)
real(kind=dp), intent(in) :: p ! probability value
integer, intent(in) :: m ! white ball count
integer, intent(in) :: n ! black ball count or sample size
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = qwilcox_vec([p], m=m, n=n)
out = tmp(1)
end function qwilcox_scalar

pure function dsignrank_scalar(x, n, log_) result(out)
real(kind=dp), intent(in) :: x ! input vector
integer, intent(in) :: n ! black ball count or sample size
logical, intent(in), optional :: log_
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = dsignrank_vec([x], n=n, log_=log_)
out = tmp(1)
end function dsignrank_scalar

pure function psignrank_scalar(x, n) result(out)
real(kind=dp), intent(in) :: x ! input vector
integer, intent(in) :: n ! black ball count or sample size
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = psignrank_vec([x], n=n)
out = tmp(1)
end function psignrank_scalar

pure function qsignrank_scalar(p, n) result(out)
real(kind=dp), intent(in) :: p ! probability value
integer, intent(in) :: n ! black ball count or sample size
real(kind=dp) :: out
real(kind=dp), allocatable :: tmp(:)
tmp = qsignrank_vec([p], n=n)
out = tmp(1)
end function qsignrank_scalar

pure function data_frame_real(names, cols) result(df)
character(len=*), intent(in) :: names(:)
real(kind=dp), intent(in) :: cols(:,:)
type(r_dataframe_t) :: df
integer :: j, name_len
name_len = 1
do j = 1, size(names)
   name_len = max(name_len, len_trim(names(j)))
end do
allocate(character(len=name_len) :: df%names(size(names)))
df%names = names
df%real_cols = cols
end function data_frame_real

function tibble_real(names, cols, row_labels, row_label_name) result(tbl)
character(len=*), intent(in) :: names(:)
real(kind=dp), intent(in) :: cols(:,:)
character(len=*), intent(in), optional :: row_labels(:)
character(len=*), intent(in), optional :: row_label_name
type(r_tibble_real_t) :: tbl
integer :: j, k, name_len, row_label_len

if (size(names) /= size(cols, 2)) &
   error stop "tibble_real: number of names must match number of columns"
do j = 1, size(names)
   if (len_trim(names(j)) == 0) error stop "tibble_real: column names must not be empty"
   do k = 1, j - 1
      if (trim(names(j)) == trim(names(k))) error stop "tibble_real: duplicate column name"
   end do
end do

name_len = 1
do j = 1, size(names)
   name_len = max(name_len, len_trim(names(j)))
end do
allocate(character(len=name_len) :: tbl%names(size(names)))
tbl%names = names
tbl%real_cols = cols
if (present(row_labels)) then
   if (size(row_labels) /= size(cols, 1)) &
      error stop "tibble_real: number of row labels must match number of rows"
   row_label_len = 1
   do j = 1, size(row_labels)
      row_label_len = max(row_label_len, len_trim(row_labels(j)))
   end do
   allocate(character(len=row_label_len) :: tbl%row_labels(size(row_labels)))
   tbl%row_labels = row_labels
   if (present(row_label_name)) then
      if (len_trim(row_label_name) == 0) &
         error stop "tibble_real: row label name must not be empty"
      tbl%row_label_name = trim(row_label_name)
   end if
else if (present(row_label_name)) then
   error stop "tibble_real: row label name requires row labels"
end if
end function tibble_real

pure integer function tibble_real_nrow(tbl) result(n)
type(r_tibble_real_t), intent(in) :: tbl

n = 0
if (allocated(tbl%real_cols)) n = size(tbl%real_cols, 1)
end function tibble_real_nrow

pure integer function tibble_real_ncol(tbl) result(n)
type(r_tibble_real_t), intent(in) :: tbl

n = 0
if (allocated(tbl%real_cols)) n = size(tbl%real_cols, 2)
end function tibble_real_ncol

function tibble_real_col(tbl, name) result(col)
type(r_tibble_real_t), intent(in) :: tbl
character(len=*), intent(in) :: name
real(kind=dp), allocatable :: col(:)
integer :: j

do j = 1, tibble_ncol(tbl)
   if (trim(tbl%names(j)) == trim(name)) then
      col = tbl%real_cols(:, j)
      return
   end if
end do
error stop "tibble_real_col: column not found"
end function tibble_real_col

function tibble_real_filter(tbl, keep) result(out)
type(r_tibble_real_t), intent(in) :: tbl
logical, intent(in) :: keep(:)
type(r_tibble_real_t) :: out
real(kind=dp), allocatable :: cols(:,:)
integer :: j

if (size(keep) /= tibble_nrow(tbl)) &
   error stop "tibble_real_filter: mask length must match number of rows"
allocate(cols(count(keep), tibble_ncol(tbl)))
do j = 1, tibble_ncol(tbl)
   cols(:, j) = pack(tbl%real_cols(:, j), keep)
end do
if (allocated(tbl%row_labels)) then
   out = tibble_real(tbl%names, cols, pack(tbl%row_labels, keep))
   if (allocated(tbl%row_label_name)) out%row_label_name = tbl%row_label_name
else
   out = tibble_real(tbl%names, cols)
end if
end function tibble_real_filter

function tibble_real_select(tbl, selected_names) result(out)
type(r_tibble_real_t), intent(in) :: tbl
character(len=*), intent(in) :: selected_names(:)
type(r_tibble_real_t) :: out
real(kind=dp), allocatable :: cols(:,:)
integer :: i, j
logical :: found

allocate(cols(tibble_nrow(tbl), size(selected_names)))
do i = 1, size(selected_names)
   found = .false.
   do j = 1, tibble_ncol(tbl)
      if (trim(tbl%names(j)) == trim(selected_names(i))) then
         cols(:, i) = tbl%real_cols(:, j)
         found = .true.
         exit
      end if
   end do
   if (.not. found) error stop "tibble_real_select: column not found"
end do
if (allocated(tbl%row_labels)) then
   out = tibble_real(selected_names, cols, tbl%row_labels)
   if (allocated(tbl%row_label_name)) out%row_label_name = tbl%row_label_name
else
   out = tibble_real(selected_names, cols)
end if
end function tibble_real_select

function tibble_real_drop(tbl, dropped_names) result(out)
type(r_tibble_real_t), intent(in) :: tbl
character(len=*), intent(in) :: dropped_names(:)
type(r_tibble_real_t) :: out
character(len=:), allocatable :: kept_names(:)
logical, allocatable :: keep(:)
integer :: i, j, name_len

allocate(keep(tibble_ncol(tbl)), source=.true.)
do i = 1, size(dropped_names)
   do j = 1, tibble_ncol(tbl)
      if (trim(tbl%names(j)) == trim(dropped_names(i))) keep(j) = .false.
   end do
end do
name_len = 1
if (allocated(tbl%names)) name_len = max(name_len, len(tbl%names))
allocate(character(len=name_len) :: kept_names(count(keep)))
kept_names = pack(tbl%names, keep)
out = tibble_real_select(tbl, kept_names)
end function tibble_real_drop

function tibble_real_mutate_vector(tbl, name, values) result(out)
type(r_tibble_real_t), intent(in) :: tbl
character(len=*), intent(in) :: name
real(kind=dp), intent(in) :: values(:)
type(r_tibble_real_t) :: out
character(len=:), allocatable :: names(:)
real(kind=dp), allocatable :: cols(:,:)
integer :: j, name_len, target

if (size(values) /= tibble_nrow(tbl)) &
   error stop "tibble_real_mutate: column length must match number of rows"
if (len_trim(name) == 0) error stop "tibble_real_mutate: column name must not be empty"

target = 0
do j = 1, tibble_ncol(tbl)
   if (trim(tbl%names(j)) == trim(name)) then
      target = j
      exit
   end if
end do
if (target > 0) then
   out = tbl
   out%real_cols(:, target) = values
   return
end if

name_len = max(1, len_trim(name))
if (allocated(tbl%names)) name_len = max(name_len, len(tbl%names))
allocate(character(len=name_len) :: names(tibble_ncol(tbl) + 1))
if (tibble_ncol(tbl) > 0) names(:tibble_ncol(tbl)) = tbl%names
names(size(names)) = name
allocate(cols(tibble_nrow(tbl), tibble_ncol(tbl) + 1))
if (tibble_ncol(tbl) > 0) cols(:, :tibble_ncol(tbl)) = tbl%real_cols
cols(:, size(cols, 2)) = values
if (allocated(tbl%row_labels)) then
   out = tibble_real(names, cols, tbl%row_labels)
   if (allocated(tbl%row_label_name)) out%row_label_name = tbl%row_label_name
else
   out = tibble_real(names, cols)
end if
end function tibble_real_mutate_vector

function tibble_real_mutate_scalar(tbl, name, value) result(out)
type(r_tibble_real_t), intent(in) :: tbl
character(len=*), intent(in) :: name
real(kind=dp), intent(in) :: value
type(r_tibble_real_t) :: out
real(kind=dp), allocatable :: values(:)

allocate(values(tibble_nrow(tbl)), source=value)
out = tibble_real_mutate_vector(tbl, name, values)
end function tibble_real_mutate_scalar

function tibble_integer(names, cols, row_labels, row_label_name) result(tbl)
character(len=*), intent(in) :: names(:)
integer, intent(in) :: cols(:,:)
character(len=*), intent(in), optional :: row_labels(:)
character(len=*), intent(in), optional :: row_label_name
type(r_tibble_integer_t) :: tbl
integer :: j, k, name_len, row_label_len

if (size(names) /= size(cols, 2)) &
   error stop "tibble_integer: number of names must match number of columns"
do j = 1, size(names)
   if (len_trim(names(j)) == 0) error stop "tibble_integer: column names must not be empty"
   do k = 1, j - 1
      if (trim(names(j)) == trim(names(k))) error stop "tibble_integer: duplicate column name"
   end do
end do
name_len = 1
do j = 1, size(names)
   name_len = max(name_len, len_trim(names(j)))
end do
allocate(character(len=name_len) :: tbl%names(size(names)))
tbl%names = names
tbl%integer_cols = cols
if (present(row_labels)) then
   if (size(row_labels) /= size(cols, 1)) &
      error stop "tibble_integer: number of row labels must match number of rows"
   row_label_len = 1
   do j = 1, size(row_labels)
      row_label_len = max(row_label_len, len_trim(row_labels(j)))
   end do
   allocate(character(len=row_label_len) :: tbl%row_labels(size(row_labels)))
   tbl%row_labels = row_labels
   if (present(row_label_name)) then
      if (len_trim(row_label_name) == 0) &
         error stop "tibble_integer: row label name must not be empty"
      tbl%row_label_name = trim(row_label_name)
   end if
else if (present(row_label_name)) then
   error stop "tibble_integer: row label name requires row labels"
end if
end function tibble_integer

pure integer function tibble_integer_nrow(tbl) result(n)
type(r_tibble_integer_t), intent(in) :: tbl
n = 0
if (allocated(tbl%integer_cols)) n = size(tbl%integer_cols, 1)
end function tibble_integer_nrow

pure integer function tibble_integer_ncol(tbl) result(n)
type(r_tibble_integer_t), intent(in) :: tbl
n = 0
if (allocated(tbl%integer_cols)) n = size(tbl%integer_cols, 2)
end function tibble_integer_ncol

function tibble_integer_col(tbl, name) result(col)
type(r_tibble_integer_t), intent(in) :: tbl
character(len=*), intent(in) :: name
integer, allocatable :: col(:)
integer :: j

do j = 1, tibble_ncol(tbl)
   if (trim(tbl%names(j)) == trim(name)) then
      col = tbl%integer_cols(:, j)
      return
   end if
end do
error stop "tibble_integer_col: column not found"
end function tibble_integer_col

function tibble_integer_filter(tbl, keep) result(out)
type(r_tibble_integer_t), intent(in) :: tbl
logical, intent(in) :: keep(:)
type(r_tibble_integer_t) :: out
integer, allocatable :: cols(:,:)
integer :: j

if (size(keep) /= tibble_nrow(tbl)) &
   error stop "tibble_integer_filter: mask length must match number of rows"
allocate(cols(count(keep), tibble_ncol(tbl)))
do j = 1, tibble_ncol(tbl)
   cols(:, j) = pack(tbl%integer_cols(:, j), keep)
end do
if (allocated(tbl%row_labels)) then
   out = tibble_integer(tbl%names, cols, pack(tbl%row_labels, keep))
   if (allocated(tbl%row_label_name)) out%row_label_name = tbl%row_label_name
else
   out = tibble_integer(tbl%names, cols)
end if
end function tibble_integer_filter

function tibble_integer_select(tbl, selected_names) result(out)
type(r_tibble_integer_t), intent(in) :: tbl
character(len=*), intent(in) :: selected_names(:)
type(r_tibble_integer_t) :: out
integer, allocatable :: cols(:,:)
integer :: i, j
logical :: found

allocate(cols(tibble_nrow(tbl), size(selected_names)))
do i = 1, size(selected_names)
   found = .false.
   do j = 1, tibble_ncol(tbl)
      if (trim(tbl%names(j)) == trim(selected_names(i))) then
         cols(:, i) = tbl%integer_cols(:, j)
         found = .true.
         exit
      end if
   end do
   if (.not. found) error stop "tibble_integer_select: column not found"
end do
if (allocated(tbl%row_labels)) then
   out = tibble_integer(selected_names, cols, tbl%row_labels)
   if (allocated(tbl%row_label_name)) out%row_label_name = tbl%row_label_name
else
   out = tibble_integer(selected_names, cols)
end if
end function tibble_integer_select

function tibble_integer_drop(tbl, dropped_names) result(out)
type(r_tibble_integer_t), intent(in) :: tbl
character(len=*), intent(in) :: dropped_names(:)
type(r_tibble_integer_t) :: out
character(len=:), allocatable :: kept_names(:)
logical, allocatable :: keep(:)
integer :: i, j, name_len

allocate(keep(tibble_ncol(tbl)), source=.true.)
do i = 1, size(dropped_names)
   do j = 1, tibble_ncol(tbl)
      if (trim(tbl%names(j)) == trim(dropped_names(i))) keep(j) = .false.
   end do
end do
name_len = 1
if (allocated(tbl%names)) name_len = max(name_len, len(tbl%names))
allocate(character(len=name_len) :: kept_names(count(keep)))
kept_names = pack(tbl%names, keep)
out = tibble_integer_select(tbl, kept_names)
end function tibble_integer_drop

function tibble_integer_mutate_vector(tbl, name, values) result(out)
type(r_tibble_integer_t), intent(in) :: tbl
character(len=*), intent(in) :: name
integer, intent(in) :: values(:)
type(r_tibble_integer_t) :: out
character(len=:), allocatable :: names(:)
integer, allocatable :: cols(:,:)
integer :: j, name_len, target

if (size(values) /= tibble_nrow(tbl)) &
   error stop "tibble_integer_mutate: column length must match number of rows"
if (len_trim(name) == 0) error stop "tibble_integer_mutate: column name must not be empty"
target = 0
do j = 1, tibble_ncol(tbl)
   if (trim(tbl%names(j)) == trim(name)) then
      target = j
      exit
   end if
end do
if (target > 0) then
   out = tbl
   out%integer_cols(:, target) = values
   return
end if
name_len = max(1, len_trim(name))
if (allocated(tbl%names)) name_len = max(name_len, len(tbl%names))
allocate(character(len=name_len) :: names(tibble_ncol(tbl) + 1))
if (tibble_ncol(tbl) > 0) names(:tibble_ncol(tbl)) = tbl%names
names(size(names)) = name
allocate(cols(tibble_nrow(tbl), tibble_ncol(tbl) + 1))
if (tibble_ncol(tbl) > 0) cols(:, :tibble_ncol(tbl)) = tbl%integer_cols
cols(:, size(cols, 2)) = values
if (allocated(tbl%row_labels)) then
   out = tibble_integer(names, cols, tbl%row_labels)
   if (allocated(tbl%row_label_name)) out%row_label_name = tbl%row_label_name
else
   out = tibble_integer(names, cols)
end if
end function tibble_integer_mutate_vector

function tibble_integer_mutate_scalar(tbl, name, value) result(out)
type(r_tibble_integer_t), intent(in) :: tbl
character(len=*), intent(in) :: name
integer, intent(in) :: value
type(r_tibble_integer_t) :: out
integer, allocatable :: values(:)

allocate(values(tibble_nrow(tbl)), source=value)
out = tibble_integer_mutate_vector(tbl, name, values)
end function tibble_integer_mutate_scalar

function tibble_integer_to_real(tbl) result(out)
type(r_tibble_integer_t), intent(in) :: tbl
type(r_tibble_real_t) :: out
real(kind=dp), allocatable :: cols(:,:)

allocate(cols(tibble_nrow(tbl), tibble_ncol(tbl)))
cols = real(tbl%integer_cols, kind=dp)
if (allocated(tbl%row_labels)) then
   out = tibble_real(tbl%names, cols, tbl%row_labels)
   if (allocated(tbl%row_label_name)) out%row_label_name = tbl%row_label_name
else
   out = tibble_real(tbl%names, cols)
end if
end function tibble_integer_to_real

function tibble_integer_mutate_real_vector(tbl, name, values) result(out)
type(r_tibble_integer_t), intent(in) :: tbl
character(len=*), intent(in) :: name
real(kind=dp), intent(in) :: values(:)
type(r_tibble_real_t) :: out

out = tibble_real_mutate(tibble_integer_to_real(tbl), name, values)
end function tibble_integer_mutate_real_vector

function tibble_integer_mutate_real_scalar(tbl, name, value) result(out)
type(r_tibble_integer_t), intent(in) :: tbl
character(len=*), intent(in) :: name
real(kind=dp), intent(in) :: value
type(r_tibble_real_t) :: out

out = tibble_real_mutate(tibble_integer_to_real(tbl), name, value)
end function tibble_integer_mutate_real_scalar

function read_csv_tibble_integer(file_path, max_rows, max_cols, index_col) result(tbl)
character(len=*), intent(in) :: file_path
integer, intent(in), optional :: max_rows, max_cols
character(len=*), intent(in), optional :: index_col
type(r_tibble_integer_t) :: tbl
type(r_tibble_real_t) :: real_tbl
integer, allocatable :: cols(:,:)
real(kind=dp) :: value
integer :: i, j

real_tbl = read_csv_tibble_real(file_path, max_rows, max_cols, index_col)
allocate(cols(tibble_nrow(real_tbl), tibble_ncol(real_tbl)))
do j = 1, tibble_ncol(real_tbl)
   do i = 1, tibble_nrow(real_tbl)
      value = real_tbl%real_cols(i, j)
      if (.not. ieee_is_finite(value) .or. value < real(-huge(0), dp) .or. &
          value > real(huge(0), dp) .or. value /= anint(value)) &
         error stop "read_csv_tibble_integer: non-integer or out-of-range value"
      cols(i, j) = int(value)
   end do
end do
if (allocated(real_tbl%row_labels)) then
   tbl = tibble_integer(real_tbl%names, cols, real_tbl%row_labels)
   if (allocated(real_tbl%row_label_name)) tbl%row_label_name = real_tbl%row_label_name
else
   tbl = tibble_integer(real_tbl%names, cols)
end if
end function read_csv_tibble_integer

function read_csv_tibble_real(file_path, max_rows, max_cols, index_col) result(tbl)
character(len=*), intent(in) :: file_path
integer, intent(in), optional :: max_rows
integer, intent(in), optional :: max_cols
character(len=*), intent(in), optional :: index_col
type(r_tibble_real_t) :: tbl
character(len=:), allocatable :: names(:), selected_names(:), row_labels(:)
real(kind=dp), allocatable :: cols(:,:), selected_cols(:,:)
integer, allocatable :: selected_indices(:)
integer :: data_limit, i, index_pos, j, ncol, required_cols

if (present(max_rows)) then
   if (max_rows < 0) error stop "read_csv_tibble_real: max_rows must be nonnegative"
end if
if (present(max_cols)) then
   if (max_cols < 0) error stop "read_csv_tibble_real: max_cols must be nonnegative"
end if
names = read_csv_header_names(file_path)
if (.not. present(index_col)) then
   call read_csv_real_matrix(file_path, cols, max_rows=max_rows, max_cols=max_cols)
   ncol = size(cols, 2)
   if (ncol > size(names)) &
      error stop "read_csv_tibble_real: header and data column counts differ"
   allocate(character(len=len(names)) :: selected_names(ncol))
   if (ncol > 0) selected_names = names(:ncol)
   tbl = tibble_real(selected_names, cols)
   return
end if

index_pos = 0
do i = 1, size(names)
   if (trim(names(i)) == trim(index_col)) then
      index_pos = i
      exit
   end if
end do
if (index_pos == 0) error stop "read_csv_tibble_real: index column not found"
data_limit = size(names) - 1
if (present(max_cols)) data_limit = min(data_limit, max_cols)
allocate(selected_indices(data_limit))
j = 0
do i = 1, size(names)
   if (i == index_pos) cycle
   if (j >= data_limit) exit
   j = j + 1
   selected_indices(j) = i
end do
required_cols = index_pos
if (data_limit > 0) required_cols = max(required_cols, maxval(selected_indices))
call read_csv_real_matrix(file_path, cols, max_rows=max_rows, max_cols=required_cols)
allocate(character(len=len(names)) :: selected_names(data_limit))
allocate(selected_cols(size(cols, 1), data_limit))
do j = 1, data_limit
   selected_names(j) = names(selected_indices(j))
   selected_cols(:, j) = cols(:, selected_indices(j))
end do
row_labels = read_csv_text_column(file_path, index_pos, max_rows=max_rows)
if (size(row_labels) /= size(selected_cols, 1)) &
   error stop "read_csv_tibble_real: index and data row counts differ"
tbl = tibble_real(selected_names, selected_cols, row_labels, row_label_name=trim(index_col))
end function read_csv_tibble_real

function tibble_real_log_returns(tbl, scale) result(out)
type(r_tibble_real_t), intent(in) :: tbl
real(kind=dp), intent(in), optional :: scale
type(r_tibble_real_t) :: out
real(kind=dp) :: multiplier
real(kind=dp), allocatable :: values(:,:)
character(len=:), allocatable :: labels(:)
logical, allocatable :: keep(:)
integer :: i, k

multiplier = 1.0_dp
if (present(scale)) multiplier = scale
if (tibble_nrow(tbl) < 2) then
   allocate(values(0, tibble_ncol(tbl)))
   if (allocated(tbl%row_labels)) then
      allocate(character(len=len(tbl%row_labels)) :: labels(0))
      out = tibble_real(tbl%names, values, labels)
      if (allocated(tbl%row_label_name)) out%row_label_name = tbl%row_label_name
   else
      out = tibble_real(tbl%names, values)
   end if
   return
end if
allocate(keep(tibble_nrow(tbl) - 1))
do i = 1, size(keep)
   keep(i) = all(ieee_is_finite(tbl%real_cols(i:i + 1, :))) .and. &
      all(tbl%real_cols(i:i + 1, :) > 0.0_dp)
end do
allocate(values(count(keep), tibble_ncol(tbl)))
if (allocated(tbl%row_labels)) &
   allocate(character(len=len(tbl%row_labels)) :: labels(count(keep)))
k = 0
do i = 1, size(keep)
   if (.not. keep(i)) cycle
   k = k + 1
   values(k, :) = multiplier * log(tbl%real_cols(i + 1, :) / tbl%real_cols(i, :))
   if (allocated(tbl%row_labels)) labels(k) = tbl%row_labels(i + 1)
end do
if (allocated(tbl%row_labels)) then
   out = tibble_real(tbl%names, values, labels)
   if (allocated(tbl%row_label_name)) out%row_label_name = tbl%row_label_name
else
   out = tibble_real(tbl%names, values)
end if
end function tibble_real_log_returns

function tibble_real_stats(tbl) result(out)
type(r_tibble_real_t), intent(in) :: tbl
type(r_tibble_real_t) :: out
real(kind=dp), allocatable :: values(:,:)
integer :: j

allocate(values(5, tibble_ncol(tbl)))
values(1, :) = real(tibble_nrow(tbl), kind=dp)
if (tibble_nrow(tbl) == 0) then
   values(2:5, :) = ieee_value(0.0_dp, ieee_quiet_nan)
else
   do j = 1, tibble_ncol(tbl)
      values(2, j) = sum(tbl%real_cols(:, j)) / real(tibble_nrow(tbl), kind=dp)
      values(3, j) = sd(tbl%real_cols(:, j))
      values(4, j) = minval(tbl%real_cols(:, j))
      values(5, j) = maxval(tbl%real_cols(:, j))
   end do
end if
out = tibble_real(tbl%names, values, &
   [character(len=7) :: "n", "mean", "sd", "minimum", "maximum"], &
   row_label_name="statistic")
end function tibble_real_stats

subroutine print_tibble_real(tbl, n, integer_row_labels, decimal_places, row_numbers, tibble_style)
type(r_tibble_real_t), intent(in) :: tbl
integer, intent(in), optional :: n
character(len=*), intent(in), optional :: integer_row_labels(:)
integer, intent(in), optional :: decimal_places
logical, intent(in), optional :: row_numbers
logical, intent(in), optional :: tibble_style
integer :: digits, field_width, i, j, k, nshow
logical :: found, print_row_numbers, show_row_labels, use_tibble_style
logical, allocatable :: integer_format(:), scientific_format(:)
real(kind=dp) :: value
character(len=32) :: fixed_fmt, header_fmt, integer_fmt, scientific_fmt

nshow = min(10, tibble_nrow(tbl))
if (present(n)) nshow = min(tibble_nrow(tbl), max(0, n))
print_row_numbers = .true.
if (present(row_numbers)) print_row_numbers = row_numbers
use_tibble_style = .true.
if (present(tibble_style)) use_tibble_style = tibble_style
digits = 6
if (present(decimal_places)) digits = decimal_places
if (digits < 0 .or. digits > 15) &
   error stop "print_tibble: decimal_places must be between 0 and 15"
field_width = max(12, digits + 8)
write(header_fmt, '("(1x,a",i0,")")') field_width
write(integer_fmt, '("(1x,i",i0,")")') field_width
write(fixed_fmt, '("(1x,f",i0,".",i0,")")') field_width, digits
write(scientific_fmt, '("(1x,es",i0,".",i0,")")') field_width, digits
show_row_labels = allocated(tbl%row_labels)
allocate(integer_format(tibble_nrow(tbl)), source=.false.)
if (present(integer_row_labels)) then
   if (.not. show_row_labels .and. size(integer_row_labels) > 0) &
      error stop "print_tibble: integer row labels require tibble row labels"
   do k = 1, size(integer_row_labels)
      found = .false.
      do i = 1, tibble_nrow(tbl)
         if (trim(tbl%row_labels(i)) == trim(integer_row_labels(k))) then
            found = .true.
            integer_format(i) = .true.
         end if
      end do
      if (.not. found) error stop "print_tibble: integer row label not found"
   end do
end if
do i = 1, tibble_nrow(tbl)
   if (.not. integer_format(i)) cycle
   do j = 1, tibble_ncol(tbl)
      value = tbl%real_cols(i, j)
      if (.not. ieee_is_finite(value) .or. &
          abs(value - anint(value)) > 100.0_dp * epsilon(value) * max(1.0_dp, abs(value)) .or. &
          abs(value) > real(huge(0), kind=dp)) &
         error stop "print_tibble: integer-formatted row contains a non-integer value"
   end do
end do
allocate(scientific_format(tibble_ncol(tbl)), source=.false.)
do j = 1, tibble_ncol(tbl)
   do i = 1, nshow
      if (integer_format(i)) cycle
      value = tbl%real_cols(i, j)
      if (ieee_is_finite(value) .and. value /= 0.0_dp) then
         if (abs(value) < 1.0e-4_dp .or. abs(value) >= 1.0e4_dp) &
            scientific_format(j) = .true.
      end if
   end do
end do
if (use_tibble_style) &
   write(*, '(a, i0, a, i0)') "# A tibble: ", tibble_nrow(tbl), " x ", tibble_ncol(tbl)
if (tibble_ncol(tbl) == 0) return
if (print_row_numbers) write(*, '(6x)', advance='no')
if (show_row_labels) then
   if (allocated(tbl%row_label_name)) then
      write(*, '(1x, a12)', advance='no') trim(tbl%row_label_name)
   else
      write(*, '(1x, a12)', advance='no') ".row"
   end if
end if
do j = 1, tibble_ncol(tbl)
   write(*, header_fmt, advance='no') trim(tbl%names(j))
end do
write(*, *)
if (use_tibble_style) then
   if (print_row_numbers) write(*, '(6x)', advance='no')
   if (show_row_labels) write(*, '(1x, a12)', advance='no') ""
   do j = 1, tibble_ncol(tbl)
      write(*, header_fmt, advance='no') '<dbl>'
   end do
   write(*, *)
end if
do i = 1, nshow
   if (print_row_numbers) write(*, '(i6)', advance='no') i
   if (show_row_labels) write(*, '(1x, a12)', advance='no') trim(tbl%row_labels(i))
   do j = 1, tibble_ncol(tbl)
      value = tbl%real_cols(i, j)
      if (integer_format(i)) then
         write(*, integer_fmt, advance='no') nint(value)
      else if (scientific_format(j)) then
         write(*, scientific_fmt, advance='no') value
      else
         write(*, fixed_fmt, advance='no') value
      end if
   end do
   write(*, *)
end do
if (use_tibble_style .and. nshow < tibble_nrow(tbl)) &
   write(*, '(a, i0, a)') "# ... with ", tibble_nrow(tbl) - nshow, " more rows"
end subroutine print_tibble_real

subroutine print_tibble_integer(tbl, n, integer_row_labels, decimal_places, row_numbers, tibble_style)
type(r_tibble_integer_t), intent(in) :: tbl
integer, intent(in), optional :: n
character(len=*), intent(in), optional :: integer_row_labels(:)
integer, intent(in), optional :: decimal_places
logical, intent(in), optional :: row_numbers
logical, intent(in), optional :: tibble_style
integer :: field_width, i, j, nshow
logical :: print_row_numbers, show_row_labels, use_tibble_style
character(len=32) :: header_fmt, integer_fmt

nshow = min(10, tibble_nrow(tbl))
if (present(n)) nshow = min(tibble_nrow(tbl), max(0, n))
print_row_numbers = .true.
if (present(row_numbers)) print_row_numbers = row_numbers
use_tibble_style = .true.
if (present(tibble_style)) use_tibble_style = tibble_style
if (present(decimal_places)) then
   if (decimal_places < 0 .or. decimal_places > 15) &
      error stop "print_tibble: decimal_places must be between 0 and 15"
end if
if (present(integer_row_labels)) then
   if (.not. allocated(tbl%row_labels) .and. size(integer_row_labels) > 0) &
      error stop "print_tibble: integer row labels require tibble row labels"
end if
field_width = 12
write(header_fmt, '("(1x,a",i0,")")') field_width
write(integer_fmt, '("(1x,i",i0,")")') field_width
show_row_labels = allocated(tbl%row_labels)
if (use_tibble_style) &
   write(*, '(a, i0, a, i0)') "# A tibble: ", tibble_nrow(tbl), " x ", tibble_ncol(tbl)
if (tibble_ncol(tbl) == 0) return
if (print_row_numbers) write(*, '(6x)', advance='no')
if (show_row_labels) then
   if (allocated(tbl%row_label_name)) then
      write(*, '(1x, a12)', advance='no') trim(tbl%row_label_name)
   else
      write(*, '(1x, a12)', advance='no') ".row"
   end if
end if
do j = 1, tibble_ncol(tbl)
   write(*, header_fmt, advance='no') trim(tbl%names(j))
end do
write(*, *)
if (use_tibble_style) then
   if (print_row_numbers) write(*, '(6x)', advance='no')
   if (show_row_labels) write(*, '(1x, a12)', advance='no') ""
   do j = 1, tibble_ncol(tbl)
      write(*, header_fmt, advance='no') '<int>'
   end do
   write(*, *)
end if
do i = 1, nshow
   if (print_row_numbers) write(*, '(i6)', advance='no') i
   if (show_row_labels) write(*, '(1x, a12)', advance='no') trim(tbl%row_labels(i))
   do j = 1, tibble_ncol(tbl)
      write(*, integer_fmt, advance='no') tbl%integer_cols(i, j)
   end do
   write(*, *)
end do
if (use_tibble_style .and. nshow < tibble_nrow(tbl)) &
   write(*, '(a, i0, a)') "# ... with ", tibble_nrow(tbl) - nshow, " more rows"
end subroutine print_tibble_integer

function dataframe_real_col(df, name) result(col)
type(r_dataframe_t), intent(in) :: df
character(len=*), intent(in) :: name
real(kind=dp), allocatable :: col(:)
integer :: j
do j = 1, size(df%names)
   if (trim(df%names(j)) == trim(name)) then
      col = df%real_cols(:, j)
      return
   end if
end do
error stop "dataframe column not found"
end function dataframe_real_col

subroutine print_dataframe(df, n)
type(r_dataframe_t), intent(in) :: df
integer, intent(in), optional :: n
integer :: i, j, nshow
if (.not. allocated(df%real_cols)) then
   write(*,"(a)") "data frame with 0 columns"
   return
end if
nshow = size(df%real_cols, 1)
if (present(n)) nshow = min(nshow, max(0, n))
write(*,"(a)", advance="no") "      "
do j = 1, size(df%real_cols, 2)
   write(*,"(1x,a12)", advance="no") trim(df%names(j))
end do
write(*,*)
do i = 1, nshow
   write(*,"(i6)", advance="no") i
   do j = 1, size(df%real_cols, 2)
      write(*,"(1x,g12.6)", advance="no") df%real_cols(i, j)
   end do
   write(*,*)
end do
end subroutine print_dataframe

subroutine print_dataframe_head(df, n)
type(r_dataframe_t), intent(in) :: df
integer, intent(in), optional :: n
if (present(n)) then
   call print_dataframe(df, n)
else
   call print_dataframe(df, 6)
end if
end subroutine print_dataframe_head

subroutine print_aggregate_result(x)
type(aggregate_result_t), intent(in) :: x
integer :: i
if (.not. allocated(x%labels) .or. .not. allocated(x%values)) then
   write(*,"(a)") "aggregate result with 0 rows"
   return
end if
write(*,"(a,1x,a)") trim(x%group_name), trim(x%value_name)
do i = 1, min(size(x%labels), size(x%values))
   write(*,"(i0,1x,a,1x,g0)") i, trim(x%labels(i)), x%values(i)
end do
end subroutine print_aggregate_result

subroutine print_by_matrix_result(x)
type(by_matrix_result_t), intent(in) :: x
integer :: i, j
if (.not. allocated(x%labels) .or. .not. allocated(x%values)) then
   write(*,"(a)") "by result with 0 groups"
   return
end if
do i = 1, min(size(x%labels), size(x%values, 1))
   write(*,"(a,1x,a)") "INDICES:", trim(x%labels(i))
   do j = 1, size(x%values, 2)
      write(*,"(a,i0,1x)", advance="no") "V", j
   end do
   write(*,*)
   do j = 1, size(x%values, 2)
      write(*,"(g0,1x)", advance="no") x%values(i, j)
   end do
   write(*,*)
   if (i < min(size(x%labels), size(x%values, 1))) write(*, &
   & "(a)") "------------------------------------------------------------"
end do
end subroutine print_by_matrix_result

end module r_mod
