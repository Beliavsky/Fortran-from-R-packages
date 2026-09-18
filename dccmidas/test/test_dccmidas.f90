program test_dccmidas
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use dccmidas
   implicit none

   integer, parameter :: nt = 60
   integer, parameter :: k = 3
   real(dp) :: ret(nt, k)
   real(dp) :: res(k, nt)
   real(dp) :: dt(k, k, nt)
   real(dp) :: a2(2, 2)
   real(dp) :: determinant_value
   real(dp), allocatable :: inverse(:, :)
   real(dp), allocatable :: ll(:)
   real(dp), allocatable :: h(:, :, :)
   real(dp), allocatable :: loss(:)
   real(dp), allocatable :: se(:)
   real(dp) :: scalar_param(8)
   real(dp) :: diag_param(12)
   real(dp) :: hessian(2, 2)
   real(dp) :: scores(3, 2)
   type(dcc_matrices) :: matrices
   type(dcc_fit_result) :: fit_result
   integer :: i
   integer :: t
   integer :: status

   do t = 1, nt
      ret(t, 1) = 0.010_dp * sin(0.17_dp * real(t, dp)) + 0.003_dp * cos(0.07_dp * real(t, dp))
      ret(t, 2) = 0.009_dp * cos(0.13_dp * real(t, dp)) + 0.002_dp * sin(0.11_dp * real(t, dp))
      ret(t, 3) = 0.008_dp * sin(0.09_dp * real(t, dp) + 0.4_dp) - 0.002_dp * cos(0.05_dp * real(t, dp))
   end do
   do i = 1, k
      res(i, :) = ret(:, i) / 0.01_dp
   end do
   dt = 0.0_dp
   do t = 1, nt
      do i = 1, k
         dt(i, i, t) = 0.01_dp
      end do
   end do

   a2 = reshape([2.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], [2, 2])
   call det_matrix(a2, determinant_value, status)
   call assert_true(status == 0 .and. abs(determinant_value - 3.0_dp) < 1.0e-12_dp, 'Det')
   call inv_matrix(a2, inverse, status)
   call assert_true(status == 0 .and. maxval(abs(matmul(a2, inverse) - identity_matrix_test(2))) < 1.0e-11_dp, 'Inv')

   call riskmetrics_mat(ret, h, status, lambda=0.94_dp)
   call assert_true(status == 0 .and. all(ieee_is_finite(h)), 'riskmetrics_mat')
   call moving_cov(ret, 10, h, status)
   call assert_true(status == 0 .and. all(ieee_is_finite(h)), 'moving_cov')

   scalar_param = [0.20_dp, 0.00_dp, 0.20_dp, 0.00_dp, 0.00_dp, 0.20_dp, 0.10_dp, 0.80_dp]
   call sBEKK_loglik(scalar_param, ret, ll, status)
   call assert_true(status == 0 .and. all(ieee_is_finite(ll)), 'sBEKK_loglik')
   call sBEKK_mat_est(scalar_param, ret, h, status)
   call assert_true(status == 0 .and. all(ieee_is_finite(h)), 'sBEKK_mat_est')

   diag_param = [0.20_dp, 0.00_dp, 0.20_dp, 0.00_dp, 0.00_dp, 0.20_dp, &
      0.10_dp, 0.10_dp, 0.10_dp, 0.80_dp, 0.80_dp, 0.80_dp]
   call dBEKK_loglik(diag_param, ret, ll, status)
   call assert_true(status == 0 .and. all(ieee_is_finite(ll)), 'dBEKK_loglik')
   call dBEKK_mat_est(diag_param, ret, h, status)
   call assert_true(status == 0 .and. all(ieee_is_finite(h)), 'dBEKK_mat_est')

   call dcc_loglik([0.02_dp, 0.90_dp], res, ll, status)
   call assert_true(status == 0 .and. all(ieee_is_finite(ll)), 'dcc_loglik')
   call dcc_mat_est([0.02_dp, 0.90_dp], res, dt, matrices, status)
   call assert_true(status == 0 .and. all(ieee_is_finite(matrices%h)), 'dcc_mat_est')

   call a_dcc_loglik([0.02_dp, 0.90_dp, 0.005_dp], res, ll, status)
   call assert_true(status == 0 .and. all(ieee_is_finite(ll)), 'a_dcc_loglik')
   call a_dcc_mat_est([0.02_dp, 0.90_dp, 0.005_dp], res, dt, matrices, status)
   call assert_true(status == 0 .and. all(ieee_is_finite(matrices%h)), 'a_dcc_mat_est')

   call dccmidas_loglik([0.02_dp, 0.90_dp, 2.0_dp], res, 'Beta', 4, 5, ll, status)
   call assert_true(status == 0 .and. all(ieee_is_finite(ll)), 'dccmidas_loglik')
   call dccmidas_mat_est([0.02_dp, 0.90_dp, 2.0_dp], res, dt, 'Beta', 4, 5, matrices, status)
   call assert_true(status == 0 .and. all(ieee_is_finite(matrices%h)), 'dccmidas_mat_est')

   call a_dccmidas_loglik([0.02_dp, 0.90_dp, 0.002_dp, 2.0_dp], res, 'Beta', 4, 5, ll, status)
   call assert_true(status == 0 .and. all(ieee_is_finite(ll)), 'a_dccmidas_loglik')
   call a_dccmidas_mat_est([0.02_dp, 0.90_dp, 0.002_dp, 2.0_dp], res, dt, 'Beta', 4, 5, matrices, status)
   call assert_true(status == 0 .and. all(ieee_is_finite(matrices%h)), 'a_dccmidas_mat_est')

   call deco_loglik([0.02_dp, 0.90_dp], res, ll, status)
   call assert_true(status == 0 .and. all(ieee_is_finite(ll)), 'deco_loglik')
   call deco_mat_est([0.02_dp, 0.90_dp], res, dt, matrices, status)
   call assert_true(status == 0 .and. all(ieee_is_finite(matrices%h)), 'deco_mat_est')

   call cov_eval(matrices%h, 'FROB', loss, status, cov_proxy=matrices%h)
   call assert_true(status == 0 .and. maxval(abs(loss)) < 1.0e-15_dp, 'cov_eval FROB')
   call cov_eval(matrices%h, 'SFROB', loss, status, cov_proxy=matrices%h)
   call assert_true(status == 0 .and. maxval(abs(loss)) < 1.0e-15_dp, 'cov_eval SFROB')
   call cov_eval(matrices%h, 'EUCL', loss, status, cov_proxy=matrices%h)
   call assert_true(status == 0 .and. maxval(abs(loss)) < 1.0e-15_dp, 'cov_eval EUCL')
   call cov_eval(matrices%h, 'RMSE', loss, status, cov_proxy=matrices%h)
   call assert_true(status == 0 .and. maxval(abs(loss)) < 1.0e-15_dp, 'cov_eval RMSE')
   call cov_eval(matrices%h, 'QLIKE', loss, status, cov_proxy=matrices%h)
   call assert_true(status == 0 .and. all(ieee_is_finite(loss)), 'cov_eval QLIKE')

   hessian = 0.0_dp
   hessian(1, 1) = -1.0_dp
   hessian(2, 2) = -1.0_dp
   scores = reshape([1.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, 2.0_dp, 0.0_dp], [3, 2])
   call qmle_sd(hessian, scores, se, status)
   call assert_true(status == 0 .and. all(ieee_is_finite(se)), 'QMLE_sd')

   call dcc_fit(ret, 'not-a-model', 'norm', 'cDCC', fit_result, status)
   call assert_true(status == DCCMIDAS_INVALID_INPUT, 'dcc_fit input validation')

   print '(a)', 'All dccmidas tests passed.'

contains

   subroutine assert_true(condition, label)
      logical, intent(in) :: condition !! Condition that must be true for the deterministic test to pass.
      character(len=*), intent(in) :: label !! Human-readable name of the tested behavior.
      if (.not. condition) then
         print '(a,1x,a)', 'FAILED:', trim(label)
         error stop 1
      end if
   end subroutine assert_true

   pure function identity_matrix_test(n) result(a)
      integer, intent(in) :: n !! Order of the identity matrix used only in this test program.
      real(dp) :: a(n, n)
      integer :: j
      a = 0.0_dp
      do j = 1, n
         a(j, j) = 1.0_dp
      end do
   end function identity_matrix_test

end program test_dccmidas
