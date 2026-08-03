program test_jll
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use multiatsm_kinds, only : dp
  use multiatsm_types, only : jll_model
  use multiatsm_jll, only : fit_jll, jll_feedback_restrictions, jll_cholesky_mask
  implicit none
  integer, parameter :: g = 1, m = 1, n = 1, c = 2, nt = 180, k = g + c * (m + n)
  real(dp) :: factors(k, nt), global, macro1, macro2, price1, price2
  real(dp), allocatable :: restrictions(:, :)
  logical, allocatable :: feedback_free(:, :), chol_free(:, :)
  type(jll_model) :: model
  integer :: t, info

  do t = 1, nt
    global = 0.03_dp * sin(0.09_dp * real(t, dp))
    macro1 = 0.6_dp * global + 0.02_dp * cos(0.17_dp * real(t, dp))
    price1 = 0.7_dp * macro1 + 0.015_dp * sin(0.23_dp * real(t, dp))
    macro2 = 0.4_dp * global + 0.3_dp * macro1 + 0.018_dp * cos(0.19_dp * real(t, dp))
    price2 = 0.5_dp * macro2 + 0.25_dp * price1 + 0.012_dp * sin(0.29_dp * real(t, dp))
    factors(:, t) = [global, macro1, price1, macro2, price2]
  end do
  call jll_feedback_restrictions(g, m, n, c, 1, restrictions, feedback_free, info)
  call check(info == 0 .and. size(restrictions, 1) == k, 'feedback restrictions')
  call jll_cholesky_mask(g, m, n, c, 1, chol_free, info)
  call check(info == 0 .and. all(.not. chol_free(3, 1:2)), 'Cholesky restrictions')

  call fit_jll(factors, g, m, n, c, 1, model, info)
  call check(info == 0, 'JLL fit status')
  call check(all(ieee_is_finite(model%k1)), 'finite JLL feedback')
  call check(maxval(abs(model%sigma_nonortho - transpose(model%sigma_nonortho))) < 1.0e-12_dp, &
    'JLL covariance symmetry')
  call check(maxval(abs(matmul(model%pi_matrix, model%orthogonal_factors) - factors)) < 0.03_dp, &
    'JLL factor reconstruction')
  call check(all(.not. model%chol_free(3, 1:2)), 'stored Cholesky mask')
  print '(a)', 'test_jll: PASS'
contains
  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*, '(a)') 'FAIL: ' // message
      error stop 1
    end if
  end subroutine check
end program test_jll
