program jll_example
  use multiatsm, only : dp, jll_model, fit_jll
  implicit none
  integer, parameter :: nt = 240
  real(dp) :: factors(5, nt), g, m1, p1, m2, p2
  type(jll_model) :: model
  integer :: t, info

  do t = 1, nt
    g = 0.02_dp * sin(0.08_dp * real(t, dp))
    m1 = 0.5_dp * g + 0.03_dp * cos(0.13_dp * real(t, dp))
    p1 = 0.8_dp * m1 + 0.02_dp * sin(0.17_dp * real(t, dp))
    m2 = 0.4_dp * g + 0.25_dp * m1 + 0.02_dp * cos(0.21_dp * real(t, dp))
    p2 = 0.6_dp * m2 + 0.2_dp * p1 + 0.015_dp * sin(0.29_dp * real(t, dp))
    factors(:, t) = [g, m1, p1, m2, p2]
  end do
  call fit_jll(factors, 1, 1, 1, 2, 1, model, info)
  if (info /= 0) error stop 'fit_jll failed'
  write(*, '(a)') 'Non-orthogonal JLL feedback matrix:'
  write(*, '(5f10.5)') model%k1
  write(*, '(a,f12.7)') 'Reconstruction max error: ', &
    maxval(abs(matmul(model%pi_matrix, model%orthogonal_factors) - factors))
end program jll_example
