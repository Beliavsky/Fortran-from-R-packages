program example_matrix_power
  use corpcor, only : dp, mpower
  implicit none
  real(dp) :: a(3, 3)
  real(dp), allocatable :: root(:, :), inverse_root(:, :)
  integer :: i

  a = reshape([4.0_dp, 1.0_dp, 0.5_dp, 1.0_dp, 3.0_dp, 0.2_dp, &
    0.5_dp, 0.2_dp, 2.0_dp], [3, 3])
  root = mpower(a, 0.5_dp)
  inverse_root = mpower(a, -0.5_dp)
  print '(a)'; print '(a)', 'A^(1/2):'
  do i = 1, 3
    print '(*(f10.6,1x))', root(i, :)
  end do
  print '(a,es12.4)', 'max |A^(1/2) A^(-1/2) - I|: ', &
    maxval(abs(matmul(root, inverse_root) - identity(3)))
contains
  pure function identity(n) result(out)
    integer, intent(in) :: n
    real(dp) :: out(n, n)
    integer :: k
    out = 0.0_dp
    do k = 1, n
      out(k, k) = 1.0_dp
    end do
  end function identity
end program example_matrix_power
