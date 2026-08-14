program basic
  use irlba
  implicit none
  real(dp) :: a(100, 30)
  type(irlba_result) :: s
  type(irlba_control) :: ctl
  integer :: i, j

  do j = 1, size(a, 2)
    do i = 1, size(a, 1)
      a(i, j) = sin(0.01_dp * real(i*j, dp)) + cos(0.03_dp * real(i + j, dp))
    end do
  end do
  ctl%tol = 1.0e-8_dp
  s = irlba_svd(a, 5, control=ctl)
  print '(a,5(1x,es12.5))', 'leading singular values:', s%d
  print '(a,i0,a,i0)', 'iterations=', s%iter, ' matrix products=', s%mprod
end program basic
