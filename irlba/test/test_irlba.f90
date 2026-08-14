program test_irlba
  use irlba
  use irlba_linalg, only : svd_real
  implicit none
  integer, parameter :: m = 40, n = 25, k = 5
  real(dp) :: a(m, n), scalev(n), cent(n), err
  real(dp), allocatable :: sref(:), uref(:, :), vref(:, :), b(:, :), x1(:, :), x2(:, :)
  type(irlba_result) :: r, rs, rr
  type(irlba_control) :: ctl
  type(csc_matrix) :: sp
  integer :: i, j, info

  do j = 1, n
    do i = 1, m
      a(i, j) = sin(12.9898_dp * real(i, dp) + 78.233_dp * real(j, dp)) + &
        0.37_dp * cos(4.123_dp * real(i*j, dp))
    end do
  end do
  ctl%work = k + 9
  ctl%tol = 1.0e-10_dp
  ctl%svtol = 1.0e-10_dp
  ctl%maxit = 1000
  r = irlba_svd(a, k, control=ctl)
  if (r%info /= 0) error stop "dense irlba did not converge"
  call svd_real(a, sref, uref, vref, info)
  if (info /= 0) error stop "reference SVD failed"
  err = maxval(abs(r%d - sref(1:k)))
  if (err > 1.0e-8_dp) error stop "dense singular values mismatch"
  err = maxval(abs(matmul(a, r%v) - spread(r%d, 1, m) * r%u))
  if (err > 2.0e-7_dp) error stop "dense residual too large"

  cent = sum(a, dim=1) / real(m, dp)
  do j = 1, n
    scalev(j) = 1.0_dp + 0.02_dp * real(j, dp)
  end do
  b = a
  do j = 1, n
    b(:, j) = (b(:, j) - cent(j)) / scalev(j)
  end do
  r = irlba_svd(a, k, control=ctl, center=cent, scale=scalev)
  call svd_real(b, sref, uref, vref, info)
  if (maxval(abs(r%d - sref(1:k))) > 1.0e-8_dp) error stop "center/scale mismatch"

  sp = csc_from_dense(a, zero_tol=0.25_dp)
  rs = irlba_svd(sp, k, control=ctl)
  b = sp%to_dense()
  call svd_real(b, sref, uref, vref, info)
  if (maxval(abs(rs%d - sref(1:k))) > 1.0e-8_dp) error stop "sparse singular values mismatch"

  rr = irlba_svd(a, 3, control=ctl)
  r = irlba_svd(a, 5, control=ctl, restart_from=rr)
  call svd_real(a, sref, uref, vref, info)
  if (maxval(abs(r%d - sref(1:5))) > 1.0e-8_dp) error stop "restart singular values mismatch"

  r = irlba_svd(a, 3, control=ctl, smallest=.true.)
  if (maxval(abs(r%d - sref(size(sref)-2:size(sref)))) > 1.0e-12_dp) error stop "smallest fallback mismatch"


  allocate(x1(30, 5), x2(20, 5))
  do j = 1, 5
    do i = 1, 30
      x1(i, j) = sin(0.07_dp * real(i*j, dp))
    end do
    do i = 1, 20
      x2(i, j) = cos(0.09_dp * real(i*(j + 1), dp))
    end do
  end do
  b = matmul(x1, transpose(x2))
  ctl%work = 15
  r = irlba_svd(b, 8, control=ctl)
  if (r%info /= 0) error stop "rank-deficient case failed"
  if (maxval(abs(r%d(6:8))) > 1.0e-10_dp) error stop "rank-deficient zeros not recovered"

  print *, "test_irlba: PASS"
end program test_irlba
