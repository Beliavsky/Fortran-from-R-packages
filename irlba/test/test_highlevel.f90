program test_highlevel
  use irlba
  use irlba_linalg, only : thin_qr
  implicit none
  integer, parameter :: m = 30, n = 12
  real(dp) :: a(m, n), q(n, n), eigval(n), target(3), z(n, n)
  real(dp), allocatable :: qtmp(:, :)
  type(eigen_result) :: e
  type(pca_result) :: p
  type(ssvd_result) :: ss
  type(svdr_result) :: sr, srs
  type(csc_matrix) :: sp
  type(irlba_control) :: ctl
  integer :: i, j, ierr

  do j = 1, n
    do i = 1, m
      a(i, j) = cos(0.07_dp * real(i * (j + 1), dp)) + 0.01_dp * real(i + 2*j, dp)
    end do
  end do
  ctl%tol = 1.0e-9_dp
  ctl%svtol = 1.0e-9_dp
  ctl%work = 10

  p = prcomp_irlba(a, 3, do_center=.true., do_scale=.true., control=ctl)
  if (p%info /= 0) error stop "PCA failed"
  if (.not. allocated(p%scores)) error stop "PCA scores missing"
  if (abs(sum(p%proportion) - p%cumulative(3)) > 1.0e-14_dp) error stop "PCA summary mismatch"

  do j = 1, n
    do i = 1, n
      z(i, j) = sin(0.19_dp * real(i*j, dp)) + cos(0.07_dp * real(i + 2*j, dp))
    end do
  end do
  call thin_qr(z, qtmp, ierr)
  if (ierr /= 0) error stop "QR setup failed"
  q = qtmp
  eigval = [12.0_dp, -11.0_dp, 10.0_dp, -9.0_dp, 8.0_dp, -7.0_dp, &
            6.0_dp, -5.0_dp, 4.0_dp, -3.0_dp, 2.0_dp, -1.0_dp]
  q = matmul(q, matmul(diag_matrix(eigval), transpose(q)))
  e = partial_eigen(q, 3, control=ctl)
  target = [12.0_dp, 10.0_dp, 8.0_dp]
  if (maxval(abs(e%values - target)) > 1.0e-8_dp) error stop "partial_eigen mismatch"

  ss = ssvd(a, 2, [4, 5], maxit=100, tol=1.0e-5_dp, control=ctl)
  if (ss%info /= 0 .and. ss%info /= -2) error stop "ssvd failed"
  if (count(abs(ss%v(:, 1)) > 1.0e-12_dp) > 4) error stop "ssvd sparsity target failed"
  if (count(abs(ss%v(:, 2)) > 1.0e-12_dp) > 5) error stop "ssvd sparsity target failed"

  sr = svdr(a, 3, tol=1.0e-8_dp, it=50, extra=7)
  if (sr%info /= 0) error stop "svdr failed"
  if (maxval(abs(matmul(a, sr%v) - spread(sr%d, 1, m) * sr%u)) > 2.0e-5_dp) error stop "svdr residual too large"


  sp = csc_from_dense(a, zero_tol=0.4_dp)
  srs = svdr(sp, 3, tol=1.0e-8_dp, it=50, extra=7)
  if (srs%info /= 0) error stop "sparse svdr failed"
  if (maxval(abs(matmul(sp%to_dense(), srs%v) - spread(srs%d, 1, m) * srs%u)) > 2.0e-5_dp) then
    error stop "sparse svdr residual too large"
  end if

  print *, "test_highlevel: PASS"
contains
  function diag_matrix(x) result(d)
    real(dp), intent(in) :: x(:)
    real(dp) :: d(size(x), size(x))
    integer :: k
    d = 0.0_dp
    do k = 1, size(x)
      d(k, k) = x(k)
    end do
  end function diag_matrix
end program test_highlevel
