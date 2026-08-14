program test_randomized
  use irlba
  use irlba_linalg, only : svd_real
  implicit none
  integer :: case_id, m, n, k, i, j, info, nseed
  integer, allocatable :: seed(:)
  real(dp), allocatable :: a(:, :), b(:, :), s(:), u(:, :), v(:, :), ctr(:), scl(:)
  real(dp), allocatable :: rtmp(:, :)
  real(dp) :: err, rv
  type(irlba_result) :: z
  type(irlba_control) :: ctl
  type(csc_matrix) :: sp

  call random_seed(size=nseed)
  allocate(seed(nseed))
  do i = 1, nseed
    seed(i) = 137 + 97 * i
  end do
  call random_seed(put=seed)

  do case_id = 1, 60
    m = 18 + mod(7 * case_id, 23)
    n = 15 + mod(11 * case_id, 19)
    k = min(4, min(m, n) - 2)
    allocate(a(m, n), rtmp(m, n), ctr(n), scl(n))
    call random_number(rtmp)
    a = 2.0_dp * rtmp - 1.0_dp
    ctr = sum(a, dim=1) / real(m, dp)
    do j = 1, n
      call random_number(rv)
      scl(j) = 0.7_dp + rv
    end do
    ctl%work = min(min(m, n), k + 10)
    ctl%maxit = 1000
    ctl%tol = 1.0e-9_dp
    ctl%svtol = 1.0e-9_dp

    z = irlba_svd(a, k, control=ctl)
    call svd_real(a, s, u, v, info)
    if (z%info /= 0) error stop "random dense nonconvergence"
    err = maxval(abs(z%d - s(1:k)) / max(1.0_dp, s(1:k)))
    if (err > 2.0e-7_dp) error stop "random dense mismatch"
    deallocate(s, u, v)

    b = a
    do j = 1, n
      b(:, j) = (b(:, j) - ctr(j)) / scl(j)
    end do
    z = irlba_svd(a, k, control=ctl, center=ctr, scale=scl)
    call svd_real(b, s, u, v, info)
    if (z%info /= 0) error stop "random center/scale nonconvergence"
    err = maxval(abs(z%d - s(1:k)) / max(1.0_dp, s(1:k)))
    if (err > 2.0e-7_dp) error stop "random center/scale mismatch"
    deallocate(s, u, v, b)

    if (mod(case_id, 3) == 0) then
      where (abs(a) < 0.65_dp) a = 0.0_dp
      sp = csc_from_dense(a)
      z = irlba_svd(sp, k, control=ctl)
      call svd_real(a, s, u, v, info)
      if (z%info /= 0) error stop "random sparse nonconvergence"
      err = maxval(abs(z%d - s(1:k)) / max(1.0_dp, s(1:k)))
      if (err > 2.0e-7_dp) error stop "random sparse mismatch"
      deallocate(s, u, v)
    end if

    deallocate(a, rtmp, ctr, scl)
  end do

  print *, "test_randomized: PASS (60 cases)"
end program test_randomized
