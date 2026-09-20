program test_rlrsim
  use rlrsim_api, only : dp, exact_test_result, lrt_sim, rlrsim_result, rlrt_sim, exact_lrt_from_design, exact_rlrt_from_design
  use rlrsim_numerics, only : singular_values_squared
  implicit none

  real(dp) :: x(8, 2)
  real(dp) :: z(8, 2)
  real(dp) :: s(2, 2)
  real(dp) :: a(2, 2)
  real(dp) :: xb(8, 1)
  real(dp) :: xrank(8, 2)
  real(dp) :: zb(8, 3)
  real(dp) :: zd(8, 3)
  real(dp) :: zt(8, 3)
  real(dp) :: s3(3, 3)
  real(dp), allocatable :: sv2(:)
  type(rlrsim_result) :: r1
  type(rlrsim_result) :: r2
  type(rlrsim_result) :: rr
  type(rlrsim_result) :: rr_general
  type(rlrsim_result) :: rr_repeat
  type(exact_test_result) :: et
  integer :: i

  a = 0.0_dp
  a(1, 1) = 3.0_dp
  a(2, 2) = 2.0_dp
  call singular_values_squared(a, sv2)
  call check(size(sv2) == 2, "singular-value count")
  call check(abs(sv2(1) - 9.0_dp) < 1.0e-10_dp, "largest squared singular value")
  call check(abs(sv2(2) - 4.0_dp) < 1.0e-10_dp, "second squared singular value")

  do i = 1, 8
    x(i, 1) = 1.0_dp
    x(i, 2) = real(i, dp) - 4.5_dp
  end do
  z = 0.0_dp
  z(1:4, 1) = 1.0_dp
  z(5:8, 2) = 1.0_dp
  s = 0.0_dp
  s(1, 1) = 1.0_dp
  s(2, 2) = 1.0_dp

  call lrt_sim(x, z, 1, s, r1, seed = 2468, nsim = 128, gridlength = 40)
  call lrt_sim(x, z, 1, s, r2, seed = 2468, nsim = 128, gridlength = 40)
  call check(size(r1%statistic) == 128, "LRTSim sample size")
  call check(all(r1%statistic >= 0.0_dp), "LRTSim nonnegative statistics")
  call check(all(r1%lambda >= 0.0_dp), "LRTSim nonnegative lambdas")
  call check(all(r1%statistic == r2%statistic), "LRTSim deterministic seeded stream")
  call check(all(r1%lambda == r2%lambda), "LRTSim deterministic seeded lambdas")

  call rlrt_sim(x, z, s, rr, seed = 1357, nsim = 128, gridlength = 40)
  call check(size(rr%statistic) == 128, "RLRTSim sample size")
  call check(all(rr%statistic >= 0.0_dp), "RLRTSim nonnegative statistics")
  call check(all(rr%lambda >= 0.0_dp), "RLRTSim nonnegative lambdas")

  call rlrt_sim(x, z, s, rr, lambda0 = 0.5_dp, seed = 1357, nsim = 32, gridlength = 40)
  call check(size(rr%statistic) == 32, "RLRTSim positive-lambda0 sample size")
  call check(all(rr%lambda >= 0.0_dp), "RLRTSim positive-lambda0 grid")

  xb(:, 1) = 1.0_dp
  xrank(:, 1) = 1.0_dp
  xrank(:, 2) = 1.0_dp
  s3 = 0.0_dp
  s3(1, 1) = 1.0_dp
  s3(2, 2) = 1.0_dp
  s3(3, 3) = 1.0_dp

  zb = 0.0_dp
  zb(1, 1) = 1.0_dp
  zb(2, 1) = -1.0_dp
  zb(3, 2) = 1.0_dp
  zb(4, 2) = -1.0_dp
  call rlrt_sim(xb, zb, s3, rr, seed = 31415, nsim = 96, use_approx = 0.9_dp, gridlength = 40)
  call rlrt_sim(xb, zb, s3, rr_repeat, seed = 31415, nsim = 96, use_approx = 0.9_dp, gridlength = 40)
  call rlrt_sim(xb, zb, s3, rr_general, seed = 31415, nsim = 96, use_approx = 0.0_dp, gridlength = 40)
  call check(all(rr%statistic >= 0.0_dp), "balanced-ANOVA approximation nonnegative statistics")
  call check(all(rr%lambda >= 0.0_dp), "balanced-ANOVA approximation nonnegative lambdas")
  call check(all(rr%statistic == rr_repeat%statistic), "balanced-ANOVA approximation deterministic")
  call check(any(rr%statistic /= rr_general%statistic), "balanced-ANOVA specialized branch differs from grid simulation")

  zd = 0.0_dp
  zd(1, 1) = 1.0_dp
  zd(2, 1) = -1.0_dp
  zd(3, 2) = sqrt(0.05_dp)
  zd(4, 2) = -sqrt(0.05_dp)
  zd(5, 3) = 0.1_dp
  zd(6, 3) = -0.1_dp
  call rlrt_sim(xb, zd, s3, rr, seed = 27182, nsim = 96, use_approx = 0.9_dp, gridlength = 40)
  call rlrt_sim(xb, zd, s3, rr_repeat, seed = 27182, nsim = 96, use_approx = 0.9_dp, gridlength = 40)
  call rlrt_sim(xb, zd, s3, rr_general, seed = 27182, nsim = 96, use_approx = 0.0_dp, gridlength = 40)
  call check(all(rr%statistic >= 0.0_dp), "dominating-eigenvalue approximation nonnegative statistics")
  call check(all(rr%lambda >= 0.0_dp), "dominating-eigenvalue approximation nonnegative lambdas")
  call check(all(rr%statistic == rr_repeat%statistic), "dominating-eigenvalue approximation deterministic")
  call check(any(rr%statistic /= rr_general%statistic), "dominating-eigenvalue branch differs from grid simulation")

  zt = 0.0_dp
  zt(1, 1) = 1.0_dp
  zt(2, 1) = -1.0_dp
  zt(3, 2) = sqrt(0.4_dp)
  zt(4, 2) = -sqrt(0.4_dp)
  zt(5, 3) = sqrt(0.1_dp)
  zt(6, 3) = -sqrt(0.1_dp)
  call rlrt_sim(xb, zt, s3, rr, seed = 16180, nsim = 48, use_approx = 0.8_dp, gridlength = 40)
  call check(size(rr%statistic) == 48, "general eigenvalue-truncation approximation sample size")
  call check(all(rr%statistic >= 0.0_dp), "general eigenvalue-truncation approximation nonnegative")

  call rlrt_sim(xrank, zt, s3, rr, seed = 4242, nsim = 48, gridlength = 40)
  call check(size(rr%statistic) == 48, "rank-deficient fixed-effect design accepted")
  call check(all(rr%statistic >= 0.0_dp), "rank-deficient fixed-effect design nonnegative")

  call exact_lrt_from_design(x, z, 1, s, 1.0_dp, et, seed = 999, nsim = 64, gridlength = 30)
  call check(et%statistic == 1.0_dp, "exactLRT observed statistic")
  call check(et%p_value >= 0.0_dp .and. et%p_value <= 1.0_dp, "exactLRT p-value range")
  call check(size(et%sample) == 64, "exactLRT retained null sample")

  call exact_rlrt_from_design(x, z, s, 0.0_dp, et, seed = 999, nsim = 64, gridlength = 30)
  call check(et%p_value == 1.0_dp, "exactRLRT zero-statistic p-value")
  call check(size(et%sample) == 0, "exactRLRT zero-statistic skips simulation")

  print '(a)', "All RLRsim tests passed."

contains

  subroutine check(condition, message)
    logical, intent(in) :: condition !! True when the tested invariant holds.
    character(len=*), intent(in) :: message !! Short diagnostic naming the tested invariant.

    if (.not. condition) then
      print '(a)', "FAIL: " // message
      error stop 1
    end if
  end subroutine check

end program test_rlrsim
