program test_probability
 use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf
 use truncated_normal
 use r_compat, only: set_seed_int
 implicit none
 real(dp) :: mu2(2), sig2(2,2), lo2(2), hi2(2)
 real(dp), allocatable :: sig15(:,:), lo15(:), hi15(:), mu15(:)
 real(dp) :: inf
 integer :: i, j
 type(prob_result) :: pr

 inf = ieee_value(1.0_dp, ieee_positive_inf)
 call set_seed_int(202)

 mu2 = 0.0_dp
 sig2 = 0.0_dp
 sig2(1,1) = 1.0_dp
 sig2(2,2) = 1.0_dp
 lo2 = 0.0_dp
 hi2 = inf
 pr = mvncdf(lo2, hi2, sig2, 4000)
 call check_close('independent normal orthant', pr%prob, 0.25_dp, 1.0e-12_dp)
 pr = mvnqmc(lo2, hi2, sig2, 4096)
 call check_close('independent normal QMC orthant', pr%prob, 0.25_dp, 1.0e-12_dp)

 mu2 = [0.2_dp, -0.1_dp]
 sig2 = reshape([1.0_dp, 0.4_dp, 0.4_dp, 1.0_dp], [2,2])
 lo2 = [0.0_dp, -1.0_dp]
 hi2 = [2.0_dp, 1.5_dp]
 pr = pmvnorm(mu2, sig2, lo2, hi2, 30000, .true.)
 call check_close('2-D normal rectangle', pr%prob, 0.44175424561021764_dp, 1.5e-3_dp)
 pr = pmvt(mu2, sig2, 5.0_dp, lo2, hi2, 36000, .true.)
 call check_close('2-D Student rectangle', pr%prob, 0.39766111423138145_dp, 5.0e-3_dp)

 allocate(sig15(15,15), lo15(15), hi15(15), mu15(15))
 sig15 = 0.5_dp
 do i = 1, 15
  sig15(i,i) = 1.0_dp
 end do
 lo15 = 0.0_dp
 hi15 = inf
 mu15 = 0.0_dp
 call set_seed_int(203)
 pr = pmvnorm(mu15, sig15, lo15, hi15, 24000, .true.)
 call check_close('15-D rho=0.5 orthant', pr%prob, 1.0_dp/16.0_dp, 4.0e-3_dp)

 print '(a)', 'test_probability: PASS'
contains
 subroutine check_close(name, got, want, tol)
  character(len=*), intent(in) :: name
  real(dp), intent(in) :: got, want, tol
  if (abs(got-want) > tol) then
   print '(a,3(1x,es24.16))', trim(name)//' failed:', got, want, tol
   error stop 1
  end if
 end subroutine check_close
end program test_probability
