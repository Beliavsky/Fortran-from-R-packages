program test_sampling
 use truncated_normal
 use r_compat, only: set_seed_int
 implicit none
 real(dp) :: mu(2), sig(2,2), lo(2), hi(2), nu
 real(dp), allocatable :: x(:,:), xt(:,:), reconstructed(:,:)
 type(tregress_result) :: tr
 integer :: i

 mu = [0.2_dp, -0.1_dp]
 sig = reshape([1.0_dp, 0.4_dp, 0.4_dp, 1.0_dp], [2,2])
 lo = [0.0_dp, -1.0_dp]
 hi = [2.0_dp, 1.5_dp]
 nu = 5.0_dp

 call set_seed_int(303)
 x = rtmvnorm(250, mu, sig, lo, hi)
 call check_bounds('rtmvnorm', x, lo, hi)

 call set_seed_int(304)
 xt = rtmvt(250, mu, sig, nu, lo, hi)
 call check_bounds('rtmvt', xt, lo, hi)

 call set_seed_int(305)
 tr = tregress(180, lo, hi, sig, nu)
 if (tr%status /= 0) then
  print '(a,i0)', 'tregress solver status: ', tr%status
  error stop 1
 end if
 allocate(reconstructed(size(tr%z,1), size(tr%z,2)))
 do i = 1, size(tr%z,1)
  reconstructed(i,:) = sqrt(nu)*tr%z(i,:)/tr%r(i)
 end do
 call check_bounds('tregress reconstruction', reconstructed, lo, hi)
 if (any(tr%r <= 0.0_dp)) error stop 'tregress returned nonpositive scale'

 ! The R wrappers condition analytically when a coordinate is degenerate.
 lo = [0.5_dp, -1.0_dp]
 hi = [0.5_dp, 1.5_dp]
 call set_seed_int(306)
 x = rtmvnorm(120, mu, sig, lo, hi)
 if (maxval(abs(x(:,1)-0.5_dp)) > 1.0e-14_dp) error stop 'degenerate rtmvnorm coordinate failed'
 call check_bounds('degenerate rtmvnorm', x, lo, hi)
 call set_seed_int(307)
 xt = rtmvt(120, mu, sig, nu, lo, hi)
 if (maxval(abs(xt(:,1)-0.5_dp)) > 1.0e-14_dp) error stop 'degenerate rtmvt coordinate failed'
 call check_bounds('degenerate rtmvt', xt, lo, hi)

 print '(a)', 'test_sampling: PASS'
contains
 subroutine check_bounds(name, y, lower, upper)
  character(len=*), intent(in) :: name
  real(dp), intent(in) :: y(:,:), lower(:), upper(:)
  integer :: k
  do k = 1, size(y,1)
   if (any(y(k,:) < lower-1.0e-12_dp) .or. any(y(k,:) > upper+1.0e-12_dp)) then
    print '(a,i0)', trim(name)//' bound failure at row ', k
    error stop 1
   end if
  end do
 end subroutine check_bounds
end program test_sampling
