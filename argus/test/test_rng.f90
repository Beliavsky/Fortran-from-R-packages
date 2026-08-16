program test_rng
   use argus, only : dp, pargus, rargus, rargus_varying, seed_argus_rng, &
      ARGUS_INVERSION, ARGUS_ROU
   implicit none
   integer, parameter :: n = 20000
   real(dp), allocatable :: x(:), u(:), chi(:)
   real(dp) :: mean_u, var_u
   integer :: i

   allocate(x(n),u(n),chi(n))

   call seed_argus_rng(12345)
   call rargus(x,0.3_dp,ARGUS_INVERSION)
   if (any(x < 0.0_dp) .or. any(x > 1.0_dp)) error stop "inversion RNG outside support"
   u = pargus(x,0.3_dp)
   call check_uniform(u,"inversion RNG")

   call seed_argus_rng(24680)
   call rargus(x,3.0_dp,ARGUS_ROU)
   if (any(x < 0.0_dp) .or. any(x > 1.0_dp)) error stop "RoU RNG outside support"
   u = pargus(x,3.0_dp)
   call check_uniform(u,"RoU RNG")

   do i = 1, n
      chi(i) = 0.02_dp + 5.98_dp*real(i-1,dp)/real(n-1,dp)
   end do
   call seed_argus_rng(98765)
   call rargus_varying(x,chi,ARGUS_INVERSION)
   do i = 1, n
      u(i) = pargus(x(i),chi(i))
   end do
   mean_u = sum(u)/real(n,dp)
   var_u = sum((u-mean_u)**2)/real(n-1,dp)
   if (abs(mean_u-0.5_dp) > 0.012_dp) error stop "varying-chi RNG mean check failed"
   if (abs(var_u-1.0_dp/12.0_dp) > 0.008_dp) error stop "varying-chi RNG variance check failed"

   print '(a)', "test_rng: PASS"

contains

   subroutine check_uniform(v,label)
      real(dp), intent(in) :: v(:)
      character(*), intent(in) :: label
      real(dp) :: m, vv

      m = sum(v)/real(size(v),dp)
      vv = sum((v-m)**2)/real(size(v)-1,dp)
      if (abs(m-0.5_dp) > 0.012_dp .or. abs(vv-1.0_dp/12.0_dp) > 0.008_dp) then
         print '(a)', "FAIL: "//label
         print '(a,f12.7)', "mean = ", m
         print '(a,f12.7)', "var  = ", vv
         error stop 1
      end if
   end subroutine check_uniform

end program test_rng
