program test_rng
   use bridgedist, only : dp, rbridge, bridge_variance, pbridge
   implicit none
   integer, parameter :: n = 120000
   real(dp), allocatable :: x(:)
   real(dp) :: meanx, varx, phi, pit_mean
   integer :: fails

   fails = 0
   phi = 1.0_dp / sqrt(1.0_dp + 3.0_dp / acos(-1.0_dp)**2)
   allocate(x(n))
   call seed_rng()
   call rbridge(x, phi)
   meanx = sum(x) / real(n, dp)
   varx = sum((x - meanx)**2) / real(n - 1, dp)
   pit_mean = sum(pbridge(x, phi)) / real(n, dp)

   if (abs(meanx) > 0.025_dp) then
      print '(a,es14.6)', 'mean FAIL ', meanx
      fails = fails + 1
   end if
   if (abs(varx - bridge_variance(phi)) > 0.04_dp) then
      print '(a,2es14.6)', 'variance FAIL ', varx, bridge_variance(phi)
      fails = fails + 1
   end if
   if (abs(pit_mean - 0.5_dp) > 0.005_dp) then
      print '(a,es14.6)', 'PIT mean FAIL ', pit_mean
      fails = fails + 1
   end if

   if (fails /= 0) error stop 1
   print '(a)', 'test_rng: PASS'

contains

   subroutine seed_rng()
      integer :: nseed, i
      integer, allocatable :: seed(:)
      call random_seed(size=nseed)
      allocate(seed(nseed))
      do i = 1, nseed
         seed(i) = 314159 + 2718 * i
      end do
      call random_seed(put=seed)
   end subroutine seed_rng

end program test_rng
