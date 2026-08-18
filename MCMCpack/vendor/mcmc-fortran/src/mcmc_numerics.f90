module mcmc_numerics
   use mcmc_kinds, only : dp
   implicit none
   private
   real(dp), parameter :: pi = acos(-1.0_dp)

   public :: set_mcmc_seed, rand_uniform, rand_normal_vec, euclid_norm

contains

   subroutine set_mcmc_seed(seed)
      integer, intent(in) :: seed
      integer :: n, i
      integer, allocatable :: state(:)
      call random_seed(size=n)
      allocate(state(n))
      do i = 1, n
         state(i) = mod(seed + 104729*i, 2147483646) + 1
      end do
      call random_seed(put=state)
   end subroutine set_mcmc_seed

   real(dp) function rand_uniform() result(u)
      call random_number(u)
      if (u <= 0.0_dp) u = tiny(1.0_dp)
      if (u >= 1.0_dp) u = 1.0_dp - spacing(1.0_dp)
   end function rand_uniform

   subroutine rand_normal_vec(z)
      real(dp), intent(out) :: z(:)
      integer :: i
      real(dp) :: u1, u2, r, theta
      i = 1
      do while (i <= size(z))
         u1 = rand_uniform()
         u2 = rand_uniform()
         r = sqrt(-2.0_dp*log(u1))
         theta = 2.0_dp*pi*u2
         z(i) = r*cos(theta)
         if (i+1 <= size(z)) z(i+1) = r*sin(theta)
         i = i + 2
      end do
   end subroutine rand_normal_vec

   pure real(dp) function euclid_norm(x) result(v)
      real(dp), intent(in) :: x(:)
      v = sqrt(sum(x*x))
   end function euclid_norm

end module mcmc_numerics
