module new_dist_numerics
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use new_dist_kinds, only : dp
   implicit none
   private
   real(dp), parameter, public :: pi = acos(-1.0_dp)
   public :: target_probability, finish_probability, rand_uniform, log1p_stable
   public :: set_new_dist_seed, sign0, nan_dp

contains

   pure elemental real(dp) function log1p_stable(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: x2
      if (abs(x) < 1.0e-5_dp) then
         x2=x*x
         y=x-0.5_dp*x2+x*x2/3.0_dp-x2*x2/4.0_dp+x2*x2*x/5.0_dp
      else
         y=log(1.0_dp+x)
      end if
   end function log1p_stable

   pure elemental real(dp) function nan_dp() result(x)
      x = ieee_value(0.0_dp,ieee_quiet_nan)
   end function nan_dp

   pure elemental real(dp) function sign0(x) result(s)
      real(dp), intent(in) :: x
      if (x > 0.0_dp) then
         s = 1.0_dp
      else if (x < 0.0_dp) then
         s = -1.0_dp
      else
         s = 0.0_dp
      end if
   end function sign0

   pure elemental real(dp) function target_probability(p,lower_tail) result(t)
      real(dp), intent(in) :: p
      logical, intent(in), optional :: lower_tail
      logical :: lt
      lt = .true.
      if (present(lower_tail)) lt = lower_tail
      if (lt) then
         t = p
      else
         t = 1.0_dp-p
      end if
   end function target_probability

   pure elemental real(dp) function finish_probability(cdf,lower_tail,log_p) result(v)
      real(dp), intent(in) :: cdf
      logical, intent(in), optional :: lower_tail,log_p
      logical :: lt,lp
      real(dp) :: p
      lt = .true.; lp = .false.
      if (present(lower_tail)) lt = lower_tail
      if (present(log_p)) lp = log_p
      if (lt) then
         p = min(1.0_dp,max(0.0_dp,cdf))
      else
         p = min(1.0_dp,max(0.0_dp,1.0_dp-cdf))
      end if
      if (lp) then
         if (p <= 0.0_dp) then
            v = -huge(1.0_dp)
         else
            v = log(p)
         end if
      else
         v = p
      end if
   end function finish_probability

   real(dp) function rand_uniform() result(u)
      call random_number(u)
      if (u <= 0.0_dp) u = tiny(1.0_dp)
      if (u >= 1.0_dp) u = 1.0_dp-spacing(1.0_dp)
   end function rand_uniform

   subroutine set_new_dist_seed(seed)
      integer, intent(in) :: seed
      integer :: n,i
      integer, allocatable :: state(:)
      call random_seed(size=n)
      allocate(state(n))
      do i=1,n
         state(i)=mod(seed+104729*i,2147483646)+1
      end do
      call random_seed(put=state)
   end subroutine set_new_dist_seed


end module new_dist_numerics
