program test_discrete
   use ceoptim
   implicit none
   type(ce_control) :: ctl
   type(ce_discrete_control) :: dc
   type(ce_result) :: res

   allocate(dc%categories(3))
   dc%categories = [5, 4, 3]
   dc%smooth_prob = 0.8_dp
   dc%prob_thr = 1.0e-4_dp

   ctl%n = 240
   ctl%rho = 0.1_dp
   ctl%iter_thr = 100
   ctl%no_improve_thr = 8
   ctl%seed = 81723_i64

   call ce_optimize(obj, res, ctl, discrete=dc)
   if (res%status /= 0) error stop res%message
   if (abs(res%optimum) > tiny(1.0_dp)) error stop 'discrete optimum is not exact'
   if (any(res%discrete /= [2, 1, 0])) error stop 'discrete optimizer incorrect'
   if (size(res%states) < 1) error stop 'state history missing'
   print *, 'test_discrete: PASS', res%optimum, res%discrete

contains

   function obj(xc, xd) result(v)
      real(dp), intent(in) :: xc(:)
      integer, intent(in) :: xd(:)
      real(dp) :: v
      if (size(xc) /= 0) error stop 'unexpected continuous argument'
      v = real((xd(1)-2)**2 + 2*(xd(2)-1)**2 + 3*xd(3)**2, dp)
   end function obj
end program test_discrete
