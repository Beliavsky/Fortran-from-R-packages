program mixed_example
   use ceoptim
   implicit none
   type(ce_control) :: ctl
   type(ce_continuous_control) :: continuous
   type(ce_discrete_control) :: discrete
   type(ce_result) :: res

   continuous%mean = [0.0_dp, 0.0_dp]
   continuous%sd = [2.0_dp, 2.0_dp]
   continuous%con_mat = reshape([-1.0_dp,0.0_dp,1.0_dp, 0.0_dp,-1.0_dp,1.0_dp], [3,2])
   continuous%con_vec = [0.0_dp, 0.0_dp, 1.2_dp]
   discrete%categories = [4]
   ctl%n = 250
   ctl%rho = 0.1_dp
   ctl%seed = 777_i64

   call ce_optimize(objective, res, ctl, continuous, discrete)
   if (res%status /= 0) error stop res%message
   print '(a,es14.6)', 'minimum = ', res%optimum
   print '(a,2f10.5)', 'continuous optimizer = ', res%continuous
   print '(a,i0)', 'discrete optimizer = ', res%discrete(1)

contains
   function objective(xc, xd) result(v)
      real(dp), intent(in) :: xc(:)
      integer, intent(in) :: xd(:)
      real(dp) :: v
      v = (xc(1)-0.4_dp)**2 + (xc(2)-0.6_dp)**2 + real((xd(1)-2)**2, dp)
   end function objective
end program mixed_example
