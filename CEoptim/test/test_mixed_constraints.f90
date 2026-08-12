program test_mixed_constraints
   use ceoptim
   implicit none
   type(ce_control) :: ctl
   type(ce_continuous_control) :: cc
   type(ce_discrete_control) :: dc
   type(ce_result) :: res

   allocate(cc%mean(2), cc%sd(2), cc%con_mat(3,2), cc%con_vec(3))
   cc%mean = [0.0_dp, 0.0_dp]
   cc%sd = [2.0_dp, 2.0_dp]
   cc%con_mat = reshape([ -1.0_dp, 0.0_dp, 1.0_dp, &
                           0.0_dp,-1.0_dp, 1.0_dp ], [3,2])
   cc%con_vec = [0.0_dp, 0.0_dp, 1.2_dp]
   cc%smooth_mean = 0.9_dp
   cc%smooth_sd = 0.8_dp
   cc%sd_thr = 2.0e-4_dp

   allocate(dc%categories(1))
   dc%categories = [4]
   dc%smooth_prob = 0.8_dp
   dc%prob_thr = 1.0e-4_dp

   ctl%n = 350
   ctl%rho = 0.1_dp
   ctl%iter_thr = 120
   ctl%no_improve_thr = 10
   ctl%seed = 777_i64

   call ce_optimize(obj, res, ctl, continuous=cc, discrete=dc)
   if (res%status /= 0) error stop res%message
   if (res%optimum > 2.0e-3_dp) error stop 'mixed constrained optimum too high'
   if (abs(res%continuous(1)-0.4_dp) > 0.08_dp) error stop 'mixed x1 inaccurate'
   if (abs(res%continuous(2)-0.6_dp) > 0.08_dp) error stop 'mixed x2 inaccurate'
   if (res%discrete(1) /= 2) error stop 'mixed discrete optimizer incorrect'
   if (res%continuous(1) < -1.0e-10_dp .or. res%continuous(2) < -1.0e-10_dp) &
      error stop 'constraint violation'
   if (sum(res%continuous) > 1.2_dp + 1.0e-10_dp) error stop 'constraint violation'
   print *, 'test_mixed_constraints: PASS', res%optimum, res%continuous, res%discrete

contains

   function obj(xc, xd) result(v)
      real(dp), intent(in) :: xc(:)
      integer, intent(in) :: xd(:)
      real(dp) :: v
      v = (xc(1)-0.4_dp)**2 + (xc(2)-0.6_dp)**2 + real((xd(1)-2)**2, dp)
   end function obj
end program test_mixed_constraints
