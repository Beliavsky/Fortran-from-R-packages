program causal_example
   use grf, only : dp, grf_options, grf_forest, grf_ate_result
   use grf, only : causal_forest, predict_causal_forest, average_treatment_effect
   implicit none

   integer, parameter :: n = 100
   real(dp) :: x(n,2)
   real(dp) :: y(n)
   real(dp) :: w(n)
   real(dp), allocatable :: tau_hat(:)
   type(grf_options) :: options
   type(grf_forest) :: forest
   type(grf_ate_result) :: ate
   integer :: i
   integer :: info

   do i = 1, n
      x(i,1) = -1.0_dp + 2.0_dp * real(i - 1,dp) / real(n - 1,dp)
      x(i,2) = cos(0.23_dp * real(i,dp))
      if (mod(i,2) == 0) then
         w(i) = 1.0_dp
      else
         w(i) = 0.0_dp
      end if
      y(i) = 0.5_dp * x(i,1) + (1.0_dp + x(i,1)) * w(i)
   end do
   options%num_trees = 50
   options%min_node_size = 4
   options%sample_fraction = 0.8_dp
   options%seed = 777
   call causal_forest(x, y, w, forest, info, options)
   if (info /= 0) error stop 'causal_forest failed'
   call predict_causal_forest(forest, x(1:5,:), tau_hat)
   call average_treatment_effect(forest, ate, info)
   if (info /= 0) error stop 'average_treatment_effect failed'
   print '(a,f10.5,a,f10.5)', 'ATE = ', ate%estimate(1), '  SE = ', ate%std_err(1)
   print '(a,5f10.5)', 'First five CATE estimates: ', tau_hat
end program causal_example
