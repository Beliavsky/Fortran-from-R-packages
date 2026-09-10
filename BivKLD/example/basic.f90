program basic
   use bivkld, only : BIVKLD_SUCCESS, biv_kld, biv_kld_normal, biv_kld_pareto2, dp
   implicit none

   real(dp) :: x(6, 2), y(6, 2), H(2, 2), estimate
   integer :: info

   x(:, 1) = [-1.2_dp, -0.7_dp, -0.1_dp, 0.4_dp, 0.9_dp, 1.3_dp]
   x(:, 2) = [-0.4_dp, 0.5_dp, -0.8_dp, 0.9_dp, 0.1_dp, 1.1_dp]
   y(:, 1) = x(:, 1) + 0.30_dp
   y(:, 2) = x(:, 2) - 0.20_dp
   H = reshape([0.35_dp, 0.0_dp, 0.0_dp, 0.45_dp], [2, 2])

   call biv_kld(x, y, estimate, Hx=H, Hy=H, grid_size=[30, 30], info=info)
   if (info /= BIVKLD_SUCCESS) error stop 'biv_kld failed'

   print '(a,f12.6)', 'Kernel D(x || y) = ', estimate
   print '(a,f12.6)', 'Pareto-II D(1 || 2) = ', biv_kld_pareto2(1.0_dp, 2.0_dp)
   print '(a,f12.6)', 'Normal identity D = ', &
      biv_kld_normal([0.0_dp, 0.0_dp], reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2]), &
                     [0.0_dp, 0.0_dp], reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2]))
end program basic
