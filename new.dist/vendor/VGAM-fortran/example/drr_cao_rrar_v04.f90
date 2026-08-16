program drr_cao_rrar_v04
   use vgam
   implicit none
   integer, parameter :: n = 36
   real(dp) :: x(n, 3), yd(n, 2), h_a(2, 1, 1), h_c(1, 1, 2)
   real(dp) :: xenv(n, 2), ya(n, 2), z, t, u
   real(dp) :: yts(80, 2), phi0(2, 2)
   integer :: i
   integer :: fam2(2)
   logical :: nr(3)
   type(drrvglm_result_t) :: dfit
   type(cao_result_t) :: afit
   type(rrar_result_t) :: rfit

   fam2 = family_gaussian
   nr = [.true., .false., .false.]
   h_a(:, 1, 1) = [1.0_dp, 2.0_dp]
   h_c = 1.0_dp
   do i = 1, n
      t = -1.0_dp + 2.0_dp*real(i - 1, dp)/real(n - 1, dp)
      x(i, :) = [1.0_dp, t, t*t - 0.3_dp]
      z = 0.8_dp*x(i, 2) - 0.35_dp*x(i, 3)
      yd(i, :) = [0.25_dp + z, -0.4_dp + 2.0_dp*z]
      u = sin(0.7_dp*real(i, dp))
      xenv(i, :) = [t, u]
      z = t + 0.55_dp*u
      ya(i, :) = [sin(1.7_dp*z) + 0.15_dp*z, cos(1.2_dp*z) - 0.1_dp*z*z]
   end do
   call fit_drrvglm(yd, x, 1, fam2, h_a, h_c, dfit, no_rrr=nr)
   print '(a,f12.6)', 'DRR residual deviance: ', dfit%deviance
   print '(a,2f10.4)', 'DRR constrained loadings: ', dfit%loadings(:, 1)

   call fit_cao_rank1(ya, xenv, fam2, afit, df=6, lambda=0.05_dp, max_iter=8)
   print '(a,f12.6)', 'CAO residual deviance: ', afit%deviance
   print '(a,2f10.4)', 'CAO canonical coefficients: ', afit%canonical_coefficients

   phi0 = reshape([0.50_dp, 0.25_dp, 0.10_dp, 0.05_dp], [2, 2])
   yts = 0.0_dp
   yts(1, :) = [0.6_dp, -0.2_dp]
   do i = 2, size(yts, 1)
      yts(i, :) = matmul(phi0, yts(i - 1, :)) + &
         [0.025_dp*sin(0.71_dp*real(i, dp)), 0.020_dp*cos(0.53_dp*real(i, dp))]
   end do
   call fit_rrar(yts, [1], rfit, compute_covariance=.false.)
   print '(a)', 'RRAR fitted lag-1 matrix:'
   do i = 1, 2
      print '(2f10.4)', rfit%phi(i, :, 1)
   end do
end program drr_cao_rrar_v04
