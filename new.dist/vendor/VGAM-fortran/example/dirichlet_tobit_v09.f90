program dirichlet_tobit_v09
   use vgam
   implicit none
   integer, parameter :: nd = 500, nt = 650
   real(dp) :: yd(nd, 3), xd(nd, 1), alpha(3), draw(3)
   real(dp) :: yt(nt), xt(nt, 2), z, mu
   type(dirichlet_regression_result_t) :: dfit
   type(tobit_result_t) :: tfit
   integer :: i, stat, nseed
   integer, allocatable :: seed(:)

   call random_seed(size=nseed)
   allocate(seed(nseed))
   seed = 909900
   call random_seed(put=seed)

   alpha = [2.0_dp, 3.0_dp, 5.0_dp]
   xd = 1.0_dp
   do i = 1, nd
      call random_dirichlet(alpha, draw, stat)
      yd(i, :) = draw
   end do
   call fit_dirichlet_regression(yd, xd, dfit, max_iter=300, tol=2.0e-6_dp)

   do i = 1, nt
      call random_number(z)
      z = 2.0_dp*z - 1.0_dp
      xt(i, :) = [1.0_dp, z]
      mu = 0.35_dp + 0.75_dp*z
      yt(i) = rtobit_v(mu, 0.85_dp, 0.0_dp, 1.5_dp)
   end do
   call fit_tobit(yt, xt, tfit, 0.0_dp, 1.5_dp, max_iter=350, tol=2.0e-6_dp)

   print '(a,3f10.4)', 'Dirichlet fitted shapes:', dfit%fitted_shape(1, :)
   print '(a,2f10.4)', 'Tobit mean coefficients:', tfit%mean_coefficients
   print '(a,f10.4)', 'Tobit fitted sd:', exp(tfit%scale_coefficients(1))
end program dirichlet_tobit_v09
