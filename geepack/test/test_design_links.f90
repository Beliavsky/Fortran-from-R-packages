program test_design_links
   use geepack
   implicit none
   integer :: cs(2)
   integer :: waves(5)
   real(dp), allocatable :: z(:, :)
   real(dp), allocatable :: zo(:, :)
   real(dp) :: r(3, 3)
   real(dp) :: d(3, 1)
   real(dp) :: rho(1)
   real(dp) :: ru(2, 2)
   real(dp) :: rhou(6)
   integer :: status

   cs = [2, 3]
   waves = [1, 2, 1, 2, 3]
   call gen_zcor(cs, waves, COR_UNSTRUCTURED, z, status)
   if (status /= GEE_OK .or. size(z, 1) /= 4 .or. size(z, 2) /= 3) error stop 1
   if (maxval(abs(z(1, :) - [1.0_dp, 0.0_dp, 0.0_dp])) > 1.0e-14_dp) error stop 2
   if (maxval(abs(z(4, :) - [0.0_dp, 0.0_dp, 1.0_dp])) > 1.0e-14_dp) error stop 3
   call gen_zodds(cs, waves, COR_UNSTRUCTURED, 2, zo, status)
   if (status /= GEE_OK .or. size(zo, 1) /= 16 .or. size(zo, 2) /= 3) error stop 4
   rho = 0.4_dp
   call working_correlation(rho, [1, 2, 3], COR_AR1, r, status)
   if (status /= GEE_OK) error stop 5
   if (abs(r(1, 3) - 0.16_dp) > 1.0e-14_dp) error stop 6
   call correlation_rho_derivative(rho, [1, 2, 3], COR_AR1, d, status)
   if (status /= GEE_OK) error stop 7
   if (maxval(abs(d(:, 1) - [1.0_dp, 0.8_dp, 1.0_dp])) > 1.0e-14_dp) error stop 8
   rhou = [0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp, 0.5_dp, 0.6_dp]
   call working_correlation(rhou, [2, 3], COR_UNSTRUCTURED, ru, status)
   if (status /= GEE_OK .or. abs(ru(1, 2) - 0.4_dp) > 1.0e-14_dp) error stop 9
   if (abs(link_inverse(0.0_dp, LINK_LOGIT) - 0.5_dp) > 1.0e-14_dp) error stop 10
   if (abs(variance_function(0.5_dp, VAR_BINOMIAL) - 0.25_dp) > 1.0e-14_dp) error stop 11
   print *, 'test_design_links: PASS'
end program test_design_links
