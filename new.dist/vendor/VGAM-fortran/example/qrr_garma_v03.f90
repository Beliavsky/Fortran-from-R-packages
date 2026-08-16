program qrr_garma_v03
   use vgam
   implicit none
   integer, parameter :: n = 48
   real(dp) :: x(n,3), y(n,2), z(n), t
   real(dp) :: xt(n,2), yt(n)
   real(dp), allocatable :: optimum(:,:), curvature(:,:,:)
   integer :: fam(2), i
   logical :: no_rrr(3)
   type(qrrvglm_result_t) :: qfit
   type(garma_result_t) :: gfit

   fam = [family_gaussian, family_gaussian]
   no_rrr = [.true., .false., .false.]
   do i = 1, n
      t = -1.0_dp + 2.0_dp*real(i-1,dp)/real(n-1,dp)
      x(i,:) = [1.0_dp, t, sin(0.4_dp*real(i,dp))]
      z(i) = 0.9_dp*x(i,2) - 0.35_dp*x(i,3)
      y(i,1) = 1.0_dp + 0.6_dp*z(i) - 0.55_dp*z(i)**2
      y(i,2) = -0.3_dp - 0.2_dp*z(i) - 0.80_dp*z(i)**2
   end do

   call fit_qrrvglm(y, x, 1, fam, qfit, no_rrr=no_rrr)
   call qfit%optima(optimum, curvature)
   print '(a,l1)', 'QRR-VGLM converged: ', qfit%converged
   print '(a,f12.6)', 'QRR residual deviance: ', qfit%deviance
   print '(a,2f10.4)', 'Latent optima: ', optimum(:,1)

   do i = 1, n
      t = real(i-1,dp)/real(n-1,dp)
      xt(i,:) = [1.0_dp, t]
   end do
   yt(1) = 0.8_dp
   do i = 2, n
      yt(i) = 0.2_dp + 0.7_dp*xt(i,2) + 0.5_dp * &
              (yt(i-1) - (0.2_dp + 0.7_dp*xt(i-1,2)))
   end do
   call fit_garma(yt, xt, 1, link_identity, gfit)
   print '(a,2f10.4)', 'GARMA regression coefficients: ', gfit%coefficients
   print '(a,f10.4)', 'GARMA AR(1) coefficient: ', gfit%ar(1)
end program qrr_garma_v03
