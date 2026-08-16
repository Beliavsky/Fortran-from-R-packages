program test_v03
   use vgam
   implicit none
   integer, parameter :: n = 64
   real(dp) :: x(n,3), yq(n,2), yq2(n,3), z(n), u, v
   real(dp) :: yg(n), xg(n,2), beta_true(2), phi_true
   real(dp) :: yl(n), xl(n,1), xm(n,2), xs(n,1), t, noise
   real(dp), allocatable :: pred(:,:), pred2(:,:), optimum(:,:), curvature(:,:,:)
   real(dp), allocatable :: fc(:), qq25(:), qq75(:)
   integer :: fam(2), fam3(3), i
   logical :: nr(3)
   type(qrrvglm_result_t) :: qfit, qfit2
   type(garma_result_t) :: gfit
   type(lms_yj_result_t) :: lfit

   fam = [family_gaussian, family_gaussian]
   fam3 = [family_gaussian, family_gaussian, family_gaussian]
   nr = [.true., .false., .false.]
   do i = 1, n
      u = -1.5_dp + 3.0_dp*real(i-1,dp)/real(n-1,dp)
      v = sin(0.31_dp*real(i,dp))
      x(i,:) = [1.0_dp, u, v]
      z(i) = 0.8_dp*u - 0.45_dp*v
      yq(i,1) = 1.2_dp + 0.7_dp*z(i) - 0.55_dp*z(i)**2
      yq(i,2) = -0.4_dp - 0.25_dp*z(i) - 0.85_dp*z(i)**2
   end do
   call fit_qrrvglm(yq, x, 1, fam, qfit, no_rrr=nr, max_iter=40, tol=1.0e-7_dp)
   call qfit%predict(x, pred)
   call qfit%optima(optimum, curvature)
   call assert_true(allocated(qfit%quadratic), 'QRR allocation')
   call assert_true(maxval(abs(pred-yq)) < 2.0e-3_dp, 'QRR exact quadratic fit')
   call assert_true(all(curvature(:,1,1) < 0.0_dp), 'QRR negative curvature')
   call assert_true(all(abs(optimum(:,1)) < 5.0_dp), 'QRR finite optima')

   do i = 1, n
      u = x(i,2)
      v = x(i,3)
      yq2(i,1) = 0.3_dp + 0.7_dp*u - 0.2_dp*v - 0.4_dp*u*u + &
                  0.35_dp*u*v - 0.6_dp*v*v
      yq2(i,2) = -0.2_dp - 0.1_dp*u + 0.8_dp*v - 0.2_dp*u*u - &
                  0.25_dp*u*v - 0.45_dp*v*v
      yq2(i,3) = 0.8_dp + 0.4_dp*u + 0.5_dp*v - 0.7_dp*u*u + &
                  0.15_dp*u*v - 0.3_dp*v*v
   end do
   call fit_qrrvglm(yq2, x, 2, fam3, qfit2, no_rrr=nr, max_iter=15, tol=1.0e-8_dp)
   call qfit2%predict(x, pred2)
   call assert_true(maxval(abs(pred2-yq2)) < 2.0e-5_dp, 'QRR full cross-quadratic fit')

   beta_true = [0.35_dp, 0.8_dp]
   phi_true = 0.55_dp
   do i = 1, n
      t = real(i-1,dp)/real(n-1,dp)
      xg(i,:) = [1.0_dp, t]
   end do
   yg(1) = dot_product(xg(1,:), beta_true) + 0.45_dp
   do i = 2, n
      yg(i) = dot_product(xg(i,:), beta_true) + phi_true * &
              (yg(i-1) - dot_product(xg(i-1,:), beta_true))
   end do
   call fit_garma(yg, xg, 1, link_identity, gfit, max_iter=250, tol=1.0e-9_dp)
   call assert_true(maxval(abs(gfit%coefficients-beta_true)) < 2.0e-3_dp, 'GARMA beta')
   call assert_true(abs(gfit%ar(1)-phi_true) < 2.0e-3_dp, 'GARMA AR')
   call gfit%forecast(xg(n:n,:), yg(n:n), xg(n:n,:), fc)
   call assert_true(size(fc) == 1 .and. abs(fc(1)) < huge(1.0_dp), 'GARMA forecast')

   do i = 1, n
      t = -1.0_dp + 2.0_dp*real(i-1,dp)/real(n-1,dp)
      xl(i,1) = 1.0_dp
      xm(i,:) = [1.0_dp, t]
      xs(i,1) = 1.0_dp
      noise = 0.30_dp*sin(1.7_dp*real(i,dp))
      yl(i) = 0.25_dp + 0.65_dp*t + noise
   end do
   call fit_lms_yj(yl, xl, xm, xs, lfit, max_iter=350, tol=2.0e-7_dp)
   call lfit%predict_quantile(xl, xm, xs, 0.25_dp, qq25)
   call lfit%predict_quantile(xl, xm, xs, 0.75_dp, qq75)
   call assert_true(size(qq25) == n, 'LMS quantile size')
   call assert_true(all(qq75 > qq25), 'LMS quantile ordering')
   call assert_true(all(lfit%sigma > 0.0_dp), 'LMS positive sigma')

   print '(a)', 'test_v03: PASS'
contains
   subroutine assert_true(ok, label)
      logical, intent(in) :: ok
      character(*), intent(in) :: label
      if (.not. ok) then
         print '(a)', 'FAIL: '//trim(label)
         error stop 1
      end if
   end subroutine assert_true
end program test_v03
