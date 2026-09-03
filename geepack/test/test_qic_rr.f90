program test_qic_rr
   use geepack
   implicit none
   real(dp) :: yq(4)
   real(dp) :: mu(4)
   real(dp) :: vr(2, 2)
   real(dp) :: vi(2, 2)
   type(qic_result) :: qic
   integer :: status
   real(dp) :: y(8)
   real(dp) :: x(8, 2)
   integer :: cs(4)
   real(dp), allocatable :: yc(:)
   real(dp), allocatable :: xc(:, :)
   real(dp), allocatable :: wc(:)
   integer, allocatable :: cc(:)
   type(gee_result) :: fit
   real(dp) :: b0(2)
   integer :: i
   real(dp) :: se(2)
   real(dp) :: wald(2)
   real(dp) :: pv(2)
   real(dp) :: stat
   real(dp) :: pval
   integer :: df
   real(dp) :: beta_test(2)
   real(dp) :: contrast(1, 2)
   real(dp) :: rhs(1)

   yq = [0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp]
   mu = [0.2_dp, 0.7_dp, 0.8_dp, 0.3_dp]
   vr = reshape([0.04_dp, 0.01_dp, 0.01_dp, 0.09_dp], [2, 2])
   vi = reshape([0.05_dp, 0.005_dp, 0.005_dp, 0.08_dp], [2, 2])
   call compute_qic(yq, mu, VAR_BINOMIAL, vr, vi, 1, 10, qic, status)
   if (status /= GEE_OK) error stop 1
   if (abs(qic%cic - 1.9119496855345912_dp) > 1.0e-12_dp) error stop 2
   if (abs(qic%qic - 6.143173352080951_dp) > 1.0e-12_dp) error stop 3
   if (abs(qic%qicc - 10.14317335208095_dp) > 1.0e-12_dp) error stop 4

   cs = 2
   y = [0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp]
   x(:, 1) = 1.0_dp
   do i = 1, 8
      x(i, 2) = merge(1.0_dp, 0.0_dp, mod(i, 2) == 0)
   end do
   call make_relative_risk_copy(y, x, cs, 100, yc, xc, cc, wc, status)
   if (status /= GEE_OK) error stop 5
   if (size(yc) /= 16 .or. size(cc) /= 8) error stop 6
   if (maxval(abs(yc(9:16) - (1.0_dp - y))) > 0.0_dp) error stop 7
   if (abs(sum(wc(1:8)) - 7.92_dp) > 1.0e-12_dp) error stop 8
   if (abs(sum(wc(9:16)) - 0.08_dp) > 1.0e-12_dp) error stop 9
   b0 = [log(0.35_dp), 0.2_dp]
   call fit_relative_risk(y, x, cs, 100, COR_INDEPENDENCE, fit, beta_initial=b0)
   if (fit%error /= GEE_OK) error stop 10
   if (any(exp(matmul(x, fit%beta)) <= 0.0_dp)) error stop 11
   beta_test = [1.0_dp, 2.0_dp]
   call coefficient_wald_summary(beta_test, reshape([0.25_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2]), &
      se, wald, pv, status)
   if (status /= GEE_OK .or. maxval(abs(wald - 4.0_dp)) > 1.0e-12_dp) error stop 12
   if (maxval(abs(pv - [0.0455002638963584_dp, 0.0455002638963584_dp])) > 2.0e-12_dp) error stop 13
   contrast = reshape([0.0_dp, 1.0_dp], [1, 2])
   rhs = 0.0_dp
   call wald_contrast(beta_test, reshape([0.25_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2]), &
      contrast, rhs, stat, df, pval, status)
   if (status /= GEE_OK .or. abs(stat - 4.0_dp) > 1.0e-12_dp .or. df /= 1) error stop 14
   print *, 'test_qic_rr: PASS'
end program test_qic_rr
