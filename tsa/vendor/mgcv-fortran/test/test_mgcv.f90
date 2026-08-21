program test_mgcv
   use mgcv
   implicit none
   integer :: failures
   failures = 0
   call test_transforms(failures)
   call test_linalg(failures)
   call test_smooths(failures)
   call test_gaussian_gam(failures)
   call test_magic_reference(failures)
   call test_glm_families(failures)
   call test_constraints(failures)
   call test_distributions(failures)
   call test_simulation(failures)
   if (failures /= 0) then
      write(*,'(a,i0)') 'FAILURES: ', failures
      error stop 1
   end if
   write(*,'(a)') 'All mgcv-fortran tests passed.'
contains

   subroutine assert_true(condition, label, failures)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      integer, intent(inout) :: failures
      if (.not. condition) then
         failures = failures + 1
         write(*,'(a,a)') 'FAIL: ', trim(label)
      end if
   end subroutine assert_true

   subroutine assert_close(value, expected, tolerance, label, failures)
      real(dp), intent(in) :: value, expected, tolerance
      character(len=*), intent(in) :: label
      integer, intent(inout) :: failures
      call assert_true(abs(value - expected) <= tolerance, label, failures)
      if (abs(value - expected) > tolerance) then
         write(*,'(a,3es16.7)') '  value expected tolerance: ', value, expected, tolerance
      end if
   end subroutine assert_close

   subroutine test_transforms(failures)
      integer, intent(inout) :: failures
      real(dp) :: x
      integer :: i
      do i = -10, 10
         x = 0.5_dp * real(i, dp)
         call assert_close(not_log(not_exp(x)), x, 2.0e-12_dp, 'notLog(notExp(x))', failures)
         call assert_close(not_log2(not_exp2(x)), x, 2.0e-12_dp, 'notLog2(notExp2(x))', failures)
      end do
      call assert_true(null_space_dimension(2, 2) == 3, 'null space dimension', failures)
   end subroutine test_transforms

   subroutine test_linalg(failures)
      integer, intent(inout) :: failures
      real(dp) :: a(3,3), ld(3), sd(2)
      real(dp), allocatable :: inv(:, :), rd(:), ro(:), root(:, :)
      integer :: status, rank
      a = reshape([4.0_dp,1.0_dp,0.5_dp, 1.0_dp,3.0_dp,0.2_dp, 0.5_dp,0.2_dp,2.0_dp],[3,3])
      call spd_inverse(a, inv, status)
      call assert_true(status == 0, 'spd inverse status', failures)
      call assert_true(maxval(abs(matmul(a, inv) - &
         reshape([1.0_dp,0.0_dp,0.0_dp,0.0_dp,1.0_dp,0.0_dp, &
                  0.0_dp,0.0_dp,1.0_dp],[3,3]))) < 1.0e-10_dp, &
         'spd inverse', failures)
      call mroot(a, root, rank, status)
      call assert_true(status == 0 .and. rank == 3, 'mroot rank', failures)
      call assert_true(maxval(abs(matmul(root, transpose(root)) - a)) < 1.0e-9_dp, 'mroot reconstruction', failures)
      ld = [4.0_dp, 5.0_dp, 6.0_dp]; sd = [1.0_dp, 0.5_dp]
      call tri_cholesky(ld, sd, rd, ro, status)
      call assert_true(status == 0, 'trichol status', failures)
      call assert_close(rd(1)**2, ld(1), 1.0e-12_dp, 'trichol first diagonal', failures)
   end subroutine test_linalg

   subroutine test_smooths(failures)
      integer, intent(inout) :: failures
      real(dp) :: x(51), y(51)
      real(dp), allocatable :: b(:, :), bp(:, :), bc(:, :), bt(:, :), pens(:, :, :)
      type(smooth_spec_t) :: cr, ps, cyc
      integer :: i, status
      do i = 1, 51
         x(i) = real(i - 1, dp) / 50.0_dp
         y(i) = 2.0_dp * x(i)
      end do
      call construct_cr_smooth(x, 8, b, cr, status)
      call assert_true(status == 0 .and. all(shape(b) == [51,8]), 'cr smooth dimensions', failures)
      call assert_true(maxval(abs(sum(b, dim=1))) < 1.0e-10_dp, 'cr smooth centered', failures)
      call predict_smooth(cr, x=x, basis=bp, status=status)
      call assert_true(maxval(abs(bp - b)) < 1.0e-11_dp, 'cr prediction reproduces basis', failures)
      call construct_ps_smooth(x, 9, bp, ps, status)
      call assert_true(status == 0 .and. all(shape(bp) == [51,9]), 'ps smooth dimensions', failures)
      call construct_cyclic_smooth(2.0_dp * pi_dp * x, 6, bc, cyc, status, period=2.0_dp*pi_dp, origin=0.0_dp)
      call assert_true(status == 0, 'cyclic smooth status', failures)
      call tensor_smooth(b(:,1:3), cr%penalties(1:3,1:3,1), bp(:,1:4), ps%penalties(1:4,1:4,1), bt, pens, status)
      call assert_true(status == 0 .and. all(shape(bt) == [51,12]) .and. size(pens,3) == 2, 'tensor smooth', failures)
   end subroutine test_smooths

   subroutine test_gaussian_gam(failures)
      integer, intent(inout) :: failures
      integer, parameter :: n = 140, k = 12
      real(dp) :: x(n), y(n), truth(n)
      real(dp), allocatable :: b(:, :), design(:, :), s(:, :, :), fulls(:, :)
      type(smooth_spec_t) :: ps
      type(gam_model_t) :: model
      type(family_t) :: fam
      type(gam_control_t) :: ctrl
      integer :: i, status
      real(dp) :: rmse
      do i = 1, n
         x(i) = real(i - 1, dp) / real(n - 1, dp)
         truth(i) = sin(2.0_dp * pi_dp * x(i)) + 0.4_dp * cos(6.0_dp * pi_dp * x(i))
         y(i) = truth(i) + 0.12_dp * sin(31.0_dp * x(i))
      end do
      call construct_ps_smooth(x, k, b, ps, status)
      design = append_columns(reshape([(1.0_dp,i=1,n)],[n,1]), b)
      allocate(s(k+1,k+1,1)); s = 0.0_dp
      fulls = embed_penalty(ps%penalties(:,:,1), 2, k+1)
      s(:,:,1) = fulls
      fam%id = family_gaussian
      ctrl%max_outer = 24; ctrl%max_irls = 40
      call gam_fit(design, y, s, model, status, family=fam, method=method_gcv, control=ctrl)
      call assert_true(status == 0 .and. model%converged, 'gaussian GAM convergence', failures)
      rmse = sqrt(sum((model%fitted - truth)**2) / real(n, dp))
      call assert_true(rmse < 0.18_dp, 'gaussian GAM fit RMSE', failures)
      call assert_true(model%edf > 2.0_dp .and. model%edf < real(k+1,dp), 'gaussian GAM EDF', failures)
      call assert_true(model%lambda(1) > 0.0_dp, 'gaussian smoothing parameter', failures)
   end subroutine test_gaussian_gam


   subroutine test_magic_reference(failures)
      integer, intent(inout) :: failures
      real(dp) :: x(5,2), y(5), s(2,2,1), lambda(1)
      type(gam_model_t) :: model
      integer :: status
      x(:,1) = 1.0_dp
      x(:,2) = [-2.0_dp,-1.0_dp,0.0_dp,1.0_dp,2.0_dp]
      y = [0.0_dp,1.0_dp,1.0_dp,3.0_dp,5.0_dp]
      s = 0.0_dp; s(2,2,1) = 1.0_dp; lambda = 2.0_dp
      call magic_fit(x,y,s,model,status,lambda=lambda,method=method_fixed)
      call assert_true(status == 0, 'magic fixed-lambda status', failures)
      call assert_close(model%coefficients(1), 2.0_dp, 1.0e-8_dp, 'magic intercept reference', failures)
      call assert_close(model%coefficients(2), 1.0_dp, 1.0e-8_dp, 'magic slope reference', failures)
   end subroutine test_magic_reference

   subroutine test_glm_families(failures)
      integer, intent(inout) :: failures
      integer, parameter :: n = 100
      real(dp) :: x(n), design(n,2), yp(n), yb(n)
      real(dp), allocatable :: s(:, :, :)
      type(gam_model_t) :: model
      type(family_t) :: fam
      integer :: i, status
      allocate(s(2,2,0))
      do i = 1, n
         x(i) = -1.0_dp + 2.0_dp * real(i - 1, dp) / real(n - 1, dp)
         design(i,:) = [1.0_dp, x(i)]
         yp(i) = real(max(0, nint(exp(0.2_dp + 0.7_dp*x(i)) + 0.3_dp*sin(real(i,dp)))), dp)
         yb(i) = merge(1.0_dp, 0.0_dp, mod(i,10) < nint(10.0_dp/(1.0_dp+exp(-1.2_dp*x(i)))))
      end do
      fam%id = family_poisson
      call gam_fit(design, yp, s, model, status, family=fam, method=method_fixed)
      call assert_true(status == 0 .and. model%converged, 'Poisson PIRLS', failures)
      call assert_true(model%coefficients(2) > 0.2_dp, 'Poisson slope sign', failures)
      fam%id = family_binomial
      call gam_fit(design, yb, s, model, status, family=fam, method=method_fixed)
      call assert_true(status == 0 .and. model%converged, 'binomial PIRLS', failures)
      call assert_true(model%coefficients(2) > 0.0_dp, 'binomial slope sign', failures)
   end subroutine test_glm_families

   subroutine test_constraints(failures)
      integer, intent(inout) :: failures
      integer, parameter :: n=8
      real(dp) :: x(n,n), y(n), w(n)
      real(dp), allocatable :: a(:, :), b(:), beta(:), s(:, :, :)
      integer :: i, status
      x = 0.0_dp
      do i=1,n; x(i,i)=1.0_dp; end do
      y = [0.0_dp, 1.0_dp, 0.5_dp, 2.0_dp, 1.8_dp, 3.0_dp, 2.7_dp, 4.0_dp]
      w = 1.0_dp
      allocate(s(n,n,0))
      call monotonicity_constraints(n, a, b)
      call pcls_fit(x,y,w,s,[real(dp)::],beta,status,a_ineq=a,b_ineq=b)
      call assert_true(status == 0, 'PCLS status', failures)
      call assert_true(minval(beta(2:) - beta(:n-1)) >= -1.0e-7_dp, 'PCLS monotonicity', failures)
   end subroutine test_constraints

   subroutine test_distributions(failures)
      integer, intent(inout) :: failures
      real(dp) :: p
      real(dp), allocatable :: z(:, :), tw(:)
      real(dp) :: mu(2), cov(2,2)
      integer :: status
      mu = [1.0_dp, -1.0_dp]
      cov = reshape([1.0_dp,0.3_dp,0.3_dp,2.0_dp],[2,2])
      call rmvn(3000,mu,cov,z,status,seed=1234)
      call assert_true(status == 0, 'rmvn status', failures)
      call assert_true(maxval(abs(sum(z,dim=1)/real(size(z,1),dp)-mu)) < 0.08_dp, 'rmvn means', failures)
      p = weighted_chisq_cdf(3.841458820694_dp,[1.0_dp],status=status,n_grid=10000)
      call assert_true(status == 0 .and. abs(p-0.95_dp) < 0.015_dp, 'weighted chi-square CDF', failures)
      call rtweedie(1000,2.0_dp,0.5_dp,1.5_dp,tw,status,seed=99)
      call assert_true(status == 0 .and. minval(tw) >= 0.0_dp, 'Tweedie simulation', failures)
   end subroutine test_distributions


   subroutine test_simulation(failures)
      integer, intent(inout) :: failures
      type(gam_sim_data_t) :: data
      integer :: status
      call gam_sim(200, example=1, distribution=sim_normal, scale=0.2_dp, data=data, status=status, seed=42)
      call assert_true(status == 0 .and. all(shape(data%x) == [200,4]), 'gamSim dimensions', failures)
      call assert_true(sqrt(sum((data%y-data%f)**2)/200.0_dp) < 0.3_dp, 'gamSim noise scale', failures)
   end subroutine test_simulation

end program test_mgcv
