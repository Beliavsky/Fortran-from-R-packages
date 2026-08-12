program test_sampling
   use ceoptim
   implicit none
   type(rng_state) :: rng
   type(tmvn_result) :: out
   real(dp), allocatable :: x(:, :), sigma(:, :), a(:, :), b(:)
   real(dp) :: alpha(3)
   integer :: status, i

   call rng_seed(rng, 24680_i64)
   alpha = [0.5_dp, 2.0_dp, 3.5_dp]
   call dirichlet_rand(alpha, 500, rng, x, status)
   if (status /= 0) error stop 'dirichlet status'
   do i = 1, size(x,1)
      if (abs(sum(x(i,:))-1.0_dp) > 1.0e-12_dp) error stop 'dirichlet row sum'
      if (any(x(i,:) <= 0.0_dp)) error stop 'dirichlet positivity'
   end do

   allocate(sigma(2,2), a(2,2), b(2))
   sigma = 0.0_dp
   sigma(1,1) = 1.0_dp
   sigma(2,2) = 1.0_dp
   a = reshape([-1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2,2])
   b = [0.0_dp, 0.5_dp]
   call rtmvnorm(600, [0.0_dp,0.0_dp], sigma, rng, out, a=a, b=b)
   if (out%status /= 0) error stop out%message
   if (out%nar /= 600) error stop 'accept/reject did not fill requested sample'
   if (out%ngibbs /= 0) error stop 'unexpected Gibbs fallback'
   do i = 1, size(out%x,1)
      if (out%x(i,1) < -1.0e-10_dp) error stop 'accept/reject lower constraint violation'
      if (out%x(i,2) > 0.5_dp + 1.0e-10_dp) error stop 'accept/reject upper constraint violation'
   end do

   call rtmvnorm(600, [0.0_dp,0.0_dp], sigma, rng, out, a=a, b=b, rho_thr=1.0_dp)
   if (out%status /= 0) error stop out%message
   if (out%nar /= 0) error stop 'forced Gibbs should skip accept/reject'
   if (out%ngibbs /= 600) error stop 'Gibbs count mismatch'
   do i = 1, size(out%x,1)
      if (out%x(i,1) < -1.0e-10_dp) error stop 'Gibbs lower constraint violation'
      if (out%x(i,2) > 0.5_dp + 1.0e-10_dp) error stop 'Gibbs upper constraint violation'
   end do
   print *, 'test_sampling: PASS', sum(out%x(:,1))/600.0_dp, sum(out%x(:,2))/600.0_dp
end program test_sampling
