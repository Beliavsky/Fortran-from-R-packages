module binary_fixture
   use trustoptim
   use trustoptim_binary
   implicit none
   type(binary_data), save :: dat
   type(binary_priors), save :: pri
contains
   subroutine init_fixture()
      allocate(dat%y(4), dat%x(2,4))
      dat%y = [7.0_dp, 3.0_dp, 8.0_dp, 2.0_dp]
      dat%x = reshape([1.0_dp,0.2_dp, -0.5_dp,1.0_dp, 0.8_dp,-1.2_dp, -1.0_dp,-0.3_dp],[2,4])
      dat%trials = 10
      allocate(pri%inv_sigma(2,2), pri%inv_omega(2,2))
      pri%inv_sigma = 0.0_dp
      pri%inv_omega = 0.0_dp
      pri%inv_sigma(1,1) = 1.5_dp
      pri%inv_sigma(2,2) = 0.8_dp
      pri%inv_omega(1,1) = 0.4_dp
      pri%inv_omega(2,2) = 0.6_dp
   end subroutine init_fixture

   function obj(p) result(f)
      real(dp), intent(in) :: p(:)
      real(dp) :: f
      f = binary_value(p, dat, pri)
   end function obj

   subroutine grad(p, g)
      real(dp), intent(in) :: p(:)
      real(dp), intent(out) :: g(:)
      call binary_gradient(p, dat, pri, g)
   end subroutine grad

   subroutine hess(p, h)
      real(dp), intent(in) :: p(:)
      type(sparse_symmetric_matrix), intent(inout) :: h
      call binary_hessian(p, dat, pri, h)
   end subroutine hess
end module binary_fixture

program test_binary
   use trustoptim
   use trustoptim_binary
   use binary_fixture
   implicit none
   real(dp) :: p(10), g(10), gn(10), pp(10), pm(10), eps
   type(sparse_symmetric_matrix) :: h
   type(trustoptim_control) :: con
   type(trustoptim_result) :: res
   integer :: i

   call init_fixture()
   p = [0.1_dp,-0.2_dp, 0.3_dp,0.1_dp, -0.1_dp,0.4_dp, 0.2_dp,-0.3_dp, 0.05_dp,-0.1_dp]
   call grad(p,g)
   eps = 1.0e-6_dp
   do i = 1, size(p)
      pp = p
      pm = p
      pp(i) = pp(i) + eps
      pm(i) = pm(i) - eps
      gn(i) = (obj(pp)-obj(pm))/(2.0_dp*eps)
   end do
   if (maxval(abs(g-gn)) > 2.0e-7_dp) error stop 'binary gradient finite-difference mismatch'

   call hess(p,h)
   if (h%nnz <= 0) error stop 'binary Hessian unexpectedly empty'

   p = 0.0_dp
   con%function_scale_factor = -1.0_dp
   con%maxit = 2000
   con%prec = 2.0e-7_dp
   con%cg_tol = 1.0e-9_dp
   con%stop_trust_radius = 1.0e-11_dp
   con%preconditioner = 1
   call trust_optim_sparse(p, obj, grad, hess, res, con)
   if (res%status /= trust_success) then
      write(*,*) trim(res%status_message()), sqrt(sum(res%gradient**2)), res%iterations
      error stop 'binary sparse optimization failed'
   end if
   if (sqrt(sum(res%gradient**2))/sqrt(real(size(p),dp)) > 2.1e-7_dp) then
      error stop 'binary optimum gradient too large'
   end if
   write(*,*) 'PASS binary objective/gradient/Hessian and sparse maximization'
end program test_binary
