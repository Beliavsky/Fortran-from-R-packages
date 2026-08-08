program test_preconditioners_scale
   use trustoptim
   implicit none
   real(dp) :: x0(3)
   type(trustoptim_control) :: con
   type(trustoptim_result) :: r

   x0 = [8.0_dp, -7.0_dp, 5.0_dp]
   con%maxit = 500
   con%prec = 1.0e-8_dp
   con%cg_tol = 1.0e-10_dp
   con%preconditioner = 1
   call trust_optim_bfgs(x0, convex, convex_grad, r, con)
   if (r%status /= trust_success) error stop 'BFGS Cholesky preconditioner failed'
   if (maxval(abs(r%solution - target())) > 1.0e-5_dp) error stop 'BFGS wrong solution'

   con%function_scale_factor = -1.0_dp
   con%preconditioner = 1
   call trust_optim_sparse(x0, concave, concave_grad, concave_hess, r, con)
   if (r%status /= trust_success) error stop 'Sparse modified Cholesky preconditioner failed'
   if (maxval(abs(r%solution - target())) > 1.0e-7_dp) error stop 'Sparse maximize wrong solution'
   if (abs(r%fval - 10.0_dp) > 1.0e-10_dp) error stop 'unscaled maximum fval wrong'

   write(*,*) 'PASS preconditioners and function scaling'
contains
   pure function target() result(t)
      real(dp) :: t(3)
      t = [1.0_dp, -2.0_dp, 0.5_dp]
   end function target

   function convex(x) result(f)
      real(dp), intent(in) :: x(:)
      real(dp) :: f, d(3)
      d = x - target()
      f = d(1)**2 + 4.0_dp*d(2)**2 + 9.0_dp*d(3)**2
   end function convex

   subroutine convex_grad(x, g)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: g(:)
      real(dp) :: d(3)
      d = x - target()
      g = [2.0_dp*d(1), 8.0_dp*d(2), 18.0_dp*d(3)]
   end subroutine convex_grad

   function concave(x) result(f)
      real(dp), intent(in) :: x(:)
      real(dp) :: f
      f = 10.0_dp - convex(x)
   end function concave

   subroutine concave_grad(x, g)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: g(:)
      call convex_grad(x, g)
      g = -g
   end subroutine concave_grad

   subroutine concave_hess(x, h)
      real(dp), intent(in) :: x(:)
      type(sparse_symmetric_matrix), intent(inout) :: h
      real(dp) :: a(size(x),size(x))
      a = 0.0_dp
      a(1,1) = -2.0_dp
      a(2,2) = -8.0_dp
      a(3,3) = -18.0_dp
      call h%set_from_dense(a)
   end subroutine concave_hess
end program test_preconditioners_scale
