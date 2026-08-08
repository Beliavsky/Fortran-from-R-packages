program test_rosenbrock
   use trustoptim
   implicit none
   real(dp) :: x0(6)
   type(trustoptim_control) :: con
   type(trustoptim_result) :: res

   x0 = [-4.0_dp, 2.0_dp, -1.5_dp, 3.0_dp, -2.5_dp, 1.0_dp]
   con%maxit = 5000
   con%prec = 1.0e-6_dp
   con%cg_tol = 1.0e-8_dp
   con%stop_trust_radius = 1.0e-10_dp
   con%report_level = 0

   call trust_optim(x0, rosen, rosen_grad, 'BFGS', res, con)
   call check(res, 'BFGS')

   call trust_optim_sr1(x0, rosen, rosen_grad, res, con)
   call check(res, 'SR1')

   con%preconditioner = 0
   call trust_optim(x0, rosen, rosen_grad, rosen_hess, res, con)
   call check(res, 'Sparse')

contains

   function rosen(x) result(f)
      real(dp), intent(in) :: x(:)
      real(dp) :: f
      integer :: i
      f = 0.0_dp
      do i = 1, size(x), 2
         f = f + 100.0_dp * (x(i)**2 - x(i+1))**2 + (x(i)-1.0_dp)**2
      end do
   end function rosen

   subroutine rosen_grad(x, g)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: g(:)
      integer :: i
      do i = 1, size(x), 2
         g(i) = 400.0_dp * (x(i)**2 - x(i+1)) * x(i) + 2.0_dp * (x(i)-1.0_dp)
         g(i+1) = -200.0_dp * (x(i)**2 - x(i+1))
      end do
   end subroutine rosen_grad

   subroutine rosen_hess(x, h)
      real(dp), intent(in) :: x(:)
      type(sparse_symmetric_matrix), intent(inout) :: h
      real(dp) :: a(size(x),size(x))
      integer :: i
      a = 0.0_dp
      do i = 1, size(x), 2
         a(i,i) = 1200.0_dp*x(i)**2 - 400.0_dp*x(i+1) + 2.0_dp
         a(i+1,i) = -400.0_dp*x(i)
         a(i,i+1) = a(i+1,i)
         a(i+1,i+1) = 200.0_dp
      end do
      call h%set_from_dense(a)
   end subroutine rosen_hess

   subroutine check(r, label)
      type(trustoptim_result), intent(in) :: r
      character(len=*), intent(in) :: label
      real(dp) :: ng
      ng = sqrt(sum(r%gradient*r%gradient))
      if (r%status /= trust_success) then
         write(*,*) trim(label), trim(r%status_message()), r%iterations, r%trust_radius, ng, r%fval
         error stop 'optimizer did not report success'
      end if
      if (ng > 5.0e-4_dp) error stop 'gradient too large'
      if (maxval(abs(r%solution - 1.0_dp)) > 1.0e-3_dp) error stop 'wrong Rosenbrock solution'
      write(*,'(a,1x,a,1x,i0,1x,es12.4)') 'PASS', trim(label), r%iterations, ng
   end subroutine check
end program test_rosenbrock
