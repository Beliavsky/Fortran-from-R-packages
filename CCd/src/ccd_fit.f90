module ccd_fit
   use ccd_kinds, only : dp
   use ccd_optimize, only : nelder_mead, golden_maximize
   implicit none
   private
   real(dp), parameter :: pi = acos(-1.0_dp)

   type, public :: cc_fit_result
      real(dp) :: mu = 0.0_dp
      real(dp) :: lambda = 1.0_dp
      real(dp) :: loglik = -huge(1.0_dp)
      integer :: iterations = 0
      integer :: status = 0
   end type

   type, public :: cc_reg_result
      real(dp) :: lambda = 1.0_dp
      real(dp), allocatable :: beta(:)
      real(dp) :: loglik = -huge(1.0_dp)
      integer :: iterations = 0
      integer :: status = 0
   end type

   type, public :: loc0_test_result
      real(dp) :: statistic = 0.0_dp
      real(dp) :: p_value = 1.0_dp
   end type

   type :: sample_context
      real(dp), allocatable :: y(:)
   end type
   type :: regression_context
      real(dp), allocatable :: y(:), x(:,:)
   end type

   public :: cc_mle, cc_mle0, cc_reg, loc0_test
contains
   function cc_mle(y) result(res)
      real(dp), intent(in) :: y(:)
      type(cc_fit_result) :: res
      type(sample_context) :: ctx
      real(dp) :: par(2), fbest, s
      integer :: n
      n = size(y)
      if (n < 1) then
         res%status = 2
         return
      end if
      allocate(ctx%y, source=y)
      s = robust_scale(y)
      par = [sum(y)/real(n,dp), log(max(s, sqrt(epsilon(1.0_dp))))]
      call nelder_mead(cc_objective, ctx, par, fbest, maxit=5000, tol=1.0e-10_dp, &
         iterations=res%iterations, status=res%status)
      res%mu = par(1)
      res%lambda = exp(par(2))
      res%loglik = -fbest - real(n,dp)*log(pi)
   end function cc_mle

   function cc_mle0(y, tol) result(res)
      real(dp), intent(in) :: y(:)
      real(dp), intent(in), optional :: tol
      type(cc_fit_result) :: res
      type(sample_context) :: ctx
      real(dp) :: fmax, eps
      integer :: n
      n = size(y)
      if (n < 1) then
         res%status = 2
         return
      end if
      allocate(ctx%y, source=y)
      eps = 1.0e-7_dp
      if (present(tol)) eps = tol
      call golden_maximize(cc0_loglik_part, ctx, 1.0e-12_dp, 1000.0_dp, res%lambda, fmax, &
         tol=eps, iterations=res%iterations)
      res%mu = 0.0_dp
      res%loglik = fmax - real(n,dp)*log(pi)
      res%status = 0
   end function cc_mle0

   function cc_reg(y, x, tol) result(res)
      real(dp), intent(in) :: y(:), x(:,:)
      real(dp), intent(in), optional :: tol
      type(cc_reg_result) :: res
      type(regression_context) :: ctx
      real(dp), allocatable :: design(:,:), par(:), beta0(:)
      real(dp) :: s, fbest, prev, eps
      integer :: n, p, it, st, total_it

      n = size(y); p = size(x,2)
      if (size(x,1) /= n .or. n < 1) then
         res%status = 2
         return
      end if
      allocate(design(n,p+1), beta0(p+1), par(p+2))
      design(:,1) = 1.0_dp
      if (p > 0) design(:,2:) = x
      call ols_fit(design, y, beta0, st)
      if (st /= 0) beta0 = 0.0_dp
      s = robust_scale(y)
      par(1) = log(max(s, sqrt(epsilon(1.0_dp))))
      par(2:) = beta0
      allocate(ctx%y, source=y)
      allocate(ctx%x, source=design)
      eps = 1.0e-6_dp
      if (present(tol)) eps = tol
      total_it = 0
      prev = huge(1.0_dp)
      do
         call nelder_mead(cc_reg_objective, ctx, par, fbest, maxit=5000, tol=1.0e-10_dp, &
            iterations=it, status=st)
         total_it = total_it + it
         if (prev - fbest <= eps) exit
         prev = fbest
         if (total_it >= 20000) exit
      end do
      res%lambda = exp(par(1))
      allocate(res%beta(size(par)-1))
      res%beta = par(2:)
      res%loglik = -fbest - real(n,dp)*log(pi)
      res%iterations = total_it
      res%status = st
   end function cc_reg

   function loc0_test(y, tol) result(ans)
      real(dp), intent(in) :: y(:)
      real(dp), intent(in), optional :: tol
      type(loc0_test_result) :: ans
      type(cc_fit_result) :: fit0, fit1
      if (present(tol)) then
         fit0 = cc_mle0(y, tol)
      else
         fit0 = cc_mle0(y)
      end if
      fit1 = cc_mle(y)
      ans%statistic = max(0.0_dp, 2.0_dp*(fit1%loglik-fit0%loglik))
      ans%p_value = erfc(sqrt(0.5_dp*ans%statistic))
   end function loc0_test

   function cc_objective(par, context) result(f)
      real(dp), intent(in) :: par(:)
      class(*), intent(in) :: context
      real(dp) :: f, lambda, mu
      integer :: n
      select type(context)
      type is (sample_context)
         n = size(context%y)
         mu = par(1); lambda = exp(par(2))
         f = -real(n,dp)*log(tanh(lambda*pi)) - real(n,dp)*log(lambda) &
            + sum(log(lambda*lambda + (context%y-mu)**2))
      class default
         f = huge(1.0_dp)
      end select
   end function cc_objective

   function cc0_loglik_part(lambda, context) result(f)
      real(dp), intent(in) :: lambda
      class(*), intent(in) :: context
      real(dp) :: f
      integer :: n
      select type(context)
      type is (sample_context)
         if (lambda <= 0.0_dp) then
            f = -huge(1.0_dp)
         else
            n = size(context%y)
            f = real(n,dp)*log(tanh(lambda*pi)) + real(n,dp)*log(lambda) &
               - sum(log(lambda*lambda + context%y**2))
         end if
      class default
         f = -huge(1.0_dp)
      end select
   end function cc0_loglik_part

   function cc_reg_objective(par, context) result(f)
      real(dp), intent(in) :: par(:)
      class(*), intent(in) :: context
      real(dp) :: f, lambda
      real(dp), allocatable :: mu(:)
      integer :: n
      select type(context)
      type is (regression_context)
         n = size(context%y)
         lambda = exp(par(1))
         mu = matmul(context%x, par(2:))
         f = -real(n,dp)*log(tanh(lambda*pi)) - real(n,dp)*log(lambda) &
            + sum(log(lambda*lambda + (context%y-mu)**2))
      class default
         f = huge(1.0_dp)
      end select
   end function cc_reg_objective

   function robust_scale(y) result(s)
      real(dp), intent(in) :: y(:)
      real(dp) :: s
      real(dp), allocatable :: z(:)
      integer :: n, i1, i3
      n = size(y); z = y
      call sort_real(z)
      i1 = max(1, min(n, ceiling(0.25_dp*real(n,dp))))
      i3 = max(1, min(n, ceiling(0.75_dp*real(n,dp))))
      s = 0.5_dp*abs(z(i3)-z(i1))
   end function robust_scale

   subroutine sort_real(x)
      real(dp), intent(inout) :: x(:)
      integer :: i, j
      real(dp) :: key
      do i = 2, size(x)
         key = x(i); j = i - 1
         do while (j >= 1)
            if (x(j) <= key) exit
            x(j+1) = x(j); j = j - 1
         end do
         x(j+1) = key
      end do
   end subroutine sort_real

   subroutine ols_fit(x, y, beta, status)
      real(dp), intent(in) :: x(:,:), y(:)
      real(dp), intent(out) :: beta(:)
      integer, intent(out) :: status
      real(dp), allocatable :: a(:,:), b(:)
      allocate(a(size(x,2),size(x,2)), b(size(x,2)))
      a = matmul(transpose(x), x)
      b = matmul(transpose(x), y)
      call solve_linear(a, b, beta, status)
   end subroutine ols_fit

   subroutine solve_linear(a, b, x, status)
      real(dp), intent(in) :: a(:,:), b(:)
      real(dp), intent(out) :: x(:)
      integer, intent(out) :: status
      real(dp), allocatable :: m(:,:), rhs(:), rowtmp(:)
      real(dp) :: piv, fac, tmp
      integer :: n, i, j, k, ip
      n = size(b); m = a; rhs = b; status = 0
      allocate(rowtmp(n))
      do k = 1, n
         ip = k
         do i = k+1, n
            if (abs(m(i,k)) > abs(m(ip,k))) ip = i
         end do
         if (abs(m(ip,k)) <= 100.0_dp*epsilon(1.0_dp)) then
            status = 1; x = 0.0_dp; return
         end if
         if (ip /= k) then
            rowtmp = m(k,:); m(k,:) = m(ip,:); m(ip,:) = rowtmp
            tmp = rhs(k); rhs(k) = rhs(ip); rhs(ip) = tmp
         end if
         piv = m(k,k)
         do i = k+1, n
            fac = m(i,k)/piv
            m(i,k:n) = m(i,k:n) - fac*m(k,k:n)
            rhs(i) = rhs(i) - fac*rhs(k)
         end do
      end do
      do i = n, 1, -1
         tmp = rhs(i)
         do j = i+1, n
            tmp = tmp - m(i,j)*x(j)
         end do
         x(i) = tmp/m(i,i)
      end do
   end subroutine solve_linear
end module ccd_fit
