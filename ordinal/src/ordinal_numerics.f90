! Numerical support for the modern Fortran translation of ordinal.
! Copyright (C) 2011-2026 R. H. B. Christensen; translation (C) 2026.
! Distributed under GPL-2.0-or-later.
module ordinal_numerics
   use ordinal_kinds, only : dp
   implicit none
   private
   type, abstract, public :: objective_type
   contains
      procedure(objective_value_interface), deferred :: value
   end type objective_type
   abstract interface
      function objective_value_interface(self, x) result(f)
         import dp, objective_type
         class(objective_type), intent(in) :: self !! Immutable numerical objective and its problem-specific data.
         real(dp), intent(in) :: x(:) !! Parameter vector at which the objective is evaluated.
         real(dp) :: f
      end function objective_value_interface
   end interface
   public :: bfgs_minimize, numerical_gradient, numerical_hessian, invert_matrix
   public :: solve_linear_system, cholesky_factor, solve_cholesky, logdet_spd
   public :: symmetric_eigenvalues, hessian_diagnostics
contains
   subroutine numerical_gradient(objective, x, grad, step)
      class(objective_type), intent(in) :: objective !! Numerical objective whose gradient is approximated.
      real(dp), intent(in) :: x(:) !! Parameter vector at which the gradient is approximated.
      real(dp), intent(out) :: grad(:) !! Central-difference gradient; same length as x.
      real(dp), intent(in), optional :: step !! Relative finite-difference step; defaults to sqrt(machine epsilon).
      real(dp), allocatable :: xp(:), xm(:)
      real(dp) :: h, base_step, fp, fm, f0, valid_limit
      integer :: j, attempt
      base_step = sqrt(epsilon(1.0_dp))
      if (present(step)) base_step = step
      allocate(xp(size(x)), xm(size(x)))
      f0 = objective%value(x)
      valid_limit = sqrt(huge(1.0_dp))
      do j = 1, size(x)
         h = base_step*max(1.0_dp, abs(x(j)))
         do attempt = 1, 12
            xp = x
            xm = x
            xp(j) = xp(j) + h
            xm(j) = xm(j) - h
            fp = objective%value(xp)
            fm = objective%value(xm)
            if (abs(fp) < valid_limit .and. abs(fm) < valid_limit) then
               grad(j) = (fp - fm)/(2.0_dp*h)
               exit
            else if (abs(fp) < valid_limit .and. abs(f0) < valid_limit) then
               grad(j) = (fp - f0)/h
               exit
            else if (abs(fm) < valid_limit .and. abs(f0) < valid_limit) then
               grad(j) = (f0 - fm)/h
               exit
            end if
            h = 0.5_dp*h
         end do
         if (attempt > 12) grad(j) = 0.0_dp
      end do
   end subroutine numerical_gradient

   subroutine numerical_hessian(objective, x, hess, step)
      class(objective_type), intent(in) :: objective !! Numerical objective whose Hessian is approximated.
      real(dp), intent(in) :: x(:) !! Parameter vector at which the Hessian is approximated.
      real(dp), intent(out) :: hess(:, :) !! Symmetric central-difference Hessian with dimensions size(x)-by-size(x).
      real(dp), intent(in), optional :: step !! Relative finite-difference step; defaults to epsilon^(1/4).
      real(dp), allocatable :: xpp(:), xpm(:), xmp(:), xmm(:), xp(:), xm(:)
      real(dp) :: hi, hj, f0, base_step
      integer :: i, j
      base_step = epsilon(1.0_dp)**0.25_dp
      if (present(step)) base_step = step
      allocate(xpp(size(x)), xpm(size(x)), xmp(size(x)), xmm(size(x)), xp(size(x)), xm(size(x)))
      f0 = objective%value(x)
      do i = 1, size(x)
         hi = base_step*max(1.0_dp, abs(x(i)))
         xp = x
         xm = x
         xp(i) = xp(i) + hi
         xm(i) = xm(i) - hi
         hess(i, i) = (objective%value(xp) - 2.0_dp*f0 + objective%value(xm))/(hi*hi)
         do j = i + 1, size(x)
            hj = base_step*max(1.0_dp, abs(x(j)))
            xpp = x
            xpm = x
            xmp = x
            xmm = x
            xpp(i) = xpp(i) + hi
            xpp(j) = xpp(j) + hj
            xpm(i) = xpm(i) + hi
            xpm(j) = xpm(j) - hj
            xmp(i) = xmp(i) - hi
            xmp(j) = xmp(j) + hj
            xmm(i) = xmm(i) - hi
            xmm(j) = xmm(j) - hj
            hess(i, j) = (objective%value(xpp) - objective%value(xpm) - objective%value(xmp) + &
                          objective%value(xmm))/(4.0_dp*hi*hj)
            hess(j, i) = hess(i, j)
         end do
      end do
   end subroutine numerical_hessian

   subroutine invert_matrix(a, ainv, status)
      real(dp), intent(in) :: a(:, :) !! Square matrix to invert.
      real(dp), intent(out) :: ainv(:, :) !! Matrix inverse when status is zero.
      integer, intent(out) :: status !! Zero on success; one for dimension errors and two for numerical singularity.
      real(dp), allocatable :: aug(:, :), rowtmp(:)
      real(dp) :: pivot, factor, scale
      integer :: n, i, j, k, pivrow
      n = size(a, 1)
      if (size(a, 2) /= n .or. size(ainv, 1) /= n .or. size(ainv, 2) /= n) then
         status = 1
         return
      end if
      if (n == 0) then
         status = 0
         return
      end if
      allocate(aug(n, 2*n), rowtmp(2*n))
      aug(:, :n) = a
      aug(:, n + 1:) = 0.0_dp
      do i = 1, n
         aug(i, n + i) = 1.0_dp
      end do
      scale = max(1.0_dp, maxval(abs(a)))
      do i = 1, n
         pivrow = i
         do k = i + 1, n
            if (abs(aug(k, i)) > abs(aug(pivrow, i))) pivrow = k
         end do
         if (abs(aug(pivrow, i)) <= 100.0_dp*epsilon(1.0_dp)*scale) then
            status = 2
            ainv = 0.0_dp
            return
         end if
         if (pivrow /= i) then
            rowtmp = aug(i, :)
            aug(i, :) = aug(pivrow, :)
            aug(pivrow, :) = rowtmp
         end if
         pivot = aug(i, i)
         aug(i, :) = aug(i, :)/pivot
         do j = 1, n
            if (j == i) cycle
            factor = aug(j, i)
            aug(j, :) = aug(j, :) - factor*aug(i, :)
         end do
      end do
      ainv = aug(:, n + 1:)
      status = 0
   end subroutine invert_matrix

   subroutine solve_linear_system(a, b, x, status)
      real(dp), intent(in) :: a(:, :) !! Square coefficient matrix.
      real(dp), intent(in) :: b(:) !! Right-hand-side vector with length size(a,1).
      real(dp), intent(out) :: x(:) !! Solution vector when status is zero.
      integer, intent(out) :: status !! Zero on success; one for dimension errors and two for numerical singularity.
      real(dp), allocatable :: aug(:, :), rowtmp(:)
      real(dp) :: factor, pivot, scale
      integer :: n, i, j, k, pivrow
      n = size(a, 1)
      if (size(a, 2) /= n .or. size(b) /= n .or. size(x) /= n) then
         status = 1
         return
      end if
      if (n == 0) then
         status = 0
         return
      end if
      allocate(aug(n, n + 1), rowtmp(n + 1))
      aug(:, :n) = a
      aug(:, n + 1) = b
      scale = max(1.0_dp, maxval(abs(a)))
      do i = 1, n
         pivrow = i
         do k = i + 1, n
            if (abs(aug(k, i)) > abs(aug(pivrow, i))) pivrow = k
         end do
         if (abs(aug(pivrow, i)) <= 100.0_dp*epsilon(1.0_dp)*scale) then
            status = 2
            x = 0.0_dp
            return
         end if
         if (pivrow /= i) then
            rowtmp = aug(i, :)
            aug(i, :) = aug(pivrow, :)
            aug(pivrow, :) = rowtmp
         end if
         pivot = aug(i, i)
         aug(i, i:n + 1) = aug(i, i:n + 1)/pivot
         do j = i + 1, n
            factor = aug(j, i)
            aug(j, i:n + 1) = aug(j, i:n + 1) - factor*aug(i, i:n + 1)
         end do
      end do
      do i = n, 1, -1
         x(i) = aug(i, n + 1)
         if (i < n) x(i) = x(i) - dot_product(aug(i, i + 1:n), x(i + 1:n))
      end do
      status = 0
   end subroutine solve_linear_system

   subroutine cholesky_factor(a, l, status)
      real(dp), intent(in) :: a(:, :) !! Symmetric matrix to factor as L*transpose(L).
      real(dp), intent(out) :: l(:, :) !! Lower-triangular Cholesky factor when status is zero.
      integer, intent(out) :: status !! Zero for positive-definite input; one for dimensions and two otherwise.
      real(dp) :: s, tol
      integer :: n, i, j, k
      n = size(a, 1)
      if (size(a, 2) /= n .or. size(l, 1) /= n .or. size(l, 2) /= n) then
         status = 1
         return
      end if
      l = 0.0_dp
      tol = 100.0_dp*epsilon(1.0_dp)*max(1.0_dp, maxval(abs(a)))
      do i = 1, n
         do j = 1, i
            s = a(i, j)
            do k = 1, j - 1
               s = s - l(i, k)*l(j, k)
            end do
            if (i == j) then
               if (s <= tol) then
                  status = 2
                  l = 0.0_dp
                  return
               end if
               l(i, j) = sqrt(s)
            else
               l(i, j) = s/l(j, j)
            end if
         end do
      end do
      status = 0
   end subroutine cholesky_factor

   subroutine solve_cholesky(l, b, x, status)
      real(dp), intent(in) :: l(:, :) !! Nonsingular lower-triangular Cholesky factor.
      real(dp), intent(in) :: b(:) !! Right-hand side in L*transpose(L)*x=b.
      real(dp), intent(out) :: x(:) !! Solution vector when status is zero.
      integer, intent(out) :: status !! Zero on success; nonzero for dimensions or a zero diagonal.
      real(dp), allocatable :: y(:)
      real(dp) :: tol
      integer :: n, i
      n = size(l, 1)
      if (size(l, 2) /= n .or. size(b) /= n .or. size(x) /= n) then
         status = 1
         return
      end if
      tol = 100.0_dp*epsilon(1.0_dp)*max(1.0_dp, maxval(abs(l)))
      if (any(abs([(l(i, i), i = 1, n)]) <= tol)) then
         status = 2
         x = 0.0_dp
         return
      end if
      allocate(y(n))
      do i = 1, n
         y(i) = b(i)
         if (i > 1) y(i) = y(i) - dot_product(l(i, :i - 1), y(:i - 1))
         y(i) = y(i)/l(i, i)
      end do
      do i = n, 1, -1
         x(i) = y(i)
         if (i < n) x(i) = x(i) - dot_product(l(i + 1:n, i), x(i + 1:n))
         x(i) = x(i)/l(i, i)
      end do
      status = 0
   end subroutine solve_cholesky

   subroutine logdet_spd(a, logdet, status)
      real(dp), intent(in) :: a(:, :) !! Symmetric positive-definite matrix.
      real(dp), intent(out) :: logdet !! Natural logarithm of determinant when status is zero.
      integer, intent(out) :: status !! Zero on success; nonzero when Cholesky factorization fails.
      real(dp), allocatable :: l(:, :)
      integer :: i, n
      n = size(a, 1)
      allocate(l(n, n))
      call cholesky_factor(a, l, status)
      if (status /= 0) then
         logdet = huge(1.0_dp)
         return
      end if
      logdet = 0.0_dp
      do i = 1, n
         logdet = logdet + 2.0_dp*log(l(i, i))
      end do
   end subroutine logdet_spd

   subroutine symmetric_eigenvalues(a, eigenvalues, status)
      real(dp), intent(in) :: a(:, :) !! Real symmetric matrix whose eigenvalues are required.
      real(dp), intent(out) :: eigenvalues(:) !! Eigenvalues in unspecified order; length must equal matrix order.
      integer, intent(out) :: status !! Zero on convergence; one for dimensions and two for Jacobi iteration limit.
      real(dp), allocatable :: work(:, :)
      real(dp) :: app, aqq, apq, tau, t, c, s, aik, akq, offmax, tol
      integer :: n, p, q, k, sweep, max_sweeps
      n = size(a, 1)
      if (size(a, 2) /= n .or. size(eigenvalues) /= n) then
         status = 1
         return
      end if
      if (n == 0) then
         status = 0
         return
      end if
      allocate(work(n, n))
      work = 0.5_dp*(a + transpose(a))
      tol = 100.0_dp*epsilon(1.0_dp)*max(1.0_dp, maxval(abs(work)))
      max_sweeps = max(20, 20*n*n)
      do sweep = 1, max_sweeps
         offmax = 0.0_dp
         p = 1
         q = min(2, n)
         do k = 1, n - 1
            if (maxval(abs(work(k, k + 1:n))) > offmax) then
               q = k + maxloc(abs(work(k, k + 1:n)), dim=1)
               p = k
               offmax = abs(work(p, q))
            end if
         end do
         if (offmax <= tol .or. n == 1) then
            do k = 1, n
               eigenvalues(k) = work(k, k)
            end do
            status = 0
            return
         end if
         app = work(p, p)
         aqq = work(q, q)
         apq = work(p, q)
         tau = (aqq - app)/(2.0_dp*apq)
         if (tau >= 0.0_dp) then
            t = 1.0_dp/(tau + sqrt(1.0_dp + tau*tau))
         else
            t = -1.0_dp/(-tau + sqrt(1.0_dp + tau*tau))
         end if
         c = 1.0_dp/sqrt(1.0_dp + t*t)
         s = t*c
         work(p, p) = app - t*apq
         work(q, q) = aqq + t*apq
         work(p, q) = 0.0_dp
         work(q, p) = 0.0_dp
         do k = 1, n
            if (k == p .or. k == q) cycle
            aik = work(k, p)
            akq = work(k, q)
            work(k, p) = c*aik - s*akq
            work(p, k) = work(k, p)
            work(k, q) = s*aik + c*akq
            work(q, k) = work(k, q)
         end do
      end do
      do k = 1, n
         eigenvalues(k) = work(k, k)
      end do
      status = 2
   end subroutine symmetric_eigenvalues

   subroutine hessian_diagnostics(hessian, min_eigenvalue, max_eigenvalue, condition_number, rank, positive_definite, status)
      real(dp), intent(in) :: hessian(:, :) !! Symmetric Hessian matrix to diagnose.
      real(dp), intent(out) :: min_eigenvalue !! Smallest estimated eigenvalue.
      real(dp), intent(out) :: max_eigenvalue !! Largest estimated eigenvalue.
      real(dp), intent(out) :: condition_number !! Absolute eigenvalue condition estimate, or huge for rank deficiency.
      integer, intent(out) :: rank !! Numerical rank using a tolerance scaled by the largest absolute eigenvalue.
      logical, intent(out) :: positive_definite !! True when all eigenvalues are strictly above the numerical tolerance.
      integer, intent(out) :: status !! Zero when diagnostics succeed; otherwise the eigensolver status.
      real(dp), allocatable :: eigenvalues(:)
      real(dp) :: scale, tol, minabs
      integer :: n
      n = size(hessian, 1)
      if (size(hessian, 2) /= n) then
         status = 1
         min_eigenvalue = 0.0_dp
         max_eigenvalue = 0.0_dp
         condition_number = huge(1.0_dp)
         rank = 0
         positive_definite = .false.
         return
      end if
      allocate(eigenvalues(n))
      call symmetric_eigenvalues(hessian, eigenvalues, status)
      if (status /= 0) then
         min_eigenvalue = minval(eigenvalues)
         max_eigenvalue = maxval(eigenvalues)
         condition_number = huge(1.0_dp)
         rank = 0
         positive_definite = .false.
         return
      end if
      min_eigenvalue = minval(eigenvalues)
      max_eigenvalue = maxval(eigenvalues)
      scale = max(1.0_dp, maxval(abs(eigenvalues)))
      tol = 1000.0_dp*epsilon(1.0_dp)*scale
      rank = count(abs(eigenvalues) > tol)
      positive_definite = min_eigenvalue > tol
      if (rank < n) then
         condition_number = huge(1.0_dp)
      else
         minabs = minval(abs(eigenvalues))
         condition_number = maxval(abs(eigenvalues))/minabs
      end if
   end subroutine hessian_diagnostics

   subroutine bfgs_minimize(objective, x, fval, iterations, status, max_iter, grad_tol)
      class(objective_type), intent(in) :: objective !! Numerical objective to minimize.
      real(dp), intent(inout) :: x(:) !! Parameter vector; starting values on input and fitted values on output.
      real(dp), intent(out) :: fval !! Objective value at the returned parameter vector.
      integer, intent(out) :: iterations !! Number of completed quasi-Newton iterations.
      integer, intent(out) :: status !! Zero for convergence, one for iteration limit, two for failed line search.
      integer, intent(in), optional :: max_iter !! Maximum quasi-Newton iterations; defaults to 300.
      real(dp), intent(in), optional :: grad_tol !! Infinity-norm gradient tolerance; defaults to 1e-7.
      real(dp), allocatable :: hinv(:, :), g(:), gnew(:), p(:), xnew(:), s(:), y(:), v(:, :), ident(:, :)
      real(dp) :: alpha, fnew, c1, ys, tol
      integer :: n, iter, ls, limit, i
      n = size(x)
      limit = 300
      if (present(max_iter)) limit = max_iter
      tol = 1.0e-7_dp
      if (present(grad_tol)) tol = grad_tol
      allocate(hinv(n, n), g(n), gnew(n), p(n), xnew(n), s(n), y(n), v(n, n), ident(n, n))
      hinv = 0.0_dp
      ident = 0.0_dp
      do i = 1, n
         hinv(i, i) = 1.0_dp
         ident(i, i) = 1.0_dp
      end do
      fval = objective%value(x)
      call numerical_gradient(objective, x, g, 1.0e-6_dp)
      c1 = 1.0e-4_dp
      status = 1
      iterations = 0
      do iter = 1, limit
         if (maxval(abs(g)) <= tol) then
            status = 0
            exit
         end if
         p = -matmul(hinv, g)
         if (dot_product(p, g) >= 0.0_dp) p = -g
         alpha = 1.0_dp
         do ls = 1, 40
            xnew = x + alpha*p
            fnew = objective%value(xnew)
            if (fnew < huge(1.0_dp)/100.0_dp .and. fnew <= fval + c1*alpha*dot_product(g, p)) exit
            alpha = 0.5_dp*alpha
         end do
         if (ls > 40) then
            p = -g
            alpha = min(1.0_dp, 1.0_dp/max(1.0_dp, norm2(g)))
            do ls = 1, 60
               xnew = x + alpha*p
               fnew = objective%value(xnew)
               if (fnew < huge(1.0_dp)/100.0_dp .and. fnew < fval) exit
               alpha = 0.5_dp*alpha
            end do
            if (ls > 60) then
               status = 2
               exit
            end if
            hinv = ident
         end if
         call numerical_gradient(objective, xnew, gnew, 1.0e-6_dp)
         s = xnew - x
         y = gnew - g
         ys = dot_product(y, s)
         if (ys > 1.0e-12_dp*max(1.0_dp, norm2(y)*norm2(s))) then
            v = ident - outer_product(s, y)/ys
            hinv = matmul(matmul(v, hinv), transpose(v)) + outer_product(s, s)/ys
         else
            hinv = ident
         end if
         x = xnew
         g = gnew
         fval = fnew
         iterations = iter
      end do
   contains
      pure function outer_product(a, b) result(c)
         real(dp), intent(in) :: a(:) !! Left vector in the dyadic product.
         real(dp), intent(in) :: b(:) !! Right vector in the dyadic product.
         real(dp) :: c(size(a), size(b))
         integer :: ii
         do ii = 1, size(a)
            c(ii, :) = a(ii)*b
         end do
      end function outer_product
   end subroutine bfgs_minimize
end module ordinal_numerics
