module countdm_optimizer
   use countdm_kinds, only: dp
   implicit none
   private
   public :: bfgs_minimize

   abstract interface
      function objective_fn(x) result(f)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp) :: f
      end function objective_fn
   end interface

contains

   subroutine numerical_gradient(fun, x, g)
      procedure(objective_fn) :: fun
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: g(:)
      real(dp), allocatable :: xp(:), xm(:)
      real(dp) :: h
      integer :: i
      allocate(xp(size(x)), xm(size(x)))
      do i = 1, size(x)
         h = 1.0e-6_dp * max(1.0_dp, abs(x(i)))
         xp = x; xm = x
         xp(i) = xp(i) + h
         xm(i) = xm(i) - h
         g(i) = (fun(xp) - fun(xm)) / (2.0_dp * h)
      end do
   end subroutine numerical_gradient

   subroutine bfgs_minimize(fun, x, f, converged, iterations, max_iter, tol)
      procedure(objective_fn) :: fun
      real(dp), intent(inout) :: x(:)
      real(dp), intent(out) :: f
      logical, intent(out) :: converged
      integer, intent(out) :: iterations
      integer, intent(in), optional :: max_iter
      real(dp), intent(in), optional :: tol
      real(dp), allocatable :: hmat(:, :), g(:), gn(:), p(:), xn(:), s(:), y(:), hy(:)
      real(dp) :: ft, fn, alpha, c1, ys, yhy, tolerance, gnorm
      integer :: n, it, imax, i
      n = size(x)
      imax = 500; if (present(max_iter)) imax = max_iter
      tolerance = 1.0e-8_dp; if (present(tol)) tolerance = tol
      allocate(hmat(n, n), g(n), gn(n), p(n), xn(n), s(n), y(n), hy(n))
      hmat = 0.0_dp
      do i = 1, n
         hmat(i, i) = 1.0_dp
      end do
      ft = fun(x)
      call numerical_gradient(fun, x, g)
      converged = .false.
      iterations = 0
      do it = 1, imax
         iterations = it
         gnorm = maxval(abs(g))
         if (gnorm <= tolerance) then
            converged = .true.
            exit
         end if
         p = -matmul(hmat, g)
         if (dot_product(p, g) >= -1.0e-12_dp * max(1.0_dp, dot_product(g, g))) then
            p = -g
            hmat = 0.0_dp
            do i = 1, n
               hmat(i, i) = 1.0_dp
            end do
         end if
         alpha = 1.0_dp
         c1 = 1.0e-4_dp
         do
            xn = x + alpha * p
            fn = fun(xn)
            if (fn <= ft + c1 * alpha * dot_product(g, p)) exit
            alpha = 0.5_dp * alpha
            if (alpha < 1.0e-10_dp) exit
         end do
         if (alpha < 1.0e-10_dp) then
            if (gnorm <= 1.0e-4_dp .or. maxval(abs(p)) <= 1.0e-7_dp) converged = .true.
            exit
         end if
         call numerical_gradient(fun, xn, gn)
         s = xn - x
         y = gn - g
         ys = dot_product(y, s)
         if (ys > 1.0e-12_dp * sqrt(max(dot_product(y, y) * dot_product(s, s), tiny(1.0_dp)))) then
            hy = matmul(hmat, y)
            yhy = dot_product(y, hy)
            hmat = hmat + ((ys + yhy) / (ys * ys)) * outer(s, s) &
               - (outer(hy, s) + outer(s, hy)) / ys
         else
            hmat = 0.0_dp
            do i = 1, n
               hmat(i, i) = 1.0_dp
            end do
         end if
         x = xn
         g = gn
         if (abs(ft - fn) <= tolerance * (1.0_dp + abs(ft))) then
            ft = fn
            converged = .true.
            exit
         end if
         ft = fn
      end do
      f = ft
   contains
      pure function outer(a, b) result(c)
         real(dp), intent(in) :: a(:), b(:)
         real(dp) :: c(size(a), size(b))
         integer :: ii, jj
         do jj = 1, size(b)
            do ii = 1, size(a)
               c(ii, jj) = a(ii) * b(jj)
            end do
         end do
      end function outer
   end subroutine bfgs_minimize

end module countdm_optimizer
