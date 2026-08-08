! SPDX-License-Identifier: LGPL-3.0-only
module pso_lbfgsb
   use pso_kinds, only : dp
   use pso_types, only : pso_objective, pso_gradient
   implicit none
   private
   public :: bounded_lbfgs_grad, bounded_lbfgs_numgrad

contains

   subroutine bounded_lbfgs_grad(fn, gr, x, lower, upper, fnscale, maxit, memory, reltol, f, nfev)
      procedure(pso_objective) :: fn
      procedure(pso_gradient) :: gr
      real(dp), intent(inout) :: x(:)
      real(dp), intent(in) :: lower(:), upper(:), fnscale, reltol
      integer, intent(in) :: maxit, memory
      real(dp), intent(out) :: f
      integer, intent(out) :: nfev

      call optimize_common(fn, x, lower, upper, fnscale, maxit, memory, reltol, f, nfev, gr)
   end subroutine bounded_lbfgs_grad

   subroutine bounded_lbfgs_numgrad(fn, x, lower, upper, fnscale, maxit, memory, reltol, f, nfev)
      procedure(pso_objective) :: fn
      real(dp), intent(inout) :: x(:)
      real(dp), intent(in) :: lower(:), upper(:), fnscale, reltol
      integer, intent(in) :: maxit, memory
      real(dp), intent(out) :: f
      integer, intent(out) :: nfev

      call optimize_common(fn, x, lower, upper, fnscale, maxit, memory, reltol, f, nfev)
   end subroutine bounded_lbfgs_numgrad

   subroutine optimize_common(fn, x, lower, upper, fnscale, maxit, memory, reltol, f, nfev, gr)
      procedure(pso_objective) :: fn
      real(dp), intent(inout) :: x(:)
      real(dp), intent(in) :: lower(:), upper(:), fnscale, reltol
      integer, intent(in) :: maxit, memory
      real(dp), intent(out) :: f
      integer, intent(out) :: nfev
      procedure(pso_gradient), optional :: gr

      integer :: n, m, it, nhist, slot, i, j
      real(dp) :: alpha, fnew, slope, ys, yy, gamma, pgmax
      real(dp), allocatable :: g(:), gnew(:), pg(:), d(:), q(:), r(:)
      real(dp), allocatable :: xnew(:), svec(:), yvec(:)
      real(dp), allocatable :: shist(:,:), yhist(:,:), rho(:), a(:)
      integer, allocatable :: order(:)

      n = size(x)
      m = max(1, memory)
      allocate(g(n), gnew(n), pg(n), d(n), q(n), r(n), xnew(n), svec(n), yvec(n))
      allocate(shist(n,m), yhist(n,m), rho(m), a(m), order(m))
      shist = 0.0_dp
      yhist = 0.0_dp
      rho = 0.0_dp
      x = max(lower, min(upper, x))
      f = fn(x) / fnscale
      nfev = 1
      call evaluate_gradient(fn, x, lower, upper, fnscale, g, gr)
      nhist = 0
      slot = 0

      do it = 1, maxit
         call projected_gradient(x, g, lower, upper, pg)
         pgmax = maxval(abs(pg))
         if (pgmax <= reltol * max(1.0_dp, abs(f))) exit

         q = pg
         do i = 1, nhist
            j = modulo(slot - i, m) + 1
            order(i) = j
            a(i) = rho(j) * dot_product(shist(:,j), q)
            q = q - a(i) * yhist(:,j)
         end do
         if (nhist > 0) then
            j = order(1)
            ys = dot_product(shist(:,j), yhist(:,j))
            yy = dot_product(yhist(:,j), yhist(:,j))
            if (yy > 0.0_dp) then
               gamma = ys / yy
            else
               gamma = 1.0_dp
            end if
         else
            gamma = 1.0_dp
         end if
         r = gamma * q
         do i = nhist, 1, -1
            j = order(i)
            r = r + shist(:,j) * (a(i) - rho(j) * dot_product(yhist(:,j), r))
         end do
         d = -r
         call enforce_active_direction(x, g, lower, upper, d)
         slope = dot_product(g, d)
         if (slope >= -epsilon(1.0_dp) * max(1.0_dp, sqrt(sum(d*d)))) then
            d = -pg
            slope = -dot_product(pg, pg)
         end if
         if (slope >= 0.0_dp) exit

         alpha = 1.0_dp
         do
            xnew = max(lower, min(upper, x + alpha * d))
            if (maxval(abs(xnew - x)) <= epsilon(1.0_dp) * max(1.0_dp, maxval(abs(x)))) then
               alpha = 0.0_dp
               exit
            end if
            fnew = fn(xnew) / fnscale
            nfev = nfev + 1
            if (fnew <= f + 1.0e-4_dp * dot_product(g, xnew - x)) exit
            alpha = 0.5_dp * alpha
            if (alpha < 1.0e-12_dp) then
               alpha = 0.0_dp
               exit
            end if
         end do
         if (alpha <= 0.0_dp) exit

         call evaluate_gradient(fn, xnew, lower, upper, fnscale, gnew, gr)
         svec = xnew - x
         yvec = gnew - g
         ys = dot_product(svec, yvec)
         if (ys > 1.0e-12_dp * sqrt(max(0.0_dp, dot_product(svec,svec) * dot_product(yvec,yvec)))) then
            slot = modulo(slot, m) + 1
            shist(:,slot) = svec
            yhist(:,slot) = yvec
            rho(slot) = 1.0_dp / ys
            nhist = min(nhist + 1, m)
         end if
         x = xnew
         g = gnew
         f = fnew
      end do
   end subroutine optimize_common

   subroutine evaluate_gradient(fn, z, lower, upper, fnscale, gz, gr)
      procedure(pso_objective) :: fn
      real(dp), intent(in) :: z(:), lower(:), upper(:), fnscale
      real(dp), intent(out) :: gz(:)
      procedure(pso_gradient), optional :: gr
      real(dp), allocatable :: zp(:), zm(:)
      real(dp) :: h, fp, fm
      integer :: i

      if (present(gr)) then
         call gr(z, gz)
         gz = gz / fnscale
         return
      end if

      allocate(zp(size(z)), zm(size(z)))
      zp = z
      zm = z
      do i = 1, size(z)
         h = sqrt(epsilon(1.0_dp)) * max(1.0_dp, abs(z(i)))
         zp(i) = min(upper(i), z(i) + h)
         zm(i) = max(lower(i), z(i) - h)
         if (abs(zp(i) - zm(i)) <= tiny(1.0_dp)) then
            gz(i) = 0.0_dp
         else
            fp = fn(zp) / fnscale
            fm = fn(zm) / fnscale
            gz(i) = (fp - fm) / (zp(i) - zm(i))
         end if
         zp(i) = z(i)
         zm(i) = z(i)
      end do
   end subroutine evaluate_gradient

   subroutine projected_gradient(x, g, lower, upper, pg)
      real(dp), intent(in) :: x(:), g(:), lower(:), upper(:)
      real(dp), intent(out) :: pg(:)
      integer :: i
      real(dp) :: tol

      pg = g
      do i = 1, size(x)
         tol = 16.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(x(i)), abs(lower(i)), abs(upper(i)))
         if (x(i) <= lower(i) + tol .and. g(i) > 0.0_dp) pg(i) = 0.0_dp
         if (x(i) >= upper(i) - tol .and. g(i) < 0.0_dp) pg(i) = 0.0_dp
      end do
   end subroutine projected_gradient

   subroutine enforce_active_direction(x, g, lower, upper, d)
      real(dp), intent(in) :: x(:), g(:), lower(:), upper(:)
      real(dp), intent(inout) :: d(:)
      integer :: i
      real(dp) :: tol

      do i = 1, size(x)
         tol = 16.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(x(i)), abs(lower(i)), abs(upper(i)))
         if (x(i) <= lower(i) + tol .and. g(i) > 0.0_dp .and. d(i) < 0.0_dp) d(i) = 0.0_dp
         if (x(i) >= upper(i) - tol .and. g(i) < 0.0_dp .and. d(i) > 0.0_dp) d(i) = 0.0_dp
      end do
   end subroutine enforce_active_direction

end module pso_lbfgsb
