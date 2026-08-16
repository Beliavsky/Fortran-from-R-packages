module ccd_optimize
   use ccd_kinds, only : dp
   implicit none
   private
   public :: nelder_mead, golden_maximize

   abstract interface
      function objective_nd(x, context) result(f)
         import dp
         real(dp), intent(in) :: x(:)
         class(*), intent(in) :: context
         real(dp) :: f
      end function objective_nd
      function objective_1d(x, context) result(f)
         import dp
         real(dp), intent(in) :: x
         class(*), intent(in) :: context
         real(dp) :: f
      end function objective_1d
   end interface
contains
   subroutine nelder_mead(fun, context, x, fbest, maxit, tol, iterations, status)
      procedure(objective_nd) :: fun
      class(*), intent(in) :: context
      real(dp), intent(inout) :: x(:)
      real(dp), intent(out) :: fbest
      integer, intent(in), optional :: maxit
      real(dp), intent(in), optional :: tol
      integer, intent(out), optional :: iterations, status
      integer :: n, m, it, max_iter, i, j, ilo, ihi, inhi
      real(dp) :: eps, fr, fe, fc, spread, scale
      real(dp), allocatable :: simplex(:,:), f(:), centroid(:), xr(:), xe(:), xc(:)

      n = size(x); m = n + 1
      max_iter = 5000
      if (present(maxit)) max_iter = maxit
      eps = 1.0e-10_dp
      if (present(tol)) eps = tol
      allocate(simplex(n,m), f(m), centroid(n), xr(n), xe(n), xc(n))
      simplex(:,1) = x
      do j = 2, m
         simplex(:,j) = x
         i = j - 1
         scale = 0.05_dp*abs(x(i))
         if (scale <= tiny(1.0_dp)) scale = 0.00025_dp
         simplex(i,j) = x(i) + scale
      end do
      do j = 1, m
         f(j) = fun(simplex(:,j), context)
      end do

      do it = 1, max_iter
         call order_indices(f, ilo, ihi, inhi)
         spread = maxval(abs(f - f(ilo)))
         if (spread <= eps*(1.0_dp + abs(f(ilo)))) exit
         centroid = 0.0_dp
         do j = 1, m
            if (j /= ihi) centroid = centroid + simplex(:,j)
         end do
         centroid = centroid/real(n,dp)
         xr = centroid + (centroid - simplex(:,ihi))
         fr = fun(xr, context)
         if (fr < f(ilo)) then
            xe = centroid + 2.0_dp*(xr-centroid)
            fe = fun(xe, context)
            if (fe < fr) then
               simplex(:,ihi) = xe; f(ihi) = fe
            else
               simplex(:,ihi) = xr; f(ihi) = fr
            end if
         else if (fr < f(inhi)) then
            simplex(:,ihi) = xr; f(ihi) = fr
         else
            if (fr < f(ihi)) then
               xc = centroid + 0.5_dp*(xr-centroid)
            else
               xc = centroid + 0.5_dp*(simplex(:,ihi)-centroid)
            end if
            fc = fun(xc, context)
            if (fc < min(fr, f(ihi))) then
               simplex(:,ihi) = xc; f(ihi) = fc
            else
               do j = 1, m
                  if (j /= ilo) then
                     simplex(:,j) = simplex(:,ilo) + 0.5_dp*(simplex(:,j)-simplex(:,ilo))
                     f(j) = fun(simplex(:,j), context)
                  end if
               end do
            end if
         end if
      end do
      call order_indices(f, ilo, ihi, inhi)
      x = simplex(:,ilo); fbest = f(ilo)
      if (present(iterations)) iterations = min(it, max_iter)
      if (present(status)) then
         if (it <= max_iter) then
            status = 0
         else
            status = 1
         end if
      end if
   end subroutine nelder_mead

   subroutine order_indices(f, ilo, ihi, inhi)
      real(dp), intent(in) :: f(:)
      integer, intent(out) :: ilo, ihi, inhi
      integer :: i
      ilo = minloc(f, dim=1)
      ihi = maxloc(f, dim=1)
      inhi = ilo
      do i = 1, size(f)
         if (i == ihi) cycle
         if (inhi == ihi .or. f(i) > f(inhi)) inhi = i
      end do
   end subroutine order_indices

   subroutine golden_maximize(fun, context, a, b, xmax, fmax, tol, maxit, iterations)
      procedure(objective_1d) :: fun
      class(*), intent(in) :: context
      real(dp), intent(in) :: a, b
      real(dp), intent(out) :: xmax, fmax
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxit
      integer, intent(out), optional :: iterations
      real(dp) :: left, right, c, d, fc, fd, eps, gr
      integer :: it, mit
      gr = (sqrt(5.0_dp)-1.0_dp)/2.0_dp
      eps = 1.0e-7_dp
      if (present(tol)) eps = tol
      mit = 10000
      if (present(maxit)) mit = maxit
      left = a; right = b
      c = right - gr*(right-left)
      d = left + gr*(right-left)
      fc = fun(c, context); fd = fun(d, context)
      do it = 1, mit
         if (abs(right-left) <= eps*(1.0_dp+abs(c)+abs(d))) exit
         if (fc > fd) then
            right = d; d = c; fd = fc
            c = right - gr*(right-left); fc = fun(c, context)
         else
            left = c; c = d; fc = fd
            d = left + gr*(right-left); fd = fun(d, context)
         end if
      end do
      if (fc > fd) then
         xmax = c; fmax = fc
      else
         xmax = d; fmax = fd
      end if
      if (present(iterations)) iterations = min(it, mit)
   end subroutine golden_maximize
end module ccd_optimize
