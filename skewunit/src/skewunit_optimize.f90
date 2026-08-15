! Modern Fortran translation of R package skewunit.
! SPDX-License-Identifier: GPL-2.0-or-later
module skewunit_optimize
   use skewunit_kinds, only : dp, eps_dp
   implicit none
   private

   abstract interface
      function objective_1d(x, context) result(f)
         import dp
         real(dp), intent(in) :: x
         class(*), intent(in) :: context
         real(dp) :: f
      end function objective_1d
      function objective_nd(x, context) result(f)
         import dp
         real(dp), intent(in) :: x(:)
         class(*), intent(in) :: context
         real(dp) :: f
      end function objective_nd
   end interface

   public :: brent_minimize, nelder_mead

contains

   subroutine brent_minimize(fun, context, lower, upper, xmin, fmin, status, &
      iterations, tol, maxit)
      procedure(objective_1d) :: fun
      class(*), intent(in) :: context
      real(dp), intent(in) :: lower, upper
      real(dp), intent(out) :: xmin, fmin
      integer, intent(out) :: status, iterations
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxit
      real(dp), parameter :: cgold = 0.3819660112501051518_dp
      real(dp) :: a, b, d, e, etemp, fu, fv, fw, fx, p, q, r
      real(dp) :: tol1, tol2, u, v, w, x, xm, tolerance
      integer :: iter, limit

      tolerance = 1.0e-10_dp
      limit = 1000
      if (present(tol)) tolerance = max(tol,32.0_dp*eps_dp)
      if (present(maxit)) limit = max(1,maxit)

      a = min(lower,upper)
      b = max(lower,upper)
      x = a+cgold*(b-a)
      w = x
      v = x
      fx = fun(x,context)
      fw = fx
      fv = fx
      d = 0.0_dp
      e = 0.0_dp
      status = 1

      do iter = 1, limit
         xm = 0.5_dp*(a+b)
         tol1 = tolerance*abs(x)+sqrt(eps_dp)
         tol2 = 2.0_dp*tol1
         if (abs(x-xm) <= tol2-0.5_dp*(b-a)) then
            status = 0
            exit
         end if

         if (abs(e) > tol1) then
            r = (x-w)*(fx-fv)
            q = (x-v)*(fx-fw)
            p = (x-v)*q-(x-w)*r
            q = 2.0_dp*(q-r)
            if (q > 0.0_dp) p = -p
            q = abs(q)
            etemp = e
            e = d
            if (abs(p) >= abs(0.5_dp*q*etemp) .or. &
                p <= q*(a-x) .or. p >= q*(b-x)) then
               if (x >= xm) then
                  e = a-x
               else
                  e = b-x
               end if
               d = cgold*e
            else
               d = p/q
               u = x+d
               if (u-a < tol2 .or. b-u < tol2) d = sign(tol1,xm-x)
            end if
         else
            if (x >= xm) then
               e = a-x
            else
               e = b-x
            end if
            d = cgold*e
         end if

         if (abs(d) >= tol1) then
            u = x+d
         else
            u = x+sign(tol1,d)
         end if
         fu = fun(u,context)

         if (fu <= fx) then
            if (u >= x) then
               a = x
            else
               b = x
            end if
            v = w
            fv = fw
            w = x
            fw = fx
            x = u
            fx = fu
         else
            if (u < x) then
               a = u
            else
               b = u
            end if
            if (fu <= fw .or. w == x) then
               v = w
               fv = fw
               w = u
               fw = fu
            else if (fu <= fv .or. v == x .or. v == w) then
               v = u
               fv = fu
            end if
         end if
      end do

      iterations = min(iter,limit)
      xmin = x
      fmin = fx
   end subroutine brent_minimize

   subroutine nelder_mead(fun, context, x0, xmin, fmin, status, iterations, &
      tol, maxit)
      procedure(objective_nd) :: fun
      class(*), intent(in) :: context
      real(dp), intent(in) :: x0(:)
      real(dp), intent(out) :: xmin(size(x0)), fmin
      integer, intent(out) :: status, iterations
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxit
      integer :: n, m, i, j, iter, limit
      real(dp) :: tolerance, fr, fe, fc, fbest, spread, xspread
      real(dp), allocatable :: simplex(:,:), fval(:), centroid(:), xr(:), xe(:), xc(:)

      n = size(x0)
      m = n+1
      tolerance = 1.0e-9_dp
      limit = 10000
      if (present(tol)) tolerance = max(tol,64.0_dp*eps_dp)
      if (present(maxit)) limit = max(1,maxit)

      allocate(simplex(n,m),fval(m),centroid(n),xr(n),xe(n),xc(n))
      simplex(:,1) = x0
      do j = 2, m
         simplex(:,j) = x0
         i = j-1
         if (x0(i) /= 0.0_dp) then
            simplex(i,j) = x0(i)*1.05_dp
         else
            simplex(i,j) = 0.00025_dp
         end if
      end do
      do j = 1, m
         fval(j) = fun(simplex(:,j),context)
      end do

      status = 1
      do iter = 1, limit
         call sort_simplex(simplex,fval)
         fbest = fval(1)
         spread = maxval(abs(fval-fbest))
         xspread = 0.0_dp
         do j = 2, m
            xspread = max(xspread,maxval(abs(simplex(:,j)-simplex(:,1))))
         end do
         if (spread <= tolerance*(1.0_dp+abs(fbest)) .and. &
             xspread <= sqrt(tolerance)*(1.0_dp+maxval(abs(simplex(:,1))))) then
            status = 0
            exit
         end if

         centroid = sum(simplex(:,1:n),dim=2)/real(n,dp)
         xr = centroid+(centroid-simplex(:,m))
         fr = fun(xr,context)

         if (fr < fval(1)) then
            xe = centroid+2.0_dp*(xr-centroid)
            fe = fun(xe,context)
            if (fe < fr) then
               simplex(:,m) = xe
               fval(m) = fe
            else
               simplex(:,m) = xr
               fval(m) = fr
            end if
         else if (fr < fval(n)) then
            simplex(:,m) = xr
            fval(m) = fr
         else
            if (fr < fval(m)) then
               xc = centroid+0.5_dp*(xr-centroid)
               fc = fun(xc,context)
               if (fc <= fr) then
                  simplex(:,m) = xc
                  fval(m) = fc
               else
                  call shrink_simplex(fun,context,simplex,fval)
               end if
            else
               xc = centroid+0.5_dp*(simplex(:,m)-centroid)
               fc = fun(xc,context)
               if (fc < fval(m)) then
                  simplex(:,m) = xc
                  fval(m) = fc
               else
                  call shrink_simplex(fun,context,simplex,fval)
               end if
            end if
         end if
      end do

      call sort_simplex(simplex,fval)
      xmin = simplex(:,1)
      fmin = fval(1)
      iterations = min(iter,limit)
   end subroutine nelder_mead

   subroutine sort_simplex(simplex, fval)
      real(dp), intent(inout) :: simplex(:,:), fval(:)
      integer :: i, j, k
      real(dp) :: tf
      real(dp), allocatable :: tx(:)
      allocate(tx(size(simplex,1)))
      do i = 1, size(fval)-1
         k = i
         do j = i+1, size(fval)
            if (fval(j) < fval(k)) k = j
         end do
         if (k /= i) then
            tf = fval(i)
            fval(i) = fval(k)
            fval(k) = tf
            tx = simplex(:,i)
            simplex(:,i) = simplex(:,k)
            simplex(:,k) = tx
         end if
      end do
   end subroutine sort_simplex

   subroutine shrink_simplex(fun, context, simplex, fval)
      procedure(objective_nd) :: fun
      class(*), intent(in) :: context
      real(dp), intent(inout) :: simplex(:,:), fval(:)
      integer :: j
      do j = 2, size(fval)
         simplex(:,j) = simplex(:,1)+0.5_dp*(simplex(:,j)-simplex(:,1))
         fval(j) = fun(simplex(:,j),context)
      end do
   end subroutine shrink_simplex

end module skewunit_optimize
