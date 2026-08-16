module discrete_weibull_numerics
   use discrete_weibull_kinds, only : dp, i64
   implicit none
   private
   real(dp), parameter :: euler_gamma = 0.577215664901532860606512090082402431_dp
   public :: log1p_stable, expm1_stable, log1mexp
   public :: logistic, logit, harmonic_number, invert_2x2
   public :: nelder_mead_2d

   abstract interface
      function objective_2d(x, data) result(value)
         import dp
         real(dp), intent(in) :: x(2)
         class(*), intent(in) :: data
         real(dp) :: value
      end function objective_2d
   end interface

contains

   pure elemental real(dp) function log1p_stable(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: x2
      if (abs(x) < 1.0e-5_dp) then
         x2 = x*x
         y = x - 0.5_dp*x2 + x*x2/3.0_dp - x2*x2/4.0_dp + x2*x2*x/5.0_dp
      else
         y = log(1.0_dp+x)
      end if
   end function log1p_stable

   pure elemental real(dp) function expm1_stable(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: x2
      if (abs(x) < 1.0e-5_dp) then
         x2 = x*x
         y = x + 0.5_dp*x2 + x*x2/6.0_dp + x2*x2/24.0_dp + x2*x2*x/120.0_dp
      else
         y = exp(x)-1.0_dp
      end if
   end function expm1_stable

   pure elemental real(dp) function log1mexp(x) result(y)
      ! log(1-exp(x)) for x <= 0.
      real(dp), intent(in) :: x
      if (x >= 0.0_dp) then
         y = -huge(1.0_dp)
      else if (x < -log(2.0_dp)) then
         y = log1p_stable(-exp(x))
      else
         y = log(-expm1_stable(x))
      end if
   end function log1mexp

   pure elemental real(dp) function logistic(x) result(p)
      real(dp), intent(in) :: x
      if (x >= 0.0_dp) then
         p = 1.0_dp/(1.0_dp+exp(-x))
      else
         p = exp(x)/(1.0_dp+exp(x))
      end if
   end function logistic

   pure elemental real(dp) function logit(p) result(x)
      real(dp), intent(in) :: p
      x = log(p)-log1p_stable(-p)
   end function logit

   pure real(dp) function harmonic_number(n) result(h)
      integer(i64), intent(in) :: n
      integer(i64) :: k
      real(dp) :: rn, r2
      if (n <= 0_i64) then
         h = 0.0_dp
      else if (n <= 100000_i64) then
         h = 0.0_dp
         do k = 1_i64, n
            h = h + 1.0_dp/real(k,dp)
         end do
      else
         rn = real(n,dp)
         r2 = rn*rn
         h = log(rn)+euler_gamma+0.5_dp/rn-1.0_dp/(12.0_dp*r2) + &
             1.0_dp/(120.0_dp*r2*r2)-1.0_dp/(252.0_dp*r2*r2*rn*rn)
      end if
   end function harmonic_number

   subroutine invert_2x2(a, ainv, status)
      real(dp), intent(in) :: a(2,2)
      real(dp), intent(out) :: ainv(2,2)
      integer, intent(out) :: status
      real(dp) :: det, scale
      scale = max(1.0_dp,maxval(abs(a)))
      det = a(1,1)*a(2,2)-a(1,2)*a(2,1)
      if (abs(det) <= 100.0_dp*epsilon(1.0_dp)*scale*scale) then
         ainv = huge(1.0_dp)
         status = 1
      else
         ainv(1,1) =  a(2,2)/det
         ainv(1,2) = -a(1,2)/det
         ainv(2,1) = -a(2,1)/det
         ainv(2,2) =  a(1,1)/det
         status = 0
      end if
   end subroutine invert_2x2

   subroutine order3(f, ibest, imid, iworst)
      real(dp), intent(in) :: f(3)
      integer, intent(out) :: ibest, imid, iworst
      integer :: idx(3), i, j, tmp
      idx = [1,2,3]
      do i = 1, 2
         do j = i+1, 3
            if (f(idx(j)) < f(idx(i))) then
               tmp = idx(i)
               idx(i) = idx(j)
               idx(j) = tmp
            end if
         end do
      end do
      ibest = idx(1)
      imid = idx(2)
      iworst = idx(3)
   end subroutine order3

   subroutine nelder_mead_2d(fn, start, xbest, fbest, iterations, status, data, &
                             step, tol, max_iter)
      procedure(objective_2d) :: fn
      real(dp), intent(in) :: start(2)
      class(*), intent(in) :: data
      real(dp), intent(out) :: xbest(2), fbest
      integer, intent(out) :: iterations, status
      real(dp), intent(in), optional :: step, tol
      integer, intent(in), optional :: max_iter
      real(dp) :: x(2,3), f(3), centroid(2), xr(2), xe(2), xc(2)
      real(dp) :: fr, fe, fc, st, eps
      integer :: imax, ibest, imid, iworst, i

      st = 0.15_dp
      eps = 1.0e-9_dp
      imax = 1000
      if (present(step)) st = step
      if (present(tol)) eps = tol
      if (present(max_iter)) imax = max_iter

      x(:,1) = start
      x(:,2) = start + [st,0.0_dp]
      x(:,3) = start + [0.0_dp,st]
      do i = 1, 3
         f(i) = fn(x(:,i),data)
      end do

      status = 1
      do iterations = 1, imax
         call order3(f,ibest,imid,iworst)
         centroid = 0.5_dp*(x(:,ibest)+x(:,imid))
         xr = centroid + (centroid-x(:,iworst))
         fr = fn(xr,data)

         if (fr < f(ibest)) then
            xe = centroid + 2.0_dp*(xr-centroid)
            fe = fn(xe,data)
            if (fe < fr) then
               x(:,iworst) = xe
               f(iworst) = fe
            else
               x(:,iworst) = xr
               f(iworst) = fr
            end if
         else if (fr < f(imid)) then
            x(:,iworst) = xr
            f(iworst) = fr
         else
            if (fr < f(iworst)) then
               xc = centroid + 0.5_dp*(xr-centroid)
            else
               xc = centroid + 0.5_dp*(x(:,iworst)-centroid)
            end if
            fc = fn(xc,data)
            if (fc < min(fr,f(iworst))) then
               x(:,iworst) = xc
               f(iworst) = fc
            else
               do i = 1, 3
                  if (i /= ibest) then
                     x(:,i) = x(:,ibest)+0.5_dp*(x(:,i)-x(:,ibest))
                     f(i) = fn(x(:,i),data)
                  end if
               end do
            end if
         end if

         if (maxval(abs(x(:,1)-x(:,2))) <= eps*(1.0_dp+maxval(abs(x))) .and. &
             maxval(abs(x(:,1)-x(:,3))) <= eps*(1.0_dp+maxval(abs(x)))) then
            status = 0
            exit
         end if
      end do

      call order3(f,ibest,imid,iworst)
      xbest = x(:,ibest)
      fbest = f(ibest)
   end subroutine nelder_mead_2d

end module discrete_weibull_numerics
