module discrete_weibull_dist
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use discrete_weibull_kinds, only : dp, i64
   use discrete_weibull_numerics, only : expm1_stable, log1mexp, harmonic_number
   implicit none
   private

   public :: ddweibull, pdweibull, qdweibull, rdweibull
   public :: ddweibull3, pdweibull3, qdweibull3, rdweibull3, hdweibull3
   public :: log_ddweibull, log_ddweibull3
   public :: Edweibull, E2dweibull, Vdweibull, ERdweibull
   public :: Edweibull3, E2dweibull3

contains

   pure elemental logical function valid_type1(q,beta) result(ok)
      real(dp), intent(in) :: q,beta
      ok = q >= 0.0_dp .and. q < 1.0_dp .and. beta > 0.0_dp .and. &
           ieee_is_finite(q) .and. ieee_is_finite(beta)
   end function valid_type1

   pure elemental real(dp) function log_ddweibull(x,q,beta,zero) result(lp)
      integer(i64), intent(in) :: x
      real(dp), intent(in) :: q,beta
      logical, intent(in), optional :: zero
      logical :: z
      integer(i64) :: lo
      real(dp) :: a,b,lq,la,lb

      z = .false.
      if (present(zero)) z = zero
      lo = merge(0_i64,1_i64,z)

      if (.not. valid_type1(q,beta) .or. x < lo) then
         lp = -huge(1.0_dp)
         return
      end if
      if (q <= 0.0_dp) then
         if (x == lo) then
            lp = 0.0_dp
         else
            lp = -huge(1.0_dp)
         end if
         return
      end if

      lq = log(q)
      if (z) then
         a = real(x,dp)**beta
         b = real(x+1_i64,dp)**beta
      else
         a = real(x-1_i64,dp)**beta
         b = real(x,dp)**beta
      end if
      la = a*lq
      lb = b*lq
      lp = la + log1mexp(lb-la)
   end function log_ddweibull

   pure elemental real(dp) function ddweibull(x,q,beta,zero) result(p)
      integer(i64), intent(in) :: x
      real(dp), intent(in) :: q,beta
      logical, intent(in), optional :: zero
      real(dp) :: lp
      lp = log_ddweibull(x,q,beta,zero)
      if (lp <= -0.5_dp*huge(1.0_dp)) then
         p = 0.0_dp
      else
         p = exp(lp)
      end if
   end function ddweibull

   pure elemental real(dp) function pdweibull(x,q,beta,zero,lower_tail) result(p)
      real(dp), intent(in) :: x,q,beta
      logical, intent(in), optional :: zero,lower_tail
      logical :: z,lt
      integer(i64) :: k,lo
      real(dp) :: exponent,surv

      z = .false.
      lt = .true.
      if (present(zero)) z = zero
      if (present(lower_tail)) lt = lower_tail
      lo = merge(0_i64,1_i64,z)

      if (.not. valid_type1(q,beta)) then
         p = 0.0_dp
         return
      end if
      if (x < real(lo,dp)) then
         p = merge(0.0_dp,1.0_dp,lt)
         return
      end if
      if (q <= 0.0_dp) then
         p = merge(1.0_dp,0.0_dp,lt)
         return
      end if

      k = floor(x,kind=i64)
      if (z) then
         exponent = real(k+1_i64,dp)**beta
      else
         exponent = real(k,dp)**beta
      end if
      surv = exp(log(q)*exponent)
      p = merge(1.0_dp-surv,surv,lt)
   end function pdweibull

   pure elemental integer(i64) function qdweibull(p,q,beta,zero) result(x)
      real(dp), intent(in) :: p,q,beta
      logical, intent(in), optional :: zero
      logical :: z
      integer(i64) :: lo
      real(dp) :: r

      z = .false.
      if (present(zero)) z = zero
      lo = merge(0_i64,1_i64,z)

      if (.not. valid_type1(q,beta) .or. p < 0.0_dp .or. p > 1.0_dp .or. &
          .not. ieee_is_finite(p)) then
         x = -huge(1_i64)
         return
      end if
      if (p <= 0.0_dp) then
         x = lo
         return
      end if
      if (p >= 1.0_dp) then
         x = huge(1_i64)
         return
      end if
      if (q <= 0.0_dp) then
         x = lo
         return
      end if

      r = (log(1.0_dp-p)/log(q))**(1.0_dp/beta)
      if (z) then
         x = max(lo,int(ceiling(r-1.0_dp),i64))
      else
         x = max(lo,int(ceiling(r),i64))
      end if
   end function qdweibull

   subroutine rdweibull(x,q,beta,zero)
      integer(i64), intent(out) :: x(:)
      real(dp), intent(in) :: q,beta
      logical, intent(in), optional :: zero
      integer :: i
      real(dp) :: u
      if (.not. valid_type1(q,beta)) error stop "rdweibull: invalid parameters"
      do i = 1, size(x)
         call random_number(u)
         x(i) = qdweibull(u,q,beta,zero)
      end do
   end subroutine rdweibull

   pure real(dp) function power_sum(n,beta) result(s)
      integer(i64), intent(in) :: n
      real(dp), intent(in) :: beta
      integer(i64) :: i
      if (n <= 0_i64) then
         s = 0.0_dp
      else if (abs(beta+1.0_dp) <= 10.0_dp*epsilon(1.0_dp)) then
         s = harmonic_number(n)
      else
         s = 0.0_dp
         do i = 1_i64, n
            s = s + real(i,dp)**beta
         end do
      end if
   end function power_sum

   pure elemental logical function valid_type3(c,beta) result(ok)
      real(dp), intent(in) :: c,beta
      ok = c > 0.0_dp .and. beta >= -1.0_dp .and. &
           ieee_is_finite(c) .and. ieee_is_finite(beta)
   end function valid_type3

   pure elemental real(dp) function log_ddweibull3(x,c,beta) result(lp)
      integer(i64), intent(in) :: x
      real(dp), intent(in) :: c,beta
      real(dp) :: h
      if (.not. valid_type3(c,beta) .or. x < 0_i64) then
         lp = -huge(1.0_dp)
         return
      end if
      if (x == 0_i64) then
         lp = log1mexp(-c)
      else
         h = power_sum(x,beta)
         lp = log1mexp(-c*real(x+1_i64,dp)**beta)-c*h
      end if
   end function log_ddweibull3

   pure elemental real(dp) function ddweibull3(x,c,beta) result(p)
      integer(i64), intent(in) :: x
      real(dp), intent(in) :: c,beta
      real(dp) :: lp
      lp = log_ddweibull3(x,c,beta)
      if (lp <= -0.5_dp*huge(1.0_dp)) then
         p = 0.0_dp
      else
         p = exp(lp)
      end if
   end function ddweibull3

   pure elemental real(dp) function pdweibull3(x,c,beta,lower_tail) result(p)
      real(dp), intent(in) :: x,c,beta
      logical, intent(in), optional :: lower_tail
      logical :: lt
      integer(i64) :: k
      real(dp) :: surv
      lt = .true.
      if (present(lower_tail)) lt = lower_tail

      if (.not. valid_type3(c,beta)) then
         p = 0.0_dp
         return
      end if
      if (x < 0.0_dp) then
         p = merge(0.0_dp,1.0_dp,lt)
         return
      end if

      k = floor(x,kind=i64)
      surv = exp(-c*power_sum(k+1_i64,beta))
      p = merge(1.0_dp-surv,surv,lt)
   end function pdweibull3

   pure elemental integer(i64) function qdweibull3(p,c,beta) result(x)
      real(dp), intent(in) :: p,c,beta
      real(dp) :: target,s
      integer(i64) :: n,lo,hi,mid
      if (.not. valid_type3(c,beta) .or. p < 0.0_dp .or. p > 1.0_dp .or. &
          .not. ieee_is_finite(p)) then
         x = -huge(1_i64)
         return
      end if
      if (p <= 0.0_dp) then
         x = 0_i64
         return
      end if
      if (p >= 1.0_dp) then
         x = huge(1_i64)
         return
      end if

      target = -log(1.0_dp-p)/c
      if (target <= 1.0_dp) then
         x = 0_i64
         return
      end if

      if (abs(beta+1.0_dp) <= 10.0_dp*epsilon(1.0_dp)) then
         ! Find smallest n=x+1 such that H_n >= target.
         if (target-euler_constant_local() > log(real(huge(1_i64),dp))) then
            x = huge(1_i64)
            return
         end if
         hi = max(2_i64,int(exp(min(target-euler_constant_local(), &
                  log(real(2305843009213693951_i64,dp)))),i64))
         do while (harmonic_number(hi) < target .and. hi < 4611686018427387903_i64)
            hi = 2_i64*hi
         end do
         lo = 1_i64
         do while (lo < hi)
            mid = lo+(hi-lo)/2_i64
            if (harmonic_number(mid) >= target) then
               hi = mid
            else
               lo = mid+1_i64
            end if
         end do
         x = lo-1_i64
         return
      end if

      s = 0.0_dp
      n = 0_i64
      do
         n = n+1_i64
         s = s+real(n,dp)**beta
         if (s >= target) exit
         if (n >= 100000000_i64) then
            x = huge(1_i64)
            return
         end if
      end do
      x = n-1_i64
   end function qdweibull3

   pure real(dp) function euler_constant_local() result(v)
      v = 0.577215664901532860606512090082402431_dp
   end function euler_constant_local

   subroutine rdweibull3(x,c,beta)
      integer(i64), intent(out) :: x(:)
      real(dp), intent(in) :: c,beta
      integer :: i
      real(dp) :: u
      if (.not. valid_type3(c,beta)) error stop "rdweibull3: invalid parameters"
      do i = 1, size(x)
         call random_number(u)
         x(i) = qdweibull3(u,c,beta)
      end do
   end subroutine rdweibull3

   pure elemental real(dp) function hdweibull3(x,c,beta) result(h)
      integer(i64), intent(in) :: x
      real(dp), intent(in) :: c,beta
      if (.not. valid_type3(c,beta) .or. x < 0_i64) then
         h = 0.0_dp
      else
         ! Correct mathematical hazard; upstream R accidentally omits "*c".
         h = -expm1_stable(-c*real(x+1_i64,dp)**beta)
      end if
   end function hdweibull3

   real(dp) function Edweibull(q,beta,eps,nmax,zero) result(e)
      real(dp), intent(in) :: q,beta
      real(dp), intent(in), optional :: eps
      integer, intent(in), optional :: nmax
      logical, intent(in), optional :: zero
      real(dp) :: tol,lambda
      integer :: mx
      integer(i64) :: xmax,x
      logical :: z
      tol = 1.0e-4_dp
      mx = 1000
      z = .false.
      if (present(eps)) tol = eps
      if (present(nmax)) mx = nmax
      if (present(zero)) z = zero

      if (.not. valid_type1(q,beta)) then
         e = huge(1.0_dp)
         return
      end if
      if (q <= 0.0_dp) then
         e = merge(0.0_dp,1.0_dp,z)
         return
      end if
      if (abs(beta-1.0_dp) <= 10.0_dp*epsilon(1.0_dp)) then
         if (z) then
            e = q/(1.0_dp-q)
         else
            e = 1.0_dp/(1.0_dp-q)
         end if
         return
      end if

      xmax = min(2_i64*qdweibull(1.0_dp-tol,q,beta,z),int(mx,i64))
      if (xmax < int(mx,i64)) then
         e = 0.0_dp
         do x = 1_i64, xmax
            e = e+ddweibull(x,q,beta,z)*real(x,dp)
         end do
      else
         lambda = (-1.0_dp/log(q))**(1.0_dp/beta)
         e = lambda*gamma(1.0_dp+1.0_dp/beta)+0.5_dp
         if (z) e = e-1.0_dp
      end if
   end function Edweibull

   real(dp) function E2dweibull(q,beta,eps,nmax,zero) result(e)
      real(dp), intent(in) :: q,beta
      real(dp), intent(in), optional :: eps
      integer, intent(in), optional :: nmax
      logical, intent(in), optional :: zero
      real(dp) :: tol,lambda,m
      integer :: mx
      integer(i64) :: xmax,x
      logical :: z
      tol = 1.0e-4_dp
      mx = 1000
      z = .false.
      if (present(eps)) tol = eps
      if (present(nmax)) mx = nmax
      if (present(zero)) z = zero

      if (.not. valid_type1(q,beta)) then
         e = huge(1.0_dp)
         return
      end if
      if (q <= 0.0_dp) then
         e = merge(0.0_dp,1.0_dp,z)
         return
      end if
      if (abs(beta-1.0_dp) <= 10.0_dp*epsilon(1.0_dp)) then
         if (z) then
            e = q*(1.0_dp+q)/(1.0_dp-q)**2
         else
            e = (1.0_dp+q)/(1.0_dp-q)**2
         end if
         return
      end if

      xmax = 2_i64*qdweibull(1.0_dp-tol,q,beta,z)
      if (xmax < int(mx,i64)) then
         e = 0.0_dp
         do x = 1_i64, xmax
            e = e+ddweibull(x,q,beta,z)*real(x*x,dp)
         end do
      else
         lambda = (-1.0_dp/log(q))**(1.0_dp/beta)
         m = Edweibull(q,beta,tol,mx,z)
         e = lambda**2*(gamma(1.0_dp+2.0_dp/beta)- &
             gamma(1.0_dp+1.0_dp/beta)**2)+m*m
         e = ceiling(e)
         if (z) e = e-1.0_dp
      end if
   end function E2dweibull

   real(dp) function Vdweibull(q,beta,eps,nmax,zero) result(v)
      real(dp), intent(in) :: q,beta
      real(dp), intent(in), optional :: eps
      integer, intent(in), optional :: nmax
      logical, intent(in), optional :: zero
      real(dp) :: m1,m2
      if (abs(beta-1.0_dp) <= 10.0_dp*epsilon(1.0_dp)) then
         v = q/(1.0_dp-q)**2
      else
         m1 = Edweibull(q,beta,eps,nmax,zero)
         m2 = E2dweibull(q,beta,eps,nmax,zero)
         v = m2-m1*m1
      end if
   end function Vdweibull

   real(dp) function ERdweibull(q,beta,eps,nmax) result(e)
      real(dp), intent(in) :: q,beta
      real(dp), intent(in), optional :: eps
      integer, intent(in), optional :: nmax
      real(dp) :: tol
      integer :: mx
      integer(i64) :: xmax,x
      tol = 1.0e-4_dp
      mx = 1000
      if (present(eps)) tol = eps
      if (present(nmax)) mx = nmax
      if (.not. valid_type1(q,beta)) then
         e = huge(1.0_dp)
      else if (q <= 0.0_dp) then
         e = 1.0_dp
      else if (abs(beta-1.0_dp) <= 10.0_dp*epsilon(1.0_dp)) then
         e = (1.0_dp-q)/q*log(1.0_dp/(1.0_dp-q))
      else
         xmax = min(2_i64*qdweibull(1.0_dp-tol,q,beta,.false.),int(mx,i64))
         e = 0.0_dp
         do x = 1_i64, xmax
            e = e+ddweibull(x,q,beta,.false.)/real(x,dp)
         end do
      end if
   end function ERdweibull

   real(dp) function Edweibull3(c,beta,eps) result(e)
      real(dp), intent(in) :: c,beta
      real(dp), intent(in), optional :: eps
      real(dp) :: tol,csum
      integer(i64) :: nmax,i
      tol = 1.0e-4_dp
      if (present(eps)) tol = eps
      if (.not. valid_type3(c,beta)) then
         e = huge(1.0_dp)
         return
      end if
      if (abs(beta) <= 10.0_dp*epsilon(1.0_dp)) then
         e = exp(-c)/(1.0_dp-exp(-c))
         return
      end if
      nmax = 2_i64*qdweibull3(1.0_dp-tol,c,beta)
      if (nmax >= 2305843009213693951_i64) then
         e = huge(1.0_dp)
         return
      end if
      e = 0.0_dp
      csum = 0.0_dp
      do i = 1_i64, nmax
         csum = csum+real(i,dp)**beta
         e = e+exp(-c*csum)
      end do
   end function Edweibull3

   real(dp) function E2dweibull3(c,beta,eps) result(e)
      real(dp), intent(in) :: c,beta
      real(dp), intent(in), optional :: eps
      real(dp) :: tol,csum,m1,s
      integer(i64) :: nmax,i
      tol = 1.0e-4_dp
      if (present(eps)) tol = eps
      if (.not. valid_type3(c,beta)) then
         e = huge(1.0_dp)
         return
      end if
      if (abs(beta) <= 10.0_dp*epsilon(1.0_dp)) then
         e = exp(-c)*(1.0_dp+exp(-c))/(1.0_dp-exp(-c))**2
         return
      end if
      nmax = 2_i64*qdweibull3(1.0_dp-tol,c,beta)
      if (nmax >= 2305843009213693951_i64) then
         e = huge(1.0_dp)
         return
      end if
      m1 = 0.0_dp
      s = 0.0_dp
      csum = 0.0_dp
      do i = 1_i64, nmax
         csum = csum+real(i,dp)**beta
         m1 = m1+exp(-c*csum)
         s = s+2.0_dp*real(i,dp)*exp(-c*csum)
      end do
      e = s-m1
   end function E2dweibull3

end module discrete_weibull_dist
