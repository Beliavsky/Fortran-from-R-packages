module compoissonreg_numerics
   use compoissonreg_kinds, only : dp, pi_dp
   implicit none
   private
   public :: logadd, logistic, log1mexp, qnorm_std, pchisq_upper
   public :: invert_matrix, numerical_gradient, numerical_hessian

   abstract interface
      function scalar_fun(x) result(f)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp) :: f
      end function scalar_fun
   end interface

contains

   pure real(dp) function logadd(a, b)
      real(dp), intent(in) :: a, b
      real(dp) :: m
      if (a <= -huge(1.0_dp)/2) then
         logadd = b
      else if (b <= -huge(1.0_dp)/2) then
         logadd = a
      else
         m = max(a, b)
         logadd = m + log(exp(a-m) + exp(b-m))
      end if
   end function logadd

   pure real(dp) function logistic(x)
      real(dp), intent(in) :: x
      if (x >= 0.0_dp) then
         logistic = 1.0_dp / (1.0_dp + exp(-x))
      else
         logistic = exp(x) / (1.0_dp + exp(x))
      end if
   end function logistic

   pure real(dp) function log1mexp(x)
      real(dp), intent(in) :: x
      if (x > 0.0_dp) then
         log1mexp = -huge(1.0_dp)
      else if (x < log(0.5_dp)) then
         log1mexp = log(1.0_dp-exp(x))
      else
         log1mexp = log(1.0_dp-exp(x))
      end if
   end function log1mexp

   pure real(dp) function qnorm_std(p)
      real(dp), intent(in) :: p
      real(dp), parameter :: a1=-3.969683028665376e+01_dp, a2=2.209460984245205e+02_dp
      real(dp), parameter :: a3=-2.759285104469687e+02_dp, a4=1.383577518672690e+02_dp
      real(dp), parameter :: a5=-3.066479806614716e+01_dp, a6=2.506628277459239e+00_dp
      real(dp), parameter :: b1=-5.447609879822406e+01_dp, b2=1.615858368580409e+02_dp
      real(dp), parameter :: b3=-1.556989798598866e+02_dp, b4=6.680131188771972e+01_dp
      real(dp), parameter :: b5=-1.328068155288572e+01_dp
      real(dp), parameter :: c1=-7.784894002430293e-03_dp, c2=-3.223964580411365e-01_dp
      real(dp), parameter :: c3=-2.400758277161838e+00_dp, c4=-2.549732539343734e+00_dp
      real(dp), parameter :: c5=4.374664141464968e+00_dp, c6=2.938163982698783e+00_dp
      real(dp), parameter :: d1=7.784695709041462e-03_dp, d2=3.224671290700398e-01_dp
      real(dp), parameter :: d3=2.445134137142996e+00_dp, d4=3.754408661907416e+00_dp
      real(dp), parameter :: plow=0.02425_dp, phigh=1.0_dp-plow
      real(dp) :: q, r
      if (p <= 0.0_dp) then
         qnorm_std = -huge(1.0_dp)
      else if (p >= 1.0_dp) then
         qnorm_std = huge(1.0_dp)
      else if (p < plow) then
         q = sqrt(-2.0_dp*log(p))
         qnorm_std = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / &
            ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
      else if (p > phigh) then
         q = sqrt(-2.0_dp*log(1.0_dp-p))
         qnorm_std = -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / &
            ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
      else
         q = p - 0.5_dp
         r = q*q
         qnorm_std = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q / &
            (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.0_dp)
      end if
   end function qnorm_std

   pure real(dp) function gammp(a, x)
      real(dp), intent(in) :: a, x
      integer :: n
      real(dp) :: ap, del, sumv, b, c, d, h, an
      real(dp), parameter :: eps=3.0e-14_dp, fpmin=1.0e-300_dp
      if (x <= 0.0_dp) then
         gammp = 0.0_dp
         return
      end if
      if (x < a + 1.0_dp) then
         ap = a
         sumv = 1.0_dp/a
         del = sumv
         do n=1,10000
            ap = ap + 1.0_dp
            del = del*x/ap
            sumv = sumv + del
            if (abs(del) <= abs(sumv)*eps) exit
         end do
         gammp = sumv*exp(-x + a*log(x) - log_gamma(a))
      else
         b = x + 1.0_dp - a
         c = 1.0_dp/fpmin
         d = 1.0_dp/b
         h = d
         do n=1,10000
            an = -real(n,dp)*(real(n,dp)-a)
            b = b + 2.0_dp
            d = an*d + b
            if (abs(d) < fpmin) d = fpmin
            c = b + an/c
            if (abs(c) < fpmin) c = fpmin
            d = 1.0_dp/d
            del = d*c
            h = h*del
            if (abs(del-1.0_dp) <= eps) exit
         end do
         gammp = 1.0_dp - exp(-x + a*log(x) - log_gamma(a))*h
      end if
      gammp = min(1.0_dp,max(0.0_dp,gammp))
   end function gammp

   pure real(dp) function pchisq_upper(x, df)
      real(dp), intent(in) :: x
      integer, intent(in) :: df
      if (x <= 0.0_dp) then
         pchisq_upper = 1.0_dp
      else if (df <= 0) then
         pchisq_upper = 0.0_dp
      else
         pchisq_upper = 1.0_dp - gammp(0.5_dp*real(df,dp), 0.5_dp*x)
      end if
   end function pchisq_upper

   subroutine invert_matrix(a, ainv, ok)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: ainv(:,:)
      logical, intent(out) :: ok
      real(dp), allocatable :: aug(:,:), tmp(:)
      real(dp) :: pivot, fac
      integer :: n, i, j, k, imax
      n = size(a,1)
      ok = .false.
      if (size(a,2) /= n .or. size(ainv,1) /= n .or. size(ainv,2) /= n) return
      allocate(aug(n,2*n), tmp(2*n))
      aug(:,1:n) = a
      aug(:,n+1:2*n) = 0.0_dp
      do i=1,n
         aug(i,n+i) = 1.0_dp
      end do
      do i=1,n
         imax = i
         do k=i+1,n
            if (abs(aug(k,i)) > abs(aug(imax,i))) imax = k
         end do
         if (abs(aug(imax,i)) < 1.0e-14_dp) return
         if (imax /= i) then
            tmp = aug(i,:); aug(i,:) = aug(imax,:); aug(imax,:) = tmp
         end if
         pivot = aug(i,i)
         aug(i,:) = aug(i,:)/pivot
         do j=1,n
            if (j == i) cycle
            fac = aug(j,i)
            aug(j,:) = aug(j,:) - fac*aug(i,:)
         end do
      end do
      ainv = aug(:,n+1:2*n)
      ok = .true.
   end subroutine invert_matrix

   subroutine numerical_gradient(f, x, g, h)
      procedure(scalar_fun) :: f
      real(dp), intent(in) :: x(:), h
      real(dp), intent(out) :: g(:)
      real(dp), allocatable :: xp(:), xm(:)
      real(dp) :: step
      integer :: j
      allocate(xp(size(x)), xm(size(x)))
      do j=1,size(x)
         step = h*max(1.0_dp,abs(x(j)))
         xp=x; xm=x
         xp(j)=xp(j)+step; xm(j)=xm(j)-step
         g(j)=(f(xp)-f(xm))/(2.0_dp*step)
      end do
   end subroutine numerical_gradient

   subroutine numerical_hessian(f, x, hess, h)
      procedure(scalar_fun) :: f
      real(dp), intent(in) :: x(:), h
      real(dp), intent(out) :: hess(:,:)
      real(dp), allocatable :: xpp(:), xpm(:), xmp(:), xmm(:), xp(:), xm(:)
      real(dp) :: hi, hj, f0
      integer :: i,j,n
      n=size(x); allocate(xpp(n),xpm(n),xmp(n),xmm(n),xp(n),xm(n))
      f0=f(x)
      do i=1,n
         hi=h*max(1.0_dp,abs(x(i)))
         xp=x; xm=x; xp(i)=xp(i)+hi; xm(i)=xm(i)-hi
         hess(i,i)=(f(xp)-2.0_dp*f0+f(xm))/(hi*hi)
         do j=i+1,n
            hj=h*max(1.0_dp,abs(x(j)))
            xpp=x; xpm=x; xmp=x; xmm=x
            xpp(i)=xpp(i)+hi; xpp(j)=xpp(j)+hj
            xpm(i)=xpm(i)+hi; xpm(j)=xpm(j)-hj
            xmp(i)=xmp(i)-hi; xmp(j)=xmp(j)+hj
            xmm(i)=xmm(i)-hi; xmm(j)=xmm(j)-hj
            hess(i,j)=(f(xpp)-f(xpm)-f(xmp)+f(xmm))/(4.0_dp*hi*hj)
            hess(j,i)=hess(i,j)
         end do
      end do
   end subroutine numerical_hessian

end module compoissonreg_numerics
