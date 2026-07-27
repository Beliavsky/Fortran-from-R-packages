! SPDX-License-Identifier: GPL-3.0-only
! Numerical translation of fHMM 1.4.3.
module fhmm_math
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use fhmm_kinds, only: dp, pi, tiny_prob
   implicit none
   private
   public :: log_sum_exp, stationary_distribution, solve_linear, invert_matrix
   public :: normal_cdf, normal_quantile, regularized_gamma_p, regularized_beta
   public :: seed_rng, random_normal, random_gamma, random_poisson
   public :: numerical_gradient, numerical_hessian, symmetrize

contains

   pure real(dp) function log_sum_exp(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: xmax
      xmax = maxval(x)
      if (.not. ieee_is_finite(xmax)) then
         value = xmax
      else
         value = xmax + log(max(sum(exp(x - xmax)), tiny_prob))
      end if
   end function log_sum_exp

   subroutine solve_linear(a, b, x, ok)
      real(dp), intent(in) :: a(:, :), b(:)
      real(dp), intent(out) :: x(:)
      logical, intent(out) :: ok
      real(dp), allocatable :: aug(:, :), row(:)
      real(dp) :: pivot, factor
      integer :: n, i, k, p
      n = size(b)
      ok = .false.
      if (size(a,1) /= n .or. size(a,2) /= n .or. size(x) /= n) return
      allocate(aug(n,n+1), row(n+1))
      aug(:,1:n) = a
      aug(:,n+1) = b
      do k = 1, n
         p = k - 1 + maxloc(abs(aug(k:n,k)), dim=1)
         if (abs(aug(p,k)) <= 100.0_dp*epsilon(1.0_dp)) return
         if (p /= k) then
            row = aug(k,:)
            aug(k,:) = aug(p,:)
            aug(p,:) = row
         end if
         pivot = aug(k,k)
         aug(k,k:n+1) = aug(k,k:n+1)/pivot
         do i = 1, n
            if (i == k) cycle
            factor = aug(i,k)
            aug(i,k:n+1) = aug(i,k:n+1) - factor*aug(k,k:n+1)
         end do
      end do
      x = aug(:,n+1)
      ok = all(ieee_is_finite(x))
   end subroutine solve_linear

   subroutine invert_matrix(a, ainv, ok)
      real(dp), intent(in) :: a(:, :)
      real(dp), intent(out) :: ainv(:, :)
      logical, intent(out) :: ok
      real(dp), allocatable :: e(:), x(:)
      integer :: n, j
      n = size(a,1)
      ok = .false.
      if (size(a,2) /= n .or. any(shape(ainv) /= [n,n])) return
      allocate(e(n),x(n))
      do j = 1, n
         e = 0.0_dp
         e(j) = 1.0_dp
         call solve_linear(a,e,x,ok)
         if (.not. ok) return
         ainv(:,j) = x
      end do
      call symmetrize(ainv)
      ok = .true.
   end subroutine invert_matrix

   function stationary_distribution(gamma) result(delta)
      real(dp), intent(in) :: gamma(:, :)
      real(dp), allocatable :: delta(:)
      real(dp), allocatable :: a(:, :), b(:)
      logical :: ok
      integer :: n
      n = size(gamma,1)
      allocate(delta(n),a(n,n),b(n))
      if (size(gamma,2) /= n) then
         delta = 1.0_dp/real(n,dp)
         return
      end if
      a = transpose(gamma) - identity_matrix(n)
      a(n,:) = 1.0_dp
      b = 0.0_dp
      b(n) = 1.0_dp
      call solve_linear(a,b,delta,ok)
      if (.not. ok .or. any(delta < 0.0_dp)) delta = 1.0_dp/real(n,dp)
      delta = max(delta,0.0_dp)
      if (sum(delta) <= tiny_prob) delta = 1.0_dp/real(n,dp)
      delta = delta/sum(delta)
   end function stationary_distribution

   pure function identity_matrix(n) result(a)
      integer, intent(in) :: n
      real(dp) :: a(n,n)
      integer :: i
      a = 0.0_dp
      do i=1,n
         a(i,i)=1.0_dp
      end do
   end function identity_matrix

   pure real(dp) function normal_cdf(x) result(p)
      real(dp), intent(in) :: x
      p = 0.5_dp*erfc(-x/sqrt(2.0_dp))
   end function normal_cdf

   pure real(dp) function normal_quantile(p) result(x)
      real(dp), intent(in) :: p
      real(dp), parameter :: a(6) = [ &
         -3.969683028665376e1_dp, 2.209460984245205e2_dp, &
         -2.759285104469687e2_dp, 1.383577518672690e2_dp, &
         -3.066479806614716e1_dp, 2.506628277459239_dp ]
      real(dp), parameter :: b(5) = [ &
         -5.447609879822406e1_dp, 1.615858368580409e2_dp, &
         -1.556989798598866e2_dp, 6.680131188771972e1_dp, &
         -1.328068155288572e1_dp ]
      real(dp), parameter :: c(6) = [ &
         -7.784894002430293e-3_dp, -3.223964580411365e-1_dp, &
         -2.400758277161838_dp, -2.549732539343734_dp, &
          4.374664141464968_dp, 2.938163982698783_dp ]
      real(dp), parameter :: d(4) = [ &
          7.784695709041462e-3_dp, 3.224671290700398e-1_dp, &
          2.445134137142996_dp, 3.754408661907416_dp ]
      real(dp), parameter :: plow=0.02425_dp, phigh=1.0_dp-plow
      real(dp) :: q, r
      if (p <= 0.0_dp) then
         x = -huge(1.0_dp)
      else if (p >= 1.0_dp) then
         x = huge(1.0_dp)
      else if (p < plow) then
         q = sqrt(-2.0_dp*log(p))
         x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
             ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
      else if (p <= phigh) then
         q = p-0.5_dp
         r = q*q
         x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q/ &
             (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
      else
         q = sqrt(-2.0_dp*log(1.0_dp-p))
         x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
              ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
      end if
   end function normal_quantile

   pure real(dp) function regularized_gamma_p(a,x) result(p)
      real(dp), intent(in) :: a,x
      integer, parameter :: itmax=500
      real(dp), parameter :: eps=3.0e-14_dp, fpmin=1.0e-300_dp
      real(dp) :: sumv, del, ap, b, c, d, h, an
      integer :: n
      if (a <= 0.0_dp .or. x <= 0.0_dp) then
         p = merge(0.0_dp, 0.0_dp, x <= 0.0_dp)
         return
      end if
      if (x < a+1.0_dp) then
         ap=a; sumv=1.0_dp/a; del=sumv
         do n=1,itmax
            ap=ap+1.0_dp
            del=del*x/ap
            sumv=sumv+del
            if (abs(del) < abs(sumv)*eps) exit
         end do
         p=sumv*exp(-x+a*log(x)-log_gamma(a))
      else
         b=x+1.0_dp-a; c=1.0_dp/fpmin; d=1.0_dp/b; h=d
         do n=1,itmax
            an=-real(n,dp)*(real(n,dp)-a)
            b=b+2.0_dp
            d=an*d+b
            if (abs(d)<fpmin) d=fpmin
            c=b+an/c
            if (abs(c)<fpmin) c=fpmin
            d=1.0_dp/d
            del=d*c
            h=h*del
            if (abs(del-1.0_dp)<eps) exit
         end do
         p=1.0_dp-exp(-x+a*log(x)-log_gamma(a))*h
      end if
      p=min(max(p,0.0_dp),1.0_dp)
   end function regularized_gamma_p

   pure real(dp) function beta_cont_frac(a,b,x) result(h)
      real(dp), intent(in) :: a,b,x
      integer, parameter :: maxit=500
      real(dp), parameter :: eps=3.0e-14_dp, fpmin=1.0e-300_dp
      integer :: m,m2
      real(dp) :: aa,c,d,del,qab,qam,qap
      qab=a+b; qap=a+1.0_dp; qam=a-1.0_dp
      c=1.0_dp
      d=1.0_dp-qab*x/qap
      if (abs(d)<fpmin) d=fpmin
      d=1.0_dp/d; h=d
      do m=1,maxit
         m2=2*m
         aa=real(m,dp)*(b-real(m,dp))*x/((qam+real(m2,dp))*(a+real(m2,dp)))
         d=1.0_dp+aa*d; if(abs(d)<fpmin)d=fpmin
         c=1.0_dp+aa/c; if(abs(c)<fpmin)c=fpmin
         d=1.0_dp/d; h=h*d*c
         aa=-(a+real(m,dp))*(qab+real(m,dp))*x/ &
            ((a+real(m2,dp))*(qap+real(m2,dp)))
         d=1.0_dp+aa*d; if(abs(d)<fpmin)d=fpmin
         c=1.0_dp+aa/c; if(abs(c)<fpmin)c=fpmin
         d=1.0_dp/d; del=d*c; h=h*del
         if(abs(del-1.0_dp)<eps) exit
      end do
   end function beta_cont_frac

   pure real(dp) function regularized_beta(x,a,b) result(p)
      real(dp), intent(in) :: x,a,b
      real(dp) :: bt
      if (x <= 0.0_dp) then
         p=0.0_dp
      else if (x >= 1.0_dp) then
         p=1.0_dp
      else
         bt=exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b)+a*log(x)+b*log(1.0_dp-x))
         if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
            p=bt*beta_cont_frac(a,b,x)/a
         else
            p=1.0_dp-bt*beta_cont_frac(b,a,1.0_dp-x)/b
         end if
      end if
      p=min(max(p,0.0_dp),1.0_dp)
   end function regularized_beta

   subroutine seed_rng(seed)
      integer, intent(in) :: seed
      integer :: n, i
      integer, allocatable :: put(:)
      call random_seed(size=n)
      allocate(put(n))
      do i=1,n
         put(i)=modulo(seed+104729*i,huge(1)-1)+1
      end do
      call random_seed(put=put)
   end subroutine seed_rng

   real(dp) function random_normal() result(z)
      real(dp) :: u1,u2
      call random_number(u1); call random_number(u2)
      u1=max(u1,tiny_prob)
      z=sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
   end function random_normal

   recursive real(dp) function random_gamma(shape,scale) result(x)
      real(dp), intent(in) :: shape,scale
      real(dp) :: d,c,z,u
      if (shape <= 0.0_dp .or. scale <= 0.0_dp) then
         x=0.0_dp; return
      end if
      if (shape < 1.0_dp) then
         call random_number(u)
         x=random_gamma(shape+1.0_dp,scale)*u**(1.0_dp/shape)
         return
      end if
      d=shape-1.0_dp/3.0_dp
      c=1.0_dp/sqrt(9.0_dp*d)
      do
         z=random_normal()
         if (1.0_dp+c*z <= 0.0_dp) cycle
         x=(1.0_dp+c*z)**3
         call random_number(u)
         if (u < 1.0_dp-0.0331_dp*z**4) exit
         if (log(u) < 0.5_dp*z*z+d*(1.0_dp-x+log(x))) exit
      end do
      x=scale*d*x
   end function random_gamma

   integer function random_poisson(lambda) result(k)
      real(dp), intent(in) :: lambda
      real(dp) :: l,p,u,z
      if (lambda <= 0.0_dp) then
         k=0
      else if (lambda < 30.0_dp) then
         l=exp(-lambda); p=1.0_dp; k=0
         do
            k=k+1; call random_number(u); p=p*u
            if (p <= l) exit
         end do
         k=k-1
      else
         z=random_normal()
         k=max(0,nint(lambda+sqrt(lambda)*z))
      end if
   end function random_poisson

   subroutine numerical_gradient(f,x,g)
      interface
         function f(z) result(value)
            import dp
            real(dp), intent(in) :: z(:)
            real(dp) :: value
         end function f
      end interface
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: g(:)
      real(dp), allocatable :: xp(:),xm(:)
      real(dp) :: h
      integer :: i
      allocate(xp(size(x)),xm(size(x)))
      do i=1,size(x)
         h=epsilon(1.0_dp)**(1.0_dp/3.0_dp)*max(1.0_dp,abs(x(i)))
         xp=x; xm=x; xp(i)=xp(i)+h; xm(i)=xm(i)-h
         g(i)=(f(xp)-f(xm))/(2.0_dp*h)
      end do
   end subroutine numerical_gradient

   subroutine numerical_hessian(f,x,hess)
      interface
         function f(z) result(value)
            import dp
            real(dp), intent(in) :: z(:)
            real(dp) :: value
         end function f
      end interface
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: hess(:, :)
      real(dp), allocatable :: xpp(:),xpm(:),xmp(:),xmm(:),xp(:),xm(:)
      real(dp) :: hi,hj,f0
      integer :: i,j,n
      n=size(x); allocate(xpp(n),xpm(n),xmp(n),xmm(n),xp(n),xm(n)); f0=f(x)
      do i=1,n
         hi=epsilon(1.0_dp)**0.25_dp*max(1.0_dp,abs(x(i)))
         xp=x; xm=x; xp(i)=xp(i)+hi; xm(i)=xm(i)-hi
         hess(i,i)=(f(xp)-2.0_dp*f0+f(xm))/(hi*hi)
         do j=i+1,n
            hj=epsilon(1.0_dp)**0.25_dp*max(1.0_dp,abs(x(j)))
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

   subroutine symmetrize(a)
      real(dp), intent(inout) :: a(:, :)
      integer :: i,j
      do j=1,size(a,2)
         do i=j+1,size(a,1)
            a(i,j)=0.5_dp*(a(i,j)+a(j,i))
            a(j,i)=a(i,j)
         end do
      end do
   end subroutine symmetrize

end module fhmm_math
