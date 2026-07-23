! SPDX-License-Identifier: GPL-3.0-only
!
! Experimental modern Fortran translation of computational methods from the
! R package rmgarch, copyright (C) 2008-2025 Alexios Galanos.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 3 only.
module rmgarch_distributions
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use rmgarch_kinds, only : dp
   use rmgarch_math, only : pi, normal_cdf, normal_quantile, cholesky_lower, &
      logdet_spd, quadratic_form_spd
   use rmgarch_rng, only : random_normal, random_gamma
   use rmgarch_types, only : dist_gaussian, dist_student, dist_laplace
   implicit none
   private

   public :: multivariate_normal_logpdf, multivariate_student_logpdf
   public :: multivariate_laplace_logpdf, multivariate_logpdf
   public :: random_multivariate_normal, random_multivariate_student
   public :: random_multivariate_laplace, random_multivariate
   public :: weighted_margin, weighted_margin_path
   public :: standardized_student_pdf, standardized_student_cdf
   public :: standardized_student_quantile
   public :: standardized_laplace_pdf, standardized_laplace_cdf
   public :: standardized_laplace_quantile
   public :: empirical_uniform_transform, normal_score_transform
   public :: student_score_transform, clamp_probabilities
   public :: modified_bessel_k_half_integer_order

contains

   function multivariate_normal_logpdf(x, mean, covariance, valid) result(value)
      real(dp), intent(in) :: x(:), mean(:), covariance(:,:)
      logical, intent(out), optional :: valid
      real(dp) :: value, centered(size(x)), ld, quad
      logical :: ok

      ok = size(mean) == size(x) .and. size(covariance,1) == size(x) .and. &
         size(covariance,2) == size(x)
      if (ok) then
         centered = x-mean
         ld = logdet_spd(covariance,ok)
      end if
      if (ok) quad = quadratic_form_spd(covariance,centered,ok)
      if (ok) then
         value = -0.5_dp*(real(size(x),dp)*log(2.0_dp*pi)+ld+quad)
         ok = ieee_is_finite(value)
      else
         value = -huge(1.0_dp)
      end if
      if (present(valid)) valid = ok
   end function multivariate_normal_logpdf

   function multivariate_student_logpdf(x, mean, covariance, nu, valid) result(value)
      real(dp), intent(in) :: x(:), mean(:), covariance(:,:), nu
      logical, intent(out), optional :: valid
      real(dp) :: value, centered(size(x)), ld, quad, m
      logical :: ok

      ok = nu > 2.0_dp .and. size(mean) == size(x) .and. &
         size(covariance,1) == size(x) .and. size(covariance,2) == size(x)
      if (ok) then
         centered = x-mean
         ld = logdet_spd(covariance,ok)
      end if
      if (ok) quad = quadratic_form_spd(covariance,centered,ok)
      m = real(size(x),dp)
      if (ok) then
         value = log_gamma(0.5_dp*(nu+m))-log_gamma(0.5_dp*nu)- &
            0.5_dp*m*log(pi*(nu-2.0_dp))-0.5_dp*ld- &
            0.5_dp*(nu+m)*log(1.0_dp+quad/(nu-2.0_dp))
         ok = ieee_is_finite(value)
      else
         value = -huge(1.0_dp)
      end if
      if (present(valid)) valid = ok
   end function multivariate_student_logpdf

   function multivariate_laplace_logpdf(x, mean, covariance, valid) result(value)
      real(dp), intent(in) :: x(:), mean(:), covariance(:,:)
      logical, intent(out), optional :: valid
      real(dp) :: value, centered(size(x)), ld, quad, order, kval, m
      logical :: ok

      ok = size(mean) == size(x) .and. size(covariance,1) == size(x) .and. &
         size(covariance,2) == size(x)
      if (ok) then
         centered = x-mean
         ld = logdet_spd(covariance,ok)
      end if
      if (ok) quad = quadratic_form_spd(covariance,centered,ok)
      m = real(size(x),dp)
      if (ok) then
         quad = max(quad,1.0e-14_dp)
         order = abs(0.5_dp*(2.0_dp-m))
         kval = modified_bessel_k_half_integer_order(order,sqrt(2.0_dp*quad),ok)
      end if
      if (ok .and. kval > 0.0_dp) then
         value = log(2.0_dp)-0.5_dp*m*log(2.0_dp*pi)-0.5_dp*ld+ &
            0.25_dp*(2.0_dp-m)*log(0.5_dp*quad)+log(kval)
         ok = ieee_is_finite(value)
      else
         value = -huge(1.0_dp)
         ok = .false.
      end if
      if (present(valid)) valid = ok
   end function multivariate_laplace_logpdf

   function multivariate_logpdf(distribution, x, mean, covariance, shape, valid) result(value)
      integer, intent(in) :: distribution
      real(dp), intent(in) :: x(:), mean(:), covariance(:,:)
      real(dp), intent(in), optional :: shape
      logical, intent(out), optional :: valid
      real(dp) :: value, nu
      logical :: ok

      nu = 8.0_dp
      if (present(shape)) nu = shape
      select case (distribution)
      case (dist_gaussian)
         value = multivariate_normal_logpdf(x,mean,covariance,ok)
      case (dist_student)
         value = multivariate_student_logpdf(x,mean,covariance,nu,ok)
      case (dist_laplace)
         value = multivariate_laplace_logpdf(x,mean,covariance,ok)
      case default
         value = -huge(1.0_dp)
         ok = .false.
      end select
      if (present(valid)) valid = ok
   end function multivariate_logpdf

   subroutine random_multivariate_normal(mean, covariance, x, valid)
      real(dp), intent(in) :: mean(:), covariance(:,:)
      real(dp), intent(out) :: x(size(mean))
      logical, intent(out), optional :: valid
      real(dp) :: lower(size(mean),size(mean)), z(size(mean))
      logical :: ok
      integer :: i

      ok = size(covariance,1) == size(mean) .and. size(covariance,2) == size(mean)
      if (ok) call cholesky_lower(covariance,lower,ok)
      if (ok) then
         do i = 1, size(mean)
            z(i) = random_normal()
         end do
         x = mean+matmul(lower,z)
      else
         x = mean
      end if
      if (present(valid)) valid = ok
   end subroutine random_multivariate_normal

   subroutine random_multivariate_student(mean, covariance, nu, x, valid)
      real(dp), intent(in) :: mean(:), covariance(:,:), nu
      real(dp), intent(out) :: x(size(mean))
      logical, intent(out), optional :: valid
      real(dp) :: centered(size(mean)), chi2
      logical :: ok

      ok = nu > 2.0_dp
      if (ok) call random_multivariate_normal(0.0_dp*mean,covariance,centered,ok)
      if (ok) then
         chi2 = 2.0_dp*random_gamma(0.5_dp*nu)
         x = mean+sqrt((nu-2.0_dp)/max(chi2,tiny(1.0_dp)))*centered
      else
         x = mean
      end if
      if (present(valid)) valid = ok
   end subroutine random_multivariate_student

   subroutine random_multivariate_laplace(mean, covariance, x, valid)
      real(dp), intent(in) :: mean(:), covariance(:,:)
      real(dp), intent(out) :: x(size(mean))
      logical, intent(out), optional :: valid
      real(dp) :: centered(size(mean)), u, exponential
      logical :: ok

      call random_multivariate_normal(0.0_dp*mean,covariance,centered,ok)
      if (ok) then
         call random_number(u)
         exponential = -log(max(u,tiny(1.0_dp)))
         x = mean+sqrt(exponential)*centered
      else
         x = mean
      end if
      if (present(valid)) valid = ok
   end subroutine random_multivariate_laplace

   subroutine random_multivariate(distribution, mean, covariance, x, shape, valid)
      integer, intent(in) :: distribution
      real(dp), intent(in) :: mean(:), covariance(:,:)
      real(dp), intent(out) :: x(size(mean))
      real(dp), intent(in), optional :: shape
      logical, intent(out), optional :: valid
      real(dp) :: nu
      logical :: ok

      nu = 8.0_dp
      if (present(shape)) nu = shape
      select case (distribution)
      case (dist_gaussian)
         call random_multivariate_normal(mean,covariance,x,ok)
      case (dist_student)
         call random_multivariate_student(mean,covariance,nu,x,ok)
      case (dist_laplace)
         call random_multivariate_laplace(mean,covariance,x,ok)
      case default
         x = mean
         ok = .false.
      end select
      if (present(valid)) valid = ok
   end subroutine random_multivariate

   subroutine weighted_margin(distribution, weights, mean, covariance, parameters, shape, valid)
      integer, intent(in) :: distribution
      real(dp), intent(in) :: weights(:), mean(:), covariance(:,:)
      real(dp), intent(out) :: parameters(4)
      real(dp), intent(in), optional :: shape
      logical, intent(out), optional :: valid
      real(dp) :: variance, nu
      logical :: ok

      ok = size(weights) == size(mean) .and. size(covariance,1) == size(weights) .and. &
         size(covariance,2) == size(weights)
      nu = 8.0_dp
      if (present(shape)) nu = shape
      parameters = 0.0_dp
      if (ok) then
         variance = dot_product(weights,matmul(covariance,weights))
         ok = variance >= 0.0_dp
      end if
      if (ok) then
         parameters(1) = dot_product(weights,mean)
         parameters(2) = sqrt(max(variance,0.0_dp))
         parameters(3) = 0.0_dp
         select case (distribution)
         case (dist_gaussian)
            parameters(4) = 0.0_dp
         case (dist_laplace)
            parameters(4) = 1.0_dp
         case (dist_student)
            parameters(4) = nu
            ok = nu > 2.0_dp
         case default
            ok = .false.
         end select
      end if
      if (present(valid)) valid = ok
   end subroutine weighted_margin

   subroutine weighted_margin_path(distribution, weights, mean, covariance, parameters, shape, valid)
      integer, intent(in) :: distribution
      real(dp), intent(in) :: weights(:,:), mean(:,:), covariance(:,:,:)
      real(dp), intent(out) :: parameters(size(covariance,3),4)
      real(dp), intent(in), optional :: shape
      logical, intent(out), optional :: valid
      integer :: t, iw, im
      real(dp) :: step_parameters(4)
      logical :: ok, step_ok

      ok = size(covariance,1) == size(covariance,2) .and. &
         size(weights,2) == size(covariance,1) .and. size(mean,2) == size(covariance,1)
      ok = ok .and. (size(weights,1) == 1 .or. size(weights,1) == size(covariance,3))
      ok = ok .and. (size(mean,1) == 1 .or. size(mean,1) == size(covariance,3))
      parameters = 0.0_dp
      if (ok) then
         do t = 1, size(covariance,3)
            iw = merge(1,t,size(weights,1) == 1)
            im = merge(1,t,size(mean,1) == 1)
            if (present(shape)) then
               call weighted_margin(distribution,weights(iw,:),mean(im,:),covariance(:,:,t), &
                  step_parameters,shape,step_ok)
            else
               call weighted_margin(distribution,weights(iw,:),mean(im,:),covariance(:,:,t), &
                  step_parameters,valid=step_ok)
            end if
            parameters(t,:) = step_parameters
            ok = ok .and. step_ok
         end do
      end if
      if (present(valid)) valid = ok
   end subroutine weighted_margin_path

   pure elemental function standardized_student_pdf(x, nu) result(value)
      real(dp), intent(in) :: x, nu
      real(dp) :: value
      if (nu <= 2.0_dp) then
         value = 0.0_dp
      else
         value = exp(log_gamma(0.5_dp*(nu+1.0_dp))-log_gamma(0.5_dp*nu)- &
            0.5_dp*log(pi*(nu-2.0_dp))-0.5_dp*(nu+1.0_dp)* &
            log(1.0_dp+x*x/(nu-2.0_dp)))
      end if
   end function standardized_student_pdf

   pure elemental function standardized_student_cdf(x, nu) result(value)
      real(dp), intent(in) :: x, nu
      real(dp) :: value, t, ib
      if (nu <= 2.0_dp) then
         value = 0.0_dp
      else if (abs(x) <= tiny(1.0_dp)) then
         value = 0.5_dp
      else
         t = x*sqrt(nu/(nu-2.0_dp))
         ib = regularized_beta(0.5_dp*nu,0.5_dp,nu/(nu+t*t))
         if (t > 0.0_dp) then
            value = 1.0_dp-0.5_dp*ib
         else
            value = 0.5_dp*ib
         end if
      end if
      value = min(1.0_dp,max(0.0_dp,value))
   end function standardized_student_cdf

   pure elemental function standardized_student_quantile(p, nu) result(value)
      real(dp), intent(in) :: p, nu
      real(dp) :: value, lo, hi, mid
      integer :: iter
      if (p <= 0.0_dp) then
         value = -huge(1.0_dp)
      else if (p >= 1.0_dp) then
         value = huge(1.0_dp)
      else if (nu <= 2.0_dp) then
         value = 0.0_dp
      else
         lo = -50.0_dp
         hi = 50.0_dp
         do iter = 1, 120
            mid = 0.5_dp*(lo+hi)
            if (standardized_student_cdf(mid,nu) < p) then
               lo = mid
            else
               hi = mid
            end if
         end do
         value = 0.5_dp*(lo+hi)
      end if
   end function standardized_student_quantile

   pure elemental function standardized_laplace_pdf(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: value
      value = exp(-sqrt(2.0_dp)*abs(x))/sqrt(2.0_dp)
   end function standardized_laplace_pdf

   pure elemental function standardized_laplace_cdf(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: value
      if (x < 0.0_dp) then
         value = 0.5_dp*exp(sqrt(2.0_dp)*x)
      else
         value = 1.0_dp-0.5_dp*exp(-sqrt(2.0_dp)*x)
      end if
   end function standardized_laplace_cdf

   pure elemental function standardized_laplace_quantile(p) result(value)
      real(dp), intent(in) :: p
      real(dp) :: value
      if (p <= 0.0_dp) then
         value = -huge(1.0_dp)
      else if (p >= 1.0_dp) then
         value = huge(1.0_dp)
      else if (p < 0.5_dp) then
         value = log(2.0_dp*p)/sqrt(2.0_dp)
      else
         value = -log(2.0_dp*(1.0_dp-p))/sqrt(2.0_dp)
      end if
   end function standardized_laplace_quantile

   pure elemental function clamp_probabilities(p) result(value)
      real(dp), intent(in) :: p
      real(dp) :: value
      value = min(1.0_dp-1.0e-10_dp,max(1.0e-10_dp,p))
   end function clamp_probabilities

   subroutine empirical_uniform_transform(x, u)
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(out) :: u(size(x,1),size(x,2))
      real(dp) :: sorted(size(x,1))
      integer :: n, m, i, j, less, equal

      n = size(x,1)
      m = size(x,2)
      if (n == 0) return
      do j = 1, m
         sorted = x(:,j)
         call sort_real(sorted)
         do i = 1, n
            less = lower_bound(sorted,x(i,j))-1
            equal = upper_bound(sorted,x(i,j))-lower_bound(sorted,x(i,j))
            u(i,j) = (real(less,dp)+0.5_dp*real(max(equal,1),dp)+0.5_dp)/real(n+1,dp)
            u(i,j) = clamp_probabilities(u(i,j))
         end do
      end do
   end subroutine empirical_uniform_transform

   subroutine normal_score_transform(u, z)
      real(dp), intent(in) :: u(:,:)
      real(dp), intent(out) :: z(size(u,1),size(u,2))
      z = normal_quantile(clamp_probabilities(u))
   end subroutine normal_score_transform

   subroutine student_score_transform(u, nu, z)
      real(dp), intent(in) :: u(:,:), nu
      real(dp), intent(out) :: z(size(u,1),size(u,2))
      z = standardized_student_quantile(clamp_probabilities(u),nu)
   end subroutine student_score_transform

   function modified_bessel_k_half_integer_order(order, x, valid) result(value)
      real(dp), intent(in) :: order, x
      logical, intent(out), optional :: valid
      real(dp) :: value, nearest_integer, nearest_half, km1, k0, kp1
      integer :: n, j
      logical :: ok

      ok = x > 0.0_dp .and. order >= 0.0_dp
      nearest_integer = real(nint(order),dp)
      nearest_half = real(nint(order-0.5_dp),dp)+0.5_dp
      if (.not. ok) then
         value = huge(1.0_dp)
      else if (abs(order-nearest_integer) < 1.0e-10_dp) then
         n = nint(order)
         if (n == 0) then
            value = bessel_k0_approx(x)
         else if (n == 1) then
            value = bessel_k1_approx(x)
         else
            km1 = bessel_k0_approx(x)
            k0 = bessel_k1_approx(x)
            do j = 1, n-1
               kp1 = km1+2.0_dp*real(j,dp)*k0/x
               km1 = k0
               k0 = kp1
            end do
            value = k0
         end if
      else if (abs(order-nearest_half) < 1.0e-10_dp) then
         n = nint(order-0.5_dp)
         value = half_integer_bessel_k(n,x)
      else
         value = bessel_k_integral(order,x)
      end if
      ok = ok .and. value > 0.0_dp .and. ieee_is_finite(value)
      if (present(valid)) valid = ok
   end function modified_bessel_k_half_integer_order

   pure function half_integer_bessel_k(n, x) result(value)
      integer, intent(in) :: n
      real(dp), intent(in) :: x
      real(dp) :: value, term, series
      integer :: k

      term = 1.0_dp
      series = 1.0_dp
      do k = 1, n
         term = term*real((n+k)*(n-k+1),dp)/(real(k,dp)*2.0_dp*x)
         series = series+term
      end do
      value = sqrt(pi/(2.0_dp*x))*exp(-x)*series
   end function half_integer_bessel_k

   pure function bessel_i0_approx(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: value, ax, y
      ax = abs(x)
      if (ax < 3.75_dp) then
         y = (x/3.75_dp)**2
         value = 1.0_dp+y*(3.5156229_dp+y*(3.0899424_dp+y*(1.2067492_dp+ &
            y*(0.2659732_dp+y*(0.0360768_dp+y*0.0045813_dp)))))
      else
         y = 3.75_dp/ax
         value = exp(ax)/sqrt(ax)*(0.39894228_dp+y*(0.01328592_dp+y*(0.00225319_dp+ &
            y*(-0.00157565_dp+y*(0.00916281_dp+y*(-0.02057706_dp+y*(0.02635537_dp+ &
            y*(-0.01647633_dp+y*0.00392377_dp))))))))
      end if
   end function bessel_i0_approx

   pure function bessel_i1_approx(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: value, ax, y
      ax = abs(x)
      if (ax < 3.75_dp) then
         y = (x/3.75_dp)**2
         value = ax*(0.5_dp+y*(0.87890594_dp+y*(0.51498869_dp+y*(0.15084934_dp+ &
            y*(0.02658733_dp+y*(0.00301532_dp+y*0.00032411_dp))))))
      else
         y = 3.75_dp/ax
         value = exp(ax)/sqrt(ax)*(0.39894228_dp+y*(-0.03988024_dp+y*(-0.00362018_dp+ &
            y*(0.00163801_dp+y*(-0.01031555_dp+y*(0.02282967_dp+y*(-0.02895312_dp+ &
            y*(0.01787654_dp-y*0.00420059_dp))))))))
      end if
      if (x < 0.0_dp) value = -value
   end function bessel_i1_approx

   pure function bessel_k0_approx(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: value, y
      if (x <= 2.0_dp) then
         y = 0.25_dp*x*x
         value = -log(0.5_dp*x)*bessel_i0_approx(x)+(-0.57721566_dp+y*(0.42278420_dp+ &
            y*(0.23069756_dp+y*(0.03488590_dp+y*(0.00262698_dp+ &
            y*(0.00010750_dp+y*0.00000740_dp))))))
      else
         y = 2.0_dp/x
         value = exp(-x)/sqrt(x)*(1.25331414_dp+y*(-0.07832358_dp+y*(0.02189568_dp+ &
            y*(-0.01062446_dp+y*(0.00587872_dp+y*(-0.00251540_dp+y*0.00053208_dp))))))
      end if
   end function bessel_k0_approx

   pure function bessel_k1_approx(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: value, y
      if (x <= 2.0_dp) then
         y = 0.25_dp*x*x
         value = log(0.5_dp*x)*bessel_i1_approx(x)+(1.0_dp/x)*(1.0_dp+y*(0.15443144_dp+ &
            y*(-0.67278579_dp+y*(-0.18156897_dp+y*(-0.01919402_dp+ &
            y*(-0.00110404_dp+y*(-0.00004686_dp)))))))
      else
         y = 2.0_dp/x
         value = exp(-x)/sqrt(x)*(1.25331414_dp+y*(0.23498619_dp+y*(-0.03655620_dp+ &
            y*(0.01504268_dp+y*(-0.00780353_dp+y*(0.00325614_dp-y*0.00068245_dp))))))
      end if
   end function bessel_k1_approx

   function bessel_k_integral(order, x) result(value)
      real(dp), intent(in) :: order, x
      real(dp) :: value, h, t, sumv
      integer, parameter :: n = 2000
      integer :: i
      h = 20.0_dp/real(n,dp)
      sumv = integrand(0.0_dp)+integrand(20.0_dp)
      do i = 1, n-1
         t = real(i,dp)*h
         sumv = sumv+merge(4.0_dp,2.0_dp,mod(i,2) == 1)*integrand(t)
      end do
      value = h*sumv/3.0_dp
   contains
      function integrand(tin) result(y)
         real(dp), intent(in) :: tin
         real(dp) :: y, exponent
         exponent = -x*cosh(tin)+order*tin
         if (exponent < -700.0_dp) then
            y = 0.0_dp
         else
            y = exp(-x*cosh(tin))*cosh(order*tin)
         end if
      end function integrand
   end function bessel_k_integral

   pure function regularized_beta(a, b, x) result(value)
      real(dp), intent(in) :: a, b, x
      real(dp) :: value, bt
      if (x <= 0.0_dp) then
         value = 0.0_dp
      else if (x >= 1.0_dp) then
         value = 1.0_dp
      else
         bt = exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b)+a*log(x)+b*log(1.0_dp-x))
         if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
            value = bt*beta_continued_fraction(a,b,x)/a
         else
            value = 1.0_dp-bt*beta_continued_fraction(b,a,1.0_dp-x)/b
         end if
      end if
      value = min(1.0_dp,max(0.0_dp,value))
   end function regularized_beta

   pure function beta_continued_fraction(a, b, x) result(value)
      real(dp), intent(in) :: a, b, x
      real(dp) :: value, qab, qap, qam, c, d, h, aa, delta
      integer :: m, m2
      integer, parameter :: maxit = 300
      real(dp), parameter :: eps = 3.0e-14_dp, fpmin = 1.0e-300_dp

      qab = a+b
      qap = a+1.0_dp
      qam = a-1.0_dp
      c = 1.0_dp
      d = 1.0_dp-qab*x/qap
      if (abs(d) < fpmin) d = fpmin
      d = 1.0_dp/d
      h = d
      do m = 1, maxit
         m2 = 2*m
         aa = real(m,dp)*(b-real(m,dp))*x/((qam+real(m2,dp))*(a+real(m2,dp)))
         d = 1.0_dp+aa*d
         if (abs(d) < fpmin) d = fpmin
         c = 1.0_dp+aa/c
         if (abs(c) < fpmin) c = fpmin
         d = 1.0_dp/d
         h = h*d*c
         aa = -(a+real(m,dp))*(qab+real(m,dp))*x/ &
            ((a+real(m2,dp))*(qap+real(m2,dp)))
         d = 1.0_dp+aa*d
         if (abs(d) < fpmin) d = fpmin
         c = 1.0_dp+aa/c
         if (abs(c) < fpmin) c = fpmin
         d = 1.0_dp/d
         delta = d*c
         h = h*delta
         if (abs(delta-1.0_dp) < eps) exit
      end do
      value = h
   end function beta_continued_fraction

   subroutine sort_real(x)
      real(dp), intent(inout) :: x(:)
      integer :: i, j
      real(dp) :: key
      do i = 2, size(x)
         key = x(i)
         j = i-1
         do while (j >= 1)
            if (x(j) <= key) exit
            x(j+1) = x(j)
            j = j-1
         end do
         x(j+1) = key
      end do
   end subroutine sort_real

   pure integer function lower_bound(x, value) result(index)
      real(dp), intent(in) :: x(:), value
      integer :: lo, hi, mid
      lo = 1
      hi = size(x)+1
      do while (lo < hi)
         mid = (lo+hi)/2
         if (mid <= size(x) .and. x(mid) < value) then
            lo = mid+1
         else
            hi = mid
         end if
      end do
      index = lo
   end function lower_bound

   pure integer function upper_bound(x, value) result(index)
      real(dp), intent(in) :: x(:), value
      integer :: lo, hi, mid
      lo = 1
      hi = size(x)+1
      do while (lo < hi)
         mid = (lo+hi)/2
         if (mid <= size(x) .and. x(mid) <= value) then
            lo = mid+1
         else
            hi = mid
         end if
      end do
      index = lo
   end function upper_bound

end module rmgarch_distributions
