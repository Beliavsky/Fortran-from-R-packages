module rebayes_math
   use rebayes_kinds, only : dp, pi
   implicit none
   private
   public :: normal_pdf, normal_cdf, student_pdf, poisson_pmf, poisson_cdf
   public :: binomial_pmf, gamma_pdf, weibull_pdf, weibull_cdf, gompertz_pdf
   public :: huber_eps, huber_pdf, linspace, normalize_prob, safe_log
   public :: solve_spd, soft_threshold, log1pexp, rng_normal, sample_discrete
   public :: log_beta, beta_pdf, traprule_values, noncentral_t_pdf
contains
   elemental real(dp) function normal_pdf(x, mu, sigma) result(y)
      real(dp), intent(in) :: x, mu, sigma
      if (sigma <= 0.0_dp) then
         y = 0.0_dp
      else
         y = exp(-0.5_dp*((x-mu)/sigma)**2)/(sqrt(2.0_dp*pi)*sigma)
      end if
   end function normal_pdf

   elemental real(dp) function normal_cdf(x) result(y)
      real(dp), intent(in) :: x
      y = 0.5_dp*erfc(-x/sqrt(2.0_dp))
   end function normal_cdf

   elemental real(dp) function student_pdf(x, nu) result(y)
      real(dp), intent(in) :: x, nu
      if (nu <= 0.0_dp) then
         y = 0.0_dp
      else
         y = exp(log_gamma(0.5_dp*(nu+1.0_dp))-log_gamma(0.5_dp*nu)) &
            /sqrt(nu*pi) * (1.0_dp+x*x/nu)**(-0.5_dp*(nu+1.0_dp))
      end if
   end function student_pdf

   elemental real(dp) function poisson_pmf(k, lambda) result(y)
      integer, intent(in) :: k
      real(dp), intent(in) :: lambda
      if (k < 0 .or. lambda < 0.0_dp) then
         y = 0.0_dp
      else if (abs(lambda) <= tiny(1.0_dp)) then
         if (k == 0) then
            y = 1.0_dp
         else
            y = 0.0_dp
         end if
      else
         y = exp(real(k,dp)*log(lambda)-lambda-log_gamma(real(k+1,dp)))
      end if
   end function poisson_pmf

   real(dp) function poisson_cdf(k, lambda) result(y)
      integer, intent(in) :: k
      real(dp), intent(in) :: lambda
      integer :: j
      real(dp) :: term
      if (k < 0) then
         y = 0.0_dp
         return
      end if
      if (lambda <= 0.0_dp) then
         y = 1.0_dp
         return
      end if
      term = exp(-lambda)
      y = term
      do j = 1, k
         term = term*lambda/real(j,dp)
         y = y + term
      end do
      y = min(1.0_dp, max(0.0_dp,y))
   end function poisson_cdf

   elemental real(dp) function binomial_pmf(x, n, p) result(y)
      integer, intent(in) :: x, n
      real(dp), intent(in) :: p
      real(dp) :: lp
      if (x < 0 .or. x > n .or. p < 0.0_dp .or. p > 1.0_dp) then
         y = 0.0_dp
      else if (p <= tiny(1.0_dp)) then
         y = merge(1.0_dp, 0.0_dp, x == 0)
      else if (p >= 1.0_dp-epsilon(1.0_dp)) then
         y = merge(1.0_dp, 0.0_dp, x == n)
      else
         lp = log_gamma(real(n+1,dp))-log_gamma(real(x+1,dp))-log_gamma(real(n-x+1,dp)) &
            + real(x,dp)*log(p) + real(n-x,dp)*log(1.0_dp-p)
         y = exp(lp)
      end if
   end function binomial_pmf

   elemental real(dp) function gamma_pdf(x, shape, scale) result(y)
      real(dp), intent(in) :: x, shape, scale
      if (x <= 0.0_dp .or. shape <= 0.0_dp .or. scale <= 0.0_dp) then
         y = 0.0_dp
      else
         y = exp((shape-1.0_dp)*log(x)-x/scale-log_gamma(shape)-shape*log(scale))
      end if
   end function gamma_pdf

   elemental real(dp) function weibull_pdf(x, shape, scale) result(y)
      real(dp), intent(in) :: x, shape, scale
      if (x < 0.0_dp .or. shape <= 0.0_dp .or. scale <= 0.0_dp) then
         y = 0.0_dp
      else if (x <= tiny(1.0_dp) .and. shape < 1.0_dp) then
         y = huge(1.0_dp)
      else
         y = shape/scale*(x/scale)**(shape-1.0_dp)*exp(-(x/scale)**shape)
      end if
   end function weibull_pdf

   elemental real(dp) function weibull_cdf(x, shape, scale) result(y)
      real(dp), intent(in) :: x, shape, scale
      if (x <= 0.0_dp) then
         y = 0.0_dp
      else if (shape <= 0.0_dp .or. scale <= 0.0_dp) then
         y = 0.0_dp
      else
         y = 1.0_dp-exp(-(x/scale)**shape)
      end if
   end function weibull_cdf

   elemental real(dp) function gompertz_pdf(x, alpha, theta) result(y)
      real(dp), intent(in) :: x, alpha, theta
      real(dp) :: u
      if (x <= 0.0_dp .or. theta <= 0.0_dp .or. abs(alpha) <= tiny(1.0_dp)) then
         y = 0.0_dp
      else
         u = alpha*x
         y = theta*exp(u)*exp((theta/alpha)*(1.0_dp-exp(u)))
      end if
   end function gompertz_pdf

   elemental real(dp) function huber_eps(k) result(eps)
      real(dp), intent(in) :: k
      real(dp) :: q
      q = 2.0_dp*normal_pdf(k,0.0_dp,1.0_dp)/k - 2.0_dp*normal_cdf(-k)
      eps = q/(1.0_dp+q)
   end function huber_eps

   elemental real(dp) function huber_pdf(x, sigma, k, heps) result(y)
      real(dp), intent(in) :: x, sigma, k, heps
      real(dp) :: z, base
      if (sigma <= 0.0_dp) then
         y = 0.0_dp
         return
      end if
      z = x/sigma
      base = (1.0_dp-heps)/sigma
      if (z > k) then
         y = base*normal_pdf(k,0.0_dp,1.0_dp)*exp(-k*(z-k))
      else if (z < -k) then
         y = base*normal_pdf(k,0.0_dp,1.0_dp)*exp(k*(z+k))
      else
         y = base*normal_pdf(z,0.0_dp,1.0_dp)
      end if
   end function huber_pdf

   function linspace(a, b, n) result(x)
      real(dp), intent(in) :: a, b
      integer, intent(in) :: n
      real(dp), allocatable :: x(:)
      integer :: i
      allocate(x(n))
      if (n <= 1) then
         if (n == 1) x(1) = a
         return
      end if
      do i = 1, n
         x(i) = a + (b-a)*real(i-1,dp)/real(n-1,dp)
      end do
   end function linspace

   subroutine normalize_prob(x)
      real(dp), intent(inout) :: x(:)
      real(dp) :: s
      x = max(x,0.0_dp)
      s = sum(x)
      if (s > 0.0_dp) then
         x = x/s
      else
         x = 1.0_dp/real(size(x),dp)
      end if
   end subroutine normalize_prob

   elemental real(dp) function safe_log(x) result(y)
      real(dp), intent(in) :: x
      y = log(max(x,tiny(1.0_dp)))
   end function safe_log

   elemental real(dp) function soft_threshold(x, t) result(y)
      real(dp), intent(in) :: x, t
      y = sign(max(abs(x)-t,0.0_dp),x)
   end function soft_threshold

   elemental real(dp) function log1pexp(x) result(y)
      real(dp), intent(in) :: x
      if (x > 0.0_dp) then
         y = x + log(1.0_dp+exp(-x))
      else
         y = log(1.0_dp+exp(x))
      end if
   end function log1pexp

   subroutine solve_spd(a, b, ok)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(inout) :: b(:)
      logical, intent(out) :: ok
      real(dp), allocatable :: l(:,:)
      real(dp) :: s
      integer :: n, i, j, k
      n = size(b)
      allocate(l(n,n)); l = 0.0_dp
      ok = .true.
      do i = 1, n
         do j = 1, i
            s = a(i,j)
            do k = 1, j-1
               s = s - l(i,k)*l(j,k)
            end do
            if (i == j) then
               if (s <= epsilon(1.0_dp)*max(1.0_dp,maxval(abs(a)))) then
                  ok = .false.; return
               end if
               l(i,j) = sqrt(s)
            else
               l(i,j) = s/l(j,j)
            end if
         end do
      end do
      do i = 1, n
         s = b(i)
         do k = 1, i-1
            s = s-l(i,k)*b(k)
         end do
         b(i) = s/l(i,i)
      end do
      do i = n, 1, -1
         s = b(i)
         do k = i+1, n
            s = s-l(k,i)*b(k)
         end do
         b(i) = s/l(i,i)
      end do
   end subroutine solve_spd

   subroutine rng_normal(z)
      real(dp), intent(out) :: z(:)
      integer :: i
      real(dp) :: u1, u2
      i = 1
      do while (i <= size(z))
         call random_number(u1); call random_number(u2)
         u1 = max(u1,tiny(1.0_dp))
         z(i) = sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
         if (i+1 <= size(z)) z(i+1) = sqrt(-2.0_dp*log(u1))*sin(2.0_dp*pi*u2)
         i = i + 2
      end do
   end subroutine rng_normal

   integer function sample_discrete(prob) result(idx)
      real(dp), intent(in) :: prob(:)
      real(dp) :: u, s
      integer :: i
      call random_number(u)
      s = 0.0_dp
      idx = size(prob)
      do i = 1, size(prob)
         s = s + prob(i)
         if (u <= s) then
            idx = i
            return
         end if
      end do
   end function sample_discrete

   elemental real(dp) function log_beta(a,b) result(y)
      real(dp), intent(in) :: a,b
      y = log_gamma(a)+log_gamma(b)-log_gamma(a+b)
   end function log_beta

   elemental real(dp) function beta_pdf(x,a,b) result(y)
      real(dp), intent(in) :: x,a,b
      if (x <= 0.0_dp .or. x >= 1.0_dp .or. a <= 0.0_dp .or. b <= 0.0_dp) then
         y = 0.0_dp
      else
         y = exp((a-1.0_dp)*log(x)+(b-1.0_dp)*log(1.0_dp-x)-log_beta(a,b))
      end if
   end function beta_pdf

   real(dp) function noncentral_t_pdf(x, nu, delta) result(y)
      real(dp), intent(in) :: x, nu, delta
      integer, parameter :: nq = 96
      real(dp), save :: nodes(nq), weights(nq)
      logical, save :: initialized = .false.
      real(dp) :: t, u, jac, a, z, lg
      integer :: i
      if (nu <= 0.0_dp) then
         y = 0.0_dp; return
      end if
      if (.not. initialized) then
         call gauss_legendre_01(nq,nodes,weights)
         initialized = .true.
      end if
      a = 0.5_dp*nu
      lg = log_gamma(a)
      y = 0.0_dp
      do i = 1, nq
         t = nodes(i)
         u = t/(1.0_dp-t)
         jac = 1.0_dp/(1.0_dp-t)**2
         z = x*sqrt(2.0_dp*u/nu)-delta
         y = y + weights(i)*sqrt(2.0_dp*u/nu)*normal_pdf(z,0.0_dp,1.0_dp) &
            * exp((a-1.0_dp)*log(max(u,tiny(1.0_dp)))-u-lg)*jac
      end do
   end function noncentral_t_pdf

   subroutine gauss_legendre_01(n,x,w)
      integer, intent(in) :: n
      real(dp), intent(out) :: x(n), w(n)
      integer :: i,j,m
      real(dp) :: z,z1,p1,p2,p3,pp,xx,ww
      m=(n+1)/2
      do i=1,m
         z=cos(pi*(real(i,dp)-0.25_dp)/(real(n,dp)+0.5_dp))
         do
            p1=1.0_dp; p2=0.0_dp
            do j=1,n
               p3=p2; p2=p1
               p1=((2.0_dp*real(j,dp)-1.0_dp)*z*p2-(real(j,dp)-1.0_dp)*p3)/real(j,dp)
            end do
            pp=real(n,dp)*(z*p1-p2)/(z*z-1.0_dp)
            z1=z; z=z1-p1/pp
            if(abs(z-z1) < 4.0_dp*epsilon(1.0_dp)) exit
         end do
         xx=0.5_dp*(1.0_dp-z); ww=1.0_dp/((1.0_dp-z*z)*pp*pp)
         x(i)=xx; x(n+1-i)=1.0_dp-xx
         w(i)=ww; w(n+1-i)=ww
      end do
   end subroutine gauss_legendre_01

   real(dp) function traprule_values(x,y) result(v)
      real(dp), intent(in) :: x(:), y(:)
      integer :: i
      v = 0.0_dp
      do i = 1, size(x)-1
         v = v + 0.5_dp*(x(i+1)-x(i))*(y(i+1)+y(i))
      end do
   end function traprule_values
end module rebayes_math
