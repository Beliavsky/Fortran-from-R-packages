module powerlaw
   use pracma_special, only : zeta_pracma => zeta
   implicit none
   private

   integer, parameter, public :: dp = kind(1.0d0)
   real(dp), parameter :: pi = acos(-1.0_dp)

   integer, parameter, public :: F_DISPL = 1, F_CONPL = 2
   integer, parameter, public :: F_DISEXP = 3, F_CONEXP = 4
   integer, parameter, public :: F_DISLNORM = 5, F_CONLNORM = 6
   integer, parameter, public :: F_DISPOIS = 7, F_CONWEIBULL = 8

   type, public :: powerlaw_dist
      integer :: family = 0
      real(dp), allocatable :: data(:)
      real(dp), allocatable :: pars(:)
      real(dp) :: xmin = 0.0_dp
   contains
      procedure :: set_data
      procedure :: set_xmin
      procedure :: set_pars
      procedure :: pdf => dist_pdf_scalar
      procedure :: cdf => dist_cdf_scalar
      procedure :: loglik => dist_loglik
      procedure :: random => dist_random
      procedure :: get_n => model_get_n
      procedure :: get_ntail => model_get_ntail
      procedure :: no_pars
   end type

   type, public :: estimate_pars_result
      real(dp), allocatable :: pars(:)
      real(dp) :: ll = -huge(1.0_dp)
      integer :: status = 0
      integer :: iterations = 0
   end type

   type, public :: estimate_xmin_result
      real(dp) :: gof = huge(1.0_dp)
      real(dp) :: xmin = 0.0_dp
      real(dp), allocatable :: pars(:)
      integer :: ntail = 0
      integer :: status = 0
   end type

   type, public :: compare_result
      real(dp) :: test_statistic = 0.0_dp
      real(dp) :: p_one_sided = 1.0_dp
      real(dp) :: p_two_sided = 1.0_dp
      real(dp), allocatable :: x(:)
      real(dp), allocatable :: ratio(:)
   end type

   type, public :: bootstrap_result
      real(dp) :: gof = huge(1.0_dp)
      real(dp) :: p = -1.0_dp
      real(dp), allocatable :: xmin(:)
      real(dp), allocatable :: pars(:,:)
      real(dp), allocatable :: sim_gof(:)
      integer, allocatable :: ntail(:)
      integer :: successful = 0
   end type

   interface displ
      module procedure make_displ
   end interface
   interface conpl
      module procedure make_conpl
   end interface
   interface disexp
      module procedure make_disexp
   end interface
   interface conexp
      module procedure make_conexp
   end interface
   interface dislnorm
      module procedure make_dislnorm
   end interface
   interface conlnorm
      module procedure make_conlnorm
   end interface
   interface dispois
      module procedure make_dispois
   end interface
   interface conweibull
      module procedure make_conweibull
   end interface

   public :: displ, conpl, disexp, conexp, dislnorm, conlnorm
   public :: dispois, conweibull
   public :: dpldis, ppldis, rpldis, dplcon, pplcon, rplcon
   public :: dist_pdf, dist_cdf, dist_all_cdf, dist_ll, dist_rand
   public :: dist_data_cdf, dist_data_all_cdf
   public :: estimate_pars, estimate_xmin
   public :: get_distance_statistic, get_KS_statistic
   public :: compare_distributions
   public :: bootstrap, bootstrap_p
   public :: get_bootstrap_sims, get_bootstrap_p_sims
   public :: get_n, get_ntail, set_rng_seed

contains

   function make_dist(family, data) result(m)
      integer, intent(in) :: family
      real(dp), intent(in), optional :: data(:)
      type(powerlaw_dist) :: m
      m%family = family
      if (present(data)) call m%set_data(data)
   end function

   function make_displ(data) result(m)
      real(dp), intent(in), optional :: data(:)
      type(powerlaw_dist) :: m
      m = make_dist(F_DISPL, data)
   end function
   function make_conpl(data) result(m)
      real(dp), intent(in), optional :: data(:)
      type(powerlaw_dist) :: m
      m = make_dist(F_CONPL, data)
   end function
   function make_disexp(data) result(m)
      real(dp), intent(in), optional :: data(:)
      type(powerlaw_dist) :: m
      m = make_dist(F_DISEXP, data)
   end function
   function make_conexp(data) result(m)
      real(dp), intent(in), optional :: data(:)
      type(powerlaw_dist) :: m
      m = make_dist(F_CONEXP, data)
   end function
   function make_dislnorm(data) result(m)
      real(dp), intent(in), optional :: data(:)
      type(powerlaw_dist) :: m
      m = make_dist(F_DISLNORM, data)
   end function
   function make_conlnorm(data) result(m)
      real(dp), intent(in), optional :: data(:)
      type(powerlaw_dist) :: m
      m = make_dist(F_CONLNORM, data)
   end function
   function make_dispois(data) result(m)
      real(dp), intent(in), optional :: data(:)
      type(powerlaw_dist) :: m
      m = make_dist(F_DISPOIS, data)
   end function
   function make_conweibull(data) result(m)
      real(dp), intent(in), optional :: data(:)
      type(powerlaw_dist) :: m
      m = make_dist(F_CONWEIBULL, data)
   end function

   logical function is_discrete(family)
      integer, intent(in) :: family
      is_discrete = family == F_DISPL .or. family == F_DISEXP .or. &
                    family == F_DISLNORM .or. family == F_DISPOIS
   end function

   subroutine set_data(self, x)
      class(powerlaw_dist), intent(inout) :: self
      real(dp), intent(in) :: x(:)
      integer :: i
      if (size(x) == 0) then
         if (allocated(self%data)) deallocate(self%data)
         return
      end if
      if (any(x <= 0.0_dp)) error stop "powerlaw: data must be positive"
      if (is_discrete(self%family)) then
         do i = 1, size(x)
            if (abs(x(i) - anint(x(i))) > 10.0_dp*epsilon(1.0_dp)) &
               error stop "powerlaw: discrete data must be integral"
         end do
      end if
      self%data = x
      call sort_real(self%data)
      self%xmin = minval(self%data)
   end subroutine

   subroutine set_xmin(self, xmin)
      class(powerlaw_dist), intent(inout) :: self
      real(dp), intent(in) :: xmin
      if (is_discrete(self%family)) then
         self%xmin = real(floor(xmin), dp)
      else
         self%xmin = xmin
      end if
   end subroutine

   subroutine set_pars(self, pars)
      class(powerlaw_dist), intent(inout) :: self
      real(dp), intent(in) :: pars(:)
      if (size(pars) /= self%no_pars()) error stop "powerlaw: wrong parameter count"
      self%pars = pars
   end subroutine

   integer function no_pars(self)
      class(powerlaw_dist), intent(in) :: self
      if (self%family == F_DISLNORM .or. self%family == F_CONLNORM .or. &
          self%family == F_CONWEIBULL) then
         no_pars = 2
      else
         no_pars = 1
      end if
   end function

   integer function model_get_n(self)
      class(powerlaw_dist), intent(in) :: self
      if (allocated(self%data)) then
         model_get_n = size(self%data)
      else
         model_get_n = 0
      end if
   end function

   integer function model_get_ntail(self)
      class(powerlaw_dist), intent(in) :: self
      if (allocated(self%data)) then
         model_get_ntail = count(self%data >= self%xmin)
      else
         model_get_ntail = 0
      end if
   end function

   subroutine sort_real(x)
      real(dp), intent(inout) :: x(:)
      integer :: i, j
      real(dp) :: v
      do i = 2, size(x)
         v = x(i)
         j = i - 1
         do while (j >= 1)
            if (x(j) <= v) exit
            x(j+1) = x(j)
            j = j - 1
         end do
         x(j+1) = v
      end do
   end subroutine

   subroutine unique_sorted(x, u)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable, intent(out) :: u(:)
      real(dp), allocatable :: s(:), tmp(:)
      integer :: i, n
      if (size(x) == 0) then
         allocate(u(0))
         return
      end if
      s = x
      call sort_real(s)
      allocate(tmp(size(s)))
      n = 1
      tmp(1) = s(1)
      do i = 2, size(s)
         if (s(i) /= tmp(n)) then
            n = n + 1
            tmp(n) = s(i)
         end if
      end do
      allocate(u(n))
      u = tmp(:n)
   end subroutine

   elemental real(dp) function normal_cdf(x)
      real(dp), intent(in) :: x
      normal_cdf = 0.5_dp*erfc(-x/sqrt(2.0_dp))
   end function

   elemental real(dp) function normal_sf(x)
      real(dp), intent(in) :: x
      normal_sf = 0.5_dp*erfc(x/sqrt(2.0_dp))
   end function

   elemental real(dp) function normal_pdf(x)
      real(dp), intent(in) :: x
      normal_pdf = exp(-0.5_dp*x*x)/sqrt(2.0_dp*pi)
   end function

   elemental real(dp) function lognormal_cdf(x, mu, sigma)
      real(dp), intent(in) :: x, mu, sigma
      if (x <= 0.0_dp) then
         lognormal_cdf = 0.0_dp
      else
         lognormal_cdf = normal_cdf((log(x)-mu)/sigma)
      end if
   end function

   elemental real(dp) function lognormal_sf(x, mu, sigma)
      real(dp), intent(in) :: x, mu, sigma
      if (x <= 0.0_dp) then
         lognormal_sf = 1.0_dp
      else
         lognormal_sf = normal_sf((log(x)-mu)/sigma)
      end if
   end function

   elemental real(dp) function lognormal_pdf(x, mu, sigma)
      real(dp), intent(in) :: x, mu, sigma
      real(dp) :: z
      if (x <= 0.0_dp .or. sigma <= 0.0_dp) then
         lognormal_pdf = 0.0_dp
      else
         z = (log(x)-mu)/sigma
         lognormal_pdf = normal_pdf(z)/(x*sigma)
      end if
   end function

   elemental real(dp) function weibull_sf(x, shape, scale)
      real(dp), intent(in) :: x, shape, scale
      if (x <= 0.0_dp) then
         weibull_sf = 1.0_dp
      else
         weibull_sf = exp(-(x/scale)**shape)
      end if
   end function

   elemental real(dp) function weibull_cdf(x, shape, scale)
      real(dp), intent(in) :: x, shape, scale
      weibull_cdf = 1.0_dp - weibull_sf(x, shape, scale)
   end function

   elemental real(dp) function weibull_pdf(x, shape, scale)
      real(dp), intent(in) :: x, shape, scale
      if (x <= 0.0_dp .or. shape <= 0.0_dp .or. scale <= 0.0_dp) then
         weibull_pdf = 0.0_dp
      else
         weibull_pdf = shape/scale*(x/scale)**(shape-1.0_dp) * &
                       exp(-(x/scale)**shape)
      end if
   end function

   real(dp) function poisson_pmf(k, lambda)
      integer, intent(in) :: k
      real(dp), intent(in) :: lambda
      if (k < 0 .or. lambda < 0.0_dp) then
         poisson_pmf = 0.0_dp
      else if (lambda == 0.0_dp) then
         if (k == 0) then
            poisson_pmf = 1.0_dp
         else
            poisson_pmf = 0.0_dp
         end if
      else
         poisson_pmf = exp(-lambda + real(k,dp)*log(lambda) - &
                           log_gamma(real(k+1,dp)))
      end if
   end function

   real(dp) function poisson_cdf(k, lambda)
      integer, intent(in) :: k
      real(dp), intent(in) :: lambda
      integer :: j
      real(dp) :: term, s
      if (k < 0) then
         poisson_cdf = 0.0_dp
         return
      end if
      if (lambda == 0.0_dp) then
         poisson_cdf = 1.0_dp
         return
      end if
      if (lambda > 700.0_dp) then
         poisson_cdf = normal_cdf((real(k,dp)+0.5_dp-lambda)/sqrt(lambda))
         return
      end if
      term = exp(-lambda)
      s = term
      do j = 1, k
         term = term*lambda/real(j,dp)
         s = s + term
      end do
      poisson_cdf = min(1.0_dp, max(0.0_dp, s))
   end function

   real(dp) function poisson_sf(k, lambda)
      integer, intent(in) :: k
      real(dp), intent(in) :: lambda
      poisson_sf = max(0.0_dp, 1.0_dp-poisson_cdf(k,lambda))
   end function

   real(dp) function rand_uniform()
      call random_number(rand_uniform)
      if (rand_uniform <= 0.0_dp) rand_uniform = spacing(1.0_dp)
      if (rand_uniform >= 1.0_dp) rand_uniform = 1.0_dp-spacing(1.0_dp)
   end function

   real(dp) function rand_normal()
      real(dp) :: u1, u2
      u1 = rand_uniform()
      u2 = rand_uniform()
      rand_normal = sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
   end function

   real(dp) function rand_lognormal(mu, sigma)
      real(dp), intent(in) :: mu, sigma
      rand_lognormal = exp(mu + sigma*rand_normal())
   end function

   real(dp) function rand_weibull(shape, scale)
      real(dp), intent(in) :: shape, scale
      rand_weibull = scale*(-log(rand_uniform()))**(1.0_dp/shape)
   end function

   integer function rand_poisson(lambda)
      real(dp), intent(in) :: lambda
      real(dp) :: u, p, cum
      integer :: k
      if (lambda <= 0.0_dp) then
         rand_poisson = 0
         return
      end if
      if (lambda > 100.0_dp) then
         do
            k = nint(lambda + sqrt(lambda)*rand_normal())
            if (k >= 0) exit
         end do
         rand_poisson = k
         return
      end if
      u = rand_uniform()
      p = exp(-lambda)
      cum = p
      k = 0
      do while (u > cum)
         k = k + 1
         p = p*lambda/real(k,dp)
         cum = cum + p
         if (p <= tiny(1.0_dp)) exit
      end do
      rand_poisson = k
   end function

   real(dp) function sample_mean(x)
      real(dp), intent(in) :: x(:)
      sample_mean = sum(x)/real(size(x),dp)
   end function

   real(dp) function sample_sd(x)
      real(dp), intent(in) :: x(:)
      real(dp) :: m
      if (size(x) < 2) then
         sample_sd = 0.0_dp
      else
         m = sample_mean(x)
         sample_sd = sqrt(sum((x-m)**2)/real(size(x)-1,dp))
      end if
   end function

   real(dp) function hurwitz_zeta_fast(s, q)
      real(dp), intent(in) :: s, q
      real(dp), parameter :: c(8) = [ &
         1.0_dp/12.0_dp, -1.0_dp/720.0_dp, 1.0_dp/30240.0_dp, &
         -1.0_dp/1209600.0_dp, 1.0_dp/47900160.0_dp, &
         -691.0_dp/1307674368000.0_dp, 1.0_dp/74724249600.0_dp, &
         -3617.0_dp/10670622842880000.0_dp ]
      real(dp) :: a, rising
      integer :: k, j, n
      if (s <= 1.0_dp .or. q <= 0.0_dp) then
         hurwitz_zeta_fast = huge(1.0_dp)
         return
      end if
      ! Euler-Maclaurin evaluation of Hurwitz zeta.  Twenty-four direct
      ! terms give double-precision accuracy throughout the fitted range.
      n = 24
      hurwitz_zeta_fast = 0.0_dp
      do k = 0, n-1
         hurwitz_zeta_fast = hurwitz_zeta_fast + (q+real(k,dp))**(-s)
      end do
      a = q + real(n,dp)
      hurwitz_zeta_fast = hurwitz_zeta_fast + a**(1.0_dp-s)/(s-1.0_dp) + &
                          0.5_dp*a**(-s)
      rising = s
      do j = 1, size(c)
         if (j > 1) then
            rising = rising*(s+real(2*j-3,dp))*(s+real(2*j-2,dp))
         end if
         hurwitz_zeta_fast = hurwitz_zeta_fast + &
                             c(j)*rising*a**(-s-real(2*j-1,dp))
      end do
   end function

   real(dp) function discrete_pl_constant(xmin, alpha)
      real(dp), intent(in) :: xmin, alpha
      real(dp) :: q
      q = real(int(floor(xmin)),dp)
      if (alpha >= 12.0_dp .and. q == 1.0_dp) then
         ! For rapidly convergent cases, reuse the supplied pracma zeta.
         discrete_pl_constant = zeta_pracma(alpha)
      else
         discrete_pl_constant = hurwitz_zeta_fast(alpha,q)
      end if
   end function

   real(dp) function dpldis(x, xmin, alpha, log_p)
      real(dp), intent(in) :: x, xmin, alpha
      logical, intent(in), optional :: log_p
      logical :: lp
      real(dp) :: lv, c
      lp = .false.
      if (present(log_p)) lp = log_p
      if (nint(x) < int(floor(xmin)) .or. x <= 0.0_dp .or. alpha <= 1.0_dp) then
         if (lp) then
            dpldis = -huge(1.0_dp)
         else
            dpldis = 0.0_dp
         end if
         return
      end if
      c = discrete_pl_constant(xmin, alpha)
      lv = -alpha*log(x) - log(c)
      if (lp) then
         dpldis = lv
      else
         dpldis = exp(lv)
      end if
   end function

   real(dp) function ppldis(q, xmin, alpha, lower_tail)
      real(dp), intent(in) :: q, xmin, alpha
      logical, intent(in), optional :: lower_tail
      logical :: lt
      real(dp) :: c, s
      integer :: k, lo, hi
      lt = .true.
      if (present(lower_tail)) lt = lower_tail
      lo = int(floor(xmin))
      hi = int(floor(q))
      if (hi < lo) then
         if (lt) then
            ppldis = 0.0_dp
         else
            ppldis = 1.0_dp
         end if
         return
      end if
      c = discrete_pl_constant(xmin, alpha)
      s = 0.0_dp
      do k = lo, hi
         s = s + real(k,dp)**(-alpha)
      end do
      s = min(1.0_dp, s/c)
      if (lt) then
         ppldis = s
      else
         ppldis = 1.0_dp-s
      end if
   end function

   subroutine rpldis(x, xmin, alpha, discrete_max)
      real(dp), intent(out) :: x(:)
      real(dp), intent(in) :: xmin, alpha
      integer, intent(in), optional :: discrete_max
      integer :: i, k, kmax, lo
      real(dp) :: u, c, s
      kmax = 10000
      if (present(discrete_max)) kmax = discrete_max
      lo = int(floor(xmin))
      c = discrete_pl_constant(xmin, alpha)
      do i = 1, size(x)
         u = rand_uniform()
         if (kmax >= lo) then
            s = 0.0_dp
            do k = lo, kmax
               s = s + real(k,dp)**(-alpha)/c
               if (u <= s) exit
            end do
            if (k <= kmax) then
               x(i) = real(k,dp)
               cycle
            end if
         end if
         x(i) = floor((real(lo,dp)-0.5_dp) * &
                (1.0_dp-u)**(-1.0_dp/(alpha-1.0_dp)) + 0.5_dp)
      end do
   end subroutine

   real(dp) function dplcon(x, xmin, alpha, log_p)
      real(dp), intent(in) :: x, xmin, alpha
      logical, intent(in), optional :: log_p
      logical :: lp
      real(dp) :: lv
      lp = .false.
      if (present(log_p)) lp = log_p
      if (x < xmin .or. xmin <= 0.0_dp .or. alpha <= 1.0_dp) then
         if (lp) then
            dplcon = -huge(1.0_dp)
         else
            dplcon = 0.0_dp
         end if
         return
      end if
      lv = log(alpha-1.0_dp) - log(xmin) - alpha*log(x/xmin)
      if (lp) then
         dplcon = lv
      else
         dplcon = exp(lv)
      end if
   end function

   real(dp) function pplcon(q, xmin, alpha, lower_tail)
      real(dp), intent(in) :: q, xmin, alpha
      logical, intent(in), optional :: lower_tail
      logical :: lt
      real(dp) :: p
      lt = .true.
      if (present(lower_tail)) lt = lower_tail
      if (q < xmin) then
         p = 0.0_dp
      else
         p = 1.0_dp - (q/xmin)**(1.0_dp-alpha)
      end if
      if (lt) then
         pplcon = p
      else
         pplcon = 1.0_dp-p
      end if
   end function

   subroutine rplcon(x, xmin, alpha)
      real(dp), intent(out) :: x(:)
      real(dp), intent(in) :: xmin, alpha
      integer :: i
      do i = 1, size(x)
         x(i) = xmin*(1.0_dp-rand_uniform())**(-1.0_dp/(alpha-1.0_dp))
      end do
   end subroutine

   real(dp) function dist_pdf_scalar(self, q, log_p)
      class(powerlaw_dist), intent(in) :: self
      real(dp), intent(in) :: q
      logical, intent(in), optional :: log_p
      logical :: lp
      real(dp) :: p, sf, mass
      integer :: k
      if (.not. allocated(self%pars)) error stop "powerlaw: parameters not set"
      lp = .false.
      if (present(log_p)) lp = log_p
      select case (self%family)
      case (F_DISPL)
         dist_pdf_scalar = dpldis(q,self%xmin,self%pars(1),lp)
         return
      case (F_CONPL)
         dist_pdf_scalar = dplcon(q,self%xmin,self%pars(1),lp)
         return
      case (F_DISEXP)
         if (q < self%xmin) then
            p = 0.0_dp
         else
            p = (1.0_dp-exp(-self%pars(1))) * &
                exp(-self%pars(1)*real(nint(q)-nint(self%xmin),dp))
         end if
      case (F_CONEXP)
         if (q < self%xmin) then
            p = 0.0_dp
         else
            p = self%pars(1)*exp(-self%pars(1)*(q-self%xmin))
         end if
      case (F_DISLNORM)
         if (q < self%xmin) then
            p = 0.0_dp
         else
            mass = lognormal_cdf(q+0.5_dp,self%pars(1),self%pars(2)) - &
                   lognormal_cdf(q-0.5_dp,self%pars(1),self%pars(2))
            sf = lognormal_sf(self%xmin-0.5_dp,self%pars(1),self%pars(2))
            p = mass/sf
         end if
      case (F_CONLNORM)
         if (q < self%xmin) then
            p = 0.0_dp
         else
            sf = lognormal_sf(self%xmin,self%pars(1),self%pars(2))
            p = lognormal_pdf(q,self%pars(1),self%pars(2))/sf
         end if
      case (F_DISPOIS)
         k = nint(q)
         if (q < self%xmin) then
            p = 0.0_dp
         else
            sf = poisson_sf(nint(self%xmin)-1,self%pars(1))
            p = poisson_pmf(k,self%pars(1))/sf
         end if
      case (F_CONWEIBULL)
         if (q < self%xmin) then
            p = 0.0_dp
         else
            sf = weibull_sf(self%xmin,self%pars(1),self%pars(2))
            p = weibull_pdf(q,self%pars(1),self%pars(2))/sf
         end if
      case default
         p = 0.0_dp
      end select
      if (lp) then
         if (p > 0.0_dp) then
            dist_pdf_scalar = log(p)
         else
            dist_pdf_scalar = -huge(1.0_dp)
         end if
      else
         dist_pdf_scalar = p
      end if
   end function

   real(dp) function dist_cdf_scalar(self, q, lower_tail)
      class(powerlaw_dist), intent(in) :: self
      real(dp), intent(in) :: q
      logical, intent(in), optional :: lower_tail
      logical :: lt
      real(dp) :: p, f0, s0
      integer :: k
      if (.not. allocated(self%pars)) error stop "powerlaw: parameters not set"
      lt = .true.
      if (present(lower_tail)) lt = lower_tail
      if (q < self%xmin) then
         if (lt) then
            dist_cdf_scalar = 0.0_dp
         else
            dist_cdf_scalar = 1.0_dp
         end if
         return
      end if
      select case (self%family)
      case (F_DISPL)
         p = ppldis(q,self%xmin,self%pars(1),.true.)
      case (F_CONPL)
         p = pplcon(q,self%xmin,self%pars(1),.true.)
      case (F_DISEXP)
         k = int(floor(q))
         p = 1.0_dp - exp(-self%pars(1)*real(k-nint(self%xmin)+1,dp))
      case (F_CONEXP)
         p = 1.0_dp - exp(-self%pars(1)*(q-self%xmin))
      case (F_DISLNORM)
         f0 = lognormal_cdf(self%xmin-0.5_dp,self%pars(1),self%pars(2))
         s0 = 1.0_dp-f0
         p = (lognormal_cdf(floor(q)+0.5_dp,self%pars(1),self%pars(2))-f0)/s0
      case (F_CONLNORM)
         f0 = lognormal_cdf(self%xmin,self%pars(1),self%pars(2))
         s0 = 1.0_dp-f0
         p = (lognormal_cdf(q,self%pars(1),self%pars(2))-f0)/s0
      case (F_DISPOIS)
         f0 = poisson_cdf(nint(self%xmin)-1,self%pars(1))
         s0 = 1.0_dp-f0
         p = (poisson_cdf(int(floor(q)),self%pars(1))-f0)/s0
      case (F_CONWEIBULL)
         f0 = weibull_cdf(self%xmin,self%pars(1),self%pars(2))
         s0 = 1.0_dp-f0
         p = (weibull_cdf(q,self%pars(1),self%pars(2))-f0)/s0
      case default
         p = 0.0_dp
      end select
      p = min(1.0_dp,max(0.0_dp,p))
      if (lt) then
         dist_cdf_scalar = p
      else
         dist_cdf_scalar = 1.0_dp-p
      end if
   end function

   real(dp) function dist_loglik(self)
      class(powerlaw_dist), intent(in) :: self
      real(dp), allocatable :: tail(:)
      real(dp) :: c
      integer :: i, n
      if (.not. allocated(self%data) .or. .not. allocated(self%pars)) then
         dist_loglik = -huge(1.0_dp)
         return
      end if
      tail = pack(self%data,self%data >= self%xmin)
      n = size(tail)
      if (n == 0) then
         dist_loglik = -huge(1.0_dp)
         return
      end if
      if (self%family == F_DISPL) then
         c = discrete_pl_constant(self%xmin,self%pars(1))
         dist_loglik = -self%pars(1)*sum(log(tail)) - real(n,dp)*log(c)
         return
      end if
      dist_loglik = 0.0_dp
      do i = 1, n
         dist_loglik = dist_loglik + self%pdf(tail(i),.true.)
      end do
   end function

   subroutine dist_random(self, x)
      class(powerlaw_dist), intent(in) :: self
      real(dp), intent(out) :: x(:)
      integer :: i, k
      real(dp) :: y
      if (.not. allocated(self%pars)) error stop "powerlaw: parameters not set"
      select case (self%family)
      case (F_DISPL)
         call rpldis(x,self%xmin,self%pars(1))
      case (F_CONPL)
         call rplcon(x,self%xmin,self%pars(1))
      case (F_DISEXP)
         do i = 1, size(x)
            y = self%xmin-0.5_dp - log(rand_uniform())/self%pars(1)
            x(i) = anint(y)
         end do
      case (F_CONEXP)
         do i = 1, size(x)
            x(i) = self%xmin - log(rand_uniform())/self%pars(1)
         end do
      case (F_DISLNORM)
         do i = 1, size(x)
            do
               y = rand_lognormal(self%pars(1),self%pars(2))
               if (y >= self%xmin-0.5_dp) exit
            end do
            x(i) = anint(y)
         end do
      case (F_CONLNORM)
         do i = 1, size(x)
            do
               y = rand_lognormal(self%pars(1),self%pars(2))
               if (y > self%xmin) exit
            end do
            x(i) = y
         end do
      case (F_DISPOIS)
         do i = 1, size(x)
            do
               k = rand_poisson(self%pars(1))
               if (real(k,dp) >= self%xmin) exit
            end do
            x(i) = real(k,dp)
         end do
      case (F_CONWEIBULL)
         do i = 1, size(x)
            do
               y = rand_weibull(self%pars(1),self%pars(2))
               if (y > self%xmin) exit
            end do
            x(i) = y
         end do
      end select
   end subroutine

   subroutine dist_pdf(m,q,p,log_p)
      type(powerlaw_dist), intent(in) :: m
      real(dp), intent(in) :: q(:)
      real(dp), intent(out) :: p(size(q))
      logical, intent(in), optional :: log_p
      integer :: i
      do i = 1, size(q)
         p(i) = m%pdf(q(i),log_p)
      end do
   end subroutine

   subroutine dist_cdf(m,q,p,lower_tail)
      type(powerlaw_dist), intent(in) :: m
      real(dp), intent(in) :: q(:)
      real(dp), intent(out) :: p(size(q))
      logical, intent(in), optional :: lower_tail
      integer :: i
      do i = 1, size(q)
         p(i) = m%cdf(q(i),lower_tail)
      end do
   end subroutine

   subroutine dist_all_cdf(m,p,lower_tail,xmax)
      type(powerlaw_dist), intent(in) :: m
      real(dp), allocatable, intent(out) :: p(:)
      logical, intent(in), optional :: lower_tail
      real(dp), intent(in), optional :: xmax
      real(dp) :: xm
      integer :: lo, hi, i
      xm = 1.0e5_dp
      if (present(xmax)) xm = xmax
      lo = int(ceiling(m%xmin))
      if (allocated(m%data)) then
         hi = int(floor(min(maxval(m%data),xm)))
      else
         hi = int(floor(xm))
      end if
      if (hi < lo) then
         allocate(p(0))
         return
      end if
      allocate(p(hi-lo+1))
      do i = lo, hi
         p(i-lo+1) = m%cdf(real(i,dp),lower_tail)
      end do
   end subroutine

   real(dp) function dist_ll(m)
      type(powerlaw_dist), intent(in) :: m
      dist_ll = m%loglik()
   end function

   subroutine dist_rand(m,x)
      type(powerlaw_dist), intent(in) :: m
      real(dp), intent(out) :: x(:)
      call m%random(x)
   end subroutine

   subroutine dist_data_cdf(m,p,lower_tail)
      type(powerlaw_dist), intent(in) :: m
      real(dp), allocatable, intent(out) :: p(:)
      logical, intent(in), optional :: lower_tail
      logical :: lt
      integer :: i, n, j
      real(dp), allocatable :: tail(:)
      lt = .true.
      if (present(lower_tail)) lt = lower_tail
      if (.not. allocated(m%data)) then
         allocate(p(0))
         return
      end if
      tail = pack(m%data,m%data >= m%xmin)
      n = size(tail)
      allocate(p(n))
      if (is_discrete(m%family)) then
         do i = 1, n
            j = count(tail <= tail(i))
            p(i) = real(j,dp)/real(n,dp)
         end do
      else
         do i = 1, n
            p(i) = real(i-1,dp)/real(n,dp)
         end do
      end if
      if (.not. lt) p = 1.0_dp-p
   end subroutine

   subroutine dist_data_all_cdf(m,p,lower_tail,xmax)
      type(powerlaw_dist), intent(in) :: m
      real(dp), allocatable, intent(out) :: p(:)
      logical, intent(in), optional :: lower_tail
      real(dp), intent(in), optional :: xmax
      logical :: lt
      real(dp) :: xm
      integer :: lo, hi, k, n
      lt = .true.
      if (present(lower_tail)) lt = lower_tail
      xm = 1.0e5_dp
      if (present(xmax)) xm = xmax
      if (.not. is_discrete(m%family)) then
         call dist_data_cdf(m,p,lt)
         return
      end if
      lo = nint(m%xmin)
      hi = int(floor(min(maxval(m%data),xm)))
      if (hi < lo) then
         allocate(p(0))
         return
      end if
      n = count(m%data >= m%xmin)
      allocate(p(hi-lo+1))
      do k = lo, hi
         p(k-lo+1) = real(count(m%data >= m%xmin .and. &
                         m%data <= real(k,dp)),dp)/real(n,dp)
      end do
      if (.not. lt) p = 1.0_dp-p
   end subroutine

   real(dp) function objective(m,u)
      type(powerlaw_dist), intent(in) :: m
      real(dp), intent(in) :: u(:)
      type(powerlaw_dist) :: w
      real(dp) :: p1(1), p2(2), ll
      w = m
      select case(m%family)
      case(F_DISPL)
         p1 = [1.0_dp+exp(u(1))]
         call w%set_pars(p1)
      case(F_DISPOIS)
         p1 = [exp(u(1))]
         call w%set_pars(p1)
      case(F_CONLNORM,F_DISLNORM)
         p2 = [u(1),exp(u(2))]
         call w%set_pars(p2)
      case(F_CONWEIBULL)
         p2 = [exp(u(1)),exp(u(2))]
         call w%set_pars(p2)
      case default
         objective = huge(1.0_dp)
         return
      end select
      ll = w%loglik()
      if (ll /= ll .or. ll < -0.9_dp*huge(1.0_dp)) then
         objective = huge(1.0_dp)/100.0_dp
      else
         objective = -ll
      end if
   end function

   subroutine golden_fit(m,a0,b0,x,iters)
      type(powerlaw_dist), intent(in) :: m
      real(dp), intent(in) :: a0,b0
      real(dp), intent(out) :: x
      integer, intent(out) :: iters
      real(dp) :: a,b,c,d,fc,fd,gr
      gr = (sqrt(5.0_dp)-1.0_dp)/2.0_dp
      a=a0; b=b0
      c=b-gr*(b-a); d=a+gr*(b-a)
      fc=objective(m,[c]); fd=objective(m,[d])
      do iters=1,250
         if (abs(b-a) < 1.0e-10_dp*(1.0_dp+abs(a)+abs(b))) exit
         if (fc < fd) then
            b=d; d=c; fd=fc
            c=b-gr*(b-a); fc=objective(m,[c])
         else
            a=c; c=d; fc=fd
            d=a+gr*(b-a); fd=objective(m,[d])
         end if
      end do
      x=0.5_dp*(a+b)
   end subroutine

   subroutine order3(f,b,m,w)
      real(dp), intent(in) :: f(3)
      integer, intent(out) :: b,m,w
      integer :: idx(3),i,j,t
      idx=[1,2,3]
      do i=1,2
         do j=i+1,3
            if (f(idx(j)) < f(idx(i))) then
               t=idx(i); idx(i)=idx(j); idx(j)=t
            end if
         end do
      end do
      b=idx(1); m=idx(2); w=idx(3)
   end subroutine

   subroutine nelder2_fit(m,start,u,iters)
      type(powerlaw_dist), intent(in) :: m
      real(dp), intent(in) :: start(2)
      real(dp), intent(out) :: u(2)
      integer, intent(out) :: iters
      real(dp) :: x(2,3), f(3), cent(2), xr(2), xe(2), xc(2)
      real(dp) :: fr,fe,fc
      integer :: i,best,mid,worst
      x(:,1)=start
      x(:,2)=start+[0.15_dp,0.0_dp]
      x(:,3)=start+[0.0_dp,0.15_dp]
      do i=1,3
         f(i)=objective(m,x(:,i))
      end do
      do iters=1,1000
         call order3(f,best,mid,worst)
         cent=0.5_dp*(x(:,best)+x(:,mid))
         xr=cent+(cent-x(:,worst))
         fr=objective(m,xr)
         if (fr < f(best)) then
            xe=cent+2.0_dp*(xr-cent)
            fe=objective(m,xe)
            if (fe < fr) then
               x(:,worst)=xe; f(worst)=fe
            else
               x(:,worst)=xr; f(worst)=fr
            end if
         else if (fr < f(mid)) then
            x(:,worst)=xr; f(worst)=fr
         else
            if (fr < f(worst)) then
               xc=cent+0.5_dp*(xr-cent)
            else
               xc=cent+0.5_dp*(x(:,worst)-cent)
            end if
            fc=objective(m,xc)
            if (fc < min(fr,f(worst))) then
               x(:,worst)=xc; f(worst)=fc
            else
               do i=1,3
                  if (i /= best) then
                     x(:,i)=x(:,best)+0.5_dp*(x(:,i)-x(:,best))
                     f(i)=objective(m,x(:,i))
                  end if
               end do
            end if
         end if
         if (maxval(abs(x(:,1)-x(:,2))) < 1.0e-8_dp .and. &
             maxval(abs(x(:,1)-x(:,3))) < 1.0e-8_dp) exit
      end do
      best=minloc(f,dim=1)
      u=x(:,best)
   end subroutine

   function estimate_pars(m,candidates) result(res)
      type(powerlaw_dist), intent(in) :: m
      real(dp), intent(in), optional :: candidates(:,:)
      type(estimate_pars_result) :: res
      type(powerlaw_dist) :: w
      real(dp), allocatable :: tail(:)
      real(dp) :: den, mean_ex, u, start(2), opt(2), ll, bestll
      integer :: n, i, best
      if (.not. allocated(m%data)) then
         res%status=1
         allocate(res%pars(m%no_pars()),source=0.0_dp)
         return
      end if
      w=m
      if (present(candidates)) then
         if (size(candidates,2) /= m%no_pars()) then
            res%status=5
            allocate(res%pars(m%no_pars()),source=0.0_dp)
            return
         end if
         best=0; bestll=-huge(1.0_dp)
         do i=1,size(candidates,1)
            call w%set_pars(candidates(i,:))
            ll=w%loglik()
            if (ll > bestll) then
               bestll=ll; best=i
            end if
         end do
         allocate(res%pars(m%no_pars()))
         if (best > 0) then
            res%pars=candidates(best,:)
            res%ll=bestll
         else
            res%status=6
         end if
         return
      end if
      tail=pack(m%data,m%data >= m%xmin)
      n=size(tail)
      if (n == 0) then
         res%status=2
         allocate(res%pars(m%no_pars()),source=0.0_dp)
         return
      end if
      select case(m%family)
      case(F_CONPL)
         den=sum(log(tail/m%xmin))
         allocate(res%pars(1))
         res%pars(1)=1.0_dp+real(n,dp)/max(den,tiny(1.0_dp))
      case(F_CONEXP)
         mean_ex=sum(tail-m%xmin)/real(n,dp)
         allocate(res%pars(1))
         res%pars(1)=1.0_dp/max(mean_ex,tiny(1.0_dp))
      case(F_DISEXP)
         mean_ex=sum(tail-m%xmin)/real(n,dp)
         allocate(res%pars(1))
         if (mean_ex <= 0.0_dp) then
            res%pars(1)=50.0_dp
         else
            res%pars(1)=log(1.0_dp+1.0_dp/mean_ex)
         end if
      case(F_DISPL)
         call golden_fit(w,log(1.0e-5_dp),log(99.0_dp),u,res%iterations)
         allocate(res%pars(1))
         res%pars(1)=1.0_dp+exp(u)
      case(F_DISPOIS)
         call golden_fit(w,log(1.0e-8_dp), &
              log(max(10.0_dp,10.0_dp*sample_mean(tail))),u,res%iterations)
         allocate(res%pars(1))
         res%pars(1)=exp(u)
      case(F_CONLNORM,F_DISLNORM)
         start=[sample_mean(log(tail)), &
                log(max(sample_sd(log(tail)),0.2_dp))]
         call nelder2_fit(w,start,opt,res%iterations)
         allocate(res%pars(2))
         res%pars=[opt(1),exp(opt(2))]
      case(F_CONWEIBULL)
         start=[0.0_dp,log(max(sample_mean(tail),tiny(1.0_dp)))]
         call nelder2_fit(w,start,opt,res%iterations)
         allocate(res%pars(2))
         res%pars=[exp(opt(1)),exp(opt(2))]
      case default
         res%status=3
         allocate(res%pars(1),source=0.0_dp)
         return
      end select
      call w%set_pars(res%pars)
      res%ll=w%loglik()
   end function

   real(dp) function get_distance_statistic(m,xmax,distance)
      type(powerlaw_dist), intent(in) :: m
      real(dp), intent(in), optional :: xmax
      character(len=*), intent(in), optional :: distance
      real(dp), allocatable :: tail(:)
      real(dp) :: xm,fit,emp,d,den
      integer :: i,n,k,lo,hi
      character(len=16) :: kind
      xm=1.0e5_dp
      if (present(xmax)) xm=xmax
      kind="ks"
      if (present(distance)) kind=trim(distance)
      get_distance_statistic=0.0_dp
      if (is_discrete(m%family)) then
         lo=nint(m%xmin)
         hi=int(floor(min(maxval(m%data),xm)))
         n=count(m%data >= m%xmin)
         do k=lo,hi
            fit=m%cdf(real(k,dp))
            emp=real(count(m%data >= m%xmin .and. &
                     m%data <= real(k,dp)),dp)/real(n,dp)
            d=abs(emp-fit)
            if (kind == "reweight") then
               den=sqrt(max(fit*(1.0_dp-fit),tiny(1.0_dp)))
               d=d/den
            end if
            get_distance_statistic=max(get_distance_statistic,d)
         end do
      else
         tail=pack(m%data,m%data >= m%xmin .and. m%data <= xm)
         n=size(tail)
         do i=1,n
            fit=m%cdf(tail(i))
            emp=real(i-1,dp)/real(n,dp)
            d=abs(emp-fit)
            if (kind == "reweight") then
               den=sqrt(max(fit*(1.0_dp-fit),tiny(1.0_dp)))
               d=d/den
            end if
            get_distance_statistic=max(get_distance_statistic,d)
         end do
      end if
   end function

   real(dp) function get_KS_statistic(m,xmax)
      type(powerlaw_dist), intent(in) :: m
      real(dp), intent(in), optional :: xmax
      get_KS_statistic=get_distance_statistic(m,xmax,"ks")
   end function

   function estimate_xmin(m,xmins,xmax,distance,candidates) result(res)
      type(powerlaw_dist), intent(in) :: m
      real(dp), intent(in), optional :: xmins(:),xmax,candidates(:,:)
      character(len=*), intent(in), optional :: distance
      type(estimate_xmin_result) :: res
      type(powerlaw_dist) :: w
      type(estimate_pars_result) :: fit
      real(dp), allocatable :: cand(:),udata(:)
      real(dp) :: xm,g
      integer :: i,np,nu
      character(len=16) :: kind
      if (.not. allocated(m%data)) then
         res%status=1
         allocate(res%pars(m%no_pars()),source=0.0_dp)
         return
      end if
      xm=1.0e5_dp
      if (present(xmax)) xm=xmax
      kind="ks"
      if (present(distance)) kind=trim(distance)
      np=m%no_pars()
      if (present(xmins)) then
         cand=pack(xmins,xmins <= xm)
      else
         call unique_sorted(m%data,cand)
         cand=pack(cand,cand <= xm)
      end if
      call unique_sorted(m%data,udata)
      nu=size(udata)
      if (nu <= np+1) then
         res%status=3
         allocate(res%pars(np),source=0.0_dp)
         return
      end if
      cand=pack(cand,cand <= udata(nu-np-1))
      if (size(cand) == 0) then
         res%status=4
         allocate(res%pars(np),source=0.0_dp)
         return
      end if
      res%gof=huge(1.0_dp)
      allocate(res%pars(np),source=0.0_dp)
      do i=1,size(cand)
         w=m
         call w%set_xmin(cand(i))
         fit=estimate_pars(w,candidates)
         if (fit%status /= 0) cycle
         call w%set_pars(fit%pars)
         g=get_distance_statistic(w,xm,kind)
         if (g < res%gof) then
            res%gof=g
            res%xmin=w%xmin
            res%pars=fit%pars
            res%ntail=w%get_ntail()
            res%status=0
         end if
      end do
      if (res%gof == huge(1.0_dp)) res%status=5
   end function

   function compare_distributions(d1,d2) result(r)
      type(powerlaw_dist), intent(in) :: d1,d2
      type(compare_result) :: r
      real(dp), allocatable :: q(:)
      real(dp) :: mu,sdv,p1
      integer :: i,n
      if (is_discrete(d1%family) .neqv. is_discrete(d2%family)) &
         error stop "powerlaw: cannot compare discrete and continuous models"
      if (abs(d1%xmin-d2%xmin) > 100.0_dp*epsilon(1.0_dp)) &
         error stop "powerlaw: compare requires equal xmin"
      q=pack(d1%data,d1%data >= d1%xmin)
      n=size(q)
      allocate(r%x(n),r%ratio(n))
      r%x=q
      do i=1,n
         r%ratio(i)=d1%pdf(q(i),.true.)-d2%pdf(q(i),.true.)
      end do
      mu=sum(r%ratio)/real(n,dp)
      sdv=sqrt(sum((r%ratio-mu)**2)/real(max(1,n-1),dp))
      if (sdv <= tiny(1.0_dp)) then
         r%test_statistic=0.0_dp
         r%p_one_sided=0.5_dp
         r%p_two_sided=1.0_dp
         return
      end if
      r%test_statistic=sqrt(real(n,dp))*mu/sdv
      p1=normal_cdf(r%test_statistic)
      r%p_one_sided=1.0_dp-p1
      r%p_two_sided=2.0_dp*min(p1,1.0_dp-p1)
   end function

   integer function get_n(m)
      type(powerlaw_dist), intent(in) :: m
      get_n=m%get_n()
   end function

   real(dp) function get_ntail(m,prop,lower)
      type(powerlaw_dist), intent(in) :: m
      logical, intent(in), optional :: prop,lower
      logical :: pr,lo
      integer :: n,nt
      pr=.false.; lo=.false.
      if (present(prop)) pr=prop
      if (present(lower)) lo=lower
      n=m%get_n(); nt=m%get_ntail()
      if (lo) nt=n-nt
      if (pr .and. n > 0) then
         get_ntail=real(nt,dp)/real(n,dp)
      else
         get_ntail=real(nt,dp)
      end if
   end function

   subroutine set_rng_seed(seed)
      integer, intent(in) :: seed
      integer :: n,i
      integer, allocatable :: s(:)
      call random_seed(size=n)
      allocate(s(n))
      do i=1,n
         s(i)=mod(seed+104729*i,2147483646)+1
      end do
      call random_seed(put=s)
   end subroutine

   subroutine get_bootstrap_sims(m,sims,seed)
      type(powerlaw_dist), intent(in) :: m
      real(dp), intent(out) :: sims(:,:)
      integer, intent(in), optional :: seed
      integer :: i,j,n,ix
      if (present(seed)) call set_rng_seed(seed)
      n=m%get_n()
      if (size(sims,1) /= n) error stop "powerlaw: bootstrap size mismatch"
      do j=1,size(sims,2)
         do i=1,n
            ix=1+int(rand_uniform()*real(n,dp))
            ix=min(n,max(1,ix))
            sims(i,j)=m%data(ix)
         end do
      end do
   end subroutine

   subroutine get_bootstrap_p_sims(m,sims,seed)
      type(powerlaw_dist), intent(in) :: m
      real(dp), intent(out) :: sims(:,:)
      integer, intent(in), optional :: seed
      real(dp), allocatable :: lower(:),one(:)
      real(dp) :: ptail
      integer :: i,j,n,ix
      if (.not. allocated(m%pars)) &
         error stop "powerlaw: parameters must be set for bootstrap_p simulations"
      if (present(seed)) call set_rng_seed(seed)
      n=m%get_n()
      if (size(sims,1) /= n) error stop "powerlaw: bootstrap size mismatch"
      lower=pack(m%data,m%data < m%xmin)
      ptail=real(m%get_ntail(),dp)/real(n,dp)
      allocate(one(1))
      do j=1,size(sims,2)
         do i=1,n
            if (rand_uniform() > ptail .and. size(lower) > 0) then
               ix=1+int(rand_uniform()*real(size(lower),dp))
               ix=min(size(lower),max(1,ix))
               sims(i,j)=lower(ix)
            else
               call m%random(one)
               sims(i,j)=one(1)
            end if
         end do
      end do
   end subroutine

   function bootstrap(m,no_of_sims,seed,xmax,distance,xmins,candidates) result(r)
      type(powerlaw_dist), intent(in) :: m
      integer, intent(in) :: no_of_sims
      integer, intent(in), optional :: seed
      real(dp), intent(in), optional :: xmax
      character(len=*), intent(in), optional :: distance
      real(dp), intent(in), optional :: xmins(:), candidates(:,:)
      type(bootstrap_result) :: r
      type(powerlaw_dist) :: w
      type(estimate_xmin_result) :: base,e
      real(dp), allocatable :: sample(:)
      integer :: i,j,n,np,ix
      if (present(seed)) call set_rng_seed(seed)
      base=estimate_xmin(m,xmins=xmins,xmax=xmax,distance=distance,candidates=candidates)
      r%gof=base%gof
      n=m%get_n(); np=m%no_pars()
      allocate(r%xmin(no_of_sims),r%pars(no_of_sims,np))
      allocate(r%sim_gof(no_of_sims),r%ntail(no_of_sims),sample(n))
      do i=1,no_of_sims
         do j=1,n
            ix=1+int(rand_uniform()*real(n,dp))
            ix=min(n,max(1,ix))
            sample(j)=m%data(ix)
         end do
         w=m
         call w%set_data(sample)
         e=estimate_xmin(w,xmins=xmins,xmax=xmax,distance=distance,candidates=candidates)
         r%xmin(i)=e%xmin
         r%pars(i,:)=e%pars
         r%sim_gof(i)=e%gof
         r%ntail(i)=e%ntail
         if (e%status == 0) r%successful=r%successful+1
      end do
   end function

   function bootstrap_p(m,no_of_sims,seed,xmax,distance,xmins,candidates) result(r)
      type(powerlaw_dist), intent(in) :: m
      integer, intent(in) :: no_of_sims
      integer, intent(in), optional :: seed
      real(dp), intent(in), optional :: xmax
      character(len=*), intent(in), optional :: distance
      real(dp), intent(in), optional :: xmins(:), candidates(:,:)
      type(bootstrap_result) :: r
      type(powerlaw_dist) :: base_model,w
      type(estimate_xmin_result) :: base,e
      real(dp), allocatable :: sample(:),lower(:),one(:)
      real(dp) :: ptail
      integer :: i,j,n,np,ix
      if (present(seed)) call set_rng_seed(seed)
      base=estimate_xmin(m,xmins=xmins,xmax=xmax,distance=distance,candidates=candidates)
      base_model=m
      if (.not. allocated(base_model%pars)) then
         call base_model%set_xmin(base%xmin)
         call base_model%set_pars(base%pars)
      end if
      r%gof=base%gof
      n=m%get_n(); np=m%no_pars()
      lower=pack(m%data,m%data < base_model%xmin)
      ptail=real(base_model%get_ntail(),dp)/real(n,dp)
      allocate(r%xmin(no_of_sims),r%pars(no_of_sims,np))
      allocate(r%sim_gof(no_of_sims),r%ntail(no_of_sims),sample(n),one(1))
      do i=1,no_of_sims
         do j=1,n
            if (rand_uniform() > ptail .and. size(lower) > 0) then
               ix=1+int(rand_uniform()*real(size(lower),dp))
               ix=min(size(lower),max(1,ix))
               sample(j)=lower(ix)
            else
               call base_model%random(one)
               sample(j)=one(1)
            end if
         end do
         w=m
         call w%set_data(sample)
         e=estimate_xmin(w,xmins=xmins,xmax=xmax,distance=distance,candidates=candidates)
         r%xmin(i)=e%xmin
         r%pars(i,:)=e%pars
         r%sim_gof(i)=e%gof
         r%ntail(i)=e%ntail
         if (e%status == 0) r%successful=r%successful+1
      end do
      if (no_of_sims > 0) then
         r%p=real(count(r%sim_gof >= r%gof),dp)/real(no_of_sims,dp)
      end if
   end function

end module powerlaw
