! distr-fortran -- computational translation of the R package distr.
! Copyright (C) 2005-2025 distr authors.
! SPDX-License-Identifier: LGPL-3.0-only
!
! The R package uses S4 classes containing r/d/p/q closures.  This port maps
! the computational semantics to a single recursive value type.  Named and
! composed distributions therefore share one Fortran API without an R runtime.
module distr_core
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use distr_kinds, only : dp, pi, sqrt2, sqrt2pi, eps_dp, nan_dp, inf_dp
   use distr_special, only : normal_pdf_std, normal_cdf_std, normal_quantile_std, &
      regularized_gamma_p, regularized_gamma_q, regularized_beta, beta_log_density, &
      gamma_log_density, log_choose, is_integer_value, clamp01, &
      poisson_weighted_beta, poisson_weighted_chisq, noncentral_t_cdf, &
      noncentral_t_density
   use distr_fft, only : real_convolution_fft, real_convolution_power_fft
   use distr_rng, only : rand_uniform, rand_normal, rand_exponential, rand_gamma, &
      rand_poisson, rand_binomial, rand_negative_binomial, rand_geometric, &
      rand_chisq, rand_noncentral_chisq, rand_beta, rand_f, rand_t, rand_cauchy, &
      rand_logistic, rand_weibull, rand_hypergeometric
   implicit none
   private

   integer, parameter :: K_INVALID=0, K_DIRAC=1, K_BINOM=2, K_HYPER=3, K_POIS=4
   integer, parameter :: K_NBINOM=5, K_GEOM=6, K_UNIF=7, K_NORM=8, K_LNORM=9
   integer, parameter :: K_CAUCHY=10, K_F=11, K_T=12, K_CHISQ=13, K_EXP=14
   integer, parameter :: K_LAPLACE=15, K_GAMMA=16, K_BETA=17, K_LOGIS=18
   integer, parameter :: K_WEIBULL=19, K_ARCSINE=20, K_DISCRETE=21, K_GRID=22
   integer, parameter :: K_MIXTURE=30, K_AFFINE=31, K_TRUNC=32, K_MIN=33
   integer, parameter :: K_MAX=34, K_CONV=35, K_EXPTRANS=36, K_ABSTRANS=37
   integer, parameter :: K_HUBER=38, K_POWERTRANS=39, K_LOGTRANS=40

   type, public :: distribution_t
      private
      integer :: kind = K_INVALID
      real(dp) :: par(8) = 0.0_dp
      integer :: ipar(4) = 0
      logical :: discrete = .false.
      real(dp), allocatable :: support(:)
      real(dp), allocatable :: probs(:)
      real(dp), allocatable :: weights(:)
      type(distribution_t), allocatable :: component(:)
   contains
      procedure, public :: density => distribution_density
      procedure, public :: cdf => distribution_cdf
      procedure, public :: sf => distribution_sf
      procedure, public :: logcdf => distribution_logcdf
      procedure, public :: logsf => distribution_logsf
      procedure, public :: cdf_left => distribution_cdf_left
      procedure, public :: quantile => distribution_quantile
      procedure, public :: random => distribution_random
      procedure, public :: mean => distribution_mean
      procedure, public :: variance => distribution_variance
      procedure, public :: sd => distribution_sd
      procedure, public :: raw_moment => distribution_raw_moment
      procedure, public :: central_moment => distribution_central_moment
      procedure, public :: skewness => distribution_skewness
      procedure, public :: excess_kurtosis => distribution_excess_kurtosis
      procedure, public :: is_discrete => distribution_is_discrete
      procedure, public :: label => distribution_label
   end type distribution_t

   public :: dirac_dist, binomial_dist, hypergeometric_dist, poisson_dist
   public :: negative_binomial_dist, geometric_dist, uniform_dist, normal_dist
   public :: lognormal_dist, cauchy_dist, f_dist, student_t_dist, chisq_dist
   public :: exponential_dist, laplace_dist, gamma_dist, beta_dist, logistic_dist
   public :: weibull_dist, arcsine_dist, discrete_dist, lattice_dist, empirical_dist, weighted_empirical_dist
   public :: mixture_dist, lebesgue_mixture_dist, affine_dist, truncate_dist
   public :: minimum_dist, maximum_dist, convolve, convolve_fft, convpow, convpow_fft, compound_dist
   public :: exp_transform, log_transform, sqrt_transform, reciprocal_transform, &
      power_transform, abs_transform, huberize_dist, kde_dist
   public :: density_vec, cdf_vec, sf_vec, logcdf_vec, logsf_vec, quantile_vec
   public :: operator(+), operator(-), operator(*), operator(/)

   interface operator(+)
      module procedure add_dist_dist, add_dist_real, add_real_dist
   end interface
   interface operator(-)
      module procedure subtract_dist_dist, subtract_dist_real, subtract_real_dist, negate_dist
   end interface
   interface operator(*)
      module procedure multiply_dist_real, multiply_real_dist
   end interface
   interface operator(/)
      module procedure divide_dist_real
   end interface

   abstract interface
      function scalar_fun(x) result(y)
         import dp
         real(dp), intent(in) :: x
         real(dp) :: y
      end function scalar_fun
   end interface

contains

   function dirac_dist(location) result(d)
      real(dp), intent(in), optional :: location
      type(distribution_t) :: d
      d%kind=K_DIRAC; d%discrete=.true.; d%par(1)=0.0_dp
      if (present(location)) d%par(1)=location
   end function dirac_dist

   function binomial_dist(size,prob) result(d)
      integer, intent(in) :: size
      real(dp), intent(in) :: prob
      type(distribution_t) :: d
      if (size < 1 .or. prob < 0.0_dp .or. prob > 1.0_dp) error stop 'invalid binomial parameters'
      d%kind=K_BINOM; d%discrete=.true.; d%ipar(1)=size; d%par(1)=prob
   end function binomial_dist

   function hypergeometric_dist(m,n,k) result(d)
      integer, intent(in) :: m,n,k
      type(distribution_t) :: d
      if (m < 0 .or. n < 0 .or. k < 0 .or. k > m+n) error stop 'invalid hypergeometric parameters'
      d%kind=K_HYPER; d%discrete=.true.; d%ipar(1)=m; d%ipar(2)=n; d%ipar(3)=k
   end function hypergeometric_dist

   function poisson_dist(lambda) result(d)
      real(dp), intent(in) :: lambda
      type(distribution_t) :: d
      if (lambda < 0.0_dp) error stop 'invalid Poisson lambda'
      d%kind=K_POIS; d%discrete=.true.; d%par(1)=lambda
   end function poisson_dist

   function negative_binomial_dist(size,prob) result(d)
      real(dp), intent(in) :: size,prob
      type(distribution_t) :: d
      if (size <= 0.0_dp .or. prob <= 0.0_dp .or. prob > 1.0_dp) error stop 'invalid negative-binomial parameters'
      d%kind=K_NBINOM; d%discrete=.true.; d%par(1)=size; d%par(2)=prob
   end function negative_binomial_dist

   function geometric_dist(prob) result(d)
      real(dp), intent(in) :: prob
      type(distribution_t) :: d
      if (prob <= 0.0_dp .or. prob > 1.0_dp) error stop 'invalid geometric probability'
      d%kind=K_GEOM; d%discrete=.true.; d%par(1)=prob
   end function geometric_dist

   function uniform_dist(lower,upper) result(d)
      real(dp), intent(in), optional :: lower,upper
      type(distribution_t) :: d
      real(dp) :: a,b
      a=0.0_dp; b=1.0_dp
      if (present(lower)) a=lower
      if (present(upper)) b=upper
      if (b <= a) error stop 'uniform upper must exceed lower'
      d%kind=K_UNIF; d%par(1)=a; d%par(2)=b
   end function uniform_dist

   function normal_dist(mean,sd) result(d)
      real(dp), intent(in), optional :: mean,sd
      type(distribution_t) :: d
      real(dp) :: m,s
      m=0.0_dp; s=1.0_dp
      if (present(mean)) m=mean
      if (present(sd)) s=sd
      if (s <= 0.0_dp) error stop 'normal sd must be positive'
      d%kind=K_NORM; d%par(1)=m; d%par(2)=s
   end function normal_dist

   function lognormal_dist(meanlog,sdlog) result(d)
      real(dp), intent(in), optional :: meanlog,sdlog
      type(distribution_t) :: d
      real(dp) :: m,s
      m=0.0_dp; s=1.0_dp
      if (present(meanlog)) m=meanlog
      if (present(sdlog)) s=sdlog
      if (s <= 0.0_dp) error stop 'lognormal sdlog must be positive'
      d%kind=K_LNORM; d%par(1)=m; d%par(2)=s
   end function lognormal_dist

   function cauchy_dist(location,scale) result(d)
      real(dp), intent(in), optional :: location,scale
      type(distribution_t) :: d
      real(dp) :: m,s
      m=0.0_dp; s=1.0_dp
      if (present(location)) m=location
      if (present(scale)) s=scale
      if (s <= 0.0_dp) error stop 'Cauchy scale must be positive'
      d%kind=K_CAUCHY; d%par(1)=m; d%par(2)=s
   end function cauchy_dist

   function f_dist(df1,df2,ncp) result(d)
      real(dp), intent(in) :: df1,df2
      real(dp), intent(in), optional :: ncp
      type(distribution_t) :: d
      real(dp) :: nc
      nc=0.0_dp; if (present(ncp)) nc=ncp
      if (df1 <= 0.0_dp .or. df2 <= 0.0_dp .or. nc < 0.0_dp) error stop 'invalid F parameters'
      d%kind=K_F; d%par(1)=df1; d%par(2)=df2; d%par(3)=nc
   end function f_dist

   function student_t_dist(df,ncp) result(d)
      real(dp), intent(in) :: df
      real(dp), intent(in), optional :: ncp
      type(distribution_t) :: d
      real(dp) :: nc
      nc=0.0_dp; if (present(ncp)) nc=ncp
      if (df <= 0.0_dp) error stop 'student-t df must be positive'
      d%kind=K_T; d%par(1)=df; d%par(2)=nc
   end function student_t_dist

   function chisq_dist(df,ncp) result(d)
      real(dp), intent(in) :: df
      real(dp), intent(in), optional :: ncp
      type(distribution_t) :: d
      real(dp) :: nc
      nc=0.0_dp; if (present(ncp)) nc=ncp
      if (df <= 0.0_dp .or. nc < 0.0_dp) error stop 'invalid chi-square parameters'
      d%kind=K_CHISQ; d%par(1)=df; d%par(2)=nc
   end function chisq_dist

   function exponential_dist(rate) result(d)
      real(dp), intent(in), optional :: rate
      type(distribution_t) :: d
      real(dp) :: r
      r=1.0_dp; if (present(rate)) r=rate
      if (r <= 0.0_dp) error stop 'exponential rate must be positive'
      d%kind=K_EXP; d%par(1)=r
   end function exponential_dist

   function laplace_dist(rate) result(d)
      real(dp), intent(in), optional :: rate
      type(distribution_t) :: d
      real(dp) :: r
      r=1.0_dp; if (present(rate)) r=rate
      if (r <= 0.0_dp) error stop 'Laplace rate must be positive'
      d%kind=K_LAPLACE; d%par(1)=r
   end function laplace_dist

   function gamma_dist(shape,scale) result(d)
      real(dp), intent(in) :: shape
      real(dp), intent(in), optional :: scale
      type(distribution_t) :: d
      real(dp) :: s
      s=1.0_dp; if (present(scale)) s=scale
      if (shape <= 0.0_dp .or. s <= 0.0_dp) error stop 'invalid gamma parameters'
      d%kind=K_GAMMA; d%par(1)=shape; d%par(2)=s
   end function gamma_dist

   function beta_dist(shape1,shape2,ncp) result(d)
      real(dp), intent(in) :: shape1,shape2
      real(dp), intent(in), optional :: ncp
      type(distribution_t) :: d
      real(dp) :: nc
      nc=0.0_dp; if (present(ncp)) nc=ncp
      if (shape1 <= 0.0_dp .or. shape2 <= 0.0_dp .or. nc < 0.0_dp) error stop 'invalid beta parameters'
      d%kind=K_BETA; d%par(1)=shape1; d%par(2)=shape2; d%par(3)=nc
   end function beta_dist

   function logistic_dist(location,scale) result(d)
      real(dp), intent(in), optional :: location,scale
      type(distribution_t) :: d
      real(dp) :: m,s
      m=0.0_dp; s=1.0_dp
      if (present(location)) m=location
      if (present(scale)) s=scale
      if (s <= 0.0_dp) error stop 'logistic scale must be positive'
      d%kind=K_LOGIS; d%par(1)=m; d%par(2)=s
   end function logistic_dist

   function weibull_dist(shape,scale) result(d)
      real(dp), intent(in) :: shape
      real(dp), intent(in), optional :: scale
      type(distribution_t) :: d
      real(dp) :: s
      s=1.0_dp; if (present(scale)) s=scale
      if (shape <= 0.0_dp .or. s <= 0.0_dp) error stop 'invalid Weibull parameters'
      d%kind=K_WEIBULL; d%par(1)=shape; d%par(2)=s
   end function weibull_dist

   function arcsine_dist() result(d)
      type(distribution_t) :: d
      d%kind=K_ARCSINE
   end function arcsine_dist

   function discrete_dist(support,prob) result(d)
      real(dp), intent(in) :: support(:)
      real(dp), intent(in), optional :: prob(:)
      type(distribution_t) :: d
      real(dp), allocatable :: x(:),p(:)
      integer :: n
      n=size(support)
      if (n < 1) error stop 'discrete support may not be empty'
      allocate(x(n),p(n)); x=support
      if (present(prob)) then
         if (size(prob) /= n) error stop 'support/probability size mismatch'
         if (any(prob < 0.0_dp)) error stop 'negative discrete probability'
         if (sum(prob) <= 0.0_dp) error stop 'zero discrete probability sum'
         p=prob/sum(prob)
      else
         p=1.0_dp/real(n,dp)
      end if
      call sort_collapse(x,p)
      d%kind=K_DISCRETE; d%discrete=.true.
      call move_alloc(x,d%support); call move_alloc(p,d%probs)
   end function discrete_dist

   function lattice_dist(offset,step,prob) result(d)
      real(dp), intent(in) :: offset,step,prob(:)
      type(distribution_t) :: d
      real(dp), allocatable :: x(:)
      integer :: i
      if (step<=0.0_dp) error stop 'lattice step must be positive'
      if (size(prob)<1) error stop 'lattice probabilities may not be empty'
      allocate(x(size(prob)))
      do i=1,size(prob)
         x(i)=offset+real(i-1,dp)*step
      end do
      d=discrete_dist(x,prob)
   end function lattice_dist

   function empirical_dist(data) result(d)
      real(dp), intent(in) :: data(:)
      type(distribution_t) :: d
      d=discrete_dist(data)
   end function empirical_dist

   function weighted_empirical_dist(data,weights) result(d)
      real(dp), intent(in) :: data(:),weights(:)
      type(distribution_t) :: d
      if (size(data) /= size(weights)) error stop 'weighted empirical size mismatch'
      d=discrete_dist(data,weights)
   end function weighted_empirical_dist

   function grid_dist(x,prob) result(d)
      real(dp), intent(in) :: x(:),prob(:)
      type(distribution_t) :: d
      real(dp) :: s
      if (size(x)<2 .or. size(prob)/=size(x)) error stop 'invalid grid distribution arrays'
      if (any(prob<0.0_dp)) error stop 'negative grid probability'
      s=sum(prob); if (s<=0.0_dp) error stop 'zero grid mass'
      d%kind=K_GRID; d%discrete=.false.
      allocate(d%support(size(x)),d%probs(size(prob)))
      d%support=x; d%probs=prob/s
   end function grid_dist

   function mixture_dist(weights,components) result(d)
      real(dp), intent(in) :: weights(:)
      type(distribution_t), intent(in) :: components(:)
      type(distribution_t) :: d
      if (size(weights)<1 .or. size(weights)/=size(components)) error stop 'invalid mixture sizes'
      if (any(weights<0.0_dp) .or. sum(weights)<=0.0_dp) error stop 'invalid mixture weights'
      d%kind=K_MIXTURE
      allocate(d%weights(size(weights)),d%component(size(components)))
      d%weights=weights/sum(weights); d%component=components
      d%discrete=all(components%discrete)
   end function mixture_dist

   function lebesgue_mixture_dist(discrete_part,ac_part,discrete_weight) result(d)
      type(distribution_t), intent(in) :: discrete_part,ac_part
      real(dp), intent(in) :: discrete_weight
      type(distribution_t) :: d
      type(distribution_t) :: c(2)
      real(dp) :: w(2)
      if (discrete_weight<0.0_dp .or. discrete_weight>1.0_dp) error stop 'invalid discrete weight'
      c(1)=discrete_part; c(2)=ac_part; w=[discrete_weight,1.0_dp-discrete_weight]
      d=mixture_dist(w,c)
      d%discrete=.false.
   end function lebesgue_mixture_dist

   recursive function affine_dist(base,scale,shift) result(d)
      type(distribution_t), intent(in) :: base
      real(dp), intent(in), optional :: scale,shift
      type(distribution_t) :: d
      real(dp) :: a,b,x1,x2
      a=1.0_dp; b=0.0_dp
      if (present(scale)) a=scale
      if (present(shift)) b=shift
      if (a==0.0_dp) then
         d=dirac_dist(b)
         return
      end if
      select case(base%kind)
      case(K_DIRAC)
         d=dirac_dist(a*base%par(1)+b); return
      case(K_NORM)
         d=normal_dist(a*base%par(1)+b,abs(a)*base%par(2)); return
      case(K_CAUCHY)
         d=cauchy_dist(a*base%par(1)+b,abs(a)*base%par(2)); return
      case(K_UNIF)
         x1=a*base%par(1)+b; x2=a*base%par(2)+b
         d=uniform_dist(min(x1,x2),max(x1,x2)); return
      case(K_LOGIS)
         d=logistic_dist(a*base%par(1)+b,abs(a)*base%par(2)); return
      case(K_LNORM)
         if (a>0.0_dp .and. b==0.0_dp) then
            d=lognormal_dist(base%par(1)+log(a),base%par(2)); return
         end if
      case(K_GAMMA)
         if (a>0.0_dp .and. b==0.0_dp) then
            d=gamma_dist(base%par(1),a*base%par(2)); return
         end if
      case(K_EXP)
         if (a>0.0_dp .and. b==0.0_dp) then
            d=exponential_dist(base%par(1)/a); return
         end if
      case(K_LAPLACE)
         if (b==0.0_dp) then
            d=laplace_dist(base%par(1)/abs(a)); return
         end if
      case(K_WEIBULL)
         if (a>0.0_dp .and. b==0.0_dp) then
            d=weibull_dist(base%par(1),a*base%par(2)); return
         end if
      case(K_AFFINE)
         d=affine_dist(base%component(1),a*base%par(1),a*base%par(2)+b); return
      end select
      d%kind=K_AFFINE; d%par(1)=a; d%par(2)=b; d%discrete=base%discrete
      allocate(d%component(1)); d%component(1)=base
   end function affine_dist

   function truncate_dist(base,lower,upper) result(d)
      type(distribution_t), intent(in) :: base
      real(dp), intent(in), optional :: lower,upper
      type(distribution_t) :: d
      real(dp) :: lo,hi,z
      lo=inf_dp(-1); hi=inf_dp()
      if (present(lower)) lo=lower
      if (present(upper)) hi=upper
      if (hi<lo) error stop 'truncation upper below lower'
      z=raw_cdf(base,hi)-cdf_left_raw(base,lo)
      if (z<=0.0_dp) error stop 'zero probability in truncation interval'
      d%kind=K_TRUNC; d%par(1)=lo; d%par(2)=hi; d%par(3)=z; d%discrete=base%discrete
      allocate(d%component(1)); d%component(1)=base
   end function truncate_dist

   function minimum_dist(a,b) result(d)
      type(distribution_t), intent(in) :: a,b
      type(distribution_t) :: d
      d%kind=K_MIN; d%discrete=a%discrete.and.b%discrete
      allocate(d%component(2)); d%component(1)=a; d%component(2)=b
   end function minimum_dist

   function maximum_dist(a,b) result(d)
      type(distribution_t), intent(in) :: a,b
      type(distribution_t) :: d
      d%kind=K_MAX; d%discrete=a%discrete.and.b%discrete
      allocate(d%component(2)); d%component(1)=a; d%component(2)=b
   end function maximum_dist

   function convolve(a,b) result(d)
      type(distribution_t), intent(in) :: a,b
      type(distribution_t) :: d
      real(dp) :: tol
      tol=64.0_dp*eps_dp
      if (a%kind==K_DIRAC) then
         d=affine_dist(b,1.0_dp,a%par(1))
      else if (b%kind==K_DIRAC) then
         d=affine_dist(a,1.0_dp,b%par(1))
      else if (a%kind==K_NORM .and. b%kind==K_NORM) then
         d=normal_dist(a%par(1)+b%par(1),sqrt(a%par(2)**2+b%par(2)**2))
      else if (a%kind==K_POIS .and. b%kind==K_POIS) then
         d=poisson_dist(a%par(1)+b%par(1))
      else if (a%kind==K_CAUCHY .and. b%kind==K_CAUCHY) then
         d=cauchy_dist(a%par(1)+b%par(1),a%par(2)+b%par(2))
      else if (a%kind==K_BINOM .and. b%kind==K_BINOM .and. &
               abs(a%par(1)-b%par(1))<=tol*max(1.0_dp,abs(a%par(1)),abs(b%par(1)))) then
         d=binomial_dist(a%ipar(1)+b%ipar(1),0.5_dp*(a%par(1)+b%par(1)))
      else if (a%kind==K_NBINOM .and. b%kind==K_NBINOM .and. &
               abs(a%par(2)-b%par(2))<=tol*max(1.0_dp,abs(a%par(2)),abs(b%par(2)))) then
         d=negative_binomial_dist(a%par(1)+b%par(1),0.5_dp*(a%par(2)+b%par(2)))
      else if (a%kind==K_GAMMA .and. b%kind==K_GAMMA .and. &
               abs(a%par(2)-b%par(2))<=tol*max(1.0_dp,abs(a%par(2)),abs(b%par(2)))) then
         d=gamma_dist(a%par(1)+b%par(1),0.5_dp*(a%par(2)+b%par(2)))
      else if (a%kind==K_EXP .and. b%kind==K_EXP .and. &
               abs(a%par(1)-b%par(1))<=tol*max(1.0_dp,abs(a%par(1)),abs(b%par(1)))) then
         d=gamma_dist(2.0_dp,1.0_dp/(0.5_dp*(a%par(1)+b%par(1))))
      else if (a%kind==K_GAMMA .and. b%kind==K_EXP .and. &
               abs(a%par(2)-1.0_dp/b%par(1))<=tol*max(1.0_dp,abs(a%par(2)))) then
         d=gamma_dist(a%par(1)+1.0_dp,a%par(2))
      else if (a%kind==K_EXP .and. b%kind==K_GAMMA .and. &
               abs(b%par(2)-1.0_dp/a%par(1))<=tol*max(1.0_dp,abs(b%par(2)))) then
         d=gamma_dist(b%par(1)+1.0_dp,b%par(2))
      else if (a%kind==K_CHISQ .and. b%kind==K_CHISQ .and. &
               a%par(2)==0.0_dp .and. b%par(2)==0.0_dp) then
         d=chisq_dist(a%par(1)+b%par(1))
      else
         d%kind=K_CONV; d%discrete=a%discrete.and.b%discrete
         allocate(d%component(2)); d%component(1)=a; d%component(2)=b
      end if
   end function convolve

   function convolve_fft(a,b,grid_points,tail_prob) result(d)
      type(distribution_t), intent(in) :: a,b
      integer, intent(in), optional :: grid_points
      real(dp), intent(in), optional :: tail_prob
      type(distribution_t) :: d,closed
      real(dp), allocatable :: pa(:),pb(:),pc(:),x(:),xa(:),xb(:)
      real(dp) :: alo,ahi,blo,bhi,h,ha,hb,eps,sa,sb
      integer :: m,na,nb,i
      logical :: lattice_a,lattice_b

      closed=convolve(a,b)
      if (closed%kind/=K_CONV) then
         d=closed
         return
      end if

      eps=1.0e-10_dp; if (present(tail_prob)) eps=max(1.0e-15_dp,min(1.0e-2_dp,tail_prob))
      if (a%discrete .and. b%discrete) then
         call support_values(a,sqrt(eps),xa,pa)
         call support_values(b,sqrt(eps),xb,pb)
         call lattice_spacing(xa,lattice_a,ha)
         call lattice_spacing(xb,lattice_b,hb)
         if (size(xa)==1) then
            d=affine_dist(b,1.0_dp,xa(1)); return
         else if (size(xb)==1) then
            d=affine_dist(a,1.0_dp,xb(1)); return
         else if (lattice_a .and. lattice_b .and. &
                  abs(ha-hb)<=128.0_dp*eps_dp*max(1.0_dp,abs(ha),abs(hb))) then
            h=0.5_dp*(ha+hb)
            call real_convolution_fft(pa,pb,pc)
            where (pc<0.0_dp) pc=0.0_dp
            pc=pc/sum(pc)
            allocate(x(size(pc)))
            do i=1,size(pc)
               x(i)=xa(1)+xb(1)+real(i-1,dp)*h
            end do
            d=discrete_dist(x,pc)
            return
         else
            d=closed
            return
         end if
      else if (a%discrete .or. b%discrete) then
         d=closed
         return
      end if

      m=2048; if (present(grid_points)) m=max(32,grid_points)
      call finite_probability_bounds(a,0.5_dp*eps,alo,ahi)
      call finite_probability_bounds(b,0.5_dp*eps,blo,bhi)
      h=max(ahi-alo,bhi-blo)/real(m,dp)
      if (.not.ieee_is_finite(h) .or. h<=0.0_dp) error stop 'invalid FFT convolution grid width'
      na=max(1,int(ceiling((ahi-alo)/h)))
      nb=max(1,int(ceiling((bhi-blo)/h)))
      allocate(pa(na),pb(nb))
      do i=1,na
         pa(i)=max(0.0_dp,raw_cdf(a,alo+real(i,dp)*h)-raw_cdf(a,alo+real(i-1,dp)*h))
      end do
      do i=1,nb
         pb(i)=max(0.0_dp,raw_cdf(b,blo+real(i,dp)*h)-raw_cdf(b,blo+real(i-1,dp)*h))
      end do
      sa=sum(pa); sb=sum(pb)
      if (sa<=0.0_dp .or. sb<=0.0_dp) error stop 'FFT convolution captured zero probability'
      pa=pa/sa; pb=pb/sb
      call real_convolution_fft(pa,pb,pc)
      where (pc<0.0_dp) pc=0.0_dp
      if (sum(pc)<=0.0_dp) error stop 'FFT convolution produced zero mass'
      pc=pc/sum(pc)
      allocate(x(size(pc)))
      do i=1,size(pc)
         x(i)=alo+blo+real(i,dp)*h
      end do
      d=grid_dist(x,pc)
   end function convolve_fft

   recursive function convpow(base,n,grid_points,tail_prob) result(d)
      type(distribution_t), intent(in) :: base
      integer, intent(in) :: n
      integer, intent(in), optional :: grid_points
      real(dp), intent(in), optional :: tail_prob
      type(distribution_t) :: d,tmp
      real(dp) :: ep
      if (n<0) error stop 'convpow exponent must be nonnegative'
      ep=1.0e-10_dp; if (present(tail_prob)) ep=tail_prob
      if (n==0) then
         d=dirac_dist(0.0_dp)
      else if (n==1) then
         d=base
      else if (base%kind==K_DIRAC) then
         d=dirac_dist(real(n,dp)*base%par(1))
      else if (base%kind==K_NORM) then
         d=normal_dist(real(n,dp)*base%par(1),sqrt(real(n,dp))*base%par(2))
      else if (base%kind==K_POIS) then
         d=poisson_dist(real(n,dp)*base%par(1))
      else if (base%kind==K_CAUCHY) then
         d=cauchy_dist(real(n,dp)*base%par(1),real(n,dp)*base%par(2))
      else if (base%kind==K_GAMMA) then
         d=gamma_dist(real(n,dp)*base%par(1),base%par(2))
      else if (base%kind==K_EXP) then
         d=gamma_dist(real(n,dp),1.0_dp/base%par(1))
      else if (base%kind==K_CHISQ .and. base%par(2)==0.0_dp) then
         d=chisq_dist(real(n,dp)*base%par(1))
      else if (base%kind==K_BINOM) then
         d=binomial_dist(n*base%ipar(1),base%par(1))
      else if (base%kind==K_NBINOM) then
         d=negative_binomial_dist(real(n,dp)*base%par(1),base%par(2))
      else if (.not.base%discrete .and. present(grid_points)) then
         d=convpow_fft(base,n,grid_points,ep)
      else if (mod(n,2)==0) then
         tmp=convpow(base,n/2)
         d=convolve(tmp,tmp)
      else
         tmp=convpow(base,n-1)
         d=convolve(tmp,base)
      end if
   end function convpow

   function convpow_fft(base,n,grid_points,tail_prob) result(d)
      type(distribution_t), intent(in) :: base
      integer, intent(in) :: n
      integer, intent(in), optional :: grid_points
      real(dp), intent(in), optional :: tail_prob
      type(distribution_t) :: d
      real(dp), allocatable :: p(:),pc(:),x(:),xb(:)
      real(dp) :: lo,hi,h,eps,s
      integer :: m,i
      logical :: is_lattice

      if (n<0) error stop 'convpow_fft exponent must be nonnegative'
      if (n==0) then
         d=dirac_dist(0.0_dp); return
      else if (n==1) then
         d=base; return
      end if

      select case(base%kind)
      case(K_DIRAC,K_NORM,K_POIS,K_CAUCHY,K_GAMMA,K_EXP,K_BINOM,K_NBINOM)
         d=convpow(base,n); return
      case(K_CHISQ)
         if (base%par(2)==0.0_dp) then; d=convpow(base,n); return; end if
      end select

      eps=1.0e-10_dp; if (present(tail_prob)) eps=max(1.0e-15_dp,min(1.0e-2_dp,tail_prob))
      if (base%discrete) then
         call support_values(base,eps,xb,p)
         call lattice_spacing(xb,is_lattice,h)
         if (.not.is_lattice .or. size(xb)<2) then
            d=convpow(base,n); return
         end if
         call real_convolution_power_fft(p,n,pc)
         where (pc<0.0_dp) pc=0.0_dp
         pc=pc/sum(pc)
         allocate(x(size(pc)))
         do i=1,size(pc)
            x(i)=real(n,dp)*xb(1)+real(i-1,dp)*h
         end do
         d=discrete_dist(x,pc)
         return
      end if

      m=2048; if (present(grid_points)) m=max(32,grid_points)
      call finite_probability_bounds(base,eps,lo,hi)
      h=(hi-lo)/real(m,dp)
      if (.not.ieee_is_finite(h) .or. h<=0.0_dp) error stop 'invalid FFT convolution-power grid width'
      allocate(p(m))
      do i=1,m
         p(i)=max(0.0_dp,raw_cdf(base,lo+real(i,dp)*h)-raw_cdf(base,lo+real(i-1,dp)*h))
      end do
      s=sum(p)
      if (s<=0.0_dp) error stop 'FFT convolution power captured zero probability'
      p=p/s
      call real_convolution_power_fft(p,n,pc)
      where (pc<0.0_dp) pc=0.0_dp
      pc=pc/sum(pc)
      allocate(x(size(pc)))
      do i=1,size(pc)
         x(i)=real(n,dp)*(lo+0.5_dp*h)+real(i-1,dp)*h
      end do
      d=grid_dist(x,pc)
   end function convpow_fft

   function compound_dist(count_dist,summand_dist,tail_prob) result(d)
      type(distribution_t), intent(in) :: count_dist,summand_dist
      real(dp), intent(in), optional :: tail_prob
      type(distribution_t) :: d
      type(distribution_t), allocatable :: comps(:)
      real(dp), allocatable :: ns(:),pr(:),w(:)
      real(dp) :: ep
      integer :: i,n
      if (.not.count_dist%discrete) error stop 'compound count distribution must be discrete'
      ep=1.0e-10_dp; if (present(tail_prob)) ep=tail_prob
      call support_values(count_dist,ep,ns,pr)
      if (size(ns)==0) error stop 'empty count support'
      if (any(ns<0.0_dp) .or. any(.not.is_integer_value(ns))) error stop 'compound counts must be nonnegative integers'
      allocate(comps(size(ns)),w(size(ns)))
      do i=1,size(ns)
         n=int(anint(ns(i)))
         comps(i)=convpow(summand_dist,n)
         w(i)=pr(i)
      end do
      d=mixture_dist(w,comps)
   end function compound_dist

   function exp_transform(base) result(d)
      type(distribution_t), intent(in) :: base
      type(distribution_t) :: d
      d%kind=K_EXPTRANS; d%discrete=base%discrete
      allocate(d%component(1)); d%component(1)=base
   end function exp_transform

   function log_transform(base) result(d)
      type(distribution_t), intent(in) :: base
      type(distribution_t) :: d
      real(dp) :: lo,hi
      call theoretical_bounds(base,lo,hi)
      if (hi<=0.0_dp) error stop 'log transform requires positive support'
      if (base%discrete .and. lo<=0.0_dp) then
         if (raw_cdf(base,0.0_dp)>0.0_dp) error stop 'log transform cannot map discrete mass at nonpositive values'
      end if
      if (lo<0.0_dp) error stop 'log transform requires nonnegative support'
      d%kind=K_LOGTRANS; d%discrete=base%discrete
      allocate(d%component(1)); d%component(1)=base
   end function log_transform

   function sqrt_transform(base) result(d)
      type(distribution_t), intent(in) :: base
      type(distribution_t) :: d
      d=power_transform(base,0.5_dp)
   end function sqrt_transform

   function reciprocal_transform(base) result(d)
      type(distribution_t), intent(in) :: base
      type(distribution_t) :: d
      d=power_transform(base,-1.0_dp)
   end function reciprocal_transform

   function power_transform(base,power) result(d)
      type(distribution_t), intent(in) :: base
      real(dp), intent(in) :: power
      type(distribution_t) :: d
      real(dp) :: lo,hi
      integer :: n
      logical :: isint

      if (.not.ieee_is_finite(power)) error stop 'power transform exponent must be finite'
      if (power==0.0_dp) then
         d=dirac_dist(1.0_dp)
         return
      end if
      call theoretical_bounds(base,lo,hi)
      isint=is_integer_value(power)
      if (power<0.0_dp) then
         if (lo<0.0_dp .or. raw_cdf(base,0.0_dp)>0.0_dp) &
            error stop 'negative power transform requires positive support with no mass at zero'
         d%ipar(1)=1
      else if (lo>=0.0_dp) then
         d%ipar(1)=1
      else if (isint) then
         n=int(anint(power))
         if (mod(abs(n),2)==1) then
            d%ipar(1)=2
         else
            d%ipar(1)=3
         end if
      else
         error stop 'noninteger power transform requires nonnegative support'
      end if
      d%kind=K_POWERTRANS; d%par(1)=power; d%discrete=base%discrete
      allocate(d%component(1)); d%component(1)=base
   end function power_transform

   function abs_transform(base) result(d)
      type(distribution_t), intent(in) :: base
      type(distribution_t) :: d
      d%kind=K_ABSTRANS; d%discrete=base%discrete
      allocate(d%component(1)); d%component(1)=base
   end function abs_transform

   function huberize_dist(base,lower,upper) result(d)
      type(distribution_t), intent(in) :: base
      real(dp), intent(in) :: lower,upper
      type(distribution_t) :: d
      if (upper<lower) error stop 'Huberize upper below lower'
      d%kind=K_HUBER; d%discrete=base%discrete; d%par(1)=lower; d%par(2)=upper
      allocate(d%component(1)); d%component(1)=base
   end function huberize_dist

   function kde_dist(sample,bandwidth,ngrid) result(d)
      real(dp), intent(in) :: sample(:)
      real(dp), intent(in), optional :: bandwidth
      integer, intent(in), optional :: ngrid
      type(distribution_t) :: d
      integer :: n,m,i,j
      real(dp) :: h,s,mu,var,lo,hi,dx,z
      real(dp), allocatable :: x(:),mass(:)
      n=size(sample); if (n<2) error stop 'KDE requires at least two observations'
      mu=sum(sample)/real(n,dp)
      var=sum((sample-mu)**2)/real(max(1,n-1),dp)
      s=sqrt(max(var,tiny(1.0_dp)))
      h=1.06_dp*s*real(n,dp)**(-0.2_dp)
      if (present(bandwidth)) h=bandwidth
      if (h<=0.0_dp) error stop 'KDE bandwidth must be positive'
      m=512; if (present(ngrid)) m=max(32,ngrid)
      lo=minval(sample)-4.0_dp*h; hi=maxval(sample)+4.0_dp*h
      dx=(hi-lo)/real(m,dp)
      allocate(x(m),mass(m)); mass=0.0_dp
      do i=1,m
         x(i)=lo+(real(i,dp)-0.5_dp)*dx
         do j=1,n
            z=(x(i)-sample(j))/h
            mass(i)=mass(i)+normal_pdf_std(z)/(h*real(n,dp))*dx
         end do
      end do
      d=grid_dist(x,mass)
   end function kde_dist

   recursive real(dp) function distribution_density(self,x,log_value) result(v)
      class(distribution_t), intent(in) :: self
      real(dp), intent(in) :: x
      logical, intent(in), optional :: log_value
      logical :: lg
      lg=.false.; if (present(log_value)) lg=log_value
      if (lg) then
         v=raw_log_density(self,x)
      else
         v=raw_density(self,x)
      end if
   end function distribution_density

   recursive real(dp) function distribution_cdf(self,x,lower_tail,log_p) result(v)
      class(distribution_t), intent(in) :: self
      real(dp), intent(in) :: x
      logical, intent(in), optional :: lower_tail,log_p
      logical :: lt,lp
      lt=.true.; lp=.false.
      if (present(lower_tail)) lt=lower_tail
      if (present(log_p)) lp=log_p
      if (lt) then
         v=raw_cdf(self,x)
      else
         v=raw_sf(self,x)
      end if
      if (lp) v=log_probability(v)
   end function distribution_cdf

   recursive real(dp) function distribution_sf(self,x,log_p) result(v)
      class(distribution_t), intent(in) :: self
      real(dp), intent(in) :: x
      logical, intent(in), optional :: log_p
      logical :: lp
      lp=.false.; if (present(log_p)) lp=log_p
      v=raw_sf(self,x)
      if (lp) v=log_probability(v)
   end function distribution_sf

   recursive real(dp) function distribution_logcdf(self,x) result(v)
      class(distribution_t), intent(in) :: self
      real(dp), intent(in) :: x
      v=log_probability(raw_cdf(self,x))
   end function distribution_logcdf

   recursive real(dp) function distribution_logsf(self,x) result(v)
      class(distribution_t), intent(in) :: self
      real(dp), intent(in) :: x
      v=log_probability(raw_sf(self,x))
   end function distribution_logsf

   recursive real(dp) function distribution_cdf_left(self,x,lower_tail,log_p) result(v)
      class(distribution_t), intent(in) :: self
      real(dp), intent(in) :: x
      logical, intent(in), optional :: lower_tail,log_p
      logical :: lt,lp
      lt=.true.; lp=.false.
      if (present(lower_tail)) lt=lower_tail
      if (present(log_p)) lp=log_p
      if (lt) then
         v=cdf_left_raw(self,x)
      else
         v=raw_sf(self,x)
         if (self%discrete) v=v+raw_density(self,x)
         v=clamp01(v)
      end if
      if (lp) v=log_probability(v)
   end function distribution_cdf_left

   recursive real(dp) function distribution_quantile(self,p,lower_tail,log_p) result(x)
      class(distribution_t), intent(in) :: self
      real(dp), intent(in) :: p
      logical, intent(in), optional :: lower_tail,log_p
      logical :: lt,lp
      real(dp) :: pp
      lt=.true.; lp=.false.
      if (present(lower_tail)) lt=lower_tail
      if (present(log_p)) lp=log_p
      if (lp) then
         if (p>0.0_dp) then; x=nan_dp(); return; end if
         pp=exp(p)
      else
         pp=p
      end if
      if (pp<0.0_dp .or. pp>1.0_dp) then; x=nan_dp(); return; end if
      if (lt) then
         x=raw_quantile(self,pp)
      else
         x=raw_quantile_upper(self,pp)
      end if
   end function distribution_quantile

   function distribution_random(self,n) result(x)
      class(distribution_t), intent(in) :: self
      integer, intent(in) :: n
      real(dp), allocatable :: x(:)
      integer :: i
      if (n<0) error stop 'random sample size must be nonnegative'
      allocate(x(n))
      do i=1,n
         x(i)=random_one(self)
      end do
   end function distribution_random

   logical function distribution_is_discrete(self) result(v)
      class(distribution_t), intent(in) :: self
      v=self%discrete
   end function distribution_is_discrete

   function distribution_label(self) result(s)
      class(distribution_t), intent(in) :: self
      character(len=:), allocatable :: s
      select case(self%kind)
      case(K_DIRAC); s='Dirac'
      case(K_BINOM); s='Binomial'
      case(K_HYPER); s='Hypergeometric'
      case(K_POIS); s='Poisson'
      case(K_NBINOM); s='Negative binomial'
      case(K_GEOM); s='Geometric'
      case(K_UNIF); s='Uniform'
      case(K_NORM); s='Normal'
      case(K_LNORM); s='Lognormal'
      case(K_CAUCHY); s='Cauchy'
      case(K_F); s='F'
      case(K_T); s='Student t'
      case(K_CHISQ); s='Chi-square'
      case(K_EXP); s='Exponential'
      case(K_LAPLACE); s='Double exponential (Laplace)'
      case(K_GAMMA); s='Gamma'
      case(K_BETA); s='Beta'
      case(K_LOGIS); s='Logistic'
      case(K_WEIBULL); s='Weibull'
      case(K_ARCSINE); s='Arcsine'
      case(K_DISCRETE); s='DiscreteDistribution'
      case(K_GRID); s='Grid continuous distribution'
      case(K_MIXTURE); s='Mixture'
      case(K_AFFINE); s='Affine transform'
      case(K_TRUNC); s='Truncated distribution'
      case(K_MIN); s='Minimum'
      case(K_MAX); s='Maximum'
      case(K_CONV); s='Convolution'
      case(K_EXPTRANS); s='Exponential transform'
      case(K_LOGTRANS); s='Logarithm transform'
      case(K_POWERTRANS); s='Power transform'
      case(K_ABSTRANS); s='Absolute-value transform'
      case(K_HUBER); s='Huberized distribution'
      case default; s='Invalid distribution'
      end select
   end function distribution_label

   recursive real(dp) function raw_density(self,x) result(v)
      class(distribution_t), intent(in) :: self
      real(dp), intent(in) :: x
      real(dp) :: a,b,c,z,den,df1,df2,nc,d1,d2,f1,f2,fl1,fl2,h
      real(dp), allocatable :: sx(:),sp(:)
      integer :: k,lo,hi,i
      select case(self%kind)
      case(K_DIRAC)
         v=merge(1.0_dp,0.0_dp,nearly_equal(x,self%par(1)))
      case(K_BINOM)
         if (.not.is_integer_value(x)) then; v=0.0_dp; return; end if
         k=int(anint(x)); if (k<0.or.k>self%ipar(1)) then; v=0.0_dp; return; end if
         a=self%par(1)
         if (a==0.0_dp) then; v=merge(1.0_dp,0.0_dp,k==0)
         else if (a==1.0_dp) then; v=merge(1.0_dp,0.0_dp,k==self%ipar(1))
         else
            v=exp(log_choose(self%ipar(1),k)+real(k,dp)*log(a)+real(self%ipar(1)-k,dp)*log(1.0_dp-a))
         end if
      case(K_HYPER)
         if (.not.is_integer_value(x)) then; v=0.0_dp; return; end if
         k=int(anint(x)); lo=max(0,self%ipar(3)-self%ipar(2)); hi=min(self%ipar(3),self%ipar(1))
         if (k<lo.or.k>hi) then; v=0.0_dp; return; end if
         v=exp(log_choose(self%ipar(1),k)+log_choose(self%ipar(2),self%ipar(3)-k)- &
               log_choose(self%ipar(1)+self%ipar(2),self%ipar(3)))
      case(K_POIS)
         if (.not.is_integer_value(x)) then; v=0.0_dp; return; end if
         k=int(anint(x)); if (k<0) then; v=0.0_dp; return; end if
         a=self%par(1)
         if (a==0.0_dp) then; v=merge(1.0_dp,0.0_dp,k==0)
         else; v=exp(-a+real(k,dp)*log(a)-log_gamma(real(k+1,dp))); end if
      case(K_NBINOM)
         if (.not.is_integer_value(x)) then; v=0.0_dp; return; end if
         k=int(anint(x)); if (k<0) then; v=0.0_dp; return; end if
         a=self%par(1); b=self%par(2)
         v=exp(log_gamma(real(k,dp)+a)-log_gamma(a)-log_gamma(real(k+1,dp)) + &
               a*log(b)+real(k,dp)*log(1.0_dp-b))
      case(K_GEOM)
         if (.not.is_integer_value(x)) then; v=0.0_dp; return; end if
         k=int(anint(x)); if (k<0) then; v=0.0_dp; return; end if
         a=self%par(1); v=a*(1.0_dp-a)**k
      case(K_UNIF)
         a=self%par(1); b=self%par(2); v=merge(1.0_dp/(b-a),0.0_dp,x>=a.and.x<=b)
      case(K_NORM)
         z=(x-self%par(1))/self%par(2); v=normal_pdf_std(z)/self%par(2)
      case(K_LNORM)
         if (x<=0.0_dp) then; v=0.0_dp; else; z=(log(x)-self%par(1))/self%par(2); v=normal_pdf_std(z)/(x*self%par(2)); end if
      case(K_CAUCHY)
         z=(x-self%par(1))/self%par(2); v=1.0_dp/(pi*self%par(2)*(1.0_dp+z*z))
      case(K_F)
         if (x<0.0_dp) then; v=0.0_dp; return; end if
         df1=self%par(1); df2=self%par(2); nc=self%par(3)
         z=df1*x/(df1*x+df2)
         den=poisson_weighted_beta(z,0.5_dp*df1,0.5_dp*df2,nc,.true.)
         v=den*df1*df2/(df1*x+df2)**2
      case(K_T)
         v=noncentral_t_density(x,self%par(1),self%par(2))
      case(K_CHISQ)
         v=poisson_weighted_chisq(x,self%par(1),self%par(2),.true.)
      case(K_EXP)
         if (x<0.0_dp) then; v=0.0_dp; else; v=self%par(1)*exp(-self%par(1)*x); end if
      case(K_LAPLACE)
         v=0.5_dp*self%par(1)*exp(-self%par(1)*abs(x))
      case(K_GAMMA)
         if (x<0.0_dp) then; v=0.0_dp; else; v=exp(gamma_log_density(x,self%par(1),self%par(2))); end if
      case(K_BETA)
         if (x < 0.0_dp .or. x > 1.0_dp) then
            v = 0.0_dp
         else
            v = poisson_weighted_beta(x, self%par(1), self%par(2), &
                 self%par(3), .true.)
         end if
      case(K_LOGIS)
         z=(x-self%par(1))/self%par(2)
         if (z>=0.0_dp) then; a=exp(-z); v=a/(self%par(2)*(1.0_dp+a)**2)
         else; a=exp(z); v=a/(self%par(2)*(1.0_dp+a)**2); end if
      case(K_WEIBULL)
         if (x<0.0_dp) then; v=0.0_dp
         else if (x==0.0_dp .and. self%par(1)<1.0_dp) then; v=inf_dp()
         else if (x==0.0_dp .and. self%par(1)>1.0_dp) then; v=0.0_dp
         else; a=self%par(1); b=self%par(2); v=(a/b)*(x/b)**(a-1.0_dp)*exp(-(x/b)**a); end if
      case(K_ARCSINE)
         if (abs(x)>1.0_dp) then; v=0.0_dp
         else if (abs(x)==1.0_dp) then; v=inf_dp()
         else; v=1.0_dp/(pi*sqrt(1.0_dp-x*x)); end if
      case(K_DISCRETE)
         v=0.0_dp
         do i=1,size(self%support)
            if (nearly_equal(x,self%support(i))) then; v=self%probs(i); return; end if
         end do
      case(K_GRID)
         v=grid_density(self,x)
      case(K_MIXTURE)
         v=0.0_dp
         do i=1,size(self%component); v=v+self%weights(i)*raw_density(self%component(i),x); end do
      case(K_AFFINE)
         a=self%par(1); b=self%par(2); z=(x-b)/a
         v=raw_density(self%component(1),z)
         if (.not.self%discrete) v=v/abs(a)
      case(K_TRUNC)
         if (x<self%par(1).or.x>self%par(2)) then; v=0.0_dp
         else; v=raw_density(self%component(1),x)/self%par(3); end if
      case(K_MIN)
         d1=raw_density(self%component(1),x); d2=raw_density(self%component(2),x)
         if (self%discrete) then
            f1=raw_cdf(self%component(1),x); f2=raw_cdf(self%component(2),x)
            fl1=f1-d1; fl2=f2-d2
            v=(1.0_dp-fl1)*(1.0_dp-fl2)-(1.0_dp-f1)*(1.0_dp-f2)
         else
            v=d1*(1.0_dp-raw_cdf(self%component(2),x))+d2*(1.0_dp-raw_cdf(self%component(1),x))
         end if
      case(K_MAX)
         d1=raw_density(self%component(1),x); d2=raw_density(self%component(2),x)
         if (self%discrete) then
            f1=raw_cdf(self%component(1),x); f2=raw_cdf(self%component(2),x)
            fl1=f1-d1; fl2=f2-d2; v=f1*f2-fl1*fl2
         else
            v=d1*raw_cdf(self%component(2),x)+d2*raw_cdf(self%component(1),x)
         end if
      case(K_CONV)
         v=convolution_density(self%component(1),self%component(2),x)
      case(K_EXPTRANS)
         if (x<=0.0_dp) then; v=0.0_dp
         else; v=raw_density(self%component(1),log(x)); if (.not.self%discrete) v=v/x; end if
      case(K_LOGTRANS)
         z=exp(x)
         if (.not.ieee_is_finite(z)) then
            v=0.0_dp
         else
            v=raw_density(self%component(1),z)
            if (.not.self%discrete) v=v*z
         end if
      case(K_POWERTRANS)
         a=self%par(1)
         select case(self%ipar(1))
         case(1)
            if (x<0.0_dp .or. (a<0.0_dp .and. x<=0.0_dp)) then
               v=0.0_dp
            else
               z=x**(1.0_dp/a)
               v=raw_density(self%component(1),z)
               if (.not.self%discrete) v=v*power_inverse_jacobian(z,x,a)
            end if
         case(2)
            if (x==0.0_dp) then
               z=0.0_dp
            else
               z=sign(abs(x)**(1.0_dp/a),x)
            end if
            v=raw_density(self%component(1),z)
            if (.not.self%discrete) v=v*power_inverse_jacobian(z,x,a)
         case(3)
            if (x<0.0_dp) then
               v=0.0_dp
            else
               z=x**(1.0_dp/a)
               v=raw_density(self%component(1),z)
               if (z>0.0_dp) v=v+raw_density(self%component(1),-z)
               if (.not.self%discrete) v=v*power_inverse_jacobian(z,x,a)
            end if
         case default
            v=nan_dp()
         end select
      case(K_ABSTRANS)
         if (x<0.0_dp) then; v=0.0_dp
         else if (x==0.0_dp) then; v=raw_density(self%component(1),0.0_dp)
         else; v=raw_density(self%component(1),x)+raw_density(self%component(1),-x); end if
      case(K_HUBER)
         if (x<self%par(1).or.x>self%par(2)) then; v=0.0_dp
         else if (nearly_equal(x,self%par(1))) then; v=raw_cdf(self%component(1),self%par(1))
         else if (nearly_equal(x,self%par(2))) then; v=1.0_dp-cdf_left_raw(self%component(1),self%par(2))
         else; v=raw_density(self%component(1),x); end if
      case default
         v=nan_dp()
      end select
      if (v<0.0_dp .and. v>-1.0e-13_dp) v=0.0_dp
   end function raw_density

   recursive real(dp) function raw_log_density(self,x) result(v)
      class(distribution_t), intent(in) :: self
      real(dp), intent(in) :: x
      real(dp) :: a,b,z,term
      integer :: k,lo,hi,i
      select case(self%kind)
      case(K_DIRAC)
         v=merge(0.0_dp,inf_dp(-1),nearly_equal(x,self%par(1)))
      case(K_BINOM)
         if (.not.is_integer_value(x)) then; v=inf_dp(-1); return; end if
         k=int(anint(x))
         if (k<0 .or. k>self%ipar(1)) then; v=inf_dp(-1); return; end if
         a=self%par(1)
         if (a==0.0_dp) then
            v=merge(0.0_dp,inf_dp(-1),k==0)
         else if (a==1.0_dp) then
            v=merge(0.0_dp,inf_dp(-1),k==self%ipar(1))
         else
            v=log_choose(self%ipar(1),k)+real(k,dp)*log(a)+ &
               real(self%ipar(1)-k,dp)*log1p_safe(-a)
         end if
      case(K_HYPER)
         if (.not.is_integer_value(x)) then; v=inf_dp(-1); return; end if
         k=int(anint(x)); lo=max(0,self%ipar(3)-self%ipar(2)); hi=min(self%ipar(3),self%ipar(1))
         if (k<lo .or. k>hi) then; v=inf_dp(-1); return; end if
         v=log_choose(self%ipar(1),k)+log_choose(self%ipar(2),self%ipar(3)-k)- &
            log_choose(self%ipar(1)+self%ipar(2),self%ipar(3))
      case(K_POIS)
         if (.not.is_integer_value(x)) then; v=inf_dp(-1); return; end if
         k=int(anint(x))
         if (k<0) then; v=inf_dp(-1); return; end if
         a=self%par(1)
         if (a==0.0_dp) then; v=merge(0.0_dp,inf_dp(-1),k==0)
         else; v=-a+real(k,dp)*log(a)-log_gamma(real(k+1,dp)); end if
      case(K_NBINOM)
         if (.not.is_integer_value(x)) then; v=inf_dp(-1); return; end if
         k=int(anint(x)); if (k<0) then; v=inf_dp(-1); return; end if
         a=self%par(1); b=self%par(2)
         v=log_gamma(real(k,dp)+a)-log_gamma(a)-log_gamma(real(k+1,dp))+ &
            a*log(b)+real(k,dp)*log1p_safe(-b)
      case(K_GEOM)
         if (.not.is_integer_value(x)) then; v=inf_dp(-1); return; end if
         k=int(anint(x)); if (k<0) then; v=inf_dp(-1); return; end if
         v=log(self%par(1))+real(k,dp)*log1p_safe(-self%par(1))
      case(K_UNIF)
         if (x>=self%par(1) .and. x<=self%par(2)) then
            v=-log(self%par(2)-self%par(1))
         else
            v=inf_dp(-1)
         end if
      case(K_NORM)
         z=(x-self%par(1))/self%par(2)
         v=-0.5_dp*z*z-log(self%par(2))-log(sqrt2pi)
      case(K_LNORM)
         if (x<=0.0_dp) then; v=inf_dp(-1)
         else
            z=(log(x)-self%par(1))/self%par(2)
            v=-0.5_dp*z*z-log(x)-log(self%par(2))-log(sqrt2pi)
         end if
      case(K_CAUCHY)
         z=(x-self%par(1))/self%par(2)
         v=-log(pi*self%par(2))-log1p_safe(z*z)
      case(K_CHISQ)
         if (self%par(2)==0.0_dp) then
            v=gamma_log_density(x,0.5_dp*self%par(1),2.0_dp)
         else
            v=log_nonnegative(raw_density(self,x))
         end if
      case(K_EXP)
         if (x<0.0_dp) then; v=inf_dp(-1)
         else; v=log(self%par(1))-self%par(1)*x; end if
      case(K_LAPLACE)
         v=log(0.5_dp*self%par(1))-self%par(1)*abs(x)
      case(K_GAMMA)
         v=gamma_log_density(x,self%par(1),self%par(2))
      case(K_BETA)
         if (self%par(3)==0.0_dp) then
            v=beta_log_density(x,self%par(1),self%par(2))
         else
            v=log_nonnegative(raw_density(self,x))
         end if
      case(K_LOGIS)
         z=(x-self%par(1))/self%par(2)
         if (z>=0.0_dp) then
            v=-log(self%par(2))-z-2.0_dp*log1p_safe(exp(-z))
         else
            v=-log(self%par(2))+z-2.0_dp*log1p_safe(exp(z))
         end if
      case(K_WEIBULL)
         if (x<0.0_dp) then
            v=inf_dp(-1)
         else if (x==0.0_dp) then
            if (self%par(1)<1.0_dp) then; v=inf_dp()
            else if (self%par(1)>1.0_dp) then; v=inf_dp(-1)
            else; v=-log(self%par(2)); end if
         else
            a=self%par(1); b=self%par(2)
            v=log(a/b)+(a-1.0_dp)*log(x/b)-(x/b)**a
         end if
      case(K_ARCSINE)
         if (abs(x)>1.0_dp) then; v=inf_dp(-1)
         else if (abs(x)==1.0_dp) then; v=inf_dp()
         else; v=-log(pi)-0.5_dp*log1p_safe(-x*x); end if
      case(K_MIXTURE)
         v=inf_dp(-1)
         do i=1,size(self%component)
            if (self%weights(i)>0.0_dp) then
               term=log(self%weights(i))+raw_log_density(self%component(i),x)
               v=logaddexp(v,term)
            end if
         end do
      case(K_AFFINE)
         a=self%par(1); b=self%par(2); z=(x-b)/a
         v=raw_log_density(self%component(1),z)
         if (.not.self%discrete) v=v-log(abs(a))
      case(K_TRUNC)
         if (x<self%par(1) .or. x>self%par(2)) then
            v=inf_dp(-1)
         else
            v=raw_log_density(self%component(1),x)-log(self%par(3))
         end if
      case(K_EXPTRANS)
         if (x<=0.0_dp) then; v=inf_dp(-1)
         else
            v=raw_log_density(self%component(1),log(x))
            if (.not.self%discrete) v=v-log(x)
         end if
      case(K_LOGTRANS)
         z=exp(x)
         if (.not.ieee_is_finite(z)) then; v=inf_dp(-1)
         else
            v=raw_log_density(self%component(1),z)
            if (.not.self%discrete) v=v+x
         end if
      case default
         v=log_nonnegative(raw_density(self,x))
      end select
   end function raw_log_density

   recursive real(dp) function raw_cdf(self,x) result(v)
      class(distribution_t), intent(in) :: self
      real(dp), intent(in) :: x
      real(dp) :: a,b,z,df1,df2,nc,flo
      integer :: k,lo,hi,i
      select case(self%kind)
      case(K_DIRAC)
         v=merge(1.0_dp,0.0_dp,x>=self%par(1)-8.0_dp*eps_dp*max(1.0_dp,abs(self%par(1))))
      case(K_BINOM)
         k=int(floor(x+8.0_dp*eps_dp)); if (k<0) then; v=0.0_dp
         else if (k>=self%ipar(1)) then; v=1.0_dp
         else if (self%par(1)==0.0_dp) then; v=1.0_dp
         else if (self%par(1)==1.0_dp) then; v=0.0_dp
         else; v=regularized_beta(1.0_dp-self%par(1),real(self%ipar(1)-k,dp),real(k+1,dp)); end if
      case(K_HYPER)
         k=int(floor(x+8.0_dp*eps_dp)); lo=max(0,self%ipar(3)-self%ipar(2)); hi=min(self%ipar(3),self%ipar(1))
         if (k<lo) then; v=0.0_dp
         else if (k>=hi) then; v=1.0_dp
         else; v=0.0_dp; do i=lo,k; v=v+raw_density(self,real(i,dp)); end do; end if
      case(K_POIS)
         k=int(floor(x+8.0_dp*eps_dp)); if (k<0) then; v=0.0_dp
         else if (self%par(1)==0.0_dp) then; v=1.0_dp
         else; v=regularized_gamma_q(real(k+1,dp),self%par(1)); end if
      case(K_NBINOM)
         k=int(floor(x+8.0_dp*eps_dp)); if (k<0) then; v=0.0_dp
         else; v=regularized_beta(self%par(2),self%par(1),real(k+1,dp)); end if
      case(K_GEOM)
         k=int(floor(x+8.0_dp*eps_dp)); if (k<0) then; v=0.0_dp
         else; v=-expm1_safe(real(k+1,dp)*log1p_safe(-self%par(1))); end if
      case(K_UNIF)
         a=self%par(1); b=self%par(2); v=clamp01((x-a)/(b-a))
      case(K_NORM)
         v=normal_cdf_std((x-self%par(1))/self%par(2))
      case(K_LNORM)
         if (x<=0.0_dp) then; v=0.0_dp; else; v=normal_cdf_std((log(x)-self%par(1))/self%par(2)); end if
      case(K_CAUCHY)
         z=(x-self%par(1))/self%par(2)
         if (z<0.0_dp) then
            v=atan(-1.0_dp/z)/pi
         else
            v=0.5_dp+atan(z)/pi
         end if
      case(K_F)
         if (x<=0.0_dp) then; v=0.0_dp; else
            df1=self%par(1); df2=self%par(2); nc=self%par(3); z=df1*x/(df1*x+df2)
            v=poisson_weighted_beta(z,0.5_dp*df1,0.5_dp*df2,nc,.false.)
         end if
      case(K_T)
         v=noncentral_t_cdf(x,self%par(1),self%par(2))
      case(K_CHISQ)
         v=poisson_weighted_chisq(x,self%par(1),self%par(2),.false.)
      case(K_EXP)
         if (x<=0.0_dp) then; v=0.0_dp; else; v=-expm1_safe(-self%par(1)*x); end if
      case(K_LAPLACE)
         if (x<=0.0_dp) then; v=0.5_dp*exp(self%par(1)*x)
         else; v=1.0_dp-0.5_dp*exp(-self%par(1)*x); end if
      case(K_GAMMA)
         if (x<=0.0_dp) then; v=0.0_dp; else; v=regularized_gamma_p(self%par(1),x/self%par(2)); end if
      case(K_BETA)
         if (x<=0.0_dp) then; v=0.0_dp
         else if (x>=1.0_dp) then; v=1.0_dp
         else; v=poisson_weighted_beta(x,self%par(1),self%par(2),self%par(3),.false.); end if
      case(K_LOGIS)
         z=(x-self%par(1))/self%par(2)
         if (z>=0.0_dp) then; v=1.0_dp/(1.0_dp+exp(-z)); else; a=exp(z); v=a/(1.0_dp+a); end if
      case(K_WEIBULL)
         if (x<=0.0_dp) then; v=0.0_dp; else; v=-expm1_safe(-(x/self%par(2))**self%par(1)); end if
      case(K_ARCSINE)
         if (x<=-1.0_dp) then; v=0.0_dp
         else if (x>=1.0_dp) then; v=1.0_dp
         else; v=asin(x)/pi+0.5_dp; end if
      case(K_DISCRETE)
         v=0.0_dp; do i=1,size(self%support); if (self%support(i)<=x+8.0_dp*eps_dp*max(1.0_dp,abs(x))) v=v+self%probs(i); end do
      case(K_GRID)
         v=grid_cdf(self,x)
      case(K_MIXTURE)
         v=0.0_dp; do i=1,size(self%component); v=v+self%weights(i)*raw_cdf(self%component(i),x); end do
      case(K_AFFINE)
         a=self%par(1); b=self%par(2); z=(x-b)/a
         if (a>0.0_dp) then; v=raw_cdf(self%component(1),z)
         else; v=1.0_dp-cdf_left_raw(self%component(1),z); end if
      case(K_TRUNC)
         if (x<self%par(1)) then; v=0.0_dp
         else if (x>=self%par(2)) then; v=1.0_dp
         else; flo=cdf_left_raw(self%component(1),self%par(1)); v=(raw_cdf(self%component(1),x)-flo)/self%par(3); end if
      case(K_MIN)
         a=raw_cdf(self%component(1),x); b=raw_cdf(self%component(2),x)
         v=a+b-a*b
      case(K_MAX)
         v=raw_cdf(self%component(1),x)*raw_cdf(self%component(2),x)
      case(K_CONV)
         v=convolution_cdf(self%component(1),self%component(2),x)
      case(K_EXPTRANS)
         if (x<=0.0_dp) then; v=0.0_dp; else; v=raw_cdf(self%component(1),log(x)); end if
      case(K_LOGTRANS)
         z=exp(x)
         if (ieee_is_finite(z)) then; v=raw_cdf(self%component(1),z); else; v=1.0_dp; end if
      case(K_POWERTRANS)
         a=self%par(1)
         select case(self%ipar(1))
         case(1)
            if (x<0.0_dp .or. (a<0.0_dp .and. x<=0.0_dp)) then
               v=0.0_dp
            else
               z=x**(1.0_dp/a)
               if (a>0.0_dp) then
                  v=raw_cdf(self%component(1),z)
               else
                  v=raw_sf(self%component(1),z)
                  if (self%component(1)%discrete) v=v+raw_density(self%component(1),z)
               end if
            end if
         case(2)
            if (x==0.0_dp) then; z=0.0_dp; else; z=sign(abs(x)**(1.0_dp/a),x); end if
            v=raw_cdf(self%component(1),z)
         case(3)
            if (x<0.0_dp) then; v=0.0_dp
            else; z=x**(1.0_dp/a); v=raw_cdf(self%component(1),z)-cdf_left_raw(self%component(1),-z); end if
         case default
            v=nan_dp()
         end select
      case(K_ABSTRANS)
         if (x<0.0_dp) then; v=0.0_dp; else; v=raw_cdf(self%component(1),x)-cdf_left_raw(self%component(1),-x); end if
      case(K_HUBER)
         if (x<self%par(1)) then; v=0.0_dp
         else if (x>=self%par(2)) then; v=1.0_dp
         else; v=raw_cdf(self%component(1),x); end if
      case default
         v=nan_dp()
      end select
      if (ieee_is_finite(v)) v=clamp01(v)
   end function raw_cdf

   recursive real(dp) function raw_sf(self,x) result(v)
      class(distribution_t), intent(in) :: self
      real(dp), intent(in) :: x
      real(dp) :: a,b,z,df1,df2,nc,flo,f1,f2
      integer :: k,lo,hi,i
      select case(self%kind)
      case(K_DIRAC)
         v=merge(1.0_dp,0.0_dp,x<self%par(1)-8.0_dp*eps_dp*max(1.0_dp,abs(self%par(1))))
      case(K_BINOM)
         k=int(floor(x+8.0_dp*eps_dp))
         if (k<0) then; v=1.0_dp
         else if (k>=self%ipar(1)) then; v=0.0_dp
         else if (self%par(1)==0.0_dp) then; v=0.0_dp
         else if (self%par(1)==1.0_dp) then; v=1.0_dp
         else; v=regularized_beta(self%par(1),real(k+1,dp),real(self%ipar(1)-k,dp)); end if
      case(K_HYPER)
         k=int(floor(x+8.0_dp*eps_dp)); lo=max(0,self%ipar(3)-self%ipar(2)); hi=min(self%ipar(3),self%ipar(1))
         if (k<lo) then; v=1.0_dp
         else if (k>=hi) then; v=0.0_dp
         else; v=0.0_dp; do i=k+1,hi; v=v+raw_density(self,real(i,dp)); end do; end if
      case(K_POIS)
         k=int(floor(x+8.0_dp*eps_dp))
         if (k<0) then; v=1.0_dp
         else if (self%par(1)==0.0_dp) then; v=0.0_dp
         else; v=regularized_gamma_p(real(k+1,dp),self%par(1)); end if
      case(K_NBINOM)
         k=int(floor(x+8.0_dp*eps_dp))
         if (k<0) then; v=1.0_dp
         else; v=regularized_beta(1.0_dp-self%par(2),real(k+1,dp),self%par(1)); end if
      case(K_GEOM)
         k=int(floor(x+8.0_dp*eps_dp))
         if (k<0) then; v=1.0_dp
         else; v=exp(real(k+1,dp)*log1p_safe(-self%par(1))); end if
      case(K_UNIF)
         a=self%par(1); b=self%par(2); v=clamp01((b-x)/(b-a))
      case(K_NORM)
         z=(x-self%par(1))/self%par(2); v=0.5_dp*erfc(z/sqrt2)
      case(K_LNORM)
         if (x<=0.0_dp) then; v=1.0_dp
         else; z=(log(x)-self%par(1))/self%par(2); v=0.5_dp*erfc(z/sqrt2); end if
      case(K_CAUCHY)
         z=(x-self%par(1))/self%par(2)
         if (z>0.0_dp) then; v=atan(1.0_dp/z)/pi
         else; v=0.5_dp-atan(z)/pi; end if
      case(K_F)
         if (x<=0.0_dp) then; v=1.0_dp
         else
            df1=self%par(1); df2=self%par(2); nc=self%par(3); z=df1*x/(df1*x+df2)
            v=poisson_weighted_beta(z,0.5_dp*df1,0.5_dp*df2,nc,.false.,upper_tail=.true.)
         end if
      case(K_T)
         v=noncentral_t_cdf(-x,self%par(1),-self%par(2))
      case(K_CHISQ)
         v=poisson_weighted_chisq(x,self%par(1),self%par(2),.false.,upper_tail=.true.)
      case(K_EXP)
         if (x<0.0_dp) then; v=1.0_dp; else; v=exp(-self%par(1)*x); end if
      case(K_LAPLACE)
         if (x<0.0_dp) then; v=1.0_dp-0.5_dp*exp(self%par(1)*x)
         else; v=0.5_dp*exp(-self%par(1)*x); end if
      case(K_GAMMA)
         if (x<=0.0_dp) then; v=1.0_dp; else; v=regularized_gamma_q(self%par(1),x/self%par(2)); end if
      case(K_BETA)
         if (x<0.0_dp) then; v=1.0_dp
         else if (x>=1.0_dp) then; v=0.0_dp
         else; v=poisson_weighted_beta(x,self%par(1),self%par(2),self%par(3),.false.,upper_tail=.true.); end if
      case(K_LOGIS)
         z=(x-self%par(1))/self%par(2)
         if (z>=0.0_dp) then; a=exp(-z); v=a/(1.0_dp+a)
         else; a=exp(z); v=1.0_dp/(1.0_dp+a); end if
      case(K_WEIBULL)
         if (x<0.0_dp) then; v=1.0_dp; else; v=exp(-(x/self%par(2))**self%par(1)); end if
      case(K_ARCSINE)
         if (x<=-1.0_dp) then; v=1.0_dp
         else if (x>=1.0_dp) then; v=0.0_dp
         else; v=acos(x)/pi; end if
      case(K_DISCRETE)
         v=0.0_dp
         do i=1,size(self%support)
            if (self%support(i)>x+8.0_dp*eps_dp*max(1.0_dp,abs(x))) v=v+self%probs(i)
         end do
      case(K_GRID)
         v=grid_sf(self,x)
      case(K_MIXTURE)
         v=0.0_dp; do i=1,size(self%component); v=v+self%weights(i)*raw_sf(self%component(i),x); end do
      case(K_AFFINE)
         a=self%par(1); b=self%par(2); z=(x-b)/a
         if (a>0.0_dp) then; v=raw_sf(self%component(1),z)
         else; v=cdf_left_raw(self%component(1),z); end if
      case(K_TRUNC)
         if (x<self%par(1)) then; v=1.0_dp
         else if (x>=self%par(2)) then; v=0.0_dp
         else; v=(raw_sf(self%component(1),x)-raw_sf(self%component(1),self%par(2)))/self%par(3); end if
      case(K_MIN)
         v=raw_sf(self%component(1),x)*raw_sf(self%component(2),x)
      case(K_MAX)
         f1=raw_sf(self%component(1),x); f2=raw_sf(self%component(2),x); v=f1+f2-f1*f2
      case(K_CONV)
         v=convolution_sf(self%component(1),self%component(2),x)
      case(K_EXPTRANS)
         if (x<=0.0_dp) then; v=1.0_dp; else; v=raw_sf(self%component(1),log(x)); end if
      case(K_LOGTRANS)
         z=exp(x)
         if (ieee_is_finite(z)) then; v=raw_sf(self%component(1),z); else; v=0.0_dp; end if
      case(K_POWERTRANS)
         a=self%par(1)
         select case(self%ipar(1))
         case(1)
            if (x<0.0_dp .or. (a<0.0_dp .and. x<=0.0_dp)) then
               v=1.0_dp
            else
               z=x**(1.0_dp/a)
               if (a>0.0_dp) then; v=raw_sf(self%component(1),z)
               else; v=cdf_left_raw(self%component(1),z); end if
            end if
         case(2)
            if (x==0.0_dp) then; z=0.0_dp; else; z=sign(abs(x)**(1.0_dp/a),x); end if
            v=raw_sf(self%component(1),z)
         case(3)
            if (x<0.0_dp) then; v=1.0_dp
            else; z=x**(1.0_dp/a); v=raw_sf(self%component(1),z)+cdf_left_raw(self%component(1),-z); end if
         case default
            v=nan_dp()
         end select
      case(K_ABSTRANS)
         if (x<0.0_dp) then; v=1.0_dp
         else; v=raw_sf(self%component(1),x)+cdf_left_raw(self%component(1),-x); end if
      case(K_HUBER)
         if (x<self%par(1)) then; v=1.0_dp
         else if (x>=self%par(2)) then; v=0.0_dp
         else; v=raw_sf(self%component(1),x); end if
      case default
         v=max(0.0_dp,1.0_dp-raw_cdf(self,x))
      end select
      if (ieee_is_finite(v)) v=clamp01(v)
   end function raw_sf

   recursive real(dp) function cdf_left_raw(self,x) result(v)
      class(distribution_t), intent(in) :: self
      real(dp), intent(in) :: x
      real(dp) :: a,b,z,fl1,fl2,flo
      if (self%discrete) then
         v=max(0.0_dp,raw_cdf(self,x)-raw_density(self,x))
         return
      end if
      select case(self%kind)
      case(K_MIXTURE)
         v=mixture_cdf_left(self,x)
      case(K_AFFINE)
         a=self%par(1); b=self%par(2); z=(x-b)/a
         if (a>0.0_dp) then; v=cdf_left_raw(self%component(1),z)
         else; v=raw_sf(self%component(1),z); end if
      case(K_TRUNC)
         if (x<=self%par(1)) then; v=0.0_dp
         else if (x>self%par(2)) then; v=1.0_dp
         else
            flo=cdf_left_raw(self%component(1),self%par(1))
            v=(cdf_left_raw(self%component(1),x)-flo)/self%par(3)
         end if
      case(K_MIN)
         fl1=cdf_left_raw(self%component(1),x); fl2=cdf_left_raw(self%component(2),x)
         v=fl1+fl2-fl1*fl2
      case(K_MAX)
         fl1=cdf_left_raw(self%component(1),x); fl2=cdf_left_raw(self%component(2),x)
         v=fl1*fl2
      case(K_CONV)
         if (self%component(1)%discrete .or. self%component(2)%discrete) then
            v=convolution_cdf_left(self%component(1),self%component(2),x)
         else
            v=raw_cdf(self,x)
         end if
      case(K_EXPTRANS)
         if (x<=0.0_dp) then; v=0.0_dp; else; v=cdf_left_raw(self%component(1),log(x)); end if
      case(K_LOGTRANS)
         z=exp(x)
         if (ieee_is_finite(z)) then; v=cdf_left_raw(self%component(1),z); else; v=1.0_dp; end if
      case(K_POWERTRANS)
         a=self%par(1)
         select case(self%ipar(1))
         case(1)
            if (x<=0.0_dp .and. a<0.0_dp) then
               v=0.0_dp
            else if (x<0.0_dp) then
               v=0.0_dp
            else
               z=x**(1.0_dp/a)
               if (a>0.0_dp) then; v=cdf_left_raw(self%component(1),z)
               else; v=raw_sf(self%component(1),z); end if
            end if
         case(2)
            if (x==0.0_dp) then; z=0.0_dp; else; z=sign(abs(x)**(1.0_dp/a),x); end if
            v=cdf_left_raw(self%component(1),z)
         case(3)
            if (x<=0.0_dp) then; v=0.0_dp
            else
               z=x**(1.0_dp/a)
               v=cdf_left_raw(self%component(1),z)-raw_cdf(self%component(1),-z)
            end if
         case default
            v=nan_dp()
         end select
      case(K_ABSTRANS)
         if (x<=0.0_dp) then; v=0.0_dp
         else; v=cdf_left_raw(self%component(1),x)-raw_cdf(self%component(1),-x); end if
      case(K_HUBER)
         if (x<=self%par(1)) then; v=0.0_dp
         else if (x>self%par(2)) then; v=1.0_dp
         else; v=cdf_left_raw(self%component(1),x); end if
      case default
         v=raw_cdf(self,x)
      end select
      if (ieee_is_finite(v)) v=clamp01(v)
   end function cdf_left_raw

   recursive real(dp) function raw_quantile(self,p) result(x)
      class(distribution_t), intent(in) :: self
      real(dp), intent(in) :: p
      real(dp) :: a,b,flo,target
      real(dp), allocatable :: sx(:),sp(:)
      integer :: i
      if (p<=0.0_dp) then; call theoretical_bounds(self,a,b); x=a; return; end if
      if (p>=1.0_dp) then; call theoretical_bounds(self,a,b); x=b; return; end if
      select case(self%kind)
      case(K_DIRAC); x=self%par(1)
      case(K_BINOM,K_HYPER,K_POIS,K_NBINOM,K_GEOM)
         x=integer_quantile(self,p)
      case(K_UNIF); x=self%par(1)+p*(self%par(2)-self%par(1))
      case(K_NORM); x=self%par(1)+self%par(2)*normal_quantile_std(p)
      case(K_LNORM); x=exp(self%par(1)+self%par(2)*normal_quantile_std(p))
      case(K_CAUCHY); x=self%par(1)+self%par(2)*tan(pi*(p-0.5_dp))
      case(K_EXP); x=-log(1.0_dp-p)/self%par(1)
      case(K_LAPLACE)
         if (p<0.5_dp) then; x=log(2.0_dp*p)/self%par(1)
         else; x=-log(2.0_dp*(1.0_dp-p))/self%par(1); end if
      case(K_LOGIS); x=self%par(1)+self%par(2)*log(p/(1.0_dp-p))
      case(K_WEIBULL); x=self%par(2)*(-log(1.0_dp-p))**(1.0_dp/self%par(1))
      case(K_ARCSINE); x=sin(pi*(p-0.5_dp))
      case(K_DISCRETE)
         target=0.0_dp; x=self%support(size(self%support))
         do i=1,size(self%support); target=target+self%probs(i); if (target>=p) then; x=self%support(i); exit; end if; end do
      case(K_GRID)
         x=grid_quantile(self,p)
      case(K_AFFINE)
         a=self%par(1); b=self%par(2)
         if (a>0.0_dp) then
            x=a*raw_quantile(self%component(1),p)+b
         else if (self%discrete) then
            x=support_quantile(self,p)
         else
            x=a*raw_quantile_upper(self%component(1),p)+b
         end if
      case(K_TRUNC)
         flo=cdf_left_raw(self%component(1),self%par(1)); target=flo+p*self%par(3)
         x=raw_quantile(self%component(1),target); x=min(self%par(2),max(self%par(1),x))
      case(K_EXPTRANS); x=exp(raw_quantile(self%component(1),p))
      case(K_LOGTRANS); x=log(raw_quantile(self%component(1),p))
      case(K_POWERTRANS)
         a=self%par(1)
         if (self%ipar(1)==1 .and. a>0.0_dp) then
            x=raw_quantile(self%component(1),p)**a
         else if (self%ipar(1)==1 .and. a<0.0_dp .and. self%discrete) then
            x=support_quantile(self,p)
         else if (self%ipar(1)==1 .and. a<0.0_dp) then
            x=raw_quantile_upper(self%component(1),p)**a
         else if (self%ipar(1)==2) then
            target=raw_quantile(self%component(1),p)
            x=sign(abs(target)**a,target)
         else if (self%discrete) then
            x=support_quantile(self,p)
         else
            x=continuous_quantile_bisect(self,p)
         end if
      case default
         if (self%discrete) then
            call support_values(self,1.0e-12_dp,sx,sp)
            target=0.0_dp; x=sx(size(sx))
            do i=1,size(sx); target=target+sp(i); if (target>=p) then; x=sx(i); exit; end if; end do
         else
            x=continuous_quantile_bisect(self,p)
         end if
      end select
   end function raw_quantile

   recursive real(dp) function raw_quantile_upper(self,p) result(x)
      class(distribution_t), intent(in) :: self
      real(dp), intent(in) :: p
      real(dp) :: a,b,target
      real(dp), allocatable :: sx(:),sp(:)
      integer :: i

      if (p<=0.0_dp) then; call theoretical_bounds(self,a,b); x=b; return; end if
      if (p>=1.0_dp) then; call theoretical_bounds(self,a,b); x=a; return; end if

      select case(self%kind)
      case(K_DIRAC)
         x=self%par(1)
      case(K_BINOM,K_HYPER,K_POIS,K_NBINOM,K_GEOM)
         x=integer_quantile_upper(self,p)
      case(K_UNIF)
         x=self%par(2)-p*(self%par(2)-self%par(1))
      case(K_NORM)
         x=self%par(1)-self%par(2)*normal_quantile_std(p)
      case(K_LNORM)
         x=exp(self%par(1)-self%par(2)*normal_quantile_std(p))
      case(K_CAUCHY)
         x=self%par(1)+self%par(2)/tan(pi*p)
      case(K_EXP)
         x=-log(p)/self%par(1)
      case(K_LAPLACE)
         if (p<0.5_dp) then; x=-log(2.0_dp*p)/self%par(1)
         else; x=log(2.0_dp*(1.0_dp-p))/self%par(1); end if
      case(K_LOGIS)
         x=self%par(1)+self%par(2)*(log1p_safe(-p)-log(p))
      case(K_WEIBULL)
         x=self%par(2)*(-log(p))**(1.0_dp/self%par(1))
      case(K_ARCSINE)
         x=cos(pi*p)
      case(K_DISCRETE)
         target=0.0_dp
         x=self%support(1)
         do i=size(self%support),1,-1
            if (target<=p) x=self%support(i)
            target=target+self%probs(i)
            if (target>p) exit
         end do
      case(K_GRID)
         x=grid_quantile_upper(self,p)
      case(K_AFFINE)
         a=self%par(1); b=self%par(2)
         if (self%discrete) then
            x=support_quantile_upper(self,p)
         else if (a>0.0_dp) then
            x=a*raw_quantile_upper(self%component(1),p)+b
         else
            x=continuous_quantile_upper_bisect(self,p)
         end if
      case(K_TRUNC,K_MIN,K_MAX,K_CONV,K_ABSTRANS,K_HUBER)
         if (self%discrete) then
            call support_values(self,1.0e-12_dp,sx,sp)
            target=0.0_dp; x=sx(1)
            do i=size(sx),1,-1
               if (target<=p) x=sx(i)
               target=target+sp(i)
               if (target>p) exit
            end do
         else
            x=continuous_quantile_upper_bisect(self,p)
         end if
      case(K_EXPTRANS)
         x=exp(raw_quantile_upper(self%component(1),p))
      case(K_LOGTRANS)
         x=log(raw_quantile_upper(self%component(1),p))
      case(K_POWERTRANS)
         a=self%par(1)
         if (self%discrete .and. (self%ipar(1)==3 .or. a<0.0_dp)) then
            x=support_quantile_upper(self,p)
         else if (self%ipar(1)==1 .and. a>0.0_dp) then
            x=raw_quantile_upper(self%component(1),p)**a
         else if (self%ipar(1)==1 .and. a<0.0_dp) then
            x=raw_quantile(self%component(1),p)**a
         else if (self%ipar(1)==2) then
            target=raw_quantile_upper(self%component(1),p)
            x=sign(abs(target)**a,target)
         else
            x=continuous_quantile_upper_bisect(self,p)
         end if
      case default
         if (self%discrete) then
            call support_values(self,1.0e-12_dp,sx,sp)
            target=0.0_dp; x=sx(1)
            do i=size(sx),1,-1
               if (target<=p) x=sx(i)
               target=target+sp(i)
               if (target>p) exit
            end do
         else
            x=continuous_quantile_upper_bisect(self,p)
         end if
      end select
   end function raw_quantile_upper

   recursive real(dp) function random_one(self) result(x)
      class(distribution_t), intent(in) :: self
      real(dp) :: u,a,b,flo,z
      integer :: i,k
      select case(self%kind)
      case(K_DIRAC); x=self%par(1)
      case(K_BINOM); x=real(rand_binomial(self%ipar(1),self%par(1)),dp)
      case(K_HYPER); x=real(rand_hypergeometric(self%ipar(1),self%ipar(2),self%ipar(3)),dp)
      case(K_POIS); x=real(rand_poisson(self%par(1)),dp)
      case(K_NBINOM); x=real(rand_negative_binomial(self%par(1),self%par(2)),dp)
      case(K_GEOM); x=real(rand_geometric(self%par(1)),dp)
      case(K_UNIF); x=self%par(1)+(self%par(2)-self%par(1))*rand_uniform()
      case(K_NORM); x=self%par(1)+self%par(2)*rand_normal()
      case(K_LNORM); x=exp(self%par(1)+self%par(2)*rand_normal())
      case(K_CAUCHY); x=rand_cauchy(self%par(1),self%par(2))
      case(K_F); x=rand_f(self%par(1),self%par(2),self%par(3))
      case(K_T); x=rand_t(self%par(1),self%par(2))
      case(K_CHISQ); x=rand_noncentral_chisq(self%par(1),self%par(2))
      case(K_EXP); x=rand_exponential(self%par(1))
      case(K_LAPLACE); x=merge(1.0_dp,-1.0_dp,rand_uniform()<0.5_dp)*rand_exponential(self%par(1))
      case(K_GAMMA); x=rand_gamma(self%par(1),self%par(2))
      case(K_BETA); x=rand_beta(self%par(1),self%par(2),self%par(3))
      case(K_LOGIS); x=rand_logistic(self%par(1),self%par(2))
      case(K_WEIBULL); x=rand_weibull(self%par(1),self%par(2))
      case(K_ARCSINE); x=sin(pi*(rand_uniform()-0.5_dp))
      case(K_DISCRETE)
         u=rand_uniform(); z=0.0_dp; x=self%support(size(self%support))
         do i=1,size(self%support); z=z+self%probs(i); if (u<=z) then; x=self%support(i); exit; end if; end do
      case(K_GRID); x=grid_quantile(self,rand_uniform())
      case(K_MIXTURE)
         u=rand_uniform(); z=0.0_dp; k=size(self%component)
         do i=1,size(self%component); z=z+self%weights(i); if (u<=z) then; k=i; exit; end if; end do
         x=random_one(self%component(k))
      case(K_AFFINE); x=self%par(1)*random_one(self%component(1))+self%par(2)
      case(K_TRUNC)
         flo=cdf_left_raw(self%component(1),self%par(1)); x=raw_quantile(self%component(1),flo+rand_uniform()*self%par(3))
         x=min(self%par(2),max(self%par(1),x))
      case(K_MIN); x=min(random_one(self%component(1)),random_one(self%component(2)))
      case(K_MAX); x=max(random_one(self%component(1)),random_one(self%component(2)))
      case(K_CONV); x=random_one(self%component(1))+random_one(self%component(2))
      case(K_EXPTRANS); x=exp(random_one(self%component(1)))
      case(K_LOGTRANS); x=log(random_one(self%component(1)))
      case(K_POWERTRANS)
         z=random_one(self%component(1))
         if (self%ipar(1)==1) then
            x=z**self%par(1)
         else
            x=z**int(anint(self%par(1)))
         end if
      case(K_ABSTRANS); x=abs(random_one(self%component(1)))
      case(K_HUBER); x=min(self%par(2),max(self%par(1),random_one(self%component(1))))
      case default; x=nan_dp()
      end select
   end function random_one

   recursive real(dp) function distribution_mean(self) result(m)
      class(distribution_t), intent(in) :: self
      integer :: i
      real(dp) :: a,b,mu
      select case(self%kind)
      case(K_DIRAC); m=self%par(1)
      case(K_BINOM); m=real(self%ipar(1),dp)*self%par(1)
      case(K_HYPER); m=real(self%ipar(3)*self%ipar(1),dp)/real(self%ipar(1)+self%ipar(2),dp)
      case(K_POIS); m=self%par(1)
      case(K_NBINOM); m=self%par(1)*(1.0_dp-self%par(2))/self%par(2)
      case(K_GEOM); m=(1.0_dp-self%par(1))/self%par(1)
      case(K_UNIF); m=0.5_dp*(self%par(1)+self%par(2))
      case(K_NORM); m=self%par(1)
      case(K_LNORM); m=exp(self%par(1)+0.5_dp*self%par(2)**2)
      case(K_CAUCHY); m=nan_dp()
      case(K_F)
         if (self%par(2)<=2.0_dp) then; m=nan_dp()
         else; m=self%par(2)*(self%par(1)+self%par(3))/(self%par(1)*(self%par(2)-2.0_dp)); end if
      case(K_T)
         if (self%par(1)<=1.0_dp) then; m=nan_dp()
         else
            m = self%par(2)*sqrt(0.5_dp*self%par(1))* &
                 exp(log_gamma(0.5_dp*(self%par(1)-1.0_dp)) - &
                 log_gamma(0.5_dp*self%par(1)))
         end if
      case(K_CHISQ); m=self%par(1)+self%par(2)
      case(K_EXP); m=1.0_dp/self%par(1)
      case(K_LAPLACE); m=0.0_dp
      case(K_GAMMA); m=self%par(1)*self%par(2)
      case(K_BETA)
         if (self%par(3)==0.0_dp) then; m=self%par(1)/(self%par(1)+self%par(2)); else; m=numeric_moment(self,1); end if
      case(K_LOGIS); m=self%par(1)
      case(K_WEIBULL); m=self%par(2)*gamma(1.0_dp+1.0_dp/self%par(1))
      case(K_ARCSINE); m=0.0_dp
      case(K_DISCRETE); m=sum(self%support*self%probs)
      case(K_GRID); m=sum(self%support*self%probs)
      case(K_MIXTURE); m=0.0_dp; do i=1,size(self%component); m=m+self%weights(i)*self%component(i)%mean(); end do
      case(K_AFFINE); m=self%par(1)*self%component(1)%mean()+self%par(2)
      case(K_CONV); m=self%component(1)%mean()+self%component(2)%mean()
      case(K_EXPTRANS,K_LOGTRANS,K_POWERTRANS,K_ABSTRANS,K_TRUNC,K_MIN,K_MAX,K_HUBER); m=numeric_moment(self,1)
      case default; m=nan_dp()
      end select
   end function distribution_mean

   recursive real(dp) function distribution_variance(self) result(v)
      class(distribution_t), intent(in) :: self
      integer :: i
      real(dp) :: m,mi,vi,p,q,n,df1,df2,nc
      select case(self%kind)
      case(K_DIRAC); v=0.0_dp
      case(K_BINOM); p=self%par(1); v=real(self%ipar(1),dp)*p*(1.0_dp-p)
      case(K_HYPER)
         n=real(self%ipar(1)+self%ipar(2),dp); p=real(self%ipar(1),dp)/n
         if (n<=1.0_dp) then; v=0.0_dp; else; v=real(self%ipar(3),dp)*p*(1.0_dp-p)*(n-real(self%ipar(3),dp))/(n-1.0_dp); end if
      case(K_POIS); v=self%par(1)
      case(K_NBINOM); v=self%par(1)*(1.0_dp-self%par(2))/self%par(2)**2
      case(K_GEOM); v=(1.0_dp-self%par(1))/self%par(1)**2
      case(K_UNIF); v=(self%par(2)-self%par(1))**2/12.0_dp
      case(K_NORM); v=self%par(2)**2
      case(K_LNORM); v=(exp(self%par(2)**2)-1.0_dp)*exp(2.0_dp*self%par(1)+self%par(2)**2)
      case(K_CAUCHY); v=nan_dp()
      case(K_F)
         df1=self%par(1); df2=self%par(2); nc=self%par(3)
         if (df2<=4.0_dp) then; v=nan_dp()
         else
            v=2.0_dp*(df2/df1)**2*((df1+nc)**2+(df1+2.0_dp*nc)*(df2-2.0_dp))/((df2-2.0_dp)**2*(df2-4.0_dp))
         end if
      case(K_T)
         if (self%par(1)<=2.0_dp) then; v=nan_dp()
         else
            m=self%mean(); v=self%par(1)*(1.0_dp+self%par(2)**2)/(self%par(1)-2.0_dp)-m*m
         end if
      case(K_CHISQ); v=2.0_dp*(self%par(1)+2.0_dp*self%par(2))
      case(K_EXP); v=1.0_dp/self%par(1)**2
      case(K_LAPLACE); v=2.0_dp/self%par(1)**2
      case(K_GAMMA); v=self%par(1)*self%par(2)**2
      case(K_BETA)
         if (self%par(3)==0.0_dp) then
            p=self%par(1); q=self%par(2); v=p*q/((p+q)**2*(p+q+1.0_dp))
         else; m=self%mean(); v=numeric_moment(self,2)-m*m; end if
      case(K_LOGIS); v=(pi*self%par(2))**2/3.0_dp
      case(K_WEIBULL); v=self%par(2)**2*(gamma(1.0_dp+2.0_dp/self%par(1))-gamma(1.0_dp+1.0_dp/self%par(1))**2)
      case(K_ARCSINE); v=0.5_dp
      case(K_DISCRETE,K_GRID)
         m=self%mean(); v=sum((self%support-m)**2*self%probs)
      case(K_MIXTURE)
         m=self%mean(); v=0.0_dp
         do i = 1, size(self%component)
            mi = self%component(i)%mean()
            vi = self%component(i)%variance()
            v = v + self%weights(i)*(vi + (mi-m)**2)
         end do
      case(K_AFFINE); v=self%par(1)**2*self%component(1)%variance()
      case(K_CONV); v=self%component(1)%variance()+self%component(2)%variance()
      case default
         m=self%mean(); v=numeric_moment(self,2)-m*m
      end select
      if (ieee_is_finite(v) .and. v<0.0_dp .and. v>-1.0e-12_dp) v=0.0_dp
   end function distribution_variance

   real(dp) function distribution_sd(self) result(v)
      class(distribution_t), intent(in) :: self
      v=self%variance()
      if (v>=0.0_dp) v=sqrt(v)
   end function distribution_sd

   real(dp) function distribution_raw_moment(self,order) result(v)
      class(distribution_t), intent(in) :: self
      integer, intent(in) :: order
      real(dp) :: m,var
      if (order<0) error stop 'raw moment order must be nonnegative'
      select case(order)
      case(0)
         v=1.0_dp
      case(1)
         v=self%mean()
      case(2)
         m=self%mean(); var=self%variance()
         if (ieee_is_finite(m) .and. ieee_is_finite(var)) then
            v=var+m*m
         else
            v=nan_dp()
         end if
      case default
         v=numeric_moment(self,order)
      end select
   end function distribution_raw_moment

   real(dp) function distribution_central_moment(self,order) result(v)
      class(distribution_t), intent(in) :: self
      integer, intent(in) :: order
      if (order<0) error stop 'central moment order must be nonnegative'
      select case(order)
      case(0)
         v=1.0_dp
      case(1)
         v=0.0_dp
      case(2)
         v=self%variance()
      case default
         v=numeric_central_moment(self,order)
      end select
   end function distribution_central_moment

   real(dp) function distribution_skewness(self) result(v)
      class(distribution_t), intent(in) :: self
      real(dp) :: s
      s=self%sd()
      if (.not.ieee_is_finite(s) .or. s<=0.0_dp) then
         v=nan_dp()
      else
         v=self%central_moment(3)/(s**3)
      end if
   end function distribution_skewness

   real(dp) function distribution_excess_kurtosis(self) result(v)
      class(distribution_t), intent(in) :: self
      real(dp) :: var
      var=self%variance()
      if (.not.ieee_is_finite(var) .or. var<=0.0_dp) then
         v=nan_dp()
      else
         v=self%central_moment(4)/(var*var)-3.0_dp
      end if
   end function distribution_excess_kurtosis

   subroutine density_vec(d,x,y,log_value)
      type(distribution_t), intent(in) :: d
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: y(size(x))
      logical, intent(in), optional :: log_value
      integer :: i
      do i=1,size(x); y(i)=d%density(x(i),log_value); end do
   end subroutine density_vec

   subroutine cdf_vec(d,x,y,lower_tail,log_p)
      type(distribution_t), intent(in) :: d
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: y(size(x))
      logical, intent(in), optional :: lower_tail,log_p
      integer :: i
      do i=1,size(x); y(i)=d%cdf(x(i),lower_tail,log_p); end do
   end subroutine cdf_vec

   subroutine sf_vec(d,x,y,log_p)
      type(distribution_t), intent(in) :: d
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: y(size(x))
      logical, intent(in), optional :: log_p
      integer :: i
      do i=1,size(x); y(i)=d%sf(x(i),log_p); end do
   end subroutine sf_vec

   subroutine logcdf_vec(d,x,y)
      type(distribution_t), intent(in) :: d
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: y(size(x))
      integer :: i
      do i=1,size(x); y(i)=d%logcdf(x(i)); end do
   end subroutine logcdf_vec

   subroutine logsf_vec(d,x,y)
      type(distribution_t), intent(in) :: d
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: y(size(x))
      integer :: i
      do i=1,size(x); y(i)=d%logsf(x(i)); end do
   end subroutine logsf_vec

   subroutine quantile_vec(d,p,x,lower_tail,log_p)
      type(distribution_t), intent(in) :: d
      real(dp), intent(in) :: p(:)
      real(dp), intent(out) :: x(size(p))
      logical, intent(in), optional :: lower_tail,log_p
      integer :: i
      do i=1,size(p); x(i)=d%quantile(p(i),lower_tail,log_p); end do
   end subroutine quantile_vec

   function add_dist_dist(a,b) result(c)
      type(distribution_t), intent(in) :: a,b
      type(distribution_t) :: c
      if (a%kind==K_NORM .and. b%kind==K_NORM) then
         c=normal_dist(a%par(1)+b%par(1),sqrt(a%par(2)**2+b%par(2)**2))
      else if (a%kind==K_POIS .and. b%kind==K_POIS) then
         c=poisson_dist(a%par(1)+b%par(1))
      else
         c=convolve(a,b)
      end if
   end function add_dist_dist

   function add_dist_real(a,b) result(c)
      type(distribution_t), intent(in) :: a
      real(dp), intent(in) :: b
      type(distribution_t) :: c
      c=affine_dist(a,1.0_dp,b)
   end function add_dist_real

   function add_real_dist(a,b) result(c)
      real(dp), intent(in) :: a
      type(distribution_t), intent(in) :: b
      type(distribution_t) :: c
      c=affine_dist(b,1.0_dp,a)
   end function add_real_dist

   function subtract_dist_dist(a,b) result(c)
      type(distribution_t), intent(in) :: a,b
      type(distribution_t) :: c, negb
      if (a%kind==K_NORM .and. b%kind==K_NORM) then
         c=normal_dist(a%par(1)-b%par(1),sqrt(a%par(2)**2+b%par(2)**2))
      else
         negb=affine_dist(b,-1.0_dp,0.0_dp)
         c=convolve(a,negb)
      end if
   end function subtract_dist_dist

   function subtract_dist_real(a,b) result(c)
      type(distribution_t), intent(in) :: a
      real(dp), intent(in) :: b
      type(distribution_t) :: c
      c=affine_dist(a,1.0_dp,-b)
   end function subtract_dist_real

   function subtract_real_dist(a,b) result(c)
      real(dp), intent(in) :: a
      type(distribution_t), intent(in) :: b
      type(distribution_t) :: c
      c=affine_dist(b,-1.0_dp,a)
   end function subtract_real_dist

   function negate_dist(a) result(c)
      type(distribution_t), intent(in) :: a
      type(distribution_t) :: c
      c=affine_dist(a,-1.0_dp,0.0_dp)
   end function negate_dist

   function multiply_dist_real(a,b) result(c)
      type(distribution_t), intent(in) :: a
      real(dp), intent(in) :: b
      type(distribution_t) :: c
      c=affine_dist(a,b,0.0_dp)
   end function multiply_dist_real

   function multiply_real_dist(a,b) result(c)
      real(dp), intent(in) :: a
      type(distribution_t), intent(in) :: b
      type(distribution_t) :: c
      c=affine_dist(b,a,0.0_dp)
   end function multiply_real_dist

   function divide_dist_real(a,b) result(c)
      type(distribution_t), intent(in) :: a
      real(dp), intent(in) :: b
      type(distribution_t) :: c
      if (b == 0.0_dp) error stop 'distribution division by zero'
      c=affine_dist(a,1.0_dp/b,0.0_dp)
   end function divide_dist_real

   real(dp) function integer_quantile(d,p) result(x)
      type(distribution_t), intent(in) :: d
      real(dp), intent(in) :: p
      integer :: lo,hi,mid
      select case(d%kind)
      case(K_BINOM); lo=0; hi=d%ipar(1)
      case(K_HYPER); lo=max(0,d%ipar(3)-d%ipar(2)); hi=min(d%ipar(3),d%ipar(1))
      case default
         lo=0; hi=1
         do
            if (raw_cdf(d,real(hi,dp)) >= p) exit
            if (hi >= huge(1)/4) exit
            hi = 2*hi + 1
         end do
      end select
      do while (lo<hi)
         mid=lo+(hi-lo)/2
         if (raw_cdf(d,real(mid,dp))>=p) then; hi=mid; else; lo=mid+1; end if
      end do
      x=real(lo,dp)
   end function integer_quantile

   real(dp) function integer_quantile_upper(d,p) result(x)
      type(distribution_t), intent(in) :: d
      real(dp), intent(in) :: p
      integer :: lo,hi,mid
      select case(d%kind)
      case(K_BINOM)
         lo=0; hi=d%ipar(1)
      case(K_HYPER)
         lo=max(0,d%ipar(3)-d%ipar(2)); hi=min(d%ipar(3),d%ipar(1))
      case default
         lo=0; hi=1
         do
            if (raw_sf(d,real(hi,dp)) <= p) exit
            if (hi >= huge(1)/4) exit
            hi=2*hi+1
         end do
      end select
      do while (lo<hi)
         mid=lo+(hi-lo)/2
         if (raw_sf(d,real(mid,dp))<=p) then; hi=mid; else; lo=mid+1; end if
      end do
      x=real(lo,dp)
   end function integer_quantile_upper

   real(dp) function support_quantile(d,p) result(x)
      type(distribution_t), intent(in) :: d
      real(dp), intent(in) :: p
      real(dp), allocatable :: ss(:),pp(:)
      real(dp) :: cum
      integer :: i
      call support_values(d,1.0e-12_dp,ss,pp)
      if (size(ss)==0) then; x=nan_dp(); return; end if
      cum=0.0_dp; x=ss(size(ss))
      do i=1,size(ss)
         cum=cum+pp(i)
         if (cum>=p) then; x=ss(i); return; end if
      end do
   end function support_quantile

   real(dp) function support_quantile_upper(d,p) result(x)
      type(distribution_t), intent(in) :: d
      real(dp), intent(in) :: p
      real(dp), allocatable :: ss(:),pp(:)
      real(dp) :: tail
      integer :: i
      call support_values(d,1.0e-12_dp,ss,pp)
      if (size(ss)==0) then; x=nan_dp(); return; end if
      tail=0.0_dp; x=ss(1)
      do i=size(ss),1,-1
         if (tail<=p) x=ss(i)
         tail=tail+pp(i)
         if (tail>p) return
      end do
   end function support_quantile_upper

   real(dp) function continuous_quantile_bisect(d,p) result(x)
      type(distribution_t), intent(in) :: d
      real(dp), intent(in) :: p
      real(dp) :: lo,hi,mid,f
      integer :: it
      call theoretical_bounds(d,lo,hi)
      if (.not.ieee_is_finite(lo)) then
         lo=-1.0_dp
         do while (raw_cdf(d,lo)>p); lo=2.0_dp*lo; if (lo < -huge(1.0_dp)/4.0_dp) exit; end do
      end if
      if (.not.ieee_is_finite(hi)) then
         hi=1.0_dp
         do while (raw_cdf(d,hi)<p); hi=2.0_dp*hi; if (hi > huge(1.0_dp)/4.0_dp) exit; end do
      end if
      do it=1,220
         mid=lo+0.5_dp*(hi-lo); f=raw_cdf(d,mid)
         if (f>=p) then; hi=mid; else; lo=mid; end if
         if (abs(hi-lo)<=64.0_dp*eps_dp*max(1.0_dp,abs(mid))) exit
      end do
      x=0.5_dp*(lo+hi)
   end function continuous_quantile_bisect

   real(dp) function continuous_quantile_upper_bisect(d,p) result(x)
      type(distribution_t), intent(in) :: d
      real(dp), intent(in) :: p
      real(dp) :: lo,hi,mid,f
      integer :: it
      call theoretical_bounds(d,lo,hi)
      if (.not.ieee_is_finite(lo)) then
         lo=-1.0_dp
         do while (raw_sf(d,lo)<p)
            lo=2.0_dp*lo
            if (lo < -huge(1.0_dp)/4.0_dp) exit
         end do
      end if
      if (.not.ieee_is_finite(hi)) then
         hi=1.0_dp
         do while (raw_sf(d,hi)>p)
            hi=2.0_dp*hi
            if (hi > huge(1.0_dp)/4.0_dp) exit
         end do
      end if
      do it=1,220
         mid=lo+0.5_dp*(hi-lo); f=raw_sf(d,mid)
         if (f<=p) then; hi=mid; else; lo=mid; end if
         if (abs(hi-lo)<=64.0_dp*eps_dp*max(1.0_dp,abs(mid))) exit
      end do
      x=0.5_dp*(lo+hi)
   end function continuous_quantile_upper_bisect

   recursive subroutine theoretical_bounds(d,lo,hi)
      type(distribution_t), intent(in) :: d
      real(dp), intent(out) :: lo,hi
      real(dp) :: a1,b1,a2,b2
      integer :: i
      select case(d%kind)
      case(K_DIRAC); lo=d%par(1); hi=d%par(1)
      case(K_BINOM); lo=0.0_dp; hi=real(d%ipar(1),dp)
      case(K_HYPER); lo=real(max(0,d%ipar(3)-d%ipar(2)),dp); hi=real(min(d%ipar(3),d%ipar(1)),dp)
      case(K_POIS,K_NBINOM,K_GEOM,K_EXP,K_GAMMA,K_CHISQ,K_F,K_WEIBULL,K_LNORM); lo=0.0_dp; hi=inf_dp()
      case(K_UNIF); lo=d%par(1); hi=d%par(2)
      case(K_BETA); lo=0.0_dp; hi=1.0_dp
      case(K_ARCSINE); lo=-1.0_dp; hi=1.0_dp
      case(K_DISCRETE,K_GRID); lo=minval(d%support); hi=maxval(d%support)
      case(K_NORM,K_CAUCHY,K_T,K_LAPLACE,K_LOGIS); lo=inf_dp(-1); hi=inf_dp()
      case(K_MIXTURE)
         lo=inf_dp(); hi=inf_dp(-1)
         do i=1,size(d%component); call theoretical_bounds(d%component(i),a1,b1); lo=min(lo,a1); hi=max(hi,b1); end do
      case(K_AFFINE)
         call theoretical_bounds(d%component(1),a1,b1)
         if (d%par(1)>0.0_dp) then; lo=d%par(1)*a1+d%par(2); hi=d%par(1)*b1+d%par(2)
         else; lo=d%par(1)*b1+d%par(2); hi=d%par(1)*a1+d%par(2); end if
      case(K_TRUNC); lo=d%par(1); hi=d%par(2)
      case(K_MIN)
         call theoretical_bounds(d%component(1),a1,b1); call theoretical_bounds(d%component(2),a2,b2)
         lo=min(a1,a2); hi=min(b1,b2)
      case(K_MAX)
         call theoretical_bounds(d%component(1),a1,b1); call theoretical_bounds(d%component(2),a2,b2)
         lo=max(a1,a2); hi=max(b1,b2)
      case(K_CONV)
         call theoretical_bounds(d%component(1),a1,b1); call theoretical_bounds(d%component(2),a2,b2)
         lo=a1+a2; hi=b1+b2
      case(K_EXPTRANS)
         call theoretical_bounds(d%component(1),a1,b1)
         if (ieee_is_finite(a1)) then; lo=exp(a1); else; lo=0.0_dp; end if
         if (ieee_is_finite(b1)) then; hi=exp(b1); else; hi=inf_dp(); end if
      case(K_LOGTRANS)
         call theoretical_bounds(d%component(1),a1,b1)
         if (.not.ieee_is_finite(a1) .or. a1<=0.0_dp) then; lo=inf_dp(-1); else; lo=log(a1); end if
         if (.not.ieee_is_finite(b1)) then; hi=inf_dp(); else; hi=log(b1); end if
      case(K_POWERTRANS)
         call theoretical_bounds(d%component(1),a1,b1)
         call power_bounds(a1,b1,d%par(1),d%ipar(1),lo,hi)
      case(K_ABSTRANS)
         call theoretical_bounds(d%component(1),a1,b1); lo=0.0_dp
         if (.not.ieee_is_finite(a1).or..not.ieee_is_finite(b1)) then; hi=inf_dp(); else; hi=max(abs(a1),abs(b1)); end if
      case(K_HUBER); lo=d%par(1); hi=d%par(2)
      case default; lo=inf_dp(-1); hi=inf_dp()
      end select
   end subroutine theoretical_bounds

   recursive subroutine support_values(d,eps,s,p)
      type(distribution_t), intent(in) :: d
      real(dp), intent(in) :: eps
      real(dp), allocatable, intent(out) :: s(:),p(:)
      real(dp), allocatable :: s1(:),p1(:),s2(:),p2(:),st(:),pt(:)
      integer :: i,j,k,n,lo,hi,pos
      real(dp) :: a,b
      select case(d%kind)
      case(K_DIRAC)
         allocate(s(1),p(1)); s=d%par(1); p=1.0_dp
      case(K_BINOM)
         lo=0; hi=d%ipar(1); allocate(s(hi-lo+1),p(hi-lo+1))
         do i=lo,hi; s(i-lo+1)=real(i,dp); p(i-lo+1)=raw_density(d,real(i,dp)); end do
      case(K_HYPER)
         lo=max(0,d%ipar(3)-d%ipar(2)); hi=min(d%ipar(3),d%ipar(1)); allocate(s(hi-lo+1),p(hi-lo+1))
         do i=lo,hi; s(i-lo+1)=real(i,dp); p(i-lo+1)=raw_density(d,real(i,dp)); end do
      case(K_POIS,K_NBINOM,K_GEOM)
         lo=0; hi=int(integer_quantile(d,1.0_dp-max(eps,1.0e-14_dp)))
         allocate(s(hi+1),p(hi+1)); do i=0,hi; s(i+1)=real(i,dp); p(i+1)=raw_density(d,real(i,dp)); end do
         p=p/sum(p)
      case(K_DISCRETE)
         allocate(s(size(d%support)),p(size(d%probs))); s=d%support; p=d%probs
      case(K_AFFINE)
         call support_values(d%component(1), eps, s1, p1)
         allocate(s(size(s1)), p(size(p1)))
         s = d%par(1)*s1 + d%par(2)
         p = p1
         call sort_collapse(s,p)
      case(K_TRUNC)
         call support_values(d%component(1),eps,s1,p1); n=count(s1>=d%par(1).and.s1<=d%par(2)); allocate(s(n),p(n)); k=0
         do i=1,size(s1); if (s1(i)>=d%par(1).and.s1(i)<=d%par(2)) then; k=k+1; s(k)=s1(i); p(k)=p1(i); end if; end do
         if (sum(p)>0.0_dp) p=p/sum(p)
      case(K_MIXTURE)
         n=0; do i=1,size(d%component); call support_values(d%component(i),eps,s1,p1); n=n+size(s1); end do
         allocate(st(n),pt(n)); pos=0
         do i = 1, size(d%component)
            call support_values(d%component(i), eps, s1, p1)
            st(pos+1:pos+size(s1)) = s1
            pt(pos+1:pos+size(p1)) = d%weights(i)*p1
            pos = pos + size(s1)
         end do
         call sort_collapse(st,pt); s=st; p=pt/sum(pt)
      case(K_CONV)
         call support_values(d%component(1),sqrt(eps),s1,p1); call support_values(d%component(2),sqrt(eps),s2,p2)
         allocate(st(size(s1)*size(s2)),pt(size(s1)*size(s2))); k=0
         do i=1,size(s1); do j=1,size(s2); k=k+1; st(k)=s1(i)+s2(j); pt(k)=p1(i)*p2(j); end do; end do
         call sort_collapse(st,pt); s=st; p=pt/sum(pt)
      case(K_MIN,K_MAX)
         call support_values(d%component(1),eps,s1,p1); call support_values(d%component(2),eps,s2,p2)
         allocate(st(size(s1)+size(s2)),pt(size(s1)+size(s2))); st=[s1,s2]; pt=1.0_dp; call sort_collapse(st,pt)
         allocate(s(size(st)),p(size(st))); s=st
         do i=1,size(s); p(i)=raw_density(d,s(i)); end do
         if (sum(p)>0.0_dp) p=p/sum(p)
      case(K_EXPTRANS)
         call support_values(d%component(1),eps,s1,p1); allocate(s(size(s1)),p(size(p1))); s=exp(s1); p=p1; call sort_collapse(s,p)
      case(K_LOGTRANS)
         call support_values(d%component(1),eps,s1,p1); allocate(s(size(s1)),p(size(p1))); s=log(s1); p=p1; call sort_collapse(s,p)
      case(K_POWERTRANS)
         call support_values(d%component(1),eps,s1,p1); allocate(s(size(s1)),p(size(p1)))
         if (d%ipar(1)==1) then
            s=s1**d%par(1)
         else
            s=s1**int(anint(d%par(1)))
         end if
         p=p1; call sort_collapse(s,p)
      case(K_ABSTRANS)
         call support_values(d%component(1),eps,s1,p1); allocate(s(size(s1)),p(size(p1))); s=abs(s1); p=p1; call sort_collapse(s,p)
      case(K_HUBER)
         call support_values(d%component(1),eps,s1,p1); allocate(s(size(s1)),p(size(p1)))
         s=min(d%par(2),max(d%par(1),s1)); p=p1; call sort_collapse(s,p)
      case default
         allocate(s(0),p(0))
      end select
   end subroutine support_values

   subroutine sort_collapse(x,p)
      real(dp), allocatable, intent(inout) :: x(:),p(:)
      real(dp), allocatable :: xx(:),pp(:)
      real(dp) :: tx,tp
      integer :: i,j,n,k
      n=size(x)
      do i=2,n
         tx=x(i); tp=p(i); j=i-1
         do while(j>=1)
            if (x(j)<=tx) exit
            x(j+1)=x(j); p(j+1)=p(j); j=j-1
         end do
         x(j+1)=tx; p(j+1)=tp
      end do
      allocate(xx(n),pp(n)); k=0
      do i=1,n
         if (k > 0) then
            if (nearly_equal(x(i),xx(k))) then
               pp(k)=pp(k)+p(i)
               cycle
            end if
         end if
         k=k+1; xx(k)=x(i); pp(k)=p(i)
      end do
      x=xx(:k); p=pp(:k)
      if (sum(p)>0.0_dp) p=p/sum(p)
   end subroutine sort_collapse


   recursive real(dp) function mixture_cdf_left(d,x) result(v)
      type(distribution_t), intent(in) :: d
      real(dp), intent(in) :: x
      integer :: i
      v=0.0_dp
      do i=1,size(d%component)
         v=v+d%weights(i)*cdf_left_raw(d%component(i),x)
      end do
      v=clamp01(v)
   end function mixture_cdf_left

   subroutine lattice_spacing(x,is_lattice,h)
      real(dp), intent(in) :: x(:)
      logical, intent(out) :: is_lattice
      real(dp), intent(out) :: h
      integer :: i
      if (size(x)<2) then
         is_lattice=.true.; h=0.0_dp; return
      end if
      h=x(2)-x(1)
      if (h<=0.0_dp) then; is_lattice=.false.; return; end if
      is_lattice=.true.
      do i=3,size(x)
         if (abs((x(i)-x(i-1))-h)>128.0_dp*eps_dp*max(1.0_dp,abs(h))) then
            is_lattice=.false.; return
         end if
      end do
   end subroutine lattice_spacing

   subroutine finite_probability_bounds(d,tail_prob,lo,hi)
      type(distribution_t), intent(in) :: d
      real(dp), intent(in) :: tail_prob
      real(dp), intent(out) :: lo,hi
      real(dp) :: tlo,thi,ep
      ep=max(1.0e-15_dp,min(1.0e-2_dp,tail_prob))
      call theoretical_bounds(d,tlo,thi)
      if (ieee_is_finite(tlo)) then
         lo=tlo
      else
         lo=raw_quantile(d,0.5_dp*ep)
      end if
      if (ieee_is_finite(thi)) then
         hi=thi
      else
         hi=raw_quantile_upper(d,0.5_dp*ep)
      end if
      if (.not.ieee_is_finite(lo) .or. .not.ieee_is_finite(hi) .or. hi<=lo) &
         error stop 'could not determine finite probability bounds'
   end subroutine finite_probability_bounds

   subroutine convolution_bounds(a,b,x,for_cdf,lo,hi)
      type(distribution_t), intent(in) :: a,b
      real(dp), intent(in) :: x
      logical, intent(in) :: for_cdf
      real(dp), intent(out) :: lo,hi
      real(dp) :: alo,ahi,blo,bhi
      call theoretical_bounds(a,alo,ahi)
      call theoretical_bounds(b,blo,bhi)
      if (ieee_is_finite(alo)) then
         lo=alo
      else
         lo=raw_quantile(a,1.0e-10_dp)
      end if
      if (ieee_is_finite(ahi)) then
         hi=ahi
      else
         hi=raw_quantile_upper(a,1.0e-10_dp)
      end if
      if (.not.for_cdf) then
         if (ieee_is_finite(bhi)) lo=max(lo,x-bhi)
      end if
      if (ieee_is_finite(blo)) hi=min(hi,x-blo)
   end subroutine convolution_bounds

   recursive real(dp) function convolution_density(a,b,x) result(v)
      type(distribution_t), intent(in) :: a,b
      real(dp), intent(in) :: x
      real(dp), allocatable :: s(:),p(:)
      real(dp) :: lo,hi
      integer :: i
      if (a%discrete) then
         call support_values(a,1.0e-11_dp,s,p); v=0.0_dp
         do i=1,size(s); v=v+p(i)*raw_density(b,x-s(i)); end do
         return
      else if (b%discrete) then
         call support_values(b,1.0e-11_dp,s,p); v=0.0_dp
         do i=1,size(s); v=v+p(i)*raw_density(a,x-s(i)); end do
         return
      end if
      call convolution_bounds(a,b,x,.false.,lo,hi)
      if (hi <= lo) then
         v=0.0_dp
      else
         v=adaptive_integral(integrand,lo,hi,1.0e-9_dp,18)
      end if
   contains
      function integrand(t) result(y)
         real(dp), intent(in) :: t
         real(dp) :: y
         y=raw_density(a,t)*raw_density(b,x-t)
      end function integrand
   end function convolution_density

   recursive real(dp) function convolution_cdf(a,b,x) result(v)
      type(distribution_t), intent(in) :: a,b
      real(dp), intent(in) :: x
      real(dp), allocatable :: s(:),p(:)
      real(dp) :: lo,hi
      integer :: i
      if (a%discrete) then
         call support_values(a,1.0e-11_dp,s,p); v=0.0_dp
         do i=1,size(s); v=v+p(i)*raw_cdf(b,x-s(i)); end do
         v=clamp01(v); return
      else if (b%discrete) then
         call support_values(b,1.0e-11_dp,s,p); v=0.0_dp
         do i=1,size(s); v=v+p(i)*raw_cdf(a,x-s(i)); end do
         v=clamp01(v); return
      end if
      call convolution_bounds(a,b,x,.true.,lo,hi)
      if (hi <= lo) then
         v=0.0_dp
      else
         v=adaptive_integral(integrand,lo,hi,1.0e-9_dp,18)
      end if
      v=clamp01(v)
   contains
      function integrand(t) result(y)
         real(dp), intent(in) :: t
         real(dp) :: y
         y=raw_density(a,t)*raw_cdf(b,x-t)
      end function integrand
   end function convolution_cdf

   recursive real(dp) function convolution_sf(a,b,x) result(v)
      type(distribution_t), intent(in) :: a,b
      real(dp), intent(in) :: x
      real(dp), allocatable :: ss(:),pp(:)
      real(dp) :: lo,hi,blo,bhi
      integer :: i
      if (a%discrete) then
         call support_values(a,1.0e-11_dp,ss,pp); v=0.0_dp
         do i=1,size(ss); v=v+pp(i)*raw_sf(b,x-ss(i)); end do
         v=clamp01(v); return
      else if (b%discrete) then
         call support_values(b,1.0e-11_dp,ss,pp); v=0.0_dp
         do i=1,size(ss); v=v+pp(i)*raw_sf(a,x-ss(i)); end do
         v=clamp01(v); return
      end if
      call finite_probability_bounds(a,1.0e-10_dp,lo,hi)
      call theoretical_bounds(b,blo,bhi)
      if (ieee_is_finite(bhi)) lo=max(lo,x-bhi)
      if (hi<=lo) then
         v=0.0_dp
      else
         v=adaptive_integral(integrand,lo,hi,1.0e-9_dp,18)
      end if
      v=clamp01(v)
   contains
      function integrand(t) result(y)
         real(dp), intent(in) :: t
         real(dp) :: y
         y=raw_density(a,t)*raw_sf(b,x-t)
      end function integrand
   end function convolution_sf

   recursive real(dp) function convolution_cdf_left(a,b,x) result(v)
      type(distribution_t), intent(in) :: a,b
      real(dp), intent(in) :: x
      real(dp), allocatable :: ss(:),pp(:)
      integer :: i
      if (a%discrete) then
         call support_values(a,1.0e-11_dp,ss,pp); v=0.0_dp
         do i=1,size(ss); v=v+pp(i)*cdf_left_raw(b,x-ss(i)); end do
      else if (b%discrete) then
         call support_values(b,1.0e-11_dp,ss,pp); v=0.0_dp
         do i=1,size(ss); v=v+pp(i)*cdf_left_raw(a,x-ss(i)); end do
      else
         v=convolution_cdf(a,b,x)
      end if
      v=clamp01(v)
   end function convolution_cdf_left

   recursive real(dp) function numeric_moment(d,order) result(v)
      type(distribution_t), intent(in) :: d
      integer, intent(in) :: order
      real(dp), allocatable :: s(:),p(:)
      real(dp) :: lo,hi
      integer :: i
      if (d%discrete) then
         call support_values(d,1.0e-11_dp,s,p); v=sum((s**order)*p); return
      end if
      lo=raw_quantile(d,1.0e-8_dp); hi=raw_quantile_upper(d,1.0e-8_dp)
      v=adaptive_integral(integrand,lo,hi,1.0e-7_dp,18)
   contains
      function integrand(x) result(y)
         real(dp), intent(in) :: x
         real(dp) :: y
         y=x**order*raw_density(d,x)
      end function integrand
   end function numeric_moment

   recursive real(dp) function numeric_central_moment(d,order) result(v)
      type(distribution_t), intent(in) :: d
      integer, intent(in) :: order
      real(dp), allocatable :: ss(:),pp(:)
      real(dp) :: lo,hi,mu
      mu=d%mean()
      if (.not.ieee_is_finite(mu)) then
         v=nan_dp(); return
      end if
      if (d%discrete) then
         call support_values(d,1.0e-11_dp,ss,pp)
         v=sum(((ss-mu)**order)*pp)
         return
      end if
      lo=raw_quantile(d,1.0e-8_dp); hi=raw_quantile_upper(d,1.0e-8_dp)
      v=adaptive_integral(integrand,lo,hi,1.0e-7_dp,18)
   contains
      function integrand(x) result(y)
         real(dp), intent(in) :: x
         real(dp) :: y
         y=(x-mu)**order*raw_density(d,x)
      end function integrand
   end function numeric_central_moment

   recursive real(dp) function adaptive_integral(f,a,b,tol,depth) result(v)
      procedure(scalar_fun) :: f
      real(dp), intent(in) :: a,b,tol
      integer, intent(in) :: depth
      real(dp) :: fa,fb,fc,s
      if (a==b) then; v=0.0_dp; return; end if
      fa=f(a); fb=f(b); fc=f(0.5_dp*(a+b)); s=(b-a)*(fa+4.0_dp*fc+fb)/6.0_dp
      v=adapt_step(f,a,b,fa,fb,fc,s,tol,depth)
   end function adaptive_integral

   recursive real(dp) function adapt_step(f,a,b,fa,fb,fc,s,tol,depth) result(v)
      procedure(scalar_fun) :: f
      real(dp), intent(in) :: a,b,fa,fb,fc,s,tol
      integer, intent(in) :: depth
      real(dp) :: c,d,e,fd,fe,sl,sr,s2
      c=0.5_dp*(a+b); d=0.5_dp*(a+c); e=0.5_dp*(c+b); fd=f(d); fe=f(e)
      sl=(c-a)*(fa+4.0_dp*fd+fc)/6.0_dp; sr=(b-c)*(fc+4.0_dp*fe+fb)/6.0_dp; s2=sl+sr
      if (depth<=0 .or. abs(s2-s)<=15.0_dp*tol) then
         v=s2+(s2-s)/15.0_dp
      else
         v=adapt_step(f,a,c,fa,fc,fd,sl,0.5_dp*tol,depth-1)+adapt_step(f,c,b,fc,fb,fe,sr,0.5_dp*tol,depth-1)
      end if
   end function adapt_step

   real(dp) function grid_density(d,x) result(v)
      type(distribution_t), intent(in) :: d
      real(dp), intent(in) :: x
      real(dp) :: h,left
      integer :: i,n
      n=size(d%support); h=d%support(2)-d%support(1); left=d%support(1)-0.5_dp*h
      i=int(floor((x-left)/h))+1
      if (i<1.or.i>n) then; v=0.0_dp; else; v=d%probs(i)/h; end if
   end function grid_density

   real(dp) function grid_cdf(d,x) result(v)
      type(distribution_t), intent(in) :: d
      real(dp), intent(in) :: x
      real(dp) :: h,left,frac
      integer :: i,n
      n=size(d%support); h=d%support(2)-d%support(1); left=d%support(1)-0.5_dp*h
      if (x<=left) then; v=0.0_dp; return; end if
      if (x>=d%support(n)+0.5_dp*h) then; v=1.0_dp; return; end if
      i=int(floor((x-left)/h))+1; frac=(x-(left+real(i-1,dp)*h))/h
      if (i>1) then; v=sum(d%probs(:i-1)); else; v=0.0_dp; end if
      v=v+frac*d%probs(i)
   end function grid_cdf

   real(dp) function grid_quantile(d,p) result(x)
      type(distribution_t), intent(in) :: d
      real(dp), intent(in) :: p
      real(dp) :: h,cum,prev,frac
      integer :: i
      h=d%support(2)-d%support(1); cum=0.0_dp
      do i=1,size(d%support)
         prev=cum; cum=cum+d%probs(i)
         if (p<=cum) then
            frac=(p-prev)/max(d%probs(i),tiny(1.0_dp))
            x=d%support(i)-0.5_dp*h+frac*h; return
         end if
      end do
      x=d%support(size(d%support))+0.5_dp*h
   end function grid_quantile

   real(dp) function grid_sf(d,x) result(v)
      type(distribution_t), intent(in) :: d
      real(dp), intent(in) :: x
      real(dp) :: h,left,right,frac
      integer :: i,n
      n=size(d%support); h=d%support(2)-d%support(1)
      left=d%support(1)-0.5_dp*h; right=d%support(n)+0.5_dp*h
      if (x<=left) then; v=1.0_dp; return; end if
      if (x>=right) then; v=0.0_dp; return; end if
      i=int(floor((x-left)/h))+1
      frac=(x-(left+real(i-1,dp)*h))/h
      if (i<n) then; v=sum(d%probs(i+1:)); else; v=0.0_dp; end if
      v=v+(1.0_dp-frac)*d%probs(i)
      v=clamp01(v)
   end function grid_sf

   real(dp) function grid_quantile_upper(d,p) result(x)
      type(distribution_t), intent(in) :: d
      real(dp), intent(in) :: p
      real(dp) :: h,tail,newtail,frac,left
      integer :: i,n
      n=size(d%support); h=d%support(2)-d%support(1)
      tail=0.0_dp
      do i=n,1,-1
         newtail=tail+d%probs(i)
         if (p<=newtail) then
            frac=(newtail-p)/max(d%probs(i),tiny(1.0_dp))
            left=d%support(i)-0.5_dp*h
            x=left+frac*h
            return
         end if
         tail=newtail
      end do
      x=d%support(1)-0.5_dp*h
   end function grid_quantile_upper

   elemental real(dp) function power_inverse_jacobian(root,y,power) result(jac)
      real(dp), intent(in) :: root,y,power
      if (y/=0.0_dp) then
         jac=abs(root/(power*y))
      else if (power<1.0_dp) then
         jac=0.0_dp
      else if (power==1.0_dp) then
         jac=1.0_dp
      else
         jac=huge(1.0_dp)
      end if
   end function power_inverse_jacobian

   subroutine power_bounds(a,b,power,mode,lo,hi)
      real(dp), intent(in) :: a,b,power
      integer, intent(in) :: mode
      real(dp), intent(out) :: lo,hi
      integer :: n
      select case(mode)
      case(1)
         if (power>0.0_dp) then
            if (a<=0.0_dp) then; lo=0.0_dp; else; lo=a**power; end if
            if (ieee_is_finite(b)) then; hi=b**power; else; hi=inf_dp(); end if
         else
            if (ieee_is_finite(b)) then; lo=b**power; else; lo=0.0_dp; end if
            if (a<=0.0_dp) then; hi=inf_dp(); else; hi=a**power; end if
         end if
      case(2)
         n=int(anint(power))
         if (ieee_is_finite(a)) then; lo=a**n; else; lo=inf_dp(-1); end if
         if (ieee_is_finite(b)) then; hi=b**n; else; hi=inf_dp(); end if
      case(3)
         n=int(anint(power))
         if (a<=0.0_dp .and. b>=0.0_dp) then
            lo=0.0_dp
         else if (ieee_is_finite(a) .and. ieee_is_finite(b)) then
            lo=min(abs(a),abs(b))**n
         else
            lo=0.0_dp
         end if
         if (.not.ieee_is_finite(a) .or. .not.ieee_is_finite(b)) then
            hi=inf_dp()
         else
            hi=max(abs(a),abs(b))**n
         end if
      case default
         lo=inf_dp(-1); hi=inf_dp()
      end select
   end subroutine power_bounds

   elemental real(dp) function log_nonnegative(x) result(v)
      real(dp), intent(in) :: x
      if (x>0.0_dp) then
         v=log(x)
      else if (x==0.0_dp) then
         v=inf_dp(-1)
      else
         v=nan_dp()
      end if
   end function log_nonnegative

   elemental real(dp) function logaddexp(a,b) result(v)
      real(dp), intent(in) :: a,b
      real(dp) :: m
      if (.not.ieee_is_finite(a)) then
         if (a>0.0_dp) then; v=a; else; v=b; end if
         return
      else if (.not.ieee_is_finite(b)) then
         if (b>0.0_dp) then; v=b; else; v=a; end if
         return
      end if
      m=max(a,b)
      v=m+log(exp(a-m)+exp(b-m))
   end function logaddexp

   elemental real(dp) function log_probability(p) result(v)
      real(dp), intent(in) :: p
      if (p>0.0_dp) then
         v=log(p)
      else if (p==0.0_dp) then
         v=inf_dp(-1)
      else
         v=nan_dp()
      end if
   end function log_probability

   elemental real(dp) function log1p_safe(x) result(y)
      real(dp), intent(in) :: x
      if (x<=-1.0_dp) then
         if (x==-1.0_dp) then; y=inf_dp(-1); else; y=nan_dp(); end if
      else if (abs(x)<1.0e-4_dp) then
         y=x*(1.0_dp+x*(-0.5_dp+x*(1.0_dp/3.0_dp+x*(-0.25_dp+x*0.2_dp))))
      else
         y=log(1.0_dp+x)
      end if
   end function log1p_safe

   elemental logical function nearly_equal(a,b) result(ok)
      real(dp), intent(in) :: a,b
      ok=abs(a-b)<=32.0_dp*eps_dp*max(1.0_dp,abs(a),abs(b))
   end function nearly_equal

   elemental real(dp) function expm1_safe(x) result(y)
      real(dp), intent(in) :: x
      if (abs(x)<1.0e-5_dp) then
         y=x*(1.0_dp+x*(0.5_dp+x*(1.0_dp/6.0_dp+x*(1.0_dp/24.0_dp+x/120.0_dp))))
      else
         y=exp(x)-1.0_dp
      end if
   end function expm1_safe

end module distr_core
