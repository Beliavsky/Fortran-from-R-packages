! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from ghyp 1.6.5 by Marc Weibel, David Luethi, and Henriette-Elise Breymann.
module ghyp_risk
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_finite
   use ghyp_kinds, only : dp
   use ghyp_special, only : gauss_legendre_rule
   use ghyp_gig, only : gig_raw_moment
   use ghyp_model, only : ghyp_model_type, ghyp_moments, moments_result, &
      model_gaussian, transform_ghyp
   use ghyp_distribution, only : dghyp, qghyp
   implicit none
   private

   type, public :: attribution_result
      real(dp) :: portfolio_es = 0.0_dp
      real(dp), allocatable :: contribution(:)
      real(dp), allocatable :: sensitivity(:)
      logical :: ok = .false.
      character(len=160) :: message = ''
   end type attribution_result

   public :: ghyp_moment, ghyp_skewness, ghyp_kurtosis
   public :: esghyp, ghyp_omega, esghyp_attribution

contains

   pure function binomial(n, k) result(value)
      integer, intent(in) :: n, k
      real(dp) :: value
      integer :: i, kk
      if (k < 0 .or. k > n) then
         value = 0.0_dp
         return
      end if
      kk = min(k,n-k)
      value = 1.0_dp
      do i = 1, kk
         value = value*real(n-kk+i,dp)/real(i,dp)
      end do
   end function binomial

   pure function odd_double_factorial(k) result(value)
      integer, intent(in) :: k
      real(dp) :: value
      integer :: i
      value = 1.0_dp
      do i = 1, k
         value = value*real(2*i-1,dp)
      end do
   end function odd_double_factorial

   function integer_moment(model, n, offset) result(value)
      type(ghyp_model_type), intent(in) :: model
      integer, intent(in) :: n
      real(dp), intent(in) :: offset
      real(dp) :: value, c, s, g, term, ew
      integer :: k, j, m
      if (n == 0) then
         value = 1.0_dp
         return
      end if
      c = model%mu(1)-offset
      s = sqrt(model%scatter(1,1))
      g = model%gamma(1)
      if (model%family == model_gaussian) then
         value = 0.0_dp
         do k = 0, n/2
            m = n-2*k
            term = binomial(n,2*k)*odd_double_factorial(k)*s**(2*k)*c**m
            value = value+term
         end do
         return
      end if
      value = 0.0_dp
      do k = 0, n/2
         m = n-2*k
         do j = 0, m
            ew = gig_raw_moment(real(k+j,dp),model%lambda,model%chi,model%psi)
            if (.not. ieee_is_finite(ew)) then
               value = ew
               return
            end if
            term = binomial(n,2*k)*odd_double_factorial(k)*s**(2*k)* &
               binomial(m,j)*c**(m-j)*g**j*ew
            value = value+term
         end do
      end do
   end function integer_moment

   function ghyp_moment(model, order, absolute, central) result(value)
      type(ghyp_model_type), intent(in) :: model
      real(dp), intent(in) :: order
      logical, intent(in), optional :: absolute, central
      real(dp) :: value, offset, t, u, c, x, base, f
      real(dp), allocatable :: nodes(:), weights(:)
      logical :: use_abs, use_central
      type(moments_result) :: mom
      integer :: n, i
      use_abs = .false.; use_central = .true.
      if (present(absolute)) use_abs = absolute
      if (present(central)) use_central = central
      if (.not. model%ok .or. model%dimension() /= 1 .or. order < 0.0_dp) then
         value = ieee_value(1.0_dp,ieee_quiet_nan)
         return
      end if
      offset = 0.0_dp
      if (use_central) then
         mom = ghyp_moments(model)
         if (.not. mom%ok) then
            value = ieee_value(1.0_dp,ieee_quiet_nan)
            return
         end if
         offset = mom%mean(1)
      end if
      n = nint(order)
      if (.not. use_abs .and. abs(order-real(n,dp)) <= 16.0_dp*epsilon(1.0_dp)) then
         value = integer_moment(model,n,offset)
      else
         call gauss_legendre_rule(256,nodes,weights)
         value=0.0_dp
         do i=1,size(nodes)
            t=0.5_dp*(nodes(i)+1.0_dp)
            u=acos(-1.0_dp)*(t-0.5_dp)
            c=cos(u);x=tan(u);base=x-offset
            if(use_abs)base=abs(base)
            f=0.0_dp
            if(.not.(base<0.0_dp.and.abs(order-real(n,dp))>16.0_dp*epsilon(1.0_dp))) then
               f=acos(-1.0_dp)*base**order*dghyp(x,model)/(c*c)
               if(.not.ieee_is_finite(f))f=0.0_dp
            end if
            value=value+weights(i)*f
         end do
         value=0.5_dp*value
      end if
   end function ghyp_moment

   function ghyp_skewness(model) result(value)
      type(ghyp_model_type), intent(in) :: model
      real(dp) :: value, m2, m3
      m2 = ghyp_moment(model,2.0_dp,central=.true.)
      m3 = ghyp_moment(model,3.0_dp,central=.true.)
      if (m2 <= 0.0_dp .or. .not. ieee_is_finite(m2)) then
         value = ieee_value(1.0_dp,ieee_quiet_nan)
      else
         value = m3/m2**1.5_dp
      end if
   end function ghyp_skewness

   function ghyp_kurtosis(model) result(value)
      type(ghyp_model_type), intent(in) :: model
      real(dp) :: value, m2, m4
      m2 = ghyp_moment(model,2.0_dp,central=.true.)
      m4 = ghyp_moment(model,4.0_dp,central=.true.)
      if (m2 <= 0.0_dp .or. .not. ieee_is_finite(m2)) then
         value = ieee_value(1.0_dp,ieee_quiet_nan)
      else
         value = m4/(m2*m2)
      end if
   end function ghyp_kurtosis

   function esghyp(alpha, model, loss) result(value)
      real(dp), intent(in) :: alpha
      type(ghyp_model_type), intent(in) :: model
      logical, intent(in), optional :: loss
      real(dp) :: value, q, denom, t, om, x, scale
      real(dp), allocatable :: nodes(:), weights(:)
      logical :: upper
      integer :: i
      upper = .false.
      if (present(loss)) upper = loss
      if (alpha <= 0.0_dp .or. alpha >= 1.0_dp .or. model%dimension() /= 1) then
         value = ieee_value(1.0_dp,ieee_quiet_nan)
         return
      end if
      q = qghyp(alpha,model)
      scale=sqrt(model%scatter(1,1))
      call gauss_legendre_rule(256,nodes,weights)
      value=0.0_dp
      do i=1,size(nodes)
         t=0.5_dp*(nodes(i)+1.0_dp);om=max(1.0_dp-t,1.0e-14_dp)
         if(upper)then;x=q+scale*t/om;else;x=q-scale*t/om;end if
         value=value+weights(i)*scale*x*dghyp(x,model)/(om*om)
      end do
      value=0.5_dp*value
      if (upper) then
         denom = 1.0_dp-alpha
      else
         denom = alpha
      end if
      value = value/denom
   end function esghyp

   function ghyp_omega(level, model) result(value)
      real(dp), intent(in) :: level
      type(ghyp_model_type), intent(in) :: model
      real(dp) :: value, numerator, denominator, t, om, xp, xm, scale
      real(dp), allocatable :: nodes(:), weights(:)
      integer :: i
      if (model%dimension() /= 1) then
         value = ieee_value(1.0_dp,ieee_quiet_nan)
         return
      end if
      scale=sqrt(model%scatter(1,1));numerator=0.0_dp;denominator=0.0_dp
      call gauss_legendre_rule(256,nodes,weights)
      do i=1,size(nodes)
         t=0.5_dp*(nodes(i)+1.0_dp);om=max(1.0_dp-t,1.0e-14_dp)
         xp=level+scale*t/om;xm=level-scale*t/om
         numerator=numerator+weights(i)*scale*(xp-level)*dghyp(xp,model)/(om*om)
         denominator=denominator+weights(i)*scale*(level-xm)*dghyp(xm,model)/(om*om)
      end do
      numerator=0.5_dp*numerator;denominator=0.5_dp*denominator
      if (denominator <= 0.0_dp) then
         value = huge(1.0_dp)
      else
         value = numerator/denominator
      end if
   end function ghyp_omega

   function esghyp_attribution(alpha, model, weights, loss) result(result)
      real(dp), intent(in) :: alpha
      type(ghyp_model_type), intent(in) :: model
      real(dp), intent(in) :: weights(:)
      logical, intent(in), optional :: loss
      type(attribution_result) :: result
      real(dp), allocatable :: a(:,:), wp(:), wm(:)
      type(ghyp_model_type) :: portfolio, plus_model, minus_model
      real(dp) :: h, ep, em
      logical :: upper
      integer :: i, d
      upper = .true.
      if (present(loss)) upper = loss
      d = model%dimension()
      if (size(weights) /= d) then
         result%message = 'weights have wrong dimension'
         return
      end if
      allocate(a(1,d),wp(d),wm(d),result%contribution(d),result%sensitivity(d))
      a(1,:) = weights
      portfolio = transform_ghyp(model,a)
      if (.not. portfolio%ok) then
         result%message = 'portfolio transformation failed'
         return
      end if
      result%portfolio_es = esghyp(alpha,portfolio,upper)
      do i = 1, d
         h = 1.0e-5_dp*max(1.0_dp,abs(weights(i)))
         wp = weights; wm = weights
         wp(i) = wp(i)+h; wm(i) = wm(i)-h
         a(1,:) = wp; plus_model = transform_ghyp(model,a)
         a(1,:) = wm; minus_model = transform_ghyp(model,a)
         ep = esghyp(alpha,plus_model,upper)
         em = esghyp(alpha,minus_model,upper)
         result%sensitivity(i) = (ep-em)/(2.0_dp*h)
         result%contribution(i) = weights(i)*result%sensitivity(i)
      end do
      result%ok = .true.
   end function esghyp_attribution

end module ghyp_risk
