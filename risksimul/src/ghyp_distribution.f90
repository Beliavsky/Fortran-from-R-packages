! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Reused under GPL-2/GPL-3 from the ghyp-fortran numerical implementation.
! Derived from ghyp 1.6.5 by Marc Weibel, David Luethi, and Henriette-Elise Breymann.
module ghyp_distribution
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_finite
   use ghyp_kinds, only : dp, i8, pi, log_two_pi
   use ghyp_special, only : normal_cdf, normal_quantile, student_cdf, &
      log_bessel_k, gauss_legendre_rule
   use ghyp_rng, only : rng_state, seed_rng, normal_rng
   use ghyp_gig, only : rgig
   use ghyp_linalg, only : inverse_spd, logdet_spd, cholesky_lower, quadratic_form
   use ghyp_model, only : ghyp_model_type, moments_result, ghyp_moments, &
      model_gaussian, model_student, model_vg
   implicit none
   private

   type, public :: probability_result
      real(dp) :: value = 0.0_dp
      real(dp) :: standard_error = 0.0_dp
      integer :: evaluations = 0
      logical :: ok = .false.
      character(len=160) :: message = ''
   end type probability_result

   public :: log_dghyp, dghyp, pghyp, pghyp_rectangle, qghyp
   public :: rghyp, rghyp_one

   interface log_dghyp
      module procedure log_dghyp_uv
      module procedure log_dghyp_mv
   end interface
   interface dghyp
      module procedure dghyp_uv
      module procedure dghyp_mv
   end interface

contains

   function log_dghyp_uv(x, model) result(value)
      real(dp), intent(in) :: x
      type(ghyp_model_type), intent(in) :: model
      real(dp) :: value, xv(1)
      xv(1) = x
      value = log_dghyp_mv(xv,model)
   end function log_dghyp_uv

   function dghyp_uv(x, model) result(value)
      real(dp), intent(in) :: x
      type(ghyp_model_type), intent(in) :: model
      real(dp) :: value, lv
      lv = log_dghyp_uv(x,model)
      if (lv <= log(tiny(1.0_dp))) then
         value = 0.0_dp
      else if (lv >= log(huge(1.0_dp))) then
         value = huge(1.0_dp)
      else
         value = exp(lv)
      end if
   end function dghyp_uv

   function log_dghyp_mv(x, model) result(value)
      real(dp), intent(in) :: x(:)
      type(ghyp_model_type), intent(in) :: model
      real(dp) :: value
      real(dp), allocatable :: inv(:,:), y(:)
      real(dp) :: q, g, a, logdet, p, aa, bb, z, logc
      real(dp) :: alpha
      integer :: d
      logical :: ok

      if (.not. model%ok .or. size(x) /= model%dimension()) then
         value = -huge(1.0_dp)
         return
      end if
      d = model%dimension()
      call inverse_spd(model%scatter,inv,ok)
      if (.not. ok) then
         value = -huge(1.0_dp)
         return
      end if
      logdet = logdet_spd(model%scatter,ok)
      allocate(y(d))
      y = x-model%mu
      q = max(0.0_dp,quadratic_form(y,inv))
      g = max(0.0_dp,quadratic_form(model%gamma,inv))
      a = dot_product(y,matmul(inv,model%gamma))

      if (model%family == model_gaussian) then
         value = -0.5_dp*(real(d,dp)*log_two_pi+logdet+q)
         return
      end if

      p = model%lambda-0.5_dp*real(d,dp)
      aa = model%chi+q
      bb = model%psi+g

      if (model%psi <= tiny(1.0_dp) .and. g <= tiny(1.0_dp)) then
         alpha = -model%lambda
         value = log_gamma(alpha+0.5_dp*real(d,dp))-log_gamma(alpha)- &
            0.5_dp*(real(d,dp)*log(pi*model%chi)+logdet)- &
            (alpha+0.5_dp*real(d,dp))*log(1.0_dp+q/model%chi)
         return
      end if

      if (model%chi <= tiny(1.0_dp)) then
         if (q <= 32.0_dp*epsilon(1.0_dp)) then
            if (p > 0.0_dp) then
               logc = model%lambda*log(0.5_dp*model%psi)-log_gamma(model%lambda)- &
                  0.5_dp*(real(d,dp)*log_two_pi+logdet)+a
               value = logc+log_gamma(p)+p*log(2.0_dp/bb)
            else
               value = log(huge(1.0_dp))
            end if
         else
            z = sqrt(q*bb)
            value = model%lambda*log(0.5_dp*model%psi)-log_gamma(model%lambda)- &
               0.5_dp*(real(d,dp)*log_two_pi+logdet)+a+log(2.0_dp)+ &
               0.5_dp*p*log(q/bb)+log_bessel_k(p,z)
         end if
         return
      end if

      if (model%psi <= tiny(1.0_dp)) then
         z = sqrt(aa*max(g,tiny(1.0_dp)))
         alpha = -model%lambda
         value = alpha*log(0.5_dp*model%chi)-log_gamma(alpha)- &
            0.5_dp*(real(d,dp)*log_two_pi+logdet)+a+log(2.0_dp)+ &
            0.5_dp*p*log(aa/max(g,tiny(1.0_dp)))+log_bessel_k(p,z)
         return
      end if

      z = sqrt(aa*bb)
      value = 0.5_dp*model%lambda*log(model%psi/model%chi)- &
         log_bessel_k(model%lambda,sqrt(model%chi*model%psi))- &
         0.5_dp*(real(d,dp)*log_two_pi+logdet)+a+ &
         0.5_dp*p*log(aa/bb)+log_bessel_k(p,z)
   end function log_dghyp_mv

   function dghyp_mv(x, model) result(value)
      real(dp), intent(in) :: x(:)
      type(ghyp_model_type), intent(in) :: model
      real(dp) :: value, lv
      lv = log_dghyp_mv(x,model)
      if (lv <= log(tiny(1.0_dp))) then
         value = 0.0_dp
      else if (lv >= log(huge(1.0_dp))) then
         value = huge(1.0_dp)
      else
         value = exp(lv)
      end if
   end function dghyp_mv

   function pghyp(q, model, lower_tail) result(value)
      real(dp), intent(in) :: q
      type(ghyp_model_type), intent(in) :: model
      logical, intent(in), optional :: lower_tail
      real(dp) :: value, mu, sd, nu, scale, t, om, xv
      real(dp), allocatable :: nodes(:), weights(:)
      logical :: lower
      type(moments_result) :: mom
      integer :: i
      lower = .true.
      if (present(lower_tail)) lower = lower_tail
      if (model%dimension() /= 1 .or. .not. model%ok) then
         value = ieee_value(1.0_dp,ieee_quiet_nan)
         return
      end if
      if (model%family == model_gaussian) then
         value = normal_cdf((q-model%mu(1))/sqrt(model%scatter(1,1)))
      else if (model%family == model_student .and. model%is_symmetric()) then
         nu = -2.0_dp*model%lambda
         scale = sqrt(model%chi/nu)*sqrt(model%scatter(1,1))
         value = student_cdf((q-model%mu(1))/scale,nu)
      else
         mom = ghyp_moments(model)
         mu = model%mu(1)
         sd = sqrt(model%scatter(1,1))
         if (mom%ok .and. ieee_is_finite(mom%mean(1))) mu = mom%mean(1)
         if (mom%ok .and. ieee_is_finite(mom%covariance(1,1))) &
            sd = sqrt(max(mom%covariance(1,1),tiny(1.0_dp)))
         call gauss_legendre_rule(192,nodes,weights)
         value=0.0_dp
         do i=1,size(nodes)
            t=0.5_dp*(nodes(i)+1.0_dp)
            om=max(1.0_dp-t,1.0e-14_dp)
            xv=q-sd*t/om
            value=value+weights(i)*sd*dghyp_uv(xv,model)/(om*om)
         end do
         value=0.5_dp*value
         value = min(1.0_dp,max(0.0_dp,value))
      end if
      if (.not. lower) value = 1.0_dp-value
   end function pghyp

   function qghyp(p, model) result(value)
      real(dp), intent(in) :: p
      type(ghyp_model_type), intent(in) :: model
      real(dp) :: value, lo, hi, mid, mu, sd, cdf
      integer :: iter
      type(moments_result) :: mom
      if (p <= 0.0_dp) then
         value = -huge(1.0_dp)
         return
      else if (p >= 1.0_dp) then
         value = huge(1.0_dp)
         return
      end if
      if (model%family == model_gaussian) then
         value = model%mu(1)+sqrt(model%scatter(1,1))*normal_quantile(p)
         return
      end if
      mom = ghyp_moments(model)
      mu = model%mu(1)
      sd = sqrt(model%scatter(1,1))
      if (mom%ok .and. ieee_is_finite(mom%mean(1))) mu = mom%mean(1)
      if (mom%ok .and. ieee_is_finite(mom%covariance(1,1))) &
         sd = sqrt(max(mom%covariance(1,1),tiny(1.0_dp)))
      lo = mu-sd
      hi = mu+sd
      do while (pghyp(lo,model) > p)
         hi = lo
         lo = lo-2.0_dp*sd
         sd = 2.0_dp*sd
      end do
      sd = max(sqrt(model%scatter(1,1)),1.0e-3_dp)
      do while (pghyp(hi,model) < p)
         lo = hi
         hi = hi+2.0_dp*sd
         sd = 2.0_dp*sd
      end do
      do iter = 1, 100
         mid = 0.5_dp*(lo+hi)
         cdf = pghyp(mid,model)
         if (cdf < p) then
            lo = mid
         else
            hi = mid
         end if
         if (abs(hi-lo) <= 1.0e-9_dp*max(1.0_dp,abs(mid))) exit
      end do
      value = 0.5_dp*(lo+hi)
   end function qghyp

   function rghyp_one(model, rng) result(x)
      type(ghyp_model_type), intent(in) :: model
      type(rng_state), intent(inout) :: rng
      real(dp), allocatable :: x(:)
      real(dp), allocatable :: l(:,:), z(:)
      real(dp) :: w
      logical :: ok
      integer :: j, d
      d = model%dimension()
      allocate(x(d),z(d))
      call cholesky_lower(model%scatter,l,ok)
      if (.not. ok) then
         x = ieee_value(1.0_dp,ieee_quiet_nan)
         return
      end if
      do j = 1, d
         z(j) = normal_rng(rng)
      end do
      if (model%family == model_gaussian) then
         x = model%mu+matmul(l,z)
      else
         w = rgig(model%lambda,model%chi,model%psi,rng)
         x = model%mu+w*model%gamma+sqrt(w)*matmul(l,z)
      end if
   end function rghyp_one

   subroutine rghyp(n, model, sample, ok, seed)
      integer, intent(in) :: n
      type(ghyp_model_type), intent(in) :: model
      real(dp), allocatable, intent(out) :: sample(:,:)
      logical, intent(out) :: ok
      integer(i8), intent(in), optional :: seed
      type(rng_state) :: rng
      integer :: i
      call seed_rng(rng,123456789_i8)
      if (present(seed)) call seed_rng(rng,seed)
      ok = model%ok .and. n >= 0
      allocate(sample(max(n,0),model%dimension()))
      if (.not. ok) then
         sample = 0.0_dp
         return
      end if
      do i = 1, n
         sample(i,:) = rghyp_one(model,rng)
      end do
   end subroutine rghyp

   function pghyp_rectangle(upper, model, n_sim, seed) result(result)
      real(dp), intent(in) :: upper(:)
      type(ghyp_model_type), intent(in) :: model
      integer, intent(in), optional :: n_sim
      integer(i8), intent(in), optional :: seed
      type(probability_result) :: result
      real(dp), allocatable :: sample(:,:)
      integer :: i, n, count
      logical :: ok
      n = 20000
      if (present(n_sim)) n = max(100,n_sim)
      if (size(upper) /= model%dimension()) then
         result%message = 'upper bound has wrong dimension'
         return
      end if
      if (present(seed)) then
         call rghyp(n,model,sample,ok,seed)
      else
         call rghyp(n,model,sample,ok,246813579_i8)
      end if
      if (.not. ok) then
         result%message = 'simulation failed'
         return
      end if
      count = 0
      do i = 1, n
         if (all(sample(i,:) <= upper)) count = count+1
      end do
      result%value = real(count,dp)/real(n,dp)
      result%standard_error = sqrt(result%value*(1.0_dp-result%value)/real(n,dp))
      result%evaluations = n
      result%ok = .true.
   end function pghyp_rectangle

end module ghyp_distribution
