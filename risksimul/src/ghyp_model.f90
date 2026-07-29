! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Reused under GPL-2/GPL-3 from the ghyp-fortran numerical implementation.
! Derived from ghyp 1.6.5 by Marc Weibel, David Luethi, and Henriette-Elise Breymann.
module ghyp_model
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use ghyp_kinds, only : dp
   use ghyp_special, only : log_bessel_k
   use ghyp_gig, only : gig_valid, gig_mean, gig_variance
   use ghyp_linalg, only : cholesky_lower, symmetrize
   implicit none
   private

   integer, parameter, public :: model_ghyp = 1
   integer, parameter, public :: model_hyp = 2
   integer, parameter, public :: model_nig = 3
   integer, parameter, public :: model_student = 4
   integer, parameter, public :: model_vg = 5
   integer, parameter, public :: model_gaussian = 6

   type, public :: ghyp_model_type
      real(dp) :: lambda = 0.5_dp
      real(dp) :: chi = 0.5_dp
      real(dp) :: psi = 2.0_dp
      real(dp), allocatable :: mu(:)
      real(dp), allocatable :: scatter(:,:)
      real(dp), allocatable :: gamma(:)
      integer :: family = model_ghyp
      logical :: ok = .false.
      character(len=160) :: message = ''
   contains
      procedure :: dimension => ghyp_dimension
      procedure :: is_symmetric => ghyp_is_symmetric
   end type ghyp_model_type

   type, public :: moments_result
      real(dp), allocatable :: mean(:)
      real(dp), allocatable :: covariance(:,:)
      logical :: ok = .false.
      character(len=160) :: message = ''
   end type moments_result

   public :: make_ghyp, ghyp_uv, ghyp_mv, hyp_uv, nig_uv, student_t_uv, vg_uv
   public :: gaussian_uv, gaussian_mv, ghyp_ad
   public :: alpha_bar_to_chi_psi, ghyp_moments, transform_ghyp, ghyp_family_name

contains

   pure function ghyp_dimension(self) result(d)
      class(ghyp_model_type), intent(in) :: self
      integer :: d
      if (allocated(self%mu)) then
         d = size(self%mu)
      else
         d = 0
      end if
   end function ghyp_dimension

   pure function ghyp_is_symmetric(self) result(value)
      class(ghyp_model_type), intent(in) :: self
      logical :: value
      if (allocated(self%gamma)) then
         value = all(abs(self%gamma) <= 32.0_dp*epsilon(1.0_dp))
      else
         value = .true.
      end if
   end function ghyp_is_symmetric

   function infer_family(lambda, chi, psi, d) result(family)
      real(dp), intent(in) :: lambda, chi, psi
      integer, intent(in) :: d
      integer :: family
      real(dp), parameter :: tol = 1.0e-10_dp
      if (.not. ieee_is_finite(chi) .and. .not. ieee_is_finite(psi)) then
         family = model_gaussian
      else if (abs(psi) <= tol) then
         family = model_student
      else if (abs(chi) <= tol) then
         family = model_vg
      else if (abs(lambda+0.5_dp) <= tol) then
         family = model_nig
      else if (abs(lambda-0.5_dp*real(d+1,dp)) <= tol) then
         family = model_hyp
      else
         family = model_ghyp
      end if
   end function infer_family

   function make_ghyp(lambda, chi, psi, mu, scatter, gamma) result(model)
      real(dp), intent(in) :: lambda, chi, psi
      real(dp), intent(in) :: mu(:), scatter(:,:), gamma(:)
      type(ghyp_model_type) :: model
      real(dp), allocatable :: l(:,:)
      logical :: chol_ok, gaussian_case
      integer :: d

      d = size(mu)
      allocate(model%mu(d),model%gamma(d),model%scatter(d,d))
      model%mu = mu
      model%gamma = gamma
      model%scatter = scatter
      model%lambda = lambda
      model%chi = chi
      model%psi = psi
      model%family = infer_family(lambda,chi,psi,d)
      gaussian_case = model%family == model_gaussian

      if (d < 1 .or. size(gamma) /= d .or. size(scatter,1) /= d .or. &
          size(scatter,2) /= d) then
         model%message = 'incompatible model dimensions'
         return
      end if
      if (.not. all(ieee_is_finite(mu)) .or. .not. all(ieee_is_finite(gamma)) .or. &
          .not. all(ieee_is_finite(scatter))) then
         model%message = 'parameters must be finite'
         return
      end if
      call cholesky_lower(scatter,l,chol_ok)
      if (.not. chol_ok) then
         model%message = 'scatter matrix must be positive definite'
         return
      end if
      if (.not. gaussian_case) then
         if (.not. gig_valid(lambda,chi,psi)) then
            model%message = 'invalid GIG mixing parameters'
            return
         end if
      end if
      model%ok = .true.
      model%message = ''
   end function make_ghyp

   function ghyp_uv(lambda, chi, psi, mu, sigma, gamma) result(model)
      real(dp), intent(in) :: lambda, chi, psi, mu, sigma, gamma
      type(ghyp_model_type) :: model
      real(dp) :: muv(1), gamv(1), sc(1,1)
      muv = mu
      gamv = gamma
      sc(1,1) = sigma*sigma
      model = make_ghyp(lambda,chi,psi,muv,sc,gamv)
   end function ghyp_uv

   function ghyp_mv(lambda, chi, psi, mu, scatter, gamma) result(model)
      real(dp), intent(in) :: lambda, chi, psi
      real(dp), intent(in) :: mu(:), scatter(:,:), gamma(:)
      type(ghyp_model_type) :: model
      model = make_ghyp(lambda,chi,psi,mu,scatter,gamma)
   end function ghyp_mv

   function hyp_uv(chi, psi, mu, sigma, gamma) result(model)
      real(dp), intent(in), optional :: chi, psi, mu, sigma, gamma
      type(ghyp_model_type) :: model
      real(dp) :: c, p, m, s, g
      c = 0.5_dp; p = 2.0_dp; m = 0.0_dp; s = 1.0_dp; g = 0.0_dp
      if (present(chi)) c = chi
      if (present(psi)) p = psi
      if (present(mu)) m = mu
      if (present(sigma)) s = sigma
      if (present(gamma)) g = gamma
      model = ghyp_uv(1.0_dp,c,p,m,s,g)
   end function hyp_uv

   function nig_uv(chi, psi, mu, sigma, gamma) result(model)
      real(dp), intent(in), optional :: chi, psi, mu, sigma, gamma
      type(ghyp_model_type) :: model
      real(dp) :: c, p, m, s, g
      c = 2.0_dp; p = 2.0_dp; m = 0.0_dp; s = 1.0_dp; g = 0.0_dp
      if (present(chi)) c = chi
      if (present(psi)) p = psi
      if (present(mu)) m = mu
      if (present(sigma)) s = sigma
      if (present(gamma)) g = gamma
      model = ghyp_uv(-0.5_dp,c,p,m,s,g)
   end function nig_uv

   function student_t_uv(nu, chi, mu, sigma, gamma) result(model)
      real(dp), intent(in) :: nu
      real(dp), intent(in), optional :: chi, mu, sigma, gamma
      type(ghyp_model_type) :: model
      real(dp) :: c, m, s, g
      c = max(nu-2.0_dp,tiny(1.0_dp)); m = 0.0_dp; s = 1.0_dp; g = 0.0_dp
      if (present(chi)) c = chi
      if (present(mu)) m = mu
      if (present(sigma)) s = sigma
      if (present(gamma)) g = gamma
      model = ghyp_uv(-0.5_dp*nu,c,0.0_dp,m,s,g)
   end function student_t_uv

   function vg_uv(lambda, psi, mu, sigma, gamma) result(model)
      real(dp), intent(in) :: lambda
      real(dp), intent(in), optional :: psi, mu, sigma, gamma
      type(ghyp_model_type) :: model
      real(dp) :: p, m, s, g
      p = 2.0_dp*lambda; m = 0.0_dp; s = 1.0_dp; g = 0.0_dp
      if (present(psi)) p = psi
      if (present(mu)) m = mu
      if (present(sigma)) s = sigma
      if (present(gamma)) g = gamma
      model = ghyp_uv(lambda,0.0_dp,p,m,s,g)
   end function vg_uv

   function gaussian_uv(mu, sigma) result(model)
      real(dp), intent(in), optional :: mu, sigma
      type(ghyp_model_type) :: model
      real(dp) :: m, s
      m = 0.0_dp; s = 1.0_dp
      if (present(mu)) m = mu
      if (present(sigma)) s = sigma
      model = ghyp_uv(0.0_dp,huge(1.0_dp),huge(1.0_dp),m,s,0.0_dp)
      if (model%ok) model%family = model_gaussian
   end function gaussian_uv

   function gaussian_mv(mu, covariance) result(model)
      real(dp), intent(in) :: mu(:), covariance(:,:)
      type(ghyp_model_type) :: model
      real(dp) :: gamma(size(mu))
      gamma = 0.0_dp
      model = make_ghyp(0.0_dp,huge(1.0_dp),huge(1.0_dp),mu,covariance,gamma)
      if (model%ok) model%family = model_gaussian
   end function gaussian_mv

   subroutine alpha_bar_to_chi_psi(alpha_bar, lambda, chi, psi, ok)
      real(dp), intent(in) :: alpha_bar, lambda
      real(dp), intent(out) :: chi, psi
      logical, intent(out) :: ok
      real(dp) :: logratio
      ok = alpha_bar >= 0.0_dp
      if (.not. ok) then
         chi = 0.0_dp; psi = 0.0_dp
         return
      end if
      if (alpha_bar > sqrt(epsilon(1.0_dp))) then
         if (lambda >= 0.0_dp) then
            logratio = log_bessel_k(lambda+1.0_dp,alpha_bar)- &
               log_bessel_k(lambda,alpha_bar)
            psi = alpha_bar*exp(logratio)
            chi = alpha_bar*alpha_bar/psi
         else
            logratio = log_bessel_k(lambda,alpha_bar)- &
               log_bessel_k(lambda+1.0_dp,alpha_bar)
            chi = alpha_bar*exp(logratio)
            psi = alpha_bar*alpha_bar/chi
         end if
      else if (lambda > 0.0_dp) then
         chi = 0.0_dp
         psi = 2.0_dp*lambda
      else if (lambda < 0.0_dp) then
         psi = 0.0_dp
         chi = -2.0_dp*(lambda+1.0_dp)
         ok = chi > 0.0_dp
      else
         chi = 0.0_dp; psi = 0.0_dp; ok = .false.
      end if
   end subroutine alpha_bar_to_chi_psi

   function ghyp_ad(lambda, alpha, delta, beta, mu, delta_matrix) result(model)
      real(dp), intent(in) :: lambda, alpha, delta
      real(dp), intent(in) :: beta(:), mu(:), delta_matrix(:,:)
      type(ghyp_model_type) :: model
      real(dp) :: chi, psi
      real(dp), allocatable :: gamma(:)
      chi = delta*delta
      psi = alpha*alpha-dot_product(beta,matmul(delta_matrix,beta))
      if (abs(psi) < 1.0e-12_dp) psi = 0.0_dp
      allocate(gamma(size(beta)))
      gamma = matmul(delta_matrix,beta)
      model = make_ghyp(lambda,chi,psi,mu,delta_matrix,gamma)
   end function ghyp_ad

   function ghyp_moments(model) result(result)
      type(ghyp_model_type), intent(in) :: model
      type(moments_result) :: result
      real(dp) :: ew, vw
      integer :: d
      if (.not. model%ok) then
         result%message = 'invalid model'
         return
      end if
      d = model%dimension()
      allocate(result%mean(d),result%covariance(d,d))
      if (model%family == model_gaussian) then
         result%mean = model%mu
         result%covariance = model%scatter
      else
         ew = gig_mean(model%lambda,model%chi,model%psi)
         vw = gig_variance(model%lambda,model%chi,model%psi)
         result%mean = model%mu+ew*model%gamma
         result%covariance = ew*model%scatter+vw* &
            spread(model%gamma,2,d)*spread(model%gamma,1,d)
      end if
      call symmetrize(result%covariance)
      result%ok = .true.
   end function ghyp_moments

   function transform_ghyp(model, multiplier, summand) result(output)
      type(ghyp_model_type), intent(in) :: model
      real(dp), intent(in) :: multiplier(:,:)
      real(dp), intent(in), optional :: summand(:)
      type(ghyp_model_type) :: output
      real(dp), allocatable :: mu(:), gamma(:), scatter(:,:), shift(:)
      integer :: m
      m = size(multiplier,1)
      allocate(mu(m),gamma(m),scatter(m,m),shift(m))
      shift = 0.0_dp
      if (present(summand)) then
         if (size(summand) /= m) then
            output%message = 'summand has wrong dimension'
            return
         end if
         shift = summand
      end if
      if (size(multiplier,2) /= model%dimension()) then
         output%message = 'multiplier has wrong dimension'
         return
      end if
      mu = shift+matmul(multiplier,model%mu)
      gamma = matmul(multiplier,model%gamma)
      scatter = matmul(multiplier,matmul(model%scatter,transpose(multiplier)))
      output = make_ghyp(model%lambda,model%chi,model%psi,mu,scatter,gamma)
      if (model%family == model_gaussian .and. output%ok) output%family = model_gaussian
   end function transform_ghyp

   pure function ghyp_family_name(model) result(name)
      type(ghyp_model_type), intent(in) :: model
      character(len=32) :: name
      select case(model%family)
      case(model_hyp); name = 'hyperbolic'
      case(model_nig); name = 'normal inverse Gaussian'
      case(model_student); name = 'skewed Student t'
      case(model_vg); name = 'variance gamma'
      case(model_gaussian); name = 'Gaussian'
      case default; name = 'generalized hyperbolic'
      end select
   end function ghyp_family_name

end module ghyp_model
