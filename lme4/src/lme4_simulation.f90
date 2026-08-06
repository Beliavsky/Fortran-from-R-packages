module lme4_simulation
   use lme4_kinds, only : dp, pi
   use lme4_types, only : random_term_t, covariance_block_t, family_binomial, &
      family_poisson, family_gamma, family_inverse_gaussian, &
      family_negative_binomial
   use lme4_covariance, only : build_random_design, total_random_effects
   use lme4_linalg, only : cholesky_lower
   implicit none
   private
   public :: simulate_lmm, simulate_glmm, set_random_seed

contains

   subroutine simulate_lmm(x, terms, beta, varcorr, sigma, y, random_effects, seed)
      real(dp), intent(in) :: x(:,:), beta(:)
      type(random_term_t), intent(in) :: terms(:)
      type(covariance_block_t), intent(in) :: varcorr(:)
      real(dp), intent(in) :: sigma
      real(dp), allocatable, intent(out) :: y(:), random_effects(:)
      integer, intent(in), optional :: seed
      real(dp), allocatable :: z(:,:), noise(:)
      integer, allocatable :: offsets(:)
      integer :: info

      if (present(seed)) call set_random_seed(seed)
      call build_random_design(terms,z,offsets)
      call draw_random_effects(terms,varcorr,random_effects,info)
      allocate(noise(size(x,1)))
      call random_normal_vector(noise)
      y = matmul(x,beta)+matmul(z,random_effects)+sigma*noise
   end subroutine simulate_lmm

   subroutine simulate_glmm(x, terms, beta, varcorr, family, y, random_effects, &
      offset, seed, dispersion)
      real(dp), intent(in) :: x(:,:), beta(:)
      type(random_term_t), intent(in) :: terms(:)
      type(covariance_block_t), intent(in) :: varcorr(:)
      integer, intent(in) :: family
      real(dp), allocatable, intent(out) :: y(:), random_effects(:)
      real(dp), intent(in), optional :: offset(:), dispersion
      integer, intent(in), optional :: seed
      real(dp), allocatable :: z(:,:), eta(:), uniforms(:)
      real(dp) :: disp, mu, lambda
      integer, allocatable :: offsets(:)
      integer :: info, i

      disp = 1.0_dp
      if (present(dispersion)) disp = dispersion
      if (present(seed)) call set_random_seed(seed)
      call build_random_design(terms,z,offsets)
      call draw_random_effects(terms,varcorr,random_effects,info)
      eta = matmul(x,beta)+matmul(z,random_effects)
      if (present(offset)) eta = eta+offset
      allocate(y(size(eta)),uniforms(size(eta)))
      select case (family)
      case (family_binomial)
         call random_number(uniforms)
         do i = 1, size(y)
            if (uniforms(i) < logistic(eta(i))) then
               y(i) = 1.0_dp
            else
               y(i) = 0.0_dp
            end if
         end do
      case (family_poisson)
         do i = 1, size(y)
            mu = exp(min(20.0_dp,max(-20.0_dp,eta(i))))
            y(i) = real(poisson_draw(mu),dp)
         end do
      case (family_gamma)
         do i = 1, size(y)
            mu = exp(min(20.0_dp,max(-20.0_dp,eta(i))))
            y(i) = gamma_draw(1.0_dp/disp,disp*mu)
         end do
      case (family_inverse_gaussian)
         do i = 1, size(y)
            mu = exp(min(20.0_dp,max(-20.0_dp,eta(i))))
            y(i) = inverse_gaussian_draw(mu,1.0_dp/disp)
         end do
      case (family_negative_binomial)
         do i = 1, size(y)
            mu = exp(min(20.0_dp,max(-20.0_dp,eta(i))))
            lambda = gamma_draw(disp,mu/disp)
            y(i) = real(poisson_draw(lambda),dp)
         end do
      case default
         y = 0.0_dp
      end select
   end subroutine simulate_glmm

   subroutine draw_random_effects(terms,varcorr,u,info)
      type(random_term_t), intent(in) :: terms(:)
      type(covariance_block_t), intent(in) :: varcorr(:)
      real(dp), allocatable, intent(out) :: u(:)
      integer, intent(out) :: info
      real(dp), allocatable :: l(:,:), z(:)
      integer :: total, k, lev, q, pos

      total = total_random_effects(terms)
      allocate(u(total))
      u = 0.0_dp
      pos = 0
      info = 0
      do k = 1, size(terms)
         q = terms(k)%n_coefficients()
         call cholesky_lower(varcorr(k)%covariance,l,info)
         if (info /= 0) return
         allocate(z(q))
         do lev = 1, terms(k)%n_levels
            call random_normal_vector(z)
            u(pos+1:pos+q) = matmul(l,z)
            pos = pos+q
         end do
         deallocate(l,z)
      end do
   end subroutine draw_random_effects

   subroutine random_normal_vector(x)
      real(dp), intent(out) :: x(:)
      real(dp) :: u1, u2, radius
      integer :: i
      i = 1
      do while (i <= size(x))
         call random_number(u1)
         call random_number(u2)
         u1 = max(u1,tiny(1.0_dp))
         radius = sqrt(-2.0_dp*log(u1))
         x(i) = radius*cos(2.0_dp*pi*u2)
         if (i+1 <= size(x)) x(i+1) = radius*sin(2.0_dp*pi*u2)
         i = i+2
      end do
   end subroutine random_normal_vector

   integer function poisson_draw(lambda) result(value)
      real(dp), intent(in) :: lambda
      real(dp) :: p, product, u, z
      integer :: k
      if (lambda < 30.0_dp) then
         p = exp(-lambda)
         product = 1.0_dp
         k = 0
         do
            k = k+1
            call random_number(u)
            product = product*u
            if (product <= p) exit
         end do
         value = k-1
      else
         call random_normal_vector_scalar(z)
         value = max(0,nint(lambda+sqrt(lambda)*z))
      end if
   end function poisson_draw

   recursive real(dp) function gamma_draw(shape, scale) result(value)
      real(dp), intent(in) :: shape, scale
      real(dp) :: d, c, x, v, u
      if (shape <= 0.0_dp .or. scale <= 0.0_dp) then
         value = 0.0_dp
         return
      end if
      if (shape < 1.0_dp) then
         call random_number(u)
         value = gamma_draw(shape+1.0_dp,scale)*max(u,tiny(1.0_dp))**(1.0_dp/shape)
         return
      end if
      d = shape-1.0_dp/3.0_dp
      c = 1.0_dp/sqrt(9.0_dp*d)
      do
         call random_normal_vector_scalar(x)
         v = 1.0_dp+c*x
         if (v <= 0.0_dp) cycle
         v = v**3
         call random_number(u)
         if (u < 1.0_dp-0.0331_dp*x**4) exit
         if (log(max(u,tiny(1.0_dp))) < 0.5_dp*x*x+d*(1.0_dp-v+log(v))) exit
      end do
      value = scale*d*v
   end function gamma_draw

   real(dp) function inverse_gaussian_draw(mu, lambda) result(value)
      real(dp), intent(in) :: mu, lambda
      real(dp) :: z, v, x, u
      call random_normal_vector_scalar(z)
      v = z*z
      x = mu+mu*mu*v/(2.0_dp*lambda) - &
         mu*sqrt(4.0_dp*mu*lambda*v+mu*mu*v*v)/(2.0_dp*lambda)
      call random_number(u)
      if (u <= mu/(mu+x)) then
         value = x
      else
         value = mu*mu/x
      end if
   end function inverse_gaussian_draw

   subroutine random_normal_vector_scalar(x)
      real(dp), intent(out) :: x
      real(dp) :: u1, u2
      call random_number(u1)
      call random_number(u2)
      x = sqrt(-2.0_dp*log(max(u1,tiny(1.0_dp))))*cos(2.0_dp*pi*u2)
   end subroutine random_normal_vector_scalar

   pure real(dp) function logistic(x) result(value)
      real(dp), intent(in) :: x
      if (x >= 0.0_dp) then
         value = 1.0_dp/(1.0_dp+exp(-min(700.0_dp,x)))
      else
         value = exp(max(-700.0_dp,x))/(1.0_dp+exp(max(-700.0_dp,x)))
      end if
   end function logistic

   subroutine set_random_seed(seed)
      integer, intent(in) :: seed
      integer :: n, i
      integer, allocatable :: values(:)
      call random_seed(size=n)
      allocate(values(n))
      do i = 1, n
         values(i) = modulo(seed+104729*i,huge(1)-1)
         if (values(i) <= 0) values(i) = i
      end do
      call random_seed(put=values)
   end subroutine set_random_seed

end module lme4_simulation
