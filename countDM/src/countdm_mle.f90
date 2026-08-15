module countdm_mle
   use countdm_kinds, only: dp
   use countdm_math, only: lambert_w0, logit, logistic, invert_matrix, numerical_hessian
   use countdm_optimizer, only: bfgs_minimize
   use countdm_distributions, only: dbell, dborel, dpoisson_count, dbellt, dzibellt, dzip, dzibell, dzoip, dzoibell
   implicit none
   private
   public :: mle_result_t, bell_closed_result_t
   public :: bell_mle, bell_mle_closed, mle_bell, mle_borel, mle_poisson, mle_bt
   public :: mle_zip, mle_zibell, mle_zibellt, mle_zoip, mle_zoibell

   integer, parameter :: MODEL_BELL = 1, MODEL_BOREL = 2, MODEL_POISSON = 3, MODEL_BT = 4
   integer, parameter :: MODEL_ZIP = 5, MODEL_ZIBELL = 6, MODEL_ZIBT = 7, MODEL_ZOIP = 8, MODEL_ZOIBELL = 9
   integer, save :: active_model = 0
   integer, allocatable, save :: active_x(:)

   type :: mle_result_t
      real(dp), allocatable :: estimate(:)
      real(dp), allocatable :: se(:)
      real(dp), allocatable :: covariance(:, :)
      real(dp) :: loglik = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      logical :: converged = .false.
      integer :: iterations = 0
   end type mle_result_t

   type :: bell_closed_result_t
      real(dp) :: theta = 0.0_dp
      real(dp) :: loglik = -huge(1.0_dp)
   end type bell_closed_result_t

contains

   function bell_mle_closed(x) result(res)
      integer, intent(in) :: x(:)
      type(bell_closed_result_t) :: res
      real(dp) :: meanx
      integer :: i
      if (size(x) == 0 .or. any(x < 0)) return
      meanx = real(sum(x), dp) / real(size(x), dp)
      res%theta = lambert_w0(meanx)
      if (res%theta <= 0.0_dp) return
      res%loglik = 0.0_dp
      do i = 1, size(x)
         res%loglik = res%loglik + dbell(x(i), res%theta, .true.)
      end do
   end function bell_mle_closed


   function bell_mle(x) result(res)
      integer, intent(in) :: x(:)
      type(bell_closed_result_t) :: res
      res = bell_mle_closed(x)
   end function bell_mle

   function mle_bell(x, theta0) result(res)
      integer, intent(in) :: x(:)
      real(dp), intent(in) :: theta0
      type(mle_result_t) :: res
      call fit_model(x, MODEL_BELL, [log(max(theta0, 1.0e-6_dp))], res)
   end function mle_bell

   function mle_borel(x, alpha0) result(res)
      integer, intent(in) :: x(:)
      real(dp), intent(in) :: alpha0
      type(mle_result_t) :: res
      call fit_model(x, MODEL_BOREL, [logit(alpha0)], res)
   end function mle_borel

   function mle_poisson(x, theta0) result(res)
      integer, intent(in) :: x(:)
      real(dp), intent(in) :: theta0
      type(mle_result_t) :: res
      call fit_model(x, MODEL_POISSON, [log(max(theta0, 1.0e-6_dp))], res)
   end function mle_poisson

   function mle_bt(x, lambda0, theta0) result(res)
      integer, intent(in) :: x(:)
      real(dp), intent(in) :: lambda0, theta0
      type(mle_result_t) :: res
      call fit_model(x, MODEL_BT, [log(max(lambda0, 1.0e-6_dp)), log(max(theta0, 1.0e-6_dp))], res)
   end function mle_bt

   function mle_zip(x, alpha0, theta0) result(res)
      integer, intent(in) :: x(:)
      real(dp), intent(in) :: alpha0, theta0
      type(mle_result_t) :: res
      call fit_model(x, MODEL_ZIP, [logit(alpha0), log(max(theta0, 1.0e-6_dp))], res)
   end function mle_zip

   function mle_zibell(x, alpha0, lambda0) result(res)
      integer, intent(in) :: x(:)
      real(dp), intent(in) :: alpha0, lambda0
      type(mle_result_t) :: res
      call fit_model(x, MODEL_ZIBELL, [logit(alpha0), log(max(lambda0, 1.0e-6_dp))], res)
   end function mle_zibell

   function mle_zibellt(x, lambda0, theta0, pi0) result(res)
      integer, intent(in) :: x(:)
      real(dp), intent(in) :: lambda0, theta0, pi0
      type(mle_result_t) :: res
      real(dp) :: z0(3)
      z0 = [log(max(lambda0, 1.0e-6_dp)), log(max(theta0, 1.0e-6_dp)), logit(pi0)]
      call fit_model(x, MODEL_ZIBT, z0, res)
   end function mle_zibellt

   function mle_zoip(x, alpha0, beta0, theta0) result(res)
      integer, intent(in) :: x(:)
      real(dp), intent(in) :: alpha0, beta0, theta0
      type(mle_result_t) :: res
      real(dp) :: base, z0(3)
      base = max(1.0_dp - alpha0 - beta0, 1.0e-8_dp)
      z0 = [log(max(alpha0, 1.0e-8_dp) / base), log(max(beta0, 1.0e-8_dp) / base), &
         log(max(theta0, 1.0e-6_dp))]
      call fit_model(x, MODEL_ZOIP, z0, res)
   end function mle_zoip

   function mle_zoibell(x, alpha0, beta0, theta0) result(res)
      integer, intent(in) :: x(:)
      real(dp), intent(in) :: alpha0, beta0, theta0
      type(mle_result_t) :: res
      real(dp) :: base, z0(3)
      base = max(1.0_dp - alpha0 - beta0, 1.0e-8_dp)
      z0 = [log(max(alpha0, 1.0e-8_dp) / base), log(max(beta0, 1.0e-8_dp) / base), &
         log(max(theta0, 1.0e-6_dp))]
      call fit_model(x, MODEL_ZOIBELL, z0, res)
   end function mle_zoibell

   subroutine fit_model(x, model, z0, res)
      integer, intent(in) :: x(:), model
      real(dp), intent(in) :: z0(:)
      type(mle_result_t), intent(out) :: res
      real(dp), allocatable :: z(:), p(:), jac(:, :), hess(:, :), covz(:, :), tmp(:, :)
      real(dp) :: nll
      logical :: ok
      integer :: k
      active_model = model
      if (allocated(active_x)) deallocate(active_x)
      allocate(active_x(size(x)))
      active_x = x
      z = z0
      call bfgs_minimize(nll_active, z, nll, res%converged, res%iterations)
      k = size(z)
      allocate(p(k), jac(k, k), hess(k, k))
      call decode(model, z, p, jac)
      allocate(res%estimate(k), res%se(k), res%covariance(k, k))
      res%estimate = p
      res%loglik = -nll
      res%aic = 2.0_dp * real(k, dp) + 2.0_dp * nll
      call numerical_hessian(nll_active, z, hess)
      call invert_matrix(hess, covz, ok)
      if (ok) then
         tmp = matmul(jac, covz)
         res%covariance = matmul(tmp, transpose(jac))
         res%se = sqrt(max(diagonal(res%covariance), 0.0_dp))
      else
         res%covariance = 0.0_dp
         res%se = huge(1.0_dp)
      end if
      active_model = 0
      deallocate(active_x)
   contains
      pure function diagonal(a) result(d)
         real(dp), intent(in) :: a(:, :)
         real(dp) :: d(min(size(a, 1), size(a, 2)))
         integer :: i
         do i = 1, size(d)
            d(i) = a(i, i)
         end do
      end function diagonal
   end subroutine fit_model

   function nll_active(z) result(v)
      real(dp), intent(in) :: z(:)
      real(dp) :: v
      real(dp), allocatable :: p(:), jac(:, :)
      real(dp) :: lp
      integer :: i
      allocate(p(size(z)), jac(size(z), size(z)))
      call decode(active_model, z, p, jac)
      v = 0.0_dp
      do i = 1, size(active_x)
         select case (active_model)
         case (MODEL_BELL)
            lp = dbell(active_x(i), p(1), .true.)
         case (MODEL_BOREL)
            lp = dborel(active_x(i), p(1), .true.)
         case (MODEL_POISSON)
            lp = dpoisson_count(active_x(i), p(1), .true.)
         case (MODEL_BT)
            lp = dbellt(active_x(i), p(1), p(2), .true.)
         case (MODEL_ZIP)
            lp = dzip(active_x(i), p(1), p(2), .true.)
         case (MODEL_ZIBELL)
            lp = dzibell(active_x(i), p(1), p(2), .true.)
         case (MODEL_ZIBT)
            lp = dzibellt(active_x(i), p(1), p(2), p(3), .true.)
         case (MODEL_ZOIP)
            lp = dzoip(active_x(i), p(1), p(2), p(3), .true.)
         case (MODEL_ZOIBELL)
            lp = dzoibell(active_x(i), p(1), p(2), p(3), .true.)
         case default
            lp = -huge(1.0_dp)
         end select
         if (lp <= -0.5_dp * huge(1.0_dp)) then
            v = huge(1.0_dp) / 10.0_dp
            return
         end if
         v = v - lp
         if (.not. (v < huge(1.0_dp) / 100.0_dp)) then
            v = huge(1.0_dp) / 10.0_dp
            return
         end if
      end do
   end function nll_active

   subroutine decode(model, z, p, jac)
      integer, intent(in) :: model
      real(dp), intent(in) :: z(:)
      real(dp), intent(out) :: p(:), jac(:, :)
      real(dp) :: m, ea, eb, e0, den
      jac = 0.0_dp
      select case (model)
      case (MODEL_BELL, MODEL_POISSON)
         p(1) = exp(z(1)); jac(1, 1) = p(1)
      case (MODEL_BOREL)
         p(1) = logistic(z(1)); jac(1, 1) = p(1) * (1.0_dp - p(1))
      case (MODEL_BT)
         p = exp(z); jac(1, 1) = p(1); jac(2, 2) = p(2)
      case (MODEL_ZIP, MODEL_ZIBELL)
         p(1) = logistic(z(1)); p(2) = exp(z(2))
         jac(1, 1) = p(1) * (1.0_dp - p(1)); jac(2, 2) = p(2)
      case (MODEL_ZIBT)
         p(1) = exp(z(1)); p(2) = exp(z(2)); p(3) = logistic(z(3))
         jac(1, 1) = p(1); jac(2, 2) = p(2); jac(3, 3) = p(3) * (1.0_dp - p(3))
      case (MODEL_ZOIP, MODEL_ZOIBELL)
         m = max(0.0_dp, max(z(1), z(2)))
         ea = exp(z(1) - m); eb = exp(z(2) - m); e0 = exp(-m)
         den = ea + eb + e0
         p(1) = ea / den; p(2) = eb / den; p(3) = exp(z(3))
         jac(1, 1) = p(1) * (1.0_dp - p(1)); jac(1, 2) = -p(1) * p(2)
         jac(2, 1) = -p(1) * p(2); jac(2, 2) = p(2) * (1.0_dp - p(2))
         jac(3, 3) = p(3)
      end select
   end subroutine decode

end module countdm_mle
