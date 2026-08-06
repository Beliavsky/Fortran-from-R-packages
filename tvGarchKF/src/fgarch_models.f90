module fgarch_models
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use fgarch_kinds, only : dp
   use fgarch_distributions, only : distribution_pdf, random_innovation
   use fgarch_types, only : garch_spec, model_garch, model_aparch, model_egarch
   implicit none
   private

   public :: garch_filter, garch_log_likelihood, simulate_garch
   public :: forecast_volatility, garch_kappa, true_persistence

contains

   subroutine garch_filter(y, spec, residuals, sigma, valid)
      real(dp), intent(in) :: y(:)
      type(garch_spec), intent(in) :: spec
      real(dp), intent(out) :: residuals(size(y)), sigma(size(y))
      logical, intent(out), optional :: valid
      real(dp), allocatable :: hpow(:), logvar(:), z(:)
      real(dp) :: v0, persistence, mean_abs
      integer :: n, i, j, start, par, qar, pg, qg
      logical :: ok

      n = size(y)
      par = size_or_zero(spec%ar)
      qar = size_or_zero(spec%ma)
      pg = size_or_zero(spec%alpha)
      qg = size_or_zero(spec%beta)
      start = max(1,max(max(par,qar),max(pg,qg))+1)
      ok = n > max(par,qar) .and. spec%omega > 0.0_dp .and. spec%delta > 0.0_dp
      if (size_or_zero(spec%gamma) /= pg) ok = .false.
      if (.not. ok) then
         residuals = 0.0_dp
         sigma = huge(1.0_dp)
         if (present(valid)) valid = .false.
         return
      end if

      v0 = sample_variance(y)
      v0 = max(v0,1.0e-12_dp)
      residuals = 0.0_dp

      do i = 1, n
         residuals(i) = y(i)-spec%mean
         do j = 1, min(par,i-1)
            residuals(i) = residuals(i)-spec%ar(j)*(y(i-j)-spec%mean)
         end do
         do j = 1, min(qar,i-1)
            residuals(i) = residuals(i)-spec%ma(j)*residuals(i-j)
         end do
      end do

      if (spec%model == model_egarch) then
         allocate(logvar(n),z(n))
         logvar = log(v0)
         sigma = sqrt(v0)
         z = residuals/sigma
         mean_abs = sqrt(2.0_dp/acos(-1.0_dp))
         do i = start, n
            logvar(i) = spec%omega
            do j = 1, pg
               logvar(i) = logvar(i)+spec%alpha(j)*(abs(z(i-j))-mean_abs)
               logvar(i) = logvar(i)+spec%gamma(j)*z(i-j)
            end do
            do j = 1, qg
               logvar(i) = logvar(i)+spec%beta(j)*logvar(i-j)
            end do
            if (.not. is_finite(logvar(i)) .or. abs(logvar(i)) > 100.0_dp) ok = .false.
            sigma(i) = exp(0.5_dp*max(-100.0_dp,min(100.0_dp,logvar(i))))
            z(i) = residuals(i)/max(sigma(i),1.0e-12_dp)
         end do
      else
         allocate(hpow(n))
         persistence = sum_if_allocated(spec%alpha)+sum_if_allocated(spec%beta)
         if (persistence < 0.999_dp .and. spec%omega > 0.0_dp) then
            hpow = max(spec%omega/(1.0_dp-persistence),v0**(0.5_dp*spec%delta))
         else
            hpow = v0**(0.5_dp*spec%delta)
         end if
         do i = start, n
            hpow(i) = spec%omega
            do j = 1, pg
               hpow(i) = hpow(i)+spec%alpha(j)* &
                  (abs(residuals(i-j))-spec%gamma(j)*residuals(i-j))**spec%delta
            end do
            do j = 1, qg
               hpow(i) = hpow(i)+spec%beta(j)*hpow(i-j)
            end do
            if (.not. is_finite(hpow(i)) .or. hpow(i) <= 0.0_dp) ok = .false.
            hpow(i) = max(hpow(i),1.0e-20_dp)
         end do
         sigma = hpow**(1.0_dp/spec%delta)
      end if

      if (any(.not. is_finite(sigma)) .or. any(sigma <= 0.0_dp)) ok = .false.
      if (present(valid)) valid = ok
   end subroutine garch_filter

   function garch_log_likelihood(y, spec, residuals, sigma) result(loglik)
      real(dp), intent(in) :: y(:)
      type(garch_spec), intent(in) :: spec
      real(dp), intent(out), optional :: residuals(size(y)), sigma(size(y))
      real(dp) :: loglik
      real(dp), allocatable :: eps(:), vol(:)
      real(dp) :: density
      integer :: i, start
      logical :: valid

      allocate(eps(size(y)),vol(size(y)))
      call garch_filter(y,spec,eps,vol,valid)
      if (.not. valid) then
         loglik = -huge(1.0_dp)
      else
         start = max(1,max(max(size_or_zero(spec%ar),size_or_zero(spec%ma)), &
                 max(size_or_zero(spec%alpha),size_or_zero(spec%beta)))+1)
         loglik = 0.0_dp
         do i = start, size(y)
            density = distribution_pdf(eps(i)/vol(i),spec%cond_dist,spec%shape,spec%skew)/vol(i)
            if (density <= tiny(1.0_dp) .or. .not. is_finite(density)) then
               loglik = -huge(1.0_dp)
               exit
            end if
            loglik = loglik+log(density)
         end do
      end if
      if (present(residuals)) residuals = eps
      if (present(sigma)) sigma = vol
   end function garch_log_likelihood

   subroutine simulate_garch(spec, n, y, sigma, residuals, burn_in)
      type(garch_spec), intent(in) :: spec
      integer, intent(in) :: n
      real(dp), intent(out) :: y(n), sigma(n), residuals(n)
      integer, intent(in), optional :: burn_in
      real(dp), allocatable :: yy(:), ss(:), eps(:), hpow(:), logvar(:), z(:)
      real(dp) :: persistence, h0, innovation, mean_abs
      integer :: burn, nt, i, j, par, qar, pg, qg, start

      burn = 500
      if (present(burn_in)) burn = max(0,burn_in)
      nt = n+burn
      par = size_or_zero(spec%ar)
      qar = size_or_zero(spec%ma)
      pg = size_or_zero(spec%alpha)
      qg = size_or_zero(spec%beta)
      start = max(2,max(max(par,qar),max(pg,qg))+1)
      allocate(yy(nt),ss(nt),eps(nt))
      yy = spec%mean
      eps = 0.0_dp

      if (spec%model == model_egarch) then
         allocate(logvar(nt),z(nt))
         if (sum_if_allocated(spec%beta) < 0.999_dp) then
            logvar = spec%omega/max(1.0e-6_dp,1.0_dp-sum_if_allocated(spec%beta))
         else
            logvar = log(max(spec%omega,1.0e-6_dp))
         end if
         ss = exp(0.5_dp*logvar)
         z = 0.0_dp
         mean_abs = sqrt(2.0_dp/acos(-1.0_dp))
         do i = 1, nt
            if (i >= start) then
               logvar(i) = spec%omega
               do j = 1, pg
                  logvar(i) = logvar(i)+spec%alpha(j)*(abs(z(i-j))-mean_abs)+spec%gamma(j)*z(i-j)
               end do
               do j = 1, qg
                  logvar(i) = logvar(i)+spec%beta(j)*logvar(i-j)
               end do
               ss(i) = exp(0.5_dp*max(-100.0_dp,min(100.0_dp,logvar(i))))
            end if
            innovation = random_innovation(spec%cond_dist,spec%shape,spec%skew)
            z(i) = innovation
            eps(i) = ss(i)*innovation
            yy(i) = spec%mean+eps(i)
            do j = 1, min(par,i-1)
               yy(i) = yy(i)+spec%ar(j)*(yy(i-j)-spec%mean)
            end do
            do j = 1, min(qar,i-1)
               yy(i) = yy(i)+spec%ma(j)*eps(i-j)
            end do
         end do
      else
         allocate(hpow(nt))
         persistence = sum_if_allocated(spec%alpha)+sum_if_allocated(spec%beta)
         if (persistence < 0.999_dp) then
            h0 = spec%omega/max(1.0e-6_dp,1.0_dp-persistence)
         else
            h0 = max(spec%omega,1.0e-6_dp)
         end if
         hpow = h0
         ss = hpow**(1.0_dp/spec%delta)
         do i = 1, nt
            if (i >= start) then
               hpow(i) = spec%omega
               do j = 1, pg
                  hpow(i) = hpow(i)+spec%alpha(j)* &
                     (abs(eps(i-j))-spec%gamma(j)*eps(i-j))**spec%delta
               end do
               do j = 1, qg
                  hpow(i) = hpow(i)+spec%beta(j)*hpow(i-j)
               end do
               ss(i) = max(hpow(i),1.0e-20_dp)**(1.0_dp/spec%delta)
            end if
            innovation = random_innovation(spec%cond_dist,spec%shape,spec%skew)
            eps(i) = ss(i)*innovation
            yy(i) = spec%mean+eps(i)
            do j = 1, min(par,i-1)
               yy(i) = yy(i)+spec%ar(j)*(yy(i-j)-spec%mean)
            end do
            do j = 1, min(qar,i-1)
               yy(i) = yy(i)+spec%ma(j)*eps(i-j)
            end do
         end do
      end if

      y = yy(burn+1:nt)
      sigma = ss(burn+1:nt)
      residuals = eps(burn+1:nt)
   end subroutine simulate_garch

   function garch_kappa(spec, leverage) result(value)
      type(garch_spec), intent(in) :: spec
      real(dp), intent(in) :: leverage
      real(dp) :: value
      integer, parameter :: ngrid = 8000
      real(dp), parameter :: limit = 50.0_dp
      real(dp) :: dx, x, fx
      integer :: i, weight

      dx = 2.0_dp*limit/real(ngrid,dp)
      value = 0.0_dp
      do i = 0, ngrid
         x = -limit+real(i,dp)*dx
         fx = (abs(x)-leverage*x)**spec%delta * &
              distribution_pdf(x,spec%cond_dist,spec%shape,spec%skew)
         if (i == 0 .or. i == ngrid) then
            weight = 1
         else if (mod(i,2) == 0) then
            weight = 2
         else
            weight = 4
         end if
         value = value+real(weight,dp)*fx
      end do
      value = value*dx/3.0_dp
   end function garch_kappa

   function true_persistence(spec) result(value)
      type(garch_spec), intent(in) :: spec
      real(dp) :: value
      integer :: i

      value = sum_if_allocated(spec%beta)
      do i = 1, size_or_zero(spec%alpha)
         value = value+spec%alpha(i)*garch_kappa(spec,spec%gamma(i))
      end do
   end function true_persistence

   subroutine forecast_volatility(spec, residuals, sigma, horizon, forecast)
      type(garch_spec), intent(in) :: spec
      real(dp), intent(in) :: residuals(:), sigma(:)
      integer, intent(in) :: horizon
      real(dp), intent(out) :: forecast(horizon)
      real(dp), allocatable :: eps_all(:), hp_all(:)
      real(dp) :: expected_shock
      integer :: n, pg, qg, h, j, idx

      n = size(residuals)
      pg = size_or_zero(spec%alpha)
      qg = size_or_zero(spec%beta)
      if (spec%model == model_egarch) then
         forecast = sigma(n)
         return
      end if
      allocate(eps_all(n+horizon),hp_all(n+horizon))
      eps_all(1:n) = residuals
      hp_all(1:n) = sigma**spec%delta
      do h = 1, horizon
         idx = n+h
         hp_all(idx) = spec%omega
         do j = 1, pg
            if (idx-j <= n) then
               expected_shock = (abs(eps_all(idx-j))-spec%gamma(j)*eps_all(idx-j))**spec%delta
            else
               expected_shock = garch_kappa(spec,spec%gamma(j))*hp_all(idx-j)
            end if
            hp_all(idx) = hp_all(idx)+spec%alpha(j)*expected_shock
         end do
         do j = 1, qg
            hp_all(idx) = hp_all(idx)+spec%beta(j)*hp_all(idx-j)
         end do
         forecast(h) = max(hp_all(idx),1.0e-20_dp)**(1.0_dp/spec%delta)
         eps_all(idx) = 0.0_dp
      end do
   end subroutine forecast_volatility

   pure integer function size_or_zero(x) result(n)
      real(dp), allocatable, intent(in) :: x(:)
      if (allocated(x)) then
         n = size(x)
      else
         n = 0
      end if
   end function size_or_zero

   pure function sum_if_allocated(x) result(value)
      real(dp), allocatable, intent(in) :: x(:)
      real(dp) :: value
      if (allocated(x)) then
         value = sum(x)
      else
         value = 0.0_dp
      end if
   end function sum_if_allocated

   pure function sample_variance(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value, meanx
      if (size(x) < 2) then
         value = 0.0_dp
      else
         meanx = sum(x)/real(size(x),dp)
         value = sum((x-meanx)**2)/real(size(x)-1,dp)
      end if
   end function sample_variance

   pure elemental logical function is_finite(x) result(value)
      real(dp), intent(in) :: x
      value = ieee_is_finite(x)
   end function is_finite

end module fgarch_models
