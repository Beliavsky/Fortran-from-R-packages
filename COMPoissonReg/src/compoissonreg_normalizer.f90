module compoissonreg_normalizer
   use compoissonreg_kinds, only : dp, pi_dp
   use compoissonreg_types, only : cmp_control_t
   use compoissonreg_numerics, only : logadd
   implicit none
   private
   public :: cmp_logz_trunc, cmp_logz_approx, cmp_logz_hybrid
   public :: cmp_z_trunc, cmp_z_approx, cmp_z_hybrid, cmp_y_trunc
   public :: z_prodj, z_prodj2, z_prodjlogj, z_prodlogj, z_prodlogj2

contains

   subroutine cmp_logz_trunc(lambda, nu, tol, ymax, logz, last_y, converged)
      real(dp), intent(in) :: lambda, nu, tol
      integer, intent(in) :: ymax
      real(dp), intent(out) :: logz
      integer, intent(out) :: last_y
      logical, intent(out) :: converged
      real(dp) :: log_tol, diff, lp, log_ratio, log_delta, den, llambda
      integer :: y

      converged = .true.
      if (lambda < 0.0_dp .or. nu < 0.0_dp .or. tol <= 0.0_dp) then
         logz = huge(1.0_dp); last_y = 0; converged = .false.; return
      end if
      if (lambda == 0.0_dp) then
         logz = 0.0_dp; last_y = 0; return
      end if
      if (nu == 0.0_dp) then
         if (lambda < 1.0_dp) then
            logz = -log(1.0_dp-lambda); last_y = min(ymax, 1000); return
         else
            logz = huge(1.0_dp); last_y = ymax; converged = .false.; return
         end if
      end if

      log_tol = log(tol)
      diff = huge(1.0_dp)
      logz = 0.0_dp
      llambda = log(lambda)
      lp = 0.0_dp
      last_y = 0
      do y=1,ymax
         lp = lp + llambda - nu*log(real(y,dp))
         logz = logadd(logz, lp)
         log_ratio = llambda + nu - nu*log(real(y+1,dp))
         if (log_ratio < 0.0_dp) then
            den = 1.0_dp - lambda*exp(nu)/real(y+1,dp)**nu
            if (den > 0.0_dp) then
               log_delta = -0.5_dp*nu*log(2.0_dp*pi_dp) &
                  - nu*(real(y,dp)+1.5_dp)*log(real(y+1,dp)) &
                  + real(y+1,dp)*(nu+llambda) - log(den)
               diff = log_delta - logz
            end if
         end if
         last_y = y
         if (diff <= log_tol) exit
      end do
      if (diff > log_tol) converged = .false.
   end subroutine cmp_logz_trunc

   pure real(dp) function cmp_logz_approx(lambda, nu)
      real(dp), intent(in) :: lambda, nu
      if (lambda <= 0.0_dp .or. nu <= 0.0_dp) then
         cmp_logz_approx = huge(1.0_dp)
      else
         cmp_logz_approx = nu*exp(log(lambda)/nu) &
            - (nu-1.0_dp)/(2.0_dp*nu)*log(lambda) &
            - 0.5_dp*(nu-1.0_dp)*log(2.0_dp*pi_dp) - 0.5_dp*log(nu)
      end if
   end function cmp_logz_approx

   function cmp_logz_hybrid(lambda, nu, control) result(logz)
      real(dp), intent(in) :: lambda, nu
      type(cmp_control_t), intent(in), optional :: control
      real(dp) :: logz
      type(cmp_control_t) :: ctrl
      integer :: last_y
      logical :: ok
      ctrl = cmp_control_t()
      if (present(control)) ctrl = control
      if (lambda == 0.0_dp) then
         logz = 0.0_dp
      else if (nu > 0.0_dp .and. -log(lambda)/nu < log(ctrl%hybrid_tol)) then
         logz = cmp_logz_approx(lambda,nu)
      else
         call cmp_logz_trunc(lambda,nu,ctrl%truncate_tol,ctrl%ymax,logz,last_y,ok)
      end if
   end function cmp_logz_hybrid

   function cmp_z_trunc(lambda, nu, control) result(z)
      real(dp), intent(in) :: lambda, nu
      type(cmp_control_t), intent(in), optional :: control
      real(dp) :: z, logz
      type(cmp_control_t) :: ctrl
      integer :: last_y
      logical :: ok
      ctrl=cmp_control_t(); if (present(control)) ctrl=control
      call cmp_logz_trunc(lambda,nu,ctrl%truncate_tol,ctrl%ymax,logz,last_y,ok)
      z=exp(logz)
   end function cmp_z_trunc

   pure real(dp) function cmp_z_approx(lambda, nu)
      real(dp), intent(in) :: lambda, nu
      cmp_z_approx=exp(cmp_logz_approx(lambda,nu))
   end function cmp_z_approx

   function cmp_z_hybrid(lambda,nu,control) result(z)
      real(dp), intent(in) :: lambda,nu
      type(cmp_control_t), intent(in), optional :: control
      real(dp) :: z
      z=exp(cmp_logz_hybrid(lambda,nu,control))
   end function cmp_z_hybrid

   function cmp_y_trunc(lambda,nu,control) result(ymax_used)
      real(dp), intent(in) :: lambda,nu
      type(cmp_control_t), intent(in), optional :: control
      integer :: ymax_used
      real(dp) :: logz
      type(cmp_control_t) :: ctrl
      logical :: ok
      ctrl=cmp_control_t(); if (present(control)) ctrl=control
      call cmp_logz_trunc(lambda,nu,ctrl%truncate_tol,ctrl%ymax,logz,ymax_used,ok)
   end function cmp_y_trunc

   pure real(dp) function term_log(lambda,nu,j)
      real(dp), intent(in) :: lambda,nu
      integer, intent(in) :: j
      if (j == 0) then
         term_log=0.0_dp
      else if (lambda <= 0.0_dp) then
         term_log=-huge(1.0_dp)
      else
         term_log=real(j,dp)*log(lambda)-nu*log_gamma(real(j+1,dp))
      end if
   end function term_log

   real(dp) function z_prodj(lambda,nu,maxj)
      real(dp), intent(in) :: lambda,nu
      integer, intent(in) :: maxj
      integer :: j
      z_prodj=0.0_dp
      do j=1,maxj
         z_prodj=z_prodj+real(j,dp)*exp(term_log(lambda,nu,j))
      end do
   end function z_prodj

   real(dp) function z_prodj2(lambda,nu,maxj)
      real(dp), intent(in) :: lambda,nu
      integer, intent(in) :: maxj
      integer :: j
      z_prodj2=0.0_dp
      do j=1,maxj
         z_prodj2=z_prodj2+real(j,dp)**2*exp(term_log(lambda,nu,j))
      end do
   end function z_prodj2

   real(dp) function z_prodjlogj(lambda,nu,maxj)
      real(dp), intent(in) :: lambda,nu
      integer, intent(in) :: maxj
      integer :: j
      real(dp) :: lf
      z_prodjlogj=0.0_dp
      do j=1,maxj
         lf=log_gamma(real(j+1,dp))
         z_prodjlogj=z_prodjlogj+real(j,dp)*lf*exp(term_log(lambda,nu,j))
      end do
   end function z_prodjlogj

   real(dp) function z_prodlogj(lambda,nu,maxj)
      real(dp), intent(in) :: lambda,nu
      integer, intent(in) :: maxj
      integer :: j
      real(dp) :: lf
      z_prodlogj=0.0_dp
      do j=1,maxj
         lf=log_gamma(real(j+1,dp))
         z_prodlogj=z_prodlogj+lf*exp(term_log(lambda,nu,j))
      end do
   end function z_prodlogj

   real(dp) function z_prodlogj2(lambda,nu,maxj)
      real(dp), intent(in) :: lambda,nu
      integer, intent(in) :: maxj
      integer :: j
      real(dp) :: lf
      z_prodlogj2=0.0_dp
      do j=1,maxj
         lf=log_gamma(real(j+1,dp))
         z_prodlogj2=z_prodlogj2+lf*lf*exp(term_log(lambda,nu,j))
      end do
   end function z_prodlogj2

end module compoissonreg_normalizer
