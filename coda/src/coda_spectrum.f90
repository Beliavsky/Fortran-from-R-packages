! coda-fortran: computational translation of the R package coda.
! Original coda license: GPL (>= 2). This translation is GPL-2.0-or-later.
module coda_spectrum
   use coda_kinds, only : dp
   use coda_types, only : mcmc_chain, mcmc_list
   use coda_math, only : mean_vec, variance_vec
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   implicit none
   private

   type, public :: spectrum_ar_result
      real(dp), allocatable :: spec(:)
      integer, allocatable :: order(:)
   end type spectrum_ar_result

   public :: spectrum0_ar, effective_size, effective_size_list, spectrum0

contains

   subroutine yule_walker_aic(x, phi, order, var_pred)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable, intent(out) :: phi(:)
      integer, intent(out) :: order
      real(dp), intent(out) :: var_pred
      real(dp), allocatable :: gamma(:), a(:), anew(:), vseq(:)
      real(dp) :: mu, kappa, v, score, best_score
      integer :: n, pmax, p, j, best_p

      n = size(x)
      pmax = min(n - 1, max(0, floor(10.0_dp * log10(real(max(n,2),dp)))))
      allocate(gamma(0:pmax), a(max(1,pmax)), anew(max(1,pmax)), vseq(0:pmax))
      mu = mean_vec(x)
      do p = 0, pmax
         gamma(p) = sum((x(1:n-p) - mu) * (x(1+p:n) - mu)) / real(n,dp)
      end do
      if (gamma(0) <= 0.0_dp) then
         order = 0
         var_pred = 0.0_dp
         allocate(phi(0))
         return
      end if

      v = gamma(0)
      vseq(0) = v
      a = 0.0_dp
      best_p = 0
      best_score = log(max(v, tiny(1.0_dp)))
      do p = 1, pmax
         kappa = gamma(p)
         do j = 1, p - 1
            kappa = kappa - a(j) * gamma(p-j)
         end do
         if (v <= tiny(1.0_dp)) exit
         kappa = kappa / v
         if (abs(kappa) >= 1.0_dp) exit
         anew = a
         anew(p) = kappa
         do j = 1, p - 1
            anew(j) = a(j) - kappa * a(p-j)
         end do
         a = anew
         v = v * (1.0_dp - kappa*kappa)
         if (v <= tiny(1.0_dp)) exit
         vseq(p) = v
         score = log(v) + 2.0_dp * real(p,dp) / real(n,dp)
         if (score < best_score) then
            best_score = score
            best_p = p
         end if
      end do

      if (best_p == 0) then
         order = 0
         var_pred = gamma(0)
         allocate(phi(0))
      else
         ! Recompute recursion only to the selected order.
         a = 0.0_dp
         v = gamma(0)
         do p = 1, best_p
            kappa = gamma(p)
            do j = 1, p - 1
               kappa = kappa - a(j) * gamma(p-j)
            end do
            kappa = kappa / v
            anew = a
            anew(p) = kappa
            do j = 1, p - 1
               anew(j) = a(j) - kappa * a(p-j)
            end do
            a = anew
            v = v * (1.0_dp - kappa*kappa)
         end do
         order = best_p
         var_pred = v
         allocate(phi(best_p))
         phi = a(1:best_p)
      end if
   end subroutine yule_walker_aic

   function spectrum0_ar(chain) result(out)
      type(mcmc_chain), intent(in) :: chain
      type(spectrum_ar_result) :: out
      real(dp), allocatable :: phi(:)
      real(dp) :: vp, denom
      integer :: j, ord

      allocate(out%spec(chain%nvar()), out%order(chain%nvar()))
      do j = 1, chain%nvar()
         if (variance_vec(chain%x(:,j)) == 0.0_dp) then
            out%spec(j) = 0.0_dp
            out%order(j) = 0
         else
            call yule_walker_aic(chain%x(:,j), phi, ord, vp)
            denom = 1.0_dp - sum(phi)
            if (abs(denom) <= sqrt(epsilon(1.0_dp))) then
               out%spec(j) = huge(1.0_dp)
            else
               out%spec(j) = vp / (denom * denom)
            end if
            out%order(j) = ord
         end if
      end do
   end function spectrum0_ar

   function effective_size(chain) result(ess)
      type(mcmc_chain), intent(in) :: chain
      real(dp), allocatable :: ess(:)
      type(spectrum_ar_result) :: sp
      integer :: j
      sp = spectrum0_ar(chain)
      allocate(ess(chain%nvar()))
      do j = 1, chain%nvar()
         if (sp%spec(j) == 0.0_dp) then
            ess(j) = 0.0_dp
         else
            ess(j) = real(chain%niter(),dp) * variance_vec(chain%x(:,j)) / sp%spec(j)
         end if
      end do
   end function effective_size

   function effective_size_list(lst) result(ess)
      type(mcmc_list), intent(in) :: lst
      real(dp), allocatable :: ess(:), tmp(:)
      integer :: k
      allocate(ess(lst%nvar()))
      ess = 0.0_dp
      do k = 1, lst%nchain()
         tmp = effective_size(lst%chain(k))
         ess = ess + tmp
      end do
   end function effective_size_list

   subroutine solve_small_normal(x, y, beta, ok)
      real(dp), intent(in) :: x(:,:), y(:)
      real(dp), intent(out) :: beta(size(x,2))
      logical, intent(out) :: ok
      real(dp), allocatable :: a(:,:), b(:)
      real(dp) :: pivot, fac, tmp
      integer :: p, i, j, k, ip

      p = size(x,2)
      allocate(a(p,p), b(p))
      a = matmul(transpose(x), x)
      b = matmul(transpose(x), y)
      ok = .true.
      do k = 1, p
         ip = k
         do i = k + 1, p
            if (abs(a(i,k)) > abs(a(ip,k))) ip = i
         end do
         if (abs(a(ip,k)) < 1.0e-14_dp) then
            ok = .false.
            beta = 0.0_dp
            return
         end if
         if (ip /= k) then
            do j = k, p
               tmp = a(k,j); a(k,j) = a(ip,j); a(ip,j) = tmp
            end do
            tmp = b(k); b(k) = b(ip); b(ip) = tmp
         end if
         pivot = a(k,k)
         do i = k + 1, p
            fac = a(i,k) / pivot
            a(i,k:p) = a(i,k:p) - fac * a(k,k:p)
            b(i) = b(i) - fac * b(k)
         end do
      end do
      beta = 0.0_dp
      do i = p, 1, -1
         beta(i) = (b(i) - sum(a(i,i+1:p) * beta(i+1:p))) / a(i,i)
      end do
   end subroutine solve_small_normal

   function spectrum0(chain, max_freq, order, max_length) result(spec_out)
      type(mcmc_chain), intent(in) :: chain
      real(dp), intent(in), optional :: max_freq
      integer, intent(in), optional :: order, max_length
      real(dp), allocatable :: spec_out(:)
      real(dp), allocatable :: z(:,:), y(:), freq(:), f1(:), f2(:), period(:), design(:,:), beta(:), eta(:), mu(:), work(:)
      real(dp) :: mf, angle, pred
      complex(dp) :: dft
      integer :: ord, ml, batch, n, nb, p, j, k, t, nf, m, it
      logical :: ok

      mf = 0.5_dp
      if (present(max_freq)) mf = max_freq
      ord = 1
      if (present(order)) ord = order
      if (ord < 0 .or. ord > 2) error stop "spectrum0: order must be 0, 1, or 2"
      ml = 200
      if (present(max_length)) ml = max_length
      if (ml < 2) error stop "spectrum0: max_length must be >= 2"

      n = chain%niter()
      p = chain%nvar()
      if (n > ml) then
         batch = ceiling(real(n,dp) / real(ml,dp))
         nb = n / batch
         if (nb < 2) nb = 2
         allocate(z(nb,p))
         do j = 1, p
            do k = 1, nb
               z(k,j) = mean_vec(chain%x((k-1)*batch+1:min(k*batch,n),j))
            end do
         end do
      else
         batch = 1
         nb = n
         z = chain%x
      end if

      nf = nb / 2
      allocate(freq(nf), f1(nf), f2(nf), period(nf), spec_out(p))
      do k = 1, nf
         freq(k) = real(k,dp) / real(nb,dp)
         f1(k) = sqrt(3.0_dp) * (4.0_dp*freq(k) - 1.0_dp)
         f2(k) = sqrt(5.0_dp) * (24.0_dp*freq(k)**2 - 12.0_dp*freq(k) + 1.0_dp)
      end do
      m = count(freq <= mf)
      if (m < ord + 1) error stop "spectrum0: too few frequencies for fit"
      allocate(design(m,ord+1), beta(ord+1), eta(m), mu(m), work(m))
      design(:,1) = 1.0_dp
      if (ord >= 1) design(:,2) = f1(1:m)
      if (ord >= 2) design(:,3) = f2(1:m)

      do j = 1, p
         y = z(:,j)
         if (variance_vec(y) == 0.0_dp) then
            spec_out(j) = 0.0_dp
            cycle
         end if
         do k = 1, nf
            dft = (0.0_dp, 0.0_dp)
            do t = 1, nb
               angle = -2.0_dp * acos(-1.0_dp) * real(k*t,dp) / real(nb,dp)
               dft = dft + cmplx(y(t)*cos(angle), y(t)*sin(angle), kind=dp)
            end do
            period(k) = real(dft * conjg(dft),dp) / real(nb,dp)
         end do
         beta = 0.0_dp
         beta(1) = log(max(mean_vec(period(1:m)), tiny(1.0_dp)))
         do it = 1, 100
            eta = matmul(design, beta)
            mu = exp(min(eta, log(huge(1.0_dp))*0.25_dp))
            work = eta + (period(1:m) - mu) / max(mu, tiny(1.0_dp))
            call solve_small_normal(design, work, beta, ok)
            if (.not. ok) exit
            if (maxval(abs(matmul(design,beta) - eta)) < 1.0e-10_dp) exit
         end do
         pred = beta(1)
         if (ord >= 1) pred = pred - sqrt(3.0_dp) * beta(2)
         if (ord >= 2) pred = pred + sqrt(5.0_dp) * beta(3)
         spec_out(j) = exp(pred) * real(batch,dp)
      end do
   end function spectrum0

end module coda_spectrum
