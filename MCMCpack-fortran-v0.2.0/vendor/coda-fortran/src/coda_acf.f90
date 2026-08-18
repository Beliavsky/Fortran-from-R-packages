! coda-fortran: computational translation of the R package coda.
! Original coda license: GPL (>= 2). This translation is GPL-2.0-or-later.
module coda_acf
   use coda_kinds, only : dp
   use coda_types, only : mcmc_chain, mcmc_list
   use coda_math, only : mean_vec, covariance_matrix
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   implicit none
   private
   public :: autocorr, autocorr_diag, autocorr_diag_list, crosscorr, rejection_rate, rejection_rate_list

contains

   function autocorr(chain, lags, relative) result(ac)
      type(mcmc_chain), intent(in) :: chain
      integer, intent(in) :: lags(:)
      logical, intent(in), optional :: relative
      real(dp), allocatable :: ac(:,:,:)
      real(dp), allocatable :: mu(:), g0(:)
      logical :: rel
      integer :: n, p, nl, h, lag_samples, i, j, t
      real(dp) :: num

      rel = .true.
      if (present(relative)) rel = relative
      n = chain%niter()
      p = chain%nvar()
      nl = size(lags)
      allocate(ac(nl,p,p), mu(p), g0(p))
      do j = 1, p
         mu(j) = mean_vec(chain%x(:,j))
         g0(j) = sum((chain%x(:,j) - mu(j))**2) / real(n,dp)
      end do
      do h = 1, nl
         if (rel) then
            lag_samples = lags(h)
         else
            if (mod(lags(h), chain%thin) /= 0) error stop "autocorr: lag incompatible with thinning"
            lag_samples = lags(h) / chain%thin
         end if
         do j = 1, p
            do i = 1, p
               if (lag_samples < 0 .or. lag_samples >= n .or. g0(i) <= 0.0_dp .or. g0(j) <= 0.0_dp) then
                  ac(h,i,j) = ieee_value(0.0_dp, ieee_quiet_nan)
               else
                  num = 0.0_dp
                  do t = 1, n - lag_samples
                     num = num + (chain%x(t + lag_samples,i) - mu(i)) * (chain%x(t,j) - mu(j))
                  end do
                  ac(h,i,j) = num / real(n,dp) / sqrt(g0(i) * g0(j))
               end if
            end do
         end do
      end do
   end function autocorr

   function autocorr_diag(chain, lags, relative) result(out)
      type(mcmc_chain), intent(in) :: chain
      integer, intent(in) :: lags(:)
      logical, intent(in), optional :: relative
      real(dp), allocatable :: out(:,:)
      real(dp), allocatable :: ac(:,:,:)
      integer :: j
      if (present(relative)) then
         ac = autocorr(chain, lags, relative)
      else
         ac = autocorr(chain, lags)
      end if
      allocate(out(size(lags), chain%nvar()))
      do j = 1, chain%nvar()
         out(:,j) = ac(:,j,j)
      end do
   end function autocorr_diag

   function autocorr_diag_list(lst, lags, relative) result(out)
      type(mcmc_list), intent(in) :: lst
      integer, intent(in) :: lags(:)
      logical, intent(in), optional :: relative
      real(dp), allocatable :: out(:,:), tmp(:,:)
      integer :: k
      allocate(out(size(lags), lst%nvar()))
      out = 0.0_dp
      do k = 1, lst%nchain()
         if (present(relative)) then
            tmp = autocorr_diag(lst%chain(k), lags, relative)
         else
            tmp = autocorr_diag(lst%chain(k), lags)
         end if
         out = out + tmp
      end do
      out = out / real(lst%nchain(),dp)
   end function autocorr_diag_list

   function crosscorr(chain) result(cor)
      type(mcmc_chain), intent(in) :: chain
      real(dp), allocatable :: cor(:,:)
      real(dp), allocatable :: cov(:,:)
      integer :: i, j, p
      cov = covariance_matrix(chain%x)
      p = chain%nvar()
      allocate(cor(p,p))
      do j = 1, p
         do i = 1, p
            if (cov(i,i) <= 0.0_dp .or. cov(j,j) <= 0.0_dp) then
               cor(i,j) = ieee_value(0.0_dp, ieee_quiet_nan)
            else
               cor(i,j) = cov(i,j) / sqrt(cov(i,i) * cov(j,j))
            end if
         end do
      end do
   end function crosscorr

   function rejection_rate(chain) result(rate)
      type(mcmc_chain), intent(in) :: chain
      real(dp), allocatable :: rate(:)
      integer :: j, n
      n = chain%niter()
      allocate(rate(chain%nvar()))
      if (n < 2) then
         rate = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      do j = 1, chain%nvar()
         rate(j) = real(count(chain%x(1:n-1,j) == chain%x(2:n,j)), dp) / real(n - 1,dp)
      end do
   end function rejection_rate

   function rejection_rate_list(lst) result(rate)
      type(mcmc_list), intent(in) :: lst
      real(dp), allocatable :: rate(:), tmp(:)
      integer :: k
      allocate(rate(lst%nvar()))
      rate = 0.0_dp
      do k = 1, lst%nchain()
         tmp = rejection_rate(lst%chain(k))
         rate = rate + tmp
      end do
      rate = rate / real(lst%nchain(),dp)
   end function rejection_rate_list

end module coda_acf
