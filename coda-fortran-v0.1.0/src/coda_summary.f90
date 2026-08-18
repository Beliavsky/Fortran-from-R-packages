! coda-fortran: computational translation of the R package coda.
! Original coda license: GPL (>= 2). This translation is GPL-2.0-or-later.
module coda_summary
   use coda_kinds, only : dp
   use coda_types, only : mcmc_chain, mcmc_list, pool_chains
   use coda_math, only : mean_vec, variance_vec, quantile_type7, sort_real
   use coda_spectrum, only : spectrum0_ar, spectrum_ar_result
   implicit none
   private

   type, public :: mcmc_summary
      real(dp), allocatable :: statistics(:,:) ! mean, sd, naive se, time-series se
      real(dp), allocatable :: quantiles(:,:)
      real(dp), allocatable :: probs(:)
      integer :: start = 1
      integer :: finish = 0
      integer :: thin = 1
      integer :: nchain = 1
   end type mcmc_summary

   public :: hpd_interval, batch_se, batch_se_list, summarize_mcmc, summarize_mcmc_list

contains

   function hpd_interval(chain, prob) result(ans)
      type(mcmc_chain), intent(in) :: chain
      real(dp), intent(in), optional :: prob
      real(dp), allocatable :: ans(:,:)
      real(dp), allocatable :: vals(:)
      real(dp) :: pr, width, best
      integer :: n, p, gap, j, i, best_i

      pr = 0.95_dp
      if (present(prob)) pr = prob
      if (pr <= 0.0_dp .or. pr >= 1.0_dp) error stop "hpd_interval: probability must be in (0,1)"
      n = chain%niter()
      p = chain%nvar()
      if (n < 2) error stop "hpd_interval: at least two samples required"
      gap = max(1, min(n - 1, nint(real(n,dp) * pr)))
      allocate(ans(p,2), vals(n))
      do j = 1, p
         vals = chain%x(:,j)
         call sort_real(vals)
         best = huge(1.0_dp)
         best_i = 1
         do i = 1, n - gap
            width = vals(i + gap) - vals(i)
            if (width < best) then
               best = width
               best_i = i
            end if
         end do
         ans(j,1) = vals(best_i)
         ans(j,2) = vals(best_i + gap)
      end do
   end function hpd_interval

   function batch_se(chain, batch_size) result(se)
      type(mcmc_chain), intent(in) :: chain
      integer, intent(in), optional :: batch_size
      real(dp), allocatable :: se(:), means(:,:)
      integer :: bs, nb, nuse, b, j
      real(dp) :: gm

      bs = 100
      if (present(batch_size)) bs = batch_size
      if (bs <= 0) error stop "batch_se: batch_size must be positive"
      nb = chain%niter() / bs
      if (nb < 2) error stop "batch_se: need at least two complete batches"
      nuse = nb * bs
      allocate(means(nb,chain%nvar()), se(chain%nvar()))
      do j = 1, chain%nvar()
         do b = 1, nb
            means(b,j) = mean_vec(chain%x((b-1)*bs+1:b*bs,j))
         end do
         gm = mean_vec(means(:,j))
         se(j) = sqrt(sum((means(:,j)-gm)**2) * real(bs,dp) / real(nb-1,dp)) / sqrt(real(chain%niter(),dp))
      end do
   end function batch_se

   function batch_se_list(lst, batch_size) result(se)
      type(mcmc_list), intent(in) :: lst
      integer, intent(in), optional :: batch_size
      real(dp), allocatable :: se(:), means(:,:)
      integer :: bs, nb, nuse, b, j, k, row
      real(dp) :: gm

      bs = 100
      if (present(batch_size)) bs = batch_size
      nb = lst%niter() / bs
      if (nb * lst%nchain() < 2) error stop "batch_se_list: insufficient batches"
      nuse = nb * bs
      allocate(means(nb*lst%nchain(),lst%nvar()), se(lst%nvar()))
      row = 0
      do k = 1, lst%nchain()
         do b = 1, nb
            row = row + 1
            do j = 1, lst%nvar()
               means(row,j) = mean_vec(lst%chain(k)%x((b-1)*bs+1:b*bs,j))
            end do
         end do
      end do
      do j = 1, lst%nvar()
         gm = mean_vec(means(:,j))
         se(j) = sqrt(sum((means(:,j)-gm)**2) * real(bs,dp) / real(size(means,1)-1,dp)) / &
                 sqrt(real(lst%niter()*lst%nchain(),dp))
      end do
   end function batch_se_list

   function summarize_mcmc(chain, probs) result(out)
      type(mcmc_chain), intent(in) :: chain
      real(dp), intent(in), optional :: probs(:)
      type(mcmc_summary) :: out
      real(dp), allocatable :: pp(:)
      type(spectrum_ar_result) :: sp
      real(dp), allocatable :: spec(:)
      integer :: j, k

      if (present(probs)) then
         pp = probs
      else
         pp = [0.025_dp, 0.25_dp, 0.5_dp, 0.75_dp, 0.975_dp]
      end if
      allocate(out%statistics(chain%nvar(),4), out%quantiles(chain%nvar(),size(pp)))
      out%probs = pp
      sp = spectrum0_ar(chain)
      spec = sp%spec
      do j = 1, chain%nvar()
         out%statistics(j,1) = mean_vec(chain%x(:,j))
         out%statistics(j,2) = sqrt(variance_vec(chain%x(:,j)))
         out%statistics(j,3) = sqrt(variance_vec(chain%x(:,j))/real(chain%niter(),dp))
         out%statistics(j,4) = sqrt(spec(j)/real(chain%niter(),dp))
         do k = 1, size(pp)
            out%quantiles(j,k) = quantile_type7(chain%x(:,j), pp(k))
         end do
      end do
      out%start = chain%start
      out%finish = chain%finish
      out%thin = chain%thin
      out%nchain = 1
   end function summarize_mcmc

   function summarize_mcmc_list(lst, probs) result(out)
      type(mcmc_list), intent(in) :: lst
      real(dp), intent(in), optional :: probs(:)
      type(mcmc_summary) :: out
      real(dp), allocatable :: pp(:), pooled(:,:), tsvar(:), s(:)
      type(spectrum_ar_result) :: sp
      integer :: j, k, c, ntot

      if (present(probs)) then
         pp = probs
      else
         pp = [0.025_dp, 0.25_dp, 0.5_dp, 0.75_dp, 0.975_dp]
      end if
      pooled = pool_chains(lst)
      ntot = size(pooled,1)
      allocate(tsvar(lst%nvar()))
      tsvar = 0.0_dp
      do c = 1, lst%nchain()
         sp = spectrum0_ar(lst%chain(c))
         s = sp%spec
         tsvar = tsvar + s
      end do
      tsvar = tsvar / real(lst%nchain(),dp)
      allocate(out%statistics(lst%nvar(),4), out%quantiles(lst%nvar(),size(pp)))
      out%probs = pp
      do j = 1, lst%nvar()
         out%statistics(j,1) = mean_vec(pooled(:,j))
         out%statistics(j,2) = sqrt(variance_vec(pooled(:,j)))
         out%statistics(j,3) = sqrt(variance_vec(pooled(:,j))/real(ntot,dp))
         out%statistics(j,4) = sqrt(tsvar(j)/real(ntot,dp))
         do k = 1, size(pp)
            out%quantiles(j,k) = quantile_type7(pooled(:,j), pp(k))
         end do
      end do
      out%start = lst%chain(1)%start
      out%finish = lst%chain(1)%finish
      out%thin = lst%chain(1)%thin
      out%nchain = lst%nchain()
   end function summarize_mcmc_list

end module coda_summary
