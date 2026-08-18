! coda-fortran: computational translation of the R package coda.
! Original coda license: GPL (>= 2). This translation is GPL-2.0-or-later.
module coda_diagnostics
   use coda_kinds, only : dp
   use coda_types, only : mcmc_chain, mcmc_list, window_mcmc, make_mcmc
   use coda_math, only : mean_vec, variance_vec, covariance_matrix, sample_covariance, &
                         quantile_type7, normal_quantile, f_quantile, cholesky_lower, &
                         symmetric_max_eigenvalue, cramer_cdf, is_finite
   use coda_spectrum, only : spectrum0_ar, spectrum_ar_result
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   implicit none
   private

   type, public :: geweke_result
      real(dp), allocatable :: z(:)
      real(dp) :: frac1 = 0.1_dp
      real(dp) :: frac2 = 0.5_dp
   end type geweke_result

   type, public :: gelman_result
      real(dp), allocatable :: psrf(:,:) ! point estimate, upper CI
      real(dp) :: mpsrf = 0.0_dp
      logical :: has_mpsrf = .false.
   end type gelman_result

   type, public :: heidel_result
      logical, allocatable :: stationarity_pass(:)
      integer, allocatable :: start(:)
      real(dp), allocatable :: pvalue(:)
      logical, allocatable :: halfwidth_pass(:)
      real(dp), allocatable :: mean(:)
      real(dp), allocatable :: halfwidth(:)
   end type heidel_result

   type, public :: raftery_result
      real(dp) :: q = 0.025_dp
      real(dp) :: r = 0.005_dp
      real(dp) :: s = 0.95_dp
      integer :: nmin = 0
      logical :: enough_samples = .true.
      integer, allocatable :: burn_in(:)
      integer, allocatable :: total(:)
      real(dp), allocatable :: dependence_factor(:)
   end type raftery_result

   public :: geweke_diag, gelman_diag, heidel_diag, raftery_diag

contains

   function geweke_diag(chain, frac1, frac2) result(out)
      type(mcmc_chain), intent(in) :: chain
      real(dp), intent(in), optional :: frac1, frac2
      type(geweke_result) :: out
      real(dp) :: f1, f2
      integer :: s1, e1, s2, e2, j
      type(mcmc_chain) :: y1, y2
      type(spectrum_ar_result) :: sp1, sp2
      real(dp) :: v1, v2

      f1 = 0.1_dp
      f2 = 0.5_dp
      if (present(frac1)) f1 = frac1
      if (present(frac2)) f2 = frac2
      if (f1 < 0.0_dp .or. f1 > 1.0_dp .or. f2 < 0.0_dp .or. f2 > 1.0_dp) &
         error stop "geweke_diag: invalid fraction"
      if (f1 + f2 > 1.0_dp) error stop "geweke_diag: windows overlap"

      s1 = chain%start
      e1 = ceiling(real(chain%start,dp) + f1 * real(chain%finish - chain%start,dp))
      s2 = floor(real(chain%finish,dp) - f2 * real(chain%finish - chain%start,dp))
      e2 = chain%finish
      y1 = window_mcmc(chain, start=s1, finish=e1)
      y2 = window_mcmc(chain, start=s2, finish=e2)
      sp1 = spectrum0_ar(y1)
      sp2 = spectrum0_ar(y2)
      allocate(out%z(chain%nvar()))
      do j = 1, chain%nvar()
         v1 = sp1%spec(j) / real(y1%niter(),dp)
         v2 = sp2%spec(j) / real(y2%niter(),dp)
         if (v1 + v2 <= 0.0_dp) then
            out%z(j) = ieee_value(0.0_dp, ieee_quiet_nan)
         else
            out%z(j) = (mean_vec(y1%x(:,j)) - mean_vec(y2%x(:,j))) / sqrt(v1 + v2)
         end if
      end do
      out%frac1 = f1
      out%frac2 = f2
   end function geweke_diag

   subroutine transform_gelman_data(lst, x)
      type(mcmc_list), intent(in) :: lst
      real(dp), allocatable, intent(out) :: x(:,:,:)
      real(dp) :: zmin, zmax
      integer :: n, p, m, i, j, k
      n = lst%niter(); p = lst%nvar(); m = lst%nchain()
      allocate(x(n,p,m))
      do k = 1, m
         x(:,:,k) = lst%chain(k)%x
      end do
      do j = 1, p
         zmin = minval(x(:,j,:))
         zmax = maxval(x(:,j,:))
         if (zmin > 0.0_dp) then
            if (zmax < 1.0_dp) then
               do k = 1, m
                  do i = 1, n
                     x(i,j,k) = log(x(i,j,k)/(1.0_dp-x(i,j,k)))
                  end do
               end do
            else
               x(:,j,:) = log(x(:,j,:))
            end if
         end if
      end do
   end subroutine transform_gelman_data

   subroutine solve_lower_vec(l, b, x)
      real(dp), intent(in) :: l(:,:), b(:)
      real(dp), intent(out) :: x(size(b))
      integer :: i, k
      real(dp) :: s
      do i = 1, size(b)
         s = b(i)
         do k = 1, i-1
            s = s - l(i,k)*x(k)
         end do
         x(i) = s/l(i,i)
      end do
   end subroutine solve_lower_vec

   subroutine solve_upper_vec(l, b, x)
      real(dp), intent(in) :: l(:,:), b(:)
      real(dp), intent(out) :: x(size(b))
      integer :: i, k, n
      real(dp) :: s
      n = size(b)
      do i = n, 1, -1
         s = b(i)
         do k = i+1, n
            s = s - l(k,i)*x(k)
         end do
         x(i) = s/l(i,i)
      end do
   end subroutine solve_upper_vec

   function gelman_diag(lst, confidence, transform, autoburnin, multivariate) result(out)
      type(mcmc_list), intent(in) :: lst
      real(dp), intent(in), optional :: confidence
      logical, intent(in), optional :: transform, autoburnin, multivariate
      type(gelman_result) :: out
      type(mcmc_list) :: work
      type(mcmc_chain), allocatable :: trimmed(:)
      real(dp), allocatable :: x(:,:,:), s2(:,:,:), wmat(:,:), bmat(:,:), xbar(:,:), &
                               w(:), b(:), vars2(:), varw(:), varb(:), covwb(:), muhat(:), &
                               v(:), varv(:), dfv(:), dfadj(:), wdf(:), r2rand(:), l(:,:), &
                               cmat(:,:), tmpv(:), tmp2(:), y(:), xx(:)
      real(dp) :: conf, r2fixed, fq, emax
      logical :: tr, ab, mv
      integer :: n, p, m, j, k, info, i

      conf = 0.95_dp; tr = .false.; ab = .true.; mv = .true.
      if (present(confidence)) conf = confidence
      if (present(transform)) tr = transform
      if (present(autoburnin)) ab = autoburnin
      if (present(multivariate)) mv = multivariate
      if (lst%nchain() < 2) error stop "gelman_diag: at least two chains required"

      work = lst
      if (ab .and. work%chain(1)%start < work%chain(1)%finish/2) then
         allocate(trimmed(work%nchain()))
         do k = 1, work%nchain()
            trimmed(k) = window_mcmc(work%chain(k), start=work%chain(k)%finish/2 + 1)
         end do
         work%chain = trimmed
      end if
      n = work%niter(); p = work%nvar(); m = work%nchain()
      if (n < 2) error stop "gelman_diag: chains too short"
      if (tr) then
         call transform_gelman_data(work, x)
      else
         allocate(x(n,p,m))
         do k = 1, m
            x(:,:,k) = work%chain(k)%x
         end do
      end if

      allocate(s2(p,p,m), wmat(p,p), bmat(p,p), xbar(p,m))
      do k = 1, m
         s2(:,:,k) = covariance_matrix(x(:,:,k))
         do j = 1, p
            xbar(j,k) = mean_vec(x(:,j,k))
         end do
      end do
      wmat = sum(s2, dim=3)/real(m,dp)
      bmat = real(n,dp) * covariance_matrix(transpose(xbar))

      if (p > 1 .and. mv) then
         allocate(l(p,p), cmat(p,p), tmpv(p), tmp2(p))
         call cholesky_lower(wmat, l, info)
         if (info == 0) then
            ! C = L^{-1} B L^{-T}, formed column by column.
            do j = 1, p
               call solve_lower_vec(l, bmat(:,j), tmpv)
               cmat(:,j) = tmpv
            end do
            do i = 1, p
               call solve_lower_vec(l, cmat(i,:), tmpv)
               cmat(i,:) = tmpv
            end do
            cmat = 0.5_dp*(cmat + transpose(cmat))
            emax = symmetric_max_eigenvalue(cmat)
            out%mpsrf = sqrt((1.0_dp - 1.0_dp/real(n,dp)) + &
                              (1.0_dp + 1.0_dp/real(p,dp))*emax/real(n,dp))
            out%has_mpsrf = .true.
         end if
      end if

      allocate(w(p), b(p), varw(p), varb(p), covwb(p), muhat(p), v(p), varv(p), dfv(p), &
               dfadj(p), wdf(p), r2rand(p), out%psrf(p,2), y(m), xx(m))
      do j = 1, p
         w(j) = wmat(j,j)
         b(j) = bmat(j,j)
         do k = 1, m
            y(k) = s2(j,j,k)
            xx(k) = xbar(j,k)
         end do
         muhat(j) = mean_vec(xx)
         varw(j) = variance_vec(y)/real(m,dp)
         varb(j) = 2.0_dp*b(j)*b(j)/real(m-1,dp)
         covwb(j) = real(n,dp)/real(m,dp) * &
                    (sample_covariance(y, xx*xx) - 2.0_dp*muhat(j)*sample_covariance(y,xx))
         v(j) = real(n-1,dp)*w(j)/real(n,dp) + (1.0_dp + 1.0_dp/real(m,dp))*b(j)/real(n,dp)
         varv(j) = (real(n-1,dp)**2*varw(j) + (1.0_dp + 1.0_dp/real(m,dp))**2*varb(j) + &
                    2.0_dp*real(n-1,dp)*(1.0_dp + 1.0_dp/real(m,dp))*covwb(j))/real(n,dp)**2
         if (varv(j) <= 0.0_dp .or. w(j) <= 0.0_dp .or. varw(j) <= 0.0_dp) then
            out%psrf(j,:) = ieee_value(0.0_dp, ieee_quiet_nan)
         else
            dfv(j) = 2.0_dp*v(j)*v(j)/varv(j)
            dfadj(j) = (dfv(j)+3.0_dp)/(dfv(j)+1.0_dp)
            wdf(j) = 2.0_dp*w(j)*w(j)/varw(j)
            r2fixed = real(n-1,dp)/real(n,dp)
            r2rand(j) = (1.0_dp + 1.0_dp/real(m,dp))/real(n,dp)*(b(j)/w(j))
            out%psrf(j,1) = sqrt(dfadj(j)*(r2fixed+r2rand(j)))
            fq = f_quantile((1.0_dp+conf)/2.0_dp, real(m-1,dp), wdf(j))
            out%psrf(j,2) = sqrt(dfadj(j)*(r2fixed+fq*r2rand(j)))
         end if
      end do
   end function gelman_diag

   function heidel_diag(chain, eps, pvalue) result(out)
      type(mcmc_chain), intent(in) :: chain
      real(dp), intent(in), optional :: eps, pvalue
      type(heidel_result) :: out
      real(dp) :: ep, pv, s0, s0ci, ybar, halfw, stat, pcr
      integer :: n, p, j, step, idx, start_idx, i
      type(mcmc_chain) :: seg, second
      type(spectrum_ar_result) :: sp
      real(dp), allocatable :: bcum(:)
      logical :: conv

      ep = 0.1_dp; pv = 0.05_dp
      if (present(eps)) ep = eps
      if (present(pvalue)) pv = pvalue
      n = chain%niter(); p = chain%nvar()
      allocate(out%stationarity_pass(p), out%start(p), out%pvalue(p), out%halfwidth_pass(p), &
               out%mean(p), out%halfwidth(p))
      step = max(1, n/10)
      do j = 1, p
         second = make_mcmc(reshape(chain%x(max(1,n/2):n,j), [n-max(1,n/2)+1,1]), &
                            start=chain%start+(max(1,n/2)-1)*chain%thin, thin=chain%thin)
         sp = spectrum0_ar(second)
         s0 = sp%spec(1)
         conv = .false.
         stat = ieee_value(0.0_dp, ieee_quiet_nan)
         start_idx = 1
         do idx = 1, max(1,n/2), step
            seg = make_mcmc(reshape(chain%x(idx:n,j), [n-idx+1,1]), &
                            start=chain%start+(idx-1)*chain%thin, thin=chain%thin)
            ybar = mean_vec(seg%x(:,1))
            allocate(bcum(seg%niter()))
            bcum(1) = seg%x(1,1)-ybar
            do i = 2, seg%niter()
               bcum(i) = bcum(i-1) + seg%x(i,1)-ybar
            end do
            if (s0 > 0.0_dp) then
               stat = sum((bcum*bcum)/(real(seg%niter(),dp)*s0))/real(seg%niter(),dp)
               pcr = cramer_cdf(stat)
               if (is_finite(stat) .and. pcr < 1.0_dp-pv) then
                  conv = .true.
                  start_idx = idx
                  deallocate(bcum)
                  exit
               end if
            end if
            deallocate(bcum)
         end do
         if (conv) then
            seg = make_mcmc(reshape(chain%x(start_idx:n,j), [n-start_idx+1,1]), &
                            start=chain%start+(start_idx-1)*chain%thin, thin=chain%thin)
            sp = spectrum0_ar(seg)
            s0ci = sp%spec(1)
            ybar = mean_vec(seg%x(:,1))
            halfw = 1.96_dp*sqrt(s0ci/real(seg%niter(),dp))
            out%stationarity_pass(j) = .true.
            out%start(j) = seg%start
            out%pvalue(j) = 1.0_dp-cramer_cdf(stat)
            out%halfwidth_pass(j) = is_finite(halfw) .and. abs(halfw/ybar) <= ep
            out%mean(j) = ybar
            out%halfwidth(j) = halfw
         else
            out%stationarity_pass(j) = .false.
            out%start(j) = -huge(1)
            out%pvalue(j) = ieee_value(0.0_dp, ieee_quiet_nan)
            out%halfwidth_pass(j) = .false.
            out%mean(j) = ieee_value(0.0_dp, ieee_quiet_nan)
            out%halfwidth(j) = ieee_value(0.0_dp, ieee_quiet_nan)
         end if
      end do
   end function heidel_diag

   function raftery_diag(chain, q, r, s, converge_eps) result(out)
      type(mcmc_chain), intent(in) :: chain
      real(dp), intent(in), optional :: q, r, s, converge_eps
      type(raftery_result) :: out
      real(dp) :: qq, rr, ss, ce, phi, quant, bic, g2, fitted, alpha, beta, tempburn, tempprec
      logical, allocatable :: dich(:), test(:)
      integer :: n, p, j, factor, newdim, i, i1, i2, i3, nburn, nkeep
      integer :: tri(2,2,2), tran(2,2), a, b, c, denom

      qq=0.025_dp; rr=0.005_dp; ss=0.95_dp; ce=0.001_dp
      if (present(q)) qq=q
      if (present(r)) rr=r
      if (present(s)) ss=s
      if (present(converge_eps)) ce=converge_eps
      out%q=qq; out%r=rr; out%s=ss
      phi = normal_quantile(0.5_dp*(1.0_dp+ss))
      out%nmin = ceiling((qq*(1.0_dp-qq)*phi*phi)/(rr*rr))
      n=chain%niter(); p=chain%nvar()
      allocate(out%burn_in(p), out%total(p), out%dependence_factor(p))
      if (out%nmin > n) then
         out%enough_samples=.false.
         out%burn_in=0; out%total=0; out%dependence_factor=ieee_value(0.0_dp,ieee_quiet_nan)
         return
      end if
      out%enough_samples=.true.
      do j=1,p
         quant=quantile_type7(chain%x(:,j),qq)
         allocate(dich(n))
         dich=chain%x(:,j) <= quant
         factor=0; bic=1.0_dp
         do while (bic >= 0.0_dp)
            factor=factor+1
            newdim=1+(n-1)/factor
            if (newdim < 4) exit
            allocate(test(newdim))
            do i=1,newdim
               test(i)=dich(1+(i-1)*factor)
            end do
            tri=0
            do i=1,newdim-2
               a=merge(2,1,test(i)); b=merge(2,1,test(i+1)); c=merge(2,1,test(i+2))
               tri(a,b,c)=tri(a,b,c)+1
            end do
            g2=0.0_dp
            do i1=1,2; do i2=1,2; do i3=1,2
               if (tri(i1,i2,i3) /= 0) then
                  denom=sum(tri(:,i2,:))
                  if (denom > 0) then
                     fitted=real(sum(tri(i1,i2,:))*sum(tri(:,i2,i3)),dp)/real(denom,dp)
                     if (fitted > 0.0_dp) g2=g2+2.0_dp*real(tri(i1,i2,i3),dp)*log(real(tri(i1,i2,i3),dp)/fitted)
                  end if
               end if
            end do; end do; end do
            bic=g2-2.0_dp*log(real(newdim-2,dp))
            if (bic >= 0.0_dp) deallocate(test)
         end do
         if (.not. allocated(test)) then
            allocate(test(newdim))
            do i=1,newdim
               test(i)=dich(1+(i-1)*factor)
            end do
         end if
         tran=0
         do i=1,newdim-1
            a=merge(2,1,test(i)); b=merge(2,1,test(i+1)); tran(a,b)=tran(a,b)+1
         end do
         if (sum(tran(1,:)) == 0 .or. sum(tran(2,:)) == 0) then
            out%burn_in(j)=0; out%total(j)=0; out%dependence_factor(j)=ieee_value(0.0_dp,ieee_quiet_nan)
         else
            alpha=real(tran(1,2),dp)/real(sum(tran(1,:)),dp)
            beta=real(tran(2,1),dp)/real(sum(tran(2,:)),dp)
            if (alpha <= 0.0_dp .or. beta <= 0.0_dp .or. abs(1.0_dp-alpha-beta) <= tiny(1.0_dp)) then
               out%burn_in(j)=0; out%total(j)=0; out%dependence_factor(j)=ieee_value(0.0_dp,ieee_quiet_nan)
            else
               tempburn=log((ce*(alpha+beta))/max(alpha,beta))/log(abs(1.0_dp-alpha-beta))
               nburn=ceiling(tempburn)*factor*chain%thin
               tempprec=((2.0_dp-alpha-beta)*alpha*beta*phi*phi)/(((alpha+beta)**3)*rr*rr)
               nkeep=ceiling(tempprec)*factor*chain%thin
               out%burn_in(j)=nburn
               out%total(j)=nburn+nkeep
               out%dependence_factor(j)=real(nburn+nkeep,dp)/real(out%nmin,dp)
            end if
         end if
         deallocate(dich,test)
      end do
   end function raftery_diag

end module coda_diagnostics
