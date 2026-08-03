! SPDX-License-Identifier: Artistic-2.0
module mts_stats
   use mts_kinds, only : dp
   use mts_linalg, only : inverse_matrix, log_determinant
   use mts_types, only : mts_success, mts_invalid_input, mts_singular
   implicit none
   private

   real(dp), parameter :: pi = acos(-1.0_dp)

   public :: column_mean, center_columns, covariance_matrix, correlation_matrix
   public :: autocovariance_matrix, autocorrelation_matrix
   public :: normal_cdf, normal_logpdf, multivariate_normal_logpdf
   public :: student_t_logpdf, multivariate_student_t_logpdf
   public :: chi_square_cdf, chi_square_survival, regularized_gamma_p
   public :: ranks, sample_variance, determinant_ic

contains

   pure function column_mean(x) result(mu)
      real(dp), intent(in) :: x(:,:)
      real(dp) :: mu(size(x,2))
      if (size(x,1) > 0) then
         mu = sum(x,dim=1)/real(size(x,1),dp)
      else
         mu = 0.0_dp
      end if
   end function column_mean

   pure function center_columns(x,mu) result(y)
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(in), optional :: mu(:)
      real(dp) :: y(size(x,1),size(x,2))
      real(dp) :: m(size(x,2))
      integer :: i
      if (present(mu)) then
         m = mu
      else
         m = column_mean(x)
      end if
      do i = 1, size(x,1)
         y(i,:) = x(i,:)-m
      end do
   end function center_columns

   function covariance_matrix(x,unbiased) result(cov)
      real(dp), intent(in) :: x(:,:)
      logical, intent(in), optional :: unbiased
      real(dp) :: cov(size(x,2),size(x,2))
      real(dp) :: y(size(x,1),size(x,2)), denom
      logical :: use_unbiased
      use_unbiased = .true.
      if (present(unbiased)) use_unbiased = unbiased
      y = center_columns(x)
      if (use_unbiased) then
         denom = real(max(1,size(x,1)-1),dp)
      else
         denom = real(max(1,size(x,1)),dp)
      end if
      cov = matmul(transpose(y),y)/denom
   end function covariance_matrix

   function correlation_matrix(x) result(corr)
      real(dp), intent(in) :: x(:,:)
      real(dp) :: corr(size(x,2),size(x,2))
      real(dp) :: cov(size(x,2),size(x,2)), sd(size(x,2))
      integer :: i, j
      cov = covariance_matrix(x)
      sd = sqrt(max(0.0_dp,[(cov(i,i),i=1,size(cov,1))]))
      corr = 0.0_dp
      do i = 1, size(corr,1)
         do j = 1, size(corr,2)
            if (sd(i) > 0.0_dp .and. sd(j) > 0.0_dp) corr(i,j) = cov(i,j)/(sd(i)*sd(j))
         end do
         corr(i,i) = 1.0_dp
      end do
   end function correlation_matrix

   function autocovariance_matrix(x,lag,biased) result(gamma)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: lag
      logical, intent(in), optional :: biased
      real(dp) :: gamma(size(x,2),size(x,2))
      real(dp) :: y(size(x,1),size(x,2)), denom
      logical :: use_biased
      integer :: n
      n = size(x,1)
      gamma = 0.0_dp
      if (lag < 0 .or. lag >= n) return
      y = center_columns(x)
      use_biased = .true.
      if (present(biased)) use_biased = biased
      if (use_biased) then
         denom = real(n,dp)
      else
         denom = real(n-lag,dp)
      end if
      gamma = matmul(transpose(y(lag+1:n,:)),y(1:n-lag,:))/max(1.0_dp,denom)
   end function autocovariance_matrix

   function autocorrelation_matrix(x,lag) result(rho)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: lag
      real(dp) :: rho(size(x,2),size(x,2))
      real(dp) :: gamma0(size(x,2),size(x,2)), gammal(size(x,2),size(x,2))
      real(dp) :: sd(size(x,2))
      integer :: i, j
      gamma0 = autocovariance_matrix(x,0)
      gammal = autocovariance_matrix(x,lag)
      sd = sqrt(max(0.0_dp,[(gamma0(i,i),i=1,size(gamma0,1))]))
      rho = 0.0_dp
      do i = 1, size(rho,1)
         do j = 1, size(rho,2)
            if (sd(i) > 0.0_dp .and. sd(j) > 0.0_dp) rho(i,j) = gammal(i,j)/(sd(i)*sd(j))
         end do
      end do
   end function autocorrelation_matrix

   pure function sample_variance(x,unbiased) result(v)
      real(dp), intent(in) :: x(:)
      logical, intent(in), optional :: unbiased
      real(dp) :: v, mu, denom
      logical :: ub
      ub = .true.
      if (present(unbiased)) ub = unbiased
      if (size(x) < 1) then
         v = 0.0_dp
         return
      end if
      mu = sum(x)/real(size(x),dp)
      denom = real(size(x),dp)
      if (ub .and. size(x) > 1) denom = real(size(x)-1,dp)
      v = sum((x-mu)**2)/max(1.0_dp,denom)
   end function sample_variance

   pure function normal_cdf(x) result(p)
      real(dp), intent(in) :: x
      real(dp) :: p
      p = 0.5_dp*erfc(-x/sqrt(2.0_dp))
   end function normal_cdf

   pure function normal_logpdf(x,mu,sigma) result(value)
      real(dp), intent(in) :: x, mu, sigma
      real(dp) :: value
      if (sigma <= 0.0_dp) then
         value = -huge(1.0_dp)
      else
         value = -0.5_dp*log(2.0_dp*pi)-log(sigma)-0.5_dp*((x-mu)/sigma)**2
      end if
   end function normal_logpdf

   function multivariate_normal_logpdf(x,mu,sigma,status) result(value)
      real(dp), intent(in) :: x(:), mu(:), sigma(:,:)
      integer, intent(out), optional :: status
      real(dp) :: value
      real(dp), allocatable :: inv(:,:), diff(:)
      real(dp) :: logdet
      integer :: sgn, istat, k
      k = size(x)
      if (size(mu) /= k .or. size(sigma,1) /= k .or. size(sigma,2) /= k) then
         value = -huge(1.0_dp)
         if (present(status)) status = mts_invalid_input
         return
      end if
      call inverse_matrix(sigma,inv,istat)
      logdet = log_determinant(sigma,sgn,istat)
      if (istat /= mts_success .or. sgn <= 0) then
         value = -huge(1.0_dp)
         if (present(status)) status = mts_singular
         return
      end if
      diff = x-mu
      value = -0.5_dp*(real(k,dp)*log(2.0_dp*pi)+logdet+dot_product(diff,matmul(inv,diff)))
      if (present(status)) status = mts_success
   end function multivariate_normal_logpdf

   pure function student_t_logpdf(x,nu,variance_standardized) result(value)
      real(dp), intent(in) :: x, nu
      logical, intent(in), optional :: variance_standardized
      real(dp) :: value, scale2
      logical :: standardized
      standardized = .false.
      if (present(variance_standardized)) standardized = variance_standardized
      if (nu <= 0.0_dp) then
         value = -huge(1.0_dp)
         return
      end if
      scale2 = 1.0_dp
      if (standardized) then
         if (nu <= 2.0_dp) then
            value = -huge(1.0_dp)
            return
         end if
         scale2 = (nu-2.0_dp)/nu
      end if
      value = log_gamma(0.5_dp*(nu+1.0_dp))-log_gamma(0.5_dp*nu) &
              -0.5_dp*log(nu*pi*scale2)-0.5_dp*(nu+1.0_dp)*log(1.0_dp+x*x/(nu*scale2))
   end function student_t_logpdf

   function multivariate_student_t_logpdf(x,mu,sigma,nu,variance_standardized,status) result(value)
      real(dp), intent(in) :: x(:), mu(:), sigma(:,:), nu
      logical, intent(in), optional :: variance_standardized
      integer, intent(out), optional :: status
      real(dp) :: value, quad, scale
      real(dp), allocatable :: inv(:,:), diff(:)
      real(dp) :: logdet
      integer :: sgn, istat, k
      logical :: standardized
      standardized = .false.
      if (present(variance_standardized)) standardized = variance_standardized
      k = size(x)
      if (nu <= 0.0_dp .or. size(mu) /= k .or. size(sigma,1) /= k .or. size(sigma,2) /= k) then
         value = -huge(1.0_dp)
         if (present(status)) status = mts_invalid_input
         return
      end if
      scale = 1.0_dp
      if (standardized) then
         if (nu <= 2.0_dp) then
            value = -huge(1.0_dp)
            if (present(status)) status = mts_invalid_input
            return
         end if
         scale = (nu-2.0_dp)/nu
      end if
      call inverse_matrix(sigma,inv,istat)
      logdet = log_determinant(sigma,sgn,istat)
      if (istat /= mts_success .or. sgn <= 0) then
         value = -huge(1.0_dp)
         if (present(status)) status = mts_singular
         return
      end if
      diff = x-mu
      quad = dot_product(diff,matmul(inv,diff))/scale
      value = log_gamma(0.5_dp*(nu+real(k,dp)))-log_gamma(0.5_dp*nu) &
              -0.5_dp*real(k,dp)*log(nu*pi)-0.5_dp*logdet-0.5_dp*real(k,dp)*log(scale) &
              -0.5_dp*(nu+real(k,dp))*log(1.0_dp+quad/nu)
      if (present(status)) status = mts_success
   end function multivariate_student_t_logpdf

   function regularized_gamma_p(a,x) result(p)
      real(dp), intent(in) :: a, x
      real(dp) :: p
      real(dp) :: sumv, del, ap, b, c, d, h, an
      integer :: n
      if (a <= 0.0_dp .or. x < 0.0_dp) then
         p = 0.0_dp
      else if (x == 0.0_dp) then
         p = 0.0_dp
      else if (x < a+1.0_dp) then
         ap = a
         sumv = 1.0_dp/a
         del = sumv
         do n = 1, 10000
            ap = ap+1.0_dp
            del = del*x/ap
            sumv = sumv+del
            if (abs(del) <= abs(sumv)*1.0e-14_dp) exit
         end do
         p = sumv*exp(-x+a*log(x)-log_gamma(a))
      else
         b = x+1.0_dp-a
         c = 1.0_dp/tiny(1.0_dp)
         d = 1.0_dp/b
         h = d
         do n = 1, 10000
            an = -real(n,dp)*(real(n,dp)-a)
            b = b+2.0_dp
            d = an*d+b
            if (abs(d) < tiny(1.0_dp)) d = tiny(1.0_dp)
            c = b+an/c
            if (abs(c) < tiny(1.0_dp)) c = tiny(1.0_dp)
            d = 1.0_dp/d
            del = d*c
            h = h*del
            if (abs(del-1.0_dp) <= 1.0e-14_dp) exit
         end do
         p = 1.0_dp-exp(-x+a*log(x)-log_gamma(a))*h
      end if
      p = max(0.0_dp,min(1.0_dp,p))
   end function regularized_gamma_p

   function chi_square_cdf(x,degrees_freedom) result(p)
      real(dp), intent(in) :: x
      integer, intent(in) :: degrees_freedom
      real(dp) :: p
      if (degrees_freedom <= 0) then
         p = 0.0_dp
      else
         p = regularized_gamma_p(0.5_dp*real(degrees_freedom,dp),0.5_dp*max(0.0_dp,x))
      end if
   end function chi_square_cdf

   function chi_square_survival(x,degrees_freedom) result(p)
      real(dp), intent(in) :: x
      integer, intent(in) :: degrees_freedom
      real(dp) :: p
      p = max(0.0_dp,1.0_dp-chi_square_cdf(x,degrees_freedom))
   end function chi_square_survival

   function ranks(x) result(r)
      real(dp), intent(in) :: x(:)
      real(dp) :: r(size(x))
      integer, allocatable :: idx(:)
      integer :: i, j, k, first, last, tmp
      allocate(idx(size(x)))
      idx = [(i,i=1,size(x))]
      do i = 2, size(x)
         tmp = idx(i)
         j = i-1
         do while (j >= 1)
            if (x(idx(j)) <= x(tmp)) exit
            idx(j+1) = idx(j)
            j = j-1
         end do
         idx(j+1) = tmp
      end do
      r = 0.0_dp
      first = 1
      do while (first <= size(x))
         last = first
         do while (last < size(x))
            if (x(idx(last+1)) /= x(idx(first))) exit
            last = last+1
         end do
         do k = first, last
            r(idx(k)) = 0.5_dp*real(first+last,dp)
         end do
         first = last+1
      end do
   end function ranks

   subroutine determinant_ic(sigma,nobs,n_parameters,aic,bic,hq,status)
      real(dp), intent(in) :: sigma(:,:)
      integer, intent(in) :: nobs, n_parameters
      real(dp), intent(out) :: aic, bic, hq
      integer, intent(out), optional :: status
      real(dp) :: logdet
      integer :: sgn, istat
      logdet = log_determinant(sigma,sgn,istat)
      if (istat /= mts_success .or. sgn <= 0 .or. nobs < 1) then
         aic = huge(1.0_dp)
         bic = huge(1.0_dp)
         hq = huge(1.0_dp)
         istat = mts_invalid_input
      else
         aic = logdet+2.0_dp*real(n_parameters,dp)/real(nobs,dp)
         bic = logdet+log(real(nobs,dp))*real(n_parameters,dp)/real(nobs,dp)
         hq = logdet+2.0_dp*log(log(real(max(3,nobs),dp)))*real(n_parameters,dp)/real(nobs,dp)
      end if
      if (present(status)) status = istat
   end subroutine determinant_ic

end module mts_stats
