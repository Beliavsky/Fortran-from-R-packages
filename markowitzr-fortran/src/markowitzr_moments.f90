! SPDX-License-Identifier: LGPL-3.0-or-later
! Based on MarkowitzR, copyright 2014-2020 Steven E. Pav.
module markowitzr_moments
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use markowitzr_kinds, only: dp
   use markowitzr_types, only: theta_result
   use markowitzr_linalg, only: column_mean, sample_covariance, covariance_of_mean
   use markowitzr_linalg, only: symmetric_ivech, symmetric_vech, invert_matrix
   use markowitzr_linalg, only: kronecker_product, duplication_matrix
   use markowitzr_linalg, only: lower_vector_indices, symmetrize_matrix
   implicit none
   private

   integer, parameter, public :: covariance_empirical = 0
   integer, parameter, public :: covariance_normal = 1
   integer, parameter, public :: covariance_hac = 2

   abstract interface
      subroutine moment_covariance_callback(z, covariance, status)
         import dp
         real(dp), intent(in) :: z(:, :)
         real(dp), intent(out) :: covariance(:, :)
         integer, intent(out) :: status
      end subroutine moment_covariance_callback
   end interface

   public :: theta_vcov, itheta_vcov, moment_covariance_callback
   public :: hac_covariance_of_mean

contains

   function theta_vcov(x, fit_intercept, covariance_method, hac_lags, &
                       covariance_callback) result(out)
      real(dp), intent(in) :: x(:, :)
      logical, intent(in), optional :: fit_intercept
      integer, intent(in), optional :: covariance_method, hac_lags
      procedure(moment_covariance_callback), optional :: covariance_callback
      type(theta_result) :: out
      real(dp), allocatable :: clean_x(:, :), y(:, :), z(:, :)
      real(dp), allocatable :: sigma(:, :), mean_x(:)
      integer, allocatable :: pair_i(:), pair_j(:)
      logical :: use_intercept
      integer :: method, n, p, pp, q, i, j, k, a, b, callback_status

      use_intercept = .true.
      if (present(fit_intercept)) use_intercept = fit_intercept
      method = covariance_empirical
      if (present(covariance_method)) method = covariance_method

      call remove_incomplete_rows(x,clean_x)
      n = size(clean_x,1)
      p = size(clean_x,2)
      pp = p+merge(1,0,use_intercept)
      q = pp*(pp+1)/2
      out%n = n
      out%pp = pp

      if (p < 1) then
         out%status = 1
         out%message = 'x must contain at least one column'
         return
      end if
      if (n < 2) then
         out%status = 2
         out%message = 'at least two complete observations are required'
         return
      end if
      if (method < covariance_empirical .or. method > covariance_hac) then
         out%status = 3
         out%message = 'unknown covariance method'
         return
      end if
      if (method == covariance_normal .and. .not. use_intercept) then
         ! The Fortran port supports this case even though upstream R rejects it.
      end if

      allocate(y(n,pp),z(n,q),pair_i(q),pair_j(q))
      if (use_intercept) then
         y(:,1) = 1.0_dp
         y(:,2:) = clean_x
      else
         y = clean_x
      end if

      k = 0
      do j = 1, pp
         do i = j, pp
            k = k+1
            pair_i(k) = i
            pair_j(k) = j
            z(:,k) = y(:,i)*y(:,j)
         end do
      end do

      allocate(out%mu(q),out%covariance(q,q))
      out%mu = column_mean(z)

      if (present(covariance_callback)) then
         call covariance_callback(z,out%covariance,callback_status)
         if (callback_status /= 0) then
            out%status = 4
            out%message = 'covariance callback failed'
            return
         end if
      else
         select case (method)
         case (covariance_empirical)
            out%covariance = covariance_of_mean(z)
         case (covariance_normal)
            mean_x = column_mean(clean_x)
            sigma = sample_covariance(clean_x)
            do b = 1, q
               do a = b, q
                  out%covariance(a,b) = normal_product_covariance( &
                     pair_i(a),pair_j(a),pair_i(b),pair_j(b), &
                     mean_x,sigma,use_intercept)/real(n,dp)
                  out%covariance(b,a) = out%covariance(a,b)
               end do
            end do
         case (covariance_hac)
            if (present(hac_lags)) then
               out%covariance = hac_covariance_of_mean(z,hac_lags)
            else
               out%covariance = hac_covariance_of_mean(z)
            end if
         end select
      end if
      call symmetrize_matrix(out%covariance)
   end function theta_vcov

   function itheta_vcov(x, fit_intercept, covariance_method, hac_lags, &
                        covariance_callback) result(out)
      real(dp), intent(in) :: x(:, :)
      logical, intent(in), optional :: fit_intercept
      integer, intent(in), optional :: covariance_method, hac_lags
      procedure(moment_covariance_callback), optional :: covariance_callback
      type(theta_result) :: out
      type(theta_result) :: sm
      real(dp), allocatable :: theta(:, :), inverse_theta(:, :)
      real(dp), allocatable :: kron(:, :), d(:, :), selected(:, :), h(:, :)
      integer, allocatable :: lower_indices(:)
      logical :: use_intercept
      integer :: method, lags, status, q, i, pp

      use_intercept = .true.
      if (present(fit_intercept)) use_intercept = fit_intercept
      method = covariance_empirical
      if (present(covariance_method)) method = covariance_method
      lags = -1
      if (present(hac_lags)) lags = hac_lags

      if (present(covariance_callback)) then
         sm = theta_vcov(x,use_intercept,method,lags,covariance_callback)
      else
         sm = theta_vcov(x,use_intercept,method,lags)
      end if
      out%n = sm%n
      out%pp = sm%pp
      if (sm%status /= 0) then
         out%status = sm%status
         out%message = sm%message
         return
      end if

      theta = symmetric_ivech(sm%mu,status)
      if (status /= 0) then
         out%status = 10
         out%message = 'invalid packed second moment'
         return
      end if
      pp = size(theta,1)
      allocate(inverse_theta(pp,pp))
      call invert_matrix(theta,inverse_theta,status)
      if (status /= 0) then
         out%status = 11
         out%message = 'second moment matrix is singular'
         return
      end if

      q = pp*(pp+1)/2
      allocate(out%mu(q),out%covariance(q,q))
      out%mu = symmetric_vech(inverse_theta)
      kron = kronecker_product(inverse_theta,inverse_theta)
      d = duplication_matrix(pp)
      lower_indices = lower_vector_indices(pp)
      allocate(selected(q,pp*pp),h(q,q))
      do i = 1, q
         selected(i,:) = kron(lower_indices(i),:)
      end do
      h = -matmul(selected,d)
      out%covariance = matmul(h,matmul(sm%covariance,transpose(h)))
      call symmetrize_matrix(out%covariance)
   end function itheta_vcov

   function hac_covariance_of_mean(z, lags) result(covariance)
      real(dp), intent(in) :: z(:, :)
      integer, intent(in), optional :: lags
      real(dp) :: covariance(size(z,2),size(z,2))
      real(dp), allocatable :: centered(:, :), gamma(:, :)
      real(dp) :: mean_z(size(z,2)), weight
      integer :: n, q, lag_count, lag, t

      n = size(z,1)
      q = size(z,2)
      covariance = 0.0_dp
      if (n < 2 .or. q < 1) return
      if (present(lags)) then
         lag_count = max(0,min(lags,n-1))
      else
         lag_count = int(4.0_dp*(real(n,dp)/100.0_dp)**(2.0_dp/9.0_dp))
         lag_count = max(0,min(lag_count,n-1))
      end if

      mean_z = column_mean(z)
      allocate(centered(n,q),gamma(q,q))
      centered = z-spread(mean_z,1,n)
      do t = 1, n
         covariance = covariance+outer_product(centered(t,:),centered(t,:))
      end do
      covariance = covariance/real(n,dp)

      do lag = 1, lag_count
         gamma = 0.0_dp
         do t = lag+1, n
            gamma = gamma+outer_product(centered(t,:),centered(t-lag,:))
         end do
         gamma = gamma/real(n,dp)
         weight = 1.0_dp-real(lag,dp)/real(lag_count+1,dp)
         covariance = covariance+weight*(gamma+transpose(gamma))
      end do
      covariance = covariance/real(n,dp)
   end function hac_covariance_of_mean

   pure function outer_product(a, b) result(c)
      real(dp), intent(in) :: a(:), b(:)
      real(dp) :: c(size(a),size(b))
      c = spread(a,2,size(b))*spread(b,1,size(a))
   end function outer_product

   pure function normal_product_covariance(i, j, k, l, mu, sigma, &
                                           fit_intercept) result(value)
      integer, intent(in) :: i, j, k, l
      real(dp), intent(in) :: mu(:), sigma(:, :)
      logical, intent(in) :: fit_intercept
      real(dp) :: value
      integer :: first_vars(2), second_vars(2), n_first, n_second
      integer :: a, b, c, d

      call product_variables(i,j,fit_intercept,first_vars,n_first)
      call product_variables(k,l,fit_intercept,second_vars,n_second)
      value = 0.0_dp
      if (n_first == 0 .or. n_second == 0) return

      if (n_first == 1 .and. n_second == 1) then
         value = sigma(first_vars(1),second_vars(1))
      else if (n_first == 1 .and. n_second == 2) then
         a = first_vars(1)
         c = second_vars(1)
         d = second_vars(2)
         value = mu(c)*sigma(a,d)+mu(d)*sigma(a,c)
      else if (n_first == 2 .and. n_second == 1) then
         a = second_vars(1)
         c = first_vars(1)
         d = first_vars(2)
         value = mu(c)*sigma(a,d)+mu(d)*sigma(a,c)
      else
         a = first_vars(1)
         b = first_vars(2)
         c = second_vars(1)
         d = second_vars(2)
         value = mu(a)*mu(c)*sigma(b,d)+mu(a)*mu(d)*sigma(b,c) &
            +mu(b)*mu(c)*sigma(a,d)+mu(b)*mu(d)*sigma(a,c) &
            +sigma(a,c)*sigma(b,d)+sigma(a,d)*sigma(b,c)
      end if
   end function normal_product_covariance

   pure subroutine product_variables(i, j, fit_intercept, variables, count)
      integer, intent(in) :: i, j
      logical, intent(in) :: fit_intercept
      integer, intent(out) :: variables(2), count
      integer :: mapped

      variables = 0
      count = 0
      if (.not. fit_intercept .or. i /= 1) then
         mapped = i-merge(1,0,fit_intercept)
         count = count+1
         variables(count) = mapped
      end if
      if (.not. fit_intercept .or. j /= 1) then
         mapped = j-merge(1,0,fit_intercept)
         count = count+1
         variables(count) = mapped
      end if
   end subroutine product_variables

   subroutine remove_incomplete_rows(x, clean_x)
      real(dp), intent(in) :: x(:, :)
      real(dp), allocatable, intent(out) :: clean_x(:, :)
      logical, allocatable :: complete(:)
      integer :: i, n_complete, destination

      allocate(complete(size(x,1)))
      do i = 1, size(x,1)
         complete(i) = all(ieee_is_finite(x(i,:)))
      end do
      n_complete = count(complete)
      allocate(clean_x(n_complete,size(x,2)))
      destination = 0
      do i = 1, size(x,1)
         if (.not. complete(i)) cycle
         destination = destination+1
         clean_x(destination,:) = x(i,:)
      end do
   end subroutine remove_incomplete_rows

end module markowitzr_moments
