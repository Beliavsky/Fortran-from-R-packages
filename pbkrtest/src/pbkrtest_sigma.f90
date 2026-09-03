! SPDX-License-Identifier: GPL-2.0-or-later
module pbkrtest_sigma
   use r_kinds, only : dp
   use pbkrtest_types, only : pbkr_invalid_argument, pbkr_invalid_shape, pbkr_success, &
      random_sigma_term_t, sigma_g_result_t
   implicit none
   private
   public :: build_sigma_g

contains

   pure subroutine build_sigma_g(terms, residual_variance, result, status)
      type(random_sigma_term_t), intent(in) :: terms(:) !! Random-effect terms with expanded `Z'` blocks and covariance matrices.
      real(dp), intent(in) :: residual_variance !! Final residual variance component; must be nonnegative.
      type(sigma_g_result_t), intent(out) :: result !! Decomposition `Sigma = sum(gamma(r) * G(:,:,r))`.
      integer, intent(out) :: status !! `pbkr_success` on success or a package error code.
      integer :: c
      integer :: i
      integer :: idx
      integer :: j
      integer :: level
      integer :: m
      integer :: n
      integer :: offset
      integer :: r
      integer :: s
      integer :: t

      status = pbkr_success
      if (residual_variance < 0.0_dp) then
         status = pbkr_invalid_argument
         return
      end if
      if (size(terms) == 0) then
         status = pbkr_invalid_argument
         return
      end if
      if (.not. allocated(terms(1)%z)) then
         status = pbkr_invalid_shape
         return
      end if
      n = size(terms(1)%z, 2)
      m = 1
      do t = 1, size(terms)
         if (.not. allocated(terms(t)%z) .or. .not. allocated(terms(t)%covariance)) then
            status = pbkr_invalid_shape
            return
         end if
         c = size(terms(t)%covariance, 1)
         if (c <= 0 .or. size(terms(t)%covariance, 2) /= c) then
            status = pbkr_invalid_shape
            return
         end if
         if (terms(t)%n_levels <= 0) then
            status = pbkr_invalid_argument
            return
         end if
         if (size(terms(t)%z, 1) /= terms(t)%n_levels * c .or. size(terms(t)%z, 2) /= n) then
            status = pbkr_invalid_shape
            return
         end if
         m = m + c * (c + 1) / 2
      end do

      allocate(result%sigma(n, n), result%g(n, n, m), result%gamma(m))
      result%sigma = 0.0_dp
      result%g = 0.0_dp
      result%gamma = 0.0_dp
      idx = 0

      do t = 1, size(terms)
         c = size(terms(t)%covariance, 1)
         do i = 1, c
            do j = i, c
               idx = idx + 1
               result%gamma(idx) = terms(t)%covariance(j, i)
               do level = 1, terms(t)%n_levels
                  offset = (level - 1) * c
                  do r = 1, n
                     do s = 1, n
                        if (i == j) then
                           result%g(r, s, idx) = result%g(r, s, idx) + &
                              terms(t)%z(offset + i, r) * terms(t)%z(offset + i, s)
                        else
                           result%g(r, s, idx) = result%g(r, s, idx) + &
                              terms(t)%z(offset + i, r) * terms(t)%z(offset + j, s) + &
                              terms(t)%z(offset + j, r) * terms(t)%z(offset + i, s)
                        end if
                     end do
                  end do
               end do
            end do
         end do
      end do

      idx = idx + 1
      result%gamma(idx) = residual_variance
      do i = 1, n
         result%g(i, i, idx) = 1.0_dp
      end do
      do i = 1, m
         result%sigma = result%sigma + result%gamma(i) * result%g(:, :, i)
      end do
   end subroutine build_sigma_g

end module pbkrtest_sigma
