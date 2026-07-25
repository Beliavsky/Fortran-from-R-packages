! SPDX-License-Identifier: GPL-2.0-or-later
!
! Computational translation of gogarch, copyright (C) 2008-2026 Bernhard Pfaff.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 2 or later.
module gogarch_core
   use gogarch_kinds, only : dp
   use gogarch_linalg, only : covariance_matrix, symmetric_eigen, symmetric_sqrt, symmetric_invsqrt
   use gogarch_linalg, only : covariance_to_correlation, outer_product
   implicit none
   private
   public :: initialize_gogarch, cora, build_covariance_path
   public :: conditional_variances, conditional_correlations
contains

   subroutine initialize_gogarch(data, covariance, eigenvectors, eigenvalues, covariance_sqrt, covariance_invsqrt, ok)
      real(dp), intent(in) :: data(:,:)
      real(dp), intent(out) :: covariance(size(data,2),size(data,2))
      real(dp), intent(out) :: eigenvectors(size(data,2),size(data,2)), eigenvalues(size(data,2))
      real(dp), intent(out) :: covariance_sqrt(size(data,2),size(data,2))
      real(dp), intent(out) :: covariance_invsqrt(size(data,2),size(data,2))
      logical, intent(out) :: ok
      logical :: ok1, ok2
      covariance = covariance_matrix(data,center=.false.)
      call symmetric_eigen(covariance,eigenvalues,eigenvectors,ok)
      if (.not. ok .or. minval(eigenvalues) <= 1.0e-14_dp) then
         covariance_sqrt = 0.0_dp
         covariance_invsqrt = 0.0_dp
         ok = .false.
         return
      end if
      covariance_sqrt = symmetric_sqrt(covariance,1.0e-14_dp,ok1)
      covariance_invsqrt = symmetric_invsqrt(covariance,1.0e-14_dp,ok2)
      ok = ok1 .and. ok2
   end subroutine initialize_gogarch

   function cora(ssi, lag, standardize, ok) result(gamma)
      real(dp), intent(in) :: ssi(:,:,:)
      integer, intent(in), optional :: lag
      logical, intent(in), optional :: standardize
      logical, intent(out), optional :: ok
      real(dp) :: gamma(size(ssi,1),size(ssi,2))
      real(dp) :: gamma0(size(ssi,1),size(ssi,2)), invroot(size(ssi,1),size(ssi,2))
      logical :: do_standardize, eig_ok
      integer :: n, lag_value, t, count
      n = size(ssi,3)
      lag_value = 1
      if (present(lag)) lag_value = abs(lag)
      do_standardize = .true.
      if (present(standardize)) do_standardize = standardize
      if (size(ssi,1) /= size(ssi,2) .or. lag_value >= n) then
         gamma = 0.0_dp
         if (present(ok)) ok = .false.
         return
      end if
      gamma0 = 0.0_dp
      do t = 1, n
         gamma0 = gamma0+matmul(ssi(:,:,t),ssi(:,:,t))
      end do
      gamma0 = gamma0/real(n,dp)
      gamma = 0.0_dp
      count = n-lag_value
      do t = 1, count
         gamma = gamma+matmul(ssi(:,:,t+lag_value),ssi(:,:,t))
      end do
      gamma = gamma/real(count,dp)
      if (do_standardize) then
         invroot = symmetric_invsqrt(gamma0,1.0e-12_dp,eig_ok)
         if (.not. eig_ok) then
            gamma = 0.0_dp
            if (present(ok)) ok = .false.
            return
         end if
         gamma = matmul(invroot,matmul(gamma,invroot))
      end if
      gamma = 0.5_dp*(gamma+transpose(gamma))
      if (present(ok)) ok = .true.
   end function cora

   subroutine build_covariance_path(mixing, factor_variance, covariance)
      real(dp), intent(in) :: mixing(:,:), factor_variance(:,:)
      real(dp), intent(out) :: covariance(size(mixing,1),size(mixing,1),size(factor_variance,1))
      real(dp) :: scaled(size(mixing,1),size(mixing,2))
      integer :: t, j
      do t = 1, size(factor_variance,1)
         scaled = mixing
         do j = 1, size(mixing,2)
            scaled(:,j) = scaled(:,j)*factor_variance(t,j)
         end do
         covariance(:,:,t) = matmul(scaled,transpose(mixing))
         covariance(:,:,t) = 0.5_dp*(covariance(:,:,t)+transpose(covariance(:,:,t)))
      end do
   end subroutine build_covariance_path

   subroutine conditional_variances(covariance, variances)
      real(dp), intent(in) :: covariance(:,:,:)
      real(dp), intent(out) :: variances(size(covariance,3),size(covariance,1))
      integer :: t, i
      do t = 1, size(covariance,3)
         do i = 1, size(covariance,1)
            variances(t,i) = covariance(i,i,t)
         end do
      end do
   end subroutine conditional_variances

   subroutine conditional_correlations(covariance, correlations)
      real(dp), intent(in) :: covariance(:,:,:)
      real(dp), intent(out) :: correlations(size(covariance,1),size(covariance,2),size(covariance,3))
      integer :: t
      do t = 1, size(covariance,3)
         correlations(:,:,t) = covariance_to_correlation(covariance(:,:,t))
      end do
   end subroutine conditional_correlations

end module gogarch_core
