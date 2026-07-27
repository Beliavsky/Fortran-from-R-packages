! SPDX-License-Identifier: GPL-3.0-only
!
! Experimental modern Fortran translation of computational methods from the
! R package rmgarch, copyright (C) 2008-2025 Alexios Galanos.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 3 only.
module rmgarch_diagnostics
   use rmgarch_kinds, only : dp
   use rmgarch_math, only : chi_square_cdf, normal_cdf, covariance_matrix, inverse_spd, &
      symmetric_eigen_jacobi, quadratic_form_spd
   implicit none
   private

   public :: correlation_distance, correlation_distance_matrix
   public :: dcc_constancy_test, mardia_test

contains

   function correlation_distance(c1, c2, method) result(value)
      real(dp), intent(in) :: c1(:,:), c2(:,:)
      character(len=*), intent(in), optional :: method
      real(dp) :: value
      real(dp), allocatable :: differences(:), vals1(:), vals2(:), vecs(:,:)
      character(len=16) :: selected
      integer :: i, j, k, m
      logical :: ok
      selected = 'ma'; if (present(method)) selected = adjustl(method)
      m = size(c1,1); allocate(differences(m*(m-1)/2)); k = 0
      do j = 2, m
         do i = 1, j-1
            k = k+1; differences(k) = c1(i,j)-c2(i,j)
         end do
      end do
      select case (trim(selected))
      case ('ma')
         value = sum(abs(differences))/real(max(k,1),dp)
      case ('ms')
         value = sum(differences*differences)/real(max(k,1),dp)
      case ('meda')
         call sort_values(differences); value = median_sorted(differences)
      case ('meds')
         differences = differences*differences; call sort_values(differences); value = median_sorted(differences)
      case ('eigen')
         allocate(vals1(m),vals2(m),vecs(m,m))
         call symmetric_eigen_jacobi(c1,vals1,vecs,ok)
         call symmetric_eigen_jacobi(c2,vals2,vecs,ok)
         value = vals1(1)-vals2(1)
      case default
         value = sum(abs(differences))/real(max(k,1),dp)
      end select
   end function correlation_distance

   function correlation_distance_matrix(r, stride, method) result(distance)
      real(dp), intent(in) :: r(:,:,:)
      integer, intent(in), optional :: stride
      character(len=*), intent(in), optional :: method
      real(dp), allocatable :: distance(:,:)
      integer, allocatable :: indices(:)
      integer :: step, count, t, i, j
      step = 25; if (present(stride)) step = max(1,stride)
      count = (size(r,3)-1)/step+1
      if (1+(count-1)*step /= size(r,3)) count = count+1
      allocate(indices(count)); t = 1
      do i = 1, count-1
         indices(i) = t; t = t+step
      end do
      indices(count) = size(r,3)
      allocate(distance(count,count)); distance = 0.0_dp
      do j = 2, count
         do i = 1, j-1
            distance(i,j) = correlation_distance(r(:,:,indices(i)),r(:,:,indices(j)),method)
            distance(j,i) = distance(i,j)
         end do
      end do
   end function correlation_distance_matrix

   subroutine dcc_constancy_test(z, lags, statistic, p_value)
      real(dp), intent(in) :: z(:,:)
      integer, intent(in) :: lags
      real(dp), intent(out) :: statistic, p_value
      real(dp), allocatable :: op(:,:), regressors(:,:), regressand(:), xtx(:,:), inv(:,:), beta(:), residual(:)
      integer :: n, m, pairs, i, j, k, t, row, cols
      logical :: ok
      real(dp) :: sigma
      n = size(z,1); m = size(z,2); pairs = m*(m-1)/2
      if (lags < 0 .or. n <= lags+1 .or. pairs == 0) then
         statistic = 0.0_dp; p_value = 1.0_dp; return
      end if
      allocate(op(n,pairs)); k = 0
      do j = 2, m
         do i = 1, j-1
            k = k+1; op(:,k) = z(:,i)*z(:,j)
         end do
      end do
      cols = lags+1
      allocate(regressors((n-lags)*pairs,cols),regressand((n-lags)*pairs))
      row = 0
      do k = 1, pairs
         do t = lags+1, n
            row = row+1
            regressors(row,1) = 1.0_dp
            do i = 1, lags
               regressors(row,i+1) = op(t-i,k)
            end do
            regressand(row) = op(t,k)
         end do
      end do
      allocate(xtx(cols,cols),inv(cols,cols),beta(cols),residual(row))
      xtx = matmul(transpose(regressors),regressors)
      inv = inverse_spd(xtx,ok)
      if (.not. ok) then
         statistic = 0.0_dp; p_value = 1.0_dp; return
      end if
      beta = matmul(inv,matmul(transpose(regressors),regressand))
      residual = regressand-matmul(regressors,beta)
      sigma = dot_product(residual,residual)/real(row,dp)
      statistic = dot_product(beta,matmul(xtx,beta))/sqrt(max(sigma,1.0e-14_dp))
      p_value = 1.0_dp-chi_square_cdf(statistic,real(cols,dp))
   end subroutine dcc_constancy_test

   subroutine mardia_test(x, skew_statistic, skew_p_value, kurtosis_statistic, kurtosis_p_value)
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(out) :: skew_statistic, skew_p_value, kurtosis_statistic, kurtosis_p_value
      real(dp), allocatable :: centered(:,:), cov(:,:), inv(:,:), d(:,:)
      real(dp) :: b1p, b2p, correction, df
      integer :: n, m, i
      logical :: ok
      n = size(x,1); m = size(x,2)
      allocate(centered(n,m),cov(m,m),inv(m,m),d(n,n))
      centered = x-spread(sum(x,dim=1)/real(n,dp),1,n)
      cov = covariance_matrix(x); inv = inverse_spd(cov,ok)
      if (.not. ok) then
         skew_statistic = 0.0_dp; skew_p_value = 1.0_dp
         kurtosis_statistic = 0.0_dp; kurtosis_p_value = 1.0_dp; return
      end if
      d = matmul(centered,matmul(inv,transpose(centered)))
      b1p = sum(d**3)/real(n*n,dp)
      b2p = 0.0_dp
      do i = 1, n
         b2p = b2p+d(i,i)**2
      end do
      b2p = b2p/real(n,dp)
      correction = real((m+1)*(n+1)*(n+3),dp)/real(n*((n+1)*(m+1)-6),dp)
      df = real(m*(m+1)*(m+2),dp)/6.0_dp
      skew_statistic = real(n,dp)*b1p*correction/6.0_dp
      skew_p_value = 1.0_dp-chi_square_cdf(skew_statistic,df)
      kurtosis_statistic = (b2p-real(m*(m+2),dp))/sqrt(8.0_dp*real(m*(m+2),dp)/real(n,dp))
      kurtosis_p_value = 2.0_dp*(1.0_dp-normal_cdf(abs(kurtosis_statistic)))
   end subroutine mardia_test

   subroutine sort_values(x)
      real(dp), intent(inout) :: x(:)
      integer :: i, j
      real(dp) :: key
      do i = 2, size(x)
         key = x(i); j = i-1
         do while (j >= 1)
            if (x(j) <= key) exit
            x(j+1) = x(j); j = j-1
         end do
         x(j+1) = key
      end do
   end subroutine sort_values

   pure function median_sorted(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value
      integer :: n
      n = size(x)
      if (n == 0) then
         value = 0.0_dp
      else if (mod(n,2) == 1) then
         value = x((n+1)/2)
      else
         value = 0.5_dp*(x(n/2)+x(n/2+1))
      end if
   end function median_sorted

end module rmgarch_diagnostics
