! SPDX-License-Identifier: GPL-2.0-only
module ppcor
   use ppcor_kinds, only : dp
   use ppcor_linalg, only : symmetric_inverse_or_pinv
   use ppcor_special, only : normal_cdf, student_t_cdf
   use ppcor_stats, only : association_matrix, method_name, method_from_name, &
                           ppcor_pearson, ppcor_kendall, ppcor_spearman
   implicit none
   private

   integer, parameter, public :: ppcor_success = 0
   integer, parameter, public :: ppcor_invalid_shape = 1
   integer, parameter, public :: ppcor_nonfinite_data = 2
   integer, parameter, public :: ppcor_invalid_method = 3
   integer, parameter, public :: ppcor_constant_variable = 4
   integer, parameter, public :: ppcor_linalg_failure = 5
   integer, parameter, public :: ppcor_insufficient_df = 6
   public :: dp, ppcor_pearson, ppcor_kendall, ppcor_spearman
   public :: ppcor_result, ppcor_test_result
   public :: partial_correlation, semi_partial_correlation
   public :: pcor, spcor, pcor_test, spcor_test
   public :: method_name, method_from_name

   type :: ppcor_result
      real(dp), allocatable :: estimate(:,:)
      real(dp), allocatable :: p_value(:,:)
      real(dp), allocatable :: statistic(:,:)
      integer :: n = 0
      integer :: gp = 0
      integer :: method = ppcor_pearson
      integer :: rank = 0
      logical :: used_pseudoinverse = .false.
      integer :: status = ppcor_success
      character(len=160) :: message = ''
   end type ppcor_result

   type :: ppcor_test_result
      real(dp) :: estimate = 0.0_dp
      real(dp) :: p_value = 1.0_dp
      real(dp) :: statistic = 0.0_dp
      integer :: n = 0
      integer :: gp = 0
      integer :: method = ppcor_pearson
      logical :: used_pseudoinverse = .false.
      integer :: status = ppcor_success
      character(len=160) :: message = ''
   end type ppcor_test_result

   interface pcor
      module procedure partial_correlation
   end interface pcor

   interface spcor
      module procedure semi_partial_correlation
   end interface spcor

   interface pcor_test
      module procedure pcor_test_vector_z
      module procedure pcor_test_matrix_z
   end interface pcor_test

   interface spcor_test
      module procedure spcor_test_vector_z
      module procedure spcor_test_matrix_z
   end interface spcor_test

contains

   subroutine partial_correlation(x, result, method)
      real(dp), intent(in) :: x(:,:)
      type(ppcor_result), intent(out) :: result
      integer, intent(in), optional :: method
      call correlation_core(x, result, method, .false.)
   end subroutine partial_correlation

   subroutine semi_partial_correlation(x, result, method)
      real(dp), intent(in) :: x(:,:)
      type(ppcor_result), intent(out) :: result
      integer, intent(in), optional :: method
      call correlation_core(x, result, method, .true.)
   end subroutine semi_partial_correlation

   subroutine correlation_core(x, result, method, semi)
      real(dp), intent(in) :: x(:,:)
      type(ppcor_result), intent(out) :: result
      integer, intent(in), optional :: method
      logical, intent(in) :: semi
      real(dp), allocatable :: assoc(:,:), inv_assoc(:,:)
      real(dp) :: denom, df, r, variance
      integer :: i, j, info, m, p

      result = ppcor_result()
      m = ppcor_pearson
      if (present(method)) m = method
      result%method = m
      result%n = size(x,1)
      p = size(x,2)
      result%gp = p-2

      if (result%n < 2 .or. p < 2) then
         call set_error(result, ppcor_invalid_shape, 'x must have at least two rows and two columns')
         return
      end if
      if (m < ppcor_pearson .or. m > ppcor_spearman) then
         call set_error(result, ppcor_invalid_method, 'unknown correlation method')
         return
      end if

      call association_matrix(x, m, assoc, info)
      if (info /= 0) then
         select case (info)
         case (1)
            call set_error(result, ppcor_invalid_shape, 'invalid data shape')
         case (2)
            call set_error(result, ppcor_nonfinite_data, 'x contains a non-finite value')
         case (4)
            call set_error(result, ppcor_constant_variable, 'a variable has zero variance')
         case default
            call set_error(result, ppcor_invalid_method, 'unknown correlation method')
         end select
         return
      end if

      call symmetric_inverse_or_pinv(assoc, inv_assoc, result%used_pseudoinverse, &
                                     result%rank, info)
      if (info /= 0) then
         call set_error(result, ppcor_linalg_failure, 'association matrix eigensolver failed')
         return
      end if

      allocate(result%estimate(p,p), result%p_value(p,p), result%statistic(p,p))
      result%estimate = 0.0_dp
      result%p_value = 0.0_dp
      result%statistic = 0.0_dp
      do j = 1, p
         do i = 1, p
            if (i == j) then
               result%estimate(i,j) = 1.0_dp
            else
               denom = sqrt(max(tiny(1.0_dp), abs(inv_assoc(i,i)*inv_assoc(j,j))))
               r = -inv_assoc(i,j)/denom
               if (semi) then
                  if (abs(inv_assoc(j,j)) <= tiny(1.0_dp)) then
                     r = 0.0_dp
                  else
                     denom = sqrt(max(tiny(1.0_dp), abs(assoc(i,i)))) * &
                        sqrt(max(tiny(1.0_dp), abs(inv_assoc(i,i) - &
                                         inv_assoc(i,j)**2/inv_assoc(j,j))))
                     r = (-inv_assoc(i,j)/sqrt(max(tiny(1.0_dp), &
                          abs(inv_assoc(i,i)*inv_assoc(j,j)))))/denom
                  end if
               end if
               result%estimate(i,j) = max(-1.0_dp, min(1.0_dp, r))
            end if
         end do
      end do

      if (m == ppcor_kendall) then
         if (result%n-result%gp <= 1) then
            call set_error(result, ppcor_insufficient_df, 'insufficient observations for Kendall test')
            return
         end if
         variance = 2.0_dp*(2.0_dp*real(result%n-result%gp,dp)+5.0_dp) / &
                    (9.0_dp*real(result%n-result%gp,dp)* &
                     real(result%n-1-result%gp,dp))
         do j = 1, p
            do i = 1, p
               if (i == j) then
                  result%statistic(i,j) = 0.0_dp
                  result%p_value(i,j) = 0.0_dp
               else
                  result%statistic(i,j) = result%estimate(i,j)/sqrt(variance)
                  result%p_value(i,j) = 2.0_dp*normal_cdf(-abs(result%statistic(i,j)))
               end if
            end do
         end do
      else
         df = real(result%n-2-result%gp,dp)
         if (df <= 0.0_dp) then
            call set_error(result, ppcor_insufficient_df, 'insufficient residual degrees of freedom')
            return
         end if
         do j = 1, p
            do i = 1, p
               if (i == j) then
                  result%statistic(i,j) = 0.0_dp
                  result%p_value(i,j) = 0.0_dp
               else
                  r = result%estimate(i,j)
                  if (abs(r) >= 1.0_dp-16.0_dp*epsilon(1.0_dp)) then
                     result%statistic(i,j) = sign(huge(1.0_dp), r)
                     result%p_value(i,j) = 0.0_dp
                  else
                     result%statistic(i,j) = r*sqrt(df/(1.0_dp-r*r))
                     result%p_value(i,j) = 2.0_dp*student_t_cdf(-abs(result%statistic(i,j)), df)
                  end if
               end if
            end do
         end do
      end if

      if (result%used_pseudoinverse) then
         result%message = 'Moore-Penrose pseudoinverse used for a rank-deficient association matrix'
      else
         result%message = 'ok'
      end if
   end subroutine correlation_core

   subroutine pcor_test_vector_z(x, y, z, result, method)
      real(dp), intent(in) :: x(:), y(:), z(:)
      type(ppcor_test_result), intent(out) :: result
      integer, intent(in), optional :: method
      real(dp), allocatable :: zm(:,:)
      allocate(zm(size(z),1))
      zm(:,1) = z
      call pcor_test_matrix_z(x, y, zm, result, method)
   end subroutine pcor_test_vector_z

   subroutine pcor_test_matrix_z(x, y, z, result, method)
      real(dp), intent(in) :: x(:), y(:), z(:,:)
      type(ppcor_test_result), intent(out) :: result
      integer, intent(in), optional :: method
      call pair_test_core(x, y, z, result, method, .false.)
   end subroutine pcor_test_matrix_z

   subroutine spcor_test_vector_z(x, y, z, result, method)
      real(dp), intent(in) :: x(:), y(:), z(:)
      type(ppcor_test_result), intent(out) :: result
      integer, intent(in), optional :: method
      real(dp), allocatable :: zm(:,:)
      allocate(zm(size(z),1))
      zm(:,1) = z
      call spcor_test_matrix_z(x, y, zm, result, method)
   end subroutine spcor_test_vector_z

   subroutine spcor_test_matrix_z(x, y, z, result, method)
      real(dp), intent(in) :: x(:), y(:), z(:,:)
      type(ppcor_test_result), intent(out) :: result
      integer, intent(in), optional :: method
      call pair_test_core(x, y, z, result, method, .true.)
   end subroutine spcor_test_matrix_z

   subroutine pair_test_core(x, y, z, result, method, semi)
      real(dp), intent(in) :: x(:), y(:), z(:,:)
      type(ppcor_test_result), intent(out) :: result
      integer, intent(in), optional :: method
      logical, intent(in) :: semi
      type(ppcor_result) :: full
      real(dp), allocatable :: xyz(:,:)
      integer :: m, n

      result = ppcor_test_result()
      m = ppcor_pearson
      if (present(method)) m = method
      result%method = m
      n = size(x)
      if (size(y) /= n .or. size(z,1) /= n) then
         result%status = ppcor_invalid_shape
         result%message = 'x, y, and z must have the same number of observations'
         return
      end if
      allocate(xyz(n,2+size(z,2)))
      xyz(:,1) = x
      xyz(:,2) = y
      if (size(z,2) > 0) xyz(:,3:) = z
      if (semi) then
         call semi_partial_correlation(xyz, full, m)
      else
         call partial_correlation(xyz, full, m)
      end if
      result%n = full%n
      result%gp = full%gp
      result%used_pseudoinverse = full%used_pseudoinverse
      result%status = full%status
      result%message = full%message
      if (full%status == ppcor_success) then
         result%estimate = full%estimate(1,2)
         result%p_value = full%p_value(1,2)
         result%statistic = full%statistic(1,2)
      end if
   end subroutine pair_test_core

   subroutine set_error(result, status, message)
      type(ppcor_result), intent(inout) :: result
      integer, intent(in) :: status
      character(len=*), intent(in) :: message
      result%status = status
      result%message = message
   end subroutine set_error

end module ppcor
