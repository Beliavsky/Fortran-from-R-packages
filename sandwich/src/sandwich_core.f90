! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module sandwich_core
   use sandwich_kinds, only : dp
   use sandwich_status, only : SANDWICH_SUCCESS, SANDWICH_INVALID_ARGUMENT, &
      SANDWICH_DIMENSION_MISMATCH
   use sandwich_linalg, only : inverse_matrix
   implicit none
   private

   public :: meat, sandwich_covariance, vcov_opg, bread_from_information

contains

   subroutine meat(scores, meat_matrix, status, adjust)
      real(dp), intent(in) :: scores(:, :)
      real(dp), allocatable, intent(out) :: meat_matrix(:, :)
      integer, intent(out), optional :: status
      logical, intent(in), optional :: adjust
      integer :: n, k
      logical :: use_adjust

      n = size(scores, 1)
      k = size(scores, 2)
      if (n <= 0 .or. k <= 0) then
         allocate(meat_matrix(0, 0))
         if (present(status)) status = SANDWICH_INVALID_ARGUMENT
         return
      end if

      use_adjust = .false.
      if (present(adjust)) use_adjust = adjust
      if (use_adjust .and. n <= k) then
         allocate(meat_matrix(0, 0))
         if (present(status)) status = SANDWICH_INVALID_ARGUMENT
         return
      end if

      meat_matrix = matmul(transpose(scores), scores) / real(n, dp)
      if (use_adjust) meat_matrix = real(n, dp) / real(n - k, dp) * meat_matrix
      meat_matrix = 0.5_dp * (meat_matrix + transpose(meat_matrix))
      if (present(status)) status = SANDWICH_SUCCESS
   end subroutine meat

   subroutine sandwich_covariance(bread, meat_matrix, nobs, covariance, status)
      real(dp), intent(in) :: bread(:, :), meat_matrix(:, :)
      integer, intent(in) :: nobs
      real(dp), allocatable, intent(out) :: covariance(:, :)
      integer, intent(out), optional :: status
      integer :: k

      k = size(bread, 1)
      if (nobs <= 0 .or. k <= 0 .or. size(bread, 2) /= k .or. &
          size(meat_matrix, 1) /= k .or. size(meat_matrix, 2) /= k) then
         allocate(covariance(0, 0))
         if (present(status)) status = SANDWICH_DIMENSION_MISMATCH
         return
      end if

      covariance = matmul(bread, matmul(meat_matrix, bread)) / real(nobs, dp)
      covariance = 0.5_dp * (covariance + transpose(covariance))
      if (present(status)) status = SANDWICH_SUCCESS
   end subroutine sandwich_covariance

   subroutine vcov_opg(scores, covariance, status, adjust)
      real(dp), intent(in) :: scores(:, :)
      real(dp), allocatable, intent(out) :: covariance(:, :)
      integer, intent(out), optional :: status
      logical, intent(in), optional :: adjust
      real(dp), allocatable :: information(:, :)
      integer :: n, k, info
      logical :: use_adjust

      n = size(scores, 1)
      k = size(scores, 2)
      use_adjust = .false.
      if (present(adjust)) use_adjust = adjust

      if (n <= 0 .or. k <= 0 .or. (use_adjust .and. n <= k)) then
         allocate(covariance(0, 0))
         if (present(status)) status = SANDWICH_INVALID_ARGUMENT
         return
      end if

      information = matmul(transpose(scores), scores)
      call inverse_matrix(information, covariance, info)
      if (info == SANDWICH_SUCCESS .and. use_adjust) then
         covariance = real(n, dp) / real(n - k, dp) * covariance
      end if
      if (present(status)) status = info
   end subroutine vcov_opg

   subroutine bread_from_information(information, nobs, bread, status)
      real(dp), intent(in) :: information(:, :)
      integer, intent(in) :: nobs
      real(dp), allocatable, intent(out) :: bread(:, :)
      integer, intent(out), optional :: status
      integer :: info

      if (nobs <= 0) then
         allocate(bread(0, 0))
         if (present(status)) status = SANDWICH_INVALID_ARGUMENT
         return
      end if
      call inverse_matrix(information, bread, info)
      if (info == SANDWICH_SUCCESS) bread = real(nobs, dp) * bread
      if (present(status)) status = info
   end subroutine bread_from_information

end module sandwich_core
