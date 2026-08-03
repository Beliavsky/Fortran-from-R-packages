! SPDX-License-Identifier: GPL-2.0-or-later
module moments_transforms
   use moments_kinds, only : dp
   implicit none
   private

   public :: raw2central, central2raw, all_cumulants

   interface raw2central
      module procedure raw2central_vector
      module procedure raw2central_matrix
   end interface raw2central

   interface central2raw
      module procedure central2raw_vector
      module procedure central2raw_matrix
   end interface central2raw

   interface all_cumulants
      module procedure all_cumulants_vector
      module procedure all_cumulants_matrix
   end interface all_cumulants

contains

   pure real(dp) function binomial(n, k) result(value)
      integer, intent(in) :: n, k
      integer :: i, kk

      if (k < 0 .or. k > n) then
         value = 0.0_dp
         return
      end if
      kk = min(k, n - k)
      value = 1.0_dp
      do i = 1, kk
         value = value * real(n - kk + i, dp) / real(i, dp)
      end do
   end function binomial

   function raw2central_vector(mu_raw) result(mu_central)
      real(dp), intent(in) :: mu_raw(:)
      real(dp), allocatable :: mu_central(:)
      integer :: k, i, maximum
      real(dp) :: eta

      maximum = size(mu_raw) - 1
      allocate(mu_central(size(mu_raw)))
      if (size(mu_raw) == 0) return
      eta = 0.0_dp
      if (maximum >= 1) eta = mu_raw(2)
      do k = 0, maximum
         mu_central(k + 1) = 0.0_dp
         do i = 0, k
            mu_central(k + 1) = mu_central(k + 1) + binomial(k, i) * &
               (-eta)**(k - i) * mu_raw(i + 1)
         end do
      end do
   end function raw2central_vector

   function raw2central_matrix(mu_raw) result(mu_central)
      real(dp), intent(in) :: mu_raw(:, :)
      real(dp), allocatable :: mu_central(:, :)
      integer :: j
      real(dp), allocatable :: column(:)

      allocate(mu_central(size(mu_raw, 1), size(mu_raw, 2)))
      do j = 1, size(mu_raw, 2)
         column = raw2central_vector(mu_raw(:, j))
         mu_central(:, j) = column
      end do
   end function raw2central_matrix

   function central2raw_vector(mu_central, eta) result(mu_raw)
      real(dp), intent(in) :: mu_central(:)
      real(dp), intent(in) :: eta
      real(dp), allocatable :: mu_raw(:)
      integer :: k, i, maximum

      maximum = size(mu_central) - 1
      allocate(mu_raw(size(mu_central)))
      do k = 0, maximum
         mu_raw(k + 1) = 0.0_dp
         do i = 0, k
            mu_raw(k + 1) = mu_raw(k + 1) + binomial(k, i) * &
               eta**(k - i) * mu_central(i + 1)
         end do
      end do
   end function central2raw_vector

   function central2raw_matrix(mu_central, eta) result(mu_raw)
      real(dp), intent(in) :: mu_central(:, :)
      real(dp), intent(in) :: eta(:)
      real(dp), allocatable :: mu_raw(:, :)
      real(dp), allocatable :: column(:)
      integer :: j

      if (size(eta) /= size(mu_central, 2)) then
         allocate(mu_raw(0, 0))
         return
      end if
      allocate(mu_raw(size(mu_central, 1), size(mu_central, 2)))
      do j = 1, size(mu_central, 2)
         column = central2raw_vector(mu_central(:, j), eta(j))
         mu_raw(:, j) = column
      end do
   end function central2raw_matrix

   function all_cumulants_vector(mu_raw, legacy) result(kappa)
      real(dp), intent(in) :: mu_raw(:)
      logical, intent(in), optional :: legacy
      real(dp), allocatable :: kappa(:)
      integer :: n, m, maximum
      logical :: reproduce_upstream
      real(dp), allocatable :: mu_central(:)

      maximum = size(mu_raw) - 1
      allocate(kappa(size(mu_raw)))
      kappa = 0.0_dp
      if (maximum < 1) return

      reproduce_upstream = .false.
      if (present(legacy)) reproduce_upstream = legacy
      if (reproduce_upstream) then
         mu_central = raw2central_vector(mu_raw)
         kappa(2) = mu_central(2)
         if (maximum >= 2) kappa(3) = mu_central(3)
      else
         kappa(2) = mu_raw(2)
      end if

      do n = 2, maximum
         if (reproduce_upstream .and. n == 2) cycle
         kappa(n + 1) = mu_raw(n + 1)
         do m = 1, n - 1
            kappa(n + 1) = kappa(n + 1) - binomial(n - 1, m - 1) * &
               kappa(m + 1) * mu_raw(n - m + 1)
         end do
      end do
   end function all_cumulants_vector

   function all_cumulants_matrix(mu_raw, legacy) result(kappa)
      real(dp), intent(in) :: mu_raw(:, :)
      logical, intent(in), optional :: legacy
      real(dp), allocatable :: kappa(:, :)
      real(dp), allocatable :: column(:)
      integer :: j

      allocate(kappa(size(mu_raw, 1), size(mu_raw, 2)))
      do j = 1, size(mu_raw, 2)
         column = all_cumulants_vector(mu_raw(:, j), legacy)
         kappa(:, j) = column
      end do
   end function all_cumulants_matrix

end module moments_transforms
