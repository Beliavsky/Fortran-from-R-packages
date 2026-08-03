! SPDX-License-Identifier: GPL-2.0-or-later
module moments_statistics
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
   use moments_kinds, only : dp
   implicit none
   private

   public :: moment, all_moments, skewness, kurtosis, geary

   interface moment
      module procedure moment_vector
      module procedure moment_matrix
   end interface moment

   interface all_moments
      module procedure all_moments_vector
      module procedure all_moments_matrix
   end interface all_moments

   interface skewness
      module procedure skewness_vector
      module procedure skewness_matrix
   end interface skewness

   interface kurtosis
      module procedure kurtosis_vector
      module procedure kurtosis_matrix
   end interface kurtosis

   interface geary
      module procedure geary_vector
      module procedure geary_matrix
   end interface geary

contains

   real(dp) function quiet_nan() result(x)
      x = ieee_value(0.0_dp, ieee_quiet_nan)
   end function quiet_nan

   real(dp) function moment_vector(x, order, central, absolute, na_rm) result(value)
      real(dp), intent(in) :: x(:)
      integer, intent(in), optional :: order
      logical, intent(in), optional :: central, absolute, na_rm
      integer :: k, i, n_used
      logical :: use_central, use_absolute, remove_na
      real(dp) :: center, term, total

      k = 1
      if (present(order)) k = order
      use_central = .false.
      if (present(central)) use_central = central
      use_absolute = .false.
      if (present(absolute)) use_absolute = absolute
      remove_na = .false.
      if (present(na_rm)) remove_na = na_rm

      if (k < 0 .or. size(x) == 0) then
         value = quiet_nan()
         return
      end if
      if (k == 0) then
         value = 1.0_dp
         return
      end if

      if (.not. remove_na) then
         do i = 1, size(x)
            if (.not. ieee_is_finite(x(i))) then
               value = quiet_nan()
               return
            end if
         end do
      end if

      n_used = 0
      center = 0.0_dp
      if (use_central) then
         do i = 1, size(x)
            if (ieee_is_finite(x(i))) then
               center = center + x(i)
               n_used = n_used + 1
            end if
         end do
         if (n_used == 0) then
            value = quiet_nan()
            return
         end if
         center = center / real(n_used, dp)
      end if

      total = 0.0_dp
      n_used = 0
      do i = 1, size(x)
         if (.not. ieee_is_finite(x(i))) cycle
         term = x(i) - center
         if (use_absolute) term = abs(term)
         total = total + term**k
         n_used = n_used + 1
      end do
      if (n_used == 0) then
         value = quiet_nan()
      else
         value = total / real(n_used, dp)
      end if
   end function moment_vector

   function moment_matrix(x, order, central, absolute, na_rm) result(values)
      real(dp), intent(in) :: x(:, :)
      integer, intent(in), optional :: order
      logical, intent(in), optional :: central, absolute, na_rm
      real(dp), allocatable :: values(:)
      integer :: j

      allocate(values(size(x, 2)))
      do j = 1, size(x, 2)
         values(j) = moment_vector(x(:, j), order, central, absolute, na_rm)
      end do
   end function moment_matrix

   function all_moments_vector(x, order_max, central, absolute, na_rm) result(values)
      real(dp), intent(in) :: x(:)
      integer, intent(in), optional :: order_max
      logical, intent(in), optional :: central, absolute, na_rm
      real(dp), allocatable :: values(:)
      integer :: k, maximum

      maximum = 2
      if (present(order_max)) maximum = order_max
      if (maximum < 0) then
         allocate(values(0))
         return
      end if
      allocate(values(maximum + 1))
      do k = 0, maximum
         values(k + 1) = moment_vector(x, k, central, absolute, na_rm)
      end do
   end function all_moments_vector

   function all_moments_matrix(x, order_max, central, absolute, na_rm) result(values)
      real(dp), intent(in) :: x(:, :)
      integer, intent(in), optional :: order_max
      logical, intent(in), optional :: central, absolute, na_rm
      real(dp), allocatable :: values(:, :)
      integer :: k, j, maximum

      maximum = 2
      if (present(order_max)) maximum = order_max
      if (maximum < 0) then
         allocate(values(0, size(x, 2)))
         return
      end if
      allocate(values(maximum + 1, size(x, 2)))
      do j = 1, size(x, 2)
         do k = 0, maximum
            values(k + 1, j) = moment_vector(x(:, j), k, central, absolute, na_rm)
         end do
      end do
   end function all_moments_matrix

   real(dp) function skewness_vector(x, na_rm) result(value)
      real(dp), intent(in) :: x(:)
      logical, intent(in), optional :: na_rm
      real(dp) :: m2, m3

      m2 = moment_vector(x, 2, central=.true., na_rm=na_rm)
      m3 = moment_vector(x, 3, central=.true., na_rm=na_rm)
      if (.not. ieee_is_finite(m2) .or. m2 <= 0.0_dp) then
         value = quiet_nan()
      else
         value = m3 / m2**1.5_dp
      end if
   end function skewness_vector

   function skewness_matrix(x, na_rm) result(values)
      real(dp), intent(in) :: x(:, :)
      logical, intent(in), optional :: na_rm
      real(dp), allocatable :: values(:)
      integer :: j

      allocate(values(size(x, 2)))
      do j = 1, size(x, 2)
         values(j) = skewness_vector(x(:, j), na_rm)
      end do
   end function skewness_matrix

   real(dp) function kurtosis_vector(x, na_rm) result(value)
      real(dp), intent(in) :: x(:)
      logical, intent(in), optional :: na_rm
      real(dp) :: m2, m4

      m2 = moment_vector(x, 2, central=.true., na_rm=na_rm)
      m4 = moment_vector(x, 4, central=.true., na_rm=na_rm)
      if (.not. ieee_is_finite(m2) .or. m2 <= 0.0_dp) then
         value = quiet_nan()
      else
         value = m4 / (m2 * m2)
      end if
   end function kurtosis_vector

   function kurtosis_matrix(x, na_rm) result(values)
      real(dp), intent(in) :: x(:, :)
      logical, intent(in), optional :: na_rm
      real(dp), allocatable :: values(:)
      integer :: j

      allocate(values(size(x, 2)))
      do j = 1, size(x, 2)
         values(j) = kurtosis_vector(x(:, j), na_rm)
      end do
   end function kurtosis_matrix

   real(dp) function geary_vector(x, na_rm) result(value)
      real(dp), intent(in) :: x(:)
      logical, intent(in), optional :: na_rm
      real(dp) :: rho, tau

      rho = sqrt(moment_vector(x, 2, central=.true., na_rm=na_rm))
      tau = moment_vector(x, 1, central=.true., absolute=.true., na_rm=na_rm)
      if (.not. ieee_is_finite(rho) .or. rho <= 0.0_dp) then
         value = quiet_nan()
      else
         value = tau / rho
      end if
   end function geary_vector

   function geary_matrix(x, na_rm) result(values)
      real(dp), intent(in) :: x(:, :)
      logical, intent(in), optional :: na_rm
      real(dp), allocatable :: values(:)
      integer :: j

      allocate(values(size(x, 2)))
      do j = 1, size(x, 2)
         values(j) = geary_vector(x(:, j), na_rm)
      end do
   end function geary_matrix

end module moments_statistics
