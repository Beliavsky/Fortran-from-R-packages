! SPDX-License-Identifier: GPL-2.0-or-later
module fracdiff_polynomial
   use fracdiff_kinds, only : dp, pi_dp
   implicit none
   private

   public :: polynomial_multiply, polynomial_roots, minimum_ar_root_modulus

contains

   subroutine polynomial_multiply(a, b, product)
      real(dp), intent(in) :: a(:), b(:)
      real(dp), intent(out) :: product(:)
      integer :: i, j

      if (size(product) /= size(a) + size(b) - 1) then
         product = 0.0_dp
         return
      end if
      product = 0.0_dp
      do j = 1, size(b)
         do i = 1, size(a)
            product(i + j - 1) = product(i + j - 1) + a(i)*b(j)
         end do
      end do
   end subroutine polynomial_multiply

   subroutine polynomial_roots(coefficients, roots, status)
      real(dp), intent(in) :: coefficients(:)
      complex(dp), intent(out) :: roots(:)
      integer, intent(out) :: status

      complex(dp), allocatable :: old_roots(:)
      complex(dp) :: denominator, correction, value
      real(dp) :: radius, change
      integer :: degree, i, j, iteration

      degree = size(coefficients) - 1
      status = 0
      if (degree < 1 .or. size(roots) /= degree .or. abs(coefficients(degree + 1)) <= tiny(1.0_dp)) then
         status = 1
         return
      end if

      radius = 1.0_dp + maxval(abs(coefficients(1:degree)/coefficients(degree + 1)))
      do i = 1, degree
         roots(i) = radius*cmplx(cos(2.0_dp*pi_dp*real(i - 1, dp)/real(degree, dp)), &
            sin(2.0_dp*pi_dp*real(i - 1, dp)/real(degree, dp)), kind=dp)
      end do
      allocate(old_roots(degree))

      do iteration = 1, 500
         old_roots = roots
         change = 0.0_dp
         do i = 1, degree
            value = cmplx(coefficients(degree + 1), 0.0_dp, kind=dp)
            do j = degree, 1, -1
               value = value*old_roots(i) + coefficients(j)
            end do
            denominator = cmplx(coefficients(degree + 1), 0.0_dp, kind=dp)
            do j = 1, degree
               if (j /= i) denominator = denominator*(old_roots(i) - old_roots(j))
            end do
            if (abs(denominator) <= tiny(1.0_dp)) denominator = cmplx(tiny(1.0_dp), 0.0_dp, kind=dp)
            correction = value/denominator
            roots(i) = old_roots(i) - correction
            change = max(change, abs(correction))
         end do
         if (change <= 100.0_dp*epsilon(1.0_dp)*(1.0_dp + maxval(abs(roots)))) return
      end do
      status = 2
   end subroutine polynomial_roots

   function minimum_ar_root_modulus(ar, status) result(minimum_modulus)
      real(dp), intent(in) :: ar(:)
      integer, intent(out), optional :: status
      real(dp) :: minimum_modulus

      real(dp), allocatable :: coefficients(:)
      complex(dp), allocatable :: roots(:)
      integer :: degree, local_status

      degree = size(ar)
      if (degree == 0) then
         minimum_modulus = huge(1.0_dp)
         if (present(status)) status = 0
         return
      end if

      allocate(coefficients(degree + 1), roots(degree))
      coefficients(1) = 1.0_dp
      coefficients(2:) = -ar
      call polynomial_roots(coefficients, roots, local_status)
      if (local_status /= 0) then
         minimum_modulus = 0.0_dp
      else
         minimum_modulus = minval(abs(roots))
      end if
      if (present(status)) status = local_status
   end function minimum_ar_root_modulus

end module fracdiff_polynomial
