! SPDX-License-Identifier: Apache-2.0
module clarabel_psd
   use clarabel_kinds, only : dp
   implicit none
   private
   public :: psd_svec_upper, psd_smat_upper, psd_matrix_order

contains

   function psd_svec_upper(x) result(v)
      real(dp), intent(in) :: x(:, :)
      real(dp), allocatable :: v(:)
      integer :: n, i, j, k
      real(dp) :: scale

      n = size(x, 1)
      if (size(x, 2) /= n) error stop "psd_svec_upper: matrix must be square"
      scale = max(1.0_dp, maxval(abs(x)))
      if (maxval(abs(x - transpose(x))) > 100.0_dp * epsilon(1.0_dp) * scale) &
         error stop "psd_svec_upper: matrix is not symmetric"
      allocate(v(n * (n + 1) / 2))
      k = 0
      do j = 1, n
         do i = 1, j
            k = k + 1
            if (i == j) then
               v(k) = x(i, j)
            else
               v(k) = sqrt(2.0_dp) * x(i, j)
            end if
         end do
      end do
   end function psd_svec_upper

   function psd_smat_upper(v) result(x)
      real(dp), intent(in) :: v(:)
      real(dp), allocatable :: x(:, :)
      integer :: n, i, j, k

      n = psd_matrix_order(size(v))
      if (n < 0) error stop "psd_smat_upper: vector length is not triangular"
      allocate(x(n, n), source=0.0_dp)
      k = 0
      do j = 1, n
         do i = 1, j
            k = k + 1
            if (i == j) then
               x(i, j) = v(k)
            else
               x(i, j) = v(k) / sqrt(2.0_dp)
               x(j, i) = x(i, j)
            end if
         end do
      end do
   end function psd_smat_upper

   pure integer function psd_matrix_order(vector_length) result(n)
      integer, intent(in) :: vector_length
      integer :: candidate
      if (vector_length < 0) then
         n = -1
         return
      end if
      candidate = int((sqrt(1.0_dp + 8.0_dp * real(vector_length, dp)) - 1.0_dp) / 2.0_dp)
      if (candidate * (candidate + 1) / 2 == vector_length) then
         n = candidate
      else
         n = -1
      end if
   end function psd_matrix_order

end module clarabel_psd
