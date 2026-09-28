module matrixcalc_special
   use r_kinds, only : dp
   implicit none
   private

   public :: creation_matrix, fibonacci_matrix, frobenius_matrix
   public :: hankel_matrix, hilbert_matrix, pascal_matrix
   public :: stirling_matrix, symmetric_pascal_matrix
   public :: toeplitz_matrix, vandermonde_matrix

contains

   pure function creation_matrix(n) result(a)
      integer, intent(in) :: n !! Order of the creation matrix; upstream requires n >= 2.
      real(dp), allocatable :: a(:, :)
      integer :: i

      allocate(a(n, n))
      a = 0.0_dp
      do i = 1, n - 1
         a(i + 1, i) = real(i, dp)
      end do
   end function creation_matrix

   pure function fibonacci_matrix(n) result(fmat)
      integer, intent(in) :: n !! Upstream order parameter; result has shape (n+1) by (n+1).
      real(dp), allocatable :: fmat(:, :)
      real(dp), allocatable :: f(:)
      integer :: i
      integer :: j
      integer :: idx

      allocate(f(n + 2))
      f = 0.0_dp
      f(1) = 1.0_dp
      f(2) = 1.0_dp
      do i = 3, n + 2
         f(i) = f(i - 1) + f(i - 2)
      end do

      allocate(fmat(n + 1, n + 1))
      fmat = 0.0_dp
      do i = 0, n
         do j = 0, n
            if (i - j + 1 >= 0) then
               idx = i - j + 2
               fmat(i + 1, j + 1) = f(idx)
            end if
         end do
      end do
   end function fibonacci_matrix

   pure function frobenius_matrix(n) result(fmat)
      integer, intent(in) :: n !! Order of the Frobenius companion-style matrix.
      real(dp), allocatable :: fmat(:, :)
      integer :: i
      integer :: j
      integer :: sign_value

      allocate(fmat(n, n))
      fmat = 0.0_dp
      do i = 1, n - 1
         fmat(i + 1, i) = 1.0_dp
      end do
      do j = 0, n - 1
         if (mod(n - 1 - j, 2) == 0) then
            sign_value = 1
         else
            sign_value = -1
         end if
         fmat(j + 1, n) = real(sign_value, dp) * binomial_real(n, j)
      end do
   end function frobenius_matrix

   pure function hankel_matrix(n, x) result(hmat)
      integer, intent(in) :: n !! Order of the square Hankel-style matrix.
      real(dp), intent(in) :: x(:) !! Source vector; upstream requires at least n values.
      real(dp), allocatable :: hmat(:, :)
      integer :: i
      integer :: j
      integer :: m
      integer :: idx

      m = size(x)
      allocate(hmat(n, n))
      do i = 1, n
         do j = 1, n
            idx = modulo(i + j - 2, m) + 1
            hmat(i, j) = x(idx)
         end do
      end do
   end function hankel_matrix

   pure function hilbert_matrix(n) result(hmat)
      integer, intent(in) :: n !! Order of the Hilbert matrix.
      real(dp), allocatable :: hmat(:, :)
      integer :: i
      integer :: j

      allocate(hmat(n, n))
      do j = 1, n
         do i = 1, n
            hmat(i, j) = 1.0_dp / real(i + j - 1, dp)
         end do
      end do
   end function hilbert_matrix

   pure function pascal_matrix(n) result(pmat)
      integer, intent(in) :: n !! Order of the lower-triangular Pascal matrix.
      real(dp), allocatable :: pmat(:, :)
      integer :: i
      integer :: j

      allocate(pmat(n, n))
      pmat = 0.0_dp
      do i = 1, n
         do j = 1, i
            pmat(i, j) = binomial_real(i - 1, j - 1)
         end do
      end do
   end function pascal_matrix

   pure function stirling_matrix(n) result(smat)
      integer, intent(in) :: n !! Maximum Stirling-number order; result has shape (n+1) by (n+1).
      real(dp), allocatable :: smat(:, :)
      integer :: i
      integer :: j

      allocate(smat(n + 1, n + 1))
      smat = 0.0_dp
      smat(1, 1) = 1.0_dp
      do i = 1, n
         do j = 1, i
            smat(i + 1, j + 1) = smat(i, j) + real(j, dp) * smat(i, j + 1)
         end do
      end do
   end function stirling_matrix

   pure function symmetric_pascal_matrix(n) result(pmat)
      integer, intent(in) :: n !! Order of the symmetric Pascal matrix.
      real(dp), allocatable :: pmat(:, :)
      integer :: i
      integer :: j

      allocate(pmat(n, n))
      do j = 1, n
         do i = 1, n
            pmat(i, j) = binomial_real(i + j - 2, j - 1)
         end do
      end do
   end function symmetric_pascal_matrix

   pure function toeplitz_matrix(n, x) result(tmat)
      integer, intent(in) :: n !! Order of the square Toeplitz matrix.
      real(dp), intent(in) :: x(:) !! Vector of length 2*n-1 ordered from lowest to highest diagonal.
      real(dp), allocatable :: tmat(:, :)
      integer :: i
      integer :: j

      allocate(tmat(n, n))
      do j = 1, n
         do i = 1, n
            tmat(i, j) = x(i - j + n)
         end do
      end do
   end function toeplitz_matrix

   pure function vandermonde_matrix(alpha, n) result(vmat)
      real(dp), intent(in) :: alpha(:) !! Vector of Vandermonde base values.
      integer, intent(in) :: n !! Number of columns/powers, beginning with power zero.
      real(dp), allocatable :: vmat(:, :)
      integer :: j

      allocate(vmat(size(alpha), n))
      if (n <= 0) return
      vmat(:, 1) = 1.0_dp
      do j = 2, n
         vmat(:, j) = alpha**(j - 1)
      end do
   end function vandermonde_matrix

   pure real(dp) function binomial_real(n, k) result(value)
      integer, intent(in) :: n !! Nonnegative upper argument of the binomial coefficient.
      integer, intent(in) :: k !! Lower argument of the binomial coefficient.
      integer :: i
      integer :: kk

      if (k < 0 .or. k > n) then
         value = 0.0_dp
         return
      end if
      kk = min(k, n - k)
      value = 1.0_dp
      do i = 1, kk
         value = value * real(n - kk + i, dp) / real(i, dp)
      end do
   end function binomial_real

end module matrixcalc_special
