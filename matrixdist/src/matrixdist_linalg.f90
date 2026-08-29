! SPDX-License-Identifier: GPL-3.0-only
module matrixdist_linalg
   use r_compat, only: dp
   use r_linalg, only: shared_inverse_matrix => inverse_matrix
   use r_linalg, only: shared_solve_system => solve_system
   implicit none
   private
   public :: inf_norm, matrix_vanloan, max_diagonal, matrix_exponential
   public :: matrix_power, matrix_inverse, solve_matrix, solve_vector
   public :: eye_matrix, kronecker, kronecker_sum, block_diag2

contains

   function eye_matrix(n) result(a)
      integer, intent(in) :: n
      real(dp) :: a(n,n)
      integer :: i
      a = 0.0_dp
      do i = 1, n
         a(i,i) = 1.0_dp
      end do
   end function eye_matrix

   function inf_norm(a) result(value)
      real(dp), intent(in) :: a(:,:)
      real(dp) :: value
      integer :: i
      value = 0.0_dp
      do i = 1, size(a,1)
         value = max(value, sum(abs(a(i,:))))
      end do
   end function inf_norm

   function max_diagonal(a) result(value)
      real(dp), intent(in) :: a(:,:)
      real(dp) :: value
      integer :: i, n
      n = min(size(a,1), size(a,2))
      if (n == 0) then
         value = -huge(1.0_dp)
         return
      end if
      value = a(1,1)
      do i = 2, n
         value = max(value, a(i,i))
      end do
   end function max_diagonal

   function matrix_vanloan(a1, a2, b1) result(v)
      real(dp), intent(in) :: a1(:,:), a2(:,:), b1(:,:)
      real(dp), allocatable :: v(:,:)
      integer :: p1, p2
      p1 = size(a1,1)
      p2 = size(a2,1)
      allocate(v(p1+p2,p1+p2))
      v = 0.0_dp
      v(1:p1,1:p1) = a1
      v(1:p1,p1+1:p1+p2) = b1
      v(p1+1:p1+p2,p1+1:p1+p2) = a2
   end function matrix_vanloan

   function solve_matrix(a, b, info_out) result(x)
      real(dp), intent(in) :: a(:,:), b(:,:)
      integer, intent(out), optional :: info_out
      real(dp), allocatable :: x(:,:)
      integer :: n, nrhs, info
      n = size(a,1)
      nrhs = size(b,2)
      if (size(a,2) /= n .or. size(b,1) /= n) error stop "solve_matrix: incompatible dimensions"
      allocate(x(n,nrhs))
      call shared_solve_system(a, b, x, info)
      if (present(info_out)) info_out = info
      if (info /= 0 .and. .not. present(info_out)) error stop "solve_matrix: singular/invalid matrix"
   end function solve_matrix

   function solve_vector(a, b, info_out) result(x)
      real(dp), intent(in) :: a(:,:), b(:)
      integer, intent(out), optional :: info_out
      real(dp), allocatable :: x(:)
      real(dp), allocatable :: bb(:,:), xx(:,:)
      integer :: info
      allocate(bb(size(b),1))
      bb(:,1) = b
      xx = solve_matrix(a, bb, info)
      allocate(x(size(b)))
      x = xx(:,1)
      if (present(info_out)) info_out = info
   end function solve_vector

   function matrix_inverse(a, info_out) result(ai)
      real(dp), intent(in) :: a(:,:)
      integer, intent(out), optional :: info_out
      real(dp), allocatable :: ai(:,:)
      integer :: info
      if (size(a,2) /= size(a,1)) error stop "matrix_inverse: matrix must be square"
      call shared_inverse_matrix(a, ai, info)
      if (present(info_out)) info_out = info
   end function matrix_inverse

   function matrix_exponential(a) result(expm)
      ! Scaling/squaring with [6/6] Pade approximation, matching matrixdist's
      ! MATLAB-derived native implementation but with a zero-norm guard.
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable :: expm(:,:)
      real(dp), allocatable :: a2(:,:), x(:,:), d(:,:)
      real(dp) :: a_norm, c, t
      integer, parameter :: q = 6
      integer :: ee, s, k, parity, n, info
      n = size(a,1)
      if (size(a,2) /= n) error stop "matrix_exponential: matrix must be square"
      if (n == 0) then
         allocate(expm(0,0))
         return
      end if
      a_norm = inf_norm(a)
      if (a_norm == 0.0_dp) then
         expm = eye_matrix(n)
         return
      end if
      ee = int(log(a_norm)/log(2.0_dp)) + 1
      s = max(0, ee + 1)
      t = 2.0_dp**(-s)
      allocate(a2(n,n), x(n,n), d(n,n))
      a2 = a*t
      x = a2
      c = 0.5_dp
      expm = eye_matrix(n) + c*a2
      d = eye_matrix(n) - c*a2
      parity = 1
      do k = 2, q
         c = c*real(q-k+1,dp)/real(k*(2*q-k+1),dp)
         x = matmul(a2,x)
         expm = expm + c*x
         if (parity == 1) then
            d = d + c*x
         else
            d = d - c*x
         end if
         parity = 1 - parity
      end do
      expm = solve_matrix(d, expm, info)
      if (info /= 0) error stop "matrix_exponential: Pade solve failed"
      do k = 1, s
         expm = matmul(expm, expm)
      end do
   end function matrix_exponential

   function matrix_power(n, a) result(p)
      integer, intent(in) :: n
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable :: p(:,:), base(:,:)
      integer :: m, nn
      if (n < 0) error stop "matrix_power: n must be nonnegative"
      nn = size(a,1)
      if (size(a,2) /= nn) error stop "matrix_power: matrix must be square"
      p = eye_matrix(nn)
      if (n == 0) return
      base = a
      m = n
      do while (m > 0)
         if (mod(m,2) == 1) p = matmul(p,base)
         m = m/2
         if (m > 0) base = matmul(base,base)
      end do
   end function matrix_power

   function kronecker(a, b) result(k)
      real(dp), intent(in) :: a(:,:), b(:,:)
      real(dp), allocatable :: k(:,:)
      integer :: i, j, m, n, p, q
      m = size(a,1)
      n = size(a,2)
      p = size(b,1)
      q = size(b,2)
      allocate(k(m*p,n*q))
      do i = 1, m
         do j = 1, n
            k((i-1)*p+1:i*p,(j-1)*q+1:j*q) = a(i,j)*b
         end do
      end do
   end function kronecker

   function kronecker_sum(a, b) result(k)
      real(dp), intent(in) :: a(:,:), b(:,:)
      real(dp), allocatable :: k(:,:)
      if (size(a,1) /= size(a,2) .or. size(b,1) /= size(b,2)) &
         error stop "kronecker_sum: square matrices required"
      k = kronecker(a, eye_matrix(size(b,1))) + &
          kronecker(eye_matrix(size(a,1)), b)
   end function kronecker_sum

   function block_diag2(a, b) result(c)
      real(dp), intent(in) :: a(:,:), b(:,:)
      real(dp), allocatable :: c(:,:)
      integer :: m, n, p, q
      m=size(a,1)
      n=size(a,2)
      p=size(b,1)
      q=size(b,2)
      allocate(c(m+p,n+q))
      c=0.0_dp
      c(1:m,1:n)=a
      c(m+1:m+p,n+1:n+q)=b
   end function block_diag2

end module matrixdist_linalg
