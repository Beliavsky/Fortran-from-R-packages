! rmutil computational translation
! Copyright (C) 1998-2001 J.K. Lindsey
! Copyright (C) 2026 OpenAI (modern Fortran translation)
! SPDX-License-Identifier: GPL-2.0-or-later
module rmutil_linalg
   use rmutil_kinds, only : dp
   implicit none
   private
   public :: matrix_exp, matrix_power_integer, lin_diff_eqn
contains
   function matrix_exp(a, t, max_terms) result(ea)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(in), optional :: t
      integer, intent(in), optional :: max_terms
      real(dp), allocatable :: ea(:,:)
      real(dp), allocatable :: x(:,:), term(:,:), ident(:,:)
      real(dp) :: tt, norm1
      integer :: n, i, k, s, nterm
      n = size(a,1)
      if (size(a,2) /= n) error stop "matrix_exp: square matrix required"
      tt = 1.0_dp
      if (present(t)) tt = t
      nterm = 80
      if (present(max_terms)) nterm = max_terms
      allocate(ea(n,n), x(n,n), term(n,n), ident(n,n))
      ident = 0.0_dp
      do i = 1, n
         ident(i,i) = 1.0_dp
      end do
      x = tt*a
      norm1 = maxval(sum(abs(x),dim=1))
      s = 0
      if (norm1 > 0.5_dp) s = max(0, ceiling(log(norm1/0.5_dp)/log(2.0_dp)))
      x = x/(2.0_dp**s)
      ea = ident
      term = ident
      do k = 1, nterm
         term = matmul(term,x)/real(k,dp)
         ea = ea + term
         if (maxval(abs(term)) <= 4.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(ea)))) exit
      end do
      do k = 1, s
         ea = matmul(ea,ea)
      end do
   end function matrix_exp

   function matrix_power_integer(a, p) result(ap)
      real(dp), intent(in) :: a(:,:)
      integer, intent(in) :: p
      real(dp), allocatable :: ap(:,:), base(:,:)
      integer :: n, i, e
      n = size(a,1)
      if (size(a,2) /= n) error stop "matrix_power_integer: square matrix required"
      if (p < 0) error stop "matrix_power_integer: negative powers are not implemented"
      allocate(ap(n,n),base(n,n))
      ap = 0.0_dp
      do i = 1, n
         ap(i,i) = 1.0_dp
      end do
      base = a
      e = p
      do while (e > 0)
         if (mod(e,2) == 1) ap = matmul(ap,base)
         e = e/2
         if (e > 0) base = matmul(base,base)
      end do
   end function matrix_power_integer

   function lin_diff_eqn(a, initial, times) result(y)
      real(dp), intent(in) :: a(:,:), initial(:), times(:)
      real(dp), allocatable :: y(:,:)
      real(dp), allocatable :: e(:,:)
      integer :: i
      if (size(a,1) /= size(a,2) .or. size(initial) /= size(a,1)) &
         error stop "lin_diff_eqn: dimension mismatch"
      allocate(y(size(times),size(initial)))
      do i = 1, size(times)
         e = matrix_exp(a,times(i))
         y(i,:) = matmul(e,initial)
      end do
   end function lin_diff_eqn
end module rmutil_linalg
