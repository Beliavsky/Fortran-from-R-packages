! SPDX-License-Identifier: GPL-2.0-or-later
module nilde_diophantine
   use nilde_kinds, only : i8
   use nilde_types, only : integer_solutions_t
   use nilde_collect, only : append_i8_column
   implicit none
   private
   public :: nlde, get_subsetsum, enumerate_equation
contains

   function nlde(a, n, m, at_most, option) result(res)
      integer(i8), intent(in) :: a(:)
      integer(i8), intent(in) :: n
      integer, intent(in), optional :: m
      logical, intent(in), optional :: at_most
      integer, intent(in), optional :: option
      type(integer_solutions_t) :: res
      integer(i8), allocatable :: aa(:), bounds(:), full(:,:)
      integer :: mm, opt, j, ns, l
      logical :: am

      l = size(a)
      if (l < 2) error stop "nlde: a must have length at least 2"
      if (any(a <= 0_i8) .or. n <= 0_i8) error stop "nlde: coefficients and n must be positive"
      am = .true.; if (present(at_most)) am = at_most
      opt = 0; if (present(option)) opt = option
      mm = int(n / minval(a)); if (present(m)) mm = m
      if (mm <= 0) error stop "nlde: M must be positive"

      if (opt > 1) then
         allocate(aa(l+1), bounds(l+1))
         aa(1:l) = a; aa(l+1) = 1_i8
         bounds(1:l) = 1_i8
         bounds(l+1) = n
         mm = int(n)
         call enumerate_equation(aa, n, mm, am, bounds, full, ns)
         if (ns > 0) then
            allocate(res%x(l,ns))
            res%x = full(1:l,1:ns)
            ! Filter the original variables to 0/1, matching the R wrapper.
            res%nsol = 0
            do j = 1, ns
               if (all(res%x(:,j) == 0_i8 .or. res%x(:,j) == 1_i8)) then
                  res%nsol = res%nsol + 1
                  res%x(:,res%nsol) = res%x(:,j)
               end if
            end do
            if (res%nsol == 0) then
               deallocate(res%x)
            else if (res%nsol < ns) then
               res%x = res%x(:,1:res%nsol)
            end if
         end if
         return
      end if

      allocate(bounds(l)); bounds = n/minval(a)
      if (opt == 1) bounds = 1_i8
      call enumerate_equation(a, n, mm, am, bounds, res%x, res%nsol)
      if (res%nsol > 0 .and. opt == 1) then
         ! bounds already enforce binary; no further action.
      end if
   end function nlde

   function get_subsetsum(a, n, m, problem, bounds) result(res)
      integer(i8), intent(in) :: a(:)
      integer(i8), intent(in) :: n
      integer, intent(in), optional :: m
      character(len=*), intent(in), optional :: problem
      integer(i8), intent(in), optional :: bounds(:)
      type(integer_solutions_t) :: res
      integer(i8), allocatable :: b(:)
      integer :: mm
      character(len=:), allocatable :: p

      if (size(a) < 2) error stop "get_subsetsum: a must have length at least 2"
      if (any(a <= 0_i8) .or. n <= 0_i8) error stop "get_subsetsum: positive integer data required"
      mm = int(n/minval(a))
      if (present(m)) then
         mm = m
         if (mm > size(a)) error stop "get_subsetsum: M must not exceed length(a)"
      end if
      if (mm <= 0) error stop "get_subsetsum: invalid M"
      p = 'subsetsum01'; if (present(problem)) p = trim(problem)
      allocate(b(size(a)))
      select case (p)
      case ('subsetsum01')
         b = 1_i8
      case ('bsubsetsum')
         if (.not. present(bounds)) error stop "get_subsetsum: bounds required"
         if (size(bounds) /= size(a)) error stop "get_subsetsum: bounds length mismatch"
         if (any(bounds < 0_i8)) error stop "get_subsetsum: nonnegative bounds required"
         b = bounds
      case default
         error stop "get_subsetsum: unknown problem"
      end select
      call enumerate_equation(a, n, mm, .true., b, res%x, res%nsol)
   end function get_subsetsum

   subroutine enumerate_equation(a, n, m, at_most, bounds, sol, nsol)
      integer(i8), intent(in) :: a(:), n, bounds(:)
      integer, intent(in) :: m
      logical, intent(in) :: at_most
      integer(i8), allocatable, intent(out) :: sol(:,:)
      integer, intent(out) :: nsol
      integer(i8), allocatable :: as(:), bs(:), x(:), store(:,:)
      integer, allocatable :: ord(:), inv(:)
      integer :: l, j

      l = size(a)
      if (size(bounds) /= l) error stop "enumerate_equation: bounds mismatch"
      allocate(ord(l), inv(l), as(l), bs(l), x(l))
      call stable_argsort_i8(a, ord)
      do j = 1, l
         as(j) = a(ord(j)); bs(j) = bounds(ord(j)); inv(ord(j)) = j
      end do
      x = 0_i8; nsol = 0
      call rec(l, n, 0)
      if (nsol > 0) then
         allocate(sol(l,nsol))
         do j = 1, l
            sol(j,:) = store(inv(j),1:nsol)
         end do
      end if

   contains
      recursive subroutine rec(idx, rem, used)
         integer, intent(in) :: idx, used
         integer(i8), intent(in) :: rem
         integer(i8) :: xmax, xv
         integer :: used2
         if (idx == 1) then
            if (mod(rem,as(1)) /= 0_i8) return
            xv = rem/as(1)
            if (xv < 0_i8 .or. xv > bs(1)) return
            used2 = used + int(xv)
            if (at_most) then
               if (used2 > m) return
            else
               if (used2 /= m) return
            end if
            x(1) = xv
            call append_i8_column(store, nsol, x)
            return
         end if
         if (rem < 0_i8 .or. used > m) return
         xmax = min(bs(idx), rem/as(idx))
         xmax = min(xmax, int(m-used,i8))
         do xv = 0_i8, xmax
            x(idx) = xv
            call rec(idx-1, rem-xv*as(idx), used+int(xv))
         end do
         x(idx) = 0_i8
      end subroutine rec
   end subroutine enumerate_equation

   subroutine stable_argsort_i8(v, ord)
      integer(i8), intent(in) :: v(:)
      integer, intent(out) :: ord(:)
      integer :: i, j, t
      do i = 1, size(v); ord(i)=i; end do
      do i = 2, size(v)
         t = ord(i); j = i-1
         do while (j >= 1)
            if (v(ord(j)) <= v(t)) exit
            ord(j+1)=ord(j); j=j-1
         end do
         ord(j+1)=t
      end do
   end subroutine stable_argsort_i8

end module nilde_diophantine
