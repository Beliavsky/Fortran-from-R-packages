program test_magic_squares
   use magic
   implicit none
   integer(ik), allocatable :: a(:, :), product_square(:, :)
   integer :: n

   do n = 3, 16
      a = magic_square_of_order(n)
      call check(is_magic(a), "magic square invariant")
      call check(is_normal_square(a), "normal square invariant")
      call check(sum(a(1, :)) == magic_constant(n), "magic constant")
   end do

   call check(is_magic(magic_2np1(4)), "odd construction")
   call check(is_magic(magic_4n(3)), "doubly-even construction")
   call check(is_magic(magic_4np2(2)), "singly-even construction")
   call check(is_magic(strachey_square(3)), "Strachey construction")
   call check(is_magic(lozenge_square(3)), "lozenge construction")

   product_square = magic_product_fast(magic_square_of_order(3), magic_square_of_order(4))
   call check(is_magic(product_square), "compound magic product")
   call check(is_normal_square(product_square), "normal compound product")

   print '(a)', "test_magic_squares: PASS"
contains
   subroutine check(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not. condition) then
         write(*, '(a)') "FAIL: " // message
         error stop 1
      end if
   end subroutine check
end program test_magic_squares
