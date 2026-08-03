program test_panmagic
   use magic
   implicit none
   integer(ik), allocatable :: a(:, :)
   integer :: m, n

   do m = 1, 4
      a = panmagic_4n(m)
      call check(is_panmagic(a), "panmagic 4n")
      call check(is_normal_square(a), "normal panmagic 4n")

      a = panmagic_6np1(m)
      call check(is_panmagic(a), "panmagic 6n+1")
      call check(is_normal_square(a), "normal panmagic 6n+1")

      a = panmagic_6nm1(m)
      call check(is_panmagic(a), "panmagic 6n-1")
      call check(is_normal_square(a), "normal panmagic 6n-1")
   end do

   do n = 5, 25, 2
      if (modulo(n, 6) == 1 .or. modulo(n, 6) == 5) then
         a = hudson_square(n)
         call check(is_panmagic(a), "Hudson panmagic construction")
         call check(is_normal_square(a), "Hudson normality")
      end if
   end do

   a = magic_prime(7)
   call check(is_magic(a), "prime-order construction")
   call check(is_normal_square(a), "prime-order normality")

   print '(a)', "test_panmagic: PASS"
contains
   subroutine check(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not. condition) then
         write(*, '(a)') "FAIL: " // message
         error stop 1
      end if
   end subroutine check
end program test_panmagic
