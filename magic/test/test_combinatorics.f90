program test_combinatorics
   use magic
   implicit none
   integer(ik), allocatable :: a(:, :), b(:, :)
   type(integer_tensor) :: inc, moved
   integer, allocatable :: seed(:)
   integer :: nseed, m, u

   call random_seed(size=nseed)
   allocate(seed(nseed))
   seed = 1729
   call random_seed(put=seed)

   a = latin_square(7)
   call check(is_latin_square(a), "cyclic Latin square")
   inc = incidence(a)
   call check(is_incidence(inc, .false.), "proper incidence tensor")
   b = unincidence(inc)
   call check(all(a == b), "incidence round trip")

   moved = another_incidence(inc, 100000)
   call check(is_incidence(moved, .false.), "incidence Markov move")
   call check(.not. tensor_equal(inc, moved), "different incidence tensor")
   b = unincidence(moved)
   call check(is_latin_square(b), "moved Latin square")

   do m = 2, 7
      do u = 1, m - 1
         a = sam_square(m, u)
         call check(is_sam(a), "sparse antimagic construction")
      end do
   end do

   a = sylvester_hadamard(5)
   call check(is_hadamard(a), "Sylvester Hadamard matrix")

   a = cilleruelo_square(2, 3)
   call check(is_multiplicative_magic(a), "Cilleruelo multiplicative magic square")

   print '(a)', "test_combinatorics: PASS"
contains
   subroutine check(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not. condition) then
         write(*, '(a)') "FAIL: " // message
         error stop 1
      end if
   end subroutine check
end program test_combinatorics
