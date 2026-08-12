program test_qaplib_formats
   use qap
   implicit none
   character(len=16), parameter :: names(8) = [character(len=16) :: &
      'esc128', 'kra32', 'lipa90b', 'sko100a', 'tai256c', 'tho150', 'wil100', 'nug30']
   type(qap_problem_t) :: p
   integer :: i, stat
   character(len=:), allocatable :: msg
   character(len=256) :: path

   do i = 1, size(names)
      path = 'data/qaplib/' // trim(names(i)) // '.dat'
      call read_qaplib(trim(path), p, stat, msg)
      if (stat /= 0) then
         write(*,*) trim(path), ': ', msg
         error stop 'QAPLIB read failed'
      end if
      if (size(p%A,1) <= 0) error stop 'empty QAPLIB problem'
      if (size(p%A,1) /= size(p%B,1)) error stop 'QAPLIB size mismatch'
      if (p%has_solution) then
         if (.not. qap_is_permutation(p%solution)) error stop 'bad QAPLIB permutation'
      end if
   end do
   print *, 'test_qaplib_formats: PASS'
end program test_qaplib_formats
