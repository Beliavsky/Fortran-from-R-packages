program test_history
  use adagio
  implicit none
  type(history_buffer) :: h
  h%store_inputs=.true.
  h%nvars=2
  call h%record([1._dp,2._dp], 5._dp)
  call h%record([3._dp,4._dp], 25._dp)
  call check(h%ncalls==2,'history count')
  call check(maxval(abs(h%values-[5._dp,25._dp])) < 1e-15_dp,'history values')
  call check(maxval(abs(h%input(2,:)-[3._dp,4._dp])) < 1e-15_dp,'history inputs')
  call h%reset()
  call check(h%ncalls==0,'history reset')
  print *, 'test_history: PASS'
contains
  subroutine check(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if(.not.ok) then
      print *, 'FAIL: ',trim(msg)
      error stop 1
    end if
  end subroutine
end program test_history
