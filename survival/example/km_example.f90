program km_example
  use survival
  implicit none
  real(dp) :: time(6)=[1._dp,2._dp,2._dp,3._dp,4._dp,5._dp]
  integer :: status(6)=[1,1,0,1,0,1]
  type(survfit_result) :: fit
  integer :: i
  call kaplan_meier(time,status,fit)
  print '(a)', ' time       risk      events    survival'
  do i=1,size(fit%time)
    print '(4f10.4)',fit%time(i),fit%n_risk(i),fit%n_event(i),fit%survival(i)
  end do
end program
