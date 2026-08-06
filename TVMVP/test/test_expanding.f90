program test_expanding
  use tvmvp, only : dp, expanding_window_result, expanding_tvmvp
  implicit none
  integer, parameter :: n=44,p=4
  real(dp) :: x(n,p),rf(1),grid(5)
  integer :: i,j
  type(expanding_window_result) :: result
  do i=1,n
    do j=1,p
      x(i,j)=0.001_dp*real(j,dp)+0.012_dp*sin(0.14_dp*real(i,dp))*(0.5_dp+0.2_dp*j)+ &
             0.004_dp*cos(real(i*j,dp))
    end do
  end do
  do i=1,5
    grid(i)=0.01_dp+0.2_dp*real(i-1,dp)
  end do
  rf=0.0001_dp
  call expanding_tvmvp(x,24,5,1,result,return_type='daily',rf=rf,m0=3,rho_grid=grid, &
                       source_compatible_expanding=.false.)
  call check(.not.result%error%failed(),'expanding-window status')
  call check(size(result%tvmvp_returns)==n-24,'return length')
  do i=1,size(result%weights,2)
    call check(abs(sum(result%weights(:,i))-1.0_dp)<1.0e-8_dp,'weight budget')
  end do
  call check(result%tvmvp_metrics%annualized_standard_deviation>=0.0_dp,'annualized risk')
  call check(sum(result%holding_lengths)==n-24,'holding lengths')
  print *, 'test_expanding: PASS'
contains
  subroutine check(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if(.not.ok) then
      print *, 'FAIL: ',msg
      error stop 1
    end if
  end subroutine check
end program test_expanding
