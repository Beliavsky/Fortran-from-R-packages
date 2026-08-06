program test_hypothesis
  use tvmvp, only : dp, hypothesis_result, hyptest, compute_sigma_0, tvmvp_error
  implicit none
  integer, parameter :: n=22,p=4
  real(dp) :: x(n,p),res(n,p)
  real(dp), allocatable :: sigma0(:,:)
  integer :: i,j
  type(hypothesis_result) :: result
  type(tvmvp_error) :: err
  do i=1,n
    do j=1,p
      x(i,j)=0.02_dp*sin(0.17_dp*real(i*j,dp))+0.01_dp*cos(0.11_dp*real(i+2*j,dp))+ &
             0.004_dp*real(j,dp)*sin(0.05_dp*real(i,dp))
      res(i,j)=sin(real(i+j,dp))
    end do
  end do
  call compute_sigma_0(res,sigma0,err)
  call check(.not.err%failed(),'sigma0 status')
  call check(maxval(abs(sigma0-transpose(sigma0)))<1.0e-13_dp,'sigma0 symmetry')
  call check(abs(sigma0(1,3))<=abs(dot_product(res(:,1),res(:,3))/real(n,dp))+1.0e-13_dp, &
             'sigma0 taper')
  call hyptest(x,1,result,n_bootstrap=2,seed=90210)
  call check(.not.result%error%failed(),'hypothesis status')
  call check(result%variance_term>0.0_dp,'positive variance term')
  call check(result%p_value>=0.0_dp .and. result%p_value<=1.0_dp,'bootstrap p-value')
  call check(size(result%bootstrap_statistics)==2,'bootstrap size')
  print *, 'test_hypothesis: PASS'
contains
  subroutine check(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if(.not.ok) then
      print *, 'FAIL: ',msg
      error stop 1
    end if
  end subroutine check
end program test_hypothesis
