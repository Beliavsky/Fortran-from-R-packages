program test_kernels_linalg
  use tvmvp, only : dp, tvmvp_error, epanechnikov_kernel, boundary_kernel, &
                    two_fold_convolution_kernel, symmetric_eigen, solve_linear
  implicit none
  real(dp), allocatable :: vals(:),vecs(:,:),x(:)
  real(dp) :: a(2,2),b(2)
  type(tvmvp_error) :: err
  call check(abs(epanechnikov_kernel(0.0_dp)-0.75_dp)<1.0e-14_dp,'Epanechnikov center')
  call check(abs(epanechnikov_kernel(1.1_dp))<1.0e-15_dp,'Epanechnikov support')
  call check(abs(two_fold_convolution_kernel(0.0_dp)-0.6_dp)<1.0e-13_dp,'convolution center')
  call check(abs(two_fold_convolution_kernel(2.0_dp))<1.0e-13_dp,'convolution endpoint')
  call check(abs(boundary_kernel(1,2,20,0.2_dp)-boundary_kernel(1,10,20,0.2_dp, &
             source_compatible=.false.))>1.0e-12_dp,'boundary compatibility switch')
  a=reshape([2.0_dp,1.0_dp,1.0_dp,2.0_dp],[2,2])
  call symmetric_eigen(a,vals,vecs,err)
  call check(.not.err%failed(),'eigen status')
  call check(maxval(abs(vals-[3.0_dp,1.0_dp]))<1.0e-10_dp,'eigen values')
  b=[3.0_dp,3.0_dp]
  call solve_linear(a,b,x,err)
  call check(.not.err%failed(),'solve status')
  call check(maxval(abs(x-[1.0_dp,1.0_dp]))<1.0e-12_dp,'linear solve')
  print *, 'test_kernels_linalg: PASS'
contains
  subroutine check(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if(.not.ok) then
      print *, 'FAIL: ',msg
      error stop 1
    end if
  end subroutine check
end program test_kernels_linalg
