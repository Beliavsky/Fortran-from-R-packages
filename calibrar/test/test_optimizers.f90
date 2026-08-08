program test_optimizers
  use calibrar, only : dp, optim_options, optim_result, optim2, quadratic_objective, quadratic_gradient
  implicit none
  type(optim_options) :: op
  type(optim_result) :: res
  real(dp) :: x0(2),lo(2),hi(2)
  character(len=12),parameter :: methods(6)=[character(len=12)::"BFGS","CG","Nelder-Mead","hjn","spg","AHR-ES"]
  integer :: i
  x0=[2.0_dp,-1.0_dp];lo=-5.0_dp;hi=5.0_dp;op%maxit=400;op%maxfeval=10000;op%reltol=1.0e-7_dp
  do i=1,size(methods)
    if(trim(methods(i))=="AHR-ES") then
      op%maxit=160
      op%seed=441
    else
      op%maxit=400
    end if
    if(trim(methods(i))=="BFGS" .or. trim(methods(i))=="CG" .or. trim(methods(i))=="spg") then
      call optim2(x0,quadratic_objective,res,trim(methods(i)),lo,hi,op,quadratic_gradient)
    else
      call optim2(x0,quadratic_objective,res,trim(methods(i)),lo,hi,op)
    end if
    if(maxval(abs(res%par))>merge(1.0e-1_dp,2.0e-3_dp,trim(methods(i))=="AHR-ES")) then
      print *, trim(methods(i)), res%par, res%value
      error stop "optimizer failed"
    end if
  end do
  print *, "PASS test_optimizers"
end program test_optimizers
