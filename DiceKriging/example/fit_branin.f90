! Modern Fortran translation of the computational core of DiceKriging 1.6.1.
! Upstream DiceKriging is distributed under GPL-2 | GPL-3.
! This translation is distributed under the same license choice; see
! LICENSE-GPL-2 and LICENSE-GPL-3 in the project root.
program fit_branin
  use dicekriging
  implicit none
  real(dp) :: x(16,2), y(16), xnew(3,2), fnew(3,1)
  real(dp), allocatable :: f(:,:)
  type(km_model) :: model
  type(km_prediction) :: pred
  type(km_control) :: control
  integer :: i, j, k

  k=0
  do j=0,3
    do i=0,3
      k=k+1
      x(k,:)=[real(i,dp)/3.0_dp,real(j,dp)/3.0_dp]
      y(k)=branin(x(k,:))
    end do
  end do
  call trend_constant(x,f)
  control%multistart=6
  control%max_iter=400
  call km_fit(model,x,y,f,'matern5_2',control=control)

  xnew=reshape([0.2_dp,0.5_dp,0.8_dp, 0.3_dp,0.5_dp,0.7_dp],[3,2])
  fnew=1.0_dp
  call km_predict(model,xnew,fnew,'UK',pred,se_compute=.true.)

  print '(a,2f12.6)', 'estimated ranges: ',model%covariance%range
  do i=1,size(xnew,1)
    print '(i3,3f14.6)',i,pred%mean(i),pred%sd(i),pred%upper95(i)-pred%lower95(i)
  end do
end program fit_branin
