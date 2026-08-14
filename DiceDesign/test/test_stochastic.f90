program test_stochastic
  use iso_fortran_env, only : int64
  use dicedesign, only : dp, dmax_result, strauss_result, wsp_result, dmax_design, strauss_design, wsp_design
  implicit none
  type(dmax_result) :: dm
  type(strauss_result) :: st
  type(wsp_result) :: ws
  real(dp) :: grid(9,2)
  integer :: i,j,k

  call dmax_design(7,2,0.0_dp,80,dm,seed=123_int64)
  if (dm%det_end+1e-13_dp<dm%det_init) error stop 'dmax determinant decreased'
  if (minval(dm%design)<0.0_dp .or. maxval(dm%design)>1.0_dp) error stop 'dmax bounds'

  call strauss_design(10,2,0.25_dp,st,alpha=0.5_dp,repulsion=0.01_dp,nmc=20, &
    constraints1d=.true.,repulsion1d=0.1_dp,seed=77_int64)
  if (minval(st%design)<0.0_dp .or. maxval(st%design)>1.0_dp) error stop 'strauss bounds'
  if (maxval(abs(st%design-st%design_init))<=0.0_dp) error stop 'strauss unchanged'

  k=0
  do i=0,2
    do j=0,2
      k=k+1
      grid(k,:)=[real(i,dp)/2.0_dp,real(j,dp)/2.0_dp]
    end do
  end do
  call wsp_design(grid,0.6_dp,ws)
  do i=1,size(ws%design,1)-1
    do j=i+1,size(ws%design,1)
      if (sqrt(sum((ws%design(i,:)-ws%design(j,:))**2))<0.6_dp-1e-14_dp) error stop 'wsp separation'
    end do
  end do
  if (size(ws%design,1)+size(ws%residual_design,1)/=9) error stop 'wsp partition'

  print *, 'test_stochastic: PASS'
end program test_stochastic
