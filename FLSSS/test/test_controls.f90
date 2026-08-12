program test_controls
  use flsss_mod
  implicit none
  type(subset_solutions) :: r
  type(integerized_search_result) :: ir
  real(dp) :: v(8), m(4,2)
  integer :: i
  v=[(real(i,dp),i=1,8)]
  r=flsss(3,v,12.0_dp,0.0_dp,solution_need=20,lb=[1,3,6],ub=[2,5,8])
  do i=1,r%size()
    if(any(r%sol(i)%idx<[1,3,6]) .or. any(r%sol(i)%idx>[2,5,8])) error stop "bounded controls"
  end do
  m=reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp, 2.0_dp,4.0_dp,6.0_dp,8.0_dp],[4,2])
  ir=mflsss_par_impose_bounds_integerized(2,m,[5.0_dp,10.0_dp],[0.1_dp,0.1_dp], &
      [1,2],[3,4],precision_level=[-1,-1],return_before_mining=.true.)
  if(ir%solution%size()/=0 .or. ir%integerized%compressed_dim/=2) error stop "return-before-mining"
  r=mflssspar(2,2,m,[5.0_dp,10.0_dp],[0.0_dp,0.0_dp],solution_need=1)
  if(r%size()/=1) error stop "compat mFLSSSpar"
  print *, "test_controls: PASS"
end program test_controls
