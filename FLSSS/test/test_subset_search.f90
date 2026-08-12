program test_subset_search
  use flsss_api
  implicit none
  type(subset_solutions) :: r
  type(multiset_solutions) :: mr
  type(real_bucket) :: b(2)
  type(integerized_search_result) :: ir
  type(mflsss_decomposition) :: dec
  real(dp) :: v(6), m(5,2)
  integer :: i

  v = [(real(i,dp), i=1,6)]
  r = flsss(3, v, 10.0_dp, 0.0_dp, solution_need=3)
  if (r%size() /= 3) error stop "FLSSS count"
  do i=1,r%size()
    if (abs(sum(v(r%sol(i)%idx))-10.0_dp)>1.0e-12_dp) error stop "FLSSS validity"
  end do

  r = flsss(0, v, 7.0_dp, 0.0_dp, solution_need=1)
  if (r%size() /= 1) error stop "variable-size FLSSS"
  if (abs(sum(v(r%sol(1)%idx))-7.0_dp)>1.0e-12_dp) error stop "variable-size validity"

  allocate(b(1)%value(3),b(2)%value(3))
  b(1)%value=[1.0_dp,2.0_dp,3.0_dp]; b(2)%value=[10.0_dp,20.0_dp,30.0_dp]
  mr=flsss_multiset([1,1],b,22.0_dp,0.0_dp,solution_need=2)
  if(mr%size()/=1) error stop "multiset count"
  if(any(mr%sol(1)%bucket(1)%idx/=[2])) error stop "multiset bucket 1"
  if(any(mr%sol(1)%bucket(2)%idx/=[2])) error stop "multiset bucket 2"

  m(:,1)=[1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp]
  m(:,2)=[1.0_dp,4.0_dp,9.0_dp,16.0_dp,25.0_dp]
  r=mflsss_par(2,m,[5.0_dp,17.0_dp],[0.0_dp,0.0_dp],solution_need=1)
  if(r%size()/=1 .or. any(r%sol(1)%idx/=[1,4])) error stop "mFLSSS"
  r=mflsss_par_impose_bounds(2,m,[5.0_dp,17.0_dp],[0.0_dp,0.0_dp],[1,4],[1,4],solution_need=1)
  if(r%size()/=1 .or. any(r%sol(1)%idx/=[1,4])) error stop "mFLSSS bounds"

  ir=mflsss_par_integerized(2,m,[5.0_dp,17.0_dp],[0.1_dp,0.1_dp],precision_level=[16,16],solution_need=1)
  if(ir%solution%size()/=1) error stop "integerized search"
  if(ir%integerized%compressed_dim/=2) error stop "integerized metadata"

  dec=decompose_mflsss(2,m,[5.0_dp,17.0_dp],[0.0_dp,0.0_dp])
  r=mflsss_obj_run(dec%object(1),solution_need=1)
  if(r%size()/=1 .or. any(r%sol(1)%idx/=[1,4])) error stop "decompose/run"
  print *, "test_subset_search: PASS"
end program test_subset_search
