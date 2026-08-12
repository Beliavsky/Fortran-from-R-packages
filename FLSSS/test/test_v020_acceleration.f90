program test_v020_acceleration
  use flsss_mod
  implicit none
  type(subset_solutions) :: r, plain, fast
  type(integerized_search_result) :: ir
  type(mflsss_decomposition) :: dec
  type(arb_flsss_decomposition) :: adec
  type(ksum_table) :: tab
  real(dp) :: x(30), mv(30,3)
  character(len=32) :: av(18,2), at(2)
  integer :: i

  do i=1,30
    x(i)=real(mod(17*i+11,97),dp)
    mv(i,1)=x(i)
    mv(i,2)=real(mod(29*i+7,113),dp)
    mv(i,3)=real(mod(41*i+3,127),dp)
  end do

  ! The exact completion envelope should reject this impossible target at the root.
  r=flsss(8,x,10000.0_dp,0.0_dp,solution_need=1)
  if(r%size()/=0) error stop "v0.2 1d impossible target"
  if(r%nodes>2_i8 .or. r%pruned<1_i8) error stop "v0.2 1d envelope not active"

  r=mflsss_par(8,mv,[10000.0_dp,10000.0_dp,10000.0_dp], &
    [0.0_dp,0.0_dp,0.0_dp],solution_need=1)
  if(r%size()/=0) error stop "v0.2 md impossible target"
  if(r%nodes>2_i8 .or. r%pruned<1_i8) error stop "v0.2 md envelope not active"

  ir=mflsss_par_integerized(8,mv,[10000.0_dp,10000.0_dp,10000.0_dp], &
    [0.01_dp,0.01_dp,0.01_dp],precision_level=[1000,1000,1000],solution_need=1)
  if(ir%solution%size()/=0) error stop "v0.2 int64 impossible target"
  if(ir%solution%pruned<1_i8) error stop "v0.2 int64 envelope not active"

  ! approx_ninstance now controls exact prefix-range decomposition.
  dec=decompose_mflsss(8,mv,[0.0_dp,0.0_dp,0.0_dp],[0.0_dp,0.0_dp,0.0_dp], &
    approx_ninstance=4)
  if(size(dec%object)/=4) error stop "mFLSSS decomposition count"
  if(dec%object(1)%prefix_lo/=1) error stop "mFLSSS decomposition first prefix"
  if(dec%object(4)%prefix_hi/=23) error stop "mFLSSS decomposition final prefix"

  do i=1,18
    write(av(i,1),'(i0)') i
    write(av(i,2),'(i0)') i*i
  end do
  at(1)='93'
  at(2)='1459'
  plain=arb_flsss(6,av,at,solution_need=1)
  if(plain%size()/=1) error stop "arb plain count"
  if(any(plain%sol(1)%idx/=[13,14,15,16,17,18])) error stop "arb plain solution"

  tab=build_ksum_hash(3,av,max_entries=5000)
  fast=arb_flsss(6,av,at,solution_need=1,given_ksum=tab)
  if(fast%size()/=1) error stop "arb hash count"
  if(any(fast%sol(1)%idx/=[13,14,15,16,17,18])) error stop "arb hash solution"
  if(fast%hash_lookups<=0_i8) error stop "arb hash not used"
  if(fast%nodes>=plain%nodes) error stop "arb hash did not reduce search nodes"

  adec=decompose_arb_flsss(6,av,at,approx_ninstance=3)
  if(size(adec%object)/=3) error stop "arb decomposition count"
  if(adec%object(1)%prefix_lo/=1 .or. adec%object(3)%prefix_hi/=13) &
    error stop "arb decomposition coverage"

  print *, "test_v020_acceleration: PASS", plain%nodes, fast%nodes, fast%hash_lookups
end program test_v020_acceleration
