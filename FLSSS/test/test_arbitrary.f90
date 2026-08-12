program test_arbitrary
  use flsss_api
  implicit none
  type(subset_solutions) :: r
  type(arb_flsss_decomposition) :: dec
  type(ksum_table) :: tab
  character(len=32) :: v(5,2), target(2), nums(3)
  character(len=:), allocatable :: s

  v='0'
  v(:,1)=[character(len=32)::'1.1','2.2','3.3','4.4','5.5']
  v(:,2)=[character(len=32)::'10','20','30','40','50']
  target=[character(len=32)::'5.5','50']
  r=arb_flsss(2,v,target,solution_need=3)
  if(r%size()/=2) error stop "arbFLSSS count"
  if(any(r%sol(1)%idx/=[1,4])) error stop "arbFLSSS first"
  if(any(r%sol(2)%idx/=[2,3])) error stop "arbFLSSS second"

  nums=[character(len=32)::'12345678901234567890.5','-20.25','3.75']
  s=add_num_strings(nums)
  if(s/='12345678901234567874') error stop "addNumStrings"

  dec=decompose_arb_flsss(2,v,target)
  r=arb_flsss_obj_run(dec%object(1),solution_need=1)
  if(r%size()/=1 .or. any(r%sol(1)%idx/=[1,4])) error stop "arb decompose/run"

  tab=build_ksum_hash(2,v(:,1:1))
  if(size(tab%index,1)/=10) error stop "ksumHash count"
  if(any(tab%index(1,:)/=[1,2])) error stop "ksumHash first combination"
  print *, "test_arbitrary: PASS"
end program test_arbitrary
