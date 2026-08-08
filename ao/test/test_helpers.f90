program test_helpers
  use ao
  implicit none
  type(ao_block),allocatable::blocks(:)
  type(ao_real_vector),allocatable::parts(:)
  integer,allocatable::seen(:)
  real(dp)::x(5)
  logical::ok
  integer::i,j
  x=[1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp]
  call split_estimate(x,[2,2,1],parts,ok)
  if(.not.ok .or. size(parts)/=3) error stop 'split failed'
  if(any(abs(parts(2)%value-[3.0_dp,4.0_dp])>0.0_dp)) error stop 'split values wrong'
  call ao_seed(123)
  call generate_random_partition(10,0.5_dp,2,blocks)
  if(size(blocks)<2) error stop 'minimum block count failed'
  allocate(seen(10)); seen=0
  do i=1,size(blocks)
    do j=1,size(blocks(i)%index)
      seen(blocks(i)%index(j))=seen(blocks(i)%index(j))+1
    end do
  end do
  if(any(seen/=1)) error stop 'partition is not a partition'
  call generate_random_partition(5,0.5_dp,5,blocks)
  if(size(blocks)/=5) error stop 'singleton partition failed'
  print *, 'PASS test_helpers'
end program
