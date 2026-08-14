program test_simulate
  use mixsqp
  implicit none
  real(dp),allocatable::x(:),s(:),L(:,:),x2(:),s2(:),logL(:,:)
  call set_seed(12345)
  call simulate_mix_data(100,10,x,s,L)
  if(any(shape(L)/=[100,10])) error stop 'wrong simulation dimensions'
  if(any(L<0._dp) .or. maxval(L)>1._dp+1e-14_dp) error stop 'bad normalized likelihoods'
  if(abs(s(1))>1e-15_dp) error stop 'first scale should be zero'
  call set_seed(12345)
  call simulate_mix_data(100,10,x2,s2,logL,log_output=.true.,normalize_rows=.false.)
  if(maxval(abs(x-x2))>1e-15_dp) error stop 'seed reproducibility failed'
  if(maxval(abs(s-s2))>1e-14_dp) error stop 'scale reproducibility failed'
  print *, 'test_simulate: PASS'
end program
