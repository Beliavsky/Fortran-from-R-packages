program test_mem_models
  use rumidas
  implicit none
  integer, parameter :: n=6,k=2
  real(dp) :: x(n), ret(n), mv(k+1,n), z(n), p(2), mu(n), meanx
  real(dp) :: pmx(6)
  real(dp),allocatable :: ll(:),pred(:),lr(:),sr(:),wrapped(:)
  type(mem_spec)::spec
  integer::i,status

  x=[1.0_dp,1.3_dp,0.8_dp,1.2_dp,0.9_dp,1.1_dp]
  ret=[0.01_dp,-0.02_dp,0.01_dp,-0.01_dp,0.02_dp,-0.015_dp]
  mv=0.2_dp
  z=[0.1_dp,0.2_dp,0.1_dp,0.0_dp,0.2_dp,0.1_dp]
  p=[0.2_dp,0.6_dp]
  spec=mem_spec(RUMIDAS_MEM,0,.false.)
  call mem_evaluate(p,spec,x,ll,pred,lr,sr,status)
  call check(status==0,'MEM status')
  meanx=sum(x)/real(n,dp)
  mu(1)=meanx
  do i=2,n
    mu(i)=(1.0_dp-p(1)-p(2))*meanx+p(1)*x(i-1)+p(2)*mu(i-1)
  end do
  call check(maxval(abs(pred-mu))<1.0e-13_dp,'MEM recursion')
  call check(maxval(abs(ll-(-log(mu)-x/mu)))<1.0e-13_dp,'MEM likelihood')
  call mem_pred_no_skew(p,x,wrapped,status)
  call check(maxval(abs(wrapped-pred))<1.0e-13_dp,'MEM wrapper')

  pmx=[0.15_dp,0.65_dp,log(1.1_dp),0.2_dp,2.0_dp,0.05_dp]
  spec=mem_spec(RUMIDAS_MEM_MIDAS_X,k,.false.)
  call mem_evaluate(pmx,spec,x,ll,pred,lr,sr,status,mv_m=mv,z_variable=z)
  call check(status==0,'MEM MIDAS X status')
  call check(all(pred>0.0_dp),'MEM MIDAS X positive')
  call check(all(abs(pred-sr*lr)<1.0e-13_dp),'MEM factorization')

  print '(a)', 'test_mem_models: PASS'
contains
  subroutine check(condition,message)
    logical,intent(in)::condition
    character(len=*),intent(in)::message
    if(.not.condition) error stop message
  end subroutine check
end program test_mem_models
