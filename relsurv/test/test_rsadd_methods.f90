program test_rsadd_methods
  use relsurv, only : dp, rsadd_result, rsadd_em_core, rsadd_glm_bin, rsadd_glm_poisson
  implicit none
  type(rsadd_result) :: fit
  real(dp) :: st(3),sp(3),x0(3,0),ph(3),pc(3),nie(3)
  integer :: status(3),cause(3)
  real(dp),allocatable :: cov(:,:),nd(:),ld(:),ps(:),kt(:),dstar(:),lny(:)
  integer,allocatable :: ii(:)
  integer :: i,n
  real(dp) :: eta,mu
  st=0.0_dp;sp=[1.0_dp,2.0_dp,3.0_dp];status=[1,1,0];cause=[1,1,2]
  ph=0.0_dp;pc=0.0_dp;nie=[1.0_dp,1.0_dp,0.0_dp]
  call rsadd_em_core(st,sp,status,x0,cause,ph,pc,nie,0.0_dp,fit,5,1.0e-10_dp)
  call assert_close(fit%lambda0_ns(1),1.0_dp/3.0_dp,1.0e-12_dp,'em lambda 1')
  call assert_close(fit%lambda0_ns(2),0.5_dp,1.0e-12_dp,'em lambda 2')
  call assert_close(fit%cumulative_lambda0(2),5.0_dp/6.0_dp,1.0e-12_dp,'em cumulative')

  n=12;allocate(cov(n,1),nd(n),ld(n),ps(n),kt(n),ii(n),dstar(n),lny(n))
  do i=1,n
    cov(i,1)=-1.0_dp+2.0_dp*real(i-1,dp)/real(n-1,dp);ii(i)=1+mod(i-1,2)
    ld(i)=20000.0_dp;ps(i)=0.95_dp;kt(i)=1.0_dp
    eta=0.35_dp*cov(i,1)+merge(-2.0_dp,-1.5_dp,ii(i)==1)
    mu=1.0_dp-exp(-exp(eta))*ps(i);nd(i)=real(nint(ld(i)*mu),dp)
  end do
  call rsadd_glm_bin(cov,ii,nd,ld,ps,kt,fit,100,1.0e-10_dp)
  call assert_close(fit%coef(1),0.35_dp,0.03_dp,'glm bin beta')

  do i=1,n
    dstar(i)=2.0_dp;lny(i)=log(100.0_dp)
    eta=0.25_dp*cov(i,1)+merge(-3.0_dp,-2.6_dp,ii(i)==1)+lny(i)
    mu=dstar(i)+exp(eta);nd(i)=real(nint(mu),dp)
  end do
  call rsadd_glm_poisson(cov,ii,nd,dstar,lny,fit,100,1.0e-10_dp)
  call assert_close(fit%coef(1),0.25_dp,0.10_dp,'glm poisson beta')
  print *, 'test_rsadd_methods: PASS'
contains
  subroutine assert_close(a,b,tol,msg)
    real(dp),intent(in)::a,b,tol;character(len=*),intent(in)::msg
    if(abs(a-b)>tol)then;print *,'FAIL ',msg,a,b;error stop 1;end if
  end subroutine assert_close
end program test_rsadd_methods
