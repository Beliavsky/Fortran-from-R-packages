program test_louis_analytic
  use flexsurv_kinds, only : dp
  use flexsurv_fit, only : flexsurv_data,flexsurv_spec,flexsurv_result,prepare_survival_data,initialize_spec
  use flexsurv_distributions, only : dist_exponential
  use flexsurv_mixture_full, only : louis_information,louis_information_numeric
  implicit none
  integer,parameter::n=5,k=2
  type(flexsurv_data)::dat,dc(k)
  type(flexsurv_spec)::sp(k)
  type(flexsurv_result)::comp(k)
  real(dp)::tt(n),px(n,1),post(n,k),alpha(1),pb(1,1),p2,f1,f2,den
  real(dp),allocatable::ca(:,:),cn(:,:)
  logical::allow(n,k)
  integer::ev(n),i,sa,sn
  tt=[0.2_dp,0.6_dp,1.0_dp,1.6_dp,2.4_dp];ev=0;allow=.true.
  call prepare_survival_data(dat,tt,[1,1,1,1,1]);dc(1)=dat;dc(2)=dat
  call initialize_spec(sp(1),dist_exponential,n,[0.5_dp]);call initialize_spec(sp(2),dist_exponential,n,[1.5_dp])
  comp(1)%theta=[log(0.5_dp)];comp(2)%theta=[log(1.5_dp)]
  alpha=[0.1_dp];pb(1,1)=0.2_dp;px(:,1)=[-1.0_dp,-0.5_dp,0.0_dp,0.5_dp,1.0_dp]
  do i=1,n
    p2=1.0_dp/(1.0_dp+exp(-(alpha(1)+pb(1,1)*px(i,1))))
    f1=0.5_dp*exp(-0.5_dp*tt(i));f2=1.5_dp*exp(-1.5_dp*tt(i));den=(1.0_dp-p2)*f1+p2*f2
    post(i,1)=(1.0_dp-p2)*f1/den;post(i,2)=p2*f2/den
  end do
  call louis_information(dat,dc,sp,ev,allow,px,alpha,pb,comp,post,ca,sa)
  call louis_information_numeric(dat,dc,sp,ev,allow,px,alpha,pb,comp,post,cn,sn)
  if(sa/=0.or.sn/=0)error stop 'louis status'
  if(maxval(abs(ca-transpose(ca)))>1.0e-8_dp)error stop 'louis symmetry'
  if(maxval(abs(ca-cn))/max(1.0_dp,maxval(abs(cn)))>3.0e-3_dp)error stop 'analytic/numeric louis'
  print *,'test_louis_analytic: PASS'
end program
