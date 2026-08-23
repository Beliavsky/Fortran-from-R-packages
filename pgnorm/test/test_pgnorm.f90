program test_pgnorm
  use pgnorm
  use pgnorm_special, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer,parameter::n=40000
  real(dp)::y(n),m,v,pi,tol
  real(dp),allocatable::x(:)
  real(dp)::u(1000,2),normp(1000)
  integer::i
  pi=acos(-1.0_dp); tol=3.0e-3_dp

  call assert_close(dpgnorm(0.0_dp,2.0_dp),1.0_dp/sqrt(2.0_dp*pi),1.0e-12_dp,'normal density')
  call assert_close(ppgnorm(0.0_dp,2.0_dp),0.5_dp,1.0e-14_dp,'cdf at zero')
  call assert_close(dpgnorm(0.0_dp,1.0_dp),0.5_dp,1.0e-12_dp,'laplace density')
  call assert_close(ppgnorm(1.0_dp,1.0_dp),1.0_dp-0.5_dp*exp(-1.0_dp),1.0e-12_dp,'laplace cdf')
  call assert_close(dpgnorm(1.0_dp,2.0_dp,1.0_dp,2.0_dp),1.0_dp/(2.0_dp*sqrt(2.0_dp*pi)),1.0e-12_dp,'scale jacobian')

  call rpgnorm(n,y,p=3.0_dp,method='nardonpianca')
  m=sum(y)/n; v=sum((y-m)**2)/(n-1)
  if(abs(m)>0.03_dp) error stop 'sample mean'
  if(abs(v-pgnorm_sd(3.0_dp)**2)>0.04_dp) error stop 'sample variance'

  call rpgnorm(n,y,p=2.0_dp,method='pgenpolar')
  m=sum(y)/n; v=sum((y-m)**2)/(n-1)
  if(abs(m)>0.03_dp .or. abs(v-1.0_dp)>0.04_dp) error stop 'pgenpolar moments'

  call rpgnorm(n,y,p=1.5_dp,method='pgenpolarrej')
  m=sum(y)/n; v=sum((y-m)**2)/(n-1)
  if(abs(m)>0.04_dp .or. abs(v-pgnorm_sd(1.5_dp)**2)>0.06_dp) error stop 'pgenpolarrej moments'

  call rpgnorm(5000,y(1:5000),p=0.5_dp,method='montypython')
  if(any(.not.ieee_is_finite(y(1:5000)))) error stop 'montypython finite'
  call rpgnorm(5000,y(1:5000),p=3.0_dp,method='ziggurat')
  if(any(.not.ieee_is_finite(y(1:5000)))) error stop 'ziggurat finite'

  call rpgunif(1000,3.0_dp,u)
  do i=1,1000
    normp(i)=abs(u(i,1))**3+abs(u(i,2))**3
  end do
  if(maxval(abs(normp-1.0_dp))>1.0e-12_dp) error stop 'p-sphere uniform norm'

  call zigsetup(2.0_dp,x,n=64)
  if(size(x)/=63) error stop 'zigsetup size'
  if(any(x(2:)<=x(:size(x)-1))) error stop 'zigsetup monotonicity'

  print '(a)','test_pgnorm: PASS'
contains
  subroutine assert_close(a,b,t,msg)
    real(dp),intent(in)::a,b,t
    character(len=*),intent(in)::msg
    if(abs(a-b)>t) then
      print *,trim(msg),a,b
      error stop 1
    end if
  end subroutine
end program test_pgnorm
