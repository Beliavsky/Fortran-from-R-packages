! SPDX-License-Identifier: MIT
program test_mi
  use bayesianou
  implicit none
  integer,parameter::t=30,s=1,d=3
  type(ou_input)::inp
  type(ou_options)::opt
  type(ou_mi_result)::mi
  real(dp)::draws(t,s,d)
  integer::i,k
  allocate(inp%y(t,s),inp%x(t,s),inp%tmg(t),inp%com(t,s),inp%capital(t,s))
  do i=1,t
    inp%x(i,1)=sin(0.1_dp*i);inp%tmg(i)=cos(0.13_dp*i);inp%com(i,1)=0.5_dp;inp%capital(i,1)=100+i
    inp%y(i,1)=0.4_dp*inp%x(i,1)+0.05_dp*sin(0.7_dp*i)
  end do
  do k=1,d;draws(:,:,k)=inp%y+0.01_dp*real(k-2,dp);end do
  opt%n_levels=1;opt%chains=2;opt%iterations=50;opt%warmup=25;opt%train_frac=0.7_dp
  call fit_ou_nested_mi(draws,inp,opt,3,.false.,mi)
  if(mi%status/=status_ok)error stop 'MI failed'
  if(.not.allocated(mi%pooled%estimate))error stop 'MI pooling missing'
  if(any(mi%pooled%total_sd<0))error stop 'MI variance invalid'
  print *, 'test_mi: PASS'
end program test_mi
