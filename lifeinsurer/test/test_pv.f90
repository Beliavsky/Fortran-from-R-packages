! SPDX-License-Identifier: GPL-2.0-or-later
program test_pv
  use lifeinsurer
  implicit none
  integer,parameter::n=10
  real(dp)::q(n),p(n),adv(n),arr(n),v
  real(dp),allocatable::pv(:),expected(:)
  type(frequency_correction)::c
  integer::i
  q=[(0.01_dp*real(i,dp),i=1,n)]; p=1.0_dp-q; adv=1.0_dp; arr=0.0_dp; v=1.0_dp/1.03_dp
  c%alpha=1.0_dp;c%beta=0.0_dp
  pv=pv_survival(adv,arr,p,v,1,c)
  allocate(expected(n)); expected(n)=1.0_dp
  do i=n-1,1,-1; expected(i)=1.0_dp+v*p(i)*expected(i+1); end do
  call check(maxval(abs(pv-expected))<1e-12_dp,'survival recursion')
  pv=pv_death(adv,q,p,v)
  expected(n)=v*q(n)
  do i=n-1,1,-1; expected(i)=v*q(i)+v*p(i)*expected(i+1); end do
  call check(maxval(abs(pv-expected))<1e-12_dp,'death recursion')
  pv=pv_guaranteed(adv,arr,v,1,c)
  expected(n)=1.0_dp
  do i=n-1,1,-1; expected(i)=1.0_dp+v*expected(i+1); end do
  call check(maxval(abs(pv-expected))<1e-12_dp,'guaranteed recursion')
  print '(a)','test_pv: PASS'
contains
  subroutine check(ok,msg); logical,intent(in)::ok; character(*),intent(in)::msg
    if(.not.ok) then; print '(a,1x,a)','FAIL:',msg; error stop 1; end if
  end subroutine
end program
