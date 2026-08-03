! SPDX-License-Identifier: Artistic-2.0
program test_options
  use ecd_api
  implicit none
  real(dp) :: c1,c2,iv,k(12),price(12),kn(5),pn(5)
  integer :: i,st
  c1=bs_call_price(0.128886_dp,2100.0_dp,2089.27_dp,1.0_dp/365.0_dp,0.0_dp,0.019_dp)
  c2=bs_call_price(0.294296_dp,2040.0_dp,2089.27_dp,1.0_dp/365.0_dp,0.0_dp,0.019_dp)
  call close(c1,1.7994517438573894_dp,2.0e-13_dp,'BS example 1')
  call close(c2,49.99998767559032_dp,2.0e-13_dp,'BS example 2')
  iv=bs_implied_volatility(c1,2100.0_dp,2089.27_dp,1.0_dp/365.0_dp, &
    0.0_dp,0.019_dp,'c',st)
  call check(st==ecd_ok,'implied volatility status')
  call close(iv,0.128886_dp,2.0e-9_dp,'implied volatility inversion')
  call close(bs_put_price(0.2_dp,100.0_dp,100.0_dp,1.0_dp,0.05_dp,0.02_dp), &
    bs_call_price(0.2_dp,100.0_dp,100.0_dp,1.0_dp,0.05_dp,0.02_dp) &
    -100.0_dp*exp(-0.02_dp)+100.0_dp*exp(-0.05_dp),2.0e-14_dp,'put-call parity')

  do i=1,size(k)
    k(i)=-0.6_dp+1.2_dp*real(i-1,dp)/real(size(k)-1,dp)
    price(i)=exp(1.0_dp+0.2_dp*k(i)+0.3_dp*k(i)**2)
  end do
  kn=[-0.4_dp,-0.1_dp,0.0_dp,0.2_dp,0.5_dp]
  call polyfit_option(k,price,0.0_dp,kn,pn,2,2,st)
  call check(st==ecd_ok,'polyfit status')
  do i=1,size(kn)
    call close(pn(i),exp(1.0_dp+0.2_dp*kn(i)+0.3_dp*kn(i)**2),2.0e-11_dp,'polyfit prediction')
  end do
  print '(a)', 'test_options: PASS'
contains
  subroutine close(x,y,tol,msg)
    real(dp),intent(in)::x,y,tol
    character(len=*),intent(in)::msg
    if(abs(x-y)>tol*max(1.0_dp,abs(y)))then
      write(*,*)trim(msg),x,y; error stop 1
    end if
  end subroutine close
  subroutine check(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if(.not.ok)then;write(*,*)trim(msg);error stop 1;end if
  end subroutine check
end program test_options
