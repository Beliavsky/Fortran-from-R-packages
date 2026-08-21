program test_transrate_parsers
  use relsurv, only : dp, ratetable_type, transrate_hmd, transrate_hld
  implicit none
  type(ratetable_type)::hmd,hld
  integer::u,y,a,s
  character(len=64)::mf='test_hmd_male.txt',ff='test_hmd_female.txt',hf='test_hld.csv'
  real(dp)::expect
  open(newunit=u,file=mf,status='replace',action='write');write(u,'(A)')'Year Age qx'
  do y=2000,2001
    do a=0,110
      if(a==110)then;write(u,'(I0,1X,I0,1X,A)')y,a,'1'
      else;write(u,'(I0,1X,I0,1X,F6.3)')y,a,0.1_dp;end if
    end do
  end do
  close(u)
  open(newunit=u,file=ff,status='replace',action='write');write(u,'(A)')'Year Age qx'
  do y=2000,2001
    do a=0,110
      write(u,'(I0,1X,I0,1X,F6.3)')y,a,0.2_dp
    end do
  end do
  close(u)
  hmd=transrate_hmd(mf,ff)
  if(any(hmd%dims/=[111,2,2]))error stop 'HMD dims'
  expect=-log(0.9_dp)/365.241_dp
  call assert_close(hmd%rate(1),expect,1.0e-14_dp,'HMD male rate')

  open(newunit=u,file=hf,status='replace',action='write')
  write(u,'(A)')'Country,Year1,Year2,TypeLT,Sex,Age,AgeInt,qx'
  do y=2000,2001
    do s=1,2
      do a=0,2
        write(u,'(A,",",I0,",",I0,",1,",I0,",",I0,",1.0,",F6.3)') &
          'X',y,y,s,a,0.05_dp*real(s+a+1,dp)
      end do
    end do
  end do
  close(u)
  hld=transrate_hld([hf])
  if(any(hld%dims/=[3,2,2]))error stop 'HLD dims'
  expect=-log(1.0_dp-0.10_dp)/365.241_dp
  call assert_close(hld%rate(1),expect,1.0e-14_dp,'HLD first rate')
  call execute_command_line('rm -f '//mf//' '//ff//' '//hf)
  print *, 'test_transrate_parsers: PASS'
contains
  subroutine assert_close(x,y,tol,msg)
    real(dp),intent(in)::x,y,tol;character(len=*),intent(in)::msg
    if(abs(x-y)>tol)then;print *,'FAIL ',msg,x,y;error stop 1;end if
  end subroutine assert_close
end program test_transrate_parsers
