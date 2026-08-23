program test_v03_repeated
   use rfast
   implicit none
   integer, parameter :: n=30
   real(dp)::x(n,1),y(n),v(n)
   integer::id(n),i,g,failures
   type(random_intercept_result)::rr,ri
   type(variance_components_result)::vc
   failures=0
   do i=1,n
      if(i<=4)then;g=1
      else if(i<=10)then;g=2
      else if(i<=19)then;g=3
      else;g=4;end if
      id(i)=g
      x(i,1)=sin(real(i,dp)*0.61_dp)
      y(i)=1.2_dp+1.35_dp*x(i,1)+0.45_dp*real(g-2,dp)+0.08_dp*cos(real(i,dp)*1.17_dp)
      v(i)=y(i)-1.35_dp*x(i,1)
   end do
   rr=rint_reg(y,x,id,1.0e-9_dp,.true.,200)
   call check(rr%status==0,'rint_reg status')
   call check(abs(rr%beta(2)-1.35_dp)<0.08_dp,'rint_reg slope')
   call check(rr%sigma_error>0.0_dp.and.rr%sigma_tau>=0.0_dp,'rint variances')
   call check(allocated(rr%ranef).and.size(rr%ranef)==4,'rint ranef')
   call check(allocated(rr%se_beta).and.all(rr%se_beta>0.0_dp),'rint standard errors')

   ri=rint_mle(v,id,.true.,1.0e-9_dp,200)
   vc=varcomps_mle(v,id,.false.,1.0e-9_dp,200)
   call check(ri%status==0.and.vc%status==0,'rint_mle/varcomps status')
   call check(abs(ri%sigma_tau-vc%info(1,1))<1.0e-9_dp,'rint_mle between consistency')
   call check(abs(ri%sigma_error-vc%info(1,2))<1.0e-9_dp,'rint_mle within consistency')

   if(failures==0)then
      print *,'test_v03_repeated: PASS'
   else
      print *,'test_v03_repeated: FAIL',failures
      error stop 1
   end if
contains
   subroutine check(ok,name)
      logical,intent(in)::ok
      character(*),intent(in)::name
      if(.not.ok)then;print *,'FAIL: ',trim(name);failures=failures+1;end if
   end subroutine check
end program test_v03_repeated
