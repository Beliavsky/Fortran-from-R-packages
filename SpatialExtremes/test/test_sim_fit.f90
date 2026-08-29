program test_sim_fit
   use SpatialExtremes
   use r_compat, only: set_seed_int
   implicit none
   real(dp)::a(2,2),x(3000,2),invmean(2),coord(2,1),cor(2,2),xb(2),fb
   integer::it,conv
   type(maxstab_fit_t)::fit
   call set_seed_int(12345)
   a=reshape([0.7_dp,0.3_dp,0.3_dp,0.7_dp],[2,2])
   x=simulate_max_linear(3000,a)
   invmean=sum(1.0_dp/x,dim=1)/3000.0_dp
   call check(maxval(abs(invmean-1.0_dp))<0.08_dp,'max-linear unit Frechet margins')
   coord(:,1)=[0.0_dp,1.0_dp]
   cor=reshape([1.0_dp,0.3_dp,0.3_dp,1.0_dp],[2,2])
   x(1:100,:)=simulate_schlather_exact(100,coord,COV_POWEREXP,0.0_dp,1.5_dp,1.0_dp)
   call check(all(x(1:100,:)>0.0_dp),'exact Schlather positivity')
   fit=fit_schlather_frechet(x(1:100,:),coord,COV_POWEREXP)
   call check(fit%convergence==0 .or. fit%iterations>0,'fit ran')
   call check(fit%par(1)>=0.0_dp .and. fit%par(1)<1.0_dp .and. fit%par(2)>0.0_dp,'fit feasible')
   call nelder_mead(quad,[3.0_dp,-4.0_dp],xb,fb,it,conv,maxit=500,tol=1e-10_dp)
   call check(maxval(abs(xb-[1.0_dp,-2.0_dp]))<2e-5_dp,'Nelder Mead')
   print *,'test_sim_fit: PASS'
contains
   function quad(z) result(v)
      real(dp),intent(in)::z(:)
      real(dp)::v
      v=(z(1)-1.0_dp)**2+2.0_dp*(z(2)+2.0_dp)**2
   end function
   subroutine check(ok,msg)
      logical,intent(in)::ok
      character(*),intent(in)::msg
      if(.not.ok)then
      print *,'FAIL: ',msg
      error stop 1
      end if
   end subroutine
end program
