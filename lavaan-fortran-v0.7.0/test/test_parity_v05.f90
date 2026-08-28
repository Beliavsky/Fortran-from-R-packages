program test_parity_v05
   use lavaan_kinds, only : dp
   use lavaan_linalg, only : vech, vech_reverse, trace_matrix
   use lavaan_hayakawa, only : hayakawa_trace_result, hayakawa_trace_corrected, hayakawa_adjusted_tests
   use lavaan_miiv_variance, only : miiv_estimate_uls, miiv_estimate_gls, miiv_estimate_2rls, miiv_estimate_rls
   use lavaan_miiv_variance, only : miiv_jacobian_uls, miiv_jacobian_gls, miiv_jacobian_2rls, miiv_jacobian_rls
   use lavaan_miiv_variance, only : miiv_vcov_from_gamma
   use lavaan_sam, only : sam_yuan_chan_test
   use lavaan_robust_tests, only : scaled_test_result
   use lavaan_mml_qmc, only : mml_mixed_loglik_qmc
   implicit none
   integer :: fails
   fails=0
   call test_hayakawa(fails)
   call test_miiv_saturated(fails)
   call test_miiv_jacobians(fails)
   call test_sam_yuan_chan(fails)
   call test_mml_qmc(fails)
   if(fails/=0) then
      write(*,'(a,i0)') 'test_parity_v05: FAIL ',fails
      error stop 1
   end if
   write(*,'(a)') 'test_parity_v05: PASS'
contains
   subroutine check(ok,msg,f)
      logical,intent(in)::ok
      character(len=*),intent(in)::msg
      integer,intent(inout)::f
      if(.not.ok) then
      f=f+1
      write(*,'(a)') 'FAIL: '//trim(msg)
      end if
   end subroutine check

   subroutine test_hayakawa(f)
      integer,intent(inout)::f
      real(dp)::x(6,2),u(3,3),mean(2),yc(6,2),z(6,3),zm(3),omega(3,3),out(2,2),v(3)
      real(dp)::manual1
      integer::i
      type(hayakawa_trace_result)::r
      x=reshape([ -1.2_dp,0.2_dp, 0.3_dp,-0.7_dp, 1.1_dp,0.5_dp, &
                   0.8_dp,1.4_dp,-0.4_dp,-1.1_dp, 1.5_dp,0.9_dp ],[6,2])
      u=0.0_dp
      u(1,1)=1.2_dp
      u(2,2)=0.8_dp
      u(3,3)=1.5_dp
      u(1,2)=0.1_dp
      u(2,1)=0.1_dp
      mean=sum(x,dim=1)/6.0_dp
      do i=1,6
         yc(i,:)=x(i,:)-mean
         out=spread(yc(i,:),2,2)*spread(yc(i,:),1,2)
         z(i,:)=vech(out)
      end do
      zm=sum(z,dim=1)/6.0_dp
      omega=0.0_dp
      do i=1,6
         v=z(i,:)-zm
         omega=omega+spread(v,2,3)*spread(v,1,3)
      end do
      omega=omega/5.0_dp
      manual1=trace_matrix(matmul(u,omega))
      call hayakawa_trace_corrected(u,x,r)
      call check(r%status==0,'Hayakawa status',f)
      call check(abs(r%trace_ugamma-manual1)<1.0e-12_dp,'Hayakawa unbiased first trace',f)
      call check(r%trace_ugamma2>0.0_dp,'Hayakawa corrected second trace positive',f)
      call hayakawa_adjusted_tests(u,x,12.0_dp,3.0_dp,r)
      call check(abs(r%mv_df-r%trace_ugamma**2/r%trace_ugamma2)<1.0e-12_dp,'Hayakawa corrected df',f)
      call check(abs(r%chisq_mv-12.0_dp/r%mv_scaling)<1.0e-12_dp,'Hayakawa corrected statistic',f)
   end subroutine test_hayakawa

   subroutine test_miiv_saturated(f)
      integer,intent(inout)::f
      real(dp)::s(2,2),d(3,3),eye(3,3),gamma(3,3)
      real(dp),allocatable::theta(:),j(:,:),vc(:,:),sv(:)
      integer::info,i
      s=reshape([1.4_dp,0.3_dp,0.3_dp,0.9_dp],[2,2])
      d=0.0_dp
      eye=0.0_dp
      gamma=0.0_dp
      do i=1,3
      d(i,i)=1.0_dp
      eye(i,i)=1.0_dp
      gamma(i,i)=real(i,dp)
      end do
      sv=vech(s)
      call miiv_estimate_uls(s,d,theta,info)
      call check(info==0 .and. maxval(abs(theta-sv))<1.0e-12_dp,'MIIV ULS saturated',f)
      call miiv_estimate_gls(s,d,theta,info)
      call check(info==0 .and. maxval(abs(theta-sv))<1.0e-10_dp,'MIIV GLS saturated',f)
      call miiv_estimate_2rls(s,d,theta,info)
      call check(info==0 .and. maxval(abs(theta-sv))<1.0e-10_dp,'MIIV 2RLS saturated',f)
      call miiv_estimate_rls(s,d,theta,info)
      call check(info==0 .and. maxval(abs(theta-sv))<1.0e-10_dp,'MIIV RLS saturated',f)
      call miiv_jacobian_uls(s,d,j,info)
      call check(info==0 .and. maxval(abs(j-eye))<1.0e-12_dp,'MIIV ULS Jacobian',f)
      call miiv_jacobian_gls(s,d,j,info)
      call check(info==0 .and. maxval(abs(j-eye))<1.0e-9_dp,'MIIV GLS Jacobian',f)
      call miiv_jacobian_2rls(s,d,j,info)
      call check(info==0 .and. maxval(abs(j-eye))<1.0e-9_dp,'MIIV 2RLS Jacobian',f)
      call miiv_jacobian_rls(s,d,j,info)
      call check(info==0 .and. maxval(abs(j-eye))<1.0e-9_dp,'MIIV RLS Jacobian',f)
      call miiv_vcov_from_gamma(j,gamma,100,vc,info)
      call check(info==0 .and. maxval(abs(vc-gamma/100.0_dp))<1.0e-12_dp,'MIIV covariance from Gamma',f)
   end subroutine test_miiv_saturated

   subroutine test_miiv_jacobians(f)
      integer,intent(inout)::f
      real(dp)::s(2,2),d(3,2),h,jfd(2,3),sp(2,2),sm(2,2)
      real(dp),allocatable::j(:,:),tp(:),tm(:),sv(:)
      integer::info,k
      s=reshape([1.5_dp,0.4_dp,0.4_dp,1.2_dp],[2,2])
      d=reshape([1.0_dp,0.2_dp,0.0_dp, 0.0_dp,0.1_dp,1.0_dp],[3,2])
      sv=vech(s)
      h=1.0e-5_dp
      call miiv_jacobian_gls(s,d,j,info)
      call check(info==0,'MIIV GLS Jacobian status',f)
      do k=1,3
         sp=vech_reverse(sv+unit3(k,h),2)
         sm=vech_reverse(sv-unit3(k,h),2)
         call miiv_estimate_gls(sp,d,tp,info)
         call miiv_estimate_gls(sm,d,tm,info)
         jfd(:,k)=(tp-tm)/(2.0_dp*h)
      end do
      call check(maxval(abs(j-jfd))<3.0e-4_dp,'MIIV GLS analytic linearization',f)
      call miiv_jacobian_2rls(s,d,j,info)
      do k=1,3
         sp=vech_reverse(sv+unit3(k,h),2)
         sm=vech_reverse(sv-unit3(k,h),2)
         call miiv_estimate_2rls(sp,d,tp,info)
         call miiv_estimate_2rls(sm,d,tm,info)
         jfd(:,k)=(tp-tm)/(2.0_dp*h)
      end do
      call check(maxval(abs(j-jfd))<5.0e-4_dp,'MIIV 2RLS analytic linearization',f)
      call miiv_jacobian_rls(s,d,j,info)
      do k=1,3
         sp=vech_reverse(sv+unit3(k,h),2)
         sm=vech_reverse(sv-unit3(k,h),2)
         call miiv_estimate_rls(sp,d,tp,info)
         call miiv_estimate_rls(sm,d,tm,info)
         jfd(:,k)=(tp-tm)/(2.0_dp*h)
      end do
      call check(maxval(abs(j-jfd))<8.0e-4_dp,'MIIV RLS analytic linearization',f)
   end subroutine test_miiv_jacobians

   pure function unit3(k,a) result(v)
      integer,intent(in)::k
      real(dp),intent(in)::a
      real(dp)::v(3)
      v=0.0_dp
      v(k)=a
   end function unit3

   subroutine test_sam_yuan_chan(f)
      integer,intent(inout)::f
      real(dp)::dt(2,1),dg(2,1),pj(1,2),gamma(2,2),w(2,2)
      type(scaled_test_result)::r
      integer::status
      dt(:,1)=[1.0_dp,0.0_dp]
      dg(:,1)=[0.0_dp,1.0_dp]
      pj(1,:)=[0.0_dp,0.5_dp]
      gamma=0.0_dp
      w=0.0_dp
      gamma(1,1)=1.0_dp
      gamma(2,2)=1.0_dp
      w(1,1)=1.0_dp
      w(2,2)=1.0_dp
      call sam_yuan_chan_test(10.0_dp,1.0_dp,dt,dg,pj,gamma,w,r,status)
      call check(status==0,'SAM Yuan-Chan status',f)
      call check(abs(r%yb_scaling-0.25_dp)<1.0e-12_dp,'SAM Yuan-Chan scaling',f)
      call check(abs(r%chisq_yb-40.0_dp)<1.0e-10_dp,'SAM Yuan-Chan statistic',f)
   end subroutine test_sam_yuan_chan

   subroutine test_mml_qmc(f)
      integer,intent(inout)::f
      real(dp)::data(1,4),load(4,4),intercept(4),rsd(4),thr(1,4),lmean(4),lcov(4,4),ll,exact
      logical::ord(4)
      integer::ncat(4),i
      data=0.0_dp
      load=0.0_dp
      intercept=0.0_dp
      rsd=1.0_dp
      thr=0.0_dp
      lmean=0.0_dp
      lcov=0.0_dp
      ord=.false.
      ncat=0
      do i=1,4
      load(i,i)=1.0_dp
      lcov(i,i)=1.0_dp
      end do
      ll=mml_mixed_loglik_qmc(data,ord,load,intercept,rsd,thr,ncat,lmean,lcov,4096)
      exact=-2.0_dp*log(4.0_dp*acos(-1.0_dp))
      call check(abs(ll-exact)<4.0e-2_dp,'4-factor QMC MML likelihood',f)
   end subroutine test_mml_qmc
end program test_parity_v05
