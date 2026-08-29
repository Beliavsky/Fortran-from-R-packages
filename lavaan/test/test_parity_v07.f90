program test_parity_v07
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use lavaan_kinds, only : dp
   use lavaan_robust_difference, only : robust_difference_result, satorra_bentler_difference_2001, &
      satorra_bentler_difference_2010
   use lavaan_sam_blocks, only : sam_block_covariance_result, sam_block_covariance, sam_second_order_bias
   use lavaan_ram, only : ram_model
   use lavaan_miiv_markers, only : miiv_marker
   use lavaan_miiv_partable, only : miiv_auto_markers, miiv_named_equation, ram_miiv_named_equations
   use lavaan_multilevel_random, only : random_effects_result
   use lavaan_multilevel_random_missing, only : random_coefficient_loglik_missing, random_effects_eb_missing
   use lavaan_muthen_mixed, only : muthen_mixed_result, muthen1984_mixed, lavaan_numeric, lavaan_ordered
   implicit none
   integer :: fails
   fails=0
   call test_difference(fails)
   call test_sam_blocks(fails)
   call test_auto_miiv(fails)
   call test_missing_random(fails)
   call test_mixed_muthen(fails)
   if(fails/=0) then
      write(*,'(a,i0)') 'test_parity_v07: FAIL ',fails
      error stop 1
   end if
   write(*,'(a)') 'test_parity_v07: PASS'
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

   subroutine test_difference(f)
      integer,intent(inout)::f
      type(robust_difference_result)::r
      call satorra_bentler_difference_2001(20.0_dp,10.0_dp,1.2_dp,10.0_dp,5.0_dp,1.0_dp,r)
      call check(r%status==0,'SB2001 status',f)
      call check(abs(r%scaling_factor-1.4_dp)<1.0e-12_dp,'SB2001 scale',f)
      call check(abs(r%statistic-10.0_dp/1.4_dp)<1.0e-12_dp,'SB2001 statistic',f)
      call satorra_bentler_difference_2010(20.0_dp,10.0_dp,1.2_dp,10.0_dp,5.0_dp,1.1_dp,r)
      call check(r%status==0 .and. abs(r%scaling_factor-1.3_dp)<1.0e-12_dp,'SB2010 scale',f)
   end subroutine test_difference

   subroutine test_sam_blocks(f)
      integer,intent(inout)::f
      real(dp)::v1(2,2),v2(1,1),j(1,2),curv(1,2,2),b1(2),b2(1)
      real(dp),allocatable::bout(:)
      type(sam_block_covariance_result)::r
      integer::status
      v1=reshape([1.0_dp,0.2_dp,0.2_dp,2.0_dp],[2,2])
      v2=0.5_dp
      j=reshape([0.4_dp,-0.1_dp],[1,2])
      call sam_block_covariance(v1,v2,j,r)
      call check(r%status==0,'SAM block status',f)
      call check(abs(r%cross_vcov(1,1)-0.38_dp)<1.0e-12_dp,'SAM cross covariance',f)
      call check(abs(r%structural_marginal_vcov(1,1)-0.664_dp)<1.0e-12_dp,'SAM marginal covariance',f)
      curv=0.0_dp
      curv(1,1,1)=0.2_dp
      curv(1,2,2)=-0.1_dp
      b1=[0.1_dp,-0.2_dp]
      b2=0.05_dp
      call sam_second_order_bias(b1,b2,j,curv,v1,bout,status)
      call check(status==0 .and. abs(bout(1)-0.11_dp)<1.0e-12_dp,'SAM second-order bias',f)
   end subroutine test_sam_blocks

   subroutine test_auto_miiv(f)
      integer,intent(inout)::f
      type(ram_model)::m
      type(miiv_marker),allocatable::markers(:)
      type(miiv_named_equation),allocatable::eq(:)
      character(len=8)::names(6)
      integer::status,k,hit
      names=['eta1    ','eta2    ','y1      ','y2      ','y3      ','y4      ']
      allocate(m%a(6,6),m%s(6,6),m%m(6),m%observed(4))
      m%a=0.0_dp
      m%s=0.0_dp
      m%m=0.0_dp
      m%a(2,1)=0.5_dp
      m%a(3,1)=1.0_dp
      m%a(4,2)=1.0_dp
      m%a(5,1)=0.8_dp
      m%a(6,2)=0.9_dp
      m%s(1,1)=1.0_dp
      m%s(2,2)=0.75_dp
      do k=3,6
      m%s(k,k)=0.2_dp
      end do
      m%observed=[3,4,5,6]
      call miiv_auto_markers(m,[1,2],markers,status)
      call check(status==0 .and. markers(1)%marker_node==3 .and. markers(2)%marker_node==4,'auto markers',f)
      call ram_miiv_named_equations(m,[1,2],names,eq,status)
      call check(status==0,'named MIIV status',f)
      hit=0
      do k=1,size(eq)
      if(eq(k)%outcome_node==2) hit=k
      end do
      call check(hit>0,'named MIIV structural equation',f)
      if(hit>0) then
         call check(trim(eq(hit)%proxy_outcome_name)=='y2','named MIIV proxy name',f)
         call check(any(eq(hit)%instrument_names=='y3'),'named MIIV instrument name',f)
      end if
   end subroutine test_auto_miiv

   subroutine test_missing_random(f)
      integer,intent(inout)::f
      real(dp)::y(4,1),ym(4,1),x(4,1),beta(1,1),z(4,1),g(1,1),r(1,1),ll0,ll1
      integer::cl(4)
      type(random_effects_result)::eb
      y(:,1)=[1.0_dp,-1.0_dp,0.5_dp,-0.5_dp]
      ym=y
      x=1.0_dp
      beta=0.0_dp
      z=1.0_dp
      g=0.5_dp
      r=1.0_dp
      cl=[1,1,2,2]
      ll0=random_coefficient_loglik_missing(y,cl,x,beta,z,g,r)
      call check(abs(ll0+5.618901313378636_dp)<1.0e-10_dp,'missing random complete-data identity',f)
      ym(2,1)=ieee_value(0.0_dp,ieee_quiet_nan)
      ll1=random_coefficient_loglik_missing(ym,cl,x,beta,z,g,r)
      call check(ll1>-huge(1.0_dp)/10.0_dp .and. abs(ll1-ll0)>1.0e-10_dp,'missing random finite likelihood',f)
      call random_effects_eb_missing(ym,cl,x,beta,z,g,r,eb)
      call check(eb%status==0 .and. size(eb%mean,1)==2,'missing random EB',f)
   end subroutine test_missing_random

   subroutine test_mixed_muthen(f)
      integer,intent(inout)::f
      integer,parameter::n=120
      real(dp)::dat(n,2),exo(n,1),e1,e2,xx
      integer::typ(2),nc(2),i
      type(muthen_mixed_result)::r
      typ=[lavaan_numeric,lavaan_ordered]
      nc=[0,2]
      do i=1,n
         xx=-1.5_dp+3.0_dp*real(i-1,dp)/real(n-1,dp)
         exo(i,1)=xx
         e1=0.15_dp*sin(1.7_dp*real(i,dp))+0.05_dp*cos(0.31_dp*real(i,dp))
         e2=0.85_dp*sin(0.91_dp*real(i,dp))+0.25_dp*cos(0.37_dp*real(i,dp))
         dat(i,1)=1.0_dp+0.5_dp*xx+e1
         if(0.7_dp*xx+e2>0.0_dp) then
         dat(i,2)=2.0_dp
         else
         dat(i,2)=1.0_dp
         end if
      end do
      call muthen1984_mixed(dat,typ,nc,r,exo=exo,compute_gamma=.false.)
      call check(r%status==0,'mixed Muthen status',f)
      call check(abs(r%intercept(1)-1.0_dp)<0.05_dp,'mixed numeric intercept',f)
      call check(abs(r%slopes(1,1)-0.5_dp)<0.04_dp,'mixed numeric conditional slope',f)
      call check(r%slopes(2,1)>0.35_dp,'mixed ordinal probit slope',f)
      call check(abs(r%correlation(1,2))<0.35_dp,'mixed residual correlation',f)
      call check(size(r%stats)>=6,'mixed statistic vector',f)
   end subroutine test_mixed_muthen
end program test_parity_v07
