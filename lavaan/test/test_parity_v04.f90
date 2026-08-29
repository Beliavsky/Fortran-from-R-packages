program test_parity_v04
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use lavaan
   implicit none
   type(ram_model) :: mt,st,ivm,wm,bm
   type(ram_free_map) :: mmap,smap,wmap,bmap
   type(sam_result) :: sam
   type(miiv_equation),allocatable :: eqs(:)
   type(categorical_stats_result) :: ca,cj
   type(twolevel_ml_result) :: tl
   type(mml_mixed_result) :: mm
   type(scaled_test_result) :: rt
   real(dp) :: scov(2,2),smean(2),data2(20,2),mix(80,2),lstart(2,1),lzero(2,1),eta_cov(1,1),eta_mean(1)
   real(dp) :: ll0,zgrid(10),nanv
   integer :: ord(80,2),cl(20),i,r,info,idx
   logical :: free(2,1),ordinal(2)

   call make_sam_models(mt,mmap,st,smap)
   scov=reshape([1.0_dp,0.6_dp,0.6_dp,1.0_dp],[2,2])
   smean=0.0_dp
   call sam_fit_cov(mt,mmap,st,smap,scov,smean,300,sam)
   call check(sam%status==0 .and. sam%uncertainty_propagated,'SAM propagated uncertainty status')
   call check(size(sam%stage1_jacobian,1)==1 .and. size(sam%stage1_jacobian,2)==1,'SAM Jacobian dimensions')
   call check(sam%structural_vcov_corrected(1,1)>=sam%structural%vcov(1,1)-1.0e-12_dp,'SAM covariance inflation')

   call make_iv_ram(ivm)
   call ram_miiv_equations(ivm,eqs,info)
   call check(info==0 .and. size(eqs)==2,'MIIV equation discovery count')
   idx=0
   do i=1,size(eqs)
   if(eqs(i)%outcome_node==2) idx=i
   end do
   call check(idx>0 .and. eqs(idx)%identified,'MIIV identified equation')
   call check(any(eqs(idx)%instrument_nodes==3),'MIIV valid instrument discovery')

   call make_binary_table(ord)
   call categorical_wls_statistics_analytic(ord,ca)
   call categorical_wls_statistics(ord,cj)
   call check(ca%status==0 .and. cj%status==0,'categorical Gamma status')
   call check(maxval(abs(ca%stats-cj%stats))<1.0e-10_dp,'categorical statistic agreement')
   call check(abs(ca%gamma(1,1)-0.5_dp*acos(-1.0_dp))<0.08_dp,'analytic binary-threshold variance')
   call check(minval([(ca%gamma(i,i),i=1,size(ca%gamma,1))])>0.0_dp,'analytic Gamma positive diagonal')
   call check(abs(ca%gamma(1,1)-cj%gamma(1,1))<0.15_dp,'analytic versus jackknife threshold Gamma')

   call make_twolevel_models(wm,wmap,bm,bmap)
   nanv=ieee_value(0.0_dp,ieee_quiet_nan)
   r=0
   do i=1,5
      call add_cluster(i,r,data2,cl)
   end do
   data2(3,2)=nanv
   data2(8,1)=nanv
   data2(17,2)=nanv
   call fit_ram_twolevel_fiml(wm,wmap,bm,bmap,data2,cl,tl,.true.)
   call check(tl%status==0 .and. all(abs(tl%par)<huge(1.0_dp)/10.0_dp),'two-level missing FIML fit')
   call check(tl%h1_loglik>=tl%loglik-1.0e-7_dp,'two-level H1 dominance')
   call check(abs(tl%df-2.0_dp)<1.0e-12_dp .and. tl%chisq>=0.0_dp,'two-level H1 likelihood-ratio df')

   zgrid=[-1.5_dp,-1.2_dp,-0.9_dp,-0.6_dp,-0.3_dp,0.3_dp,0.6_dp,0.9_dp,1.2_dp,1.5_dp]
   do i=1,80
      r=mod(i-1,10)+1
      mix(i,1)=1.0_dp+0.8_dp*zgrid(r)+0.05_dp*sin(real(i,dp))
      mix(i,2)=merge(1.0_dp,2.0_dp,zgrid(r)+0.08_dp*cos(real(i,dp))<0.0_dp)
   end do
   mix(7,1)=nanv
   mix(23,2)=nanv
   ordinal=[.false.,.true.]
   lstart=reshape([0.55_dp,0.55_dp],[2,1])
   free=.true.
   eta_mean=0.0_dp
   eta_cov=1.0_dp
   call fit_mml_mixed_factor(mix,ordinal,lstart,free,eta_mean,eta_cov,mm,7)
   lzero=0.0_dp
   ll0=mml_mixed_loglik(mix,ordinal,lzero,mm%intercept,mm%residual_sd,mm%thresholds,mm%ncat,eta_mean,eta_cov,7)
   if(.not.mm%converged) print '(a,i0,a,es12.4)', 'mixed MML did not converge after ', mm%iterations, &
      ' iterations; log likelihood = ', mm%loglik
   call check(mm%status==0 .and. mm%converged,'mixed MML convergence')
   call check(mm%loglik>ll0+8.0_dp,'mixed MML improves independence')
   call check(minval(abs(mm%loadings(:,1)))>0.15_dp,'mixed MML loading recovery')

   call scaled_tests_from_ugamma(10.0_dp,2.0_dp,reshape([2.0_dp,0.0_dp,0.0_dp,1.0_dp],[2,2]),rt)
   call check(rt%status==0,'extended robust-test status')
   call check(abs(rt%ss_scaling-sqrt(2.5_dp))<1.0e-12_dp,'scaled-shifted scaling')
   call check(abs(rt%ss_shift-(2.0_dp-3.0_dp*sqrt(0.4_dp)))<1.0e-12_dp,'scaled-shifted shift')
   call yuan_bentler_from_traces(10.0_dp,2.0_dp,3.0_dp,8.0_dp,4.0_dp,3.0_dp,2.0_dp,rt)
   call check(abs(rt%yb_scaling-1.5_dp)<1.0e-12_dp .and. abs(rt%yb_h1_scaling-2.0_dp)<1.0e-12_dp, &
      'Yuan-Bentler trace scaling')

   print '(a)', 'test_parity_v04: PASS'
contains
   subroutine check(cond,label)
      logical,intent(in)::cond
      character(len=*),intent(in)::label
      if(.not.cond) then
      write(*,'(a,1x,a)') 'FAIL:',trim(label)
      error stop 1
      end if
   end subroutine check

   subroutine make_sam_models(meas,mmap,struct,smap)
      type(ram_model),intent(out)::meas,struct
      type(ram_free_map),intent(out)::mmap,smap
      allocate(meas%a(2,2),meas%s(2,2),meas%observed(2))
      meas%a=0.0_dp
      meas%a(2,1)=0.5_dp
      meas%s=0.0_dp
      meas%s(1,1)=1.0_dp
      meas%s(2,2)=0.64_dp
      meas%observed=[1,2]
      allocate(mmap%matrix_id(1),mmap%row(1),mmap%col(1))
      mmap%matrix_id=ram_a
      mmap%row=2
      mmap%col=1
      allocate(struct%a(2,2),struct%s(2,2),struct%observed(2))
      struct%a=0.0_dp
      struct%a(2,1)=0.5_dp
      struct%s=0.0_dp
      struct%s(1,1)=0.9_dp
      struct%s(2,2)=0.64_dp
      struct%observed=[1,2]
      allocate(smap%matrix_id(1),smap%row(1),smap%col(1))
      smap%matrix_id=ram_s
      smap%row=1
      smap%col=1
   end subroutine make_sam_models

   subroutine make_iv_ram(model)
      type(ram_model),intent(out)::model
      allocate(model%a(3,3),model%s(3,3),model%observed(3))
      model%a=0.0_dp
      model%s=0.0_dp
      model%a(1,3)=0.8_dp
      model%a(2,1)=2.0_dp
      model%s(1,1)=1.0_dp
      model%s(2,2)=1.0_dp
      model%s(3,3)=1.0_dp
      model%observed=[1,2,3]
   end subroutine make_iv_ram

   subroutine make_binary_table(x)
      integer,intent(out)::x(:,:)
      integer::a
      do a=1,30
      x(a,:)=[1,1]
      end do
      do a=31,40
      x(a,:)=[1,2]
      end do
      do a=41,50
      x(a,:)=[2,1]
      end do
      do a=51,80
      x(a,:)=[2,2]
      end do
   end subroutine make_binary_table

   subroutine make_twolevel_models(wmod,wmp,bmod,bmp)
      type(ram_model),intent(out)::wmod,bmod
      type(ram_free_map),intent(out)::wmp,bmp
      allocate(wmod%a(2,2),wmod%s(2,2),wmod%observed(2))
      wmod%a=0.0_dp
      wmod%s=0.0_dp
      wmod%s(1,1)=0.15_dp
      wmod%s(2,2)=0.20_dp
      wmod%observed=[1,2]
      allocate(wmp%matrix_id(2),wmp%row(2),wmp%col(2))
      wmp%matrix_id=[ram_s,ram_s]
      wmp%row=[1,2]
      wmp%col=[1,2]
      allocate(bmod%a(2,2),bmod%s(2,2),bmod%m(2),bmod%observed(2))
      bmod%a=0.0_dp
      bmod%s=0.0_dp
      bmod%s(1,1)=0.45_dp
      bmod%s(2,2)=0.50_dp
      bmod%m=0.0_dp
      bmod%observed=[1,2]
      allocate(bmp%matrix_id(4),bmp%row(4),bmp%col(4))
      bmp%matrix_id=[ram_s,ram_s,ram_m,ram_m]
      bmp%row=[1,2,1,2]
      bmp%col=[1,2,1,1]
   end subroutine make_twolevel_models

   subroutine add_cluster(g,pos,x,group)
      integer,intent(in)::g
      integer,intent(inout)::pos
      real(dp),intent(inout)::x(:,:)
      integer,intent(inout)::group(:)
      integer::a
      real(dp)::b1,b2,w
      b1=0.55_dp*real(g-3,dp)
      b2=0.75_dp*real(g-3,dp)
      do a=1,4
         pos=pos+1
         group(pos)=g
         w=real(a-2,dp)*0.16_dp
         x(pos,1)=b1+w
         x(pos,2)=b2+0.7_dp*w+0.05_dp*real(mod(a,2),dp)
      end do
   end subroutine add_cluster
end program test_parity_v04
