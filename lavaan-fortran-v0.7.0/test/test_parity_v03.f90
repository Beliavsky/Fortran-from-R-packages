program test_parity_v03
   use lavaan
   implicit none
   type(efa_result) :: efa,efaml
   type(categorical_stats_result) :: cat
   type(miiv_result) :: iv
   type(twolevel_ml_result) :: tl
   type(mml_ordinal_result) :: mm
   type(scaled_test_result) :: rt
   type(ram_model) :: wm,bm,sm,sm2,ivmodel,pmix
   type(ram_free_map) :: wmap,bmap,pmap
   real(dp) :: s(6,6), ltrue(6,2), psi(6), x(60,1), z(60,1), y(60)
   real(dp) :: data2(24,1), zero_load(3,1), ll0, ug(2,2), mixed(80,2), grid(10)
   integer :: ord(80,2), ord3(100,3), cl(24), i,r,info
   integer, allocatable :: inst(:)
   type(sem_fit_result) :: pmfit
   logical :: mask(3,1)

   ! EFA extraction + GPA varimax rotation on an exact two-factor covariance.
   ltrue=0.0_dp
   ltrue(1:3,1)=[0.80_dp,0.75_dp,0.70_dp]
   ltrue(4:6,2)=[0.82_dp,0.76_dp,0.68_dp]
   do i=1,6
   psi(i)=1.0_dp-sum(ltrue(i,:)**2)
   end do
   s=matmul(ltrue,transpose(ltrue))
   do i=1,6
   s(i,i)=s(i,i)+psi(i)
   end do
   call efa_fit_cov(s,2,efa,'PAF','varimax',.false.)
   call check(efa%status==0 .and. efa%converged,'EFA PAF convergence')
   call check(maxval(abs(efa%residual))<5.0e-3_dp,'EFA covariance reproduction')
   call check(abs(sum(efa%loadings(1:3,1)**2)+sum(efa%loadings(1:3,2)**2)- &
                  sum(ltrue(1:3,:)**2))<5.0e-3_dp,'EFA communalities')

   call efa_fit_cov(s,2,efaml,'ML','none',.false.)
   call check(efaml%status==0 .and. efaml%converged,'EFA ML convergence')
   call check(maxval(abs(efaml%residual))<2.0e-2_dp,'EFA ML covariance reproduction')

   ! Threshold-inclusive categorical Gamma matrix: 2 thresholds + one correlation.
   do i=1,30
   ord(i,:)=[1,1]
   end do
   do i=31,40
   ord(i,:)=[1,2]
   end do
   do i=41,50
   ord(i,:)=[2,1]
   end do
   do i=51,80
   ord(i,:)=[2,2]
   end do
   call categorical_wls_statistics(ord,cat)
   call check(cat%status==0,'categorical full WLS status')
   call check(size(cat%stats)==3 .and. all(shape(cat%gamma)==[3,3]),'categorical Gamma dimensions')
   call check(abs(cat%correlation(1,2)-sqrt(0.5_dp))<3.0e-3_dp,'categorical polychoric reference')
   call check(all(cat%dwls_weight>0.0_dp),'categorical DWLS weights')

   ! Raw-data MIIV/2SLS with a strong instrument and deterministic small disturbances.
   do i=1,60
      z(i,1)=real(i-30,dp)/10.0_dp
      x(i,1)=z(i,1)+0.15_dp*sin(real(i,dp))
      y(i)=1.0_dp+2.0_dp*x(i,1)+0.10_dp*cos(2.0_dp*real(i,dp))
   end do
   call miiv_2sls(y,x,z,iv)
   call check(iv%status==0,'MIIV 2SLS status')
   call check(abs(iv%beta(1)-1.0_dp)<0.08_dp .and. abs(iv%beta(2)-2.0_dp)<0.04_dp,'MIIV coefficient recovery')
   call check(iv%first_stage_f(1)>100.0_dp,'MIIV strong-instrument diagnostic')

   call make_iv_ram(ivmodel)
   call ram_miiv_candidates(ivmodel,2,[1],[3],inst,info)
   call check(info==0 .and. size(inst)==1 .and. inst(1)==3,'RAM model-implied IV candidate')

   ! Exact two-level Gaussian ML: one manifest variable, free within/between variances and grand mean.
   call variance_mean_models(wm,wmap,bm,bmap)
   do r=1,6
      do i=1,4
         cl((i-1)*6+r)=i
         data2((i-1)*6+r,1)=real(2*i-5,dp)*0.55_dp + real(r-3,dp)*0.18_dp
      end do
   end do
   call fit_ram_twolevel_ml(wm,wmap,bm,bmap,data2,cl,tl)
   call check(tl%status==0 .and. tl%converged,'two-level exact ML convergence')
   call check(tl%sigma_within(1,1)>0.01_dp .and. tl%sigma_between(1,1)>0.10_dp,'two-level variance components')
   call check(abs(tl%mu(1))<0.15_dp,'two-level grand mean')

   ! Ordinal marginal ML via Gauss-Hermite quadrature.
   r=0
   do i=1,30
   r=r+1
   ord3(r,:)=[1,1,1]
   end do
   do i=1,30
   r=r+1
   ord3(r,:)=[2,2,2]
   end do
   do i=1,8
   r=r+1
   ord3(r,:)=[1,1,2]
   end do
   do i=1,8
   r=r+1
   ord3(r,:)=[2,2,1]
   end do
   do i=1,6
   r=r+1
   ord3(r,:)=[1,2,1]
   end do
   do i=1,6
   r=r+1
   ord3(r,:)=[2,1,2]
   end do
   do i=1,6
   r=r+1
   ord3(r,:)=[1,2,2]
   end do
   do i=1,6
   r=r+1
   ord3(r,:)=[2,1,1]
   end do
   zero_load=0.0_dp
   mask=.true.
   call fit_mml_ordinal_factor(ord3,reshape([0.45_dp,0.45_dp,0.45_dp],[3,1]),mask,mm,7)
   ll0=mml_ordinal_loglik(ord3,zero_load,mm%thresholds,mm%ncat,7)
   call check(mm%status==0 .and. mm%converged,'ordinal MML convergence')
   call check(mm%loglik>ll0+2.0_dp,'ordinal MML improves independence')
   call check(minval(abs(mm%loadings(:,1)))>0.2_dp,'ordinal MML nonzero loadings')


   grid=[-1.5_dp,-1.2_dp,-0.9_dp,-0.6_dp,-0.3_dp,0.3_dp,0.6_dp,0.9_dp,1.2_dp,1.5_dp]
   do i=1,80
      mixed(i,1)=grid(mod(i-1,10)+1)
      mixed(i,2)=merge(1.0_dp,2.0_dp,mod(i-1,10)<5)
      if(mod(i,17)==0) mixed(i,2)=3.0_dp-mixed(i,2)
   end do
   call make_mixed_model(pmix,pmap)
   call fit_ram_pml_mixed(pmix,pmap,mixed,[.false.,.true.],pmfit)
   call check(pmfit%status==0 .and. pmfit%converged,'mixed continuous-ordinal PML convergence')
   call check(pmfit%par(1)>0.3_dp,'mixed PML positive latent correlation')

   ! Mean and mean/variance adjusted robust test formulas.
   ug=0.0_dp
   ug(1,1)=2.0_dp
   ug(2,2)=1.0_dp
   call scaled_tests_from_ugamma(10.0_dp,2.0_dp,ug,rt)
   call check(rt%status==0,'robust scaled tests status')
   call check(abs(rt%sb_scaling-1.5_dp)<1e-12_dp .and. abs(rt%chisq_sb-20.0_dp/3.0_dp)<1e-12_dp, &
              'Satorra-Bentler formula')
   call check(abs(rt%df_mv-1.8_dp)<1e-12_dp .and. abs(rt%chisq_mv-6.0_dp)<1e-12_dp,'mean-variance adjusted formula')

   ! SAM stage-freezing helper.
   sm=bm
   call sam_fix_measurement(sm,bmap,[1.25_dp,0.30_dp],sm2)
   call check(abs(sm2%s(1,1)-1.25_dp)<1e-12_dp .and. abs(sm2%m(1)-0.30_dp)<1e-12_dp,'SAM parameter freezing')

   print '(a)', 'test_parity_v03: PASS'
contains
   subroutine check(cond,label)
      logical,intent(in)::cond
      character(len=*),intent(in)::label
      if(.not.cond) then
      write(*,'(a,1x,a)') 'FAIL:',trim(label)
      error stop 1
      end if
   end subroutine check

   subroutine variance_mean_models(wmod,wmp,bmod,bmp)
      type(ram_model),intent(out)::wmod,bmod
      type(ram_free_map),intent(out)::wmp,bmp
      allocate(wmod%a(1,1),wmod%s(1,1),wmod%observed(1))
      wmod%a=0.0_dp
      wmod%s=0.12_dp
      wmod%observed=1
      allocate(wmp%matrix_id(1),wmp%row(1),wmp%col(1))
      wmp%matrix_id=ram_s
      wmp%row=1
      wmp%col=1
      allocate(bmod%a(1,1),bmod%s(1,1),bmod%m(1),bmod%observed(1))
      bmod%a=0.0_dp
      bmod%s=0.6_dp
      bmod%m=0.0_dp
      bmod%observed=1
      allocate(bmp%matrix_id(2),bmp%row(2),bmp%col(2))
      bmp%matrix_id=[ram_s,ram_m]
      bmp%row=[1,1]
      bmp%col=[1,1]
   end subroutine variance_mean_models

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

   subroutine make_mixed_model(model,map)
      type(ram_model),intent(out)::model
      type(ram_free_map),intent(out)::map
      allocate(model%a(2,2),model%s(2,2),model%m(2),model%observed(2))
      model%a=0.0_dp
      model%s=0.0_dp
      model%s(1,1)=1.0_dp
      model%s(2,2)=1.0_dp
      model%s(1,2)=0.2_dp
      model%s(2,1)=0.2_dp
      model%m=0.0_dp
      model%observed=[1,2]
      allocate(map%matrix_id(1),map%row(1),map%col(1))
      map%matrix_id=ram_s
      map%row=2
      map%col=1
   end subroutine make_mixed_model
end program test_parity_v03
