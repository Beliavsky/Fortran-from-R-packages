program test_inference_v02
   use lavaan
   implicit none
   type(ram_model) :: model, vmodel, wmodel, bmodel
   type(ram_free_map) :: map, vmap, wmap, bmap
   type(sem_fit_result) :: fit
   type(robust_sem_result) :: rob, crob
   type(bootstrap_sem_result) :: boot
   type(twolevel_sem_result) :: tl
   real(dp) :: data(120,2), one(60,1), ml(100,1), x1, noise
   integer :: cluster(120), cl2(100), i,g,j

   call make_path_model(model,map)
   do i=1,120
      x1=sin(0.31_dp*real(i,dp))+0.5_dp*cos(0.13_dp*real(i,dp))
      noise=0.35_dp*sin(1.7_dp*real(i,dp))*(1.0_dp+0.3_dp*abs(x1))
      data(i,1)=x1
      data(i,2)=0.55_dp*x1+noise
      cluster(i)=(i-1)/5+1
   end do
   data(:,1)=data(:,1)-sum(data(:,1))/120.0_dp
   data(:,2)=data(:,2)-sum(data(:,2))/120.0_dp
   call fit_ram_data(model,map,data,fit)
   call check(fit%converged,'robust base fit')
   call robust_ml_inference(model,map,data,fit,rob)
   call check(rob%status==0 .and. all(rob%se<huge(1.0_dp)/10.0_dp),'sandwich SE')
   call check(rob%sb_scaling>0.0_dp .and. rob%chisq_scaled>=0.0_dp,'SB scaling')
   call robust_ml_inference(model,map,data,fit,crob,cluster)
   call check(crob%status==0 .and. crob%ncluster==24,'cluster sandwich')

   call make_variance_model(vmodel,vmap,1.0_dp)
   do i=1,60
   one(i,1)=0.8_dp*sin(0.7_dp*real(i,dp))+0.4_dp*cos(0.17_dp*real(i,dp))
   end do
   one(:,1)=one(:,1)-sum(one(:,1))/60.0_dp
   call bootstrap_ram_data(vmodel,vmap,one,12,boot,seed=12345)
   call check(boot%status==0 .and. boot%n_success>=8,'bootstrap success')
   call check(boot%se(1)>0.0_dp .and. boot%ci_upper(1)>boot%ci_lower(1),'bootstrap interval')

   call make_variance_model(wmodel,wmap,0.5_dp)
   call make_variance_model(bmodel,bmap,1.0_dp)
   i=0
   do g=1,20
      do j=1,5
         i=i+1
         cl2(i)=g
         ml(i,1)=0.7_dp*sin(0.6_dp*real(g,dp))+0.25_dp*real(j-3,dp)
      end do
   end do
   call fit_ram_twolevel(wmodel,wmap,bmodel,bmap,ml,cl2,tl)
   call check(tl%status==0,'two-level fit')
   call check(tl%icc(1)>0.2_dp .and. tl%icc(1)<0.95_dp,'two-level ICC')

   print '(a)', 'test_inference_v02: PASS'
contains
   subroutine check(cond,label)
      logical,intent(in)::cond
      character(len=*),intent(in)::label
      if(.not.cond) then
      write(*,'(a,1x,a)') 'FAIL:',trim(label)
      error stop 1
      end if
   end subroutine check
   subroutine make_path_model(m,map)
      type(ram_model),intent(out)::m
      type(ram_free_map),intent(out)::map
      allocate(m%a(2,2),m%s(2,2),m%observed(2))
      m%a=0.0_dp
      m%s=0.0_dp
      m%observed=[1,2]
      m%a(2,1)=0.3_dp
      m%s(1,1)=0.8_dp
      m%s(2,2)=0.3_dp
      allocate(map%matrix_id(3),map%row(3),map%col(3))
      map%matrix_id=[ram_a,ram_s,ram_s]
      map%row=[2,1,2]
      map%col=[1,1,2]
   end subroutine make_path_model
   subroutine make_variance_model(m,map,start)
      type(ram_model),intent(out)::m
      type(ram_free_map),intent(out)::map
      real(dp),intent(in)::start
      allocate(m%a(1,1),m%s(1,1),m%observed(1))
      m%a=0.0_dp
      m%s=start
      m%observed=1
      allocate(map%matrix_id(1),map%row(1),map%col(1))
      map%matrix_id=ram_s
      map%row=1
      map%col=1
   end subroutine make_variance_model
end program test_inference_v02
