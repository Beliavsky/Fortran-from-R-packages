program test_parity_v02
   use lavaan
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   implicit none
   type(ram_group_spec), allocatable :: groups(:)
   type(sem_multigroup_result) :: mg
   type(ram_group_data), allocatable :: gd(:)
   type(ram_model) :: m1, m2
   type(ram_free_map) :: map1, map2, free2, cand
   type(sem_fit_result) :: cfit, pfit, wfit
   type(modification_result) :: mires
   type(ordinal_wls_result) :: ow
   real(dp) :: covg(1,1,2), meang(1,2), cov2(2,2), mean2(2)
   integer :: nobs(2), ord(80,2), i

   call make_variance_model(m1,map1,1.5_dp)
   allocate(groups(2))
   groups(1)%model=m1
   groups(1)%map=map1
   groups(2)%model=m1
   groups(2)%map=map1
   allocate(groups(1)%link(1),groups(2)%link(1))
   groups(1)%link=1
   groups(2)%link=1
   covg=0.0_dp
   covg(1,1,1)=2.0_dp
   covg(1,1,2)=4.0_dp
   meang=0.0_dp
   nobs=[200,200]
   call fit_ram_multigroup_cov(groups,covg,meang,nobs,mg)
   call check(mg%converged,'multigroup equality convergence')
   call check(abs(mg%par(1)-3.0_dp)<2.0e-3_dp,'multigroup shared variance')
   call independent_group_links(groups)
   call fit_ram_multigroup_cov(groups,covg,meang,nobs,mg)
   call check(mg%converged,'multigroup free convergence')
   call check(abs(mg%par(1)-2.0_dp)<2.0e-3_dp .and. abs(mg%par(2)-4.0_dp)<2.0e-3_dp,'multigroup free variances')
   groups(1)%link=1
   groups(2)%link=1
   allocate(gd(2))
   allocate(gd(1)%x(40,1))
   allocate(gd(2)%x(30,1))
   do i=1,40
   gd(1)%x(i,1)=merge(sqrt(2.0_dp),-sqrt(2.0_dp),mod(i,2)==0)
   end do
   do i=1,30
   gd(2)%x(i,1)=merge(sqrt(2.0_dp),-sqrt(2.0_dp),mod(i,2)==0)
   end do
   gd(2)%x(30,1)=ieee_value(0.0_dp,ieee_quiet_nan)
   call fit_ram_multigroup_fiml(groups,gd,mg)
   call check(mg%converged .and. abs(mg%par(1)-2.0_dp)<3.0e-3_dp,'multigroup FIML with missingness')

   call make_two_variance_model(m2,free2)
   cov2=0.0_dp
   cov2(1,1)=2.0_dp
   cov2(2,2)=4.0_dp
   mean2=0.0_dp
   call fit_ram_cov_constrained(m2,free2,cov2,mean2,500,equal_variances,cfit)
   call check(cfit%converged,'nonlinear constraint convergence')
   call check(abs(cfit%par(1)-cfit%par(2))<2.0e-4_dp,'equality constraint')
   call check(abs(sum(cfit%par)/2.0_dp-3.0_dp)<3.0e-3_dp,'constrained common variance')

   cov2(1,2)=0.6_dp
   cov2(2,1)=0.6_dp
   allocate(cand%matrix_id(1),cand%row(1),cand%col(1))
   cand%matrix_id=ram_s
   cand%row=2
   cand%col=1
   call modification_indices_cov(m2,free2,cand,cov2,mean2,500,mires)
   call check(mires%status==0 .and. mires%mi(1)>20.0_dp,'modification index detects omitted covariance')
   call check(abs(mires%epc(1))>0.3_dp,'modification EPC')

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
   call ordinal_wls_correlation_weights(ord,ow)
   call check(ow%status==0,'ordinal WLS status')
   call check(abs(ow%correlation(1,2)-sqrt(0.5_dp))<2.0e-3_dp,'polychoric correlation reference')
   call check(maxval(ow%dwls_weight)>0.0_dp,'ordinal DWLS weight')

   call make_correlation_model(m2,map2)
   call fit_ram_cov(m2,map2,ow%correlation,[0.0_dp,0.0_dp],80,wfit,'DWLS',wls_vd=ow%dwls_weight)
   call check(wfit%converged .and. abs(wfit%par(1)-sqrt(0.5_dp))<3.0e-3_dp,'automatic ordinal DWLS fit')

   call fit_ram_pml_ordinal(m2,map2,ord,pfit)
   call check(pfit%converged,'ordinal PML convergence')
   call check(abs(pfit%par(1)-sqrt(0.5_dp))<3.0e-3_dp,'ordinal PML correlation')

   print '(a)', 'test_parity_v02: PASS'
contains
   subroutine check(cond,label)
      logical,intent(in)::cond
      character(len=*),intent(in)::label
      if(.not.cond) then
      write(*,'(a,1x,a)') 'FAIL:',trim(label)
      error stop 1
      end if
   end subroutine check
   subroutine make_variance_model(model,map,start)
      type(ram_model),intent(out)::model
      type(ram_free_map),intent(out)::map
      real(dp),intent(in)::start
      allocate(model%a(1,1),model%s(1,1),model%observed(1))
      model%a=0.0_dp
      model%s=start
      model%observed=1
      allocate(map%matrix_id(1),map%row(1),map%col(1))
      map%matrix_id=ram_s
      map%row=1
      map%col=1
   end subroutine make_variance_model
   subroutine make_two_variance_model(model,map)
      type(ram_model),intent(out)::model
      type(ram_free_map),intent(out)::map
      allocate(model%a(2,2),model%s(2,2),model%observed(2))
      model%a=0.0_dp
      model%s=0.0_dp
      model%observed=[1,2]
      model%s(1,1)=1.5_dp
      model%s(2,2)=2.5_dp
      allocate(map%matrix_id(2),map%row(2),map%col(2))
      map%matrix_id=[ram_s,ram_s]
      map%row=[1,2]
      map%col=[1,2]
   end subroutine make_two_variance_model
   subroutine make_correlation_model(model,map)
      type(ram_model),intent(out)::model
      type(ram_free_map),intent(out)::map
      allocate(model%a(2,2),model%s(2,2),model%observed(2))
      model%a=0.0_dp
      model%s=0.0_dp
      model%observed=[1,2]
      model%s(1,1)=1.0_dp
      model%s(2,2)=1.0_dp
      model%s(1,2)=0.2_dp
      model%s(2,1)=0.2_dp
      allocate(map%matrix_id(1),map%row(1),map%col(1))
      map%matrix_id=ram_s
      map%row=2
      map%col=1
   end subroutine make_correlation_model
   subroutine equal_variances(x,equality,inequality)
      real(dp),intent(in)::x(:)
      real(dp),allocatable,intent(out)::equality(:),inequality(:)
      allocate(equality(1),inequality(0))
      equality(1)=x(1)-x(2)
   end subroutine equal_variances
end program test_parity_v02
