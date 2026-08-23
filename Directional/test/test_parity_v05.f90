program test_parity_v05
   use directional
   implicit none
   real(dp) :: y3(8,3), y2(8,2), x(8,2), a
   integer :: grp(8), fails, i
   type(directional_regression_result) :: r
   type(test_result) :: tr
   fails=0
   x(:,1)=1.0_dp
   x(:,2)=[-1.0_dp,-0.7_dp,-0.4_dp,-0.1_dp,0.1_dp,0.4_dp,0.7_dp,1.0_dp]
   do i=1,8
      a=0.35_dp*x(i,2)
      y3(i,:)=[cos(a)*cos(0.25_dp),sin(a)*cos(0.25_dp),sin(0.25_dp)]
      y2(i,:)=[cos(a),sin(a)]
   end do
   grp=[1,1,1,1,2,2,2,2]

   r=iag_reg(y3,x,maxit=120); call check_reg(r,8,3,fails)
   r=sipc_reg(y3,x,maxit=120); call check_reg(r,8,3,fails)
   r=spcauchy_reg(y3,x,maxit=120); call check_reg(r,8,3,fails)
   r=pkbd_reg(y3,x,maxit=120); call check_reg(r,8,3,fails)
   r=spml_reg(y2,x,maxit=120); call check_reg(r,8,2,fails)
   r=cipc_reg(y2,x,maxit=120); call check_reg(r,8,2,fails)
   r=gcpc_reg(y2,x,maxit=120); call check_reg(r,8,2,fails)
   r=esag_reg(y3,x(:,1:1),maxit=80); call check_reg(r,8,3,fails)
   r=sespc_reg(y3,x(:,1:1),maxit=80); call check_reg(r,8,3,fails)
   r=spher_reg(y3,y3); call check_reg(r,8,3,fails)
   if(r%fit<7.5_dp) fails=fails+1

   tr=embed_aov(y3,grp,mc_reps=19); call check_test(tr,fails)
   tr=hcf_aov(y3,grp,mc_reps=19); call check_test(tr,fails)
   tr=het_aov(y3,grp,mc_reps=19); call check_test(tr,fails)
   tr=lr_aov(y3,grp,mc_reps=19); call check_test(tr,fails)
   tr=embed_boot(y3(1:4,:),y3(5:8,:),b=19); call check_test(tr,fails)
   tr=hcf_boot(y3(1:4,:),y3(5:8,:),b=19); call check_test(tr,fails)
   tr=het_boot(y3(1:4,:),y3(5:8,:),b=19); call check_test(tr,fails)

   if(fails==0)then
      print *, 'test_parity_v05: PASS'
   else
      print *, 'test_parity_v05: FAIL',fails
      error stop 1
   end if
contains
   subroutine check_reg(q,n,d,f)
      type(directional_regression_result),intent(in)::q
      integer,intent(in)::n,d
      integer,intent(inout)::f
      if(.not.allocated(q%beta) .or. .not.allocated(q%fitted))then;f=f+1;return;end if
      if(size(q%fitted,1)/=n .or. size(q%fitted,2)/=d) f=f+1
      if(q%loglik/=q%loglik) f=f+1
      if(any(q%fitted/=q%fitted)) f=f+1
   end subroutine
   subroutine check_test(q,f)
      type(test_result),intent(in)::q;integer,intent(inout)::f
      if(q%statistic/=q%statistic) f=f+1
      if(q%p_value<0.0_dp .or. q%p_value>1.0_dp) f=f+1
   end subroutine
end program
