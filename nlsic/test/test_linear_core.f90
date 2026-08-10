program test_linear_core
   use nlsic
   implicit none
   type(lsi_result) :: r
   real(dp) :: u1(1,3),c1(1),u2(2,3),c2(2),a(2,3),b(2)
   real(dp) :: ad(5,3),bd(5),xe(3),mn(1,3),bm(5,2),xm(3,2)
   integer :: i,rank,status

   u1=0.0_dp; u1(1,1)=1.0_dp; c1=1.0_dp
   call ldp(u1,c1,r)
   call assert_true(r%succeeded(),'ldp feasible status')
   call assert_vec(r%x,[1.0_dp,0.0_dp,0.0_dp],1.0e-10_dp,'ldp solution')

   u2=0.0_dp; u2(1,1)=1.0_dp; u2(2,1)=-1.0_dp; c2=[1.0_dp,1.0_dp]
   call ldp(u2,c2,r)
   call assert_true(r%status==LSI_INFEASIBLE,'ldp infeasible status')

   a=0.0_dp; a(1,1)=1.0_dp; a(2,2)=1.0_dp; b=0.0_dp
   u1=1.0_dp; c1=1.0_dp
   call lsi_ln(a,b,r,u1,c1)
   call assert_vec(r%x,[0.0_dp,0.0_dp,1.0_dp],1.0e-9_dp,'rank deficient constrained LS')

   u2=0.0_dp; u2(1,:)=1.0_dp; u2(2,3)=-1.0_dp; c2=[1.0_dp,-0.5_dp]
   call lsi_ln(a,b,r,u2,c2)
   call assert_vec(r%x,[0.25_dp,0.25_dp,0.5_dp],1.0e-9_dp,'second constrained LS')

   do i=1,5
      ad(i,1)=real(i,dp)
      ad(i,2)=1.0_dp+0.3_dp*real(i,dp)
      ad(i,3)=ad(i,2)
   end do
   xe=[0.8_dp,-0.4_dp,1.2_dp]; bd=matmul(ad,xe)
   call ls_ln(ad,bd,r)
   call assert_close(r%x(2),r%x(3),1.0e-9_dp,'least norm duplicate columns')
   call assert_true(vecnorm(matmul(ad,r%x)-bd)<1.0e-9_dp,'least norm residual')
   bm(:,1)=bd; bm(:,2)=2.0_dp*bd
   call ls_ln_multi(ad,bm,xm,rank,status)
   call assert_true(status==0 .and. rank==2,'multiple RHS status/rank')
   call assert_true(maxval(abs(matmul(ad,xm)-bm))<1.0e-9_dp,'multiple RHS fit')
   call assert_close(xm(2,1),xm(3,1),1.0e-9_dp,'multiple RHS least norm')

   mn=reshape([0.0_dp,2.0_dp,1.0_dp],[1,3])
   call ls_ln(ad,bd,r,mnorm=mn)
   call assert_true(abs(dot_product(mn(1,:),r%x))<1.0e-9_dp,'custom least norm')

   print *, 'test_linear_core: PASS'
contains
   subroutine assert_close(x,y,tol,msg)
      real(dp),intent(in)::x,y,tol; character(*),intent(in)::msg
      if(abs(x-y)>tol) then; print *, 'FAIL: ',msg,x,y; error stop 1; end if
   end subroutine
   subroutine assert_vec(x,y,tol,msg)
      real(dp),intent(in)::x(:),y(:),tol; character(*),intent(in)::msg
      if(size(x)/=size(y) .or. maxval(abs(x-y))>tol) then
         print *, 'FAIL: ',msg; print *,x; print *,y; error stop 1
      end if
   end subroutine
   subroutine assert_true(ok,msg)
      logical,intent(in)::ok; character(*),intent(in)::msg
      if(.not.ok) then; print *, 'FAIL: ',msg; error stop 1; end if
   end subroutine
end program test_linear_core
