program test_equalities_utils
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf
   use nlsic
   implicit none
   type(lsi_result) :: r
   real(dp) :: a(5,4),b(5),e(1,4),ce(1),u(1,4),co(1),mn(1,4)
   real(dp) :: atls(4,2),btls(4),xtls(2),mvec(5,1),xp(3),ap(2,3),bp(2)
   real(dp) :: lower(3),upper(3),inf
   real(dp),allocatable :: basis(:,:),ub(:,:),cb(:)
   integer :: i,rank,status

   do i=1,5
      a(i,1)=real(i,dp); a(i,2)=1.0_dp; a(i,3)=a(i,1); a(i,4)=a(i,1)
      b(i)=2.0_dp*a(i,1)+0.3_dp; mvec(i,1)=real(i,dp)
   end do
   e=0.0_dp; e(1,2)=1.0_dp; ce=0.8_dp
   u=0.0_dp; u(1,1)=-1.0_dp; co=1.0_dp
   call lsie_ln(a,b,r,u,co,e,ce)
   call assert_true(r%succeeded(),'lsie status')
   call assert_close(dot_product(e(1,:),r%x),ce(1),1.0e-9_dp,'equality enforced')
   call assert_true(dot_product(u(1,:),r%x)>=co(1)-1.0e-9_dp,'inequality enforced')

   mn=reshape([0.0_dp,0.0_dp,1.0_dp,2.0_dp],[1,4])
   call lsie_ln(a,b,r,u,co,e,ce,mnorm=mn)
   call assert_close(dot_product(mn(1,:),r%x),0.0_dp,2.0e-8_dp,'mnorm with equality')

   call nulla(mvec,basis,rank)
   call assert_true(rank==1 .and. size(basis,2)==4,'Nulla dimensions')
   call assert_true(maxval(abs(matmul(transpose(mvec),basis)))<1.0e-10_dp,'Nulla orthogonality')

   ap=0.0_dp; ap(1,1)=1.0_dp; ap(2,2)=1.0_dp; bp=[1.0_dp,2.0_dp]
   call pnull(ap,bp,xp,basis,rank,status=status)
   call assert_true(status==0,'pnull status')
   call assert_vec(xp,[1.0_dp,2.0_dp,0.0_dp],1.0e-10_dp,'pnull particular')
   call assert_true(maxval(abs(matmul(ap,basis)))<1.0e-10_dp,'pnull basis')

   atls=reshape([1.0_dp,0.0_dp,1.0_dp,2.0_dp, &
                 0.0_dp,1.0_dp,1.0_dp,-1.0_dp],[4,2])
   btls=matmul(atls,[1.5_dp,-0.75_dp])
   call tls(atls,btls,xtls,status)
   call assert_true(status==0,'tls status')
   call assert_vec(xtls,[1.5_dp,-0.75_dp],2.0e-7_dp,'tls exact relation')

   inf=ieee_value(0.0_dp,ieee_positive_inf)
   lower=[1.0_dp,-inf,0.0_dp]; upper=[inf,2.0_dp,3.0_dp]
   call uplo_to_uco(lower,upper,ub,cb)
   call assert_true(size(ub,1)==4 .and. size(ub,2)==3,'uplo dimensions')
   call assert_true(all(matmul(ub,[1.5_dp,1.5_dp,2.0_dp])>=cb),'uplo feasible point')

   print *, 'test_equalities_utils: PASS'
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
end program test_equalities_utils
