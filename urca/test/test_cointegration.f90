program test_cointegration
   use urca_kinds, only : dp
   use urca_types, only : po_result, johansen_result, vecm_result
   use urca_cointegration
   implicit none
   integer, parameter :: n = 240
   real(dp) :: d(n,4), x(n,3)
   integer :: i, iu
   type(po_result) :: po
   type(johansen_result) :: jo, je
   type(vecm_result) :: v

   open(newunit=iu,file='test/reference_data.txt',status='old',action='read')
   do i=1,n
      read(iu,*) d(i,:)
   end do
   close(iu)
   x=d(:,2:4)

   po=phillips_ouliaris(x(:,1:2),PO_CONST,1,PO_PU)
   call assert_zero(po%info,'PO status')
   call assert_close(po%statistic,109.9797706037426_dp,2.0e-9_dp,'PO Pu')

   jo=johansen_test(x,JO_TRACE,JO_NONE,2,JO_LONGRUN)
   call assert_zero(jo%info,'Johansen trace status')
   call assert_close(jo%lambda(1),0.3155669409245894_dp,3.0e-12_dp,'Johansen lambda1')
   call assert_close(jo%lambda(2),0.21808937961951297_dp,3.0e-12_dp,'Johansen lambda2')
   call assert_close(jo%lambda(3),0.007898095658131796_dp,3.0e-12_dp,'Johansen lambda3')
   call assert_close(jo%teststat(1),1.887209295673415_dp,2.0e-10_dp,'Johansen trace r=2')
   call assert_close(jo%teststat(2),60.43874149122232_dp,3.0e-9_dp,'Johansen trace r=1')
   call assert_close(jo%teststat(3),150.67987687768203_dp,8.0e-9_dp,'Johansen trace r=0')

   je=johansen_test(x,JO_EIGEN,JO_CONST,2,JO_TRANSITORY)
   call assert_zero(je%info,'Johansen eigen/transitory status')
   call assert_true(all(je%lambda>=0.0_dp).and.all(je%lambda<1.0_dp),'Johansen eigen range')

   v=cajools_fit(jo)
   call assert_zero(v%info,'cajools status')
   call assert_true(size(v%coefficients,1)>0.and.size(v%coefficients,2)==3,'cajools dimensions')
   v=alphaols_fit(jo,1)
   call assert_zero(v%info,'alphaols status')
   call assert_int(v%rank,1,'alphaols rank')
   v=cajorls_fit(jo,1)
   call assert_zero(v%info,'cajorls status')
   call assert_int(v%rank,1,'cajorls rank')

   print '(a)', 'test_cointegration: PASS'
contains
   subroutine assert_close(actual,expected,tol,name)
      real(dp), intent(in) :: actual, expected, tol
      character(len=*), intent(in) :: name
      if(abs(actual-expected)>tol)then
         write(*,'(a,2(1x,es24.16))') 'FAIL '//trim(name),actual,expected
         error stop 1
      end if
   end subroutine
   subroutine assert_zero(actual,name)
      integer, intent(in) :: actual
      character(len=*), intent(in) :: name
      if(actual/=0)then
         write(*,'(a,1x,i0)') 'FAIL '//trim(name),actual
         error stop 1
      end if
   end subroutine
   subroutine assert_int(actual,expected,name)
      integer, intent(in) :: actual, expected
      character(len=*), intent(in) :: name
      if(actual/=expected)then
         write(*,'(a,2(1x,i0))') 'FAIL '//trim(name),actual,expected
         error stop 1
      end if
   end subroutine
   subroutine assert_true(ok,name)
      logical, intent(in) :: ok
      character(len=*), intent(in) :: name
      if(.not.ok)then
         write(*,'(a)') 'FAIL '//trim(name)
         error stop 1
      end if
   end subroutine
end program test_cointegration
