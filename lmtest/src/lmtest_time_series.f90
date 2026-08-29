module lmtest_time_series
   use lmtest_kinds, only : dp
   use lmtest_types, only : test_result
   use lmtest_inference, only : nested_linear_test
   implicit none
   private
   public :: granger_test, lag_matrix

contains

   function lag_matrix(x, order) result(lags)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: order
      real(dp), allocatable :: lags(:,:)
      integer :: n, i, j
      n=size(x)
      if(order<1 .or. n<=order)then
         allocate(lags(0,0))
         return
      end if
      allocate(lags(n-order,order))
      do i=1,n-order
         do j=1,order
            lags(i,j)=x(order+i-j)
         end do
      end do
   end function lag_matrix

   function granger_test(x, y, order, use_f) result(out)
      ! Tests whether lags of x add explanatory power for y after lags of y.
      real(dp),intent(in)::x(:),y(:)
      integer,intent(in),optional::order
      logical,intent(in),optional::use_f
      type(test_result)::out
      real(dp),allocatable::lx(:,:),ly(:,:),full(:,:),reduced(:,:),yt(:)
      integer::p,n,m
      logical::as_f
      p=1
      if(present(order))p=order
      n=size(y)
      if(size(x)/=n .or. n<=p .or. p<1)return
      lx=lag_matrix(x,p)
      ly=lag_matrix(y,p)
      m=n-p
      allocate(full(m,1+2*p),reduced(m,1+p),yt(m))
      full(:,1)=1.0_dp
      full(:,2:1+p)=ly
      full(:,2+p:1+2*p)=lx
      reduced(:,1)=1.0_dp
      reduced(:,2:1+p)=ly
      yt=y(p+1:n)
      as_f=.true.
      if(present(use_f))as_f=use_f
      out=nested_linear_test(yt,full,reduced,as_f)
   end function granger_test

end module lmtest_time_series
