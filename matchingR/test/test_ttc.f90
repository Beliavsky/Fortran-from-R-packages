program test_ttc
   use matchingr
   implicit none
   integer :: pref(4,4)
   integer, allocatable :: m(:)
   pref=reshape([2,3,4,1, 4,3,2,1, 3,4,2,1, 4,2,1,3],[4,4])
   m=top_trading_cycles_preferences(pref)
   if(.not.top_trading_stable(pref,m)) error stop "TTC unstable"
   if(any(sortcopy(m) /= [1,2,3,4])) error stop "TTC not permutation"
   print *, "test_ttc: PASS"
contains
   function sortcopy(x) result(y)
      integer,intent(in)::x(:)
      integer::y(size(x)),i,j,t
      y=x
      do i=2,size(y)
         t=y(i);j=i-1
         do while(j>=1)
            if(y(j)<=t) exit
            y(j+1)=y(j);j=j-1
         end do
         y(j+1)=t
      end do
   end function
end program
