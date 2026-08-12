module matchingmarkets_plp
   use matchingmarkets_kinds, only : dp
   use matchingmarkets_types, only : plp_result_t
   use lpsolve, only : lp_result, solve_lp, LP_MAX, LP_EQ, LP_OPTIMAL
   implicit none
   private
   public :: plp
contains
   function plp(v) result(out)
      real(dp),intent(in)::v(:,:)
      type(plp_result_t)::out
      integer::n,m,i,j,k,r
      real(dp),allocatable::obj(:),a(:,:),rhs(:)
      integer,allocatable::sense(:),bin(:),ii(:),jj(:)
      type(lp_result)::sol
      n=size(v,1);if(size(v,2)/=n)error stop 'plp: valuation matrix must be square'
      if(mod(n,2)/=0)error stop 'plp: market size must be even'
      m=n*(n-1)/2;allocate(obj(m),a(n,m),rhs(n),sense(n),bin(m),ii(m),jj(m))
      a=0.0_dp;rhs=1.0_dp;sense=LP_EQ;bin=[(k,k=1,m)];k=0
      do j=2,n
         do i=1,j-1
            k=k+1;ii(k)=i;jj(k)=j;obj(k)=v(i,j)+v(j,i);a(i,k)=1.0_dp;a(j,k)=1.0_dp
         end do
      end do
      call solve_lp(LP_MAX,obj,a,sense,rhs,sol,binary_variables=bin)
      out%status=sol%status;out%objective=sol%objective
      allocate(out%assignment(n,n));out%assignment=0
      if(sol%status==LP_OPTIMAL)then
         do k=1,m;if(sol%solution(k)>0.5_dp)then;out%assignment(ii(k),jj(k))=1;out%assignment(jj(k),ii(k))=1;end if;end do
      end if
      allocate(out%pairs(2,count([(sol%solution(k)>0.5_dp,k=1,m)])))
      r=0
      do k=1,m
         if(sol%solution(k)>0.5_dp)then;r=r+1;out%pairs(:,r)=[ii(k),jj(k)];end if
      end do
   end function plp
end module matchingmarkets_plp
