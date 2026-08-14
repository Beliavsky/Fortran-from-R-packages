program test_kw_random
   use rebayes_kinds, only : dp
   use rebayes_kw
   implicit none
   integer, parameter :: n=12, ncase=100
   real(dp) :: a(n,2),w(n),d(2),p,pexact,lo,hi,mid,der0,der1,der
   type(kw_result)::r
   type(kw_control)::ctl
   integer :: c,k,seed_size
   integer,allocatable::seed(:)
   call random_seed(size=seed_size);allocate(seed(seed_size));seed=314159+37*[(k,k=1,seed_size)];call random_seed(put=seed)
   d=1.0_dp;ctl%tol=1.0e-10_dp;ctl%max_iter=30000;ctl%vertex_every=10
   do c=1,ncase
      call random_number(a);a=0.02_dp+a
      call random_number(w);w=w/sum(w)
      der0=sum(w*(a(:,1)-a(:,2))/a(:,2))
      der1=sum(w*(a(:,1)-a(:,2))/a(:,1))
      if(der0<=0.0_dp)then
         pexact=0.0_dp
      else if(der1>=0.0_dp)then
         pexact=1.0_dp
      else
         lo=0.0_dp;hi=1.0_dp
         do k=1,80
            mid=0.5_dp*(lo+hi)
            der=sum(w*(a(:,1)-a(:,2))/(mid*a(:,1)+(1.0_dp-mid)*a(:,2)))
            if(der>0.0_dp)then;lo=mid;else;hi=mid;end if
         end do
         pexact=0.5_dp*(lo+hi)
      end if
      call kw_fit(a,d,w,r,ctl)
      p=r%f(1)
      if(abs(p-pexact)>2.0e-6_dp)then
         print *,"case",c,"p",p,"exact",pexact,"gap",r%kkt_gap
         error stop "random KW mismatch"
      end if
   end do
   print *,"test_kw_random: PASS (100 cases)"
end program test_kw_random
