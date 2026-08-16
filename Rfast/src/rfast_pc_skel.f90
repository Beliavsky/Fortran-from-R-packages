module rfast_pc_skel
   use rfast_special, only : dp, student_t_cdf
   use rfast_arrays, only : rank_average
   use rfast_linalg, only : correlation_matrix, solve_linear
   implicit none
   private

   integer, parameter, public :: PC_PEARSON=1, PC_SPEARMAN=2

   type, public :: pc_skeleton_result
      integer, allocatable :: graph(:,:)
      real(dp), allocatable :: statistic(:,:), pvalue(:,:)
      integer, allocatable :: sep_size(:,:), sepset(:,:,:)
      integer, allocatable :: tests_by_order(:)
      integer :: kappa = 0
      integer :: status = 0
   end type pc_skeleton_result

   public :: pc_skeleton, partial_correlation_test

contains

   function pc_skeleton(data,method,alpha,max_order) result(res)
      real(dp),intent(in)::data(:,:)
      integer,intent(in),optional::method,max_order
      real(dp),intent(in),optional::alpha
      type(pc_skeleton_result)::res
      real(dp),allocatable::z(:,:),cor(:,:)
      integer,allocatable::neighbors(:),comb(:)
      real(dp)::sig,stat,pv
      integer::meth,n,p,i,j,k,m,mx,nt,df
      logical::changed,indep
      if(size(data,1)<4.or.size(data,2)<2)then;res%status=1;return;end if
      n=size(data,1);p=size(data,2);meth=PC_PEARSON;if(present(method))meth=method
      if(meth/=PC_PEARSON.and.meth/=PC_SPEARMAN)then;res%status=2;return;end if
      sig=0.01_dp;if(present(alpha))sig=alpha
      if(sig<=0.0_dp.or.sig>=1.0_dp)then;res%status=3;return;end if
      mx=p-2;if(present(max_order))mx=min(mx,max(0,max_order));mx=min(mx,n-4)
      allocate(z(n,p));z=data
      if(meth==PC_SPEARMAN)then
         do j=1,p;z(:,j)=rank_average(z(:,j));end do
      end if
      cor=correlation_matrix(z)
      allocate(res%graph(p,p),res%statistic(p,p),res%pvalue(p,p),res%sep_size(p,p),res%sepset(p,p,max(1,mx)))
      allocate(res%tests_by_order(0:mx));res%graph=0;res%statistic=0.0_dp;res%pvalue=1.0_dp
      res%sep_size=0;res%sepset=0;res%tests_by_order=0
      do i=1,p-1
         do j=i+1,p
            call partial_correlation_test(cor,n,i,j,[integer::],meth,stat,pv,df)
            res%tests_by_order(0)=res%tests_by_order(0)+1;res%statistic(i,j)=stat;res%statistic(j,i)=stat
            res%pvalue(i,j)=pv;res%pvalue(j,i)=pv
            if(pv<=sig)then;res%graph(i,j)=1;res%graph(j,i)=1;end if
         end do
      end do
      do k=1,mx
         changed=.false.;nt=0
         do i=1,p-1
            do j=i+1,p
               if(res%graph(i,j)==0)cycle
               neighbors=neighbor_union(res%graph,i,j)
               if(size(neighbors)<k)cycle
               allocate(comb(k));comb=[(m,m=1,k)];indep=.false.
               do
                  call partial_correlation_test(cor,n,i,j,neighbors(comb),meth,stat,pv,df);nt=nt+1
                  if(pv>res%pvalue(i,j))then;res%pvalue(i,j)=pv;res%pvalue(j,i)=pv;end if
                  if(pv>sig)then
                     res%graph(i,j)=0;res%graph(j,i)=0;res%statistic(i,j)=stat;res%statistic(j,i)=stat
                     res%sep_size(i,j)=k;res%sep_size(j,i)=k;res%sepset(i,j,1:k)=neighbors(comb)
                     res%sepset(j,i,1:k)=neighbors(comb);indep=.true.;changed=.true.;exit
                  end if
                  if(.not.next_combination(comb,size(neighbors)))exit
               end do
               deallocate(comb)
               if(indep)cycle
            end do
         end do
         res%tests_by_order(k)=nt;res%kappa=k
         if(nt==0)exit
         if(.not.any_degree_at_least(res%graph,k+1))exit
      end do
      do i=1,p;res%graph(i,i)=0;res%pvalue(i,i)=0.0_dp;end do
   end function pc_skeleton

   subroutine partial_correlation_test(cor,n,i,j,cond,method,stat,pvalue,df)
      real(dp),intent(in)::cor(:,:)
      integer,intent(in)::n,i,j,cond(:),method
      real(dp),intent(out)::stat,pvalue
      integer,intent(out)::df
      real(dp)::r,den,div
      real(dp),allocatable::rss(:,:),rhs(:),sol(:),sub(:,:)
      integer::k,a,b,info
      k=size(cond);df=n-k-3
      if(df<=0)then;stat=0.0_dp;pvalue=1.0_dp;return;end if
      if(k==0)then
         r=cor(i,j)
      else
         allocate(rss(k,k),rhs(k),sol(k),sub(2,2));rss=cor(cond,cond);sub=0.0_dp
         sub(1,1)=cor(i,i);sub(1,2)=cor(i,j);sub(2,1)=cor(j,i);sub(2,2)=cor(j,j)
         do a=1,2
            if(a==1)then;rhs=cor(cond,i);else;rhs=cor(cond,j);end if
            call solve_linear(rss,rhs,sol,info)
            if(info/=0)then;stat=0.0_dp;pvalue=1.0_dp;return;end if
            do b=1,2
               if(b==1)then;sub(b,a)=sub(b,a)-dot_product(cor(i,cond),sol)
               else;sub(b,a)=sub(b,a)-dot_product(cor(j,cond),sol);end if
            end do
         end do
         den=sqrt(max(tiny(1.0_dp),sub(1,1)*sub(2,2)));r=sub(1,2)/den
      end if
      r=max(-1.0_dp+1.0e-14_dp,min(1.0_dp-1.0e-14_dp,r));div=1.0_dp
      if(method==PC_SPEARMAN)div=1.029563_dp
      stat=abs(0.5_dp*log((1.0_dp+r)/(1.0_dp-r))*sqrt(real(df,dp)))/div
      pvalue=min(1.0_dp,max(0.0_dp,2.0_dp*(1.0_dp-student_t_cdf(stat,real(df,dp)))))
   end subroutine partial_correlation_test

   function neighbor_union(g,i,j) result(v)
      integer,intent(in)::g(:,:),i,j
      integer,allocatable::v(:)
      integer::k,n
      logical,allocatable::seen(:)
      allocate(seen(size(g,1)));seen=.false.;seen(i)=.true.;seen(j)=.true.
      do k=1,size(g,1)
         if((g(i,k)/=0.or.g(j,k)/=0).and..not.seen(k))seen(k)=.true.
      end do
      n=count(seen)-2;allocate(v(max(0,n)));n=0
      do k=1,size(g,1)
         if(k/=i.and.k/=j.and.seen(k))then;n=n+1;v(n)=k;end if
      end do
   end function neighbor_union

   logical function next_combination(c,n) result(ok)
      integer,intent(inout)::c(:)
      integer,intent(in)::n
      integer::i,j,k
      k=size(c);ok=.false.
      do i=k,1,-1
         if(c(i)<n-k+i)then
            c(i)=c(i)+1
            do j=i+1,k;c(j)=c(j-1)+1;end do
            ok=.true.;return
         end if
      end do
   end function next_combination

   logical function any_degree_at_least(g,d) result(ok)
      integer,intent(in)::g(:,:),d
      integer::i
      ok=.false.
      do i=1,size(g,1)
         if(count(g(i,:)/=0)>=d)then;ok=.true.;return;end if
      end do
   end function any_degree_at_least

end module rfast_pc_skel
