module matchingmarkets_helpers
   use matchingmarkets_kinds, only : dp
   implicit none
   private
   public :: pair_combinations, coalition_partitions, consensus_mc

   interface
      subroutine dgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
         integer,intent(in)::n,nrhs,lda,ldb
         integer,intent(out)::ipiv(*)
         integer,intent(out)::info
         double precision,intent(inout)::a(lda,*),b(ldb,*)
      end subroutine dgesv
   end interface
contains
   function pair_combinations(x) result(pairs)
      integer,intent(in)::x(:)
      integer,allocatable::pairs(:,:)
      integer::i,j,k,n
      n=size(x);allocate(pairs(2,n*(n-1)/2));k=0
      do i=1,n-1
         do j=i+1,n
            k=k+1;pairs(:,k)=[x(i),x(j)]
         end do
      end do
   end function pair_combinations

   function coalition_partitions(n1,n2,drop_observed) result(parts)
      integer,intent(in)::n1,n2
      logical,intent(in),optional::drop_observed
      integer,allocatable::parts(:,:)
      logical::drop
      integer::n,k,countall
      integer,allocatable::comb(:),tmp(:,:)
      drop=.false.;if(present(drop_observed))drop=drop_observed
      n=n1+n2
      if(n1<0.or.n2<0)error stop 'coalition_partitions: invalid sizes'
      countall=choose_int(n,n1)
      allocate(tmp(n,countall));allocate(comb(n1));comb=0;k=0
      call rec(1,1)
      allocate(parts(n,k));if(k>0)parts=tmp(:,:k)
   contains
      recursive subroutine rec(pos,start)
         integer,intent(in)::pos,start
         integer::v,j,q
         logical::isobs
         integer::col(n)
         if(pos>n1)then
            col=0
            do j=1,n1;col(j)=comb(j);end do
            q=n1
            do v=1,n
               if(.not.any(comb==v))then;q=q+1;col(q)=v;end if
            end do
            isobs=all(comb==[(j,j=1,n1)]) .or. all(comb==[(j,j=n2+1,n)])
            if(drop.and.isobs)return
            k=k+1;tmp(:,k)=col;return
         end if
         do v=start,n-(n1-pos)
            comb(pos)=v;call rec(pos+1,v+1)
         end do
      end subroutine rec
   end function coalition_partitions

   function consensus_mc(subchain) result(theta)
      real(dp),intent(in)::subchain(:,:,:)
      real(dp),allocatable::theta(:,:)
      integer::d,ns,m,k,i,info,j
      real(dp),allocatable::covm(:,:,:),invm(:,:,:),sumw(:,:),rhs(:,:),mean(:),center(:,:),work(:,:)
      integer,allocatable::ipiv(:)
      d=size(subchain,1);ns=size(subchain,2);m=size(subchain,3)
      allocate(theta(d,ns))
      if(m==1)then;theta=subchain(:,:,1);return;end if
      allocate(covm(d,d,m),invm(d,d,m),sumw(d,d),mean(d),center(d,ns),ipiv(d))
      sumw=0.0_dp
      do k=1,m
         mean=sum(subchain(:,:,k),dim=2)/real(ns,dp)
         do i=1,ns;center(:,i)=subchain(:,i,k)-mean;end do
         covm(:,:,k)=matmul(center,transpose(center))/real(max(1,ns-1),dp)
         allocate(work(d,d));work=0.0_dp;do j=1,d;work(j,j)=1.0_dp;end do
         call dgesv(d,d,covm(:,:,k),d,ipiv,work,d,info)
         if(info/=0)error stop 'consensus_mc: singular subchain covariance'
         invm(:,:,k)=work;sumw=sumw+work;deallocate(work)
      end do
      allocate(work(d,d));work=0.0_dp;do j=1,d;work(j,j)=1.0_dp;end do
      call dgesv(d,d,sumw,d,ipiv,work,d,info)
      if(info/=0)error stop 'consensus_mc: singular weight sum'
      allocate(rhs(d,1))
      do i=1,ns
         rhs=0.0_dp
         do k=1,m;rhs(:,1)=rhs(:,1)+matmul(invm(:,:,k),subchain(:,i,k));end do
         theta(:,i)=matmul(work,rhs(:,1))
      end do
   end function consensus_mc

   pure integer function choose_int(n,k) result(v)
      integer,intent(in)::n,k
      integer::i,kk
      if(k<0.or.k>n)then;v=0;return;end if
      kk=min(k,n-k);v=1
      do i=1,kk;v=v*(n-kk+i)/i;end do
   end function choose_int
end module matchingmarkets_helpers
