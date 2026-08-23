module rfast_algorithms
   use rfast_special, only : dp
   use rfast_distances, only : vector_distance, DIST_EUCLIDEAN
   use rfast_arrays, only : order_real
   implicit none
   private
   public :: knn_classify, knn_regress, floyd_warshall, topological_sort
   public :: next_permutation, prev_permutation, combination_count, spatial_median

contains

   function knn_classify(train_x,train_y,test_x,k) result(pred)
      real(dp),intent(in)::train_x(:,:),test_x(:,:);integer,intent(in)::train_y(:),k;integer::pred(size(test_x,1))
      real(dp)::d(size(train_x,1));integer,allocatable::ord(:);integer::i,j,kk,cmin,cmax,c,best,bestn,nc
      cmin=minval(train_y);cmax=maxval(train_y);nc=cmax-cmin+1;kk=min(k,size(train_x,1))
      do i=1,size(test_x,1)
         do j=1,size(train_x,1);d(j)=vector_distance(test_x(i,:),train_x(j,:),DIST_EUCLIDEAN);end do;ord=order_real(d)
         best=cmin;bestn=-1
         do c=cmin,cmax
            j=count(train_y(ord(1:kk))==c);if(j>bestn)then;bestn=j;best=c;end if
         end do;pred(i)=best
      end do
   end function knn_classify

   function knn_regress(train_x,train_y,test_x,k,weighted) result(pred)
      real(dp),intent(in)::train_x(:,:),train_y(:),test_x(:,:);integer,intent(in)::k;logical,intent(in),optional::weighted
      real(dp)::pred(size(test_x,1)),d(size(train_x,1)),w;integer,allocatable::ord(:);integer::i,j,kk;logical::wt
      wt=.false.;if(present(weighted))wt=weighted;kk=min(k,size(train_x,1))
      do i=1,size(test_x,1)
         do j=1,size(train_x,1);d(j)=vector_distance(test_x(i,:),train_x(j,:),DIST_EUCLIDEAN);end do;ord=order_real(d)
         if(.not.wt)then;pred(i)=sum(train_y(ord(1:kk)))/kk
         else
            pred(i)=0.0_dp;w=0.0_dp
            do j=1,kk
               if(d(ord(j))<=tiny(1.0_dp))then;pred(i)=train_y(ord(j));w=-1.0_dp;exit;end if
               pred(i)=pred(i)+train_y(ord(j))/d(ord(j));w=w+1.0_dp/d(ord(j))
            end do
            if(w>0)pred(i)=pred(i)/w
         end if
      end do
   end function knn_regress

   function floyd_warshall(weight) result(d)
      real(dp),intent(in)::weight(:,:);real(dp)::d(size(weight,1),size(weight,2));integer::i,j,k,n
      n=size(weight,1);d=weight
      do i=1,n;d(i,i)=0.0_dp;end do
      do k=1,n;do j=1,n;do i=1,n
         if(d(i,k)<huge(1.0_dp)/4.and.d(k,j)<huge(1.0_dp)/4)d(i,j)=min(d(i,j),d(i,k)+d(k,j))
      end do;end do;end do
   end function floyd_warshall

   function topological_sort(adj,ok) result(order)
      logical,intent(in)::adj(:,:);logical,intent(out),optional::ok;integer,allocatable::order(:)
      integer::n,indeg(size(adj,1)),q(size(adj,1)),head,tail,u,v,k
      n=size(adj,1);allocate(order(n));indeg=0
      do v=1,n;indeg(v)=count(adj(:,v));end do
      head=1;tail=0;do v=1,n;if(indeg(v)==0)then;tail=tail+1;q(tail)=v;end if;end do;k=0
      do while(head<=tail)
         u=q(head);head=head+1;k=k+1;order(k)=u
         do v=1,n;if(adj(u,v))then;indeg(v)=indeg(v)-1;if(indeg(v)==0)then;tail=tail+1;q(tail)=v;end if;end if;end do
      end do
      if(present(ok))ok=(k==n);if(k<n)order(k+1:n)=0
   end function topological_sort

   logical function next_permutation(a) result(ok)
      integer,intent(inout)::a(:);integer::i,j,tmp,l,r
      i=size(a)-1;do while(i>=1);if(a(i)<a(i+1))exit;i=i-1;end do
      if(i<1)then;ok=.false.;return;end if;j=size(a);do while(a(j)<=a(i));j=j-1;end do
      tmp=a(i);a(i)=a(j);a(j)=tmp;l=i+1;r=size(a)
      do while(l<r);tmp=a(l);a(l)=a(r);a(r)=tmp;l=l+1;r=r-1;end do;ok=.true.
   end function next_permutation

   logical function prev_permutation(a) result(ok)
      integer,intent(inout)::a(:);integer::i,j,tmp,l,r
      i=size(a)-1;do while(i>=1);if(a(i)>a(i+1))exit;i=i-1;end do
      if(i<1)then;ok=.false.;return;end if;j=size(a);do while(a(j)>=a(i));j=j-1;end do
      tmp=a(i);a(i)=a(j);a(j)=tmp;l=i+1;r=size(a)
      do while(l<r);tmp=a(l);a(l)=a(r);a(r)=tmp;l=l+1;r=r-1;end do;ok=.true.
   end function prev_permutation

   pure integer(kind=8) function combination_count(n,k) result(c)
      integer,intent(in)::n,k;integer::i,kk
      if(k<0.or.k>n)then;c=0;return;end if;kk=min(k,n-k);c=1_8
      do i=1,kk;c=c*int(n-kk+i,8)/int(i,8);end do
   end function combination_count

   function spatial_median(x,tol,maxiter) result(m)
      real(dp),intent(in)::x(:,:);real(dp),intent(in),optional::tol;integer,intent(in),optional::maxiter
      real(dp)::m(size(x,2)),old(size(x,2)),num(size(x,2)),den,d,eps;integer::i,it,mi
      m=sum(x,dim=1)/real(size(x,1),dp);eps=1e-8_dp;if(present(tol))eps=tol;mi=1000;if(present(maxiter))mi=maxiter
      do it=1,mi
         old=m;num=0.0_dp;den=0.0_dp
         do i=1,size(x,1);d=sqrt(sum((x(i,:)-old)**2));if(d<=eps)cycle;num=num+x(i,:)/d;den=den+1.0_dp/d;end do
         if(den>0)m=num/den;if(sqrt(sum((m-old)**2))<eps)exit
      end do
   end function spatial_median

end module rfast_algorithms
