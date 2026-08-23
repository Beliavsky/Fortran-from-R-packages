module directional_regression
   use directional_kinds, only : dp
   implicit none
   private
   public :: knn_reg, knnreg_tune
contains
   subroutine knn_reg(xnew,y,x,k,est,spherical,harmonic)
      real(dp),intent(in)::xnew(:,:),y(:,:),x(:,:);integer,intent(in)::k
      real(dp),intent(out)::est(size(xnew,1),size(y,2))
      logical,intent(in),optional::spherical,harmonic
      real(dp)::dist(size(x,1)),nr,den
      integer::idx(size(x,1)),i,j,l,t,tmp,kk
      logical::sph,har
      sph=.false.;if(present(spherical))sph=spherical;har=.false.;if(present(harmonic))har=harmonic
      kk=max(1,min(k,size(x,1)))
      do i=1,size(xnew,1)
         do j=1,size(x,1);dist(j)=sqrt(sum((xnew(i,:)-x(j,:))**2));idx(j)=j;end do
         do j=1,kk
            t=j;do l=j+1,size(x,1);if(dist(idx(l))<dist(idx(t)))t=l;end do
            tmp=idx(j);idx(j)=idx(t);idx(t)=tmp
         end do
         do l=1,size(y,2)
            if(har)then
               den=0.0_dp
               do j=1,kk
                  if(abs(y(idx(j),l))>sqrt(tiny(1.0_dp)))den=den+1.0_dp/y(idx(j),l)
               end do
               if(abs(den)>sqrt(tiny(1.0_dp)))then;est(i,l)=real(kk,dp)/den;else;est(i,l)=0.0_dp;end if
            else
               est(i,l)=sum(y(idx(1:kk),l))/real(kk,dp)
            end if
         end do
         if(sph)then;nr=sqrt(sum(est(i,:)*est(i,:)));if(nr>tiny(1.0_dp))est(i,:)=est(i,:)/nr;end if
      end do
   end subroutine

   subroutine knnreg_tune(y,x,a,best_k,crit,nfolds,spherical,harmonic)
      real(dp),intent(in)::y(:,:),x(:,:);integer,intent(in)::a
      integer,intent(out)::best_k;real(dp),intent(out)::crit(max(0,a-1))
      integer,intent(in),optional::nfolds;logical,intent(in),optional::spherical,harmonic
      integer::nf,n,k,fold,i,j,nt,nr,ix,it
      real(dp),allocatable::xt(:,:),yt(:,:),xv(:,:),yv(:,:),pred(:,:)
      real(dp)::err
      logical::sph,har
      n=size(x,1);nf=min(10,n);if(present(nfolds))nf=max(2,min(n,nfolds));sph=.false.;if(present(spherical))sph=spherical;har=.false.;if(present(harmonic))har=harmonic
      if(a<2)then;best_k=1;return;end if
      crit=0.0_dp
      do fold=1,nf
         nt=count([(mod(i-1,nf)/=fold-1,i=1,n)]);nr=n-nt
         allocate(xt(nt,size(x,2)),yt(nt,size(y,2)),xv(nr,size(x,2)),yv(nr,size(y,2)))
         ix=0;it=0
         do i=1,n
            if(mod(i-1,nf)==fold-1)then;ix=ix+1;xv(ix,:)=x(i,:);yv(ix,:)=y(i,:)
            else;it=it+1;xt(it,:)=x(i,:);yt(it,:)=y(i,:);end if
         end do
         do k=2,a
            allocate(pred(nr,size(y,2)));call knn_reg(xv,yt,xt,k,pred,sph,har)
            if(sph)then
               err=0.0_dp;do i=1,nr;err=err+1.0_dp-dot_product(pred(i,:),yv(i,:));end do;err=err/real(nr,dp)
            else
               err=sum((pred-yv)**2)/real(size(pred),dp)
            end if
            crit(k-1)=crit(k-1)+err/real(nf,dp);deallocate(pred)
         end do
         deallocate(xt,yt,xv,yv)
      end do
      best_k=minloc(crit,dim=1)+1
   end subroutine
end module directional_regression
