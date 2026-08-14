! SPDX-License-Identifier: GPL-2.0-only
module ks_kda
   use ks_kinds, only: dp
   use ks_kde, only: kde_model, fit_kde, kde_pdf
   use ks_bandwidth, only: hns_matrix
   implicit none
   private
   public :: kda_model, fit_kda, predict_kda, confusion_matrix, classification_error

   type :: kda_model
      integer :: d=0, ngroup=0
      integer, allocatable :: labels(:)
      real(dp), allocatable :: prior(:)
      type(kde_model), allocatable :: density(:)
   end type kda_model
contains
   subroutine unique_sorted_int(x,u)
      integer,intent(in)::x(:)
      integer,allocatable,intent(out)::u(:)
      integer,allocatable::tmp(:)
      integer::i,j,k,t,n
      allocate(tmp(size(x)));tmp=x;n=size(x)
      do i=2,n
         t=tmp(i);j=i-1
         do while(j>=1)
            if(tmp(j)<=t)exit
            tmp(j+1)=tmp(j);j=j-1
         end do
         tmp(j+1)=t
      end do
      k=1
      do i=2,n
         if(tmp(i)/=tmp(k))then;k=k+1;tmp(k)=tmp(i);end if
      end do
      allocate(u(k));u=tmp(1:k)
   end subroutine unique_sorted_int

   subroutine fit_kda(x,group,model,Hs,prior,weights,info)
      real(dp),intent(in)::x(:,:)
      integer,intent(in)::group(:)
      type(kda_model),intent(out)::model
      real(dp),intent(in),optional::Hs(:,:,:),prior(:),weights(:)
      integer,intent(out),optional::info
      integer::g,i,nj,ierr
      integer,allocatable::ind(:)
      real(dp),allocatable::xx(:,:),ww(:)
      real(dp)::H(size(x,2),size(x,2))
      call unique_sorted_int(group,model%labels)
      model%ngroup=size(model%labels);model%d=size(x,2)
      allocate(model%prior(model%ngroup),model%density(model%ngroup))
      if(present(prior))then
         model%prior=prior/sum(prior)
      else
         do g=1,model%ngroup;model%prior(g)=real(count(group==model%labels(g)),dp)/real(size(group),dp);end do
      end if
      ierr=0
      do g=1,model%ngroup
         nj=count(group==model%labels(g));allocate(ind(nj));ind=pack([(i,i=1,size(group))],group==model%labels(g))
         allocate(xx(nj,model%d));xx=x(ind,:)
         if(present(weights))then;allocate(ww(nj));ww=weights(ind);end if
         if(present(Hs))then;H=Hs(:,:,g);else;call hns_matrix(xx,H);end if
         if(present(weights))then
            call fit_kde(xx,model%density(g),H,ww,ierr)
         else
            call fit_kde(xx,model%density(g),H,info=ierr)
         end if
         deallocate(ind,xx);if(allocated(ww))deallocate(ww)
         if(ierr/=0)exit
      end do
      if(present(info))info=ierr
   end subroutine fit_kda

   subroutine predict_kda(model,x,group,posterior)
      type(kda_model),intent(in)::model
      real(dp),intent(in)::x(:,:)
      integer,intent(out)::group(size(x,1))
      real(dp),intent(out),optional::posterior(size(x,1),model%ngroup)
      real(dp),allocatable::score(:,:),dens(:)
      real(dp)::s
      integer::g,i,imax
      allocate(score(size(x,1),model%ngroup),dens(size(x,1)))
      do g=1,model%ngroup
         call kde_pdf(model%density(g),x,dens);score(:,g)=model%prior(g)*dens
      end do
      do i=1,size(x,1)
         imax=maxloc(score(i,:),dim=1);group(i)=model%labels(imax)
         if(present(posterior))then
            s=sum(score(i,:));if(s>0.0_dp)then;posterior(i,:)=score(i,:)/s;else;posterior(i,:)=1.0_dp/real(model%ngroup,dp);end if
         end if
      end do
   end subroutine predict_kda

   subroutine confusion_matrix(truth,pred,labels,cm)
      integer,intent(in)::truth(:),pred(:),labels(:)
      integer,intent(out)::cm(size(labels),size(labels))
      integer::i,j,k
      cm=0
      do k=1,size(truth)
         i=findloc(labels,truth(k),dim=1);j=findloc(labels,pred(k),dim=1)
         if(i>0 .and. j>0)cm(i,j)=cm(i,j)+1
      end do
   end subroutine confusion_matrix

   function classification_error(truth,pred) result(e)
      integer,intent(in)::truth(:),pred(:)
      real(dp)::e
      e=real(count(truth/=pred),dp)/real(size(truth),dp)
   end function classification_error
end module ks_kda
