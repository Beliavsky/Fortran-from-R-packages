! SPDX-License-Identifier: MIT
module cgnm_kmeans
   use cgnm_kinds, only : dp
   use cgnm_utils, only : column_mean_sd, elbow_index, uniform01
   implicit none
   private
   public :: optimal_kmeans_labels
contains
   subroutine kmeans_lloyd(x, k, labels, within, max_iter)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: k, max_iter
      integer, intent(out) :: labels(size(x,1))
      real(dp), intent(out) :: within(k)
      real(dp), allocatable :: centers(:,:), newc(:,:), dist(:)
      integer, allocatable :: counts(:), old(:), chosen(:)
      integer :: n,p,i,c,it,ic,tries
      real(dp) :: best, d, u
      n=size(x,1); p=size(x,2)
      allocate(centers(k,p),newc(k,p),dist(k),counts(k),old(n),chosen(k))
      chosen=0
      do c=1,k
         tries=0
         do
            u=uniform01(); ic=1+min(n-1,int(u*real(n,dp)))
            tries=tries+1
            if (c==1 .or. all(chosen(1:c-1)/=ic) .or. tries>10*n) exit
         end do
         chosen(c)=ic; centers(c,:)=x(ic,:)
      end do
      labels=1
      do it=1,max_iter
         old=labels
         do i=1,n
            best=huge(1.0_dp); ic=1
            do c=1,k
               d=sum((x(i,:)-centers(c,:))**2)
               if (d<best) then; best=d; ic=c; end if
            end do
            labels(i)=ic
         end do
         counts=0; newc=0.0_dp
         do i=1,n
            counts(labels(i))=counts(labels(i))+1
            newc(labels(i),:)=newc(labels(i),:)+x(i,:)
         end do
         do c=1,k
            if (counts(c)>0) then
               centers(c,:)=newc(c,:)/real(counts(c),dp)
            else
               u=uniform01(); ic=1+min(n-1,int(u*real(n,dp)))
               centers(c,:)=x(ic,:)
            end if
         end do
         if (all(labels==old)) exit
      end do
      within=0.0_dp
      do i=1,n
         within(labels(i))=within(labels(i))+sum((x(i,:)-centers(labels(i),:))**2)
      end do
   end subroutine kmeans_lloyd

   subroutine optimal_kmeans_labels(x, labels, max_iter)
      real(dp), intent(in) :: x(:,:)
      integer, intent(out) :: labels(size(x,1))
      integer, intent(in) :: max_iter
      real(dp), allocatable :: xn(:,:), mu(:), sd(:), residual(:), within(:)
      integer, allocatable :: tmp(:), counts(:)
      integer :: n,p,k,kmax,el,kk,c
      n=size(x,1); p=size(x,2)
      allocate(xn(n,p),mu(p),sd(p))
      call column_mean_sd(x,mu,sd)
      xn=x
      do c=1,p
         if (sd(c)>0.0_dp) xn(:,c)=(x(:,c)-mu(c))/sd(c)
      end do
      kmax=min(nint(real(n,dp)/real(max(1,p),dp)/2.0_dp),10)
      kmax=max(1,min(kmax,n))
      if (kmax<=1) then
         labels=1; return
      end if
      allocate(residual(kmax-1))
      do k=2,kmax
         allocate(tmp(n),within(k))
         call kmeans_lloyd(xn,k,tmp,within,max_iter)
         residual(k-1)=sum(within)
         deallocate(tmp,within)
      end do
      el=elbow_index(residual)+1
      kk=min(max(1,el),kmax)
      do
         allocate(within(kk),counts(kk))
         call kmeans_lloyd(xn,kk,labels,within,max_iter)
         counts=0
         do c=1,n
            counts(labels(c))=counts(labels(c))+1
         end do
         if (kk<=1) then
            deallocate(within,counts); exit
         end if
         if (minval(counts)>=p .and. minval(within)>0.0_dp) then
            deallocate(within,counts); exit
         end if
         deallocate(within,counts)
         kk=kk-1
      end do
   end subroutine optimal_kmeans_labels
end module cgnm_kmeans
