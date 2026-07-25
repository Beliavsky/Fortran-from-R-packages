! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of robustbase.
! It may be redistributed and/or modified under GPL version 2 or later.
module robustbase_utilities
   use robustbase_kinds, only: dp
   use robustbase_sort, only: median
   use robustbase_scale, only: mad_scale
   use robustbase_linalg, only: matrix_rank
   implicit none
   private
   public :: row_medians, column_medians, robust_standardize, independent_columns, full_rank_matrix
contains
   subroutine row_medians(x,m)
      real(dp),intent(in)::x(:,:);real(dp),intent(out)::m(:)
      integer::i
      if(size(m)/=size(x,1)) error stop "row_medians: size mismatch"
      do i=1,size(x,1);m(i)=median(x(i,:));end do
   end subroutine
   subroutine column_medians(x,m)
      real(dp),intent(in)::x(:,:);real(dp),intent(out)::m(:)
      integer::j
      if(size(m)/=size(x,2)) error stop "column_medians: size mismatch"
      do j=1,size(x,2);m(j)=median(x(:,j));end do
   end subroutine
   subroutine robust_standardize(x,z,center,scale)
      real(dp),intent(in)::x(:,:);real(dp),intent(out)::z(:,:),center(:),scale(:)
      integer::j
      if(any(shape(z)/=shape(x)) .or. size(center)/=size(x,2) .or. size(scale)/=size(x,2)) error stop "robust_standardize: size mismatch"
      do j=1,size(x,2)
         center(j)=median(x(:,j));scale(j)=mad_scale(x(:,j),center(j));if(scale(j)<=1.0e-14_dp)scale(j)=1.0_dp;z(:,j)=(x(:,j)-center(j))/scale(j)
      end do
   end subroutine
   subroutine independent_columns(x,keep,rank)
      real(dp),intent(in)::x(:,:)
      logical,intent(out)::keep(:)
      integer,intent(out)::rank
      real(dp),allocatable::a(:,:)
      integer,allocatable::accepted(:)
      integer::j,k,r0
      if(size(keep)/=size(x,2)) error stop "independent_columns: size mismatch"
      allocate(accepted(size(x,2)))
      keep=.false.;rank=0;accepted=0
      do j=1,size(x,2)
         allocate(a(size(x,1),rank+1))
         do k=1,rank
            a(:,k)=x(:,accepted(k))
         end do
         a(:,rank+1)=x(:,j)
         r0=matrix_rank(a)
         deallocate(a)
         if(r0>rank) then
            rank=r0;accepted(rank)=j;keep(j)=.true.
         end if
      end do
   end subroutine

   subroutine full_rank_matrix(x,x_full,kept,rank)
      real(dp),intent(in)::x(:,:)
      real(dp),allocatable,intent(out)::x_full(:,:)
      integer,allocatable,intent(out)::kept(:)
      integer,intent(out)::rank
      logical,allocatable::mask(:)
      real(dp),allocatable::xt(:,:),work(:,:)
      integer::j,k
      if(size(x,1)>=size(x,2))then
         allocate(mask(size(x,2)))
         call independent_columns(x,mask,rank)
         allocate(kept(rank),x_full(size(x,1),rank));k=0
         do j=1,size(x,2)
            if(mask(j))then;k=k+1;kept(k)=j;x_full(:,k)=x(:,j);end if
         end do
      else
         xt=transpose(x);allocate(mask(size(xt,2)))
         call independent_columns(xt,mask,rank)
         allocate(kept(rank),work(size(xt,1),rank));k=0
         do j=1,size(xt,2)
            if(mask(j))then;k=k+1;kept(k)=j;work(:,k)=xt(:,j);end if
         end do
         x_full=transpose(work)
      end if
   end subroutine full_rank_matrix
end module robustbase_utilities
