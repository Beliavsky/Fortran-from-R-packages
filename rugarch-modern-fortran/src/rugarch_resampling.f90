! Part of the experimental modern Fortran translation of rugarch 1.5-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original rugarch authors retain copyright; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-3.0-only

module rugarch_resampling
   use rugarch_kinds, only : dp
   implicit none
   private

   integer,parameter,public::bootstrap_stationary=1
   integer,parameter,public::bootstrap_block=2

   public :: bootstrap_indices, bootstrap_sample

contains

   subroutine bootstrap_indices(n,nboot,mean_block,method,indices)
      integer,intent(in)::n,nboot,mean_block,method
      integer,intent(out)::indices(n,nboot)
      real(dp)::u
      integer::b,i,start,offset,block_length

      if(n<=0 .or. nboot<=0)error stop 'bootstrap_indices: invalid dimensions'
      block_length=max(1,mean_block)
      select case(method)
      case(bootstrap_stationary)
         do b=1,nboot
            call random_number(u)
            indices(1,b)=1+int(u*real(n,dp))
            do i=2,n
               call random_number(u)
               if(u<1.0_dp/real(block_length,dp))then
                  call random_number(u)
                  indices(i,b)=1+int(u*real(n,dp))
               else
                  indices(i,b)=1+mod(indices(i-1,b),n)
               end if
            end do
         end do
      case(bootstrap_block)
         do b=1,nboot
            i=1
            do while(i<=n)
               call random_number(u)
               start=1+int(u*real(n,dp))
               do offset=0,block_length-1
                  if(i>n)exit
                  indices(i,b)=1+mod(start-1+offset,n)
                  i=i+1
               end do
            end do
         end do
      case default
         error stop 'bootstrap_indices: unknown method'
      end select
   end subroutine bootstrap_indices

   subroutine bootstrap_sample(data,nboot,mean_block,method,samples,indices)
      real(dp),intent(in)::data(:)
      integer,intent(in)::nboot,mean_block,method
      real(dp),intent(out)::samples(size(data),nboot)
      integer,intent(out),optional::indices(size(data),nboot)
      integer,allocatable::idx(:,:)
      integer::b,i,n
      n=size(data)
      allocate(idx(n,nboot))
      call bootstrap_indices(n,nboot,mean_block,method,idx)
      do b=1,nboot
         do i=1,n
            samples(i,b)=data(idx(i,b))
         end do
      end do
      if(present(indices))indices=idx
   end subroutine bootstrap_sample

end module rugarch_resampling
