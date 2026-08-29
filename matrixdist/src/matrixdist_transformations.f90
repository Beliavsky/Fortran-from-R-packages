! SPDX-License-Identifier: GPL-3.0-only
module matrixdist_transformations
   use r_compat, only: dp, runif1
   use matrixdist_linalg
   use matrixdist_types, only: ph_type, dph_type
   implicit none
   private
   public :: tvr_ph, tvr_dph, linear_combination, merge_matrices
   public :: random_reward, reward_sanity_check, n_pos, plus_states

contains

   integer function n_pos(r) result(n)
      real(dp),intent(in)::r(:)
      n=count(r>0.0_dp)
   end function n_pos

   function plus_states(r) result(idx)
      real(dp),intent(in)::r(:)
      integer,allocatable::idx(:)
      integer::i,j
      allocate(idx(count(r>0.0_dp)))
      j=0
      do i=1,size(r)
      if(r(i)>0.0_dp)then
      j=j+1
      idx(j)=i
      end if
      end do
   end function plus_states

   function random_reward(p,d) result(r)
      integer,intent(in)::p,d
      real(dp)::r(p,d),sm
      integer::i,j
      do i=1,p
         do j=1,d
         r(i,j)=runif1()
         end do
         sm=sum(r(i,:))
         if(sm<=0.0_dp)then
         r(i,:)=1.0_dp/real(d,dp)
         else
         r(i,:)=r(i,:)/sm
         end if
      end do
   end function random_reward

   subroutine reward_sanity_check(r,tol)
      real(dp),intent(inout)::r(:,:)
      real(dp),intent(in)::tol
      real(dp)::sm,miss,w
      integer::i,j
      where(r<tol)r=0.0_dp
      do i=1,size(r,1)
         sm=sum(r(i,:))
         if(sm<=0.0_dp)cycle
         if(sm/=1.0_dp)then
            miss=1.0_dp-sm
            do j=1,size(r,2)
               if(r(i,j)>0.0_dp)then
               w=r(i,j)/sm
               r(i,j)=r(i,j)+miss*w
               end if
            end do
         end if
      end do
   end subroutine reward_sanity_check

   function merge_matrices(s11,s12,s22) result(s)
      real(dp),intent(in)::s11(:,:),s12(:,:),s22(:,:)
      real(dp),allocatable::s(:,:)
      integer::p1,p2
      p1=size(s11,1)
      p2=size(s22,1)
      allocate(s(p1+p2,p1+p2))
      s=0.0_dp
      s(1:p1,1:p1)=s11
      s(1:p1,p1+1:)=s12
      s(p1+1:,p1+1:)=s22
   end function merge_matrices

   subroutine partition_indices(r,keep,zero)
      real(dp),intent(in)::r(:)
      integer,allocatable,intent(out)::keep(:),zero(:)
      integer::i,j,k
      allocate(keep(count(r/=0.0_dp)),zero(count(r==0.0_dp)))
      j=0
      k=0
      do i=1,size(r)
         if(r(i)==0.0_dp)then
         k=k+1
         zero(k)=i
         else
         j=j+1
         keep(j)=i
         end if
      end do
   end subroutine partition_indices

   function submatrix(a,rows,cols) result(b)
      real(dp),intent(in)::a(:,:)
      integer,intent(in)::rows(:),cols(:)
      real(dp)::b(size(rows),size(cols))
      integer::i,j
      do i=1,size(rows)
      do j=1,size(cols)
      b(i,j)=a(rows(i),cols(j))
      end do
      end do
   end function submatrix

   function subvector(a,idx) result(b)
      real(dp),intent(in)::a(:)
      integer,intent(in)::idx(:)
      real(dp)::b(size(idx))
      integer::i
      do i=1,size(idx)
      b(i)=a(idx(i))
      end do
   end function subvector

   function tvr_ph(alpha,s,reward) result(z)
      real(dp),intent(in)::alpha(:),s(:,:),reward(:)
      type(ph_type)::z
      integer,allocatable::kp(:),z0(:)
      real(dp),allocatable::spp(:,:),sp0(:,:),s0p(:,:),s00(:,:),inv0(:,:)
      real(dp),allocatable::ap(:),a0(:),saux(:,:)
      integer::i
      call partition_indices(reward,kp,z0)
      if(size(kp)==0) error stop "tvr_ph: all rewards are zero"
      ap=subvector(alpha,kp)
      if(size(z0)==0)then
         z%alpha=ap
         z%s=submatrix(s,kp,kp)
      else
         a0=subvector(alpha,z0)
         spp=submatrix(s,kp,kp)
         sp0=submatrix(s,kp,z0)
         s0p=submatrix(s,z0,kp)
         s00=submatrix(s,z0,z0)
         inv0=matrix_inverse(-s00)
         z%alpha=ap+matmul(a0,matmul(inv0,s0p))
         z%s=spp+matmul(sp0,matmul(inv0,s0p))
      end if
      ! Reward time change: holding rates are divided by positive reward.
      do i=1,size(kp)
      z%s(i,:)=z%s(i,:)/reward(kp(i))
      end do
   end function tvr_ph

   function tvr_dph(alpha,s,reward) result(z)
      real(dp),intent(in)::alpha(:),s(:,:),reward(:)
      type(dph_type)::z
      integer,allocatable::kp(:),z0(:)
      real(dp),allocatable::spp(:,:),sp0(:,:),s0p(:,:),s00(:,:),inv0(:,:),ap(:),a0(:)
      call partition_indices(reward,kp,z0)
      if(size(kp)==0) error stop "tvr_dph: all rewards are zero"
      ap=subvector(alpha,kp)
      if(size(z0)==0)then
      z%alpha=ap
      z%s=submatrix(s,kp,kp)
      return
      end if
      a0=subvector(alpha,z0)
      spp=submatrix(s,kp,kp)
      sp0=submatrix(s,kp,z0)
      s0p=submatrix(s,z0,kp)
      s00=submatrix(s,z0,z0)
      inv0=matrix_inverse(eye_matrix(size(z0))-s00)
      z%alpha=ap+matmul(a0,matmul(inv0,s0p))
      z%s=spp+matmul(sp0,matmul(inv0,s0p))
   end function tvr_dph

   function linear_combination(w,alpha,s,reward) result(z)
      real(dp),intent(in)::w(:),alpha(:),s(:,:),reward(:,:)
      type(ph_type)::z
      real(dp),allocatable::rw(:)
      rw=matmul(reward,w)
      z=tvr_ph(alpha,s,rw)
   end function linear_combination

end module matrixdist_transformations
