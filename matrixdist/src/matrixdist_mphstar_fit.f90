! SPDX-License-Identifier: GPL-3.0-only
module matrixdist_mphstar_fit
   use r_compat, only: dp
   use matrixdist_types, only: ph_type
   use matrixdist_linalg, only: matrix_exponential, matrix_vanloan
   use matrixdist_ph, only: ph_exit_rates
   use matrixdist_em, only: emstep_ph
   use matrixdist_transformations, only: tvr_ph, plus_states, reward_sanity_check
   implicit none
   private
   public :: marginal_sojourn_expectation, emstep_mphstar
contains
   function outer_product(a,b) result(c)
      real(dp),intent(in)::a(:),b(:)
      real(dp)::c(size(a),size(b))
      integer::i,j
      do i=1,size(a)
      do j=1,size(b)
      c(i,j)=a(i)*b(j)
      end do
      end do
   end function outer_product

   function marginal_sojourn_expectation(alpha,s,obs,weight) result(z)
      ! Sum E[Z_i | X=x] over observations for a PH marginal.
      real(dp),intent(in)::alpha(:),s(:,:),obs(:)
      real(dp),intent(in),optional::weight(:)
      real(dp)::z(size(alpha)),w,density
      real(dp),allocatable::exitv(:),b(:,:),vl(:,:),c(:,:)
      integer::p,k,i
      p=size(alpha)
      z=0.0_dp
      exitv=ph_exit_rates(s)
      b=outer_product(exitv,alpha)
      do k=1,size(obs)
         w=1.0_dp
         if(present(weight))w=weight(k)
         vl=matrix_exponential(matrix_vanloan(s,s,b)*obs(k))
         density=dot_product(alpha,matmul(vl(1:p,1:p),exitv))
         if(density<=tiny(1.0_dp))cycle
         c=vl(1:p,p+1:2*p)
         do i=1,p
         z(i)=z(i)+w*c(i,i)/density
         end do
      end do
   end function marginal_sojourn_expectation

   subroutine emstep_mphstar(alpha,s,reward,y,weight,reward_tol)
      ! Algorithm-2 reward update used by matrixdist's MPHstar fit.  Because
      ! each reward row sums to one, rowSums(Y) follows the originating PH;
      ! that sufficient statistic updates alpha/S, while marginal time-change
      ! models update the reward proportions.
      real(dp),intent(inout)::alpha(:),s(:,:),reward(:,:)
      real(dp),intent(in)::y(:,:)
      real(dp),intent(in),optional::weight(:,:)
      real(dp),intent(in),optional::reward_tol
      real(dp)::tol
      real(dp),allocatable::sumy(:),ztot(:),zj(:,:),wj(:)
      integer,allocatable::pos(:)
      type(ph_type)::mar
      integer::n,p,d,j,k
      tol=1.0e-4_dp
      if(present(reward_tol))tol=reward_tol
      n=size(y,1)
      d=size(y,2)
      p=size(alpha)
      if(size(reward,1)/=p .or. size(reward,2)/=d)error stop 'emstep_mphstar: dimension mismatch'
      call reward_sanity_check(reward,tol)
      allocate(zj(p,d))
      zj=0.0_dp
      ! Marginal conditional sojourn expectations use the pre-update model.
      do j=1,d
         pos=plus_states(reward(:,j))
         if(size(pos)==0)cycle
         mar=tvr_ph(alpha,s,reward(:,j))
         if(present(weight))then
            allocate(wj(n))
            wj=weight(:,j)
            ztot=marginal_sojourn_expectation(mar%alpha,mar%s,y(:,j),wj)
            deallocate(wj)
         else
            ztot=marginal_sojourn_expectation(mar%alpha,mar%s,y(:,j))
         end if
         do k=1,size(pos)
         zj(pos(k),j)=ztot(k)
         end do
      end do
      ! The originating path duration is the sum of rewards along the path.
      allocate(sumy(n))
      sumy=sum(y,dim=2)
      call emstep_ph(alpha,s,sumy)
      do k=1,p
         if(sum(zj(k,:))<=tiny(1.0_dp))then
            reward(k,:)=0.0_dp
         else
            reward(k,:)=max(0.0_dp,zj(k,:)/sum(zj(k,:)))
         end if
      end do
   end subroutine emstep_mphstar
end module matrixdist_mphstar_fit
