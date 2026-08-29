! SPDX-License-Identifier: GPL-3.0-only
module matrixdist_numerics
   use r_compat, only: dp
   use matrixdist_linalg, only: eye_matrix, max_diagonal
   implicit none
   private
   public :: default_step_length, find_n, uniformized_exponential, pow2_matrix
   public :: reverse_gev_transform
contains
   function default_step_length(s) result(h)
      real(dp),intent(in)::s(:,:)
      real(dp)::h
      integer::i
      h=-0.1_dp/s(1,1)
      do i=2,size(s,1)
      h=min(h,-0.1_dp/s(i,i))
      end do
   end function default_step_length

   integer function find_n(h,lambda) result(n)
      real(dp),intent(in)::h,lambda
      real(dp)::term,cum
      integer::k
      if(lambda<0.0_dp)error stop "find_n: lambda must be nonnegative"
      if(lambda==0.0_dp)then
      n=0
      return
      end if
      term=exp(-lambda)
      cum=term
      k=0
      do while(cum<1.0_dp-h .and. k<1000000)
         k=k+1
         term=term*lambda/real(k,dp)
         cum=cum+term
      end do
      n=k
   end function find_n

   function uniformized_exponential(s,x,epsilon) result(ea)
      real(dp),intent(in)::s(:,:),x
      real(dp),intent(in),optional::epsilon
      real(dp),allocatable::ea(:,:),p(:,:),term(:,:)
      real(dp)::eps,a,coef
      integer::n,k,scale
      eps=1.0e-12_dp
      if(present(epsilon))eps=epsilon
      a=maxval(-[(s(k,k),k=1,size(s,1))])
      if(a<=0.0_dp)then
      ea=eye_matrix(size(s,1))
      return
      end if
      scale=0
      if(a*x>1.0_dp)scale=max(0,int(log(a*x)/log(2.0_dp))+1)
      p=eye_matrix(size(s,1))+s/a
      n=find_n(eps,1.0_dp)
      ea=eye_matrix(size(s,1))
      term=eye_matrix(size(s,1))
      coef=1.0_dp
      do k=1,n
         term=matmul(p,term)/real(k,dp)
         coef=(a*x/(2.0_dp**scale))**k
         ea=ea+coef*term
      end do
      ea=ea*exp(-a*x/(2.0_dp**scale))
      do k=1,scale
      ea=matmul(ea,ea)
      end do
   end function uniformized_exponential

   subroutine pow2_matrix(n,a)
      integer,intent(in)::n
      real(dp),intent(inout)::a(:,:)
      integer::i
      do i=1,n
      a=matmul(a,a)
      end do
   end subroutine pow2_matrix

   subroutine reverse_gev_transform(obs,weights,beta,trans_obs,trans_weights)
      real(dp),intent(in)::obs(:),weights(:),beta(3)
      real(dp),intent(out)::trans_obs(size(obs)),trans_weights(size(obs))
      integer::i,n
      n=size(obs)
      do i=1,n
         if(beta(3)==0.0_dp)then
            trans_obs(i)=exp(-(obs(n-i+1)-beta(1))/beta(2))
         else
            trans_obs(i)=(1.0_dp+beta(3)*(obs(n-i+1)-beta(1))/beta(2))**(-1.0_dp/beta(3))
         end if
         trans_weights(i)=weights(n-i+1)
      end do
   end subroutine reverse_gev_transform
end module matrixdist_numerics
