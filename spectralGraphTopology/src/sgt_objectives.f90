! SPDX-License-Identifier: GPL-3.0-only
module sgt_objectives
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf
   use sgt_kinds, only : dp
   use sgt_linalg, only : symmetric_eigen_jacobi, frobenius_norm, trace_matrix
   use sgt_updates, only : spectral_reconstruction
   implicit none
   private
   public :: laplacian_negative_log_likelihood, laplacian_prior, laplacian_objective
   public :: bipartite_negative_log_likelihood, bipartite_prior, bipartite_objective
   public :: joint_prior, joint_objective, vanilla_objective
contains
   function laplacian_negative_log_likelihood(lw,lambda,kmat) result(value)
      real(dp), intent(in) :: lw(:,:),lambda(:),kmat(:,:)
      real(dp) :: value
      if (size(lambda)==0 .or. minval(lambda)<=0.0_dp) then
         value=ieee_value(0.0_dp,ieee_positive_inf)
      else
         value=-sum(log(lambda))+sum(kmat*transpose(lw))
      end if
   end function laplacian_negative_log_likelihood

   function laplacian_prior(beta,lw,lambda,u) result(value)
      real(dp), intent(in) :: beta,lw(:,:),lambda(:),u(:,:)
      real(dp) :: value
      real(dp), allocatable :: target(:,:)
      target=spectral_reconstruction(u,lambda)
      value=0.5_dp*beta*frobenius_norm(lw-target)**2
   end function laplacian_prior

   function laplacian_objective(lw,u,lambda,kmat,beta) result(value)
      real(dp), intent(in) :: lw(:,:),u(:,:),lambda(:),kmat(:,:),beta
      real(dp) :: value
      value=laplacian_negative_log_likelihood(lw,lambda,kmat)+laplacian_prior(beta,lw,lambda,u)
   end function laplacian_objective

   function bipartite_negative_log_likelihood(lw,kmat,jmat) result(value)
      real(dp), intent(in) :: lw(:,:),kmat(:,:),jmat(:,:)
      real(dp) :: value
      real(dp), allocatable :: eig(:),vec(:,:)
      integer :: status
      call symmetric_eigen_jacobi(lw+jmat,eig,vec,status)
      if (size(eig)==0 .or. minval(eig)<=0.0_dp) then
         value=ieee_value(0.0_dp,ieee_positive_inf)
      else
         value=-sum(log(eig))+sum(kmat*transpose(lw))
      end if
   end function bipartite_negative_log_likelihood

   function bipartite_prior(nu,aw,psi,v) result(value)
      real(dp), intent(in) :: nu,aw(:,:),psi(:),v(:,:)
      real(dp) :: value
      real(dp), allocatable :: target(:,:)
      target=spectral_reconstruction(v,psi)
      value=0.5_dp*nu*frobenius_norm(aw-target)**2
   end function bipartite_prior

   function bipartite_objective(aw,lw,v,psi,kmat,jmat,nu) result(value)
      real(dp), intent(in) :: aw(:,:),lw(:,:),v(:,:),psi(:),kmat(:,:),jmat(:,:),nu
      real(dp) :: value
      value=bipartite_negative_log_likelihood(lw,kmat,jmat)+bipartite_prior(nu,aw,psi,v)
   end function bipartite_objective

   function joint_prior(beta,nu,lw,aw,u,v,lambda,psi) result(value)
      real(dp), intent(in) :: beta,nu,lw(:,:),aw(:,:),u(:,:),v(:,:),lambda(:),psi(:)
      real(dp) :: value
      value=laplacian_prior(beta,lw,lambda,u)+bipartite_prior(nu,aw,psi,v)
   end function joint_prior

   function joint_objective(lw,aw,u,v,lambda,psi,beta,nu,kmat) result(value)
      real(dp), intent(in) :: lw(:,:),aw(:,:),u(:,:),v(:,:),lambda(:),psi(:),beta,nu,kmat(:,:)
      real(dp) :: value
      value=laplacian_negative_log_likelihood(lw,lambda,kmat)+ &
         joint_prior(beta,nu,lw,aw,u,v,lambda,psi)
   end function joint_objective

   function vanilla_objective(theta,kmat) result(value)
      real(dp), intent(in) :: theta(:,:),kmat(:,:)
      real(dp) :: value
      real(dp), allocatable :: eig(:),vec(:,:)
      integer :: p,status
      p=size(theta,1)
      call symmetric_eigen_jacobi(theta,eig,vec,status)
      if (p<2 .or. size(eig)/=p .or. minval(eig(2:p))<=0.0_dp) then
         value=ieee_value(0.0_dp,ieee_positive_inf)
      else
         value=sum(theta*transpose(kmat))-sum(log(eig(2:p)))
      end if
   end function vanilla_objective
end module sgt_objectives
