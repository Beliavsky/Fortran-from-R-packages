! SPDX-License-Identifier: MIT
module cgnm_extensions
   use cgnm_kinds, only : dp
   use cgnm_types, only : cgnm_problem,cgnm_options,cgnm_result
   use cgnm_core, only : cgnm_fit
   use cgnm_postprocess, only : accepted_indices, top_indices
   implicit none
   private
   public :: cgnm_bootstrap, cgnm_ebe
contains
   subroutine cgnm_bootstrap(prob,opt,base,num_bootstrap,bootstrap_type,out,seed)
      type(cgnm_problem), intent(in) :: prob
      type(cgnm_options), intent(in) :: opt
      type(cgnm_result), intent(in) :: base
      integer, intent(in) :: num_bootstrap,bootstrap_type
      type(cgnm_result), intent(out) :: out
      integer, intent(in), optional :: seed
      type(cgnm_options) :: ob
      real(dp), allocatable :: init(:,:),tm(:,:),wm(:,:),residual(:)
      integer, allocatable :: idx(:)
      integer :: i,j,k,best,nacc,sd
      real(dp) :: u
      sd=opt%seed+7919; if (present(seed)) sd=seed
      call random_seed_from_int(sd)
      call accepted_indices(base,idx)
      if (size(idx)<1) call top_indices(base,min(100,size(base%x,1)),idx)
      if (size(idx)<100 .and. size(base%x,1)>=1) call top_indices(base,min(100,size(base%x,1)),idx)
      nacc=size(idx); best=idx(1)
      allocate(init(num_bootstrap,prob%npar),tm(num_bootstrap,prob%nobs_model))
      allocate(wm(num_bootstrap,prob%nobs_total),residual(prob%nobs_model))
      residual=prob%target-base%y(best,1:prob%nobs_model)
      do i=1,num_bootstrap
         call random_number(u); k=1+min(nacc-1,int(u*real(nacc,dp)))
         init(i,:)=base%theta(idx(k),:)
         tm(i,:)=prob%target
         wm(i,:)=1.0_dp
         select case(bootstrap_type)
         case(1)
            do j=1,prob%nobs_model
               call random_number(u); k=1+min(prob%nobs_model-1,int(u*real(prob%nobs_model,dp)))
               tm(i,j)=prob%target(j)+residual(k)
            end do
         case(2)
            wm(i,1:prob%nobs_model)=0.0_dp
            do j=1,prob%nobs_model
               call random_number(u); k=1+min(prob%nobs_model-1,int(u*real(prob%nobs_model,dp)))
               wm(i,k)=wm(i,k)+1.0_dp
            end do
         case(3)
            do j=1,prob%nobs_model; call random_number(wm(i,j)); end do
         end select
      end do
      ob=opt; ob%num_minimizers=num_bootstrap; ob%seed=sd
      if (bootstrap_type==1) then
         call cgnm_fit(prob,ob,out,initial_iterates=init,target_matrix=tm)
      else
         call cgnm_fit(prob,ob,out,initial_iterates=init,weight_matrix=wm)
      end if
   end subroutine cgnm_bootstrap

   subroutine cgnm_ebe(prob,opt,base,individual_indices,num_repeat,ebe_weight,out,seed)
      type(cgnm_problem), intent(in) :: prob
      type(cgnm_options), intent(in) :: opt
      type(cgnm_result), intent(in) :: base
      integer, intent(in) :: individual_indices(:),num_repeat
      real(dp), intent(in) :: ebe_weight
      type(cgnm_result), intent(out) :: out
      integer, intent(in), optional :: seed
      type(cgnm_options) :: oe
      integer, allocatable :: uniq(:),idx(:)
      real(dp), allocatable :: init(:,:),wm(:,:)
      integer :: i,j,r,k,ng,nrun,sd,countg
      if (size(individual_indices)/=prob%nobs_model) then
         out%status=20; out%message='individual_indices must match model observations'; return
      end if
      call unique_int(individual_indices,uniq); ng=size(uniq); nrun=ng*max(1,num_repeat)
      call accepted_indices(base,idx)
      if (size(idx)<1) call top_indices(base,size(base%x,1),idx)
      allocate(init(nrun,prob%npar),wm(nrun,prob%nobs_total)); wm=1.0_dp
      k=0
      do r=1,max(1,num_repeat)
         do i=1,ng
            k=k+1; init(k,:)=base%theta(idx(1+mod(k-1,size(idx))),:)
            countg=count(individual_indices==uniq(i))
            do j=1,prob%nobs_model
               if (individual_indices(j)==uniq(i)) then
                  wm(k,j)=wm(k,j)+ebe_weight*real(prob%nobs_model,dp)/real(countg,dp)
               end if
            end do
            if (prob%nobs_total>prob%nobs_model) wm(k,prob%nobs_model+1:)=ebe_weight
         end do
      end do
      oe=opt; oe%num_minimizers=nrun; sd=opt%seed+1543; if (present(seed)) sd=seed; oe%seed=sd
      call cgnm_fit(prob,oe,out,initial_iterates=init,weight_matrix=wm)
   end subroutine cgnm_ebe

   subroutine random_seed_from_int(seed)
      integer,intent(in)::seed
      integer :: n,i
      integer,allocatable::put(:)
      call random_seed(size=n); allocate(put(n))
      do i=1,n; put(i)=modulo(seed+65537*i,huge(1)-1); if(put(i)<=0) put(i)=i; end do
      call random_seed(put=put)
   end subroutine random_seed_from_int

   subroutine unique_int(x,u)
      integer,intent(in)::x(:)
      integer,allocatable,intent(out)::u(:)
      integer,allocatable::tmp(:)
      integer :: i,k
      allocate(tmp(size(x))); k=0
      do i=1,size(x)
         if (k==0 .or. .not.any(tmp(1:k)==x(i))) then; k=k+1; tmp(k)=x(i); end if
      end do
      allocate(u(k)); u=tmp(1:k)
   end subroutine unique_int
end module cgnm_extensions
