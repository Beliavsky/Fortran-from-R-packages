! SPDX-License-Identifier: MIT
module cgnm_core
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use cgnm_kinds, only : dp
   use cgnm_types, only : cgnm_problem, cgnm_options, cgnm_result
   use cgnm_utils, only : seed_rng, median_value
   use cgnm_linalg, only : weighted_local_linear, cgnr_reg, cgnr_reg_weight
   use cgnm_kmeans, only : optimal_kmeans_labels
   implicit none
   private
   public :: cgnm_fit, cgnm_internal_to_physical, cgnm_physical_to_internal
   public :: cgnm_evaluate, cgnr_ata_atb, cgnr_ata_atb_reg

contains

   subroutine cgnm_physical_to_internal(prob, theta, z, ierr)
      type(cgnm_problem), intent(in) :: prob
      real(dp), intent(in) :: theta(:)
      real(dp), intent(out) :: z(size(theta))
      integer, intent(out) :: ierr
      integer :: j
      real(dp) :: lo, hi
      ierr=0; z=theta
      do j=1,prob%npar
         if (prob%has_lower(j) .and. prob%has_upper(j)) then
            lo=prob%lower_bound(j); hi=prob%upper_bound(j)
            if (theta(j)<=lo .or. theta(j)>=hi) then; ierr=1; return; end if
            z(j)=log((theta(j)-lo)/(hi-theta(j)))
         else if (prob%has_lower(j)) then
            lo=prob%lower_bound(j)
            if (theta(j)<=lo) then; ierr=1; return; end if
            z(j)=log(theta(j)-lo)
         else if (prob%has_upper(j)) then
            hi=prob%upper_bound(j)
            if (theta(j)>=hi) then; ierr=1; return; end if
            z(j)=log(hi-theta(j))
         end if
      end do
   end subroutine cgnm_physical_to_internal

   subroutine cgnm_internal_to_physical(prob, z, theta)
      type(cgnm_problem), intent(in) :: prob
      real(dp), intent(in) :: z(:)
      real(dp), intent(out) :: theta(size(z))
      integer :: j
      real(dp) :: lo,hi,e
      theta=z
      do j=1,prob%npar
         if (prob%has_lower(j) .and. prob%has_upper(j)) then
            lo=prob%lower_bound(j); hi=prob%upper_bound(j)
            if (z(j)>=0.0_dp) then
               e=exp(-min(z(j),700.0_dp))
               theta(j)=lo+(hi-lo)/(1.0_dp+e)
            else
               e=exp(max(z(j),-700.0_dp))
               theta(j)=lo+(hi-lo)*e/(1.0_dp+e)
            end if
         else if (prob%has_lower(j)) then
            theta(j)=prob%lower_bound(j)+exp(min(z(j),700.0_dp))
         else if (prob%has_upper(j)) then
            theta(j)=prob%upper_bound(j)-exp(min(z(j),700.0_dp))
         end if
      end do
   end subroutine cgnm_internal_to_physical

   subroutine cgnm_evaluate(prob, z, y, ierr)
      type(cgnm_problem), intent(in) :: prob
      real(dp), intent(in) :: z(:)
      real(dp), intent(out) :: y(prob%nobs_total)
      integer, intent(out) :: ierr
      real(dp), allocatable :: theta(:), ym(:)
      integer :: j,k
      allocate(theta(prob%npar),ym(prob%nobs_model))
      call cgnm_internal_to_physical(prob,z,theta)
      call prob%model(theta,ym,ierr)
      if (ierr/=0) return
      if (.not.all(ieee_is_finite(ym))) then; ierr=2; return; end if
      y(1:prob%nobs_model)=ym
      k=prob%nobs_model
      do j=1,prob%npar
         if (abs(prob%mo_weights(j))>0.0_dp) then
            k=k+1; y(k)=z(j)*prob%mo_weights(j)
         end if
      end do
   end subroutine cgnm_evaluate

   subroutine transformed_ranges(prob,zlo,zhi,zmo,ierr)
      type(cgnm_problem), intent(in) :: prob
      real(dp), intent(out) :: zlo(prob%npar), zhi(prob%npar), zmo(prob%npar)
      integer, intent(out) :: ierr
      integer :: e
      call cgnm_physical_to_internal(prob,prob%initial_lower,zlo,e)
      if (e/=0) then; ierr=1; return; end if
      call cgnm_physical_to_internal(prob,prob%initial_upper,zhi,e)
      if (e/=0) then; ierr=2; return; end if
      call cgnm_physical_to_internal(prob,prob%mo_values,zmo,e)
      if (e/=0) then
         zmo=0.5_dp*(zlo+zhi)
      end if
      ierr=0
   end subroutine transformed_ranges

   real(dp) function ssr(y,target,w) result(v)
      real(dp), intent(in) :: y(:),target(:),w(:)
      v=sum((w*(y-target))**2)
   end function ssr

   subroutine cgnm_fit(prob,opt,res,initial_iterates,target_matrix,weight_matrix,algorithm_version)
      type(cgnm_problem), intent(in) :: prob
      type(cgnm_options), intent(in) :: opt
      type(cgnm_result), intent(out) :: res
      real(dp), intent(in), optional :: initial_iterates(:,:)
      real(dp), intent(in), optional :: target_matrix(:,:),weight_matrix(:,:)
      integer, intent(in), optional :: algorithm_version
      real(dp), allocatable :: x(:,:),y(:,:),xnew(:,:),ynew(:,:),dxstep(:,:)
      real(dp), allocatable :: zlo(:),zhi(:),zmo(:),resid(:),residnew(:),lambda(:)
      real(dp), allocatable :: target(:,:),weights(:,:), theta(:), tmpz(:), ytmp(:)
      integer :: n,p,m,i,j,k,it,ierr,tries,ver,accepted
      real(dp) :: u,gamma
      logical :: ok

      res%status=0; res%message='success'; res%iterations=0
      if (.not.associated(prob%model)) then
         res%status=1; res%message='model callback is not associated'; return
      end if
      n=opt%num_minimizers; p=prob%npar; m=prob%nobs_total
      if (n<1 .or. p<1 .or. m<1 .or. opt%num_iterations<0) then
         res%status=2; res%message='invalid dimensions/options'; return
      end if
      allocate(zlo(p),zhi(p),zmo(p))
      call transformed_ranges(prob,zlo,zhi,zmo,ierr)
      if (ierr/=0) then
         res%status=3; res%message='initial range is incompatible with bounds'; return
      end if
      if (any(zhi<=zlo)) then
         res%status=4; res%message='initial upper range must exceed lower range'; return
      end if
      allocate(x(n,p),y(n,m),xnew(n,p),ynew(n,m),dxstep(n,p))
      allocate(resid(n),residnew(n),lambda(n),target(n,m),weights(n,m),theta(p),tmpz(p),ytmp(m))
      target=0.0_dp
      target(:,1:prob%nobs_model)=spread(prob%target,1,n)
      k=prob%nobs_model
      do j=1,p
         if (abs(prob%mo_weights(j))>0.0_dp) then
            k=k+1; target(:,k)=zmo(j)*prob%mo_weights(j)
         end if
      end do
      if (present(target_matrix)) then
         if (size(target_matrix,1)==n .and. size(target_matrix,2)==m) then
            target=target_matrix
         else if (size(target_matrix,1)==n .and. size(target_matrix,2)==prob%nobs_model) then
            target(:,1:prob%nobs_model)=target_matrix
         else
            res%status=5; res%message='target_matrix has wrong dimensions'; return
         end if
      end if
      weights=1.0_dp
      if (present(weight_matrix)) then
         if (size(weight_matrix,1)/=n .or. size(weight_matrix,2)/=m) then
            res%status=6; res%message='weight_matrix has wrong dimensions'; return
         end if
         weights=weight_matrix
      end if
      call seed_rng(opt%seed)
      x=0.0_dp; y=0.0_dp
      do i=1,n
         ok=.false.; tries=0
         if (present(initial_iterates)) then
            if (i<=size(initial_iterates,1) .and. size(initial_iterates,2)==p) then
               call cgnm_physical_to_internal(prob,initial_iterates(i,:),tmpz,ierr)
               if (ierr==0) then
                  call cgnm_evaluate(prob,tmpz,ytmp,ierr)
                  if (ierr==0) then; x(i,:)=tmpz; y(i,:)=ytmp; ok=.true.; end if
               end if
            end if
         end if
         do while (.not.ok .and. tries<opt%max_initial_draws)
            do j=1,p
               call random_number(u); tmpz(j)=zlo(j)+u*(zhi(j)-zlo(j))
            end do
            call cgnm_evaluate(prob,tmpz,ytmp,ierr)
            tries=tries+1
            if (ierr==0) then; x(i,:)=tmpz; y(i,:)=ytmp; ok=.true.; end if
         end do
         if (.not.ok) then
            res%status=7; res%message='unable to generate a valid initial cluster'; return
         end if
         resid(i)=ssr(y(i,:),target(i,:),weights(i,:))
      end do
      lambda=opt%initial_lambda
      allocate(res%residual_history(n,opt%num_iterations+1))
      allocate(res%lambda_history(n,opt%num_iterations+1))
      res%residual_history=0.0_dp; res%lambda_history=0.0_dp
      res%residual_history(:,1)=resid; res%lambda_history(:,1)=lambda
      res%initial_x=x; res%initial_y=y
      allocate(res%initial_theta(n,p))
      do i=1,n
         tmpz=x(i,:); call cgnm_internal_to_physical(prob,tmpz,theta); res%initial_theta(i,:)=theta
      end do
      ver=3; if (present(algorithm_version)) ver=algorithm_version
      do it=1,opt%num_iterations
         gamma=opt%gamma
         if (opt%gamma_auto) gamma=2.0_dp**real(nint(real(min(it,50),dp)/15.0_dp),dp)
         call cluster_iteration(x,y,lambda,zlo,zhi,target,weights,gamma,opt,ver, &
                                prob%keep_initial_distribution,dxstep)
         xnew=x+dxstep; ynew=y; residnew=resid
         accepted=0
         do i=1,n
            if (lambda(i)>=opt%initial_lambda*1.0e10_dp) cycle
            tmpz=xnew(i,:); call cgnm_evaluate(prob,tmpz,ytmp,ierr)
            if (ierr==0) then
               ynew(i,:)=ytmp
               residnew(i)=ssr(ynew(i,:),target(i,:),weights(i,:))
            else
               residnew(i)=huge(1.0_dp)
            end if
            if (residnew(i)<resid(i)) then
               x(i,:)=xnew(i,:); y(i,:)=ynew(i,:); resid(i)=residnew(i)
               lambda(i)=lambda(i)/10.0_dp; accepted=accepted+1
            else
               lambda(i)=min(lambda(i)*10.0_dp,opt%initial_lambda*1.0e10_dp)
            end if
         end do
         res%residual_history(:,it+1)=resid; res%lambda_history(:,it+1)=lambda
         res%iterations=it
         if (median_value(resid)<=tiny(1.0_dp)) exit
         if (accepted==0 .and. all(lambda>=opt%initial_lambda*1.0e10_dp)) exit
      end do
      res%x=x; res%y=y; res%target_matrix=target; res%weight_matrix=weights
      allocate(res%theta(n,p))
      do i=1,n
         tmpz=x(i,:); call cgnm_internal_to_physical(prob,tmpz,theta); res%theta(i,:)=theta
      end do
   end subroutine cgnm_fit

   subroutine cluster_iteration(x,y,lambda,zlo,zhi,target,weights,gamma,opt,version,keep_mask,step)
      real(dp), intent(in) :: x(:,:),y(:,:),lambda(:),zlo(:),zhi(:),target(:,:),weights(:,:)
      real(dp), intent(in) :: gamma
      type(cgnm_options), intent(in) :: opt
      integer, intent(in) :: version
      logical, intent(in) :: keep_mask(:)
      real(dp), intent(out) :: step(size(x,1),size(x,2))
      integer, allocatable :: labels(:),members(:)
      real(dp), allocatable :: sdv(:),dx(:,:),dy(:,:),wloc(:),amat(:,:),rhs(:),d(:)
      integer :: n,p,m,i,j,k,c,nc,ierr
      real(dp) :: mu,dist,rec,ratio
      logical, allocatable :: active(:)
      n=size(x,1); p=size(x,2); m=size(y,2); step=0.0_dp
      allocate(labels(n),active(p)); active=.true.
      do j=1,p
         if (maxval(x(:,j))-minval(x(:,j))<=1.0e-14_dp*max(1.0_dp,maxval(abs(x(:,j))))) active(j)=.false.
      end do
      if (version==1) then
         labels=1
      else
         call optimal_kmeans_labels(pack_columns(x,active),labels,opt%kmeans_max_iter)
      end if
      do i=1,n
         if (lambda(i)>=opt%initial_lambda*1.0e10_dp) cycle
         c=labels(i); nc=count(labels==c)
         if (nc<max(2,count(active))) cycle
         allocate(members(nc)); members=pack([(j,j=1,n)],labels==c)
         allocate(sdv(p)); sdv=1.0e-9_dp
         if (version==1) then
            do j=1,p
               if (active(j)) sdv(j)=max(1.0e-9_dp,zhi(j)-zlo(j))
            end do
         else
            do j=1,p
               if (active(j)) then
                  mu=sum(x(members,j))/real(nc,dp)
                  if (nc>1) sdv(j)=max(1.0e-9_dp, &
                     sqrt(sum((x(members,j)-mu)**2)/real(nc-1,dp)))
               end if
            end do
         end if
         allocate(dx(n,p),dy(n,m),wloc(n),amat(m,p),rhs(m),d(p))
         do j=1,n
            dx(j,:)=x(j,:)-x(i,:); dy(j,:)=y(j,:)-y(i,:); wloc(j)=0.0_dp
            if (labels(j)==c .and. j/=i) then
               dist=0.0_dp
               do k=1,p
                  if (active(k)) dist=dist+(dx(j,k)/sdv(k))**2
               end do
               if (dist>tiny(1.0_dp)) then
                  rec=min(1.0e10_dp,1.0_dp/dist)
                  wloc(j)=min(1.0e10_dp,rec**gamma)
               end if
            end if
         end do
         call weighted_local_linear(dx,dy,wloc,amat,ierr)
         if (ierr==0) then
            do j=1,p
               if (.not.active(j) .or. keep_mask(j)) amat(:,j)=0.0_dp
            end do
            rhs=target(i,:)-y(i,:)
            if (any(abs(weights(i,:)-1.0_dp)>0.0_dp)) then
               call cgnr_reg_weight(amat,rhs,lambda(i),weights(i,:),d)
            else
               call cgnr_reg(amat,rhs,lambda(i),d)
            end if
            if (opt%stay_in_initial_range .and. version/=1) then
               do j=1,p
                  if (d(j)<0.0_dp .and. x(i,j)+d(j)<zlo(j)) then
                     ratio=d(j)/(zlo(j)-x(i,j))
                     if (ratio>1.0_dp) d(j)=d(j)/(2.0_dp*ratio)
                  else if (d(j)>0.0_dp .and. x(i,j)+d(j)>zhi(j)) then
                     ratio=d(j)/(zhi(j)-x(i,j))
                     if (ratio>1.0_dp) d(j)=d(j)/(2.0_dp*ratio)
                  end if
               end do
            end if
            step(i,:)=d
         end if
         deallocate(members,sdv,dx,dy,wloc,amat,rhs,d)
      end do
   contains
      function pack_columns(a,mask) result(out)
         real(dp), intent(in) :: a(:,:)
         logical, intent(in) :: mask(:)
         real(dp), allocatable :: out(:,:)
         integer :: q,j,k
         q=count(mask); allocate(out(size(a,1),max(1,q)))
         if (q==0) then; out(:,1)=0.0_dp; return; end if
         k=0
         do j=1,size(a,2)
            if (mask(j)) then; k=k+1; out(:,k)=a(:,j); end if
         end do
      end function pack_columns
   end subroutine cluster_iteration

   subroutine cgnr_ata_atb(a,b,x)
      real(dp), intent(in) :: a(:,:),b(:)
      real(dp), intent(out) :: x(size(a,2))
      call cgnr_reg(a,b,0.0_dp,x)
   end subroutine cgnr_ata_atb

   subroutine cgnr_ata_atb_reg(a,b,lambda,x)
      real(dp), intent(in) :: a(:,:),b(:),lambda
      real(dp), intent(out) :: x(size(a,2))
      call cgnr_reg(a,b,lambda,x)
   end subroutine cgnr_ata_atb_reg
end module cgnm_core
