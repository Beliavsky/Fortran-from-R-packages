! SPDX-License-Identifier: Artistic-2.0
module mts_ecm_missing
   use mts_kinds, only : dp
   use mts_types, only : vecm_model, mts_success, mts_invalid_input
   use mts_linalg, only : least_squares, inverse_matrix, matrix_sqrt_symmetric, symmetric_eigen
   use mts_stats, only : determinant_ic
   implicit none
   private

   public :: fit_vecm_known_beta, fit_vecm_johansen
   public :: estimate_missing_observation, estimate_partial_missing

contains

   subroutine fit_vecm_known_beta(x,p,beta,model,include_constant)
      real(dp),intent(in)::x(:,:),beta(:,:)
      integer,intent(in)::p
      type(vecm_model),intent(out)::model
      logical,intent(in),optional::include_constant
      real(dp),allocatable::dx(:,:),design(:,:),y(:,:),coef(:,:),res(:,:),covb(:,:)
      logical::constant
      integer::n,k,r,ne,d,offset,lag,istat,npar
      real(dp)::hqtmp
      n=size(x,1);k=size(x,2);r=size(beta,2);constant=.false.;if(present(include_constant))constant=include_constant
      if(p<1.or.n<=p.or.size(beta,1)/=k.or.r<1)then;model%status=mts_invalid_input;return;end if
      allocate(dx(n-1,k));dx=x(2:n,:)-x(1:n-1,:);ne=n-p
      d=merge(1,0,constant)+r+k*(p-1);allocate(design(ne,d),y(ne,k));y=dx(p:n-1,:);offset=0
      if(constant)then;design(:,1)=1.0_dp;offset=1;end if
      design(:,offset+1:offset+r)=matmul(x(p:n-1,:),beta);offset=offset+r
      do lag=1,p-1
         design(:,offset+(lag-1)*k+1:offset+lag*k)=dx(p-lag:n-1-lag,:)
      end do
      call least_squares(design,y,coef,res,covb,istat)
      model%lag_order=p;model%rank=r;model%include_constant=constant;model%beta=beta
      allocate(model%intercept(k),model%alpha(k,r),model%gamma(k,k,max(0,p-1)))
      model%intercept=0.0_dp;model%gamma=0.0_dp;offset=0
      if(constant)then;model%intercept=coef(1,:);offset=1;end if
      model%alpha=transpose(coef(offset+1:offset+r,:));offset=offset+r
      do lag=1,p-1;model%gamma(:,:,lag)=transpose(coef(offset+(lag-1)*k+1:offset+lag*k,:));end do
      model%residuals=res;model%sigma=matmul(transpose(res),res)/real(ne,dp)
      npar=size(coef);call determinant_ic(model%sigma,ne,npar,model%aic,model%bic,hqtmp,istat)
      model%status=istat
   end subroutine fit_vecm_known_beta

   subroutine fit_vecm_johansen(x,p,rank,model,include_constant)
      real(dp),intent(in)::x(:,:)
      integer,intent(in)::p,rank
      type(vecm_model),intent(out)::model
      logical,intent(in),optional::include_constant
      real(dp),allocatable::dx(:,:),z2(:,:),r0(:,:),r1(:,:),coef(:,:),tmpres(:,:),covb(:,:)
      real(dp),allocatable::s00(:,:),s11(:,:),s01(:,:),s00inv(:,:),s11invroot(:,:),m(:,:),values(:),vectors(:,:),beta(:,:)
      logical::constant
      integer::n,k,ne,d,lag,offset,istat,r
      n=size(x,1);k=size(x,2);r=min(max(1,rank),k);constant=.false.;if(present(include_constant))constant=include_constant
      if(p<1.or.n<=p+2)then;model%status=mts_invalid_input;return;end if
      allocate(dx(n-1,k));dx=x(2:n,:)-x(1:n-1,:);ne=n-p
      d=merge(1,0,constant)+k*(p-1);allocate(z2(ne,d));offset=0
      if(constant)then;z2(:,1)=1.0_dp;offset=1;end if
      do lag=1,p-1;z2(:,offset+(lag-1)*k+1:offset+lag*k)=dx(p-lag:n-1-lag,:);end do
      if(d>0)then
         call least_squares(z2,dx(p:n-1,:),coef,r0,covb,istat)
         call least_squares(z2,x(p:n-1,:),coef,r1,covb,istat)
      else
         r0=dx(p:n-1,:);r1=x(p:n-1,:)
      end if
      s00=matmul(transpose(r0),r0)/real(ne,dp);s11=matmul(transpose(r1),r1)/real(ne,dp);s01=matmul(transpose(r0),r1)/real(ne,dp)
      call inverse_matrix(s00,s00inv,istat);call matrix_sqrt_symmetric(s11,s11invroot,istat,inverse=.true.)
      m=matmul(s11invroot,matmul(transpose(s01),matmul(s00inv,matmul(s01,s11invroot))))
      call symmetric_eigen(m,values,vectors,istat)
      beta=matmul(s11invroot,vectors(:,1:r))
      call normalize_beta(beta)
      call fit_vecm_known_beta(x,p,beta,model,constant)
      model%eigenvalues=values
   end subroutine fit_vecm_johansen

   subroutine normalize_beta(beta)
      real(dp),intent(inout)::beta(:,:)
      real(dp)::scale
      integer::j,imax
      do j=1,size(beta,2)
         imax=maxloc(abs(beta(:,j)),dim=1);scale=beta(imax,j)
         if(abs(scale)>1.0e-14_dp)beta(:,j)=beta(:,j)/scale
      end do
   end subroutine normalize_beta

   subroutine estimate_missing_observation(series,pi_weights,sigma,time_index,estimate,constant,status)
      real(dp),intent(in)::series(:,:),pi_weights(:,:),sigma(:,:)
      integer,intent(in)::time_index
      real(dp),allocatable,intent(out)::estimate(:)
      real(dp),intent(in),optional::constant(:)
      integer,intent(out),optional::status
      logical,allocatable::missing(:)
      allocate(missing(size(series,2)));missing=.true.
      call estimate_partial_missing(series,pi_weights,sigma,time_index,missing,estimate,constant,status)
   end subroutine estimate_missing_observation

   subroutine estimate_partial_missing(series,pi_weights,sigma,time_index,is_missing,estimate,constant,status)
      real(dp),intent(in)::series(:,:),pi_weights(:,:),sigma(:,:)
      integer,intent(in)::time_index
      logical,intent(in)::is_missing(:)
      real(dp),allocatable,intent(out)::estimate(:)
      real(dp),intent(in),optional::constant(:)
      integer,intent(out),optional::status
      real(dp),allocatable::root_inv(:,:),work(:,:),c(:),design(:,:),target(:,:),ainv(:,:),beta(:,:),res(:,:)
      real(dp),allocatable::known(:),row(:),pred(:),block(:,:)
      integer,allocatable::miss_idx(:),obs_idx(:)
      integer::n,k,lags,nm,no,i,j,t,nrows,idx,istat,rstart,rend
      n=size(series,1);k=size(series,2);lags=size(pi_weights,2)/k;nm=count(is_missing);no=k-nm
      if(time_index<1.or.time_index>n.or.size(is_missing)/=k.or.size(pi_weights,1)/=k.or.nm<1)then
         allocate(estimate(0));if(present(status))status=mts_invalid_input;return
      end if
      allocate(miss_idx(nm),obs_idx(no));miss_idx=pack([(i,i=1,k)],is_missing);obs_idx=pack([(i,i=1,k)],.not.is_missing)
      allocate(c(k));c=0.0_dp;if(present(constant))c=constant
      call matrix_sqrt_symmetric(sigma,root_inv,istat,inverse=.true.)
      work=series;work(time_index,miss_idx)=0.0_dp
      nrows=k*(1+min(lags,n-time_index));allocate(design(nrows,nm),target(nrows,1));design=0.0_dp;target=0.0_dp
      pred=c
      do j=1,min(lags,time_index-1)
         block=pi_weights(:,(j-1)*k+1:j*k);pred=pred+matmul(block,work(time_index-j,:))
      end do
      rstart=1;rend=k
      design(rstart:rend,:)=root_inv(:,miss_idx)
      target(rstart:rend,1)=matmul(root_inv,pred)
      if(no>0)target(rstart:rend,1)=target(rstart:rend,1)-matmul(root_inv(:,obs_idx),series(time_index,obs_idx))
      do i=1,min(lags,n-time_index)
         t=time_index+i;pred=series(t,:)-c
         do j=1,min(lags,t-1)
            if(t-j==time_index)cycle
            block=pi_weights(:,(j-1)*k+1:j*k);pred=pred-matmul(block,work(t-j,:))
         end do
         block=pi_weights(:,(i-1)*k+1:i*k)
         rstart=i*k+1;rend=(i+1)*k
         design(rstart:rend,:)=-matmul(root_inv,block(:,miss_idx))
         target(rstart:rend,1)=matmul(root_inv,pred)
         if(no>0)target(rstart:rend,1)=target(rstart:rend,1)+matmul(root_inv,matmul(block(:,obs_idx),series(time_index,obs_idx)))
      end do
      call least_squares(design,target,beta,res,status=istat,ridge=1.0e-10_dp)
      estimate=beta(:,1)
      if(present(status))status=istat
   end subroutine estimate_partial_missing

end module mts_ecm_missing
