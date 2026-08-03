! SPDX-License-Identifier: Artistic-2.0
module mts_structural
   use mts_kinds, only : dp
   use mts_types, only : mts_success, mts_invalid_input
   use mts_linalg, only : eye, inverse_matrix, least_squares, symmetric_eigen
   use mts_stats, only : center_columns, covariance_matrix, correlation_matrix, chi_square_survival
   implicit none
   private

   public :: matrix_polynomial_product, matrix_polynomial_product_seasonal
   public :: pi_weight_matrices, apply_matrix_filter
   public :: corner_table, extended_cross_correlation
   public :: approximate_kronecker_indices, kronecker_specification

contains

   subroutine matrix_polynomial_product(a,b,c)
      !! Product C(B)=A(B)B(B).  Coefficient zero is stored in section 1.
      real(dp), intent(in) :: a(:,:,:), b(:,:,:)
      real(dp), allocatable, intent(out) :: c(:,:,:)
      integer :: k,na,nb,i,j
      k=size(a,1);na=size(a,3);nb=size(b,3)
      if(size(a,2)/=k .or. size(b,1)/=k .or. size(b,2)/=k .or. na<1 .or. nb<1)then
         allocate(c(0,0,0));return
      end if
      allocate(c(k,k,na+nb-1));c=0.0_dp
      do i=1,na
         do j=1,nb
            c(:,:,i+j-1)=c(:,:,i+j-1)+matmul(a(:,:,i),b(:,:,j))
         end do
      end do
   end subroutine matrix_polynomial_product

   subroutine matrix_polynomial_product_seasonal(regular,seasonal,season,c)
      !! Product with seasonal coefficients placed at lags season,2*season,...
      real(dp), intent(in) :: regular(:,:,:),seasonal(:,:,:)
      integer,intent(in)::season
      real(dp),allocatable,intent(out)::c(:,:,:)
      real(dp),allocatable::expanded(:,:,:)
      integer::k,j,maxlag
      k=size(regular,1)
      if(season<1.or.size(regular,2)/=k.or.size(seasonal,1)/=k.or.size(seasonal,2)/=k)then
         allocate(c(0,0,0));return
      end if
      maxlag=1+season*(size(seasonal,3)-1)
      allocate(expanded(k,k,maxlag));expanded=0.0_dp
      do j=1,size(seasonal,3)
         expanded(:,:,1+(j-1)*season)=seasonal(:,:,j)
      end do
      call matrix_polynomial_product(regular,expanded,c)
   end subroutine matrix_polynomial_product_seasonal

   subroutine pi_weight_matrices(phi,theta,max_lag,pi,status)
      !! Compute Pi(B)=Theta(B)^(-1) Phi(B), with coefficient zero equal to I.
      !! phi(:,:,j) and theta(:,:,j) are positive-lag matrices in
      !! Phi(B)=I-sum Phi_j B^j and Theta(B)=I-sum Theta_j B^j.
      real(dp),intent(in)::phi(:,:,:),theta(:,:,:)
      integer,intent(in)::max_lag
      real(dp),allocatable,intent(out)::pi(:,:,:)
      integer,intent(out),optional::status
      integer::k,p,q,h,j
      real(dp),allocatable::phicoef(:,:,:)
      k=max(size(phi,1),size(theta,1));p=size(phi,3);q=size(theta,3)
      if(k<1.or.max_lag<0.or.(p>0.and.(size(phi,2)/=k)).or.(q>0.and.(size(theta,2)/=k)))then
         allocate(pi(0,0,0));if(present(status))status=mts_invalid_input;return
      end if
      allocate(pi(k,k,0:max_lag),phicoef(k,k,0:max_lag));pi=0.0_dp;phicoef=0.0_dp
      pi(:,:,0)=eye(k);phicoef(:,:,0)=eye(k)
      do j=1,min(p,max_lag);phicoef(:,:,j)=-phi(:,:,j);end do
      do h=1,max_lag
         pi(:,:,h)=phicoef(:,:,h)
         do j=1,min(q,h)
            pi(:,:,h)=pi(:,:,h)+matmul(theta(:,:,j),pi(:,:,h-j))
         end do
      end do
      if(present(status))status=mts_success
   end subroutine pi_weight_matrices

   subroutine apply_matrix_filter(x,weights,y,constant)
      !! y_t = x_t - c + sum_j W_j x_(t-j).  Weights are packed k by k*p.
      real(dp),intent(in)::x(:,:),weights(:,:)
      real(dp),allocatable,intent(out)::y(:,:)
      real(dp),intent(in),optional::constant(:)
      integer::n,k,p,t,j
      real(dp),allocatable::c(:)
      n=size(x,1);k=size(x,2)
      if(size(weights,1)/=k.or.mod(size(weights,2),max(1,k))/=0)then
         allocate(y(0,0));return
      end if
      p=size(weights,2)/k;allocate(y(n,k),c(k));c=0.0_dp;if(present(constant))c=constant
      do t=1,n
         y(t,:)=x(t,:)-c
         do j=1,min(p,t-1)
            y(t,:)=y(t,:)+matmul(weights(:,(j-1)*k+1:j*k),x(t-j,:))
         end do
      end do
   end subroutine apply_matrix_filter

   subroutine corner_table(y,x,n_rows,n_cols,table,status)
      !! Sample Corner table based on lagged input/output cross-correlations.
      !! Rows are denominator orders and columns are numerator orders.
      real(dp),intent(in)::y(:),x(:)
      integer,intent(in)::n_rows,n_cols
      real(dp),allocatable,intent(out)::table(:,:)
      integer,intent(out),optional::status
      real(dp),allocatable::yc(:),xc(:),design(:,:),target(:,:),beta(:,:),res(:,:),covb(:,:)
      integer::n,r,s,t0,ne,j,istat
      n=min(size(y),size(x))
      if(n_rows<1.or.n_cols<1.or.n<=n_rows+n_cols+2)then
         allocate(table(0,0));if(present(status))status=mts_invalid_input;return
      end if
      allocate(yc(n),xc(n));yc=y(1:n)-sum(y(1:n))/real(n,dp);xc=x(1:n)-sum(x(1:n))/real(n,dp)
      allocate(table(0:n_rows-1,0:n_cols-1));table=0.0_dp
      do r=0,n_rows-1
         do s=0,n_cols-1
            t0=max(r,s)+1;ne=n-t0+1
            allocate(design(ne,1+r+s+1),target(ne,1));design(:,1)=1.0_dp;target(:,1)=yc(t0:n)
            do j=1,r;design(:,1+j)=yc(t0-j:n-j);end do
            do j=0,s;design(:,1+r+j+1)=xc(t0-j:n-j);end do
            call least_squares(design,target,beta,res,covb,istat,ridge=1.0e-10_dp)
            if(istat==mts_success)then
               table(r,s)=sqrt(sum(res(:,1)**2)/real(max(1,ne),dp))
            else
               table(r,s)=huge(1.0_dp)
            end if
            deallocate(design,target,beta,res,covb)
         end do
      end do
      if(present(status))status=mts_success
   end subroutine corner_table

   subroutine extended_cross_correlation(x,max_ar,max_ma,statistics,p_values,status)
      !! Portmanteau-style extended cross-correlation table after VAR prefilters.
      real(dp),intent(in)::x(:,:)
      integer,intent(in)::max_ar,max_ma
      real(dp),allocatable,intent(out)::statistics(:,:),p_values(:,:)
      integer,intent(out),optional::status
      real(dp),allocatable::z(:,:),design(:,:),target(:,:),beta(:,:),res(:,:),covb(:,:),c0(:,:),ch(:,:),c0inv(:,:)
      integer::n,k,p,q,t0,ne,j,h,istat
      n=size(x,1);k=size(x,2)
      if(max_ar<0.or.max_ma<1.or.n<=max_ar+max_ma+2.or.k<1)then
         allocate(statistics(0,0),p_values(0,0));if(present(status))status=mts_invalid_input;return
      end if
      allocate(statistics(0:max_ar,1:max_ma),p_values(0:max_ar,1:max_ma));statistics=0.0_dp;p_values=1.0_dp
      do p=0,max_ar
         t0=p+1;ne=n-p
         if(p==0)then
            z=x-spread(sum(x,dim=1)/real(n,dp),1,n)
         else
            allocate(design(ne,1+k*p),target(ne,k));design(:,1)=1.0_dp;target=x(t0:n,:)
            do j=1,p;design(:,1+(j-1)*k+1:1+j*k)=x(t0-j:n-j,:);end do
            call least_squares(design,target,beta,res,covb,istat,ridge=1.0e-10_dp);z=res
            deallocate(design,target,beta,res,covb)
         end if
         c0=covariance_matrix(z);call inverse_matrix(c0,c0inv,istat)
         if(istat==mts_success)then
            do h=1,max_ma
               if(size(z,1)<=h)cycle
               ch=matmul(transpose(z(h+1:,:)),z(1:size(z,1)-h,:))/real(size(z,1),dp)
               statistics(p,h)=real(size(z,1),dp)*sum(matmul(c0inv,ch)*transpose(matmul(c0inv,ch)))
               statistics(p,h)=max(0.0_dp,statistics(p,h))
               p_values(p,h)=chi_square_survival(statistics(p,h),k*k)
            end do
         end if
      end do
      if(present(status))status=mts_success
   end subroutine extended_cross_correlation

   subroutine approximate_kronecker_indices(x,max_lag,critical_value,indices,canonical_correlations,status)
      !! Approximate Kronecker-index identification from canonical-correlation ranks.
      real(dp),intent(in)::x(:,:)
      integer,intent(in)::max_lag
      real(dp),intent(in),optional::critical_value
      integer,allocatable,intent(out)::indices(:)
      real(dp),allocatable,intent(out)::canonical_correlations(:,:)
      integer,intent(out),optional::status
      real(dp),allocatable::past(:,:),future(:,:),spp(:,:),sff(:,:),spf(:,:),spinv(:,:),sfinv(:,:),m(:,:),vals(:),vecs(:,:)
      real(dp)::crit
      integer::n,k,h,ne,j,istat,r
      n=size(x,1);k=size(x,2);crit=0.05_dp;if(present(critical_value))crit=max(0.0_dp,critical_value)
      if(max_lag<1.or.n<=2*max_lag+2.or.k<1)then
         allocate(indices(0),canonical_correlations(0,0));if(present(status))status=mts_invalid_input;return
      end if
      allocate(indices(k),canonical_correlations(k,max_lag));indices=0;canonical_correlations=0.0_dp
      do h=1,max_lag
         ne=n-2*h+1;allocate(past(ne,k*h),future(ne,k*h))
         do j=1,h
            past(:,(j-1)*k+1:j*k)=x(h+1-j:h+ne-j,:)
            future(:,(j-1)*k+1:j*k)=x(h+j:h+ne+j-1,:)
         end do
         past=center_columns(past);future=center_columns(future)
         spp=matmul(transpose(past),past)/real(ne,dp);sff=matmul(transpose(future),future)/real(ne,dp)
         spf=matmul(transpose(past),future)/real(ne,dp)
         call inverse_matrix(spp+1.0e-10_dp*eye(size(spp,1)),spinv,istat);call inverse_matrix(sff+1.0e-10_dp*eye(size(sff,1)),sfinv,istat)
         m=matmul(spinv,matmul(spf,matmul(sfinv,transpose(spf))))
         call symmetric_eigen(0.5_dp*(m+transpose(m)),vals,vecs,istat)
         do j=1,k;canonical_correlations(j,h)=sqrt(max(0.0_dp,min(1.0_dp,vals(j))));end do
         deallocate(past,future,spp,sff,spf,spinv,sfinv,m,vals,vecs)
      end do
      do j=1,k
         r=0
         do h=1,max_lag
            if(canonical_correlations(j,h)>crit)r=h
         end do
         indices(j)=r
      end do
      if(present(status))status=mts_success
   end subroutine approximate_kronecker_indices

   subroutine kronecker_specification(indices,phi_id,theta_id,status)
      !! Construct the upstream Kronspec 0/1/2 parameter indicator arrays.
      !! Lag zero is represented by the third-index lower bound zero.
      integer,intent(in)::indices(:)
      integer,allocatable,intent(out)::phi_id(:,:,:),theta_id(:,:,:)
      integer,intent(out),optional::status
      integer::k,p,i,j,h,diff
      k=size(indices)
      if(k<1.or.any(indices<0))then
         allocate(phi_id(0,0,0),theta_id(0,0,0));if(present(status))status=mts_invalid_input;return
      end if
      p=maxval(indices);allocate(phi_id(k,k,0:p),theta_id(k,k,0:p));theta_id=2
      do i=1,k
         theta_id(i,i,0)=1
         if(indices(i)<p)theta_id(i,:,indices(i)+1:p)=0
      end do
      if(k>1)then
         do i=1,k-1;theta_id(i,i+1:k,0)=0;end do
      end if
      phi_id=theta_id
      if(k>1)then
         do i=2,k
            do j=1,i-1
               if(indices(j)<=indices(i))phi_id(i,j,0)=0
            end do
         end do
      end if
      theta_id(:,:,0)=phi_id(:,:,0)
      do i=1,k
         do j=1,k
            diff=indices(i)-indices(j)
            if(diff>0)then
               do h=1,min(diff,p);phi_id(i,j,h)=0;end do
            end if
         end do
      end do
      if(present(status))status=mts_success
   end subroutine kronecker_specification

end module mts_structural
