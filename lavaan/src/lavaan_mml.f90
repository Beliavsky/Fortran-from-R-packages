module lavaan_mml
   use lavaan_kinds, only : dp
   use lavaan_ordinal, only : ordinal_thresholds
   use lavaan_linalg, only : sym_eigen_jacobi, inverse_general
   use lavaan_optimizer, only : bfgs_minimize
   use numderiv, only : hessian, nd_success
   implicit none
   private

   type, public :: mml_ordinal_result
      real(dp), allocatable :: loadings(:, :), thresholds(:, :), residual_sd(:)
      real(dp), allocatable :: par(:), se(:), vcov(:, :)
      integer, allocatable :: ncat(:)
      real(dp) :: loglik=-huge(1.0_dp), aic=huge(1.0_dp), bic=huge(1.0_dp)
      integer :: iterations=0, status=0, nquad=0, nfactor=0
      logical :: converged=.false.
   end type mml_ordinal_result

   public :: fit_mml_ordinal_factor, mml_ordinal_loglik, gauss_hermite_normal

contains

   subroutine fit_mml_ordinal_factor(data,loading_start,free_mask,result,nquad)
      integer,intent(in)::data(:,:)
      real(dp),intent(in)::loading_start(:,:)
      logical,intent(in)::free_mask(:,:)
      type(mml_ordinal_result),intent(out)::result
      integer,intent(in),optional::nquad
      real(dp),allocatable::x(:),hess(:,:),hi(:,:)
      integer::p,q,k,i,j,idx,nq,info,status
      real(dp)::fval
      p=size(data,2)
      q=size(loading_start,2)
      nq=9
      if(present(nquad)) nq=max(3,nquad)
      if(size(loading_start,1)/=p .or. any(shape(free_mask)/=shape(loading_start)) .or. q<1 .or. q>3) then
         result%status=-1
         return
      end if
      k=count(free_mask)
      allocate(x(k))
      idx=0
      do j=1,q
      do i=1,p
      if(free_mask(i,j)) then
      idx=idx+1
      x(idx)=loading_start(i,j)
      end if
      end do
      end do
      call build_thresholds(data,result%ncat,result%thresholds,info)
      if(info/=0) then
      result%status=info
      return
      end if
      call bfgs_minimize(nll,x,fval,result%converged,result%iterations,maxiter=1200,tol=2.0e-7_dp)
      result%par=x
      result%loadings=loading_start
      call unpack_loadings(x,result%loadings)
      allocate(result%residual_sd(p))
      do i=1,p
         result%residual_sd(i)=sqrt(max(1.0e-5_dp,1.0_dp-sum(result%loadings(i,:)**2)))
      end do
      result%loglik=-fval
      result%aic=2.0_dp*fval+2.0_dp*real(k,dp)
      result%bic=2.0_dp*fval+log(real(size(data,1),dp))*real(k,dp)
      result%nquad=nq
      result%nfactor=q
      call hessian(nll,x,hess,status=status)
      allocate(result%vcov(k,k),result%se(k))
      result%vcov=0.0_dp
      result%se=huge(1.0_dp)
      if(status==nd_success) then
         call inverse_general(hess,hi,info)
         if(info==0) then
         result%vcov=hi
            do i=1,k
            if(hi(i,i)>=0.0_dp) result%se(i)=sqrt(hi(i,i))
            end do
         end if
      end if
      result%status=0
   contains
      function nll(z) result(v)
         real(dp),intent(in)::z(:)
         real(dp)::v
         real(dp),allocatable::ll(:,:)
         integer::ii
         ll=loading_start
         call unpack_loadings(z,ll)
         do ii=1,p
            if(sum(ll(ii,:)**2)>=0.995_dp) then
            v=huge(1.0_dp)/100.0_dp
            return
            end if
         end do
         v=-mml_ordinal_loglik(data,ll,result%thresholds,result%ncat,nq)
         if(.not.(v<huge(1.0_dp)/10.0_dp)) v=huge(1.0_dp)/100.0_dp
      end function nll
      subroutine unpack_loadings(z,ll)
         real(dp),intent(in)::z(:)
         real(dp),intent(inout)::ll(:,:)
         integer::ii,jj,kk
         kk=0
         do jj=1,q
         do ii=1,p
         if(free_mask(ii,jj)) then
         kk=kk+1
         ll(ii,jj)=z(kk)
         end if
         end do
         end do
      end subroutine unpack_loadings
   end subroutine fit_mml_ordinal_factor

   function mml_ordinal_loglik(data,loadings,thresholds,ncat,nquad) result(ll)
      integer,intent(in)::data(:,:),ncat(:),nquad
      real(dp),intent(in)::loadings(:,:),thresholds(:,:)
      real(dp)::ll
      real(dp),allocatable::nodes(:),weights(:),eta(:),sd(:)
      integer::n,p,q,r,total,code,j,c,iq,tmp,info
      real(dp)::pr,cond,wt,lo,hi,mu
      n=size(data,1)
      p=size(data,2)
      q=size(loadings,2)
      ll=0.0_dp
      if(size(loadings,1)/=p .or. size(ncat)/=p .or. q<1 .or. q>3) then
      ll=-huge(1.0_dp)
      return
      end if
      call gauss_hermite_normal(nquad,nodes,weights,info)
      if(info/=0) then
      ll=-huge(1.0_dp)
      return
      end if
      allocate(eta(q),sd(p))
      do j=1,p
      sd(j)=sqrt(max(1.0e-6_dp,1.0_dp-sum(loadings(j,:)**2)))
      end do
      total=nquad**q
      do r=1,n
         pr=0.0_dp
         do code=0,total-1
            tmp=code
            wt=1.0_dp
            do iq=1,q
               c=mod(tmp,nquad)+1
               tmp=tmp/nquad
               eta(iq)=nodes(c)
               wt=wt*weights(c)
            end do
            cond=1.0_dp
            do j=1,p
               c=data(r,j)
               if(c<1 .or. c>ncat(j)) then
               ll=-huge(1.0_dp)
               return
               end if
               mu=dot_product(loadings(j,:),eta)
               if(c==1) then
               lo=-huge(1.0_dp)
               else
               lo=(thresholds(c-1,j)-mu)/sd(j)
               end if
               if(c==ncat(j)) then
               hi=huge(1.0_dp)
               else
               hi=(thresholds(c,j)-mu)/sd(j)
               end if
               cond=cond*max(normal_cdf(hi)-normal_cdf(lo),1.0e-300_dp)
            end do
            pr=pr+wt*cond
         end do
         ll=ll+log(max(pr,1.0e-300_dp))
      end do
   end function mml_ordinal_loglik

   subroutine gauss_hermite_normal(n,x,w,info)
      integer,intent(in)::n
      real(dp),allocatable,intent(out)::x(:),w(:)
      integer,intent(out)::info
      real(dp),allocatable::jmat(:,:),ev(:),vec(:,:)
      integer::i
      if(n<1) then
      info=-1
      allocate(x(0),w(0))
      return
      end if
      allocate(jmat(n,n))
      jmat=0.0_dp
      do i=1,n-1
      jmat(i,i+1)=sqrt(real(i,dp)/2.0_dp)
      jmat(i+1,i)=jmat(i,i+1)
      end do
      call sym_eigen_jacobi(jmat,ev,vec,info)
      if(info/=0) return
      ! sym_eigen_jacobi sorts ascending; convert physicists' GH to N(0,1) quadrature.
      allocate(x(n),w(n))
      x=sqrt(2.0_dp)*ev
      w=vec(1,:)**2
      w=w/sum(w)
   end subroutine gauss_hermite_normal

   subroutine build_thresholds(data,ncat,thresholds,info)
      integer,intent(in)::data(:,:)
      integer,allocatable,intent(out)::ncat(:)
      real(dp),allocatable,intent(out)::thresholds(:,:)
      integer,intent(out)::info
      integer::p,j,r,maxc
      integer,allocatable::cnt(:)
      real(dp),allocatable::th(:)
      p=size(data,2)
      allocate(ncat(p))
      do j=1,p
      ncat(j)=maxval(data(:,j))
      end do
      if(minval(data)<1 .or. minval(ncat)<2) then
      info=-1
      return
      end if
      maxc=maxval(ncat)
      allocate(thresholds(maxc-1,p))
      thresholds=0.0_dp
      do j=1,p
         allocate(cnt(ncat(j)))
         cnt=0
         do r=1,size(data,1)
         cnt(data(r,j))=cnt(data(r,j))+1
         end do
         if(any(cnt==0)) then
         info=j
         return
         end if
         th=ordinal_thresholds(cnt)
         thresholds(1:size(th),j)=th
         deallocate(cnt)
      end do
      info=0
   end subroutine build_thresholds

   pure function normal_cdf(z) result(p)
      real(dp),intent(in)::z
      real(dp)::p
      if(z>8.0_dp) then
      p=1.0_dp
      else if(z< -8.0_dp) then
      p=0.0_dp
      else
      p=0.5_dp*erfc(-z/sqrt(2.0_dp))
      end if
   end function normal_cdf
end module lavaan_mml
