! SPDX-License-Identifier: GPL-2.0-only
module kernlab_core
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use kernlab_kinds
  use kernlab_types
  use kernlab_kernels, only: kernel_value, kernel_matrix
  use kernlab_linalg, only: quantile_value, solve_linear, vector_ranks, vec_norm
  implicit none
  private
  public :: sigest, inchol, ipop, couple, couple_vote, couple_pkpd, couple_minpair

contains

  subroutine sigest(x, estimates, status, frac, scaled)
    real(dp),intent(in)::x(:,:)
    real(dp),allocatable,intent(out)::estimates(:)
    integer,intent(out)::status
    real(dp),intent(in),optional::frac
    logical,intent(in),optional::scaled
    real(dp),allocatable::z(:,:),mu(:),sd(:),dist(:)
    real(dp)::f
    integer::n,p,i,j,k,npairs
    logical::doscale
    status=KL_INVALID_ARGUMENT;allocate(estimates(0));n=size(x,1);p=size(x,2)
    if(n<2.or.p<1)return
    doscale=.true.;if(present(scaled))doscale=scaled
    allocate(z(n,p));z=x
    if(doscale)then
      allocate(mu(p),sd(p));mu=sum(x,dim=1)/real(n,dp)
      do j=1,p
        sd(j)=sqrt(sum((x(:,j)-mu(j))**2)/real(max(1,n-1),dp));if(sd(j)<=sqrt(epsilon(1.0_dp)))sd(j)=1.0_dp
        z(:,j)=(x(:,j)-mu(j))/sd(j)
      end do
    end if
    f=0.5_dp;if(present(frac))f=max(0.05_dp,min(1.0_dp,frac))
    npairs=min(n*(n-1)/2,max(10,int(f*real(n*n,dp))))
    allocate(dist(npairs));k=0
    do i=1,n-1
      do j=i+1,n
        if(k>=npairs)exit
        k=k+1;dist(k)=sum((z(i,:)-z(j,:))**2)
      end do
      if(k>=npairs)exit
    end do
    if(k<1)return
    dist=dist(1:k)
    deallocate(estimates)
    allocate(estimates(3))
    estimates(1)=1.0_dp/max(quantile_value(dist,0.9_dp),tiny(1.0_dp))
    estimates(2)=1.0_dp/max(quantile_value(dist,0.5_dp),tiny(1.0_dp))
    estimates(3)=1.0_dp/max(quantile_value(dist,0.1_dp),tiny(1.0_dp))
    status=KL_SUCCESS
  end subroutine sigest

  subroutine inchol(x,kernel,result,tol,maxiter)
    real(dp),intent(in)::x(:,:)
    type(kernel_spec),intent(in)::kernel
    type(inchol_result),intent(out)::result
    real(dp),intent(in),optional::tol
    integer,intent(in),optional::maxiter
    real(dp),allocatable::g(:,:),diagres(:),col(:),temp(:),maxres(:)
    integer,allocatable::piv(:)
    real(dp)::threshold,residue,tau
    integer::n,i,rank,pivot,limit
    result%status=KL_INVALID_ARGUMENT;n=size(x,1);if(n<1)return
    threshold=1.0e-3_dp;if(present(tol))threshold=tol
    limit=n;if(present(maxiter))limit=min(n,maxiter)
    allocate(g(n,limit),diagres(n),piv(limit),maxres(limit),col(n),temp(limit));g=0.0_dp
    do i=1,n;diagres(i)=kernel_value(kernel,x(i,:),x(i,:));end do
    rank=0
    do while(rank<limit)
      pivot=maxloc(diagres,dim=1);residue=diagres(pivot)
      if(residue<=threshold)exit
      rank=rank+1;piv(rank)=pivot;maxres(rank)=residue
      do i=1,n;col(i)=kernel_value(kernel,x(i,:),x(pivot,:));end do
      if(rank>1)then
        temp(1:rank-1)=g(pivot,1:rank-1)
        col=col-matmul(g(:,1:rank-1),temp(1:rank-1))
      end if
      tau=sqrt(max(col(pivot),tiny(1.0_dp)))
      g(:,rank)=col/tau
      diagres=max(0.0_dp,diagres-g(:,rank)**2)
      diagres(pivot)=0.0_dp
    end do
    allocate(result%factor(n,rank),result%pivots(rank),result%diag_residues(n),result%max_residuals(rank))
    result%factor=g(:,1:rank);result%pivots=piv(1:rank);result%diag_residues=diagres
    result%max_residuals=maxres(1:rank);result%rank=rank;result%status=KL_SUCCESS
  end subroutine inchol

  subroutine ipop(c,h,a,b,l,u,r,result,maxiter,tol,penalty)
    real(dp),intent(in)::c(:),h(:,:),a(:,:),b(:),l(:),u(:),r(:)
    type(ipop_result),intent(out)::result
    integer,intent(in),optional::maxiter
    real(dp),intent(in),optional::tol,penalty
    real(dp),allocatable::x(:),grad(:),ax(:),vlow(:),vhigh(:),trial(:),gtrial(:),dual(:)
    real(dp)::rho,eps,step,obj,newobj,gnorm
    integer::n,m,it,limit,bt
    result%status=KL_INVALID_ARGUMENT;n=size(c);m=size(a,1)
    if(size(h,1)/=n.or.size(h,2)/=n.or.size(a,2)/=n.or.size(b)/=m.or.size(r)/=m.or.size(l)/=n.or.size(u)/=n)return
    if(any(l>u).or.any(r<0.0_dp))return
    limit=2000;if(present(maxiter))limit=maxiter;eps=1.0e-7_dp;if(present(tol))eps=tol
    rho=100.0_dp;if(present(penalty))rho=penalty
    allocate(x(n),grad(n),ax(m),vlow(m),vhigh(m),trial(n),gtrial(n),dual(m))
    x=min(u,max(l,0.5_dp*(l+u)))
    do it=1,limit
      ax=matmul(a,x);vlow=max(0.0_dp,b-ax);vhigh=max(0.0_dp,ax-(b+r))
      grad=matmul(h,x)+c+rho*matmul(transpose(a),vhigh-vlow)
      gtrial=x-min(u,max(l,x-grad));gnorm=vec_norm(gtrial)
      obj=0.5_dp*dot_product(x,matmul(h,x))+dot_product(c,x)+0.5_dp*rho*(dot_product(vlow,vlow)+dot_product(vhigh,vhigh))
      if(gnorm<=eps.and.max(maxval(vlow),maxval(vhigh))<=sqrt(eps))exit
      step=1.0_dp
      do bt=1,30
        trial=min(u,max(l,x-step*grad))
        ax=matmul(a,trial);vlow=max(0.0_dp,b-ax);vhigh=max(0.0_dp,ax-(b+r))
        newobj = 0.5_dp*dot_product(trial,matmul(h,trial)) + dot_product(c,trial) &
          + 0.5_dp*rho*(dot_product(vlow,vlow)+dot_product(vhigh,vhigh))
        if(newobj<=obj-1.0e-4_dp*step*dot_product(grad,x-trial))exit
        step=0.5_dp*step
      end do
      x=trial
      if(mod(it,200)==0)rho=min(1.0e8_dp,10.0_dp*rho)
    end do
    ax=matmul(a,x);vlow=max(0.0_dp,b-ax);vhigh=max(0.0_dp,ax-(b+r));dual=rho*(vhigh-vlow)
    allocate(result%primal(n),result%dual(m));result%primal=x;result%dual=dual
    result%objective=0.5_dp*dot_product(x,matmul(h,x))+dot_product(c,x)
    result%primal_infeasibility=max(maxval(vlow),maxval(vhigh));result%iterations=min(it,limit)
    if(result%primal_infeasibility<=1.0e-5_dp)then;result%status=KL_SUCCESS;else;result%status=KL_NOT_CONVERGED;end if
  end subroutine ipop

  subroutine couple(probin, probabilities, status, method)
    real(dp),intent(in)::probin(:)
    real(dp),allocatable,intent(out)::probabilities(:)
    integer,intent(out)::status
    character(len=*),intent(in),optional::method
    character(len=16)::m
    m='minpair';if(present(method))m=adjustl(method)
    select case(trim(m))
    case('vote');call couple_vote(probin,probabilities,status)
    case('pkpd');call couple_pkpd(probin,probabilities,status)
    case default;call couple_minpair(probin,probabilities,status)
    end select
  end subroutine couple

  subroutine pair_matrix(probin,p,status)
    real(dp),intent(in)::probin(:)
    real(dp),allocatable,intent(out)::p(:,:)
    integer,intent(out)::status
    integer::i,j,q,nc
    nc=nint(0.5_dp*(1.0_dp+sqrt(1.0_dp+8.0_dp*real(size(probin),dp))))
    status=KL_INVALID_ARGUMENT;allocate(p(0,0));if(nc*(nc-1)/2/=size(probin))return
    deallocate(p);allocate(p(nc,nc));p=0.0_dp;q=0
    do i=1,nc-1;do j=i+1,nc;q=q+1;p(i,j)=min(1.0_dp,max(0.0_dp,probin(q)));p(j,i)=1.0_dp-p(i,j);end do;end do
    status=KL_SUCCESS
  end subroutine pair_matrix

  subroutine couple_vote(probin,probabilities,status)
    real(dp),intent(in)::probin(:)
    real(dp),allocatable,intent(out)::probabilities(:)
    integer,intent(out)::status
    real(dp),allocatable::p(:,:)
    integer::i,j,n
    call pair_matrix(probin,p,status);if(status/=KL_SUCCESS)then;allocate(probabilities(0));return;end if
    n=size(p,1);allocate(probabilities(n));probabilities=0.0_dp
    do i=1,n-1
      do j=i+1,n
        if(p(i,j)>=0.5_dp)then
          probabilities(i)=probabilities(i)+1.0_dp
        else
          probabilities(j)=probabilities(j)+1.0_dp
        end if
      end do
    end do
    if(sum(probabilities)>0.0_dp)probabilities=probabilities/sum(probabilities)
  end subroutine couple_vote

  subroutine couple_pkpd(probin,probabilities,status)
    real(dp),intent(in)::probin(:)
    real(dp),allocatable,intent(out)::probabilities(:)
    integer,intent(out)::status
    real(dp),allocatable::p(:,:)
    integer::i,j,n
    call pair_matrix(probin,p,status);if(status/=KL_SUCCESS)then;allocate(probabilities(0));return;end if
    n=size(p,1);allocate(probabilities(n))
    do i=1,n
      probabilities(i)=1.0_dp/(sum(1.0_dp/max(p(i,:),1.0e-12_dp),mask=[(j/=i,j=1,n)]) &
        -real(n-2,dp))
    end do
    probabilities=max(probabilities,0.0_dp);probabilities=probabilities/max(sum(probabilities),tiny(1.0_dp))
  end subroutine couple_pkpd

  subroutine couple_minpair(probin,probabilities,status)
    real(dp),intent(in)::probin(:)
    real(dp),allocatable,intent(out)::probabilities(:)
    integer,intent(out)::status
    real(dp),allocatable::p(:,:),q(:,:),rhs(:,:),sol(:,:)
    integer::i,j,n,st
    call pair_matrix(probin,p,status);if(status/=KL_SUCCESS)then;allocate(probabilities(0));return;end if
    n=size(p,1);allocate(q(n+1,n+1),rhs(n+1,1));q=0.0_dp;rhs=0.0_dp
    do i=1,n
      do j=1,n
        if(i==j)cycle
        q(i,i)=q(i,i)+p(j,i)**2
        q(i,j)=q(i,j)-p(i,j)*p(j,i)
      end do
      q(i,n+1)=1.0_dp;q(n+1,i)=1.0_dp
    end do
    rhs(n+1,1)=1.0_dp;call solve_linear(q,rhs,sol,st)
    allocate(probabilities(n))
    if(st==KL_SUCCESS)then;probabilities=max(0.0_dp,sol(1:n,1));else;probabilities=1.0_dp/real(n,dp);end if
    probabilities=probabilities/max(sum(probabilities),tiny(1.0_dp));status=KL_SUCCESS
  end subroutine couple_minpair

end module kernlab_core
