! SPDX-License-Identifier: GPL-3.0-only
module pa_views
  use pa_kinds, only: dp
  use pa_types, only: entropy_result
  use pa_linalg, only: solve_linear, inverse_matrix, normal_quantile
  use pa_statistics, only: weighted_moments, sample_moments
  implicit none
  private
  public :: black_litterman, entropy_pool, meucci_moments, meucci_ranking
  public :: centroid, ac_ranking

contains

  subroutine black_litterman(mu,sigma,pick,views,post_mu,post_sigma,omega,info)
    real(dp), intent(in) :: mu(:),sigma(:,:),pick(:,:),views(:)
    real(dp), intent(out) :: post_mu(:),post_sigma(:,:)
    real(dp), intent(in), optional :: omega(:,:)
    integer, intent(out), optional :: info
    real(dp), allocatable :: om(:,:), middle(:,:), rhs(:), x(:), tmp(:,:), invmid(:,:)
    integer :: k,n,istat
    n = size(mu)
    k = size(pick,1)
    if (present(info)) info = 0
    if (size(pick,2) /= n .or. size(views) /= k) then
      if (present(info)) info = -1
      return
    end if
    allocate(om(k,k),middle(k,k),rhs(k),x(k),tmp(n,k),invmid(k,k))
    if (present(omega)) then
      om = omega
    else
      om = matmul(matmul(pick,sigma),transpose(pick))
    end if
    middle = matmul(matmul(pick,sigma),transpose(pick))+om
    rhs = views-matmul(pick,mu)
    call solve_linear(middle,rhs,x,istat)
    if (istat /= 0) then
      if (present(info)) info = istat
      return
    end if
    tmp = matmul(sigma,transpose(pick))
    post_mu = mu+matmul(tmp,x)
    call inverse_matrix(middle,invmid,istat)
    if (istat /= 0) then
      if (present(info)) info = istat
      return
    end if
    post_sigma = sigma-matmul(matmul(tmp,invmid),matmul(pick,sigma))
    post_sigma = 0.5_dp*(post_sigma+transpose(post_sigma))
  end subroutine black_litterman

  subroutine meucci_moments(returns,probabilities,mu,sigma)
    real(dp), intent(in) :: returns(:,:),probabilities(:)
    real(dp), intent(out) :: mu(:),sigma(:,:)
    call weighted_moments(returns,probabilities,mu,sigma)
  end subroutine meucci_moments

  subroutine entropy_pool(prior,result,a,b,aeq,beq,tolerance,max_iterations)
    real(dp), intent(in) :: prior(:)
    type(entropy_result), intent(out) :: result
    real(dp), intent(in), optional :: a(:,:),b(:),aeq(:,:),beq(:)
    real(dp), intent(in), optional :: tolerance
    integer, intent(in), optional :: max_iterations
    real(dp), allocatable :: eqmat(:,:),eqrhs(:),cmat(:,:),crhs(:),lambda(:),q(:)
    real(dp), allocatable :: qtrial(:),res(:),jac(:,:),delta(:),ineq(:),newlambda(:)
    integer, allocatable :: active(:)
    logical, allocatable :: is_active(:)
    real(dp) :: tol,norm0,norm1,step,maxviol,obj
    integer :: j,neq,nineq,neq0,nactive,maxit,outer,iter,istat,idx,m,total_iter
    logical :: has_sum,changed

    j = size(prior)
    tol = 1.0e-10_dp
    if (present(tolerance)) tol = tolerance
    maxit = 200
    if (present(max_iterations)) maxit = max_iterations
    allocate(q(j))
    result%converged = .false.
    result%iterations = 0
    if (j == 0 .or. any(prior <= 0.0_dp) .or. sum(prior) <= 0.0_dp) return

    neq0 = 0
    if (present(aeq)) neq0 = size(aeq,1)
    has_sum = .false.
    if (present(aeq) .and. present(beq)) then
      do m=1,neq0
        if (maxval(abs(aeq(m,:)-1.0_dp)) < 1.0e-12_dp .and. abs(beq(m)-1.0_dp) < 1.0e-12_dp) has_sum=.true.
      end do
    end if
    neq = neq0+merge(0,1,has_sum)
    allocate(eqmat(neq,j),eqrhs(neq))
    m=0
    if (.not. has_sum) then
      m=1
      eqmat(1,:)=1.0_dp
      eqrhs(1)=1.0_dp
    end if
    if (neq0>0) then
      eqmat(m+1:m+neq0,:)=aeq
      eqrhs(m+1:m+neq0)=beq
    end if
    nineq=0
    if (present(a)) nineq=size(a,1)
    allocate(active(max(1,nineq)),is_active(max(1,nineq)))
    active=0
    is_active=.false.
    nactive=0
    total_iter=0

    do outer=1,max(10,nineq*3+5)
      allocate(cmat(neq+nactive,j),crhs(neq+nactive),lambda(neq+nactive), &
               res(neq+nactive),jac(neq+nactive,neq+nactive),delta(neq+nactive), &
               newlambda(neq+nactive),qtrial(j))
      cmat(1:neq,:)=eqmat
      crhs(1:neq)=eqrhs
      do m=1,nactive
        cmat(neq+m,:)=a(active(m),:)
        crhs(neq+m)=b(active(m))
      end do
      lambda=0.0_dp
      do m=1,neq
        if (maxval(abs(cmat(m,:)-1.0_dp))<1.0e-12_dp .and. abs(crhs(m)-1.0_dp)<1.0e-12_dp) lambda(m)=-1.0_dp
      end do

      do iter=1,maxit
        call entropy_probabilities(prior,cmat,lambda,q)
        res=matmul(cmat,q)-crhs
        norm0=maxval(abs(res))
        if (norm0<tol) exit
        call entropy_jacobian(cmat,q,jac)
        call solve_linear(jac,-res,delta,istat)
        if (istat/=0) then
          do m=1,size(jac,1)
            jac(m,m)=jac(m,m)-1.0e-10_dp
          end do
          call solve_linear(jac,-res,delta,istat)
          if (istat/=0) exit
        end if
        step=1.0_dp
        do
          newlambda=lambda+step*delta
          call entropy_probabilities(prior,cmat,newlambda,qtrial)
          norm1=maxval(abs(matmul(cmat,qtrial)-crhs))
          if (norm1<norm0 .or. step<1.0e-8_dp) exit
          step=0.5_dp*step
        end do
        lambda=newlambda
      end do
      total_iter=total_iter+min(iter,maxit)
      call entropy_probabilities(prior,cmat,lambda,q)
      changed=.false.
      if (nactive>0) then
        idx=0
        do m=1,nactive
          if (lambda(neq+m)<-sqrt(tol)) then
            if (idx==0) then
              idx=m
            else if (lambda(neq+m)<lambda(neq+idx)) then
              idx=m
            end if
          end if
        end do
        if (idx>0) then
          is_active(active(idx))=.false.
          if (idx<nactive) active(idx:nactive-1)=active(idx+1:nactive)
          nactive=nactive-1
          changed=.true.
        end if
      end if
      if (.not.changed .and. nineq>0) then
        allocate(ineq(nineq))
        ineq=matmul(a,q)-b
        maxviol=maxval(ineq)
        if (maxviol>tol) then
          idx=maxloc(ineq,dim=1)
          if (.not.is_active(idx)) then
            nactive=nactive+1
            active(nactive)=idx
            is_active(idx)=.true.
            changed=.true.
          end if
        end if
        deallocate(ineq)
      end if
      deallocate(cmat,crhs,res,jac,delta,newlambda,qtrial)
      if (.not.changed) exit
      deallocate(lambda)
    end do

    allocate(result%probabilities(j))
    result%probabilities=q
    if (allocated(lambda)) then
      allocate(result%dual(size(lambda)))
      result%dual=lambda
    end if
    obj=sum(q*log(q/(prior/sum(prior))))
    result%objective=obj
    result%iterations=total_iter
    result%max_violation=maxval(abs(matmul(eqmat,q)-eqrhs))
    if (nineq>0) result%max_violation=max(result%max_violation,max(0.0_dp,maxval(matmul(a,q)-b)))
    result%converged=result%max_violation<sqrt(tol) .and. abs(sum(q)-1.0_dp)<sqrt(tol)
  end subroutine entropy_pool

  subroutine entropy_probabilities(prior,c,lambda,q)
    real(dp),intent(in)::prior(:),c(:,:),lambda(:)
    real(dp),intent(out)::q(:)
    real(dp),allocatable::z(:)
    allocate(z(size(prior)))
    z=-1.0_dp-matmul(transpose(c),lambda)
    z=min(max(z,-700.0_dp),700.0_dp)
    q=(prior/sum(prior))*exp(z)
  end subroutine entropy_probabilities

  subroutine entropy_jacobian(c,q,jac)
    real(dp),intent(in)::c(:,:),q(:)
    real(dp),intent(out)::jac(:,:)
    integer::i
    jac=0.0_dp
    do i=1,size(q)
      jac=jac-q(i)*outer_product(c(:,i),c(:,i))
    end do
  end subroutine entropy_jacobian

  pure function outer_product(x,y) result(a)
    real(dp),intent(in)::x(:),y(:)
    real(dp)::a(size(x),size(y))
    integer::i
    do i=1,size(x)
      a(i,:)=x(i)*y
    end do
  end function outer_product

  subroutine meucci_ranking(returns,prior,order,mu,sigma,posterior,converged)
    real(dp),intent(in)::returns(:,:),prior(:)
    integer,intent(in)::order(:)
    real(dp),intent(out)::mu(:),sigma(:,:),posterior(:)
    logical,intent(out)::converged
    real(dp),allocatable::a(:,:),b(:),aeq(:,:),beq(:)
    type(entropy_result)::ep
    integer::k,i
    k=size(order)
    allocate(a(max(0,k-1),size(returns,1)),b(max(0,k-1)),aeq(1,size(returns,1)),beq(1))
    do i=1,k-1
      a(i,:)=returns(:,order(i))-returns(:,order(i+1))
    end do
    b=0.0_dp
    aeq=1.0_dp
    beq=1.0_dp
    call entropy_pool(prior,ep,a,b,aeq,beq)
    posterior=ep%probabilities
    converged=ep%converged
    call weighted_moments(returns,posterior,mu,sigma)
  end subroutine meucci_ranking

  subroutine centroid(n,values)
    integer,intent(in)::n
    real(dp),intent(out)::values(:)
    real(dp),parameter::aa=0.4424_dp,bb=0.1185_dp,beta=0.21_dp
    real(dp)::alpha,p
    integer::j
    alpha=aa-bb*real(n,dp)**(-beta)
    do j=1,n
      p=(real(n+1-j,dp)-alpha)/(real(n+1,dp)-2.0_dp*alpha)
      values(j)=normal_quantile(p)
    end do
  end subroutine centroid

  subroutine ac_ranking(returns,order,expected_returns)
    real(dp),intent(in)::returns(:,:)
    integer,intent(in)::order(:)
    real(dp),intent(out)::expected_returns(:)
    real(dp),allocatable::c(:)
    real(dp)::xmin,xmax
    integer::n,j
    n=size(order)
    if (size(returns,2) /= n) then
      expected_returns=0.0_dp
      return
    end if
    allocate(c(n))
    call centroid(n,c)
    xmin=minval(c)
    xmax=maxval(c)
    if (xmax>xmin) c=(c-xmin)*0.1_dp/(xmax-xmin)-0.05_dp
    expected_returns=0.0_dp
    do j=1,n
      expected_returns(order(n+1-j))=c(j)
    end do
  end subroutine ac_ranking

end module pa_views
