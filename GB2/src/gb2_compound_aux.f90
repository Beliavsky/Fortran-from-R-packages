! SPDX-License-Identifier: GPL-2.0-or-later
module gb2_compound_aux
  use gb2_kinds, only : dp
  use gb2_compound, only : vofp_cgb2
  use gb2_optimizer, only : optimization_result, bfgs_minimize
  implicit none
  private
  public :: pkl_cavgb2, lambda0_cavgb2, loglik_cavgb2, scoreu_cavgb2, scores_cavgb2
  public :: scorez_cavgb2, hess_cavgb2, fit_cavgb2
  real(dp), allocatable, save :: ctx_fac(:,:),ctx_z(:,:),ctx_w(:)
  integer, save :: ctx_i=0,ctx_l1=0
contains
  subroutine pkl_cavgb2(z,lambda,p)
    real(dp), intent(in) :: z(:,:),lambda(:,:)
    real(dp), intent(out) :: p(:,:)
    real(dp), allocatable :: eta(:)
    real(dp) :: m,s
    integer :: n,j,l1
    l1=size(lambda,2)
    if(size(lambda,1)/=size(z,2) .or. any(shape(p)/=[size(z,1),l1+1])) error stop 'pkl_cavgb2: shape mismatch'
    allocate(eta(l1))
    do n=1,size(z,1)
      eta=matmul(z(n,:),lambda)
      m=max(0.0_dp,maxval(eta))
      s=exp(-m)+sum(exp(eta-m))
      do j=1,l1
      p(n,j)=exp(eta(j)-m)/s
      end do
      p(n,l1+1)=exp(-m)/s
    end do
  end subroutine pkl_cavgb2

  subroutine lambda0_cavgb2(pl0,z,lambda0,w)
    real(dp), intent(in) :: pl0(:),z(:,:)
    real(dp), intent(out) :: lambda0(:,:)
    real(dp), intent(in), optional :: w(:)
    real(dp), allocatable :: vl0(:),means(:)
    real(dp) :: sw,wi
    integer :: i,n,j
    if(any(shape(lambda0)/=[size(z,2),size(pl0)-1])) error stop 'lambda0_cavgb2: shape mismatch'
    allocate(vl0(size(pl0)-1),means(size(z,2)))
    call vofp_cgb2(pl0,vl0)
    means=0.0_dp
    sw=0.0_dp
    do n=1,size(z,1)
    wi=1.0_dp
    if(present(w)) wi=w(n)
    sw=sw+wi
    means=means+wi*z(n,:)
    end do
    means=means/sw
    do j=1,size(vl0)
    do i=1,size(means)
    lambda0(i,j)=vl0(j)/(real(size(z,2),dp)*means(i))
    end do
    end do
  end subroutine lambda0_cavgb2

  real(dp) function loglik_cavgb2(fac,z,lambda,w) result(v)
    real(dp), intent(in) :: fac(:,:),z(:,:),lambda(:,:)
    real(dp), intent(in), optional :: w(:)
    real(dp), allocatable :: p(:,:)
    real(dp) :: sw,wi,mix
    integer :: n
    allocate(p(size(fac,1),size(fac,2)))
    call pkl_cavgb2(z,lambda,p)
    v=0.0_dp
    sw=0.0_dp
    do n=1,size(fac,1)
      wi=1.0_dp
      if(present(w)) wi=w(n)
      mix=dot_product(p(n,:),fac(n,:))
      v=v+wi*log(mix)
      sw=sw+wi
    end do
    v=v/sw
  end function loglik_cavgb2

  subroutine scoreu_cavgb2(fac,z,lambda,u,p)
    real(dp), intent(in) :: fac(:,:),z(:,:),lambda(:,:)
    real(dp), intent(out) :: u(:,:)
    real(dp), intent(out), optional :: p(:,:)
    real(dp), allocatable :: pp(:,:)
    real(dp) :: denom
    integer :: n,j,l1
    l1=size(lambda,2)
    allocate(pp(size(fac,1),l1+1))
    call pkl_cavgb2(z,lambda,pp)
    do n=1,size(fac,1)
    denom=dot_product(pp(n,:),fac(n,:))
    do j=1,l1
    u(n,j)=pp(n,j)*(fac(n,j)/denom-1.0_dp)
    end do
    end do
    if(present(p)) p=pp
  end subroutine scoreu_cavgb2

  subroutine scores_cavgb2(fac,z,lambda,g,w)
    real(dp), intent(in) :: fac(:,:),z(:,:),lambda(:,:)
    real(dp), intent(out) :: g(:,:)
    real(dp), intent(in), optional :: w(:)
    real(dp), allocatable :: u(:,:)
    real(dp) :: wi,sw
    integer :: n,j
    allocate(u(size(fac,1),size(lambda,2)))
    call scoreu_cavgb2(fac,z,lambda,u)
    g=0.0_dp
    sw=0.0_dp
    do n=1,size(fac,1)
      wi=1.0_dp
      if(present(w)) wi=w(n)
      sw=sw+wi
      do j=1,size(lambda,2)
        g(:,j)=g(:,j)+wi*u(n,j)*z(n,:)
      end do
    end do
    g=g/sw
  end subroutine scores_cavgb2

  subroutine scorez_cavgb2(u,z,sc)
    real(dp), intent(in) :: u(:,:),z(:,:)
    real(dp), intent(out) :: sc(:,:)
    integer :: n,j,i,k
    if(any(shape(sc)/=[size(z,1),size(z,2)*size(u,2)])) error stop 'scorez_cavgb2: shape mismatch'
    do n=1,size(z,1)
    k=0
    do j=1,size(u,2)
    do i=1,size(z,2)
    k=k+1
    sc(n,k)=z(n,i)*u(n,j)
    end do
    end do
    end do
  end subroutine scorez_cavgb2

  subroutine hess_cavgb2(u,p,z,h,w)
    real(dp), intent(in) :: u(:,:),p(:,:),z(:,:)
    real(dp), intent(out) :: h(:,:)
    real(dp), intent(in), optional :: w(:)
    real(dp) :: a2,wi
    integer :: i,j,r,s,n,ii,jj,naux,l1
    naux=size(z,2)
    l1=size(u,2)
    if(any(shape(h)/=[naux*l1,naux*l1])) error stop 'hess_cavgb2: shape mismatch'
    h=0.0_dp
    do n=1,size(z,1)
      wi=1.0_dp
      if(present(w)) wi=w(n)
      do j=1,l1
      do i=1,l1
        a2=-p(n,i)*u(n,j)-p(n,j)*u(n,i)-u(n,i)*u(n,j)
        if(i==j) a2=a2+u(n,i)
        do s=1,naux
        jj=(j-1)*naux+s
        do r=1,naux
        ii=(i-1)*naux+r
        h(ii,jj)=h(ii,jj)+wi*z(n,r)*z(n,s)*a2
        end do
        end do
      end do
      end do
    end do
  end subroutine hess_cavgb2

  subroutine fit_cavgb2(fac,z,lambda0,result,w,maxiter,tol)
    real(dp), intent(in) :: fac(:,:),z(:,:),lambda0(:,:)
    type(optimization_result), intent(out) :: result
    real(dp), intent(in), optional :: w(:)
    integer, intent(in), optional :: maxiter
    real(dp), intent(in), optional :: tol
    type(optimization_result) :: raw
    integer :: nm
    real(dp) :: tt
    if(allocated(ctx_fac)) deallocate(ctx_fac)
    if(allocated(ctx_z)) deallocate(ctx_z)
    if(allocated(ctx_w)) deallocate(ctx_w)
    allocate(ctx_fac,source=fac)
    allocate(ctx_z,source=z)
    allocate(ctx_w(size(fac,1)))
    ctx_w=1.0_dp
    if(present(w)) ctx_w=w
    ctx_i=size(lambda0,1)
    ctx_l1=size(lambda0,2)
    nm=300
    if(present(maxiter)) nm=maxiter
    tt=1.0e-9_dp
    if(present(tol)) tt=tol
    call bfgs_minimize(aux_obj,aux_grad,reshape(lambda0,[size(lambda0)]),raw,nm,tt)
    allocate(result%par(size(raw%par)))
    result%par=raw%par
    result%value=raw%value
    result%iterations=raw%iterations
    result%converged=raw%converged
  end subroutine fit_cavgb2
  real(dp) function aux_obj(v) result(f)
    real(dp), intent(in) :: v(:)
    real(dp) :: lam(ctx_i,ctx_l1)
    lam=reshape(v,[ctx_i,ctx_l1])
    f=-loglik_cavgb2(ctx_fac,ctx_z,lam,ctx_w)
  end function aux_obj
  subroutine aux_grad(v,g)
    real(dp), intent(in) :: v(:)
    real(dp), intent(out) :: g(:)
    real(dp) :: lam(ctx_i,ctx_l1),gg(ctx_i,ctx_l1)
    lam=reshape(v,[ctx_i,ctx_l1])
    call scores_cavgb2(ctx_fac,ctx_z,lam,gg,ctx_w)
    g=-reshape(gg,[size(g)])
  end subroutine aux_grad
end module gb2_compound_aux
