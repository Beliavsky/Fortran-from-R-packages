module mixsqp_solver
  use mixsqp_kinds, only : dp
  use mixsqp_lapack, only : dpotrf, dpotrs
  use mixsqp_utils, only : mixobjective, compute_grad_hessian
  use mixsqp_em, only : mixem_update
  implicit none
  private
  public :: sqp_core, active_set_qp
contains
  subroutine feasible_stepsize(x,p,jblock,a)
    real(dp), intent(in) :: x(:),p(:)
    integer, intent(out) :: jblock
    real(dp), intent(out) :: a
    integer :: i
    real(dp) :: t
    a=1.0_dp; jblock=0
    do i=1,size(x)
      if (p(i)<0.0_dp) then
        t=-x(i)/p(i)
        if (t<a) then
          a=t; jblock=i
        end if
      end if
    end do
  end subroutine feasible_stepsize

  function init_hessian_correction(H,a0) result(a)
    real(dp), intent(in) :: H(:,:), a0
    real(dp) :: a,d
    integer :: i
    d=huge(1.0_dp)
    do i=1,size(H,1)
      d=min(d,H(i,i))
    end do
    if (d>a0) then
      a=0.0_dp
    else
      a=a0-d
    end if
  end function init_hessian_correction

  subroutine compute_searchdir(H,b,p,ainc)
    real(dp), intent(in) :: H(:,:),b(:),ainc
    real(dp), intent(out) :: p(size(b))
    real(dp), allocatable :: Bmat(:,:),rhs(:,:)
    real(dp) :: a
    integer :: n,i,info
    real(dp), parameter :: a0=1.0e-15_dp, amax=1.0e15_dp
    n=size(b)
    allocate(Bmat(n,n),rhs(n,1))
    a=init_hessian_correction(H,a0)
    do
      Bmat=H
      do i=1,n
        Bmat(i,i)=Bmat(i,i)+a
      end do
      call dpotrf('U',n,Bmat,n,info)
      if (info==0) exit
      if (a*ainc>amax) exit
      if (a<=0.0_dp) then
        a=a0
      else
        a=a*ainc
      end if
    end do
    if (info/=0) then
      p=0.0_dp
      return
    end if
    rhs(:,1)=-b
    call dpotrs('U',n,1,Bmat,n,rhs,n,info)
    if (info==0) then
      p=rhs(:,1)
    else
      p=0.0_dp
    end if
  end subroutine compute_searchdir

  function active_set_qp(H,g,y,maxiter,zero_search,tol,ainc) result(niter)
    real(dp), intent(in) :: H(:,:),g(:),zero_search,tol,ainc
    real(dp), intent(inout) :: y(:)
    integer, intent(in) :: maxiter
    integer :: niter,m,ni,nj,k,jblock,ii,jj
    integer, allocatable :: ia(:),ja(:)
    logical, allocatable :: freec(:)
    logical :: add_work
    real(dp), allocatable :: Hs(:,:),bs(:),ps(:),p(:),b(:)
    real(dp) :: a,minb
    m=size(g)
    allocate(freec(m),p(m),b(m),ia(m),ja(m))
    freec=(y>0.0_dp)
    niter=0
    do while(niter<maxiter)
      niter=niter+1
      ni=0; nj=0
      do ii=1,m
        if (freec(ii)) then
          ni=ni+1; ia(ni)=ii
        else
          nj=nj+1; ja(nj)=ii; y(ii)=0.0_dp
        end if
      end do
      if (ni==0) exit
      allocate(Hs(ni,ni),bs(ni),ps(ni))
      do ii=1,ni
        bs(ii)=g(ia(ii))
        do jj=1,ni
          Hs(ii,jj)=H(ia(ii),ia(jj))
          bs(ii)=bs(ii)+Hs(ii,jj)*y(ia(jj))
        end do
      end do
      p=0.0_dp
      call compute_searchdir(Hs,bs,ps,ainc)
      do ii=1,ni
        p(ia(ii))=ps(ii)
      end do
      deallocate(Hs,bs,ps)
      a=1.0_dp; add_work=.false.
      if (maxval(abs(p))<=zero_search) then
        b=g+matmul(H,y)
        if (nj==0) exit
        minb=huge(1.0_dp); k=0
        do jj=1,nj
          if (b(ja(jj))<minb) then
            minb=b(ja(jj)); k=ja(jj)
          end if
        end do
        if (minb>=-tol) exit
        freec(k)=.true.
      else
        call feasible_stepsize(y,p,jblock,a)
        if (jblock>0 .and. a<1.0_dp .and. ni>1) add_work=.true.
        y=y+a*p
        where(y<0.0_dp) y=0.0_dp
        do jj=1,nj
          y(ja(jj))=0.0_dp
        end do
        if (add_work) then
          freec(jblock)=.false.
          y(jblock)=0.0_dp
        end if
      end if
    end do
  end function active_set_qp

  function line_search(f,L,U,V,w,z,g,x,y,e,usesvd,suffdecr,beta,amin,a,xnew) result(nls)
    real(dp), intent(in) :: f,L(:,:),U(:,:),V(:,:),w(:),z(:),g(:),x(:),y(:),e(:)
    logical, intent(in) :: usesvd
    real(dp), intent(in) :: suffdecr,beta,amin
    real(dp), intent(out) :: a,xnew(:)
    integer :: nls,jblock
    real(dp), allocatable :: p(:),Luse(:,:)
    real(dp) :: afeas,fnew
    p=y-x
    call feasible_stepsize(x,p,jblock,afeas)
    nls=0
    if (afeas<=amin) then
      a=afeas; xnew=afeas*y+(1.0_dp-afeas)*x
      return
    end if
    a=min(1.0_dp,afeas)
    if (usesvd) then
      allocate(Luse(size(U,1),size(V,1)))
      Luse=matmul(U,transpose(V))
    else
      allocate(Luse(size(L,1),size(L,2)))
      Luse=L
    end if
    do
      xnew=a*y+(1.0_dp-a)*x
      fnew=mixobjective(Luse,xnew,w,z,e)
      nls=nls+1
      if (minval(xnew)>=0.0_dp .and. &
          fnew+sum(xnew)<=f+sum(x)+suffdecr*a*dot_product(y-x,g+1.0_dp)) exit
      if (a*beta<amin) then
        a=amin; xnew=a*y+(1.0_dp-a)*x
        if (minval(xnew)<0.0_dp) then
          a=0.0_dp; xnew=x
        end if
        exit
      end if
      a=a*beta
    end do
  end function line_search

  subroutine sqp_core(L,U,V,w,z,x,usesvd,runem,convtol_sqp,convtol_as, &
      zero_solution,zero_search,suffdecr,beta,amin,ainc,e,maxiter_sqp, &
      maxiter_as,status,obj,rdual,nnz,stepsize,dmax,nqp,nls,nout)
    real(dp), intent(in) :: L(:,:),U(:,:),V(:,:),w(:),z(:),e(:)
    real(dp), intent(inout) :: x(:)
    logical, intent(in) :: usesvd,runem
    real(dp), intent(in) :: convtol_sqp,convtol_as,zero_solution,zero_search
    real(dp), intent(in) :: suffdecr,beta,amin,ainc
    integer, intent(in) :: maxiter_sqp,maxiter_as
    integer, intent(out) :: status,nout
    real(dp), intent(out) :: obj(maxiter_sqp),rdual(maxiter_sqp)
    integer, intent(out) :: nnz(maxiter_sqp),nqp(maxiter_sqp),nls(maxiter_sqp)
    real(dp), intent(out) :: stepsize(maxiter_sqp),dmax(maxiter_sqp)
    real(dp), allocatable :: xold(:),g(:),ghat(:),H(:,:),y(:),xnew(:)
    integer :: i
    allocate(xold(size(x)),g(size(x)),ghat(size(x)),H(size(x),size(x)), &
             y(size(x)),xnew(size(x)))
    status=1; nout=0
    stepsize=-1.0_dp; dmax=-1.0_dp; nqp=-1; nls=-1
    do i=1,maxiter_sqp
      xold=x
      if (runem) call mixem_update(L,w,x)
      where(x<=zero_solution) x=0.0_dp
      obj(i)=mixobjective(L,x,w,z,e)
      call compute_grad_hessian(L,w,x,e,g,H,U,V,usesvd)
      if (count(x>0.0_dp)>0) then
        rdual(i)=-(1.0_dp+minval(g,mask=x>0.0_dp))
      else
        rdual(i)=huge(1.0_dp)
      end if
      nnz(i)=count(x>0.0_dp)
      nout=i
      if (-rdual(i)>=-convtol_sqp) then
        status=0
        exit
      end if
      ghat=g-matmul(H,x)+1.0_dp
      y=x
      nqp(i)=active_set_qp(H,ghat,y,maxiter_as,zero_search,convtol_as,ainc)
      nls(i)=line_search(obj(i),L,U,V,w,z,g,x,y,e,usesvd,suffdecr,beta,amin,stepsize(i),xnew)
      dmax(i)=maxval(abs(xnew-xold))
      x=xnew
    end do
  end subroutine sqp_core
end module mixsqp_solver
