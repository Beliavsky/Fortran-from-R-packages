module mixsqp_highlevel
  use mixsqp_kinds, only : dp
  use mixsqp_types, only : mixsqp_control, mixsqp_result
  use mixsqp_utils, only : mixobjective, normalize_loglikelihoods, &
    normalize_rows_with_logscale, compute_grad_hessian, truncated_svd
  use mixsqp_em, only : run_mixem
  use mixsqp_solver, only : sqp_core
  implicit none
  private
  public :: fit_mixsqp, mixsqp_default_control
contains
  function mixsqp_default_control() result(c)
    type(mixsqp_control) :: c
    c = mixsqp_control()
  end function mixsqp_default_control

  subroutine fit_mixsqp(Lin,result,w_in,x0_in,log_input,control)
    real(dp), intent(in) :: Lin(:,:)
    type(mixsqp_result), intent(out) :: result
    real(dp), intent(in), optional :: w_in(:), x0_in(:)
    logical, intent(in), optional :: log_input
    type(mixsqp_control), intent(in), optional :: control
    type(mixsqp_control) :: ctl
    real(dp), allocatable :: Lall(:,:),L(:,:),w(:),x0all(:),x(:),z(:),e(:)
    real(dp), allocatable :: U(:,:),V(:,:),obj(:),rdual(:),steps(:),dmax(:)
    real(dp), allocatable :: emobj(:),emdmax(:),g(:),H(:,:),xfinal(:)
    integer, allocatable :: nnz(:),nqp(:),nls(:),emnnz(:),cols(:)
    integer :: n,m,mred,i,j,k,status,nout,maxas,total,rank
    logical :: islog,ok,usesvd
    real(dp) :: sw,sx,emin

    ctl=mixsqp_default_control()
    if (present(control)) ctl=control
    islog=.false.; if (present(log_input)) islog=log_input
    n=size(Lin,1); m=size(Lin,2)
    if (n<1 .or. m<1) error stop "L must have at least one row and one column"
    if (islog) then
      allocate(Lall(n,m)); call normalize_loglikelihoods(Lin,Lall)
    else
      if (any(Lin<0.0_dp)) error stop "L contains negative entries"
      allocate(Lall(n,m)); Lall=Lin
    end if
    if (.not. any(Lall>0.0_dp)) error stop "L has no positive entries"

    allocate(w(n))
    if (present(w_in)) then
      if (size(w_in)/=n .or. any(w_in<0.0_dp)) error stop "invalid weights"
      w=w_in
    else
      w=1.0_dp
    end if
    sw=sum(w); if (sw<=0.0_dp) error stop "weights must have positive sum"
    w=w/sw

    allocate(x0all(m))
    if (present(x0_in)) then
      if (size(x0_in)/=m .or. any(x0_in<0.0_dp)) error stop "invalid x0"
      x0all=x0_in
    else
      x0all=1.0_dp
    end if
    sx=sum(x0all); if (sx<=0.0_dp) error stop "x0 must have positive sum"
    x0all=x0all/sx
    if (mixobjective(Lall,x0all,w)>=huge(1.0_dp)/2.0_dp) error stop "L*x0 must be positive"

    allocate(cols(m)); k=0
    do j=1,m
      if (maxval(Lall(:,j))>0.0_dp) then
        k=k+1; cols(k)=j
      end if
    end do
    if (k==1) then
      allocate(result%x(m),result%grad(m),result%hessian(m,m))
      result%x=0.0_dp; result%x(cols(1))=1.0_dp
      result%value=mixobjective(Lall,result%x,w)
      result%grad=0.0_dp; result%hessian=0.0_dp
      result%status=2; result%status_message="SQP algorithm was not run"
      result%iterations=0
      return
    end if
    mred=k
    allocate(L(n,mred),x(mred))
    do j=1,mred
      L(:,j)=Lall(:,cols(j)); x(j)=x0all(cols(j))
    end do
    x=x/sum(x)

    allocate(z(n)); z=0.0_dp
    if (ctl%normalize_rows .and. .not. islog) call normalize_rows_with_logscale(L,z)
    usesvd=.false.; rank=0
    allocate(U(n,1),V(mred,1)); U=0.0_dp; V=0.0_dp
    if (ctl%tol_svd>0.0_dp .and. mred>4) then
      call truncated_svd(L,ctl%tol_svd,U,V,rank,ok)
      if (ok .and. rank<mred) then
        L=matmul(U,transpose(V)); usesvd=.true.
      else
        if (allocated(U)) deallocate(U)
        if (allocated(V)) deallocate(V)
        allocate(U(n,1),V(mred,1)); U=0.0_dp; V=0.0_dp
      end if
    end if

    allocate(e(n))
    emin=minval(L)
    e=ctl%eps-min(0.0_dp,emin)

    total=ctl%numiter_em+ctl%maxiter_sqp
    allocate(result%objective(total),result%max_rdual(total),result%nnz(total), &
      result%stepsize(total),result%max_diff(total),result%nqp(total),result%nls(total))
    result%objective=0.0_dp; result%max_rdual=huge(1.0_dp); result%nnz=0
    result%stepsize=-1.0_dp; result%max_diff=-1.0_dp; result%nqp=-1; result%nls=-1
    i=0
    if (ctl%numiter_em>0) then
      allocate(emobj(ctl%numiter_em),emdmax(ctl%numiter_em),emnnz(ctl%numiter_em))
      call run_mixem(L,w,z,x,e,ctl%numiter_em,ctl%zero_threshold_solution,emobj,emnnz,emdmax)
      do j=1,ctl%numiter_em
        i=i+1; result%objective(i)=emobj(j); result%nnz(i)=emnnz(j)
        result%stepsize(i)=1.0_dp; result%max_diff(i)=emdmax(j)
      end do
    end if

    maxas=ctl%maxiter_activeset
    if (maxas<=0) maxas=min(20,mred+1)
    allocate(obj(ctl%maxiter_sqp),rdual(ctl%maxiter_sqp),nnz(ctl%maxiter_sqp), &
      steps(ctl%maxiter_sqp),dmax(ctl%maxiter_sqp),nqp(ctl%maxiter_sqp),nls(ctl%maxiter_sqp))
    call sqp_core(L,U,V,w,z,x,usesvd,.true.,ctl%convtol_sqp,ctl%convtol_activeset, &
      ctl%zero_threshold_solution,ctl%zero_threshold_searchdir,ctl%suffdecr_linesearch, &
      ctl%stepsize_reduce,ctl%min_stepsize,ctl%identity_contrib_increase,e,ctl%maxiter_sqp, &
      maxas,status,obj,rdual,nnz,steps,dmax,nqp,nls,nout)
    do j=1,nout
      i=i+1; result%objective(i)=obj(j); result%max_rdual(i)=rdual(j); result%nnz(i)=nnz(j)
      result%stepsize(i)=steps(j); result%max_diff(i)=dmax(j); result%nqp(i)=nqp(j); result%nls(i)=nls(j)
    end do
    result%iterations=i
    if (i<total) then
      call trim_progress(result,i)
    end if

    x=x/sum(x)
    allocate(g(mred),H(mred,mred))
    e=ctl%eps
    call compute_grad_hessian(L,w,x,e,g,H,U,V,usesvd)
    result%value=mixobjective(L,x,w,z)
    allocate(xfinal(m)); xfinal=0.0_dp
    do j=1,mred
      xfinal(cols(j))=x(j)
    end do
    call move_alloc(xfinal,result%x)
    allocate(result%grad(mred),result%hessian(mred,mred))
    result%grad=g; result%hessian=H
    result%status=status
    if (status==0) then
      result%status_message="converged to optimal solution"
    else
      result%status_message="exceeded maximum number of iterations"
    end if
    result%used_svd=usesvd; result%svd_rank=rank
  end subroutine fit_mixsqp

  subroutine trim_progress(r,n)
    type(mixsqp_result), intent(inout) :: r
    integer, intent(in) :: n
    real(dp), allocatable :: a(:)
    integer, allocatable :: ia(:)
    allocate(a(n)); a=r%objective(:n); call move_alloc(a,r%objective)
    allocate(a(n)); a=r%max_rdual(:n); call move_alloc(a,r%max_rdual)
    allocate(ia(n)); ia=r%nnz(:n); call move_alloc(ia,r%nnz)
    allocate(a(n)); a=r%stepsize(:n); call move_alloc(a,r%stepsize)
    allocate(a(n)); a=r%max_diff(:n); call move_alloc(a,r%max_diff)
    allocate(ia(n)); ia=r%nqp(:n); call move_alloc(ia,r%nqp)
    allocate(ia(n)); ia=r%nls(:n); call move_alloc(ia,r%nls)
  end subroutine trim_progress
end module mixsqp_highlevel
