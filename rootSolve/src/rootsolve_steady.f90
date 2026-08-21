! SPDX-License-Identifier: GPL-2.0-or-later
module rootsolve_steady
  use rootsolve_kinds, only : dp
  use rootsolve_types, only : steady_rhs, steady_jac, steady_options, sparse_options, steady_result, runsteady_options
  use rootsolve_derivatives, only : jacobian_full, jacobian_band
  use rootsolve_linalg, only : solve_linear
  use rootsolve_sparse, only : build_grid_pattern, discover_pattern, sparse_jacobian, solve_csr
  use rootsolve_runsteady, only : runsteady
  implicit none
  private
  public :: stode, stodes, steady
contains

  function stode(func, y0, time, options, jacfunc) result(res)
    procedure(steady_rhs) :: func
    real(dp), intent(in) :: y0(:)
    real(dp), intent(in), optional :: time
    type(steady_options), intent(in), optional :: options
    procedure(steady_jac), optional :: jacfunc
    type(steady_result) :: res
    type(steady_options) :: opt
    real(dp), allocatable :: y(:), f(:), jac(:,:), band(:,:), step(:), ewt(:)
    real(dp) :: t
    integer :: n, iter, info, i, j, lo, hi

    if (present(options)) opt=options
    t=0.0_dp
    if (present(time)) t=time
    n=size(y0)
    allocate(y(n),f(n),jac(n,n),step(n),ewt(n),res%precision(opt%maxiter))
    y=y0
    res%precision=0.0_dp

    do iter=1,opt%maxiter
      res%iterations=iter
      call func(t,y,f)
      res%precision(iter)=sum(abs(f))/real(max(1,n),dp)
      call make_weights(y,opt,ewt)
      if(maxval(abs(f)/ewt)<=1.0_dp)then
        res%steady=.true.
        exit
      end if

      select case(trim(opt%jactype))
      case('fullusr')
        if(.not.present(jacfunc))then
          res%status=-3
          exit
        end if
        call jacfunc(t,y,jac)
      case('bandusr')
        if(.not.present(jacfunc))then
          res%status=-3
          exit
        end if
        call jacfunc(t,y,jac)
        do j=1,n
          lo=max(1,j-opt%bandup)
          hi=min(n,j+opt%banddown)
          if(lo>1)jac(:lo-1,j)=0.0_dp
          if(hi<n)jac(hi+1:,j)=0.0_dp
        end do
      case('bandint','1D','1Dint')
        allocate(band(opt%bandup+opt%banddown+1,n))
        call jacobian_band(y,func,opt%bandup,opt%banddown,band,time=t,dy=f,pert=opt%pert)
        jac=0.0_dp
        do j=1,n
          lo=max(1,j-opt%bandup)
          hi=min(n,j+opt%banddown)
          do i=lo,hi
            jac(i,j)=band(i-j+opt%bandup+1,j)
          end do
        end do
        deallocate(band)
      case default
        call jacobian_full(y,func,jac,time=t,dy=f,pert=opt%pert)
      end select

      call solve_linear(jac,-f,step,info)
      if(info/=0)then
        res%status=-2
        exit
      end if
      y=y+step
      call enforce_positive(y,opt)
      if(maxval(abs(step))<=opt%ctol)then
        call func(t,y,f)
        res%steady=.true.
        if(iter<opt%maxiter)then
          res%iterations=iter+1
          res%precision(iter+1)=sum(abs(f))/real(max(1,n),dp)
        end if
        exit
      end if
    end do

    allocate(res%y(n),res%f(n))
    res%y=y
    call func(t,y,res%f)
    res%time=t
    if(res%iterations>0)res%estimated_precision=res%precision(res%iterations)
    if(res%iterations<size(res%precision))res%precision=res%precision(:res%iterations)
    if(res%steady)res%status=1
    if(.not.res%steady.and.res%status==0)res%status=-1
  end function stode

  function stodes(func, y0, time, options) result(res)
    procedure(steady_rhs) :: func
    real(dp), intent(in) :: y0(:)
    real(dp), intent(in), optional :: time
    type(sparse_options), intent(in), optional :: options
    type(steady_result) :: res
    type(sparse_options) :: opt
    real(dp), allocatable :: y(:),f(:),values(:),step(:),ewt(:)
    integer, allocatable :: rowptr(:),colind(:)
    real(dp)::t
    integer::n,iter,info,ndims

    if(present(options))opt=options
    t=0.0_dp
    if(present(time))t=time
    n=size(y0)
    allocate(y(n),f(n),step(n),ewt(n),res%precision(opt%base%maxiter))
    y=y0
    res%precision=0.0_dp
    call func(t,y,f)

    if(allocated(opt%rowptr).and.allocated(opt%colind))then
      rowptr=opt%rowptr
      colind=opt%colind
    else
      select case(trim(opt%sparsetype))
      case('1D')
        ndims=1
        if(opt%dims(1)<=0)opt%dims(1)=n/max(1,opt%nspec)
        call build_grid_pattern(opt%nspec,opt%dims,ndims,opt%cyclic,rowptr,colind)
      case('2D','2Dmap')
        ndims=2
        call build_grid_pattern(opt%nspec,opt%dims,ndims,opt%cyclic,rowptr,colind)
      case('3D','3Dmap')
        ndims=3
        call build_grid_pattern(opt%nspec,opt%dims,ndims,opt%cyclic,rowptr,colind)
      case default
        call discover_pattern(func,t,y,opt%base%pert,opt%drop_tol,rowptr,colind)
      end select
    end if
    if(size(rowptr)/=n+1)error stop 'stodes: sparse pattern dimension mismatch'
    allocate(values(size(colind)))

    do iter=1,opt%base%maxiter
      res%iterations=iter
      call func(t,y,f)
      res%precision(iter)=sum(abs(f))/real(max(1,n),dp)
      call make_weights(y,opt%base,ewt)
      if(maxval(abs(f)/ewt)<=1.0_dp)then
        res%steady=.true.
        exit
      end if
      call sparse_jacobian(func,t,y,f,opt%base%pert,rowptr,colind,values)
      call solve_csr(rowptr,colind,values,-f,step,info)
      if(info/=0)then
        res%status=-2
        exit
      end if
      y=y+step
      call enforce_positive(y,opt%base)
      if(maxval(abs(step))<=opt%base%ctol)then
        call func(t,y,f)
        res%steady=.true.
        if(iter<opt%base%maxiter)then
          res%iterations=iter+1
          res%precision(iter+1)=sum(abs(f))/real(max(1,n),dp)
        end if
        exit
      end if
    end do

    allocate(res%y(n),res%f(n),res%rowptr(size(rowptr)),res%colind(size(colind)))
    res%y=y
    call func(t,y,res%f)
    res%rowptr=rowptr
    res%colind=colind
    res%time=t
    if(res%iterations>0)res%estimated_precision=res%precision(res%iterations)
    if(res%iterations<size(res%precision))res%precision=res%precision(:res%iterations)
    if(res%steady)res%status=1
    if(.not.res%steady.and.res%status==0)res%status=-1
  end function stodes

  function steady(func,y0,time,method,options,sparse_opts,jacfunc) result(res)
    procedure(steady_rhs)::func
    real(dp),intent(in)::y0(:)
    real(dp),intent(in),optional::time
    character(len=*),intent(in),optional::method
    type(steady_options),intent(in),optional::options
    type(sparse_options),intent(in),optional::sparse_opts
    procedure(steady_jac),optional::jacfunc
    type(steady_result)::res
    type(runsteady_options)::rop
    character(len=16)::meth
    real(dp)::t0,times2(2)
    meth='stode'
    if(present(method))meth=method
    select case(trim(meth))
    case('runsteady')
      t0=0.0_dp
      if(present(time))t0=time
      times2=[t0,huge(1.0_dp)/1024.0_dp]
      if(present(options))then
        rop%rtol=options%rtol
        rop%atol=options%atol
        rop%positive=options%positive
      end if
      res=runsteady(func,y0,times2,rop)
    case('stodes')
      if(present(sparse_opts))then
        res=stodes(func,y0,time,sparse_opts)
      else
        res=stodes(func,y0,time)
      end if
    case default
      if(present(options))then
        if(present(jacfunc))then
          res=stode(func,y0,time,options,jacfunc)
        else
          res=stode(func,y0,time,options)
        end if
      else
        if(present(jacfunc))then
          res=stode(func,y0,time,jacfunc=jacfunc)
        else
          res=stode(func,y0,time)
        end if
      end if
    end select
  end function steady

  subroutine make_weights(y,opt,ewt)
    real(dp),intent(in)::y(:)
    type(steady_options),intent(in)::opt
    real(dp),intent(out)::ewt(:)
    if(allocated(opt%rtol_vec))then
      if(size(opt%rtol_vec)/=size(y))error stop 'rtol_vec size mismatch'
      ewt=opt%rtol_vec*abs(y)
    else
      ewt=opt%rtol*abs(y)
    end if
    if(allocated(opt%atol_vec))then
      if(size(opt%atol_vec)/=size(y))error stop 'atol_vec size mismatch'
      ewt=ewt+opt%atol_vec
    else
      ewt=ewt+opt%atol
    end if
    ewt=max(ewt,tiny(1.0_dp))
  end subroutine make_weights

  subroutine enforce_positive(y,opt)
    real(dp),intent(inout)::y(:)
    type(steady_options),intent(in)::opt
    integer::i
    if(opt%positive)y=max(y,0.0_dp)
    if(allocated(opt%positive_index))then
      do i=1,size(opt%positive_index)
        if(opt%positive_index(i)>=1.and.opt%positive_index(i)<=size(y))then
          y(opt%positive_index(i))=max(0.0_dp,y(opt%positive_index(i)))
        end if
      end do
    end if
  end subroutine enforce_positive
end module rootsolve_steady
