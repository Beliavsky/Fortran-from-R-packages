! SPDX-License-Identifier: GPL-2.0-or-later
module rootsolve_pde
  use rootsolve_kinds, only : dp
  use rootsolve_types, only : steady_rhs, root_func, steady_options, sparse_options, &
      runsteady_options, steady_result, root_result
  use rootsolve_steady, only : stode, stodes
  use rootsolve_runsteady, only : runsteady
  use rootsolve_derivatives, only : perturb_value
  use rootsolve_linalg, only : solve_linear
  implicit none
  private
  public :: steady_1d, steady_2d, steady_3d, steady_band, multiroot_1d
contains

  function steady_1d(func,y0,nspec,dimens,time,method,bandwidth,cyclic,options) result(res)
    procedure(steady_rhs)::func
    real(dp),intent(in)::y0(:)
    integer,intent(in)::nspec
    integer,intent(in),optional::dimens,bandwidth
    real(dp),intent(in),optional::time
    character(len=*),intent(in),optional::method
    logical,intent(in),optional::cyclic
    type(steady_options),intent(in),optional::options
    type(steady_result)::res
    type(steady_options)::op
    type(sparse_options)::sp
    type(runsteady_options)::rop
    character(len=16)::meth
    integer::nx,bw
    real(dp)::t,times(2)

    if(present(options))op=options
    if(nspec<1.or.mod(size(y0),nspec)/=0)error stop 'steady_1d: invalid nspec'
    nx=size(y0)/nspec
    if(present(dimens))nx=dimens
    if(nx*nspec/=size(y0))error stop 'steady_1d: dimens*nspec mismatch'
    bw=1
    if(present(bandwidth))bw=bandwidth
    t=0.0_dp
    if(present(time))t=time
    meth='stodes'
    if(present(method))meth=method
    select case(trim(meth))
    case('stode')
      op%jactype='bandint'
      op%bandup=nspec*bw
      op%banddown=nspec*bw
      res=stode(func,y0,t,op)
    case('runsteady')
      rop%rtol=op%rtol
      rop%atol=op%atol
      rop%positive=op%positive
      times=[t,huge(1.0_dp)/1024.0_dp]
      res=runsteady(func,y0,times,rop)
    case default
      sp%base=op
      sp%sparsetype='1D'
      sp%nspec=nspec
      sp%dims=[nx,0,0]
      if(present(cyclic))sp%cyclic(1)=cyclic
      res=stodes(func,y0,t,sp)
    end select
  end function steady_1d

  function steady_2d(func,y0,nspec,dims,time,method,cyclic,options) result(res)
    procedure(steady_rhs)::func
    real(dp),intent(in)::y0(:)
    integer,intent(in)::nspec,dims(2)
    real(dp),intent(in),optional::time
    character(len=*),intent(in),optional::method
    logical,intent(in),optional::cyclic(2)
    type(steady_options),intent(in),optional::options
    type(steady_result)::res
    type(steady_options)::op
    type(sparse_options)::sp
    type(runsteady_options)::rop
    character(len=16)::meth
    real(dp)::t,times(2)

    if(present(options))op=options
    if(nspec*product(dims)/=size(y0))error stop 'steady_2d: dimension mismatch'
    t=0.0_dp
    if(present(time))t=time
    meth='stodes'
    if(present(method))meth=method
    if(trim(meth)=='runsteady')then
      rop%rtol=op%rtol
      rop%atol=op%atol
      rop%positive=op%positive
      times=[t,huge(1.0_dp)/1024.0_dp]
      res=runsteady(func,y0,times,rop)
    else
      sp%base=op
      sp%sparsetype='2D'
      sp%nspec=nspec
      sp%dims=[dims(1),dims(2),0]
      if(present(cyclic))sp%cyclic(1:2)=cyclic
      res=stodes(func,y0,t,sp)
    end if
  end function steady_2d

  function steady_3d(func,y0,nspec,dims,time,method,cyclic,options) result(res)
    procedure(steady_rhs)::func
    real(dp),intent(in)::y0(:)
    integer,intent(in)::nspec,dims(3)
    real(dp),intent(in),optional::time
    character(len=*),intent(in),optional::method
    logical,intent(in),optional::cyclic(3)
    type(steady_options),intent(in),optional::options
    type(steady_result)::res
    type(steady_options)::op
    type(sparse_options)::sp
    type(runsteady_options)::rop
    character(len=16)::meth
    real(dp)::t,times(2)

    if(present(options))op=options
    if(nspec*product(dims)/=size(y0))error stop 'steady_3d: dimension mismatch'
    t=0.0_dp
    if(present(time))t=time
    meth='stodes'
    if(present(method))meth=method
    if(trim(meth)=='runsteady')then
      rop%rtol=op%rtol
      rop%atol=op%atol
      rop%positive=op%positive
      times=[t,huge(1.0_dp)/1024.0_dp]
      res=runsteady(func,y0,times,rop)
    else
      sp%base=op
      sp%sparsetype='3D'
      sp%nspec=nspec
      sp%dims=dims
      if(present(cyclic))sp%cyclic=cyclic
      res=stodes(func,y0,t,sp)
    end if
  end function steady_3d

  function steady_band(func,y0,bandup,banddown,time,options) result(res)
    procedure(steady_rhs)::func
    real(dp),intent(in)::y0(:)
    integer,intent(in)::bandup,banddown
    real(dp),intent(in),optional::time
    type(steady_options),intent(in),optional::options
    type(steady_result)::res
    type(steady_options)::op
    if(present(options))op=options
    op%jactype='bandint'
    op%bandup=bandup
    op%banddown=banddown
    res=stode(func,y0,time,op)
  end function steady_band

  function multiroot_1d(f,start,nspec,bandwidth,options) result(out)
    procedure(root_func)::f
    real(dp),intent(in)::start(:)
    integer,intent(in),optional::nspec,bandwidth
    type(steady_options),intent(in),optional::options
    type(root_result)::out
    type(steady_options)::op
    real(dp),allocatable::x(:),fx(:),fp(:),xp(:),jac(:,:),step(:),ewt(:),delta(:)
    integer::n,ns,bw,band,niter,j,k,i,lo,hi,info
    real(dp)::prec

    if(present(options))op=options
    n=size(start)
    ns=1
    if(present(nspec))ns=max(1,nspec)
    bw=1
    if(present(bandwidth))bw=max(1,bandwidth)
    band=min(n,2*ns*bw+1)
    allocate(x(n),fx(n),fp(n),xp(n),jac(n,n),step(n),ewt(n),delta(n))
    allocate(out%root(n),out%f_root(n),out%precision(op%maxiter))
    x=start
    out%precision=0.0_dp
    call f(x,fx)
    do niter=1,op%maxiter
      out%iterations=niter
      prec=sum(abs(fx))/real(max(1,n),dp)
      out%precision(niter)=prec
      call root_weights(x,op,ewt)
      if(maxval(abs(fx)/ewt)<1.0_dp)then
        out%converged=.true.
        exit
      end if
      do i=1,n
        delta(i)=perturb_value(x(i),op%pert)
      end do
      jac=0.0_dp
      do j=1,min(n,band)
        xp=x
        do k=j,n,band
          xp(k)=xp(k)+delta(k)
        end do
        call f(xp,fp)
        do k=j,n,band
          lo=max(1,k-ns*bw)
          hi=min(n,k+ns*bw)
          do i=lo,hi
            jac(i,k)=(fp(i)-fx(i))/delta(k)
          end do
        end do
      end do
      call solve_linear(jac,-fx,step,info)
      if(info/=0)then
        out%status=-2
        exit
      end if
      x=x+step
      call root_positive(x,op)
      call f(x,fx)
      if(maxval(abs(step))<=op%ctol)then
        out%converged=.true.
        exit
      end if
    end do
    out%root=x
    out%f_root=fx
    if(out%iterations>0)out%estimated_precision=out%precision(out%iterations)
    if(out%iterations<size(out%precision))out%precision=out%precision(:out%iterations)
    if(out%converged)out%status=1
    if(.not.out%converged.and.out%status==0)out%status=-1
  end function multiroot_1d

  subroutine root_weights(x,op,ewt)
    real(dp),intent(in)::x(:)
    type(steady_options),intent(in)::op
    real(dp),intent(out)::ewt(:)
    if(allocated(op%rtol_vec))then
      if(size(op%rtol_vec)/=size(x))error stop 'multiroot_1d: rtol size mismatch'
      ewt=op%rtol_vec*abs(x)
    else
      ewt=op%rtol*abs(x)
    end if
    if(allocated(op%atol_vec))then
      if(size(op%atol_vec)/=size(x))error stop 'multiroot_1d: atol size mismatch'
      ewt=ewt+op%atol_vec
    else
      ewt=ewt+op%atol
    end if
    ewt=max(ewt,tiny(1.0_dp))
  end subroutine root_weights

  subroutine root_positive(x,op)
    real(dp),intent(inout)::x(:)
    type(steady_options),intent(in)::op
    integer::i
    if(op%positive)x=max(x,0.0_dp)
    if(allocated(op%positive_index))then
      do i=1,size(op%positive_index)
        if(op%positive_index(i)>=1.and.op%positive_index(i)<=size(x))then
          x(op%positive_index(i))=max(0.0_dp,x(op%positive_index(i)))
        end if
      end do
    end if
  end subroutine root_positive
end module rootsolve_pde
