! SPDX-License-Identifier: GPL-2.0-only
module calibrar_optimization
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use calibrar_kinds, only : dp
  use calibrar_interfaces, only : scalar_objective, gradient_callback, vector_objective
  use calibrar_random, only : rand_uniform, rand_normal, rtnorm_matrix, set_seed
  use calibrar_utils, only : sort_indices, clamp_vector
  use calibrar_stopping, only : stop2 => smooth_stop2, stop3 => smooth_stop3, &
    stop4 => smooth_stop4, n_stop
  implicit none
  private
  public :: optim_options, optim_result, ahr_options, ahr_result
  public :: optim2, optimh, calibrate, calibrate_multi, ahres, ahres_scalar
  public :: minimize_bfgs, minimize_cg, minimize_nelder_mead, minimize_hooke_jeeves
  public :: minimize_spg, minimize_sann

  type :: optim_options
    integer :: maxit = 500
    integer :: maxfeval = 10000
    real(dp) :: reltol = sqrt(epsilon(1.0_dp))
    real(dp) :: abstol = -huge(1.0_dp)
    real(dp) :: fnscale = 1.0_dp
    real(dp) :: initial_step = 1.0_dp
    real(dp) :: temperature = 10.0_dp
    integer :: tmax = 10
    integer :: replicates = 1
    integer :: seed = 12345
    logical :: maximize = .false.
    logical :: trace = .false.
    character(len=16) :: gradient_method = "forward"
  end type optim_options

  type :: optim_result
    real(dp), allocatable :: par(:)
    real(dp) :: value = huge(1.0_dp)
    integer :: function_count = 0
    integer :: gradient_count = 0
    integer :: iterations = 0
    integer :: convergence = 1
    character(len=128) :: message = ""
    real(dp), allocatable :: hessian(:,:)
  end type optim_result

  type :: ahr_options
    integer :: maxgen = 1000
    integer :: popsize = 0
    integer :: replicates = 1
    integer :: seed = 12345
    real(dp) :: selection = 0.5_dp
    real(dp) :: alpha = 0.05_dp
    real(dp) :: step = 0.5_dp
    real(dp) :: reltol = sqrt(epsilon(1.0_dp))
    real(dp) :: convergence = 1.0e-6_dp
    integer :: termination = 2
    integer :: max_no_improvement = 10
    integer :: fn_smoothing = 5
    logical :: use_cv = .true.
    logical :: trace = .false.
  end type ahr_options

  type :: ahr_result
    real(dp), allocatable :: par(:)
    real(dp), allocatable :: partial(:)
    real(dp) :: value = huge(1.0_dp)
    integer :: function_count = 0
    integer :: generations = 0
    integer :: convergence = 1
    character(len=128) :: message = ""
    real(dp), allocatable :: value_trace(:)
    real(dp), allocatable :: best_trace(:)
  end type ahr_result

contains

  subroutine optim2(par, fn, result, method, lower, upper, options, gr, active)
    real(dp), intent(in) :: par(:)
    procedure(scalar_objective) :: fn
    type(optim_result), intent(out) :: result
    character(len=*), intent(in), optional :: method
    real(dp), intent(in), optional :: lower(:), upper(:)
    type(optim_options), intent(in), optional :: options
    procedure(gradient_callback), optional :: gr
    logical, intent(in), optional :: active(:)
    type(optim_options) :: op
    real(dp), allocatable :: lo(:), hi(:)
    character(len=24) :: meth
    integer :: n
    n=size(par); op=optim_options(); if(present(options)) op=options
    allocate(lo(n),hi(n)); lo=-huge(1.0_dp); hi=huge(1.0_dp)
    if(present(lower)) then
      if(size(lower)/=n) error stop "optim2: lower size mismatch"
      lo=lower
    end if
    if(present(upper)) then
      if(size(upper)/=n) error stop "optim2: upper size mismatch"
      hi=upper
    end if
    meth="BFGS"; if(present(method)) meth=adjustl(method)
    select case(trim(meth))
    case("BFGS","Rvmmin","nlm")
      call minimize_bfgs(par,fn,lo,hi,result,op,gr,active,.false.)
    case("L-BFGS-B","LBFGSB3","nlminb","bobyqa")
      call minimize_bfgs(par,fn,lo,hi,result,op,gr,active,.true.)
    case("CG","Rcgmin")
      call minimize_cg(par,fn,lo,hi,result,op,gr,active)
    case("Nelder-Mead","nmk","nmkb")
      call minimize_nelder_mead(par,fn,lo,hi,result,op,active)
    case("hjn","hjk","hjkb","mads")
      call minimize_hooke_jeeves(par,fn,lo,hi,result,op,active)
    case("spg")
      call minimize_spg(par,fn,lo,hi,result,op,gr,active)
    case("SANN","genSA")
      call minimize_sann(par,fn,lo,hi,result,op,active)
    case("AHR-ES")
      call ahres_scalar_as_optim(par,fn,lo,hi,result,op,active)
    case default
      error stop "optim2: method belongs to an external package; use an integration adapter"
    end select
  end subroutine optim2

  subroutine optimh(par, fn, result, method, lower, upper, options, gr, active)
    real(dp), intent(in) :: par(:)
    procedure(scalar_objective) :: fn
    type(optim_result), intent(out) :: result
    character(len=*), intent(in), optional :: method
    real(dp), intent(in), optional :: lower(:), upper(:)
    type(optim_options), intent(in), optional :: options
    procedure(gradient_callback), optional :: gr
    logical, intent(in), optional :: active(:)
    call optim2(par,fn,result,method,lower,upper,options,gr,active)
  end subroutine optimh

  subroutine calibrate(par, fn, phases, result, method, lower, upper, options, gr, replicates)
    real(dp), intent(in) :: par(:)
    procedure(scalar_objective) :: fn
    integer, intent(in), optional :: phases(:)
    type(optim_result), intent(out) :: result
    character(len=*), intent(in), optional :: method
    real(dp), intent(in), optional :: lower(:), upper(:)
    type(optim_options), intent(in), optional :: options
    procedure(gradient_callback), optional :: gr
    integer, intent(in), optional :: replicates(:)
    type(optim_options) :: op
    type(optim_result) :: tmp
    real(dp), allocatable :: x(:), lo(:), hi(:)
    integer, allocatable :: ph(:)
    logical, allocatable :: active(:)
    integer :: n,np,p
    n=size(par); allocate(x(n),lo(n),hi(n),ph(n),active(n)); x=par
    lo=-huge(1.0_dp); hi=huge(1.0_dp)
    if(present(lower)) lo=lower
    if(present(upper)) hi=upper
    if(present(phases)) then
      if(size(phases)/=n) error stop "calibrate: phases size mismatch"
      ph=phases
    else
      ph=1
    end if
    if(maxval(ph)<1 .or. minval(ph,mask=ph>0)/=1) error stop "calibrate: phases must start at 1"
    np=maxval(ph); op=optim_options(); if(present(options)) op=options
    do p=1,np
      active=(ph>0 .and. ph<=p)
      if(present(replicates)) then
        if(size(replicates)==1) then
          op%replicates=replicates(1)
        else if(size(replicates)==np) then
          op%replicates=replicates(p)
        else
          error stop "calibrate: replicates must have length 1 or number of phases"
        end if
      end if
      call optim2(x,fn,tmp,method,lo,hi,op,gr,active)
      x=tmp%par
    end do
    result=tmp
  end subroutine calibrate

  subroutine calibrate_multi(par,fn,nvar,phases,result,lower,upper,options,weights,replicates)
    real(dp), intent(in) :: par(:)
    procedure(vector_objective) :: fn
    integer, intent(in) :: nvar
    integer, intent(in), optional :: phases(:)
    type(ahr_result), intent(out) :: result
    real(dp), intent(in), optional :: lower(:),upper(:),weights(:)
    type(ahr_options), intent(in), optional :: options
    integer, intent(in), optional :: replicates(:)
    type(ahr_options) :: op
    type(ahr_result) :: tmp
    real(dp),allocatable::x(:),lo(:),hi(:),w(:)
    integer,allocatable::ph(:)
    logical,allocatable::active(:)
    integer::n,np,p
    n=size(par)
    allocate(x(n),lo(n),hi(n),ph(n),active(n),w(nvar))
    x=par;lo=-huge(1.0_dp);hi=huge(1.0_dp);w=1.0_dp
    if(present(lower))lo=lower
    if(present(upper))hi=upper
    if(present(weights))w=weights
    if(present(phases))then
      if(size(phases)/=n)error stop "calibrate_multi: phases size mismatch"
      ph=phases
    else
      ph=1
    end if
    if(maxval(ph)<1 .or. minval(ph,mask=ph>0)/=1)error stop "calibrate_multi: phases must start at 1"
    np=maxval(ph);op=ahr_options();if(present(options))op=options
    do p=1,np
      active=(ph>0 .and. ph<=p)
      if(present(replicates))then
        if(size(replicates)==1)then
          op%replicates=replicates(1)
        else if(size(replicates)==np)then
          op%replicates=replicates(p)
        else
          error stop "calibrate_multi: replicates must have length 1 or number of phases"
        end if
      end if
      call ahres(x,fn,nvar,tmp,lo,hi,op,w,active)
      x=tmp%par
    end do
    result=tmp
  end subroutine calibrate_multi

  subroutine minimize_bfgs(par,fn,lower,upper,result,options,gr,active,projected)
    real(dp), intent(in) :: par(:),lower(:),upper(:)
    procedure(scalar_objective) :: fn
    type(optim_result), intent(out) :: result
    type(optim_options), intent(in) :: options
    procedure(gradient_callback), optional :: gr
    logical, intent(in), optional :: active(:)
    logical, intent(in) :: projected
    real(dp), allocatable :: x(:),xn(:),g(:),gn(:),p(:),s(:),y(:),h(:,:),a(:),hy(:)
    real(dp) :: f,fnv,alpha,ys,rho,gnorm
    integer :: n,i,it,ls,fc,gc
    logical :: use_active(size(par))
    n=size(par); allocate(x(n),xn(n),g(n),gn(n),p(n),s(n),y(n),h(n,n),a(n),hy(n))
    use_active=.true.; if(present(active)) use_active=active
    x=par; call clamp_vector(x,lower,upper); h=0.0_dp
    do i=1,n; h(i,i)=1.0_dp; end do
    fc=0; gc=0; f=evaluate(fn,x,options,fc)
    call get_gradient(fn,x,g,options,gr,fc,gc,use_active)
    do it=1,options%maxit
      gnorm=maxval(abs(pack(g,use_active)))
      if(gnorm<=options%reltol) exit
      p=-matmul(h,g); where(.not.use_active) p=0.0_dp
      if(dot_product(p,g)>=0.0_dp) p=-g
      alpha=1.0_dp
      do ls=1,40
        xn=x+alpha*p
        if(projected) call clamp_vector(xn,lower,upper)
        where(.not.use_active) xn=x
        fnv=evaluate(fn,xn,options,fc)
        if(fnv<=f+1.0e-4_dp*dot_product(g,xn-x)) exit
        alpha=0.5_dp*alpha
      end do
      s=xn-x
      if(maxval(abs(s))<=options%reltol*(maxval(abs(x))+options%reltol)) then
        x=xn; f=fnv; exit
      end if
      call get_gradient(fn,xn,gn,options,gr,fc,gc,use_active)
      y=gn-g; ys=dot_product(y,s)
      if(ys>sqrt(epsilon(1.0_dp))*sqrt(max(dot_product(y,y)*dot_product(s,s),tiny(1.0_dp)))) then
        rho=1.0_dp/ys; hy=matmul(h,y)
        h=h-rho*outer(hy,s)-rho*outer(s,hy)+ &
          (rho+rho*rho*dot_product(y,hy))*outer(s,s)
      else
        h=0.0_dp; do i=1,n; h(i,i)=1.0_dp; end do
      end if
      x=xn; f=fnv; g=gn
      if(fc>=options%maxfeval) exit
    end do
    allocate(result%par(n)); result%par=x; result%value=unscale_value(f,options)
    result%function_count=fc; result%gradient_count=gc; result%iterations=min(it,options%maxit)
    if(maxval(abs(pack(g,use_active)))<=options%reltol) then
      result%convergence=0; result%message="converged"
    else
      result%convergence=1; result%message="iteration/evaluation limit reached"
    end if
  end subroutine minimize_bfgs

  subroutine minimize_cg(par,fn,lower,upper,result,options,gr,active)
    real(dp), intent(in) :: par(:),lower(:),upper(:)
    procedure(scalar_objective) :: fn
    type(optim_result), intent(out) :: result
    type(optim_options), intent(in) :: options
    procedure(gradient_callback), optional :: gr
    logical, intent(in), optional :: active(:)
    real(dp), allocatable :: x(:),xn(:),g(:),gn(:),d(:),y(:)
    real(dp) :: f,fnv,alpha,beta,den
    integer :: n,it,ls,fc,gc
    logical :: use_active(size(par))
    n=size(par); allocate(x(n),xn(n),g(n),gn(n),d(n),y(n)); use_active=.true.
    if(present(active)) use_active=active
    x=par; call clamp_vector(x,lower,upper); fc=0; gc=0; f=evaluate(fn,x,options,fc)
    call get_gradient(fn,x,g,options,gr,fc,gc,use_active); d=-g; where(.not.use_active)d=0.0_dp
    do it=1,options%maxit
      if(maxval(abs(pack(g,use_active)))<=options%reltol) exit
      alpha=1.0_dp
      do ls=1,40
        xn=x+alpha*d; call clamp_vector(xn,lower,upper); where(.not.use_active)xn=x
        fnv=evaluate(fn,xn,options,fc)
        if(fnv<=f+1.0e-4_dp*dot_product(g,xn-x)) exit
        alpha=0.5_dp*alpha
      end do
      call get_gradient(fn,xn,gn,options,gr,fc,gc,use_active)
      y=gn-g; den=max(dot_product(g,g),tiny(1.0_dp))
      beta=max(0.0_dp,dot_product(gn,y)/den)
      d=-gn+beta*d; where(.not.use_active)d=0.0_dp
      if(dot_product(d,gn)>=0.0_dp) d=-gn
      x=xn; g=gn; f=fnv
      if(fc>=options%maxfeval) exit
    end do
    allocate(result%par(n)); result%par=x; result%value=unscale_value(f,options)
    result%function_count=fc; result%gradient_count=gc; result%iterations=min(it,options%maxit)
    result%convergence=merge(0,1,maxval(abs(pack(g,use_active)))<=options%reltol)
    if (result%convergence == 0) then
      result%message="converged"
    else
      result%message="iteration/evaluation limit reached"
    end if
  end subroutine minimize_cg

  subroutine minimize_nelder_mead(par,fn,lower,upper,result,options,active)
    real(dp), intent(in) :: par(:),lower(:),upper(:)
    procedure(scalar_objective) :: fn
    type(optim_result), intent(out) :: result
    type(optim_options), intent(in) :: options
    logical, intent(in), optional :: active(:)
    real(dp), allocatable :: simp(:,:),fv(:),cent(:),xr(:),xe(:),xc(:)
    logical :: use_active(size(par))
    real(dp) :: fr,fe,fcv,scale
    integer :: n,i,j,it,fc,ilo,ihi,inh
    n=size(par); allocate(simp(n,n+1),fv(n+1),cent(n),xr(n),xe(n),xc(n)); use_active=.true.
    if(present(active)) use_active=active
    simp(:,1)=par; call clamp_vector(simp(:,1),lower,upper)
    do j=2,n+1
      simp(:,j)=simp(:,1); i=j-1
      if(use_active(i)) then
        scale=0.05_dp*max(abs(par(i)),1.0_dp); simp(i,j)=simp(i,j)+scale
      end if
      call clamp_vector(simp(:,j),lower,upper)
    end do
    fc=0; do j=1,n+1; fv(j)=evaluate(fn,simp(:,j),options,fc); end do
    do it=1,options%maxit
      call simplex_order(fv,ilo,ihi,inh)
      if(maxval(abs(fv-fv(ilo)))<=options%reltol*(abs(fv(ilo))+options%reltol)) exit
      cent=(sum(simp,dim=2)-simp(:,ihi))/real(n,dp)
      xr=cent+(cent-simp(:,ihi)); call clamp_vector(xr,lower,upper); where(.not.use_active)xr=par
      fr=evaluate(fn,xr,options,fc)
      if(fr<fv(ilo)) then
        xe=cent+2.0_dp*(xr-cent); call clamp_vector(xe,lower,upper); where(.not.use_active)xe=par
        fe=evaluate(fn,xe,options,fc)
        if(fe<fr) then; simp(:,ihi)=xe; fv(ihi)=fe; else; simp(:,ihi)=xr; fv(ihi)=fr; end if
      else if(fr<fv(inh)) then
        simp(:,ihi)=xr; fv(ihi)=fr
      else
        xc=cent+0.5_dp*(simp(:,ihi)-cent); call clamp_vector(xc,lower,upper); where(.not.use_active)xc=par
        fcv=evaluate(fn,xc,options,fc)
        if(fcv<fv(ihi)) then
          simp(:,ihi)=xc; fv(ihi)=fcv
        else
          do j=1,n+1
            if(j==ilo) cycle
            simp(:,j)=simp(:,ilo)+0.5_dp*(simp(:,j)-simp(:,ilo)); call clamp_vector(simp(:,j),lower,upper)
            where(.not.use_active)simp(:,j)=par
            fv(j)=evaluate(fn,simp(:,j),options,fc)
          end do
        end if
      end if
      if(fc>=options%maxfeval) exit
    end do
    call simplex_order(fv,ilo,ihi,inh); allocate(result%par(n)); result%par=simp(:,ilo)
    result%value=unscale_value(fv(ilo),options); result%function_count=fc; result%iterations=min(it,options%maxit)
    result%convergence=merge(0,1,it<=options%maxit .and. fc<options%maxfeval); result%message="Nelder-Mead finished"
  end subroutine minimize_nelder_mead

  subroutine minimize_hooke_jeeves(par,fn,lower,upper,result,options,active)
    real(dp), intent(in) :: par(:),lower(:),upper(:)
    procedure(scalar_objective) :: fn
    type(optim_result), intent(out) :: result
    type(optim_options), intent(in) :: options
    logical, intent(in), optional :: active(:)
    real(dp), allocatable :: x(:),xt(:),step(:)
    real(dp) :: f,ft
    logical :: improved,use_active(size(par))
    integer :: n,i,it,fc
    n=size(par); allocate(x(n),xt(n),step(n)); use_active=.true.; if(present(active))use_active=active
    x=par; call clamp_vector(x,lower,upper); step=options%initial_step; where(.not.use_active)step=0.0_dp
    fc=0; f=evaluate(fn,x,options,fc)
    do it=1,options%maxit
      improved=.false.
      do i=1,n
        if(.not.use_active(i)) cycle
        xt=x; xt(i)=min(upper(i),x(i)+step(i)); ft=evaluate(fn,xt,options,fc)
        if(ft<f) then
          x=xt; f=ft; improved=.true.
        else
          xt=x; xt(i)=max(lower(i),x(i)-step(i)); ft=evaluate(fn,xt,options,fc)
          if(ft<f) then; x=xt; f=ft; improved=.true.; end if
        end if
      end do
      if(.not.improved) step=0.5_dp*step
      if(maxval(step)<=options%reltol) exit
      if(fc>=options%maxfeval) exit
    end do
    allocate(result%par(n)); result%par=x; result%value=unscale_value(f,options); result%function_count=fc
    result%iterations=min(it,options%maxit); result%convergence=merge(0,1,maxval(step)<=options%reltol)
    result%message="Hooke-Jeeves finished"
  end subroutine minimize_hooke_jeeves

  subroutine minimize_spg(par,fn,lower,upper,result,options,gr,active)
    real(dp), intent(in) :: par(:),lower(:),upper(:)
    procedure(scalar_objective) :: fn
    type(optim_result), intent(out) :: result
    type(optim_options), intent(in) :: options
    procedure(gradient_callback), optional :: gr
    logical, intent(in), optional :: active(:)
    real(dp), allocatable :: x(:),xn(:),g(:),gn(:),s(:),y(:),d(:)
    real(dp) :: f,fnv,lambda,alpha,sy,ss
    logical :: use_active(size(par))
    integer :: n,it,ls,fc,gc
    n=size(par); allocate(x(n),xn(n),g(n),gn(n),s(n),y(n),d(n)); use_active=.true.
    if(present(active))use_active=active
    x=par; call clamp_vector(x,lower,upper); fc=0;gc=0; f=evaluate(fn,x,options,fc)
    call get_gradient(fn,x,g,options,gr,fc,gc,use_active); lambda=1.0_dp
    do it=1,options%maxit
      xn=x-lambda*g; call clamp_vector(xn,lower,upper); where(.not.use_active)xn=x
      d=xn-x
      if(maxval(abs(d))<=options%reltol) exit
      alpha=1.0_dp
      do ls=1,30
        xn=x+alpha*d; call clamp_vector(xn,lower,upper)
        fnv=evaluate(fn,xn,options,fc)
        if(fnv<=f+1.0e-4_dp*alpha*dot_product(g,d)) exit
        alpha=0.5_dp*alpha
      end do
      call get_gradient(fn,xn,gn,options,gr,fc,gc,use_active)
      s=xn-x; y=gn-g; sy=dot_product(s,y); ss=dot_product(s,s)
      if(sy>0.0_dp) lambda=min(1.0e5_dp,max(1.0e-10_dp,ss/sy))
      x=xn;g=gn;f=fnv
      if(fc>=options%maxfeval)exit
    end do
    allocate(result%par(n)); result%par=x; result%value=unscale_value(f,options); result%function_count=fc
    result%gradient_count=gc
    result%iterations=min(it,options%maxit)
    result%convergence=0
    result%message="SPG finished"
  end subroutine minimize_spg

  subroutine minimize_sann(par,fn,lower,upper,result,options,active)
    real(dp), intent(in) :: par(:),lower(:),upper(:)
    procedure(scalar_objective) :: fn
    type(optim_result), intent(out) :: result
    type(optim_options), intent(in) :: options
    logical, intent(in), optional :: active(:)
    real(dp), allocatable :: x(:),xn(:),best(:)
    real(dp) :: f,fnv,bestf,temp
    logical :: use_active(size(par))
    integer :: n,it,i,fc
    n=size(par);allocate(x(n),xn(n),best(n));use_active=.true.;if(present(active))use_active=active
    call set_seed(options%seed);x=par;call clamp_vector(x,lower,upper);fc=0;f=evaluate(fn,x,options,fc);best=x;bestf=f
    do it=1,options%maxit
      temp=options%temperature/log(real(it+2,dp))
      xn=x
      do i=1,n
        if(use_active(i)) xn(i)=xn(i)+temp*rand_normal()
      end do
      call clamp_vector(xn,lower,upper);fnv=evaluate(fn,xn,options,fc)
      if(fnv<f .or. rand_uniform()<exp(min(0.0_dp,(f-fnv)/max(temp,tiny(1.0_dp))))) then
        x=xn;f=fnv
        if(f<bestf) then;best=x;bestf=f;end if
      end if
      if(fc>=options%maxfeval)exit
    end do
    allocate(result%par(n));result%par=best;result%value=unscale_value(bestf,options);result%function_count=fc
    result%iterations=min(it,options%maxit);result%convergence=0;result%message="SANN finished"
  end subroutine minimize_sann

  subroutine ahres_scalar_as_optim(par,fn,lower,upper,result,options,active)
    real(dp), intent(in) :: par(:),lower(:),upper(:)
    procedure(scalar_objective) :: fn
    type(optim_result), intent(out) :: result
    type(optim_options), intent(in) :: options
    logical, intent(in), optional :: active(:)
    type(ahr_options) :: aop
    type(ahr_result) :: ares
    aop%maxgen=options%maxit
    aop%replicates=options%replicates
    aop%seed=options%seed
    aop%reltol=options%reltol
    call ahres_scalar(par,fn,ares,lower,upper,aop,active)
    allocate(result%par(size(par)))
    result%par=ares%par
    result%value=ares%value
    if(options%maximize) result%value=-result%value
    result%function_count=ares%function_count
    result%iterations=ares%generations
    result%convergence=ares%convergence
    result%message=ares%message
  end subroutine ahres_scalar_as_optim

  subroutine ahres_scalar(par,fn,result,lower,upper,options,active)
    real(dp), intent(in) :: par(:)
    procedure(scalar_objective) :: fn
    type(ahr_result), intent(out) :: result
    real(dp), intent(in), optional :: lower(:),upper(:)
    type(ahr_options), intent(in), optional :: options
    logical, intent(in), optional :: active(:)
    type(ahr_options) :: op
    real(dp),allocatable::lo(:),hi(:),mu(:),sig(:),sd(:),pc(:),ps(:),pop(:,:),fit(:),wrec(:)
    real(dp),allocatable::selected(:,:),wmat(:,:),mueffrow(:),values(:),besttrace(:)
    integer,allocatable::idx(:)
    logical::use_active(size(par)),stopflag
    logical,allocatable::sstop(:)
    real(dp)::step,selection,mueff,cc,chin,cs,damp,mucov,ccov,oldmu(size(par)),bestv
    integer::n,popsize,parents,gen,i,fc,rep
    n=size(par)
    op=ahr_options()
    if(present(options))op=options
    call set_seed(op%seed)
    allocate(lo(n),hi(n),mu(n),sig(n),sd(n),pc(n),ps(n))
    lo=-huge(1.0_dp);hi=huge(1.0_dp)
    if(present(lower))lo=lower
    if(present(upper))hi=upper
    use_active=.true.
    if(present(active))use_active=active
    selection=op%selection
    popsize=op%popsize
    if(popsize<=0)popsize=max(4,int(0.5_dp*real(4+int(3.0_dp*log(real(max(n,2),dp))),dp)/selection))
    parents=max(1,ceiling(selection*real(popsize,dp)))
    allocate(pop(n,popsize),fit(popsize),wrec(parents),selected(n,parents),wmat(n,parents),mueffrow(n))
    allocate(idx(popsize),values(op%maxgen),besttrace(op%maxgen),sstop(op%maxgen))
    mu=par
    sig=(hi-lo)**2/12.0_dp
    where(.not.ieee_is_finite(sig))sig=1.0_dp/12.0_dp
    where(.not.use_active)sig=0.0_dp
    step=op%step
    sd=step*sqrt(max(sig,0.0_dp))
    pc=0.0_dp;ps=0.0_dp
    do i=1,parents
      wrec(i)=log(real(parents+1,dp))-log(real(i,dp))
    end do
    wrec=wrec/sum(wrec)
    mueff=1.0_dp/sum(wrec*wrec)
    cc=4.0_dp/real(n+4,dp)
    chin=sqrt(real(n,dp))*(1.0_dp-1.0_dp/(4.0_dp*n)+1.0_dp/(21.0_dp*n*n))
    fc=0;values=huge(1.0_dp);besttrace=huge(1.0_dp);sstop=.false.
    do gen=1,op%maxgen
      if(all(sd<=tiny(1.0_dp)))exit
      call rtnorm_matrix(popsize,mu,sd,lo,hi,pop)
      pop(:,1)=mu
      do i=1,popsize
        fit(i)=0.0_dp
        do rep=1,max(1,op%replicates)
          fit(i)=fit(i)+eval_scalar(fn,pop(:,i))
          fc=fc+1
        end do
        fit(i)=fit(i)/real(max(1,op%replicates),dp)
      end do
      call sort_indices(fit,idx)
      selected=pop(:,idx(1:parents))
      bestv=fit(idx(1))
      wmat=spread(wrec,1,n)
      mueffrow=1.0_dp/sum(wmat*wmat,dim=2)
      mueff=max(sum(mueffrow)/real(n,dp),1.0_dp)
      cs=(mueff+2.0_dp)/(real(n,dp)+mueff+3.0_dp)
      damp=1.0_dp+2.0_dp*max(0.0_dp,sqrt((mueff-1.0_dp)/real(n+1,dp))-1.0_dp)+cs
      mucov=mueff
      ccov=(1.0_dp/mucov)*(2.0_dp/(real(n,dp)+sqrt(2.0_dp))**2)+ &
        (1.0_dp-1.0_dp/mucov)*min(1.0_dp,(2.0_dp*mueff-1.0_dp)/(real(n+2,dp)**2+mueff))
      oldmu=mu
      mu=sum(wmat*selected,dim=2)
      pc=(1.0_dp-cc)*pc+sqrt(cc*(2.0_dp-cc))*sqrt(mueffrow)*(mu-oldmu)/step
      where(sig>tiny(1.0_dp))
        ps=(1.0_dp-cs)*ps+sqrt(cs*(2.0_dp-cs))*sqrt(mueffrow)*((mu-oldmu)/sqrt(sig))/step
      elsewhere
        ps=0.0_dp
      end where
      sig=(1.0_dp-ccov)*sig+(ccov/mucov)*pc*pc+ &
        ccov*(1.0_dp-1.0_dp/mucov)*sum(wmat*((selected-spread(oldmu,2,parents))/step)**2,dim=2)
      step=step*exp(cs*(sqrt(sum(ps*ps))/chin-1.0_dp)/damp)
      sig=max(sig,0.0_dp)
      where(.not.use_active)sig=0.0_dp
      sd=step*sqrt(sig)
      values(gen)=fit(1)
      besttrace(gen)=bestv
      stopflag=.false.
      select case(op%termination)
      case(0)
      case(1)
        stopflag=step<op%convergence
      case(2)
        sstop(gen)=stop2(values(1:gen),op%reltol,op%fn_smoothing)
        stopflag=n_stop(sstop(1:gen),op%max_no_improvement)
      case(3)
        sstop(gen)=stop3(values(1:gen),op%reltol,op%fn_smoothing)
        stopflag=n_stop(sstop(1:gen),op%max_no_improvement)
      case(4)
        sstop(gen)=stop4(values(1:gen),op%reltol,op%fn_smoothing)
        stopflag=n_stop(sstop(1:gen),op%fn_smoothing)
      end select
      if(stopflag)exit
    end do
    allocate(result%par(n),result%partial(1))
    result%par=mu
    result%partial(1)=eval_scalar(fn,mu)
    fc=fc+1
    result%value=result%partial(1)
    result%function_count=fc
    result%generations=min(gen,op%maxgen)
    result%convergence=merge(0,1,gen<op%maxgen)
    if(result%convergence==0)then
      result%message="stopping criterion reached"
    else
      result%message="maximum generations reached"
    end if
    allocate(result%value_trace(result%generations),result%best_trace(result%generations))
    result%value_trace=values(1:result%generations)
    result%best_trace=besttrace(1:result%generations)
  end subroutine ahres_scalar

  subroutine ahres(par, fn, nvar, result, lower, upper, options, weights, active)
    real(dp), intent(in) :: par(:)
    procedure(vector_objective) :: fn
    integer, intent(in) :: nvar
    type(ahr_result), intent(out) :: result
    real(dp), intent(in), optional :: lower(:), upper(:), weights(:)
    type(ahr_options), intent(in), optional :: options
    logical, intent(in), optional :: active(:)
    type(ahr_options) :: op
    real(dp), allocatable :: lo(:),hi(:),wgt(:),mu(:),sigma_diag(:),sd(:),pc(:),ps(:)
    real(dp), allocatable :: pop(:,:),fit(:,:),fg(:),wrec(:),optind(:,:),optsd(:,:)
    real(dp), allocatable :: mu_state(:,:),sigma_state(:,:),wnew(:,:),wmat(:,:),xcv(:,:)
    real(dp), allocatable :: selected(:,:),values(:),besttrace(:),fit_tmp(:),rangev(:),mu_eff_row(:)
    real(dp), allocatable :: cvmins(:),cvmaxs(:)
    integer, allocatable :: idx(:),sups(:),lidx(:)
    logical :: use_active(size(par))
    real(dp) :: step,selection,mu_eff,cc,chi_n,cs,damp,ccov,mucov
    real(dp) :: bestv,mu_eff_global,old_mu(size(par))
    integer :: n,popsize,parents,gen,i,j,k,fc,rep
    logical :: stopflag
    logical, allocatable :: sstop(:)
    n=size(par);op=ahr_options();if(present(options))op=options;call set_seed(op%seed)
    allocate(lo(n),hi(n),wgt(nvar),mu(n),sigma_diag(n),sd(n),pc(n),ps(n),rangev(n));
    lo=-huge(1.0_dp);hi=huge(1.0_dp);if(present(lower))lo=lower;if(present(upper))hi=upper
    wgt=1.0_dp;if(present(weights))wgt=weights
    use_active=.true.;if(present(active))use_active=active
    selection=op%selection;popsize=op%popsize
    if(popsize<=0)popsize=max(4,int(0.5_dp*real(4+int(3.0_dp*log(real(max(n,2),dp))),dp)/selection))
    parents=max(1,ceiling(selection*real(popsize,dp)))
    allocate(pop(n,popsize),fit(popsize,nvar),fg(popsize),wrec(parents),selected(n,parents))
    allocate(optind(n,nvar),optsd(n,nvar),mu_state(n,nvar),sigma_state(n,nvar))
    allocate(wnew(n,nvar),wmat(n,parents),xcv(n,nvar),mu_eff_row(n))
    allocate(cvmins(nvar),cvmaxs(nvar),sstop(op%maxgen))
    allocate(idx(popsize),sups(parents),lidx(parents),values(op%maxgen),besttrace(op%maxgen),fit_tmp(nvar))
    mu=par
    rangev=hi-lo
    where(.not.ieee_is_finite(rangev)) rangev=1.0_dp
    sigma_diag=rangev*rangev/12.0_dp
    where(.not.use_active) sigma_diag=0.0_dp
    mu_state=spread(mu,2,nvar)
    sigma_state=0.0_dp
    step=op%step;sd=step*sqrt(max(sigma_diag,0.0_dp));pc=0.0_dp;ps=0.0_dp
    do i=1,parents;wrec(i)=log(real(parents+1,dp))-log(real(i,dp));end do;wrec=wrec/sum(wrec)
    mu_eff=1.0_dp/sum(wrec*wrec);cc=4.0_dp/real(n+4,dp)
    chi_n=sqrt(real(n,dp))*(1.0_dp-1.0_dp/(4.0_dp*n)+1.0_dp/(21.0_dp*n*n))
    cs=(mu_eff+2.0_dp)/(real(n,dp)+mu_eff+3.0_dp)
    damp=1.0_dp+2.0_dp*max(0.0_dp,sqrt((mu_eff-1.0_dp)/real(n+1,dp))-1.0_dp)+cs
    mucov=mu_eff;ccov=(1.0_dp/mucov)*(2.0_dp/(real(n,dp)+sqrt(2.0_dp))**2)+ &
      (1.0_dp-1.0_dp/mucov)*min(1.0_dp,(2.0_dp*mu_eff-1.0_dp)/(real(n+2,dp)**2+mu_eff))
    fc=0
    values=huge(1.0_dp)
    besttrace=huge(1.0_dp)
    sstop=.false.
    do gen=1,op%maxgen
      if(all(sd<=tiny(1.0_dp))) exit
      call rtnorm_matrix(popsize,mu,sd,lo,hi,pop);pop(:,1)=mu
      do i=1,popsize
        fit(i,:)=0.0_dp
        do rep=1,max(1,op%replicates)
          call eval_vector(fn,pop(:,i),fit_tmp)
          fit(i,:)=fit(i,:)+fit_tmp
          fc=fc+1
        end do
        fit(i,:)=fit(i,:)/real(max(1,op%replicates),dp)
        fg(i)=sum(fit(i,:)*wgt)
      end do
      call sort_indices(fg,idx);sups=idx(1:parents);selected=pop(:,sups);bestv=fg(idx(1))
      do j=1,nvar
        call sort_indices(fit(sups,j),lidx)
        do i=1,n
          optind(i,j)=sum(wrec*selected(i,lidx))
          optsd(i,j)=sqrt(max(sum((wrec*selected(i,lidx))**2)-optind(i,j)**2,0.0_dp))
        end do
      end do
      optind=(1.0_dp-op%alpha)*mu_state+op%alpha*optind
      optsd=sqrt(max((1.0_dp-op%alpha)*(mu_state**2+sigma_state**2)+ &
        op%alpha*(optind**2+optsd**2)-optind**2,0.0_dp))
      mu_state=optind
      sigma_state=optsd
      if(nvar==1) then
        wnew(:,1)=1.0_dp
      else
        xcv=sigma_state/spread(rangev,2,nvar)
        do j=1,nvar
          cvmins(j)=0.9_dp*(minval(xcv(:,j))+1.0e-20_dp)
          cvmaxs(j)=1.1_dp*(maxval(xcv(:,j))+1.0e-20_dp)
          do i=1,n
            wnew(i,j)=((cvmaxs(j)-xcv(i,j))/max(cvmaxs(j)-cvmins(j),tiny(1.0_dp)))**4
            if(wnew(i,j)<=0.0_dp) wnew(i,j)=1.0e-20_dp
          end do
          wnew(:,j)=wnew(:,j)/sum(wnew(:,j))
        end do
        do i=1,n
          wnew(i,:)=wnew(i,:)/sum(wnew(i,:))
        end do
      end if
      wmat=0.0_dp
      do j=1,nvar
        call sort_indices(fit(sups,j),lidx)
        do k=1,parents
          do i=1,n
            wmat(i,k)=wmat(i,k)+wnew(i,j)*wrec(lidx(k))
          end do
        end do
      end do
      mu_eff_row=1.0_dp/max(sum(wmat*wmat,dim=2),tiny(1.0_dp));mu_eff_global=max(sum(mu_eff_row)/real(n,dp),1.0_dp)
      cs=(mu_eff_global+2.0_dp)/(real(n,dp)+mu_eff_global+3.0_dp)
      damp=1.0_dp+2.0_dp*max(0.0_dp,sqrt((mu_eff_global-1.0_dp)/real(n+1,dp))-1.0_dp)+cs
      mucov=mu_eff_global;ccov=(1.0_dp/mucov)*(2.0_dp/(real(n,dp)+sqrt(2.0_dp))**2)+ &
        (1.0_dp-1.0_dp/mucov)*min(1.0_dp,(2.0_dp*mu_eff_global-1.0_dp)/(real(n+2,dp)**2+mu_eff_global))
      old_mu=mu;mu=sum(wmat*selected,dim=2)
      pc=(1.0_dp-cc)*pc+sqrt(cc*(2.0_dp-cc))*sqrt(mu_eff_row)*(mu-old_mu)/step
      where(sigma_diag>tiny(1.0_dp))
        ps=(1.0_dp-cs)*ps+sqrt(cs*(2.0_dp-cs))*sqrt(mu_eff_row)*((mu-old_mu)/sqrt(sigma_diag))/step
      elsewhere
        ps=0.0_dp
      end where
      sigma_diag=(1.0_dp-ccov)*sigma_diag+(ccov/mucov)*pc*pc+ &
        ccov*(1.0_dp-1.0_dp/mucov)*sum(wmat*((selected-spread(old_mu,2,parents))/step)**2,dim=2)
      step=step*exp(cs*(sqrt(sum(ps*ps))/chi_n-1.0_dp)/damp)
      sigma_diag=max(sigma_diag,0.0_dp);where(.not.use_active)sigma_diag=0.0_dp
      sd=step*sqrt(sigma_diag)
      values(gen)=fg(1)
      besttrace(gen)=bestv
      stopflag=.false.
      select case(op%termination)
      case(0)
        stopflag=.false.
      case(1)
        stopflag=step<op%convergence
      case(2)
        sstop(gen)=stop2(values(1:gen),op%reltol,op%fn_smoothing)
        stopflag = n_stop(sstop(1:gen),op%max_no_improvement)
      case(3)
        sstop(gen)=stop3(values(1:gen),op%reltol,op%fn_smoothing)
        stopflag = n_stop(sstop(1:gen),op%max_no_improvement)
      case(4)
        sstop(gen)=stop4(values(1:gen),op%reltol,op%fn_smoothing)
        stopflag = n_stop(sstop(1:gen),op%fn_smoothing)
      case default
        error stop "ahres: invalid termination criterion"
      end select
      if(stopflag) exit
    end do
    call eval_vector(fn,mu,fit_tmp);fc=fc+1
    allocate(result%par(n),result%partial(nvar));result%par=mu;result%partial=fit_tmp;result%value=sum(fit_tmp*wgt)
    result%function_count=fc;result%generations=min(gen,op%maxgen);result%convergence=merge(0,1,gen<op%maxgen)
    if (result%convergence == 0) then
      result%message="stopping criterion reached"
    else
      result%message="maximum generations reached"
    end if
    allocate(result%value_trace(result%generations),result%best_trace(result%generations))
    result%value_trace=values(1:result%generations);result%best_trace=besttrace(1:result%generations)
  end subroutine ahres

  function evaluate(fn,x,options,count) result(f)
    procedure(scalar_objective) :: fn
    real(dp),intent(in)::x(:)
    type(optim_options),intent(in)::options
    integer,intent(inout)::count
    real(dp)::f,tmp
    integer::r,nr
    nr=max(1,options%replicates);f=0.0_dp
    do r=1,nr
      tmp=eval_scalar(fn,x);f=f+tmp;count=count+1
    end do
    f=f/real(nr,dp)
    if(options%maximize)f=-f
    f=f/options%fnscale
  end function evaluate

  function unscale_value(f,options) result(v)
    real(dp),intent(in)::f
    type(optim_options),intent(in)::options
    real(dp)::v
    v=f*options%fnscale
    if(options%maximize)v=-v
  end function unscale_value

  subroutine get_gradient(fn,x,g,options,gr,fc,gc,active)
    procedure(scalar_objective)::fn
    real(dp),intent(in)::x(:)
    real(dp),intent(out)::g(:)
    type(optim_options),intent(in)::options
    procedure(gradient_callback),optional::gr
    integer,intent(inout)::fc,gc
    logical,intent(in)::active(:)
    if(present(gr))then
      call dispatch_gradient(gr,x,g);gc=gc+1
      if(options%maximize)g=-g
      g=g/options%fnscale
    else
      call numerical_gradient_local(fn,x,g,options,fc,active);gc=gc+1
    end if
    where(.not.active)g=0.0_dp
  end subroutine get_gradient

  subroutine numerical_gradient_local(fn,x,g,options,fc,active)
    procedure(scalar_objective)::fn
    real(dp),intent(in)::x(:)
    real(dp),intent(out)::g(:)
    type(optim_options),intent(in)::options
    integer,intent(inout)::fc
    logical,intent(in)::active(:)
    real(dp),allocatable::xp(:),xm(:),h(:),df(:,:)
    real(dp)::fx,fac
    integer::i,k,m,r
    allocate(xp(size(x)),xm(size(x)),h(size(x)));g=0.0_dp
    select case(trim(options%gradient_method))
    case("central")
      do i=1,size(x)
        if(.not.active(i)) cycle
        h(i)=1.0e-8_dp*max(abs(x(i)),1.0_dp)
        xp=x
        xm=x
        xp(i)=xp(i)+h(i)
        xm(i)=xm(i)-h(i)
        g(i)=(evaluate(fn,xp,options,fc)-evaluate(fn,xm,options,fc))/(2.0_dp*h(i))
      end do
    case("backward")
      fx=evaluate(fn,x,options,fc)
      do i=1,size(x)
        if(.not.active(i)) cycle
        h(i)=1.0e-8_dp*max(abs(x(i)),1.0_dp)
        xm=x
        xm(i)=xm(i)-h(i)
        g(i)=(fx-evaluate(fn,xm,options,fc))/h(i)
      end do
    case("richardson")
      r=4;allocate(df(r,size(x)));df=0.0_dp
      do i=1,size(x);h(i)=1.0e-4_dp*max(abs(x(i)),1.0_dp);end do
      do k=1,r
        do i=1,size(x)
          if(.not.active(i)) cycle
          xp=x
          xm=x
          xp(i)=xp(i)+h(i)
          xm(i)=xm(i)-h(i)
          df(k,i)=(evaluate(fn,xp,options,fc)-evaluate(fn,xm,options,fc))/(2.0_dp*h(i))
        end do
        h=h/2.0_dp
      end do
      do m=1,r-1;fac=4.0_dp**m;do k=1,r-m;df(k,:)=(fac*df(k+1,:)-df(k,:))/(fac-1.0_dp);end do;end do
      g=df(1,:)
    case default
      fx=evaluate(fn,x,options,fc)
      do i=1,size(x)
        if(.not.active(i)) cycle
        h(i)=1.0e-8_dp*max(abs(x(i)),1.0_dp)
        xp=x
        xp(i)=xp(i)+h(i)
        g(i)=(evaluate(fn,xp,options,fc)-fx)/h(i)
      end do
    end select
  end subroutine numerical_gradient_local


  function outer(a,b) result(c)
    real(dp),intent(in)::a(:),b(:)
    real(dp)::c(size(a),size(b))
    integer::i
    do i=1,size(a);c(i,:)=a(i)*b;end do
  end function outer

  subroutine simplex_order(f,ilo,ihi,inh)
    real(dp),intent(in)::f(:)
    integer,intent(out)::ilo,ihi,inh
    integer::i
    ilo=1;ihi=1
    do i=2,size(f);if(f(i)<f(ilo))ilo=i;if(f(i)>f(ihi))ihi=i;end do
    inh=ilo
    do i=1,size(f)
      if(i==ihi)cycle
      if(inh==ihi .or. f(i)>f(inh))inh=i
    end do
  end subroutine simplex_order

  function eval_scalar(fn,x) result(f)
    procedure(scalar_objective)::fn
    real(dp),intent(in)::x(:)
    real(dp)::f
    f=fn(x)
  end function eval_scalar

  subroutine eval_vector(fn,x,f)
    procedure(vector_objective)::fn
    real(dp),intent(in)::x(:)
    real(dp),intent(out)::f(:)
    call fn(x,f)
  end subroutine eval_vector

  subroutine dispatch_gradient(gr,x,g)
    procedure(gradient_callback)::gr
    real(dp),intent(in)::x(:)
    real(dp),intent(out)::g(:)
    call gr(x,g)
  end subroutine dispatch_gradient
end module calibrar_optimization
