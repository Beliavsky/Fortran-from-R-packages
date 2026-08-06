module arfima_fit
  use arfima_kinds, only : dp
  use arfima_status, only : arfima_ok, arfima_invalid_input, arfima_no_convergence
  use arfima_types, only : arfima_spec, arfima_parameters, arfima_fit_result, arfima_mode_set, arfima_error, &
    transfer_spec, set_error, long_memory_none, long_memory_fd, long_memory_fgn, long_memory_pla
  use arfima_polynomial, only : ar_to_pacf, pacf_to_ar, difference_series
  use arfima_autocov, only : tacvf_arfima
  use arfima_durbin, only : durbin_levinson
  use arfima_transfer, only : apply_static_regression, apply_transfer_function
  use arfima_optimizer, only : objective_function, nelder_mead, numerical_hessian
  use arfima_linalg, only : invert_matrix, solve_linear
  use arfima_random, only : random_normal_vector
  implicit none
  private
  public :: fit_arfima, fit_arfima_modes, arfima0_fit, mode_distance, weed_modes, remove_mode, best_modes

  type :: fit_context
    type(arfima_spec) :: spec
    real(dp), allocatable :: y(:)
    real(dp), allocatable :: xreg(:,:)
    type(transfer_spec) :: transfer
  end type fit_context

contains

  subroutine fit_arfima(spec,z,result,start,max_iterations,tolerance)
    type(arfima_spec),intent(in)::spec
    real(dp),intent(in)::z(:)
    type(arfima_fit_result),intent(out)::result
    type(arfima_parameters),intent(in),optional::start
    integer,intent(in),optional::max_iterations
    real(dp),intent(in),optional::tolerance
    type(fit_context)::context
    type(arfima_parameters)::initial,fitpar
    type(arfima_error)::err
    real(dp),allocatable::x(:),hess(:,:),hinv(:,:),jac(:,:),natp(:),natm(:),xp(:),xm(:)
    type(arfima_fit_result)::tmp
    type(arfima_parameters)::pp,pm
    real(dp)::fval,tol,hstep
    integer::iter,eval,info,maxit,npar,knat,i,j

    result%spec=spec
    call prepare_context(spec,z,context,err)
    if(err%code/=arfima_ok) then; result%error=err; return; end if
    call default_parameters(spec,context,initial)
    if(present(start)) initial=start
    call pack_transformed(spec,initial,x,err)
    if(err%code/=arfima_ok) then; result%error=err; return; end if
    maxit=2500; if(present(max_iterations)) maxit=max_iterations
    tol=1.0e-7_dp; if(present(tolerance)) tol=tolerance
    call nelder_mead(fit_objective,context,x,fval,iter,eval,info,maxit,tol,0.15_dp)
    call unpack_transformed(spec,x,fitpar,err)
    if(err%code/=arfima_ok) then; result%error=err; return; end if
    call evaluate_fit(context,fitpar,tmp)
    result=tmp
    result%spec=spec
    result%iterations=iter
    result%evaluations=eval
    result%converged=(info==arfima_ok .and. tmp%error%code==arfima_ok)
    if(.not.result%converged .and. result%error%code==arfima_ok) call set_error(result%error,arfima_no_convergence,'optimizer did not satisfy convergence tolerances')

    npar=size(x)
    if(npar>0 .and. tmp%error%code==arfima_ok) then
      call numerical_hessian(fit_objective,context,x,hess,info)
      if(info==arfima_ok) then
        call invert_matrix(hess,hinv,info)
        if(info==0) then
          knat=natural_count(spec)
          allocate(jac(knat,npar),xp(npar),xm(npar))
          jac=0.0_dp
          do j=1,npar
            hstep=epsilon(1.0_dp)**0.25_dp*max(1.0_dp,abs(x(j)))
            xp=x; xm=x; xp(j)=xp(j)+hstep; xm(j)=xm(j)-hstep
            call unpack_transformed(spec,xp,pp,err); natp=natural_vector(spec,pp)
            call unpack_transformed(spec,xm,pm,err); natm=natural_vector(spec,pm)
            jac(:,j)=(natp-natm)/(2.0_dp*hstep)
          end do
          result%covariance=matmul(jac,matmul(hinv,transpose(jac)))
          allocate(result%standard_error(knat))
          do i=1,knat; result%standard_error(i)=sqrt(max(0.0_dp,result%covariance(i,i))); end do
        end if
      end if
    end if
  end subroutine fit_arfima

  function fit_objective(x,context_any) result(value)
    real(dp),intent(in)::x(:)
    class(*),intent(inout)::context_any
    real(dp)::value
    type(arfima_parameters)::p
    type(arfima_error)::err
    type(arfima_fit_result)::res
    select type(context_any)
    type is(fit_context)
      call unpack_transformed(context_any%spec,x,p,err)
      if(err%code/=arfima_ok) then; value=1.0e100_dp; return; end if
      call evaluate_fit(context_any,p,res)
      if(res%error%code/=arfima_ok .or. .not.(res%loglik>-huge(1.0_dp))) then
        value=1.0e100_dp
      else
        value=-res%loglik
      end if
    class default
      value=1.0e100_dp
    end select
  end function fit_objective

  subroutine evaluate_fit(context,p,result)
    type(fit_context),intent(in)::context
    type(arfima_parameters),intent(in)::p
    type(arfima_fit_result),intent(out)::result
    real(dp),allocatable::series(:),r(:),effect(:)
    type(arfima_error)::err
    type(arfima_parameters)::pcov
    type(transfer_spec)::tlocal
    integer::n,k
    type_dl: block
      use arfima_types, only : dl_result
      type(dl_result)::dl
      result%spec=context%spec
      result%parameters=p
      n=size(context%y)
      if(context%spec%use_regression) then
        call apply_static_regression(context%y,context%xreg,p%beta,p%mean,series,err)
      else if(context%spec%use_transfer) then
        tlocal=context%transfer
        tlocal%delta=p%delta
        tlocal%omega=p%omega
        call apply_transfer_function(context%y,tlocal,p%mean,series,effect,err)
      else
        allocate(series(n)); series=context%y-p%mean; call set_error(err,arfima_ok,'')
      end if
      if(err%code/=arfima_ok) then; result%error=err; exit type_dl; end if
      pcov=p; pcov%mean=0.0_dp
      call tacvf_arfima(context%spec,pcov,n-1,1.0_dp,r,err)
      if(err%code/=arfima_ok) then; result%error=err; exit type_dl; end if
      call durbin_levinson(r,series,dl)
      if(dl%error%code/=arfima_ok) then; result%error=dl%error; exit type_dl; end if
      result%loglik=dl%loglik
      result%sigma2=dl%sigma2_mle
      result%residuals=dl%residuals
      allocate(result%fitted(n)); result%fitted=context%y-dl%residuals
      k=natural_count(context%spec)+1
      result%aic=-2.0_dp*result%loglik+2.0_dp*real(k,dp)
      result%bic=-2.0_dp*result%loglik+log(real(n,dp))*real(k,dp)
      call set_error(result%error,arfima_ok,'')
    end block type_dl
  end subroutine evaluate_fit

  subroutine prepare_context(spec,z,context,error)
    type(arfima_spec),intent(in)::spec
    real(dp),intent(in)::z(:)
    type(fit_context),intent(out)::context
    type(arfima_error),intent(out)::error
    real(dp),allocatable::col(:)
    type(arfima_error)::err
    integer::j,nout
    call set_error(error,arfima_ok,'')
    context%spec=spec
    if(spec%use_regression .and. spec%use_transfer) then
      call set_error(error,arfima_invalid_input,'static regression and dynamic transfer modes are mutually exclusive'); return
    end if
    call difference_series(z,spec%dint,spec%dseas,spec%period,context%y,err)
    if(spec%dint+spec%dseas*spec%period==0) then
      if(allocated(context%y)) deallocate(context%y)
      allocate(context%y(size(z))); context%y=z; call set_error(err,arfima_ok,'')
    end if
    if(err%code/=arfima_ok) then; error=err; return; end if
    nout=size(context%y)
    if(spec%use_regression) then
      if(.not.allocated(spec%xreg) .or. size(spec%xreg,1)/=size(z)) then
        call set_error(error,arfima_invalid_input,'xreg must have the same number of rows as z'); return
      end if
      allocate(context%xreg(nout,size(spec%xreg,2)))
      do j=1,size(spec%xreg,2)
        if(spec%dint+spec%dseas*spec%period>0) then
          call difference_series(spec%xreg(:,j),spec%dint,spec%dseas,spec%period,col,err)
          if(err%code/=arfima_ok) then; error=err; return; end if
          context%xreg(:,j)=col
        else
          context%xreg(:,j)=spec%xreg(:,j)
        end if
      end do
      context%spec%xreg=context%xreg
    end if
    if(spec%use_transfer) then
      context%transfer=spec%transfer
      if(allocated(context%transfer%x)) deallocate(context%transfer%x)
      if(.not.allocated(spec%transfer%x) .or. size(spec%transfer%x,1)/=size(z)) then
        call set_error(error,arfima_invalid_input,'transfer x must have the same number of rows as z'); return
      end if
      allocate(context%transfer%x(nout,size(spec%transfer%x,2)))
      do j=1,size(spec%transfer%x,2)
        if(spec%dint+spec%dseas*spec%period>0) then
          call difference_series(spec%transfer%x(:,j),spec%dint,spec%dseas,spec%period,col,err)
          if(err%code/=arfima_ok) then; error=err; return; end if
          context%transfer%x(:,j)=col
        else
          context%transfer%x(:,j)=spec%transfer%x(:,j)
        end if
      end do
      context%spec%transfer=context%transfer
    end if
  end subroutine prepare_context

  subroutine default_parameters(spec,context,p)
    type(arfima_spec),intent(in)::spec
    type(fit_context),intent(in)::context
    type(arfima_parameters),intent(out)::p
    real(dp),allocatable::xtx(:,:),xty(:),sol(:)
    integer::info
    allocate(p%phi(spec%p),p%theta(spec%q),p%phiseas(spec%pseas),p%thetaseas(spec%qseas))
    p%phi=0.0_dp; p%theta=0.0_dp; p%phiseas=0.0_dp; p%thetaseas=0.0_dp
    p%dfrac=0.1_dp; p%dfs=0.1_dp; p%hurst=0.7_dp; p%hurst_seasonal=0.7_dp; p%alpha=0.7_dp; p%alpha_seasonal=0.7_dp
    p%mean=sum(context%y)/real(size(context%y),dp)
    if(spec%use_regression) then
      allocate(p%beta(size(context%xreg,2))); p%beta=0.0_dp
      xtx=matmul(transpose(context%xreg),context%xreg); xty=matmul(transpose(context%xreg),context%y-p%mean)
      call solve_linear(xtx,xty,sol,info); if(info==0) p%beta=sol
    else
      allocate(p%beta(0))
    end if
    if(spec%use_transfer) then
      allocate(p%delta(sum(context%transfer%r)),p%omega(sum(context%transfer%s)))
      p%delta=0.0_dp; p%omega=0.0_dp
    else
      allocate(p%delta(0),p%omega(0))
    end if
  end subroutine default_parameters

  subroutine pack_transformed(spec,p,x,error)
    type(arfima_spec),intent(in)::spec
    type(arfima_parameters),intent(in)::p
    real(dp),allocatable,intent(out)::x(:)
    type(arfima_error),intent(out)::error
    real(dp),allocatable::pacf(:)
    integer::i,pos,k,rr
    call set_error(error,arfima_ok,'')
    allocate(x(transformed_count(spec))); x=0.0_dp; pos=0
    call pack_group(p%phi,spec%p); call pack_group(p%theta,spec%q); call pack_group(p%phiseas,spec%pseas); call pack_group(p%thetaseas,spec%qseas)
    if(spec%lmodel/=long_memory_none) then; pos=pos+1; x(pos)=long_inverse(spec%lmodel,long_value(spec%lmodel,p,.false.)); end if
    if(spec%slmodel/=long_memory_none) then; pos=pos+1; x(pos)=long_inverse(spec%slmodel,long_value(spec%slmodel,p,.true.)); end if
    if(spec%use_regression) then
      if(.not.allocated(p%beta) .or. size(p%beta)/=size(spec%xreg,2)) then; call set_error(error,arfima_invalid_input,'start beta has wrong length'); return; end if
      x(pos+1:pos+size(p%beta))=p%beta; pos=pos+size(p%beta)
    end if
    if(spec%use_transfer) then
      rr=0
      do k=1,size(spec%transfer%r)
        if(spec%transfer%r(k)>0) then
          pacf=ar_to_pacf(p%delta(rr+1:rr+spec%transfer%r(k)))
          do i=1,size(pacf); pos=pos+1; x(pos)=atanh(max(-0.999999_dp,min(0.999999_dp,pacf(i)))); end do
          rr=rr+spec%transfer%r(k)
        end if
      end do
      x(pos+1:pos+size(p%omega))=p%omega; pos=pos+size(p%omega)
    end if
    if(spec%estimate_mean) then; pos=pos+1; x(pos)=p%mean; end if
  contains
    subroutine pack_group(group,n)
      real(dp),allocatable,intent(in)::group(:)
      integer,intent(in)::n
      integer::j
      if(n==0) return
      if(.not.allocated(group) .or. size(group)/=n) then; call set_error(error,arfima_invalid_input,'start AR/MA group has wrong length'); return; end if
      pacf=ar_to_pacf(group)
      do j=1,n; pos=pos+1; x(pos)=atanh(max(-0.999999_dp,min(0.999999_dp,pacf(j)))); end do
    end subroutine pack_group
  end subroutine pack_transformed

  subroutine unpack_transformed(spec,x,p,error)
    type(arfima_spec),intent(in)::spec
    real(dp),intent(in)::x(:)
    type(arfima_parameters),intent(out)::p
    type(arfima_error),intent(out)::error
    real(dp),allocatable::pacf(:),group(:)
    integer::pos,k,nr,rr
    call set_error(error,arfima_ok,'')
    if(size(x)/=transformed_count(spec)) then; call set_error(error,arfima_invalid_input,'transformed parameter vector has wrong length'); return; end if
    pos=0
    call unpack_group(spec%p,p%phi); call unpack_group(spec%q,p%theta); call unpack_group(spec%pseas,p%phiseas); call unpack_group(spec%qseas,p%thetaseas)
    p%dfrac=0.0_dp; p%dfs=0.0_dp; p%hurst=0.5_dp; p%hurst_seasonal=0.5_dp; p%alpha=1.0_dp; p%alpha_seasonal=1.0_dp
    if(spec%lmodel/=long_memory_none) then; pos=pos+1; call set_long(spec%lmodel,p,.false.,long_forward(spec%lmodel,x(pos))); end if
    if(spec%slmodel/=long_memory_none) then; pos=pos+1; call set_long(spec%slmodel,p,.true.,long_forward(spec%slmodel,x(pos))); end if
    if(spec%use_regression) then
      allocate(p%beta(size(spec%xreg,2))); p%beta=x(pos+1:pos+size(p%beta)); pos=pos+size(p%beta)
    else; allocate(p%beta(0)); end if
    if(spec%use_transfer) then
      allocate(p%delta(sum(spec%transfer%r)),p%omega(sum(spec%transfer%s))); rr=0
      do k=1,size(spec%transfer%r)
        nr=spec%transfer%r(k)
        if(nr>0) then
          allocate(pacf(nr)); pacf=tanh(x(pos+1:pos+nr)); group=pacf_to_ar(pacf); p%delta(rr+1:rr+nr)=group
          pos=pos+nr; rr=rr+nr; deallocate(pacf)
        end if
      end do
      p%omega=x(pos+1:pos+size(p%omega)); pos=pos+size(p%omega)
    else; allocate(p%delta(0),p%omega(0)); end if
    if(spec%estimate_mean) then; pos=pos+1; p%mean=x(pos); else; p%mean=0.0_dp; end if
  contains
    subroutine unpack_group(n,out)
      integer,intent(in)::n
      real(dp),allocatable,intent(out)::out(:)
      if(n==0) then; allocate(out(0)); return; end if
      allocate(pacf(n)); pacf=tanh(x(pos+1:pos+n)); out=pacf_to_ar(pacf); pos=pos+n; deallocate(pacf)
    end subroutine unpack_group
  end subroutine unpack_transformed

  integer function transformed_count(spec) result(n)
    type(arfima_spec),intent(in)::spec
    n=spec%p+spec%q+spec%pseas+spec%qseas
    if(spec%lmodel/=long_memory_none) n=n+1
    if(spec%slmodel/=long_memory_none) n=n+1
    if(spec%use_regression) n=n+size(spec%xreg,2)
    if(spec%use_transfer) n=n+sum(spec%transfer%r)+sum(spec%transfer%s)
    if(spec%estimate_mean) n=n+1
  end function transformed_count

  integer function natural_count(spec) result(n)
    type(arfima_spec),intent(in)::spec
    n=transformed_count(spec)
  end function natural_count

  function natural_vector(spec,p) result(v)
    type(arfima_spec),intent(in)::spec
    type(arfima_parameters),intent(in)::p
    real(dp),allocatable::v(:)
    integer::pos
    allocate(v(natural_count(spec))); pos=0
    call add(p%phi); call add(p%theta); call add(p%phiseas); call add(p%thetaseas)
    if(spec%lmodel/=long_memory_none) then; pos=pos+1; v(pos)=long_value(spec%lmodel,p,.false.); end if
    if(spec%slmodel/=long_memory_none) then; pos=pos+1; v(pos)=long_value(spec%slmodel,p,.true.); end if
    if(spec%use_regression) call add(p%beta)
    if(spec%use_transfer) then; call add(p%delta); call add(p%omega); end if
    if(spec%estimate_mean) then; pos=pos+1; v(pos)=p%mean; end if
  contains
    subroutine add(a)
      real(dp),allocatable,intent(in)::a(:)
      if(allocated(a)) then; v(pos+1:pos+size(a))=a; pos=pos+size(a); end if
    end subroutine add
  end function natural_vector

  real(dp) function long_forward(model,z) result(v)
    integer,intent(in)::model
    real(dp),intent(in)::z
    real(dp)::s
    s=1.0_dp/(1.0_dp+exp(-max(-40.0_dp,min(40.0_dp,z))))
    select case(model)
    case(long_memory_fd); v=-0.99_dp+1.48_dp*s
    case(long_memory_fgn); v=0.001_dp+0.998_dp*s
    case(long_memory_pla); v=0.001_dp+2.998_dp*s
    case default; v=0.0_dp
    end select
  end function long_forward

  real(dp) function long_inverse(model,v) result(z)
    integer,intent(in)::model
    real(dp),intent(in)::v
    real(dp)::s
    select case(model)
    case(long_memory_fd); s=(v+0.99_dp)/1.48_dp
    case(long_memory_fgn); s=(v-0.001_dp)/0.998_dp
    case(long_memory_pla); s=(v-0.001_dp)/2.998_dp
    case default; s=0.5_dp
    end select
    s=max(1.0e-8_dp,min(1.0_dp-1.0e-8_dp,s)); z=log(s/(1.0_dp-s))
  end function long_inverse

  real(dp) function long_value(model,p,seasonal) result(v)
    integer,intent(in)::model
    type(arfima_parameters),intent(in)::p
    logical,intent(in)::seasonal
    v=0.0_dp
    if(.not.seasonal) then
      select case(model); case(long_memory_fd); v=p%dfrac; case(long_memory_fgn); v=p%hurst; case(long_memory_pla); v=p%alpha; end select
    else
      select case(model); case(long_memory_fd); v=p%dfs; case(long_memory_fgn); v=p%hurst_seasonal; case(long_memory_pla); v=p%alpha_seasonal; end select
    end if
  end function long_value

  subroutine set_long(model,p,seasonal,v)
    integer,intent(in)::model
    type(arfima_parameters),intent(inout)::p
    logical,intent(in)::seasonal
    real(dp),intent(in)::v
    if(.not.seasonal) then
      select case(model); case(long_memory_fd); p%dfrac=v; case(long_memory_fgn); p%hurst=v; case(long_memory_pla); p%alpha=v; end select
    else
      select case(model); case(long_memory_fd); p%dfs=v; case(long_memory_fgn); p%hurst_seasonal=v; case(long_memory_pla); p%alpha_seasonal=v; end select
    end if
  end subroutine set_long

  subroutine fit_arfima_modes(spec,z,n_starts,modes,seed,max_iterations)
    type(arfima_spec),intent(in)::spec
    real(dp),intent(in)::z(:)
    integer,intent(in)::n_starts
    type(arfima_mode_set),intent(out)::modes
    integer,intent(in),optional::seed,max_iterations
    type(arfima_parameters)::start
    type(fit_context)::context
    type(arfima_error)::err
    real(dp),allocatable::noise(:),x(:)
    integer::i,maxit
    if(present(seed)) then
      use_seed:block
        use arfima_random,only:set_random_seed
        call set_random_seed(seed)
      end block use_seed
    end if
    maxit=2000; if(present(max_iterations)) maxit=max_iterations
    allocate(modes%modes(max(1,n_starts)))
    call prepare_context(spec,z,context,err)
    if(err%code/=arfima_ok) then; modes%modes(1)%error=err; return; end if
    call default_parameters(spec,context,start)
    call fit_arfima(spec,z,modes%modes(1),start,maxit)
    do i=2,size(modes%modes)
      call pack_transformed(spec,start,x,err); allocate(noise(size(x))); call random_normal_vector(noise)
      x=x+0.8_dp*noise
      call unpack_transformed(spec,x,start,err)
      call fit_arfima(spec,z,modes%modes(i),start,maxit)
      deallocate(noise,x)
    end do
    call sort_modes(modes)
  end subroutine fit_arfima_modes

  real(dp) function mode_distance(a,b,pnorm) result(d)
    type(arfima_fit_result),intent(in)::a,b
    real(dp),intent(in),optional::pnorm
    real(dp)::p
    real(dp),allocatable::x(:),y(:)
    p=2.0_dp; if(present(pnorm)) p=pnorm
    x=natural_vector(a%spec,a%parameters); y=natural_vector(b%spec,b%parameters)
    if(size(x)/=size(y)) then; d=huge(1.0_dp); else; d=sum(abs(x-y)**p)**(1.0_dp/p); end if
  end function mode_distance

  subroutine weed_modes(modes,epsilon,pnorm)
    type(arfima_mode_set),intent(inout)::modes
    real(dp),intent(in),optional::epsilon,pnorm
    logical,allocatable::keep(:)
    type(arfima_fit_result),allocatable::tmp(:)
    real(dp)::eps,p
    integer::i,j,n
    eps=0.01_dp; if(present(epsilon)) eps=epsilon
    p=2.0_dp; if(present(pnorm)) p=pnorm
    n=size(modes%modes); allocate(keep(n)); keep=.true.
    do i=1,n
      if(.not.keep(i)) cycle
      do j=i+1,n
        if(keep(j)) then
          if(mode_distance(modes%modes(i),modes%modes(j),p)<eps) keep(j)=.false.
        end if
      end do
    end do
    allocate(tmp(count(keep))); j=0
    do i=1,n; if(keep(i)) then; j=j+1; tmp(j)=modes%modes(i); end if; end do
    call move_alloc(tmp,modes%modes)
  end subroutine weed_modes

  subroutine arfima0_fit(z,p,dint,q,lmodel,result)
    real(dp),intent(in)::z(:)
    integer,intent(in)::p,dint,q,lmodel
    type(arfima_fit_result),intent(out)::result
    type(arfima_spec)::spec
    spec%p=p; spec%q=q; spec%dint=dint; spec%lmodel=lmodel; spec%estimate_mean=.true.
    call fit_arfima(spec,z,result)
  end subroutine arfima0_fit

  subroutine remove_mode(modes,index)
    type(arfima_mode_set),intent(inout)::modes
    integer,intent(in)::index
    type(arfima_fit_result),allocatable::tmp(:)
    integer::n
    n=size(modes%modes)
    if(index<1 .or. index>n) return
    allocate(tmp(n-1))
    if(index>1) tmp(1:index-1)=modes%modes(1:index-1)
    if(index<n) tmp(index:n-1)=modes%modes(index+1:n)
    call move_alloc(tmp,modes%modes)
  end subroutine remove_mode

  subroutine best_modes(modes,bestn)
    type(arfima_mode_set),intent(inout)::modes
    integer,intent(in)::bestn
    type(arfima_fit_result),allocatable::tmp(:)
    integer::n
    call sort_modes(modes)
    n=max(0,min(bestn,size(modes%modes)))
    allocate(tmp(n))
    if(n>0) tmp=modes%modes(1:n)
    call move_alloc(tmp,modes%modes)
  end subroutine best_modes

  subroutine sort_modes(modes)
    type(arfima_mode_set),intent(inout)::modes
    type(arfima_fit_result)::temp
    integer::i,j
    do i=1,size(modes%modes)-1
      do j=i+1,size(modes%modes)
        if(modes%modes(j)%loglik>modes%modes(i)%loglik) then; temp=modes%modes(i); modes%modes(i)=modes%modes(j); modes%modes(j)=temp; end if
      end do
    end do
  end subroutine sort_modes

end module arfima_fit
