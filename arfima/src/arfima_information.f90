module arfima_information_mod
  use arfima_kinds, only : dp, pi_dp
  use arfima_status, only : arfima_ok, arfima_invalid_input, arfima_not_stationary
  use arfima_types, only : arfima_spec, arfima_parameters, arfima_error, set_error, &
    long_memory_none, long_memory_fd, long_memory_fgn, long_memory_pla
  use arfima_polynomial, only : is_stationary_polynomial
  use arfima_autocov, only : tacvf_arfima
  use arfima_linalg, only : is_positive_definite
  implicit none
  private
  public :: arfima_information, identifiable_invertible
contains

  subroutine arfima_information(spec,params,information,error,nfreq,maxlag)
    type(arfima_spec),intent(in)::spec
    type(arfima_parameters),intent(in)::params
    real(dp),allocatable,intent(out)::information(:,:)
    type(arfima_error),intent(out)::error
    integer,intent(in),optional::nfreq,maxlag
    integer::m,l,np,i,j,k
    real(dp)::omega,eps,fp,fm
    real(dp),allocatable::scores(:,:),rp(:),rm(:)
    type(arfima_parameters)::pp,pm
    type(arfima_error)::ep,em

    call set_error(error,arfima_ok,'')
    m=512; if(present(nfreq)) m=nfreq
    l=512; if(present(maxlag)) l=maxlag
    np=parameter_count(spec)
    allocate(information(np,np))
    information=0.0_dp
    if(np==0) return
    allocate(scores(np,m))
    scores=0.0_dp
    do i=1,np
      pp=params; pm=params
      eps=1.0e-5_dp*max(1.0_dp,abs(parameter_value(spec,params,i)))
      call set_parameter(spec,pp,i,parameter_value(spec,params,i)+eps)
      call set_parameter(spec,pm,i,parameter_value(spec,params,i)-eps)
      call tacvf_arfima(spec,pp,l,1.0_dp,rp,ep)
      call tacvf_arfima(spec,pm,l,1.0_dp,rm,em)
      if(ep%code/=arfima_ok .or. em%code/=arfima_ok) then
        call set_error(error,arfima_invalid_input,'unable to perturb parameters for information matrix')
        return
      end if
      do k=1,m
        omega=pi_dp*(real(k,dp)-0.5_dp)/real(m,dp)
        fp=spectral_from_acvf(rp,omega)
        fm=spectral_from_acvf(rm,omega)
        scores(i,k)=(log(max(fp,tiny(1.0_dp)))-log(max(fm,tiny(1.0_dp))))/(2.0_dp*eps)
      end do
    end do
    do i=1,np
      do j=i,np
        information(i,j)=0.5_dp*dot_product(scores(i,:),scores(j,:))/real(m,dp)
        information(j,i)=information(i,j)
      end do
    end do
  end subroutine arfima_information

  logical function identifiable_invertible(spec,params,check_information,error) result(ok)
    type(arfima_spec),intent(in)::spec
    type(arfima_parameters),intent(in)::params
    logical,intent(in),optional::check_information
    type(arfima_error),intent(out),optional::error
    logical::check
    real(dp),allocatable::info(:,:)
    type(arfima_error)::err
    ok=.false.; check=.true.; if(present(check_information)) check=check_information
    if(.not.array_stable(params%phi) .or. .not.array_stable(params%theta) .or. .not.array_stable(params%phiseas) .or. &
       .not.array_stable(params%thetaseas)) then
      call set_error(err,arfima_not_stationary,'AR or MA coefficients are outside the PACF stability region')
      if(present(error)) error=err
      return
    end if
    if(.not.valid_long(spec%lmodel,params,.false.) .or. .not.valid_long(spec%slmodel,params,.true.)) then
      call set_error(err,arfima_not_stationary,'long-memory parameter is outside its admissible range')
      if(present(error)) error=err
      return
    end if
    if(spec%period==0 .and. (array_size(params%phiseas)>0 .or. array_size(params%thetaseas)>0 .or. &
       spec%slmodel/=long_memory_none)) then
      call set_error(err,arfima_invalid_input,'seasonal terms require period >= 2')
      if(present(error)) error=err
      return
    end if
    if(check .and. parameter_count(spec)>0) then
      call arfima_information(spec,params,info,err,nfreq=256,maxlag=256)
      if(err%code/=arfima_ok) then
        call set_error(err,arfima_invalid_input,'unable to compute numerical information matrix')
        if(present(error)) error=err
        return
      end if
      if(.not.is_positive_definite(info)) then
        call set_error(err,arfima_invalid_input,'numerical information matrix is not positive definite')
        if(present(error)) error=err
        return
      end if
    end if
    call set_error(err,arfima_ok,''); if(present(error)) error=err; ok=.true.
  contains
    logical function array_stable(x)
      real(dp),allocatable,intent(in)::x(:)
      if(allocated(x)) then; array_stable=is_stationary_polynomial(x)
      else; array_stable=.true.; end if
    end function array_stable
    integer function array_size(x)
      real(dp),allocatable,intent(in)::x(:)
      if(allocated(x)) then; array_size=size(x); else; array_size=0; end if
    end function array_size
  end function identifiable_invertible

  integer function parameter_count(spec) result(n)
    type(arfima_spec),intent(in)::spec
    n=spec%p+spec%q+spec%pseas+spec%qseas
    if(spec%lmodel/=long_memory_none) n=n+1
    if(spec%slmodel/=long_memory_none) n=n+1
  end function parameter_count

  real(dp) function parameter_value(spec,p,index) result(v)
    type(arfima_spec),intent(in)::spec
    type(arfima_parameters),intent(in)::p
    integer,intent(in)::index
    integer::i
    i=index
    if(i<=spec%p) then; v=p%phi(i); return; else; i=i-spec%p; end if
    if(i<=spec%q) then; v=p%theta(i); return; else; i=i-spec%q; end if
    if(i<=spec%pseas) then; v=p%phiseas(i); return; else; i=i-spec%pseas; end if
    if(i<=spec%qseas) then; v=p%thetaseas(i); return; else; i=i-spec%qseas; end if
    if(spec%lmodel/=long_memory_none) then
      if(i==1) then; v=long_value(spec%lmodel,p,.false.); return; else; i=i-1; end if
    end if
    v=long_value(spec%slmodel,p,.true.)
  end function parameter_value

  subroutine set_parameter(spec,p,index,v)
    type(arfima_spec),intent(in)::spec
    type(arfima_parameters),intent(inout)::p
    integer,intent(in)::index
    real(dp),intent(in)::v
    integer::i
    i=index
    if(i<=spec%p) then; p%phi(i)=v; return; else; i=i-spec%p; end if
    if(i<=spec%q) then; p%theta(i)=v; return; else; i=i-spec%q; end if
    if(i<=spec%pseas) then; p%phiseas(i)=v; return; else; i=i-spec%pseas; end if
    if(i<=spec%qseas) then; p%thetaseas(i)=v; return; else; i=i-spec%qseas; end if
    if(spec%lmodel/=long_memory_none) then
      if(i==1) then; call set_long_value(spec%lmodel,p,.false.,v); return; else; i=i-1; end if
    end if
    call set_long_value(spec%slmodel,p,.true.,v)
  end subroutine set_parameter

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

  subroutine set_long_value(model,p,seasonal,v)
    integer,intent(in)::model
    type(arfima_parameters),intent(inout)::p
    logical,intent(in)::seasonal
    real(dp),intent(in)::v
    if(.not.seasonal) then
      select case(model); case(long_memory_fd); p%dfrac=v; case(long_memory_fgn); p%hurst=v; case(long_memory_pla); p%alpha=v; end select
    else
      select case(model); case(long_memory_fd); p%dfs=v; case(long_memory_fgn); p%hurst_seasonal=v; case(long_memory_pla); p%alpha_seasonal=v; end select
    end if
  end subroutine set_long_value

  logical function valid_long(model,p,seasonal)
    integer,intent(in)::model
    type(arfima_parameters),intent(in)::p
    logical,intent(in)::seasonal
    real(dp)::v
    v=long_value(model,p,seasonal)
    select case(model)
    case(long_memory_none); valid_long=.true.
    case(long_memory_fd); valid_long=v>-1.0_dp .and. v<0.5_dp
    case(long_memory_fgn); valid_long=v>0.0_dp .and. v<1.0_dp
    case(long_memory_pla); valid_long=v>0.0_dp .and. v<3.0_dp
    case default; valid_long=.false.
    end select
  end function valid_long

  real(dp) function spectral_from_acvf(r,omega) result(f)
    real(dp),intent(in)::r(:),omega
    integer::k
    f=r(1)
    do k=1,size(r)-1; f=f+2.0_dp*r(k+1)*cos(real(k,dp)*omega); end do
    f=max(f,tiny(1.0_dp))
  end function spectral_from_acvf

end module arfima_information_mod
