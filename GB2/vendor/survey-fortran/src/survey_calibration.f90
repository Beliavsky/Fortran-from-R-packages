! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module survey_calibration
  use survey_kinds, only : dp
  use survey_linalg, only : sym_pinv
  implicit none
  private
  integer, parameter, public :: CAL_LINEAR=1, CAL_RAKING=2, CAL_LOGIT=3, CAL_SINH=4
  public :: grake, calibrate_weights, poststratify_weights, rake_margins, trim_weights
contains

  subroutine grake(mm,ww,population,g,calfun,lower,upper,epsilon,maxit,sigma2,eta,converged,iterations,achieved)
    real(dp), intent(in) :: mm(:,:),ww(:),population(:)
    real(dp), intent(out) :: g(:)
    integer, intent(in), optional :: calfun,maxit
    real(dp), intent(in), optional :: lower,upper,epsilon,sigma2(:)
    real(dp), intent(inout), optional :: eta(:)
    logical, intent(out), optional :: converged
    integer, intent(out), optional :: iterations
    real(dp), intent(out), optional :: achieved
    real(dp), allocatable :: e(:),sig(:),xeta(:),deriv(:),fm1(:),sample_total(:),misfit(:),deta(:),pinv(:,:),T(:,:)
    real(dp) :: lo,hi,eps,scale,maxerr
    integer :: n,p,cf,mi,iter,i,j,k,rank,info,half
    logical :: conv
    n=size(mm,1); p=size(mm,2)
    if(size(ww)/=n .or. size(population)/=p .or. size(g)/=n) error stop 'grake: shape mismatch'
    cf=CAL_LINEAR; if(present(calfun)) cf=calfun
    lo=-huge(1.0_dp); hi=huge(1.0_dp); if(present(lower)) lo=lower; if(present(upper)) hi=upper
    eps=1e-7_dp; if(present(epsilon)) eps=epsilon
    mi=50; if(present(maxit)) mi=maxit
    if((cf==CAL_LOGIT .or. cf==CAL_SINH) .and. (abs(lo)>=huge(1.0_dp)/4 .or. abs(hi)>=huge(1.0_dp)/4)) &
      error stop 'grake: logit/sinh calibration requires finite bounds'
    allocate(e(p),sig(n),xeta(n),deriv(n),fm1(n),sample_total(p),misfit(p),deta(p),pinv(p,p),T(p,p))
    e=0; if(present(eta)) then; if(size(eta)/=p) error stop 'grake: eta size'; e=eta; end if
    sig=1; if(present(sigma2)) then; if(size(sigma2)/=n) error stop 'grake: sigma2 size'; sig=sigma2; end if
    if(any(sig<=0)) error stop 'grake: sigma2 must be positive'
    sample_total=matmul(transpose(mm),ww)
    scale=1
    if(all(abs(sample_total)>tiny(1.0_dp))) then
      if(minval(abs(population/sample_total))>20.0_dp) then
        scale=sum(population/sample_total)/real(p,dp); sample_total=sample_total*scale
      end if
    end if
    do iter=1,mi
      xeta=matmul(mm,e)/sig
      call calibration_functions(xeta,cf,lo,hi,fm1,deriv)
      g=1+fm1
      misfit=population-sample_total-matmul(transpose(mm),ww*scale*fm1)
      maxerr=maxval(abs(misfit)/(1+abs(population)))
      if(maxerr<eps) exit
      T=0
      do i=1,n
        do j=1,p
          do k=1,p
            T(j,k)=T(j,k)+mm(i,j)*mm(i,k)*ww(i)*scale*deriv(i)/sig(i)
          end do
        end do
      end do
      call sym_pinv(T,pinv,rank,info=info)
      deta=matmul(pinv,misfit); e=e+deta
      ! upstream step-halving when calibration map overflows
      do half=1,20
        xeta=matmul(mm,e)/sig; call calibration_functions(xeta,cf,lo,hi,fm1,deriv)
        if(all(ieee_finite_vec(1+fm1)) .and. all(ieee_finite_vec(deriv))) exit
        deta=deta/2; e=e-deta
      end do
    end do
    xeta=matmul(mm,e)/sig; call calibration_functions(xeta,cf,lo,hi,fm1,deriv); g=(1+fm1)*scale
    misfit=population-matmul(transpose(mm),ww*g); maxerr=maxval(abs(misfit)/(1+abs(population)))
    conv=maxerr<eps
    if(present(eta)) eta=e
    if(present(converged)) converged=conv
    if(present(iterations)) iterations=min(iter,mi)
    if(present(achieved)) achieved=maxerr
  end subroutine grake

  subroutine calibrate_weights(mm,weight,population,new_weight,calfun,lower,upper,epsilon,maxit,converged)
    real(dp), intent(in) :: mm(:,:),weight(:),population(:)
    real(dp), intent(out) :: new_weight(:)
    integer, intent(in), optional :: calfun,maxit
    real(dp), intent(in), optional :: lower,upper,epsilon
    logical, intent(out), optional :: converged
    real(dp), allocatable :: g(:)
    logical :: ok
    allocate(g(size(weight)))
    call grake(mm,weight,population,g,calfun,lower,upper,epsilon,maxit,converged=ok)
    new_weight=weight*g; if(present(converged)) converged=ok
  end subroutine calibrate_weights

  subroutine poststratify_weights(category,nlevels,population,weight,new_weight)
    integer, intent(in) :: category(:),nlevels
    real(dp), intent(in) :: population(:),weight(:)
    real(dp), intent(out) :: new_weight(:)
    real(dp) :: sw,fac
    integer :: l
    if(size(category)/=size(weight) .or. size(population)/=nlevels .or. size(new_weight)/=size(weight)) &
      error stop 'poststratify_weights: shape mismatch'
    new_weight=weight
    do l=1,nlevels
      sw=sum(weight,mask=category==l)
      if(sw>0) then
        fac=population(l)/sw; where(category==l) new_weight=weight*fac
      else if(abs(population(l))>tiny(1.0_dp)) then
        error stop 'poststratify_weights: population cell absent from sample'
      end if
    end do
  end subroutine poststratify_weights

  subroutine rake_margins(categories,nlevels,targets,weight,new_weight,maxit,epsilon,converged,iterations)
    integer, intent(in) :: categories(:,:),nlevels(:)
    real(dp), intent(in) :: targets(:,:),weight(:)
    real(dp), intent(out) :: new_weight(:)
    integer, intent(in), optional :: maxit
    real(dp), intent(in), optional :: epsilon
    logical, intent(out), optional :: converged
    integer, intent(out), optional :: iterations
    real(dp), allocatable :: old(:),tmp(:)
    real(dp) :: eps,err
    integer :: mi,it,m
    if(size(categories,1)/=size(weight) .or. size(categories,2)/=size(nlevels) .or. size(targets,2)/=size(nlevels)) &
      error stop 'rake_margins: shape mismatch'
    mi=50; if(present(maxit)) mi=maxit; eps=1e-7_dp; if(present(epsilon)) eps=epsilon
    allocate(old(size(weight)),tmp(size(weight))); new_weight=weight; err=huge(1.0_dp)
    do it=1,mi
      old=new_weight
      do m=1,size(nlevels)
        call poststratify_weights(categories(:,m),nlevels(m),targets(1:nlevels(m),m),new_weight,tmp)
        new_weight=tmp
      end do
      err=maxval(abs(new_weight-old)/(1+abs(old)))
      if(err<eps) exit
    end do
    if(present(converged)) converged=(err<eps)
    if(present(iterations)) iterations=min(it,mi)
  end subroutine rake_margins

  subroutine trim_weights(weight,lower,upper,new_weight,strict,converged)
    real(dp), intent(in) :: weight(:),lower,upper
    real(dp), intent(out) :: new_weight(:)
    logical, intent(in), optional :: strict
    logical, intent(out), optional :: converged
    logical, allocatable :: trimmed(:),outside(:),cantrim(:)
    real(dp) :: trimmings
    integer :: it
    allocate(trimmed(size(weight)),outside(size(weight)),cantrim(size(weight))); trimmed=.false.; new_weight=weight
    do it=1,size(weight)+2
      outside=(new_weight<lower .or. new_weight>upper)
      if(.not.any(outside)) exit
      trimmings=sum(new_weight-merge(max(lower,min(upper,new_weight)),new_weight,outside))
      where(outside) new_weight=max(lower,min(upper,new_weight))
      cantrim=.not.outside .and. .not.trimmed
      if(.not.any(cantrim)) exit
      where(cantrim) new_weight=new_weight+trimmings/real(count(cantrim),dp)
      trimmed=trimmed .or. outside
    end do
    if(present(converged)) converged=.not.any(new_weight<lower .or. new_weight>upper)
    if(present(strict)) then
      if(strict .and. any(new_weight<lower .or. new_weight>upper)) error stop 'trim_weights: strict trimming failed'
    end if
  end subroutine trim_weights

  subroutine calibration_functions(u,cf,lo,hi,fm1,deriv)
    real(dp),intent(in)::u(:),lo,hi;integer,intent(in)::cf;real(dp),intent(out)::fm1(:),deriv(:)
    real(dp)::blo,bup,acoef,eau,pp,mm
    integer::i
    select case(cf)
    case(CAL_LINEAR)
      do i=1,size(u); mm=max(lo,min(hi,u(i)+1)); fm1(i)=mm-1; deriv(i)=merge(1.0_dp,0.0_dp,u(i)<hi-1 .and. u(i)>lo-1); end do
    case(CAL_RAKING)
      do i=1,size(u); mm=exp(max(-700.0_dp,min(700.0_dp,u(i)))); mm=max(lo,min(hi,mm)); fm1(i)=mm-1; deriv(i)=merge(exp(max(-700.0_dp,min(700.0_dp,u(i)))),0.0_dp,u(i)<hi-1 .and. u(i)>lo-1); end do
    case(CAL_LOGIT)
      blo=lo; bup=hi; acoef=(bup-blo)/((bup-1)*(1-blo))
      do i=1,size(u)
        eau=exp(max(-700.0_dp,min(700.0_dp,acoef*u(i))))
        fm1(i)=(blo*(bup-1)+bup*(1-blo)*eau)/(bup-1+(1-blo)*eau)-1
        deriv(i)=bup*(1-blo)*eau*acoef/(bup-1+(1-blo)*eau) - &
          ((blo*(bup-1)+bup*(1-blo)*eau)*((1-blo)*eau*acoef))/(bup-1+(1-blo)*eau)**2
      end do
    case(CAL_SINH)
      do i=1,size(u)
        pp=asinh(2*u(i)); mm=(pp+sqrt(pp*pp+4.0_dp))/2; mm=max(lo,min(hi,mm)); fm1(i)=mm-1
        if(u(i)<hi-1 .and. u(i)>lo-1) then
          deriv(i)=(1.0_dp/sqrt(1.0_dp+(2*u(i))**2))*(1.0_dp+pp/sqrt(pp*pp+4.0_dp))
        else; deriv(i)=0; end if
      end do
    case default
      error stop 'unknown calibration function'
    end select
  end subroutine calibration_functions

  elemental logical function ieee_finite_vec(x) result(ok)
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    real(dp),intent(in)::x; ok=ieee_is_finite(x)
  end function ieee_finite_vec
end module survey_calibration
