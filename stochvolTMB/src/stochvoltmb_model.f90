! SPDX-License-Identifier: GPL-3.0-only
module stochvoltmb_model
  use stochvoltmb_kinds, only : dp, pi, log_two_pi, tiny_dp, huge_dp
  use stochvoltmb_status, only : sv_ok, sv_invalid_argument, sv_no_convergence, &
                                sv_singular, sv_numerical_failure
  use stochvoltmb_math, only : normal_logpdf, student_t_logpdf, skew_normal_logpdf, &
                              inv_logit_pm1, logit_pm1, clamp, finite_real
  use stochvoltmb_linalg, only : solve_tridiagonal, tridiagonal_logdet, &
                                tridiagonal_inverse_diag, invert_matrix
  use stochvoltmb_types
  implicit none
  private

  public :: parameter_count, parameters_to_theta, theta_to_parameters
  public :: joint_nll, laplace_nll, latent_mode, estimate_parameters, get_nll
  public :: model_name, model_from_name

contains

  pure integer function parameter_count(model) result(p)
    integer, intent(in) :: model
    select case(model)
    case(sv_gaussian)
      p=3
    case(sv_student_t,sv_skew_gaussian,sv_leverage)
      p=4
    case(sv_skew_gaussian_leverage)
      p=5
    case default
      p=0
    end select
  end function parameter_count

  pure function model_name(model) result(name)
    integer, intent(in) :: model
    character(len=32) :: name
    select case(model)
    case(sv_gaussian); name='gaussian'
    case(sv_student_t); name='t'
    case(sv_skew_gaussian); name='skew_gaussian'
    case(sv_leverage); name='leverage'
    case(sv_skew_gaussian_leverage); name='skew_gaussian_leverage'
    case default; name='unknown'
    end select
  end function model_name

  pure integer function model_from_name(name) result(model)
    character(len=*), intent(in) :: name
    select case(trim(adjustl(name)))
    case('gaussian'); model=sv_gaussian
    case('t','student_t'); model=sv_student_t
    case('skew_gaussian'); model=sv_skew_gaussian
    case('leverage'); model=sv_leverage
    case('skew_gaussian_leverage'); model=sv_skew_gaussian_leverage
    case default; model=-1
    end select
  end function model_from_name

  subroutine parameters_to_theta(params, model, theta, info)
    type(sv_parameters), intent(in) :: params
    integer, intent(in) :: model
    real(dp), intent(out) :: theta(:)
    integer, intent(out) :: info
    integer :: p
    p=parameter_count(model); info=sv_ok
    if (p==0 .or. size(theta)/=p .or. params%sigma_y<=0.0_dp .or. &
        params%sigma_h<=0.0_dp .or. abs(params%phi)>=1.0_dp) then
      info=sv_invalid_argument; return
    end if
    theta(1)=log(params%sigma_y)
    theta(2)=log(params%sigma_h)
    theta(3)=logit_pm1(params%phi)
    select case(model)
    case(sv_student_t)
      if (params%df<=2.0_dp) then
        info=sv_invalid_argument; return
      end if
      theta(4)=log(params%df-2.0_dp)
    case(sv_skew_gaussian)
      theta(4)=params%alpha
    case(sv_leverage)
      if (abs(params%rho)>=1.0_dp) then
        info=sv_invalid_argument; return
      end if
      theta(4)=logit_pm1(params%rho)
    case(sv_skew_gaussian_leverage)
      if (abs(params%rho)>=1.0_dp) then
        info=sv_invalid_argument; return
      end if
      theta(4)=params%alpha
      theta(5)=logit_pm1(params%rho)
    end select
  end subroutine parameters_to_theta

  subroutine theta_to_parameters(theta, model, params, info)
    real(dp), intent(in) :: theta(:)
    integer, intent(in) :: model
    type(sv_parameters), intent(out) :: params
    integer, intent(out) :: info
    integer :: p
    p=parameter_count(model); info=sv_ok
    if (p==0 .or. size(theta)/=p .or. any(.not. finite_real(theta))) then
      info=sv_invalid_argument; return
    end if
    params%sigma_y=exp(theta(1))
    params%sigma_h=exp(theta(2))
    params%phi=inv_logit_pm1(theta(3))
    params%df=8.0_dp; params%alpha=0.0_dp; params%rho=0.0_dp
    select case(model)
    case(sv_student_t)
      params%df=exp(theta(4))+2.0_dp
    case(sv_skew_gaussian)
      params%alpha=theta(4)
    case(sv_leverage)
      params%rho=inv_logit_pm1(theta(4))
    case(sv_skew_gaussian_leverage)
      params%alpha=theta(4)
      params%rho=inv_logit_pm1(theta(5))
    end select
    if (.not. finite_real(params%sigma_y) .or. .not. finite_real(params%sigma_h) .or. &
        params%sigma_y<=tiny_dp .or. params%sigma_h<=tiny_dp) info=sv_invalid_argument
  end subroutine theta_to_parameters

  real(dp) function observation_nll(y, h1, h2, params, model, has_next) result(q)
    real(dp), intent(in) :: y, h1, h2
    type(sv_parameters), intent(in) :: params
    integer, intent(in) :: model
    logical, intent(in) :: has_next
    real(dp) :: scale, z, norm_scale, delta, omega, xi, eta, mu, sd
    select case(model)
    case(sv_gaussian)
      scale=params%sigma_y*exp(0.5_dp*h1)
      z=y/scale
      q=log(scale)-normal_logpdf(z)
    case(sv_student_t)
      norm_scale=params%sigma_y*exp(0.5_dp*h1)*sqrt((params%df-2.0_dp)/params%df)
      z=y/norm_scale
      q=log(norm_scale)-student_t_logpdf(z,params%df)
    case(sv_skew_gaussian)
      scale=params%sigma_y*exp(0.5_dp*h1)
      delta=params%alpha/sqrt(1.0_dp+params%alpha*params%alpha)
      omega=scale/sqrt(max(tiny_dp,1.0_dp-2.0_dp*delta*delta/pi))
      xi=-omega*delta*sqrt(2.0_dp/pi)
      q=-skew_normal_logpdf(y,xi,omega,params%alpha)
    case(sv_leverage)
      if (.not. has_next) then
        q=0.0_dp
      else
        scale=params%sigma_y*exp(0.5_dp*h1)
        eta=(h2-params%phi*h1)/params%sigma_h
        mu=scale*params%rho*eta
        sd=scale*sqrt(max(tiny_dp,1.0_dp-params%rho*params%rho))
        q=log(sd)-normal_logpdf((y-mu)/sd)
      end if
    case(sv_skew_gaussian_leverage)
      if (.not. has_next) then
        q=0.0_dp
      else
        scale=params%sigma_y*exp(0.5_dp*h1)
        eta=(h2-params%phi*h1)/params%sigma_h
        mu=scale*params%rho*eta
        sd=scale*sqrt(max(tiny_dp,1.0_dp-params%rho*params%rho))
        delta=params%alpha/sqrt(1.0_dp+params%alpha*params%alpha)
        omega=sd/sqrt(max(tiny_dp,1.0_dp-2.0_dp*delta*delta/pi))
        xi=mu-omega*delta*sqrt(2.0_dp/pi)
        q=-skew_normal_logpdf(y,xi,omega,params%alpha)
      end if
    case default
      q=huge_dp
    end select
    if (.not. finite_real(q)) q=huge_dp
  end function observation_nll

  real(dp) function joint_nll(y, h, params, model) result(nll)
    real(dp), intent(in) :: y(:), h(:)
    type(sv_parameters), intent(in) :: params
    integer, intent(in) :: model
    real(dp) :: var0, e
    integer :: n, i
    n=size(y)
    if (size(h)/=n .or. n<2 .or. params%sigma_y<=0.0_dp .or. params%sigma_h<=0.0_dp .or. &
        abs(params%phi)>=1.0_dp) then
      nll=huge_dp; return
    end if
    var0=params%sigma_h*params%sigma_h/(1.0_dp-params%phi*params%phi)
    nll=0.5_dp*(log_two_pi+log(var0)+h(1)*h(1)/var0)
    do i=2,n
      e=h(i)-params%phi*h(i-1)
      nll=nll+0.5_dp*(log_two_pi+2.0_dp*log(params%sigma_h)+(e/params%sigma_h)**2)
    end do
    do i=1,n-1
      nll=nll+observation_nll(y(i),h(i),h(i+1),params,model,.true.)
    end do
    nll=nll+observation_nll(y(n),h(n),h(n),params,model,.false.)
  end function joint_nll

  subroutine joint_gradient_hessian(y,h,params,model,g,diag,off,nll)
    real(dp), intent(in) :: y(:), h(:)
    type(sv_parameters), intent(in) :: params
    integer, intent(in) :: model
    real(dp), intent(out) :: g(:), diag(:), off(:)
    real(dp), intent(out) :: nll
    real(dp) :: sh2, var0, e, z2, df, eps1, eps2, f0, fp, fm
    real(dp) :: fpp, fpm, fmp, fmm, gx, gy, hxx, hyy, hxy
    integer :: n, i
    n=size(y); g=0.0_dp; diag=0.0_dp; off=0.0_dp
    sh2=params%sigma_h*params%sigma_h
    var0=sh2/(1.0_dp-params%phi*params%phi)
    nll=0.5_dp*(log_two_pi+log(var0)+h(1)*h(1)/var0)
    g(1)=h(1)/var0
    diag(1)=1.0_dp/var0
    do i=2,n
      e=h(i)-params%phi*h(i-1)
      nll=nll+0.5_dp*(log_two_pi+log(sh2)+e*e/sh2)
      g(i)=g(i)+e/sh2
      g(i-1)=g(i-1)-params%phi*e/sh2
      diag(i)=diag(i)+1.0_dp/sh2
      diag(i-1)=diag(i-1)+params%phi*params%phi/sh2
      off(i-1)=off(i-1)-params%phi/sh2
    end do

    select case(model)
    case(sv_gaussian)
      do i=1,n
        z2=(y(i)/(params%sigma_y*exp(0.5_dp*h(i))))**2
        nll=nll+0.5_dp*(log_two_pi+2.0_dp*log(params%sigma_y)+h(i)+z2)
        g(i)=g(i)+0.5_dp*(1.0_dp-z2)
        diag(i)=diag(i)+0.5_dp*z2
      end do
    case(sv_student_t)
      df=params%df
      do i=1,n
        z2=(y(i)/(params%sigma_y*exp(0.5_dp*h(i))*sqrt((df-2.0_dp)/df)))**2
        nll=nll+observation_nll(y(i),h(i),h(i),params,model,.false.)
        g(i)=g(i)+0.5_dp-0.5_dp*(df+1.0_dp)*z2/(df+z2)
        diag(i)=diag(i)+0.5_dp*(df+1.0_dp)*df*z2/(df+z2)**2
      end do
    case(sv_skew_gaussian)
      do i=1,n
        eps1=1.0e-4_dp*max(1.0_dp,abs(h(i)))
        f0=observation_nll(y(i),h(i),h(i),params,model,.false.)
        fp=observation_nll(y(i),h(i)+eps1,h(i),params,model,.false.)
        fm=observation_nll(y(i),h(i)-eps1,h(i),params,model,.false.)
        nll=nll+f0
        g(i)=g(i)+(fp-fm)/(2.0_dp*eps1)
        diag(i)=diag(i)+(fp-2.0_dp*f0+fm)/(eps1*eps1)
      end do
    case(sv_leverage,sv_skew_gaussian_leverage)
      do i=1,n-1
        eps1=1.0e-4_dp*max(1.0_dp,abs(h(i)))
        eps2=1.0e-4_dp*max(1.0_dp,abs(h(i+1)))
        f0=observation_nll(y(i),h(i),h(i+1),params,model,.true.)
        fp=observation_nll(y(i),h(i)+eps1,h(i+1),params,model,.true.)
        fm=observation_nll(y(i),h(i)-eps1,h(i+1),params,model,.true.)
        gx=(fp-fm)/(2.0_dp*eps1)
        hxx=(fp-2.0_dp*f0+fm)/(eps1*eps1)
        fp=observation_nll(y(i),h(i),h(i+1)+eps2,params,model,.true.)
        fm=observation_nll(y(i),h(i),h(i+1)-eps2,params,model,.true.)
        gy=(fp-fm)/(2.0_dp*eps2)
        hyy=(fp-2.0_dp*f0+fm)/(eps2*eps2)
        fpp=observation_nll(y(i),h(i)+eps1,h(i+1)+eps2,params,model,.true.)
        fpm=observation_nll(y(i),h(i)+eps1,h(i+1)-eps2,params,model,.true.)
        fmp=observation_nll(y(i),h(i)-eps1,h(i+1)+eps2,params,model,.true.)
        fmm=observation_nll(y(i),h(i)-eps1,h(i+1)-eps2,params,model,.true.)
        hxy=(fpp-fpm-fmp+fmm)/(4.0_dp*eps1*eps2)
        nll=nll+f0
        g(i)=g(i)+gx; g(i+1)=g(i+1)+gy
        diag(i)=diag(i)+hxx; diag(i+1)=diag(i+1)+hyy
        off(i)=off(i)+hxy
      end do
    end select
  end subroutine joint_gradient_hessian

  subroutine initialize_h(y, params, h)
    real(dp), intent(in) :: y(:)
    type(sv_parameters), intent(in) :: params
    real(dp), intent(out) :: h(:)
    real(dp), allocatable :: raw(:)
    integer :: n, i, lo, hi
    n=size(y); allocate(raw(n))
    raw=log((y/params%sigma_y)**2+0.1_dp)
    raw=clamp(raw,-8.0_dp,8.0_dp)
    do i=1,n
      lo=max(1,i-2); hi=min(n,i+2)
      h(i)=sum(raw(lo:hi))/real(hi-lo+1,dp)
    end do
    h=h-sum(h)/real(n,dp)
  end subroutine initialize_h

  subroutine latent_mode(y, params, model, control, h, nll, diag, off, converged, iterations, info)
    real(dp), intent(in) :: y(:)
    type(sv_parameters), intent(in) :: params
    integer, intent(in) :: model
    type(sv_control), intent(in) :: control
    real(dp), intent(inout) :: h(:)
    real(dp), intent(out) :: nll
    real(dp), intent(out) :: diag(:), off(:)
    logical, intent(out) :: converged
    integer, intent(out) :: iterations, info
    real(dp), allocatable :: g(:), step(:), trial(:), work_diag(:)
    real(dp) :: current, candidate, scale, damping, maxg
    integer :: n, iter, stat, ls
    n=size(y); converged=.false.; info=sv_ok; iterations=0
    allocate(g(n),step(n),trial(n),work_diag(n))
    if (size(h)/=n .or. size(diag)/=n .or. size(off)/=n-1) then
      info=sv_invalid_argument; nll=huge_dp; return
    end if
    if (all(abs(h)<tiny_dp)) call initialize_h(y,params,h)
    do iter=1,control%max_inner_iter
      call joint_gradient_hessian(y,h,params,model,g,diag,off,current)
      maxg=maxval(abs(g))
      if (.not. finite_real(current) .or. any(.not. finite_real(g))) then
        info=sv_numerical_failure; nll=huge_dp; return
      end if
      if (maxg<control%inner_tolerance) then
        converged=.true.; exit
      end if
      damping=0.0_dp
      do
        work_diag=diag+damping
        call solve_tridiagonal(work_diag,off,-g,step,stat)
        if (stat==sv_ok .and. all(finite_real(step))) exit
        if (damping < tiny_dp) then
          damping=1.0e-6_dp
        else
          damping=10.0_dp*damping
        end if
        if (damping>1.0e8_dp) then
          info=sv_singular; nll=huge_dp; return
        end if
      end do
      scale=1.0_dp
      do ls=1,25
        trial=clamp(h+scale*step,-30.0_dp,30.0_dp)
        candidate=joint_nll(y,trial,params,model)
        if (candidate<current) exit
        scale=0.5_dp*scale
      end do
      if (candidate>=current) then
        trial=clamp(h-1.0e-3_dp*g/max(1.0_dp,maxg),-30.0_dp,30.0_dp)
        candidate=joint_nll(y,trial,params,model)
      end if
      h=trial
      if (maxval(abs(scale*step))<control%inner_tolerance*(1.0_dp+maxval(abs(h)))) then
        converged=.true.; exit
      end if
    end do
    iterations=min(iter,control%max_inner_iter)
    call joint_gradient_hessian(y,h,params,model,g,diag,off,nll)
    if (.not. converged) info=sv_no_convergence
  end subroutine latent_mode

  real(dp) function theta_penalty(theta, model) result(pen)
    real(dp), intent(in) :: theta(:)
    integer, intent(in) :: model
    pen=0.0_dp
    if (theta(1)<-20.0_dp .or. theta(1)>10.0_dp) pen=pen+1.0e8_dp*(1.0_dp+abs(theta(1)))
    if (theta(2)<-10.0_dp .or. theta(2)>3.0_dp) pen=pen+1.0e8_dp*(1.0_dp+abs(theta(2)))
    if (abs(theta(3))>12.0_dp) pen=pen+1.0e8_dp*(1.0_dp+abs(theta(3)))
    select case(model)
    case(sv_student_t)
      if (theta(4)<-6.0_dp .or. theta(4)>8.0_dp) pen=pen+1.0e8_dp*(1.0_dp+abs(theta(4)))
    case(sv_skew_gaussian)
      if (abs(theta(4))>30.0_dp) pen=pen+1.0e8_dp*(1.0_dp+abs(theta(4)))
    case(sv_leverage)
      if (abs(theta(4))>12.0_dp) pen=pen+1.0e8_dp*(1.0_dp+abs(theta(4)))
    case(sv_skew_gaussian_leverage)
      if (abs(theta(4))>30.0_dp) pen=pen+1.0e8_dp*(1.0_dp+abs(theta(4)))
      if (abs(theta(5))>12.0_dp) pen=pen+1.0e8_dp*(1.0_dp+abs(theta(5)))
    end select
  end function theta_penalty

  real(dp) function laplace_nll(y, theta, model, control, h_mode, h_diag, h_off, inner_info) result(value)
    real(dp), intent(in) :: y(:), theta(:)
    integer, intent(in) :: model
    type(sv_control), intent(in) :: control
    real(dp), intent(out), optional :: h_mode(:), h_diag(:), h_off(:)
    integer, intent(out), optional :: inner_info
    type(sv_parameters) :: params
    real(dp), allocatable :: h(:), diag(:), off(:)
    real(dp) :: nll, ld, pen
    integer :: n, stat, it
    logical :: conv
    n=size(y); pen=theta_penalty(theta,model)
    if (pen>0.0_dp .or. n<3 .or. any(.not. finite_real(y))) then
      value=1.0e100_dp+pen
      if (present(inner_info)) inner_info=sv_invalid_argument
      return
    end if
    call theta_to_parameters(theta,model,params,stat)
    if (stat/=sv_ok) then
      value=1.0e100_dp
      if (present(inner_info)) inner_info=stat
      return
    end if
    allocate(h(n),diag(n),off(n-1)); h=0.0_dp
    call latent_mode(y,params,model,control,h,nll,diag,off,conv,it,stat)
    if (stat==sv_singular .or. stat==sv_numerical_failure) then
      value=1.0e100_dp
      if (present(inner_info)) inner_info=stat
      return
    end if
    ld=tridiagonal_logdet(diag,off,stat)
    if (stat/=sv_ok .or. .not. finite_real(ld)) then
      value=1.0e100_dp
      if (present(inner_info)) inner_info=stat
      return
    end if
    value=nll+0.5_dp*ld-0.5_dp*real(n,dp)*log_two_pi
    if (present(h_mode)) h_mode=h
    if (present(h_diag)) h_diag=diag
    if (present(h_off)) h_off=off
    if (present(inner_info)) inner_info=merge(sv_ok,sv_no_convergence,conv)
  end function laplace_nll

  real(dp) function get_nll(y, theta, model, control) result(value)
    real(dp), intent(in) :: y(:), theta(:)
    integer, intent(in) :: model
    type(sv_control), intent(in), optional :: control
    type(sv_control) :: ctl
    ctl=sv_control(); if (present(control)) ctl=control
    value=laplace_nll(y,theta,model,ctl)
  end function get_nll

  subroutine order_vertices(f,best,worst,second)
    real(dp), intent(in) :: f(:)
    integer, intent(out) :: best,worst,second
    integer :: j
    best=minloc(f,dim=1)
    worst=maxloc(f,dim=1)
    second=0
    do j=1,size(f)
      if (j==worst) cycle
      if (second==0) then
        second=j
      else if (f(j)>f(second)) then
        second=j
      end if
    end do
  end subroutine order_vertices

  subroutine nelder_mead(y,model,control,x0,xbest,fbest,ok,iters,neval)
    real(dp), intent(in) :: y(:), x0(:)
    integer, intent(in) :: model
    type(sv_control), intent(in) :: control
    real(dp), intent(out) :: xbest(:), fbest
    logical, intent(out) :: ok
    integer, intent(out) :: iters, neval
    integer :: p, j, k, best, worst, second
    real(dp), allocatable :: simplex(:,:), fv(:), centroid(:), xr(:), xe(:), xc(:)
    real(dp) :: fr, fe, fc, spread, step, diameter
    p=size(x0)
    allocate(simplex(p,p+1),fv(p+1),centroid(p),xr(p),xe(p),xc(p))
    simplex(:,1)=x0
    do j=1,p
      simplex(:,j+1)=x0
      step=0.15_dp
      if (j==1 .or. j==2) step=0.25_dp
      if (j==3) step=0.5_dp
      if (abs(x0(j))>1.0_dp) step=0.15_dp*abs(x0(j))
      simplex(j,j+1)=simplex(j,j+1)+step
    end do
    neval=0
    do j=1,p+1
      fv(j)=laplace_nll(y,simplex(:,j),model,control)
      neval=neval+1
    end do
    ok=.false.
    do k=1,control%max_outer_iter
      call order_vertices(fv,best,worst,second)
      spread=maxval(abs(fv-fv(best)))/(1.0_dp+abs(fv(best)))
      diameter=0.0_dp
      do j=1,p+1
        diameter=max(diameter,maxval(abs(simplex(:,j)-simplex(:,best))))
      end do
      if (spread<control%outer_tolerance .and. diameter<5.0e-4_dp) then
        ok=.true.
        exit
      end if
      centroid=0.0_dp
      do j=1,p+1
        if (j/=worst) centroid=centroid+simplex(:,j)
      end do
      centroid=centroid/real(p,dp)
      xr=centroid+(centroid-simplex(:,worst))
      fr=laplace_nll(y,xr,model,control); neval=neval+1
      if (fr<fv(best)) then
        xe=centroid+2.0_dp*(xr-centroid)
        fe=laplace_nll(y,xe,model,control); neval=neval+1
        if (fe<fr) then
          simplex(:,worst)=xe; fv(worst)=fe
        else
          simplex(:,worst)=xr; fv(worst)=fr
        end if
      else if (fr<fv(second)) then
        simplex(:,worst)=xr; fv(worst)=fr
      else
        if (fr<fv(worst)) then
          xc=centroid+0.5_dp*(xr-centroid)
        else
          xc=centroid+0.5_dp*(simplex(:,worst)-centroid)
        end if
        fc=laplace_nll(y,xc,model,control); neval=neval+1
        if (fc<min(fr,fv(worst))) then
          simplex(:,worst)=xc; fv(worst)=fc
        else
          do j=1,p+1
            if (j/=best) then
              simplex(:,j)=simplex(:,best)+0.5_dp*(simplex(:,j)-simplex(:,best))
              fv(j)=laplace_nll(y,simplex(:,j),model,control)
              neval=neval+1
            end if
          end do
        end if
      end if
      if (control%verbose .and. mod(k,10)==0) then
        write(*,'(a,i0,a,es13.5)') 'outer iteration ',k,', objective ',minval(fv)
      end if
    end do
    call order_vertices(fv,best,worst,second)
    xbest=simplex(:,best)
    fbest=fv(best)
    iters=min(k,control%max_outer_iter)
  end subroutine nelder_mead

  subroutine numerical_hessian(y,theta,model,control,hess,info)
    real(dp), intent(in) :: y(:), theta(:)
    integer, intent(in) :: model
    type(sv_control), intent(in) :: control
    real(dp), intent(out) :: hess(:,:)
    integer, intent(out) :: info
    real(dp), allocatable :: xp(:), xm(:), xpp(:), xpm(:), xmp(:), xmm(:)
    real(dp) :: f0, fp, fm, fpp, fpm, fmp, fmm, ei, ej
    integer :: p, i, j
    p=size(theta); info=sv_ok; hess=0.0_dp
    allocate(xp(p),xm(p),xpp(p),xpm(p),xmp(p),xmm(p))
    f0=laplace_nll(y,theta,model,control)
    if (.not. finite_real(f0) .or. f0>1.0e90_dp) then
      info=sv_numerical_failure; return
    end if
    do i=1,p
      ei=2.0e-3_dp*max(1.0_dp,abs(theta(i)))
      xp=theta; xm=theta; xp(i)=xp(i)+ei; xm(i)=xm(i)-ei
      fp=laplace_nll(y,xp,model,control)
      fm=laplace_nll(y,xm,model,control)
      hess(i,i)=(fp-2.0_dp*f0+fm)/(ei*ei)
      do j=i+1,p
        ej=2.0e-3_dp*max(1.0_dp,abs(theta(j)))
        xpp=theta; xpm=theta; xmp=theta; xmm=theta
        xpp(i)=xpp(i)+ei; xpp(j)=xpp(j)+ej
        xpm(i)=xpm(i)+ei; xpm(j)=xpm(j)-ej
        xmp(i)=xmp(i)-ei; xmp(j)=xmp(j)+ej
        xmm(i)=xmm(i)-ei; xmm(j)=xmm(j)-ej
        fpp=laplace_nll(y,xpp,model,control)
        fpm=laplace_nll(y,xpm,model,control)
        fmp=laplace_nll(y,xmp,model,control)
        fmm=laplace_nll(y,xmm,model,control)
        hess(i,j)=(fpp-fpm-fmp+fmm)/(4.0_dp*ei*ej)
        hess(j,i)=hess(i,j)
      end do
    end do
    if (any(.not. finite_real(hess))) info=sv_numerical_failure
  end subroutine numerical_hessian

  subroutine transformed_standard_errors(theta,model,theta_se,param_se)
    real(dp), intent(in) :: theta(:), theta_se(:)
    integer, intent(in) :: model
    real(dp), intent(out) :: param_se(:)
    type(sv_parameters) :: par
    real(dp) :: deriv
    integer :: stat
    call theta_to_parameters(theta,model,par,stat)
    param_se=0.0_dp
    param_se(1)=par%sigma_y*theta_se(1)
    param_se(2)=par%sigma_h*theta_se(2)
    deriv=0.5_dp*(1.0_dp-par%phi*par%phi)
    param_se(3)=abs(deriv)*theta_se(3)
    select case(model)
    case(sv_student_t)
      param_se(4)=(par%df-2.0_dp)*theta_se(4)
    case(sv_skew_gaussian)
      param_se(4)=theta_se(4)
    case(sv_leverage)
      deriv=0.5_dp*(1.0_dp-par%rho*par%rho)
      param_se(4)=abs(deriv)*theta_se(4)
    case(sv_skew_gaussian_leverage)
      param_se(4)=theta_se(4)
      deriv=0.5_dp*(1.0_dp-par%rho*par%rho)
      param_se(5)=abs(deriv)*theta_se(5)
    end select
  end subroutine transformed_standard_errors

  subroutine estimate_parameters(y, model, fit, start, control)
    real(dp), intent(in) :: y(:)
    integer, intent(in) :: model
    type(sv_fit_result), intent(out) :: fit
    type(sv_parameters), intent(in), optional :: start
    type(sv_control), intent(in), optional :: control
    type(sv_control) :: ctl
    type(sv_parameters) :: initial
    real(dp), allocatable :: theta0(:), theta(:), hdiag(:), hoff(:), invhdiag(:)
    real(dp), allocatable :: hess(:,:), cov(:,:)
    real(dp) :: sy, obj
    integer :: p, stat, n, evals, iter, inner_stat, i
    logical :: conv
    ctl=sv_control()
    if (present(control)) ctl=control
    n=size(y); p=parameter_count(model)
    fit%nobs=n; fit%model=model
    if (n<10 .or. p==0 .or. any(.not. finite_real(y))) then
      fit%status=sv_invalid_argument; fit%message='invalid data or model'; return
    end if
    sy=sqrt(sum((y-sum(y)/real(n,dp))**2)/real(max(1,n-1),dp))
    initial%sigma_y=max(1.0e-6_dp,sy)
    initial%sigma_h=0.25_dp
    initial%phi=0.95_dp
    initial%df=8.0_dp
    initial%alpha=0.0_dp
    initial%rho=-0.4_dp
    if (present(start)) initial=start
    allocate(theta0(p),theta(p))
    call parameters_to_theta(initial,model,theta0,stat)
    if (stat/=sv_ok) then
      fit%status=stat; fit%message='invalid starting parameters'; return
    end if
    call nelder_mead(y,model,ctl,theta0,theta,obj,conv,iter,evals)
    allocate(fit%theta(p)); fit%theta=theta
    fit%objective=obj; fit%log_likelihood=-obj; fit%iterations=iter
    fit%function_evaluations=evals; fit%converged=conv
    call theta_to_parameters(theta,model,fit%params,stat)
    allocate(fit%h(n),hdiag(n),hoff(n-1),fit%h_se(n),invhdiag(n))
    obj=laplace_nll(y,theta,model,ctl,fit%h,hdiag,hoff,inner_stat)
    fit%objective=obj; fit%log_likelihood=-obj
    call tridiagonal_inverse_diag(hdiag,hoff,invhdiag,stat)
    if (stat==sv_ok) then
      fit%h_se=sqrt(max(0.0_dp,invhdiag))
    else
      fit%h_se=0.0_dp
    end if
    allocate(fit%theta_cov(p,p),fit%theta_se(p),fit%param_se(p))
    fit%theta_cov=0.0_dp; fit%theta_se=0.0_dp; fit%param_se=0.0_dp
    if (ctl%compute_covariance) then
      allocate(hess(p,p),cov(p,p))
      call numerical_hessian(y,theta,model,ctl,hess,stat)
      if (stat==sv_ok) then
        do i=1,p
          hess(i,i)=hess(i,i)+1.0e-8_dp
        end do
        call invert_matrix(hess,cov,stat)
      end if
      if (stat==sv_ok) then
        fit%theta_cov=cov
        do i=1,p
          fit%theta_se(i)=sqrt(max(0.0_dp,cov(i,i)))
        end do
        call transformed_standard_errors(theta,model,fit%theta_se,fit%param_se)
      end if
    end if
    if (conv) then
      fit%status=sv_ok; fit%message='converged'
    else
      fit%status=sv_no_convergence; fit%message='outer optimizer reached its iteration limit'
    end if
    if (inner_stat/=sv_ok .and. fit%status==sv_ok) then
      fit%message='outer optimizer converged; latent mode reached tolerance approximately'
    end if
  end subroutine estimate_parameters

end module stochvoltmb_model
