! SPDX-License-Identifier: GPL-3.0-only
! Native Fortran translation of the computational core of rvinecopulib.
module rvine_fit
  use rvine_kinds, only : dp, pi, clamp_prob
  use rvine_bicop, only : bicop_model, make_bicop, family_bounds, &
                          family_parameter_count, bicop_indep, bicop_gaussian, &
                          bicop_student, bicop_clayton, bicop_gumbel, &
                          bicop_frank, bicop_joe, bicop_bb1, bicop_bb6, &
                          bicop_bb7, bicop_bb8, bicop_tawn
  implicit none
  private
  public :: fit_bicop, select_bicop, empirical_kendall_tau

contains

  real(dp) function empirical_kendall_tau(data) result(tau)
    real(dp), intent(in) :: data(:,:)
    integer :: i, j, n
    real(dp) :: s, dx, dy
    n = size(data,2)
    s = 0.0_dp
    do i = 1, n-1
      do j = i+1, n
        dx = data(1,i)-data(1,j)
        dy = data(2,i)-data(2,j)
        if (dx*dy > 0.0_dp) s = s + 1.0_dp
        if (dx*dy < 0.0_dp) s = s - 1.0_dp
      end do
    end do
    if (n > 1) then
      tau = 2.0_dp*s/real(n*(n-1),dp)
    else
      tau = 0.0_dp
    end if
  end function empirical_kendall_tau

  subroutine fit_bicop(data, family, rotation, model, max_iter, tolerance)
    real(dp), intent(in) :: data(:,:)
    integer, intent(in) :: family
    integer, intent(in), optional :: rotation, max_iter
    real(dp), intent(in), optional :: tolerance
    type(bicop_model), intent(out) :: model
    integer :: rot, np, iter_max
    real(dp) :: tol, lower(3), upper(3), start(3), best(3), tau

    rot = 0
    if (present(rotation)) rot = rotation
    model = make_bicop(family, rot)
    model%nobs = size(data,2)
    np = family_parameter_count(family)
    if (np == 0) then
      model%loglik = 0.0_dp
      return
    end if

    iter_max = 600
    if (present(max_iter)) iter_max = max_iter
    tol = 1.0e-8_dp
    if (present(tolerance)) tol = tolerance
    call family_bounds(family,lower,upper)
    tau = empirical_kendall_tau(data)
    if (rot == 90 .or. rot == 270) tau = -tau
    call starting_parameters(family,tau,start)
    start(1:np) = min(upper(1:np),max(lower(1:np),start(1:np)))
    call nelder_mead_bicop(data,family,rot,start(1:np),lower(1:np),upper(1:np), &
                             best(1:np),iter_max,tol)
    model%parameters(1:np) = best(1:np)
    model%loglik = -bicop_objective(data,family,rot,best(1:np))
  end subroutine fit_bicop

  real(dp) function bicop_objective(data,family,rot,par) result(value)
    real(dp), intent(in) :: data(:,:), par(:)
    integer, intent(in) :: family, rot
    type(bicop_model) :: trial
    integer :: i
    real(dp) :: den
    trial = make_bicop(family,rot,par)
    value = 0.0_dp
    do i = 1, size(data,2)
      den = trial%pdf(clamp_prob(data(1,i)),clamp_prob(data(2,i)))
      if (.not. (den > 0.0_dp) .or. den >= huge(1.0_dp)) then
        value = huge(1.0_dp)/100.0_dp
        return
      end if
      value = value - log(den)
      if (value >= huge(1.0_dp)/1000.0_dp) return
    end do
  end function bicop_objective

  subroutine select_bicop(data, model, families, criterion, allow_rotations)
    real(dp), intent(in) :: data(:,:)
    type(bicop_model), intent(out) :: model
    integer, intent(in), optional :: families(:)
    character(len=*), intent(in), optional :: criterion
    logical, intent(in), optional :: allow_rotations
    integer, allocatable :: fams(:)
    integer :: i, j, nf, nr, rotations(4)
    logical :: do_rot
    character(len=8) :: crit
    real(dp) :: score, best_score
    type(bicop_model) :: trial

    if (present(families)) then
      allocate(fams(size(families)))
      fams = families
    else
      allocate(fams(7))
      fams = [bicop_indep,bicop_gaussian,bicop_student,bicop_clayton, &
              bicop_gumbel,bicop_frank,bicop_joe]
    end if
    crit = 'aic'
    if (present(criterion)) crit = adjustl(criterion)
    do_rot = .true.
    if (present(allow_rotations)) do_rot = allow_rotations
    rotations = [0,90,180,270]
    best_score = huge(1.0_dp)
    model = make_bicop(bicop_indep)
    nf = size(fams)
    do i = 1, nf
      nr = 1
      if (do_rot .and. is_positive_family(fams(i))) nr = 4
      do j = 1, nr
        call fit_bicop(data,fams(i),rotations(j),trial)
        select case (trim(crit))
        case ('bic','BIC')
          score = trial%bic()
        case ('loglik','LOGLIK')
          score = -trial%loglik
        case default
          score = trial%aic()
        end select
        if (score < best_score) then
          best_score = score
          model = trial
        end if
      end do
    end do
  end subroutine select_bicop

  pure logical function is_positive_family(family) result(answer)
    integer, intent(in) :: family
    select case (family)
    case (bicop_clayton,bicop_gumbel,bicop_joe,bicop_bb1,bicop_bb6, &
          bicop_bb7,bicop_bb8,bicop_tawn)
      answer = .true.
    case default
      answer = .false.
    end select
  end function is_positive_family

  subroutine starting_parameters(family,tau,start)
    integer, intent(in) :: family
    real(dp), intent(in) :: tau
    real(dp), intent(out) :: start(3)
    real(dp) :: t
    start = 0.0_dp
    t = min(0.95_dp,max(-0.95_dp,tau))
    select case (family)
    case (bicop_gaussian)
      start(1)=sin(0.5_dp*pi*t)
    case (bicop_student)
      start(1:2)=[sin(0.5_dp*pi*t),5.0_dp]
    case (bicop_clayton)
      t=max(0.01_dp,t)
      start(1)=2.0_dp*t/(1.0_dp-t)
    case (bicop_gumbel)
      start(1)=1.0_dp/(1.0_dp-max(0.0_dp,t))
    case (bicop_frank)
      start(1)=6.0_dp*t/(1.0_dp-abs(t)+0.2_dp)
    case (bicop_joe)
      start(1)=1.0_dp/(1.0_dp-max(0.0_dp,t))
    case (bicop_bb1)
      start(1:2)=[max(0.2_dp,2.0_dp*max(t,0.05_dp)/(1.0_dp-max(t,0.05_dp))),1.2_dp]
    case (bicop_bb6)
      start(1:2)=[1.2_dp,1.2_dp]
    case (bicop_bb7)
      start(1:2)=[1.2_dp,0.8_dp]
    case (bicop_bb8)
      start(1:2)=[1.2_dp,0.8_dp]
    case (bicop_tawn)
      start(1:3)=[1.0_dp,1.0_dp,1.0_dp/(1.0_dp-max(0.0_dp,t))]
    end select
  end subroutine starting_parameters

  subroutine nelder_mead_bicop(data,family,rot,x0,lower,upper,xbest,max_iter,tol)
    real(dp), intent(in) :: data(:,:),x0(:),lower(:),upper(:),tol
    integer, intent(in) :: family,rot
    real(dp), intent(out) :: xbest(:)
    integer, intent(in) :: max_iter
    integer :: n,j,iter,ilo,ihi,inhi
    real(dp), allocatable :: simplex(:,:),vals(:),centroid(:),xr(:),xe(:),xc(:)
    real(dp) :: fr,fe,fc,fspread,step
    real(dp), parameter :: alpha=1.0_dp,gamma=2.0_dp,rho=0.5_dp,sigma=0.5_dp
    n=size(x0)
    allocate(simplex(n,n+1),vals(n+1),centroid(n),xr(n),xe(n),xc(n))
    simplex(:,1)=min(upper,max(lower,x0))
    do j=1,n
      simplex(:,j+1)=simplex(:,1)
      step=0.08_dp*(upper(j)-lower(j))
      if (step <= 0.0_dp) step=0.05_dp*max(1.0_dp,abs(x0(j)))
      simplex(j,j+1)=min(upper(j),max(lower(j),simplex(j,j+1)+step))
      if (abs(simplex(j,j+1)-simplex(j,1)) <= epsilon(1.0_dp)) &
        simplex(j,j+1)=min(upper(j),max(lower(j),simplex(j,j+1)-step))
    end do
    do j=1,n+1
      vals(j)=bicop_objective(data,family,rot,simplex(:,j))
    end do
    do iter=1,max_iter
      ilo=minloc(vals,dim=1)
      ihi=maxloc(vals,dim=1)
      inhi=ilo
      do j=1,n+1
        if (j /= ihi) then
          if (inhi==ilo .or. vals(j)>vals(inhi)) inhi=j
        end if
      end do
      fspread=maxval(abs(vals-vals(ilo)))/(1.0_dp+abs(vals(ilo)))
      if (fspread < tol .and. maxval(abs(simplex-spread(simplex(:,ilo),2,n+1))) < sqrt(tol)) exit
      centroid=0.0_dp
      do j=1,n+1
        if (j /= ihi) centroid=centroid+simplex(:,j)
      end do
      centroid=centroid/real(n,dp)
      xr=min(upper,max(lower,centroid+alpha*(centroid-simplex(:,ihi))))
      fr=bicop_objective(data,family,rot,xr)
      if (fr < vals(ilo)) then
        xe=min(upper,max(lower,centroid+gamma*(xr-centroid)))
        fe=bicop_objective(data,family,rot,xe)
        if (fe < fr) then
          simplex(:,ihi)=xe; vals(ihi)=fe
        else
          simplex(:,ihi)=xr; vals(ihi)=fr
        end if
      else if (fr < vals(inhi)) then
        simplex(:,ihi)=xr; vals(ihi)=fr
      else
        if (fr < vals(ihi)) then
          xc=min(upper,max(lower,centroid+rho*(xr-centroid)))
        else
          xc=min(upper,max(lower,centroid-rho*(centroid-simplex(:,ihi))))
        end if
        fc=bicop_objective(data,family,rot,xc)
        if (fc < min(fr,vals(ihi))) then
          simplex(:,ihi)=xc; vals(ihi)=fc
        else
          do j=1,n+1
            if (j /= ilo) then
              simplex(:,j)=min(upper,max(lower,simplex(:,ilo)+sigma*(simplex(:,j)-simplex(:,ilo))))
              vals(j)=bicop_objective(data,family,rot,simplex(:,j))
            end if
          end do
        end if
      end if
    end do
    ilo=minloc(vals,dim=1)
    xbest=simplex(:,ilo)
  end subroutine nelder_mead_bicop

end module rvine_fit
