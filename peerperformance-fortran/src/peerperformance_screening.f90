! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from PeerPerformance 2.4.0, copyright 2012-2023 David Ardia and Kris Boudt.
module peerperformance_screening
  use peerperformance_kinds, only: dp
  use peerperformance_math, only: finite_value, missing_value
  use peerperformance_types, only: peer_control, test_result, screening_result, &
                                   rolling_result, valid_control
  use peerperformance_stats, only: alpha_coefficients, alpha_testing, sharpe, &
                                   modified_sharpe, sharpe_testing_asymptotic, &
                                   modified_sharpe_testing_asymptotic
  use peerperformance_bootstrap, only: sharpe_testing_bootstrap, &
                                       modified_sharpe_testing_bootstrap
  use peerperformance_pi, only: compute_peer_ratios
  implicit none
  private
  public :: alpha_screening, sharpe_screening, modified_sharpe_screening
  public :: target_peer_performance, roll_screening, exposure_heterogeneity

contains

  subroutine allocate_screening(result, ncoef, nfocal, npeers, cross_group)
    type(screening_result), intent(out) :: result
    integer, intent(in) :: ncoef, nfocal, npeers
    logical, intent(in) :: cross_group
    result%status = 0
    result%message = ''
    result%ncoef = ncoef
    result%n_focal = nfocal
    result%n_peer_group = npeers
    result%cross_group = cross_group
    allocate(result%nobs(nfocal),result%npeer(ncoef,nfocal), &
             result%estimate(ncoef,nfocal), &
             result%difference(ncoef,nfocal,npeers), &
             result%standard_error(ncoef,nfocal,npeers), &
             result%tstat(ncoef,nfocal,npeers),result%pvalue(ncoef,nfocal,npeers), &
             result%lambda(ncoef,nfocal),result%pizero(ncoef,nfocal), &
             result%pipos(ncoef,nfocal),result%pineg(ncoef,nfocal))
    result%nobs = 0
    result%npeer = 0
    result%estimate = missing_value()
    result%difference = missing_value()
    result%standard_error = missing_value()
    result%tstat = missing_value()
    result%pvalue = missing_value()
    result%lambda = missing_value()
    result%pizero = missing_value()
    result%pipos = missing_value()
    result%pineg = missing_value()
  end subroutine allocate_screening

  subroutine apply_peer_ratios(control, result)
    type(peer_control), intent(in) :: control
    type(screening_result), intent(inout) :: result
    real(dp) :: lam(1)
    real(dp), allocatable :: pz(:), pp(:), pn(:), lu(:)
    integer :: c, i
    allocate(pz(result%n_focal),pp(result%n_focal),pn(result%n_focal),lu(result%n_focal))
    do c = 1, result%ncoef
      if (control%has_lambda) then
        lam(1) = control%lambda
        call compute_peer_ratios(result%pvalue(c,:,:),result%difference(c,:,:), &
             result%tstat(c,:,:),pz,pp,pn,lu,lambda=lam,n_boot=control%n_boot, &
             gamma_pos=control%gamma_pos,gamma_neg=control%gamma_neg, &
             seed=control%seed+1000*c,fast=control%fast_adjust)
      else
        call compute_peer_ratios(result%pvalue(c,:,:),result%difference(c,:,:), &
             result%tstat(c,:,:),pz,pp,pn,lu,n_boot=control%n_boot, &
             gamma_pos=control%gamma_pos,gamma_neg=control%gamma_neg, &
             seed=control%seed+1000*c,fast=control%fast_adjust)
      end if
      do i = 1, result%n_focal
        result%npeer(c,i) = count(finite_value(result%pvalue(c,i,:)))
        if (result%npeer(c,i) >= control%min_obs_pi) then
          result%pizero(c,i) = pz(i)
          result%pipos(c,i) = pp(i)
          result%pineg(c,i) = pn(i)
          result%lambda(c,i) = lu(i)
        end if
      end do
    end do
  end subroutine apply_peer_ratios

  subroutine alpha_screening(x, control, result, factors, screen_beta, peers)
    real(dp), intent(in) :: x(:,:)
    type(peer_control), intent(in) :: control
    type(screening_result), intent(out) :: result
    real(dp), intent(in), optional :: factors(:,:), peers(:,:)
    logical, intent(in), optional :: screen_beta
    real(dp), allocatable :: estimates(:,:)
    integer, allocatable :: nobs(:)
    type(test_result) :: test = test_result()
    logical :: all_coef, cross
    integer :: ncoef, nx, ny, i, j, c, status

    if (.not. valid_control(control)) then
      call allocate_screening(result,0,0,0,.false.)
      result%status = 1
      result%message = 'invalid peer-control parameters'
      return
    end if
    nx = size(x,2)
    cross = present(peers)
    if (cross) then
      ny = size(peers,2)
      if (size(peers,1) /= size(x,1)) then
        call allocate_screening(result,0,0,0,.true.)
        result%status = 2
        result%message = 'x and peers must have the same row count'
        return
      end if
    else
      ny = nx
      if (nx < 2) then
        call allocate_screening(result,0,0,0,.false.)
        result%status = 3
        result%message = 'within-group screening requires at least two funds'
        return
      end if
    end if
    if (present(factors)) then
      if (size(factors,1) /= size(x,1)) then
        call allocate_screening(result,0,0,0,cross)
        result%status = 4
        result%message = 'factors must have one row per observation'
        return
      end if
    end if
    all_coef = control%screen_beta
    if (present(screen_beta)) all_coef = screen_beta
    if (.not. present(factors)) all_coef = .false.
    if (all_coef) then
      ncoef = 1+size(factors,2)
    else
      ncoef = 1
    end if
    call allocate_screening(result,ncoef,nx,ny,cross)
    if (present(factors)) then
      call alpha_coefficients(x,estimates,nobs,factors,status)
    else
      call alpha_coefficients(x,estimates,nobs,status=status)
    end if
    result%nobs = nobs
    result%estimate = estimates(1:ncoef,:)

    if (cross) then
      do i = 1, nx
        do j = 1, ny
          if (present(factors)) then
            call alpha_testing(x(:,i),peers(:,j),test,factors,control%hac,all_coef,control%min_obs)
          else
            call alpha_testing(x(:,i),peers(:,j),test,hac=control%hac, &
                               screen_beta=.false.,min_obs=control%min_obs)
          end if
          if (test%status /= 0) cycle
          do c = 1, ncoef
            result%difference(c,i,j) = test%difference(c)
            result%standard_error(c,i,j) = test%standard_error(c)
            result%tstat(c,i,j) = test%tstat(c)
            result%pvalue(c,i,j) = test%pvalue(c)
          end do
        end do
      end do
    else
      do i = 1, nx-1
        do j = i+1, nx
          if (present(factors)) then
            call alpha_testing(x(:,i),x(:,j),test,factors,control%hac,all_coef,control%min_obs)
          else
            call alpha_testing(x(:,i),x(:,j),test,hac=control%hac, &
                               screen_beta=.false.,min_obs=control%min_obs)
          end if
          if (test%status /= 0) cycle
          do c = 1, ncoef
            result%difference(c,i,j) = test%difference(c)
            result%difference(c,j,i) = -test%difference(c)
            result%standard_error(c,i,j) = test%standard_error(c)
            result%standard_error(c,j,i) = test%standard_error(c)
            result%tstat(c,i,j) = test%tstat(c)
            result%tstat(c,j,i) = -test%tstat(c)
            result%pvalue(c,i,j) = test%pvalue(c)
            result%pvalue(c,j,i) = test%pvalue(c)
          end do
        end do
      end do
    end if
    call apply_peer_ratios(control,result)
  end subroutine alpha_screening

  subroutine sharpe_screening(x, control, result, peers)
    real(dp), intent(in) :: x(:,:)
    type(peer_control), intent(in) :: control
    type(screening_result), intent(out) :: result
    real(dp), intent(in), optional :: peers(:,:)
    type(test_result) :: test = test_result()
    real(dp), allocatable :: estimates(:)
    integer, allocatable :: nobs(:)
    logical :: cross
    integer :: nx, ny, i, j, status

    if (.not. valid_control(control)) then
      call allocate_screening(result,0,0,0,.false.)
      result%status=1; result%message='invalid peer-control parameters'; return
    end if
    if (control%test_type == 2 .and. control%block_length == 0) then
      call allocate_screening(result,0,0,0,.false.)
      result%status=2
      result%message='data-driven block length is unavailable in screening; set block_length >= 1'
      return
    end if
    nx=size(x,2); cross=present(peers)
    if (cross) then
      ny=size(peers,2)
      if (size(peers,1)/=size(x,1)) then
        call allocate_screening(result,0,0,0,.true.)
        result%status=3; result%message='x and peers must have the same row count'; return
      end if
    else
      ny=nx
      if (nx<2) then
        call allocate_screening(result,0,0,0,.false.)
        result%status=4; result%message='within-group screening requires at least two funds'; return
      end if
    end if
    call allocate_screening(result,1,nx,ny,cross)
    allocate(estimates(nx),nobs(nx))
    call sharpe(x,estimates,nobs,status)
    result%nobs=nobs; result%estimate(1,:)=estimates
    if (cross) then
      do i=1,nx
        do j=1,ny
          if (control%test_type==1) then
            call sharpe_testing_asymptotic(x(:,i),peers(:,j),test,control%hac, &
                                           control%ttype,control%min_obs)
          else
            call sharpe_testing_bootstrap(x(:,i),peers(:,j),test,control%n_boot, &
                 control%block_length,control%ttype,control%p_boot, &
                 control%seed+10000*i+j,control%min_obs)
          end if
          if (test%status/=0) cycle
          result%difference(1,i,j)=test%difference(1)
          result%standard_error(1,i,j)=test%standard_error(1)
          result%tstat(1,i,j)=test%tstat(1)
          result%pvalue(1,i,j)=test%pvalue(1)
        end do
      end do
    else
      do i=1,nx-1
        do j=i+1,nx
          if (control%test_type==1) then
            call sharpe_testing_asymptotic(x(:,i),x(:,j),test,control%hac, &
                                           control%ttype,control%min_obs)
          else
            call sharpe_testing_bootstrap(x(:,i),x(:,j),test,control%n_boot, &
                 control%block_length,control%ttype,control%p_boot, &
                 control%seed+10000*i+j,control%min_obs)
          end if
          if (test%status/=0) cycle
          result%difference(1,i,j)=test%difference(1)
          result%difference(1,j,i)=-test%difference(1)
          result%standard_error(1,i,j)=test%standard_error(1)
          result%standard_error(1,j,i)=test%standard_error(1)
          result%tstat(1,i,j)=test%tstat(1)
          result%tstat(1,j,i)=-test%tstat(1)
          result%pvalue(1,i,j)=test%pvalue(1)
          result%pvalue(1,j,i)=test%pvalue(1)
        end do
      end do
    end if
    call apply_peer_ratios(control,result)
  end subroutine sharpe_screening

  subroutine modified_sharpe_screening(x, level, control, result, na_negative, peers)
    real(dp), intent(in) :: x(:,:), level
    type(peer_control), intent(in) :: control
    type(screening_result), intent(out) :: result
    logical, intent(in), optional :: na_negative
    real(dp), intent(in), optional :: peers(:,:)
    type(test_result) :: test = test_result()
    real(dp), allocatable :: estimates(:)
    integer, allocatable :: nobs(:)
    logical :: cross, reject_negative
    integer :: nx, ny, i, j, status

    reject_negative=.true.; if (present(na_negative)) reject_negative=na_negative
    if (.not. valid_control(control) .or. level<=0.0_dp .or. level>=1.0_dp) then
      call allocate_screening(result,0,0,0,.false.)
      result%status=1; result%message='invalid parameters'; return
    end if
    if (control%test_type==2 .and. control%block_length==0) then
      call allocate_screening(result,0,0,0,.false.)
      result%status=2
      result%message='data-driven block length is unavailable in screening; set block_length >= 1'
      return
    end if
    nx=size(x,2); cross=present(peers)
    if (cross) then
      ny=size(peers,2)
      if (size(peers,1)/=size(x,1)) then
        call allocate_screening(result,0,0,0,.true.)
        result%status=3; result%message='x and peers must have the same row count'; return
      end if
    else
      ny=nx
      if (nx<2) then
        call allocate_screening(result,0,0,0,.false.)
        result%status=4; result%message='within-group screening requires at least two funds'; return
      end if
    end if
    call allocate_screening(result,1,nx,ny,cross)
    allocate(estimates(nx),nobs(nx))
    call modified_sharpe(x,level,estimates,nobs,reject_negative,status)
    result%nobs=nobs; result%estimate(1,:)=estimates
    if (cross) then
      do i=1,nx
        do j=1,ny
          if (control%test_type==1) then
            call modified_sharpe_testing_asymptotic(x(:,i),peers(:,j),level,test, &
                 reject_negative,control%hac,control%ttype,control%min_obs)
          else
            call modified_sharpe_testing_bootstrap(x(:,i),peers(:,j),level,test, &
                 reject_negative,control%n_boot,control%block_length,control%ttype, &
                 control%p_boot,control%seed+10000*i+j,control%min_obs)
          end if
          if (test%status/=0) cycle
          result%difference(1,i,j)=test%difference(1)
          result%standard_error(1,i,j)=test%standard_error(1)
          result%tstat(1,i,j)=test%tstat(1)
          result%pvalue(1,i,j)=test%pvalue(1)
        end do
      end do
    else
      do i=1,nx-1
        do j=i+1,nx
          if (control%test_type==1) then
            call modified_sharpe_testing_asymptotic(x(:,i),x(:,j),level,test, &
                 reject_negative,control%hac,control%ttype,control%min_obs)
          else
            call modified_sharpe_testing_bootstrap(x(:,i),x(:,j),level,test, &
                 reject_negative,control%n_boot,control%block_length,control%ttype, &
                 control%p_boot,control%seed+10000*i+j,control%min_obs)
          end if
          if (test%status/=0) cycle
          result%difference(1,i,j)=test%difference(1)
          result%difference(1,j,i)=-test%difference(1)
          result%standard_error(1,i,j)=test%standard_error(1)
          result%standard_error(1,j,i)=test%standard_error(1)
          result%tstat(1,i,j)=test%tstat(1)
          result%tstat(1,j,i)=-test%tstat(1)
          result%pvalue(1,i,j)=test%pvalue(1)
          result%pvalue(1,j,i)=test%pvalue(1)
        end do
      end do
    end if
    call apply_peer_ratios(control,result)
  end subroutine modified_sharpe_screening

  subroutine target_peer_performance(x, fund_indices, method, control, result, &
                                     factors, level, na_negative)
    real(dp), intent(in) :: x(:,:)
    integer, intent(in) :: fund_indices(:)
    character(len=*), intent(in) :: method
    type(peer_control), intent(in) :: control
    type(screening_result), intent(out) :: result
    real(dp), intent(in), optional :: factors(:,:), level
    logical, intent(in), optional :: na_negative
    real(dp), allocatable :: focal(:,:)
    real(dp) :: lev
    logical :: reject_negative
    integer :: i
    if (size(fund_indices)<1 .or. any(fund_indices<1) .or. any(fund_indices>size(x,2))) then
      call allocate_screening(result,0,0,0,.true.)
      result%status=1; result%message='fund index out of range'; return
    end if
    do i=1,size(fund_indices)-1
      if (any(fund_indices(i+1:)==fund_indices(i))) then
        call allocate_screening(result,0,0,0,.true.)
        result%status=2; result%message='duplicate fund indices are not allowed'; return
      end if
    end do
    focal=x(:,fund_indices)
    select case(trim(method))
    case('alpha')
      if (present(factors)) then
        call alpha_screening(focal,control,result,factors=factors,peers=x)
      else
        call alpha_screening(focal,control,result,peers=x)
      end if
    case('sharpe')
      call sharpe_screening(focal,control,result,peers=x)
    case('msharpe','modified_sharpe')
      lev=0.9_dp; if (present(level)) lev=level
      reject_negative=.true.; if (present(na_negative)) reject_negative=na_negative
      call modified_sharpe_screening(focal,lev,control,result,reject_negative,peers=x)
    case default
      call allocate_screening(result,0,0,0,.true.)
      result%status=3; result%message='unknown screening method'
    end select
  end subroutine target_peer_performance

  subroutine roll_screening(x, method, width, step, control, result, factors, &
                            peers, level, na_negative)
    real(dp), intent(in) :: x(:,:)
    character(len=*), intent(in) :: method
    integer, intent(in) :: width, step
    type(peer_control), intent(in) :: control
    type(rolling_result), intent(out) :: result
    real(dp), intent(in), optional :: factors(:,:), peers(:,:), level
    logical, intent(in), optional :: na_negative
    type(screening_result) :: screen
    integer :: nw, ncoef, w, start, finish, c, nvalid
    real(dp) :: lev
    logical :: reject_negative
    if (width<2 .or. width>size(x,1) .or. step<1) then
      result%status=1; result%message='invalid rolling width or step'; return
    end if
    nw=1+(size(x,1)-width)/step
    if (trim(method)=='alpha' .and. control%screen_beta .and. present(factors)) then
      ncoef=1+size(factors,2)
    else
      ncoef=1
    end if
    result%nwindow=nw; result%ncoef=ncoef
    allocate(result%window(nw),result%end_index(nw),result%pizero(ncoef,nw), &
             result%pipos(ncoef,nw),result%pineg(ncoef,nw),result%heterogeneity(ncoef,nw))
    result%pizero=missing_value(); result%pipos=missing_value(); result%pineg=missing_value()
    result%heterogeneity=missing_value(); result%status=0; result%message=''
    lev=0.9_dp; if (present(level)) lev=level
    reject_negative=.true.; if (present(na_negative)) reject_negative=na_negative
    do w=1,nw
      start=1+(w-1)*step; finish=start+width-1
      select case(trim(method))
      case('alpha')
        if (present(factors) .and. present(peers)) then
          call alpha_screening(x(start:finish,:),control,screen, &
               factors(start:finish,:),peers=peers(start:finish,:))
        else if (present(factors)) then
          call alpha_screening(x(start:finish,:),control,screen,factors(start:finish,:))
        else if (present(peers)) then
          call alpha_screening(x(start:finish,:),control,screen,peers=peers(start:finish,:))
        else
          call alpha_screening(x(start:finish,:),control,screen)
        end if
      case('sharpe')
        if (present(peers)) then
          call sharpe_screening(x(start:finish,:),control,screen,peers(start:finish,:))
        else
          call sharpe_screening(x(start:finish,:),control,screen)
        end if
      case('msharpe','modified_sharpe')
        if (present(peers)) then
          call modified_sharpe_screening(x(start:finish,:),lev,control,screen, &
                                         reject_negative,peers(start:finish,:))
        else
          call modified_sharpe_screening(x(start:finish,:),lev,control,screen,reject_negative)
        end if
      case default
        result%status=2; result%message='unknown screening method'; return
      end select
      if (screen%status/=0) then
        result%status=screen%status; result%message=screen%message; return
      end if
      result%window(w)=w; result%end_index(w)=finish
      do c=1,ncoef
        nvalid=count(finite_value(screen%pizero(c,:)))
        if (nvalid>0) then
          result%pizero(c,w)=sum(pack(screen%pizero(c,:),finite_value(screen%pizero(c,:))))/real(nvalid,dp)
          result%pipos(c,w)=sum(pack(screen%pipos(c,:),finite_value(screen%pipos(c,:))))/real(nvalid,dp)
          result%pineg(c,w)=sum(pack(screen%pineg(c,:),finite_value(screen%pineg(c,:))))/real(nvalid,dp)
          result%heterogeneity(c,w)=1.0_dp-result%pizero(c,w)
        end if
      end do
    end do
  end subroutine roll_screening

  subroutine exposure_heterogeneity(screen, equal_exposure, heterogeneity, status)
    type(screening_result), intent(in) :: screen
    real(dp), allocatable, intent(out) :: equal_exposure(:), heterogeneity(:)
    integer, intent(out), optional :: status
    integer :: c, nvalid
    allocate(equal_exposure(screen%ncoef),heterogeneity(screen%ncoef))
    equal_exposure=missing_value(); heterogeneity=missing_value()
    if (present(status)) status=0
    do c=1,screen%ncoef
      nvalid=count(finite_value(screen%pizero(c,:)))
      if (nvalid>0) then
        equal_exposure(c)=sum(pack(screen%pizero(c,:),finite_value(screen%pizero(c,:))))/real(nvalid,dp)
        heterogeneity(c)=1.0_dp-equal_exposure(c)
      end if
    end do
  end subroutine exposure_heterogeneity

end module peerperformance_screening
