! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from tsDyn, copyright its original authors.
! This file may be redistributed and/or modified under GPL version 2 or later.
module tsdyn_rolling
  use tsdyn_kinds, only: dp
  use tsdyn_ar, only: ar_model, fit_ar, forecast_ar
  use tsdyn_setar, only: setar_model, fit_setar, forecast_setar
  use tsdyn_var, only: var_model, fit_var, forecast_var
  implicit none
  private
  public :: rolling_forecast_ar, rolling_forecast_setar, rolling_forecast_var
contains
  subroutine rolling_forecast_ar(x, initial, h, p, include, forecasts, actual, origins, info, window)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: initial,h,p,include
    real(dp), allocatable, intent(out) :: forecasts(:,:),actual(:,:)
    integer, allocatable, intent(out) :: origins(:)
    integer, intent(out) :: info
    integer, intent(in), optional :: window
    integer :: norg,o,start,istat
    type(ar_model)::m
    real(dp),allocatable::f(:)
    if(initial<=p.or.h<1.or.initial+h>size(x))then;info=-1;allocate(forecasts(0,0),actual(0,0),origins(0));return;end if
    norg=size(x)-initial-h+1;allocate(forecasts(h,norg),actual(h,norg),origins(norg))
    do o=1,norg
      origins(o)=initial+o-1;start=1
      if(present(window))start=max(1,origins(o)-max(window,p+2)+1)
      call fit_ar(x(start:origins(o)),p,include,'level',m,istat);if(istat/=0)then;info=istat;return;end if
      call forecast_ar(m,x(start:origins(o)),h,f,istat);if(istat/=0)then;info=istat;return;end if
      forecasts(:,o)=f;actual(:,o)=x(origins(o)+1:origins(o)+h)
    end do
    info=0
  end subroutine rolling_forecast_ar

  subroutine rolling_forecast_setar(x, initial, h, orders, include, nthresh, forecasts, actual, origins, info, thresholds, window)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: initial,h,orders(:),include,nthresh
    real(dp), allocatable, intent(out) :: forecasts(:,:),actual(:,:)
    integer, allocatable, intent(out) :: origins(:)
    integer, intent(out) :: info
    real(dp), intent(in), optional :: thresholds(:)
    integer, intent(in), optional :: window
    integer :: norg,o,start,istat
    type(setar_model)::m
    real(dp),allocatable::f(:)
    if(initial<=maxval(orders)+2.or.h<1.or.initial+h>size(x))then;info=-1;allocate(forecasts(0,0),actual(0,0),origins(0));return;end if
    norg=size(x)-initial-h+1;allocate(forecasts(h,norg),actual(h,norg),origins(norg))
    do o=1,norg
      origins(o)=initial+o-1;start=1
      if(present(window))start=max(1,origins(o)-max(window,maxval(orders)+10)+1)
      if(present(thresholds))then
        call fit_setar(x(start:origins(o)),orders,include,nthresh,m,istat,thresholds=thresholds)
      else
        call fit_setar(x(start:origins(o)),orders,include,nthresh,m,istat)
      end if
      if(istat/=0)then;info=istat;return;end if
      call forecast_setar(m,x(start:origins(o)),h,f,info=istat);if(istat/=0)then;info=istat;return;end if
      forecasts(:,o)=f;actual(:,o)=x(origins(o)+1:origins(o)+h)
    end do
    info=0
  end subroutine rolling_forecast_setar

  subroutine rolling_forecast_var(y, initial, h, p, include, forecasts, actual, origins, info, window)
    real(dp), intent(in) :: y(:,:)
    integer, intent(in) :: initial,h,p,include
    real(dp), allocatable, intent(out) :: forecasts(:,:,:),actual(:,:,:)
    integer, allocatable, intent(out) :: origins(:)
    integer, intent(out) :: info
    integer, intent(in), optional :: window
    integer :: norg,o,start,istat,k
    type(var_model)::m
    real(dp),allocatable::f(:,:)
    k=size(y,2)
    if(initial<=p+2.or.h<1.or.initial+h>size(y,1).or.k<1)then;info=-1;allocate(forecasts(0,0,0),actual(0,0,0),origins(0));return;end if
    norg=size(y,1)-initial-h+1;allocate(forecasts(h,k,norg),actual(h,k,norg),origins(norg))
    do o=1,norg
      origins(o)=initial+o-1;start=1
      if(present(window))start=max(1,origins(o)-max(window,p+5)+1)
      call fit_var(y(start:origins(o),:),p,include,m,istat);if(istat/=0)then;info=istat;return;end if
      call forecast_var(m,y(start:origins(o),:),h,f,istat);if(istat/=0)then;info=istat;return;end if
      forecasts(:,:,o)=f;actual(:,:,o)=y(origins(o)+1:origins(o)+h,:)
    end do
    info=0
  end subroutine rolling_forecast_var
end module tsdyn_rolling
