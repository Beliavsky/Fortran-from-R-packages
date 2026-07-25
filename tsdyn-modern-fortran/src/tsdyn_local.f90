! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from tsDyn, copyright its original authors.
! This file may be redistributed and/or modified under GPL version 2 or later.
module tsdyn_local
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use tsdyn_kinds, only: dp
  use tsdyn_linalg, only: ols_fit
  use tsdyn_utils, only: lag_embed_univariate
  implicit none
  private
  public :: llar_result, llar_fit_curve, llar_fitted, llar_predict

  type :: llar_result
    real(dp), allocatable :: eps(:), normalized_error(:), average_neighbors(:)
    integer, allocatable :: usable_points(:)
  end type llar_result
contains
  pure real(dp) function embedded_distance(xx,i,j) result(d)
    real(dp),intent(in)::xx(:,:)
    integer,intent(in)::i,j
    d=sqrt(sum((xx(i,:)-xx(j,:))**2))
  end function embedded_distance

  subroutine radius_neighbors_direct(xx,point,eps,theiler,idx,nfound)
    real(dp),intent(in)::xx(:,:),eps
    integer,intent(in)::point,theiler
    integer,allocatable,intent(out)::idx(:)
    integer,intent(out)::nfound
    integer::j,n
    integer,allocatable::tmp(:)
    allocate(tmp(size(xx,1)));n=0
    do j=1,size(xx,1)
      if(abs(j-point)<=theiler)cycle
      if(embedded_distance(xx,point,j)<eps)then;n=n+1;tmp(n)=j;end if
    end do
    allocate(idx(n));if(n>0)idx=tmp(1:n);nfound=n
  end subroutine radius_neighbors_direct

  subroutine radius_neighbors_box(xx,point,eps,theiler,idx,nfound)
    real(dp),intent(in)::xx(:,:),eps
    integer,intent(in)::point,theiler
    integer,allocatable,intent(out)::idx(:)
    integer,intent(out)::nfound
    integer,parameter::maxbox=100
    integer::nb,n,i,j,bx,by,cx,cy,ix,iy,b,lo,hi
    integer,allocatable::counts(:),starts(:),cursor(:),members(:),tmp(:)
    real(dp)::xmin,xmax,ymin,ymax,h
    n=size(xx,1);h=max(eps,tiny(1.0_dp));nb=min(maxbox,max(1,int(1.0_dp/h)+1))
    xmin=minval(xx(:,1));xmax=maxval(xx(:,1));ymin=minval(xx(:,size(xx,2)));ymax=maxval(xx(:,size(xx,2)))
    allocate(counts(nb*nb),starts(nb*nb+1),cursor(nb*nb),members(n),tmp(n));counts=0
    do i=1,n
      bx=min(nb,max(1,1+int((xx(i,1)-xmin)/max(xmax-xmin,tiny(1.0_dp))*real(nb-1,dp))))
      by=min(nb,max(1,1+int((xx(i,size(xx,2))-ymin)/max(ymax-ymin,tiny(1.0_dp))*real(nb-1,dp))))
      counts((bx-1)*nb+by)=counts((bx-1)*nb+by)+1
    end do
    starts(1)=1
    do b=1,nb*nb;starts(b+1)=starts(b)+counts(b);end do
    cursor=starts(1:nb*nb)
    do i=1,n
      bx=min(nb,max(1,1+int((xx(i,1)-xmin)/max(xmax-xmin,tiny(1.0_dp))*real(nb-1,dp))))
      by=min(nb,max(1,1+int((xx(i,size(xx,2))-ymin)/max(ymax-ymin,tiny(1.0_dp))*real(nb-1,dp))))
      b=(bx-1)*nb+by;members(cursor(b))=i;cursor(b)=cursor(b)+1
    end do
    cx=min(nb,max(1,1+int((xx(point,1)-xmin)/max(xmax-xmin,tiny(1.0_dp))*real(nb-1,dp))))
    cy=min(nb,max(1,1+int((xx(point,size(xx,2))-ymin)/max(ymax-ymin,tiny(1.0_dp))*real(nb-1,dp))))
    nfound=0
    do ix=max(1,cx-2),min(nb,cx+2)
      do iy=max(1,cy-2),min(nb,cy+2)
        b=(ix-1)*nb+iy;lo=starts(b);hi=starts(b+1)-1
        do j=lo,hi
          i=members(j);if(abs(i-point)<=theiler)cycle
          if(embedded_distance(xx,point,i)<eps)then;nfound=nfound+1;tmp(nfound)=i;end if
        end do
      end do
    end do
    allocate(idx(nfound));if(nfound>0)idx=tmp(1:nfound)
  end subroutine radius_neighbors_box

  subroutine get_neighbors(xx,point,eps,theiler,method,idx,nfound)
    real(dp),intent(in)::xx(:,:),eps
    integer,intent(in)::point,theiler
    character(len=*),intent(in)::method
    integer,allocatable,intent(out)::idx(:)
    integer,intent(out)::nfound
    character(len=8)::m
    m=adjustl(method)
    if(trim(m)=='auto')then
      if(size(xx,1)>=300.and.size(xx,2)<=8)then;m='box';else;m='direct';end if
    end if
    if(trim(m)=='box')then
      call radius_neighbors_box(xx,point,eps,theiler,idx,nfound)
    else
      call radius_neighbors_direct(xx,point,eps,theiler,idx,nfound)
    end if
  end subroutine get_neighbors

  subroutine local_prediction(xx,yy,point,eps,theiler,method,pred,nfound,info)
    real(dp),intent(in)::xx(:,:),yy(:),eps
    integer,intent(in)::point,theiler
    character(len=*),intent(in)::method
    real(dp),intent(out)::pred
    integer,intent(out)::nfound,info
    integer,allocatable::idx(:)
    real(dp),allocatable::design(:,:),ym(:,:),beta(:,:),fit(:,:),res(:,:)
    integer::rank
    real(dp)::ssr
    call get_neighbors(xx,point,eps,theiler,method,idx,nfound)
    if(nfound<=2*(size(xx,2)+1))then;info=1;pred=0.0_dp;return;end if
    allocate(design(nfound,size(xx,2)+1),ym(nfound,1));design(:,1)=1.0_dp;design(:,2:)=xx(idx,:);ym(:,1)=yy(idx)
    call ols_fit(design,ym,beta,fit,res,rank,ssr,info)
    if(info==0)pred=beta(1,1)+dot_product(xx(point,:),beta(2:,1))
  end subroutine local_prediction

  subroutine scale01(x,z,xmin,xrange)
    real(dp),intent(in)::x(:)
    real(dp),allocatable,intent(out)::z(:)
    real(dp),intent(out)::xmin,xrange
    xmin=minval(x);xrange=maxval(x)-xmin;allocate(z(size(x)))
    if(xrange>0.0_dp)then;z=(x-xmin)/xrange;else;z=0.0_dp;xrange=1.0_dp;end if
  end subroutine scale01

  subroutine llar_fit_curve(x,m,d,steps,eps_values,result,info,search_method)
    real(dp),intent(in)::x(:),eps_values(:)
    integer,intent(in)::m,d,steps
    type(llar_result),intent(out)::result
    integer,intent(out)::info
    character(len=*),intent(in),optional::search_method
    real(dp),allocatable::z(:),xx(:,:),yy(:)
    real(dp)::xmin,xrange,pred,sdy
    integer::istat,i,j,nfound,nok
    character(len=8)::method
    method='auto';if(present(search_method))method=adjustl(search_method)
    call scale01(x,z,xmin,xrange);call lag_embed_univariate(z,m,d,steps,xx,yy,istat)
    if(istat/=0)then;info=istat;return;end if
    allocate(result%eps(size(eps_values)),result%normalized_error(size(eps_values)),result%average_neighbors(size(eps_values)),result%usable_points(size(eps_values)))
    result%eps=eps_values;result%normalized_error=0.0_dp;result%average_neighbors=0.0_dp;result%usable_points=0
    sdy=sqrt(sum((yy-sum(yy)/real(size(yy),dp))**2)/real(max(1,size(yy)-1),dp))
    do i=1,size(eps_values)
      nok=0
      do j=1,size(yy)
        call local_prediction(xx,yy,j,eps_values(i),steps,method,pred,nfound,istat)
        if(istat==0)then
          nok=nok+1;result%normalized_error(i)=result%normalized_error(i)+(pred-yy(j))**2
          result%average_neighbors(i)=result%average_neighbors(i)+real(nfound,dp)
        end if
      end do
      result%usable_points(i)=nok
      if(nok>0)then
        result%normalized_error(i)=sqrt(result%normalized_error(i)/real(nok,dp))/max(sdy,tiny(1.0_dp))
        result%average_neighbors(i)=result%average_neighbors(i)/real(nok,dp)
      end if
    end do
    info=0
  end subroutine llar_fit_curve

  subroutine llar_fitted(x,m,d,steps,eps,fitted,n_neighbors,info,search_method)
    real(dp),intent(in)::x(:),eps
    integer,intent(in)::m,d,steps
    real(dp),allocatable,intent(out)::fitted(:)
    integer,allocatable,intent(out)::n_neighbors(:)
    integer,intent(out)::info
    character(len=*),intent(in),optional::search_method
    real(dp),allocatable::z(:),xx(:,:),yy(:)
    real(dp)::xmin,xrange,pred
    integer::istat,j,nfound
    character(len=8)::method
    method='auto';if(present(search_method))method=adjustl(search_method)
    call scale01(x,z,xmin,xrange);call lag_embed_univariate(z,m,d,steps,xx,yy,istat)
    if(istat/=0)then;info=istat;return;end if
    allocate(fitted(size(yy)),n_neighbors(size(yy)));fitted=ieee_value(0.0_dp,ieee_quiet_nan);n_neighbors=0
    do j=1,size(yy)
      call local_prediction(xx,yy,j,eps,steps,method,pred,nfound,istat);n_neighbors(j)=nfound
      if(istat==0)fitted(j)=pred*xrange+xmin
    end do
    info=0
  end subroutine llar_fitted

  subroutine llar_predict(x,m,d,steps,eps,h,forecast,info,search_method,enlarge_factor)
    real(dp),intent(in)::x(:),eps
    integer,intent(in)::m,d,steps,h
    real(dp),allocatable,intent(out)::forecast(:)
    integer,intent(out)::info
    character(len=*),intent(in),optional::search_method
    real(dp),intent(in),optional::enlarge_factor
    real(dp),allocatable::work(:),z(:),xx(:,:),yy(:),design(:,:),ym(:,:),beta(:,:),fit(:,:),res(:,:)
    real(dp)::xmin,xrange,e,fac,dist
    integer::s,i,j,nfound,rank,istat,nobs
    integer,allocatable::idx(:),tmp(:)
    character(len=8)::method
    real(dp)::ssr
    method='auto';if(present(search_method))method=adjustl(search_method);fac=2.0_dp;if(present(enlarge_factor))fac=max(1.1_dp,enlarge_factor)
    if(h<1)then;info=-1;allocate(forecast(0));return;end if
    allocate(work(size(x)+h));work(1:size(x))=x;allocate(forecast(h))
    do s=1,h
      call scale01(work(1:size(x)+s-1),z,xmin,xrange);call lag_embed_univariate(z,m,d,steps,xx,yy,istat)
      if(istat/=0)then;info=istat;return;end if
      ! Query is the latest m-dimensional state; compare to training rows.
      nobs=size(xx,1);e=eps
      do
        allocate(tmp(nobs));nfound=0
        do i=1,nobs
          dist=0.0_dp
          do j=1,m;dist=dist+(z(size(z)-(j-1)*d)-xx(i,j))**2;end do
          if(sqrt(dist)<=e.and.abs(i-nobs)>=steps)then;nfound=nfound+1;tmp(nfound)=i;end if
        end do
        if(nfound>2*(m+1).or.e>10.0_dp)exit
        e=e*fac
      end do
      if(nfound<=m+1)then;info=2;return;end if
      allocate(idx(nfound));idx=tmp(1:nfound);deallocate(tmp)
      allocate(design(nfound,m+1),ym(nfound,1));design(:,1)=1.0_dp;design(:,2:)=xx(idx,:);ym(:,1)=yy(idx)
      call ols_fit(design,ym,beta,fit,res,rank,ssr,istat);if(istat/=0)then;info=istat;return;end if
      forecast(s)=(beta(1,1)+dot_product(z(size(z):size(z)-(m-1)*d:-d),beta(2:,1)))*xrange+xmin
      work(size(x)+s)=forecast(s)
      deallocate(z,xx,yy,idx,design,ym,beta,fit,res)
    end do
    info=0
  end subroutine llar_predict
end module tsdyn_local
