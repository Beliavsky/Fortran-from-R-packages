! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) fMultivar authors and translation contributors.
! This file is part of fmultivar-modern-fortran and may be redistributed
! and/or modified under the GNU General Public License, version 2 or later.
module fmultivar_grid
  use fmultivar_kinds, only : dp, pi
  implicit none
  private
  public :: grid_coordinates, grid_data, binning_result
  public :: grid2d, make_grid_data, density2d, hist2d, square_binning, hex_binning

  type :: grid_coordinates
    real(dp), allocatable :: x(:)
    real(dp), allocatable :: y(:)
  end type grid_coordinates

  type :: grid_data
    real(dp), allocatable :: x(:)
    real(dp), allocatable :: y(:)
    real(dp), allocatable :: z(:,:)
  end type grid_data

  type :: binning_result
    real(dp), allocatable :: x(:)
    real(dp), allocatable :: y(:)
    integer, allocatable :: count(:)
    real(dp), allocatable :: x_center_of_mass(:)
    real(dp), allocatable :: y_center_of_mass(:)
    integer :: bins_x = 0
    integer :: bins_y = 0
  end type binning_result
contains
  function grid2d(x,y) result(grid)
    real(dp),intent(in)::x(:),y(:)
    type(grid_coordinates)::grid
    integer::nx,ny,i,j,k
    nx=size(x);ny=size(y);allocate(grid%x(nx*ny),grid%y(nx*ny));k=0
    do j=1,ny
      do i=1,nx
        k=k+1;grid%x(k)=x(i);grid%y(k)=y(j)
      end do
    end do
  end function grid2d

  function make_grid_data(x,y,z) result(data)
    real(dp),intent(in)::x(:),y(:),z(:,:)
    type(grid_data)::data
    allocate(data%x(size(x)),data%y(size(y)),data%z(size(z,1),size(z,2)))
    data%x=x;data%y=y;data%z=z
  end function make_grid_data

  function density2d(x,y,n,h,limits) result(data)
    real(dp),intent(in)::x(:),y(:)
    integer,intent(in),optional::n
    real(dp),intent(in),optional::h(2),limits(4)
    type(grid_data)::data
    integer::ng,nobs,i,j,k
    real(dp)::band(2),lims(4),dx,dy,ax,ay,norm
    ng=20;if(present(n))ng=n;nobs=size(x)
    lims=[minval(x),maxval(x),minval(y),maxval(y)];if(present(limits))lims=limits
    if(present(h))then;band=h
    else;band=[bandwidth_nrd(x),bandwidth_nrd(y)];end if
    band=max(band/4.0_dp,sqrt(epsilon(1.0_dp)))
    allocate(data%x(ng),data%y(ng),data%z(ng,ng))
    if(ng==1)then;data%x=0.5_dp*(lims(1)+lims(2));data%y=0.5_dp*(lims(3)+lims(4))
    else
      dx=(lims(2)-lims(1))/real(ng-1,dp);dy=(lims(4)-lims(3))/real(ng-1,dp)
      data%x=[(lims(1)+real(i-1,dp)*dx,i=1,ng)]
      data%y=[(lims(3)+real(i-1,dp)*dy,i=1,ng)]
    end if
    norm=1.0_dp/(real(nobs,dp)*band(1)*band(2)*2.0_dp*pi)
    do i=1,ng
      do j=1,ng
        data%z(i,j)=0.0_dp
        do k=1,nobs
          ax=(data%x(i)-x(k))/band(1);ay=(data%y(j)-y(k))/band(2)
          data%z(i,j)=data%z(i,j)+exp(-0.5_dp*(ax*ax+ay*ay))
        end do
        data%z(i,j)=norm*data%z(i,j)
      end do
    end do
  end function density2d

  function hist2d(x,y,bins) result(data)
    real(dp),intent(in)::x(:),y(:)
    integer,intent(in),optional::bins(2)
    type(grid_data)::data
    integer::nb(2),i,ix,iy,nobs
    real(dp)::lo,hi,wx,wy
    nb=[20,20];if(present(bins))nb=bins;nobs=size(x)
    lo=min(minval(x),minval(y));hi=max(maxval(x),maxval(y))
    wx=(hi-lo)/real(nb(1),dp);wy=(hi-lo)/real(nb(2),dp)
    if(wx<=0.0_dp)wx=1.0_dp;if(wy<=0.0_dp)wy=1.0_dp
    allocate(data%x(nb(1)),data%y(nb(2)),data%z(nb(1),nb(2)))
    data%x=[(lo+real(i-1,dp)*wx,i=1,nb(1))]
    data%y=[(lo+real(i-1,dp)*wy,i=1,nb(2))]
    data%z=0.0_dp
    do i=1,nobs
      ix=min(nb(1),max(1,int((x(i)-lo)/wx)+1))
      iy=min(nb(2),max(1,int((y(i)-lo)/wy)+1))
      data%z(ix,iy)=data%z(ix,iy)+1.0_dp
    end do
  end function hist2d

  function square_binning(x,y,bins_x,bins_y) result(out)
    real(dp),intent(in)::x(:),y(:)
    integer,intent(in),optional::bins_x,bins_y
    type(binning_result)::out
    integer::bx,by,n,i,ix,iy,k,nc
    integer,allocatable::cnt(:,:)
    real(dp),allocatable::sx(:,:),sy(:,:)
    real(dp)::xmin,xmax,ymin,ymax,wx,wy
    bx=30;if(present(bins_x))bx=bins_x;by=bx;if(present(bins_y))by=bins_y
    n=size(x);xmin=minval(x);xmax=maxval(x);ymin=minval(y);ymax=maxval(y)
    wx=(xmax-xmin)/real(bx,dp);wy=(ymax-ymin)/real(by,dp)
    if(wx<=0.0_dp)wx=1.0_dp;if(wy<=0.0_dp)wy=1.0_dp
    allocate(cnt(bx,by),sx(bx,by),sy(bx,by));cnt=0;sx=0.0_dp;sy=0.0_dp
    do i=1,n
      ix=min(bx,max(1,int((x(i)-xmin)/wx)+1));iy=min(by,max(1,int((y(i)-ymin)/wy)+1))
      cnt(ix,iy)=cnt(ix,iy)+1;sx(ix,iy)=sx(ix,iy)+x(i);sy(ix,iy)=sy(ix,iy)+y(i)
    end do
    nc=count(cnt>0);allocate(out%x(nc),out%y(nc),out%count(nc), &
      out%x_center_of_mass(nc),out%y_center_of_mass(nc));k=0
    do ix=1,bx
      do iy=1,by
        if(cnt(ix,iy)>0)then
          k=k+1;out%x(k)=xmin+(real(ix,dp)-0.5_dp)*wx
          out%y(k)=ymin+(real(iy,dp)-0.5_dp)*wy;out%count(k)=cnt(ix,iy)
          out%x_center_of_mass(k)=sx(ix,iy)/real(cnt(ix,iy),dp)
          out%y_center_of_mass(k)=sy(ix,iy)/real(cnt(ix,iy),dp)
        end if
      end do
    end do
    out%bins_x=bx;out%bins_y=by
  end function square_binning

  function hex_binning(x,y,bins) result(out)
    real(dp),intent(in)::x(:),y(:)
    integer,intent(in),optional::bins
    type(binning_result)::out
    integer::b,n,jmax,imax,lmax,i,j1,i1,j2,i2,idx,lat,iinc,jinc,nc,k,ii,jj
    integer,allocatable::cnt(:),cell(:)
    real(dp),allocatable::xcm(:),ycm(:)
    real(dp)::xmin,xmax,ymin,ymax,xr,yr,c1,c2,sxv,syv,dist1,test,c3,c4
    b=30;if(present(bins))b=bins;n=size(x)
    xmin=minval(x);xmax=maxval(x);ymin=minval(y);ymax=maxval(y)
    xr=max(xmax-xmin,sqrt(epsilon(1.0_dp)));yr=max(ymax-ymin,sqrt(epsilon(1.0_dp)))
    jmax=floor(real(b,dp)+1.5001_dp)
    c1=2.0_dp*floor(real(b,dp)/sqrt(3.0_dp)+1.5001_dp)
    imax=int((real(jmax,dp)*c1-1.0_dp)/real(jmax,dp)+1.0_dp)
    lmax=jmax*imax;allocate(cnt(lmax),xcm(lmax),ycm(lmax));cnt=0;xcm=0.0_dp;ycm=0.0_dp
    c1=real(b,dp)/xr;c2=real(b,dp)/(yr*sqrt(3.0_dp));jinc=jmax;lat=jinc+1;iinc=2*jinc
    do i=1,n
      sxv=c1*(x(i)-xmin);syv=c2*(y(i)-ymin);j1=floor(sxv+0.5_dp);i1=floor(syv+0.5_dp)
      dist1=(sxv-real(j1,dp))**2+3.0_dp*(syv-real(i1,dp))**2
      if(dist1<0.25_dp)then
        idx=i1*iinc+j1+1
      else if(dist1>1.0_dp/3.0_dp)then
        idx=floor(syv)*iinc+floor(sxv)+lat
      else
        j2=floor(sxv);i2=floor(syv)
        test=(sxv-real(j2,dp)-0.5_dp)**2+3.0_dp*(syv-real(i2,dp)-0.5_dp)**2
        if(dist1<=test)then;idx=i1*iinc+j1+1;else;idx=i2*iinc+j2+lat;end if
      end if
      idx=min(lmax,max(1,idx));cnt(idx)=cnt(idx)+1
      xcm(idx)=xcm(idx)+(x(i)-xcm(idx))/real(cnt(idx),dp)
      ycm(idx)=ycm(idx)+(y(i)-ycm(idx))/real(cnt(idx),dp)
    end do
    nc=count(cnt>0);allocate(cell(nc));k=0
    do i=1,lmax;if(cnt(i)>0)then;k=k+1;cell(k)=i;end if;end do
    allocate(out%x(nc),out%y(nc),out%count(nc),out%x_center_of_mass(nc),out%y_center_of_mass(nc))
    c3=xr/real(b,dp);c4=yr*sqrt(3.0_dp)/(2.0_dp*real(b,dp))
    do k=1,nc
      i=cell(k)-1;ii=i/jmax;jj=mod(i,jmax)
      out%y(k)=c4*real(ii,dp)+ymin
      if(mod(ii,2)==0)then;out%x(k)=c3*real(jj,dp)+xmin
      else;out%x(k)=c3*(real(jj,dp)+0.5_dp)+xmin;end if
      out%count(k)=cnt(cell(k));out%x_center_of_mass(k)=xcm(cell(k));out%y_center_of_mass(k)=ycm(cell(k))
    end do
    out%bins_x=b;out%bins_y=b
  end function hex_binning

  function bandwidth_nrd(x) result(h)
    real(dp),intent(in)::x(:)
    real(dp)::h,q1,q3,sd
    q1=quantile(x,0.25_dp);q3=quantile(x,0.75_dp)
    sd=sqrt(sum((x-sum(x)/real(size(x),dp))**2)/real(max(1,size(x)-1),dp))
    h=4.0_dp*1.06_dp*min(sd,(q3-q1)/1.34_dp)*real(size(x),dp)**(-0.2_dp)
    if(h<=0.0_dp)h=max(sd,1.0_dp)
  end function bandwidth_nrd

  function quantile(x,p) result(q)
    real(dp),intent(in)::x(:),p
    real(dp)::q,pos,w,tmp
    real(dp),allocatable::v(:)
    integer::i,j,k,n
    n=size(x);allocate(v(n));v=x
    do i=2,n
      tmp=v(i);j=i-1
      do while(j>=1)
        if(v(j)<=tmp)exit
        v(j+1)=v(j);j=j-1
      end do
      v(j+1)=tmp
    end do
    pos=1.0_dp+p*real(n-1,dp);k=floor(pos);w=pos-real(k,dp)
    if(k>=n)then;q=v(n);else;q=(1.0_dp-w)*v(k)+w*v(k+1);end if
  end function quantile
end module fmultivar_grid
