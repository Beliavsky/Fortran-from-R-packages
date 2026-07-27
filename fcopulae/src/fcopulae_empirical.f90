! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) fCopulae authors and translation contributors.
! This file is part of fcopulae-modern-fortran and may be redistributed
! and/or modified under the GNU General Public License, version 2 or later.
module fcopulae_empirical
  use fcopulae_kinds, only : dp
  use fcopulae_utils, only : clamp01, debye_function, pseudo_observations
  implicit none
  private
  public :: copula_grid, empirical_copula_grid, empirical_density_grid
  public :: empirical_copula_cdf, frechet_copula_cdf, marshall_olkin_cdf
  public :: debye_function
  type :: copula_grid
    real(dp),allocatable::x(:),y(:),z(:,:)
  end type copula_grid
contains
  function empirical_copula_cdf(u,v,sample_u,sample_v,use_ranks) result(c)
    real(dp),intent(in)::u,v,sample_u(:),sample_v(:)
    logical,intent(in),optional::use_ranks
    real(dp)::c
    real(dp),allocatable::a(:),b(:)
    logical::ranks
    integer::i,n,count
    ranks=.false.;if(present(use_ranks))ranks=use_ranks
    if(ranks)then;call pseudo_observations(sample_u,a);call pseudo_observations(sample_v,b)
    else;allocate(a(size(sample_u)),b(size(sample_v)));a=sample_u;b=sample_v;end if
    n=min(size(a),size(b));count=0
    do i=1,n;if(a(i)<=u.and.b(i)<=v)count=count+1;end do
    if(n>0)then;c=real(count,dp)/real(n,dp);else;c=0.0_dp;end if
  end function empirical_copula_cdf

  function empirical_copula_grid(sample_u,sample_v,n_grid,use_ranks) result(g)
    real(dp),intent(in)::sample_u(:),sample_v(:)
    integer,intent(in),optional::n_grid
    logical,intent(in),optional::use_ranks
    type(copula_grid)::g
    integer::n,i,j
    n=10;if(present(n_grid))n=max(1,n_grid)
    allocate(g%x(n+1),g%y(n+1),g%z(n+1,n+1))
    do i=0,n;g%x(i+1)=real(i,dp)/real(n,dp);g%y(i+1)=g%x(i+1);end do
    do j=1,n+1;do i=1,n+1;g%z(i,j)=empirical_copula_cdf(g%x(i),g%y(j),sample_u,sample_v,use_ranks);end do;end do
  end function empirical_copula_grid

  function empirical_density_grid(sample_u,sample_v,n_grid,use_ranks) result(g)
    real(dp),intent(in)::sample_u(:),sample_v(:)
    integer,intent(in),optional::n_grid
    logical,intent(in),optional::use_ranks
    type(copula_grid)::g
    type(copula_grid)::c
    integer::n,i,j
    n=10;if(present(n_grid))n=max(1,n_grid);c=empirical_copula_grid(sample_u,sample_v,n,use_ranks)
    allocate(g%x(n),g%y(n),g%z(n,n))
    do i=1,n;g%x(i)=0.5_dp*(c%x(i)+c%x(i+1));g%y(i)=0.5_dp*(c%y(i)+c%y(i+1));end do
    do j=1,n;do i=1,n
      g%z(i,j)=real(n*n,dp)*(c%z(i+1,j+1)+c%z(i,j)-c%z(i+1,j)-c%z(i,j+1))
      g%z(i,j)=max(0.0_dp,g%z(i,j))
    end do;end do
  end function empirical_density_grid

  elemental function frechet_copula_cdf(u,v,type_name) result(c)
    real(dp),intent(in)::u,v
    character(len=*),intent(in)::type_name
    real(dp)::c,a,b
    a=clamp01(u);b=clamp01(v)
    select case(trim(adjustl(type_name)))
    case('m','upper');c=min(a,b)
    case('pi','independence');c=a*b
    case('w','lower');c=max(a+b-1.0_dp,0.0_dp)
    case('psp');if(a+b-a*b<=0.0_dp)then;c=0.0_dp;else;c=a*b/(a+b-a*b);end if
    case default;c=a*b
    end select
  end function frechet_copula_cdf

  elemental function marshall_olkin_cdf(u,v,alpha1,alpha2) result(c)
    real(dp),intent(in)::u,v,alpha1,alpha2
    real(dp)::c,a,b
    a=clamp01(u);b=clamp01(v)
    c=min(a**(1.0_dp-alpha1)*b,a*b**(1.0_dp-alpha2))
    c=max(max(0.0_dp,a+b-1.0_dp),min(min(a,b),c))
  end function marshall_olkin_cdf
end module fcopulae_empirical
