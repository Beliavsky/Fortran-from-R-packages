! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) fMultivar authors and translation contributors.
! This file is part of fmultivar-modern-fortran and may be redistributed
! and/or modified under the GNU General Public License, version 2 or later.
program test_utilities
  use fmultivar, only : dp, i8, integration_result, integrate2d_rule, adapt_integrate2d, adapt_integrate_nd, &
    grid_coordinates, grid_data, binning_result, grid2d, make_grid_data, density2d, hist2d, &
    square_binning, hex_binning, rng_state, seed_rng, normal_rng
  implicit none
  type(integration_result)::res
  type(grid_coordinates)::g
  type(grid_data)::d,h,gd
  type(binning_result)::sb,hb
  type(rng_state)::state
  real(dp),allocatable::x(:),y(:),z(:,:)
  real(dp)::dx,dy,total
  integer::i,n,failures
  failures=0
  res=integrate2d_rule(fun_xy,1.0e-7_dp)
  call check_close(res%value,0.25_dp,2.0e-6_dp,'integrate2d rule',failures)
  res=adapt_integrate2d(fun_exp,[0.0_dp,0.0_dp],[1.0_dp,1.0_dp],1.0e-9_dp)
  call check_close(res%value,(exp(1.0_dp)-1.0_dp)**2,1.0e-8_dp,'adaptive 2d',failures)
  res=adapt_integrate_nd(fun_sum,[0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp],2.0e-3_dp,65536)
  call check_close(res%value,1.5_dp,2.0e-3_dp,'adaptive nd',failures)

  g=grid2d([1.0_dp,2.0_dp,3.0_dp],[10.0_dp,20.0_dp])
  call check_true(all(abs(g%x-[1.0_dp,2.0_dp,3.0_dp,1.0_dp,2.0_dp,3.0_dp])<1.0e-14_dp),'grid x',failures)
  call check_true(all(abs(g%y-[10.0_dp,10.0_dp,10.0_dp,20.0_dp,20.0_dp,20.0_dp])<1.0e-14_dp),'grid y',failures)
  allocate(z(2,3));z=reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp],[2,3])
  gd=make_grid_data([1.0_dp,2.0_dp],[3.0_dp,4.0_dp,5.0_dp],z)
  call check_true(all(abs(gd%z-z)<1.0e-14_dp),'grid data',failures)

  n=3000;allocate(x(n),y(n));call seed_rng(state,777_i8)
  do i=1,n;x(i)=normal_rng(state);y(i)=0.5_dp*x(i)+sqrt(0.75_dp)*normal_rng(state);end do
  d=density2d(x,y,35,limits=[-5.0_dp,5.0_dp,-5.0_dp,5.0_dp])
  dx=d%x(2)-d%x(1);dy=d%y(2)-d%y(1);total=sum(d%z)*dx*dy
  call check_close(total,1.0_dp,0.035_dp,'density integral',failures)
  h=hist2d(x,y,[25,30]);call check_close(sum(h%z),real(n,dp),1.0e-12_dp,'hist count',failures)
  sb=square_binning(x,y,24,20);call check_true(sum(sb%count)==n,'square count',failures)
  hb=hex_binning(x,y,25);call check_true(sum(hb%count)==n,'hex count',failures)
  call check_true(all(sb%count>0).and.all(hb%count>0),'positive bin counts',failures)

  if(failures>0)then;write(*,'(a,i0)')'Utility test failures: ',failures;error stop 1;end if
  write(*,'(a)')'Integration, grid, density, histogram, and binning tests passed.'
contains
  function fun_xy(a,b) result(v)
    real(dp),intent(in)::a,b
    real(dp)::v
    v=a*b
  end function fun_xy
  function fun_exp(a,b) result(v)
    real(dp),intent(in)::a,b
    real(dp)::v
    v=exp(a+b)
  end function fun_exp
  function fun_sum(a) result(v)
    real(dp),intent(in)::a(:)
    real(dp)::v
    v=sum(a)
  end function fun_sum
  subroutine check_close(actual,expected,tol,name,nfail)
    real(dp),intent(in)::actual,expected,tol
    character(len=*),intent(in)::name
    integer,intent(inout)::nfail
    if(.not.(abs(actual-expected)<=tol))then
      write(*,'(a,2es18.8)')trim(name)//' failed: ',actual,expected;nfail=nfail+1
    end if
  end subroutine check_close
  subroutine check_true(cond,name,nfail)
    logical,intent(in)::cond
    character(len=*),intent(in)::name
    integer,intent(inout)::nfail
    if(.not.cond)then;write(*,'(a)')trim(name)//' failed';nfail=nfail+1;end if
  end subroutine check_true
end program test_utilities
