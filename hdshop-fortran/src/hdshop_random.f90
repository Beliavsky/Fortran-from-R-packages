! SPDX-License-Identifier: GPL-3.0-only
! Derived from HDShOP 0.1.7.
module hdshop_random
  use, intrinsic :: iso_fortran_env, only: int64
  use hdshop_kinds, only: dp
  use hdshop_linalg, only: symmetric_eigen, symmetrize
  implicit none
  private
  public :: random_covariance_matrix

contains

  function random_covariance_matrix(p,eigenvalues,seed) result(covariance)
    integer,intent(in)::p
    real(dp),intent(in),optional::eigenvalues(:)
    integer(int64),intent(in),optional::seed
    real(dp),allocatable::covariance(:,:),eval(:),z(:,:),wish(:,:),we(:),u(:,:)
    integer::i,j,m
    logical::ok
    real(dp)::u1,u2,spare
    logical::has_spare
    integer(int64)::state
    allocate(covariance(p,p),eval(p));covariance=0.0_dp
    if(p<=0)return
    if(present(eigenvalues))then
      if(size(eigenvalues)/=p .or. any(eigenvalues<=0.0_dp))return
      eval=eigenvalues
    else
      do i=1,p;eval(i)=0.1_dp*exp(5.0_dp*real(i,dp)/real(p,dp));end do
    end if
    m=p*p;allocate(z(p,m),wish(p,p));state=104729_int64;if(present(seed))state=seed
    has_spare=.false.;spare=0.0_dp
    do j=1,m
      do i=1,p
        if(has_spare)then
          z(i,j)=spare;has_spare=.false.
        else
          u1=max(uniform(state),tiny(1.0_dp));u2=uniform(state)
          z(i,j)=sqrt(-2.0_dp*log(u1))*cos(2.0_dp*acos(-1.0_dp)*u2)
          spare=sqrt(-2.0_dp*log(u1))*sin(2.0_dp*acos(-1.0_dp)*u2);has_spare=.true.
        end if
      end do
    end do
    wish=matmul(z,transpose(z));call symmetric_eigen(wish,we,u,ok)
    if(.not.ok)return
    covariance=matmul(u*spread(eval,1,p),transpose(u));call symmetrize(covariance)
  end function random_covariance_matrix

  real(dp) function uniform(state) result(value)
    integer(int64),intent(inout)::state
    integer(int64),parameter::a=2862933555777941757_int64,c=3037000493_int64
    state=a*state+c
    value=real(iand(shiftr(state,11),int(z'001FFFFFFFFFFFFF',int64)),dp)/9007199254740992.0_dp
  end function uniform

end module hdshop_random
