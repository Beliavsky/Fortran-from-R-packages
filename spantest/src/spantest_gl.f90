! SPDX-License-Identifier: GPL-3.0-only
module spantest_gl
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use spantest_kinds, only : dp
  use spantest_types, only : gl_result, span_ok, span_invalid_input, span_singular
  use spantest_linalg, only : inverse_matrix, sum_of_squares_columns
  use spantest_random, only : rng_state, rng_seed, rng_uniform
  implicit none
  private
  public :: span_gl_a, span_gl_ad

contains

  pure real(dp) function qnan() result(x)
    x = ieee_value(0.0_dp,ieee_quiet_nan)
  end function qnan

  subroutine fail_gl(res,h0,status,message)
    type(gl_result), intent(out) :: res
    character(len=*), intent(in) :: h0,message
    integer, intent(in) :: status
    res%pval_lmc=qnan(); res%pval_bmc=qnan(); res%stat=qnan()
    res%decision=-1; res%decision_string='Inconclusive'
    res%h0=h0; res%status=status; res%message=message
  end subroutine fail_gl

  integer function rank_last_lex(x,u) result(rank)
    real(dp), intent(in) :: x(:),u(:)
    integer :: i,n
    n=size(x); rank=1
    do i=1,n-1
      if (x(n)>x(i) .or. ((.not. (x(n)>x(i))) .and. (.not. (x(n)<x(i))) .and. u(n)>u(i))) rank=rank+1
    end do
  end function rank_last_lex

  function span_gl_a(r1,r2,totsim,pval_thresh,seed,do_trace) result(res)
    real(dp), intent(in) :: r1(:,:),r2(:,:)
    integer, intent(in), optional :: totsim,seed
    real(dp), intent(in), optional :: pval_thresh
    logical, intent(in), optional :: do_trace
    type(gl_result) :: res
    res=span_gl_core(r1,r2,.false.,totsim,pval_thresh,seed,do_trace)
  end function span_gl_a

  function span_gl_ad(r1,r2,totsim,pval_thresh,seed,do_trace) result(res)
    real(dp), intent(in) :: r1(:,:),r2(:,:)
    integer, intent(in), optional :: totsim,seed
    real(dp), intent(in), optional :: pval_thresh
    logical, intent(in), optional :: do_trace
    type(gl_result) :: res
    res=span_gl_core(r1,r2,.true.,totsim,pval_thresh,seed,do_trace)
  end function span_gl_ad

  function span_gl_core(r1,r2,joint,totsim_in,threshold_in,seed_in,trace_in) result(res)
    real(dp), intent(in) :: r1(:,:),r2(:,:)
    logical, intent(in) :: joint
    integer, intent(in), optional :: totsim_in,seed_in
    real(dp), intent(in), optional :: threshold_in
    logical, intent(in), optional :: trace_in
    type(gl_result) :: res
    real(dp), allocatable :: xx(:,:),xtx(:,:),xinv(:,:),b1(:,:),e1(:,:),ssru(:)
    real(dp), allocatable :: h(:,:),c(:,:),hxht(:,:),hxht_inv(:,:),premult(:,:)
    real(dp), allocatable :: b0(:,:),e0(:,:),ssrr(:),xxb0(:,:),xinv_xt(:,:)
    real(dp), allocatable :: lmc(:),bmc(:),uu(:),esim(:,:),ysim(:,:),b1s(:,:),e1s(:,:)
    real(dp), allocatable :: e0s(:,:),signv(:),temp(:,:),constraint(:,:)
    real(dp) :: threshold, fmax
    integer :: n,k,m,nh,totsim,seed,info,s,i
    logical :: trace
    type(rng_state) :: rng
    character(len=40) :: h0

    if (joint) then
      h0='alpha = 0 and delta = 0'
      nh=2
    else
      h0='alpha = 0'
      nh=1
    end if
    totsim=500; if (present(totsim_in)) totsim=totsim_in
    seed=123; if (present(seed_in)) seed=seed_in
    threshold=0.05_dp; if (present(threshold_in)) threshold=threshold_in
    trace=.false.; if (present(trace_in)) trace=trace_in
    if (size(r1,1)/=size(r2,1) .or. size(r1,1)<1 .or. size(r1,2)<1 .or. &
        size(r2,2)<1 .or. totsim<2 .or. threshold<0.0_dp .or. threshold>1.0_dp) then
      call fail_gl(res,h0,span_invalid_input,'invalid dimensions, simulation count, or threshold')
      return
    end if
    n=size(r1,1); k=size(r1,2); m=size(r2,2)
    allocate(xx(n,k+1),xtx(k+1,k+1),xinv(k+1,k+1),b1(k+1,m),e1(n,m),ssru(m))
    xx(:,1)=1.0_dp; xx(:,2:k+1)=r1
    xtx=matmul(transpose(xx),xx)
    call inverse_matrix(xtx,xinv,info)
    if (info/=0) then
      call fail_gl(res,h0,span_singular,'singular benchmark design')
      return
    end if
    b1=matmul(matmul(xinv,transpose(xx)),r2)
    e1=r2-matmul(xx,b1)
    ssru=sum_of_squares_columns(e1)
    if (any(ssru<=0.0_dp)) then
      call fail_gl(res,h0,span_singular,'nonpositive unrestricted residual sum of squares')
      return
    end if

    allocate(h(nh,k+1),c(nh,m),hxht(nh,nh),hxht_inv(nh,nh),premult(k+1,nh))
    h=0.0_dp; c=0.0_dp; h(1,1)=1.0_dp
    if (joint) then
      h(2,2:k+1)=1.0_dp
      c(2,:)=1.0_dp
    end if
    hxht=matmul(matmul(h,xinv),transpose(h))
    call inverse_matrix(hxht,hxht_inv,info)
    if (info/=0) then
      call fail_gl(res,h0,span_singular,'singular restriction covariance')
      return
    end if
    premult=matmul(matmul(xinv,transpose(h)),hxht_inv)
    allocate(b0(k+1,m),e0(n,m),ssrr(m),xxb0(n,m),xinv_xt(k+1,n))
    b0=b1-matmul(premult,matmul(h,b1)-c)
    xxb0=matmul(xx,b0)
    e0=r2-xxb0
    ssrr=sum_of_squares_columns(e0)
    fmax=maxval((ssrr-ssru)/ssru)

    allocate(lmc(totsim),bmc(totsim),uu(totsim),esim(n,m),ysim(n,m),b1s(k+1,m), &
             e1s(n,m),e0s(n,m),signv(n),temp(n,nh),constraint(nh,m))
    lmc(totsim)=fmax; bmc(totsim)=fmax
    xinv_xt=matmul(xinv,transpose(xx))
    temp=matmul(xx,premult)
    call rng_seed(rng,seed)
    do s=1,totsim-1
      do i=1,n
        if (rng_uniform(rng)<0.5_dp) then
          signv(i)=-1.0_dp
        else
          signv(i)=1.0_dp
        end if
      end do
      esim=e0*spread(signv,2,m)
      ysim=xxb0+esim
      b1s=matmul(xinv_xt,ysim)
      e1s=ysim-matmul(xx,b1s)
      ssru=sum_of_squares_columns(e1s)
      constraint=matmul(h,b1s)-c
      e0s=e1s+matmul(temp,constraint)
      lmc(s)=maxval((sum_of_squares_columns(e0s)-ssru)/ssru)
      bmc(s)=maxval((ssrr-ssru)/ssru)
    end do
    do i=1,totsim
      uu(i)=rng_uniform(rng)
    end do
    res%pval_lmc=real(totsim-rank_last_lex(lmc,uu)+1,dp)/real(totsim,dp)
    res%pval_bmc=real(totsim-rank_last_lex(bmc,uu)+1,dp)/real(totsim,dp)
    res%stat=fmax
    if (res%pval_lmc>threshold) then
      res%decision=1; res%decision_string='Accept'
    else if (res%pval_bmc<=threshold) then
      res%decision=0; res%decision_string='Reject'
    else
      res%decision=-1; res%decision_string='Inconclusive'
    end if
    res%h0=h0; res%status=span_ok
    if (trace) then
      write(*,'(a,es14.6)') 'F-max: ',res%stat
      write(*,'(a,es14.6)') 'LMC p-value: ',res%pval_lmc
      write(*,'(a,es14.6)') 'BMC p-value: ',res%pval_bmc
      write(*,'(a,a)') 'Decision: ',trim(res%decision_string)
    end if
  end function span_gl_core

end module spantest_gl
