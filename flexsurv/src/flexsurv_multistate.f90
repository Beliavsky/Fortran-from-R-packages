! SPDX-License-Identifier: GPL-2.0-or-later
module flexsurv_multistate
  use flexsurv_kinds, only : dp
  use flexsurv_fit, only : flexsurv_spec, predict_hazard, predict_cumhaz
  use flexsurv_spline_fit, only : flexsurvspline_result, spline_predict_hazard
  use desolve_types, only : ode_result
  use desolve_rk, only : rk4
  implicit none
  private

  integer, parameter, public :: msm_parametric = 1
  integer, parameter, public :: msm_spline = 2

  type, public :: flexsurv_transition
    integer :: from = 0
    integer :: to = 0
    integer :: model_kind = msm_parametric
    integer :: row = 1
    type(flexsurv_spec) :: spec
    real(dp), allocatable :: theta(:)
    real(dp), allocatable :: covariance(:,:)
    type(flexsurvspline_result) :: spline
    integer, allocatable :: tcov_parameter(:)
    integer, allocatable :: tcov_column(:)
    real(dp), allocatable :: tcov_rate(:)
    real(dp), allocatable :: spline_gamma_tcov_rate(:)
  end type flexsurv_transition

  type, public :: msm_path
    real(dp), allocatable :: time(:)
    integer, allocatable :: state(:)
    integer :: n = 0
  end type msm_path

  public :: transition_intensity, intensity_matrix, pmatrix_flexsurv
  public :: state_occupancy, expected_length_of_stay, simulate_msm_path

contains

  real(dp) function transition_intensity(tr,t) result(h)
    type(flexsurv_transition),intent(in)::tr
    real(dp),intent(in)::t
    select case(tr%model_kind)
    case(msm_parametric)
      h=predict_hazard(tr%spec,tr%theta,tr%row,t)
    case(msm_spline)
      h=spline_predict_hazard(tr%spline,tr%row,t)
    case default
      h=0.0_dp
    end select
    if(h<0.0_dp.or.h/=h)h=0.0_dp
  end function transition_intensity

  subroutine intensity_matrix(trans,nstate,t,q)
    type(flexsurv_transition),intent(in)::trans(:)
    integer,intent(in)::nstate
    real(dp),intent(in)::t
    real(dp),intent(out)::q(nstate,nstate)
    integer::r,i
    real(dp)::h
    q=0.0_dp
    do r=1,size(trans)
      if(trans(r)%from<1.or.trans(r)%from>nstate.or.trans(r)%to<1.or.trans(r)%to>nstate)cycle
      h=transition_intensity(trans(r),t)
      q(trans(r)%from,trans(r)%to)=q(trans(r)%from,trans(r)%to)+h
    end do
    do i=1,nstate;q(i,i)=-sum(q(i,:))+q(i,i);end do
  end subroutine intensity_matrix

  subroutine pmatrix_flexsurv(trans,nstate,times,p,status)
    type(flexsurv_transition),intent(in)::trans(:)
    integer,intent(in)::nstate
    real(dp),intent(in)::times(:)
    real(dp),allocatable,intent(out)::p(:,:,:)
    integer,intent(out),optional::status
    real(dp),allocatable::y0(:)
    type(ode_result)::sol
    integer::i,j,k,idx
    if(size(times)<2)then
      allocate(p(nstate,nstate,size(times)));p=0.0_dp
      if(size(times)==1)then;do i=1,nstate;p(i,i,1)=1.0_dp;end do;end if
      if(present(status))status=0;return
    end if
    allocate(y0(nstate*nstate));y0=0.0_dp
    do i=1,nstate;y0((i-1)*nstate+i)=1.0_dp;end do
    sol=rk4(rhs,y0,times)
    allocate(p(nstate,nstate,size(times)));p=0.0_dp
    do k=1,size(times)
      do i=1,nstate
        do j=1,nstate
          idx=(i-1)*nstate+j;p(i,j,k)=sol%y(idx,k)
        end do
      end do
    end do
    if(present(status))status=sol%status
  contains
    subroutine rhs(t,y,dydt)
      real(dp),intent(in)::t,y(:)
      real(dp),intent(out)::dydt(:)
      real(dp)::q(nstate,nstate),v
      integer::a,b,c,ia,ib
      call intensity_matrix(trans,nstate,t,q);dydt=0.0_dp
      do a=1,nstate
        do b=1,nstate
          v=0.0_dp
          do c=1,nstate
            ia=(a-1)*nstate+c;v=v+y(ia)*q(c,b)
          end do
          ib=(a-1)*nstate+b;dydt(ib)=v
        end do
      end do
    end subroutine rhs
  end subroutine pmatrix_flexsurv

  subroutine state_occupancy(trans,nstate,start_state,times,prob,status)
    type(flexsurv_transition),intent(in)::trans(:)
    integer,intent(in)::nstate,start_state
    real(dp),intent(in)::times(:)
    real(dp),allocatable,intent(out)::prob(:,:)
    integer,intent(out),optional::status
    real(dp),allocatable::p(:,:,:)
    integer::st
    call pmatrix_flexsurv(trans,nstate,times,p,st)
    allocate(prob(nstate,size(times)));prob=p(start_state,:,:)
    if(present(status))status=st
  end subroutine state_occupancy

  subroutine expected_length_of_stay(trans,nstate,start_state,t0,t1,elos,ngrid)
    type(flexsurv_transition),intent(in)::trans(:)
    integer,intent(in)::nstate,start_state
    real(dp),intent(in)::t0,t1
    real(dp),intent(out)::elos(nstate)
    integer,intent(in),optional::ngrid
    real(dp),allocatable::times(:),prob(:,:)
    integer::n,i,j
    n=401;if(present(ngrid))n=max(3,ngrid)
    allocate(times(n));do i=1,n;times(i)=t0+(t1-t0)*real(i-1,dp)/real(n-1,dp);end do
    call state_occupancy(trans,nstate,start_state,times,prob)
    elos=0.0_dp
    do j=1,nstate
      do i=1,n-1
        elos(j)=elos(j)+0.5_dp*(prob(j,i)+prob(j,i+1))*(times(i+1)-times(i))
      end do
    end do
  end subroutine expected_length_of_stay

  function simulate_msm_path(trans,nstate,start_state,t0,tmax,max_events) result(path)
    type(flexsurv_transition),intent(in)::trans(:)
    integer,intent(in)::nstate,start_state
    real(dp),intent(in)::t0,tmax
    integer,intent(in),optional::max_events
    type(msm_path)::path
    real(dp),allocatable::tt(:)
    integer,allocatable::ss(:),outidx(:)
    integer::m,nev,state,nout,j,chosen
    real(dp)::cur,target,u,total,cum
    m=1000;if(present(max_events))m=max_events
    allocate(tt(m+1),ss(m+1));nev=1;tt(1)=t0;ss(1)=start_state;cur=t0;state=start_state
    do while(cur<tmax.and.nev<=m)
      nout=count([(trans(j)%from==state,j=1,size(trans))])
      if(nout==0)exit
      allocate(outidx(nout));nout=0
      do j=1,size(trans);if(trans(j)%from==state)then;nout=nout+1;outidx(nout)=j;end if;end do
      call random_number(u);u=max(u,tiny(1.0_dp));target=-log(u)
      total=cum_total(cur,tmax,outidx)
      if(total<target)then;deallocate(outidx);exit;end if
      cur=find_event_time(cur,tmax,target,outidx)
      total=0.0_dp
      do j=1,nout;total=total+transition_intensity(trans(outidx(j)),cur);end do
      if(total<=0.0_dp)then;deallocate(outidx);exit;end if
      call random_number(u);cum=0.0_dp;chosen=outidx(nout)
      do j=1,nout
        cum=cum+transition_intensity(trans(outidx(j)),cur)/total
        if(u<=cum)then;chosen=outidx(j);exit;end if
      end do
      state=trans(chosen)%to;nev=nev+1;tt(nev)=cur;ss(nev)=state;deallocate(outidx)
    end do
    path%n=nev;allocate(path%time(nev),path%state(nev));path%time=tt(1:nev);path%state=ss(1:nev)
  contains
    real(dp) function cum_total(a,b,idx) result(v)
      real(dp),intent(in)::a,b
      integer,intent(in)::idx(:)
      integer::r,k,nstep
      real(dp)::x0,x1,h0,h1,dx
      nstep=128;dx=(b-a)/real(nstep,dp);v=0.0_dp
      do k=0,nstep-1
        x0=a+real(k,dp)*dx;x1=x0+dx;h0=0.0_dp;h1=0.0_dp
        do r=1,size(idx)
          h0=h0+transition_intensity(trans(idx(r)),x0)
          h1=h1+transition_intensity(trans(idx(r)),x1)
        end do
        v=v+0.5_dp*(h0+h1)*dx
      end do
    end function cum_total
    real(dp) function find_event_time(a,b,z,idx) result(root)
      real(dp),intent(in)::a,b,z
      integer,intent(in)::idx(:)
      real(dp)::lo,hi,mid,val
      integer::it
      lo=a;hi=b
      do it=1,80
        mid=0.5_dp*(lo+hi);val=cum_total(a,mid,idx)
        if(val<z)then;lo=mid;else;hi=mid;end if
      end do
      root=0.5_dp*(lo+hi)
    end function find_event_time
  end function simulate_msm_path

end module flexsurv_multistate
