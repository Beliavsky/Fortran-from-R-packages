module mev_taildep_kj
  use mev_kinds, only: dp
  use mev_math, only: pattern_minimize, finite_diff_hessian, inverse_matrix, empirical_quantile
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  implicit none
  private
  public :: kjtail_result, kjtail, kjtail_uniform

  type :: kjtail_result
    real(dp), allocatable :: p(:)
    real(dp), allocatable :: eta(:)
    real(dp), allocatable :: eta_sd(:)
    real(dp), allocatable :: k1(:)
    real(dp), allocatable :: pat(:)
    real(dp), allocatable :: lambda(:)
    integer, allocatable :: n_tail(:)
    integer, allocatable :: convergence(:)
  end type kjtail_result

contains

  subroutine kjtail(xdat,qlev,result)
    real(dp), intent(in) :: xdat(:,:),qlev(:)
    type(kjtail_result), intent(out) :: result
    real(dp), allocatable :: u(:,:)
    integer :: n,d,j
    n=size(xdat,1); d=size(xdat,2)
    if(n<3 .or. d<2) then
      call empty_result(result,size(qlev)); return
    end if
    allocate(u(n,d))
    do j=1,d
      call rank_uniform(xdat(:,j),u(:,j))
    end do
    call kjtail_uniform(u,qlev,result)
  end subroutine kjtail

  subroutine kjtail_uniform(unif,qlev,result)
    real(dp), intent(in) :: unif(:,:),qlev(:)
    type(kjtail_result), intent(out) :: result
    real(dp), allocatable :: ud(:),sud(:),current_x(:)
    real(dp) :: tp,xopt(2),par(2),fval,h(2,2),hi(2,2),var_eta,tailfrac,current_tp
    integer :: n,d,j,m,nt,info
    n=size(unif,1); d=size(unif,2)
    call empty_result(result,size(qlev))
    if(n<3 .or. d<2 .or. any(unif<=0.0_dp) .or. any(unif>=1.0_dp)) return
    if(any(qlev<=0.0_dp) .or. any(qlev>1.0_dp)) return
    allocate(ud(n))
    do m=1,n
      ud(m)=1.0_dp-minval(unif(m,:))
    end do
    par=[1.0_dp,min(0.999_dp,max(1.0e-4_dp,1.0_dp/(real(d,dp)-0.5_dp)))]
    do j=1,size(qlev)
      result%p(j)=qlev(j)
      tp=empirical_quantile(ud,1.0_dp-qlev(j))
      nt=count(ud<tp); result%n_tail(j)=nt
      if(nt<3 .or. tp<=0.0_dp) cycle
      allocate(sud(nt)); m=0
      do info=1,n
        if(ud(info)<tp) then;m=m+1;sud(m)=ud(info);end if
      end do
      current_tp=tp; current_x=sud
      xopt=[log(max(par(1),1.0e-6_dp)),logit(min(0.999999_dp,max(1.0e-6_dp,par(2))))]
      call pattern_minimize(obj_trans,xopt,fval,info,2500,1.0e-8_dp,0.2_dp)
      call decode(xopt,par)
      result%k1(j)=par(1); result%eta(j)=par(2); result%convergence(j)=info
      tailfrac=real(nt,dp)/real(n,dp); result%pat(j)=tailfrac; result%lambda(j)=tailfrac*par(1)
      call finite_diff_hessian(obj_actual,par,h,2.0e-4_dp)
      call inverse_matrix(h,hi,info)
      if(info==0 .and. hi(2,2)>0.0_dp .and. ieee_is_finite(hi(2,2))) then
        var_eta=hi(2,2);result%eta_sd(j)=sqrt(var_eta)
      end if
      deallocate(sud,current_x)
    end do
  contains
    real(dp) function obj_trans(y) result(v)
      real(dp), intent(in) :: y(:)
      real(dp) :: p(2)
      call decode(y,p); v=obj_actual(p)
    end function obj_trans
    real(dp) function obj_actual(p) result(v)
      real(dp), intent(in) :: p(:)
      real(dp), allocatable :: dens(:),cprob(:)
      real(dp) :: k1,eta
      k1=p(1);eta=p(2)
      if(k1<=0.0_dp .or. eta<=0.0_dp .or. eta>1.0_dp) then;v=1.0e100_dp;return;end if
      allocate(cprob(size(current_x)),dens(size(current_x)))
      cprob=k1*current_x+(1.0_dp-k1*current_tp)*(current_x/current_tp)**(1.0_dp/eta)
      if(minval(cprob)<0.0_dp .or. maxval(cprob)>1.0_dp) then
        v=1.0e8_dp+1.0e6_dp*(max(0.0_dp,-minval(cprob))+max(0.0_dp,maxval(cprob)-1.0_dp));return
      end if
      dens=k1+(1.0_dp-k1*current_tp)/eta*exp((1.0_dp/eta-1.0_dp)*log(current_x)-log(current_tp)/eta)
      if(any(dens<=0.0_dp) .or. any(.not.ieee_is_finite(dens))) then;v=1.0e100_dp;return;end if
      v=-sum(log(max(1.0e-20_dp,dens)))
    end function obj_actual
    subroutine decode(y,p)
      real(dp), intent(in) :: y(:)
      real(dp), intent(out) :: p(2)
      p(1)=exp(min(30.0_dp,max(-30.0_dp,y(1))))
      p(2)=logistic(y(2))
    end subroutine decode
  end subroutine kjtail_uniform

  subroutine rank_uniform(x,u)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: u(size(x))
    real(dp) :: rank
    integer :: i,j,less,equal
    do i=1,size(x)
      less=0;equal=0
      do j=1,size(x)
        if(x(j)<x(i)) less=less+1
        if(abs(x(j)-x(i))<=8.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(x(i)),abs(x(j)))) equal=equal+1
      end do
      rank=real(less,dp)+0.5_dp*real(equal+1,dp)
      u(i)=rank/real(size(x)+1,dp)
    end do
  end subroutine rank_uniform

  subroutine empty_result(r,nq)
    type(kjtail_result), intent(out) :: r
    integer, intent(in) :: nq
    allocate(r%p(nq),r%eta(nq),r%eta_sd(nq),r%k1(nq),r%pat(nq),r%lambda(nq),r%n_tail(nq),r%convergence(nq))
    r%p=ieee_value(0.0_dp,ieee_quiet_nan);r%eta=r%p;r%eta_sd=r%p;r%k1=r%p;r%pat=r%p;r%lambda=r%p
    r%n_tail=0;r%convergence=1
  end subroutine empty_result

  pure real(dp) function logistic(x) result(v)
    real(dp), intent(in) :: x
    if(x>=0.0_dp) then;v=1.0_dp/(1.0_dp+exp(-min(50.0_dp,x)))
    else;v=exp(max(-50.0_dp,x))/(1.0_dp+exp(max(-50.0_dp,x)));end if
  end function logistic

  pure real(dp) function logit(p) result(v)
    real(dp), intent(in) :: p
    real(dp) :: q
    q=min(1.0_dp-1.0e-12_dp,max(1.0e-12_dp,p));v=log(q/(1.0_dp-q))
  end function logit

end module mev_taildep_kj
