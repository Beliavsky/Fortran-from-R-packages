module mev_erm
  use mev_kinds, only: dp
  use mev_math, only: sort_descending, pattern_minimize
  implicit none
  private
  public :: erm_result, shape_erm

  type :: erm_result
    integer, allocatable :: k(:)
    real(dp), allocatable :: shape(:)
    real(dp), allocatable :: scale(:)
    real(dp), allocatable :: rho(:)
    integer, allocatable :: convergence(:)
  end type erm_result

contains

  subroutine shape_erm(xdat,kvals,result,method,bounds)
    real(dp), intent(in) :: xdat(:)
    integer, intent(in) :: kvals(:)
    type(erm_result), intent(out) :: result
    character(len=*), intent(in), optional :: method
    real(dp), intent(in), optional :: bounds(2)
    real(dp), allocatable :: xf(:),xs(:),ld(:),z(:)
    real(dp) :: rb(2),x(3),p(3),fval
    integer :: n,kmax,i,j,info,nf,current_k
    character(len=8) :: met
    met='bdgm'; if(present(method)) met=trim(adjustl(method))
    nf=count(xdat>0.0_dp)
    allocate(result%k(size(kvals)),result%shape(size(kvals)),result%scale(size(kvals)), &
             result%rho(size(kvals)),result%convergence(size(kvals)))
    result%k=kvals; result%shape=0.0_dp; result%scale=0.0_dp; result%rho=0.0_dp; result%convergence=1
    if(nf<3 .or. size(kvals)<1) return
    allocate(xf(nf)); j=0
    do i=1,size(xdat)
      if(xdat(i)>0.0_dp) then;j=j+1;xf(j)=xdat(i);end if
    end do
    allocate(xs(nf)); call sort_descending(xf,xs); n=nf
    kmax=maxval(kvals)+1
    if(minval(kvals)<2 .or. kmax>n) return
    allocate(ld(kmax),z(kmax-1)); ld=log(xs(1:kmax))
    do j=1,kmax-1
      z(j)=real(j,dp)*(ld(j)-ld(j+1))
    end do
    if(present(bounds)) then
      rb=sort2(bounds)
    else
      rb=[-5.0_dp,merge(-0.5_dp,-0.25_dp,kmax<5000)]
    end if
    if(rb(2)>=0.0_dp) return
    do i=1,size(kvals)
      current_k=kvals(i)
      x(1)=logit_clamped(min(1.999_dp,max(1.0e-8_dp,sum(z(1:current_k))/real(current_k,dp)))/2.0_dp)
      x(2)=-2.5_dp
      x(3)=0.0_dp
      call pattern_minimize(obj,x,fval,info,2500,1.0e-8_dp,0.18_dp)
      call decode(x,p)
      result%shape(i)=p(1); result%rho(i)=min(0.0_dp,p(3))
      if(met=='fh') then
        result%scale(i)=exp(p(2))/max(p(1),1.0e-12_dp)
      else
        result%scale(i)=exp(p(2))
      end if
      result%convergence(i)=info
    end do
  contains
    real(dp) function obj(y) result(v)
      real(dp), intent(in) :: y(:)
      real(dp) :: par(3),bn,sc
      integer :: m
      call decode(y,par); bn=exp(par(2)); v=0.0_dp
      do m=1,current_k
        if(met=='fh') then
          sc=par(1)*exp(bn*(real(m,dp)/real(current_k+1,dp))**(-par(3)))
        else
          sc=par(1)+bn*(real(m,dp)/real(current_k+1,dp))**(-par(3))
        end if
        if(sc<=0.0_dp .or. sc>1.0e100_dp) then;v=1.0e100_dp;return;end if
        v=v+log(sc)+z(m)/sc
      end do
    end function obj
    subroutine decode(y,par)
      real(dp), intent(in) :: y(:)
      real(dp), intent(out) :: par(3)
      real(dp) :: s
      s=logistic(y(1)); par(1)=2.0_dp*s
      par(2)=-50.0_dp+53.0_dp*logistic(y(2))
      par(3)=rb(1)+(rb(2)-rb(1))*logistic(y(3))
    end subroutine decode
  end subroutine shape_erm

  pure real(dp) function logistic(x) result(v)
    real(dp), intent(in) :: x
    if(x>=0.0_dp) then;v=1.0_dp/(1.0_dp+exp(-min(x,50.0_dp)))
    else;v=exp(max(x,-50.0_dp))/(1.0_dp+exp(max(x,-50.0_dp)));end if
  end function logistic

  pure real(dp) function logit_clamped(p) result(v)
    real(dp), intent(in) :: p
    real(dp) :: q
    q=min(1.0_dp-1.0e-10_dp,max(1.0e-10_dp,p)); v=log(q/(1.0_dp-q))
  end function logit_clamped

  pure function sort2(x) result(y)
    real(dp), intent(in) :: x(2)
    real(dp) :: y(2)
    if(x(1)<=x(2)) then;y=x;else;y=[x(2),x(1)];end if
  end function sort2

end module mev_erm
