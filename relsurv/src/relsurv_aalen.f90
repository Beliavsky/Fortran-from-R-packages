module relsurv_aalen
  use relsurv_kinds, only : dp
  use relsurv_linalg, only : inverse_matrix, solve_linear
  implicit none
  private
  type, public :: aalen_result
    real(dp), allocatable :: time(:), increment(:,:), cumulative(:,:)
    real(dp), allocatable :: variance(:,:,:)
    real(dp), allocatable :: gamma(:)
  end type aalen_result
  public :: at_risk_vector, event_increment_vector, prepare_x
  public :: fit_ols2, aalen_fit, aalen_fit_relative, aalen_fit_const
contains
  pure subroutine at_risk_vector(start,stop,t,y)
    real(dp),intent(in)::start(:),stop(:),t
    integer,intent(out)::y(:)
    integer::i
    do i=1,size(start); y(i)=merge(1,0,t<=stop(i) .and. start(i)<t); end do
  end subroutine at_risk_vector

  pure subroutine event_increment_vector(stop,status,t,dn)
    real(dp),intent(in)::stop(:),t
    integer,intent(in)::status(:)
    real(dp),intent(out)::dn(:)
    integer::i
    do i=1,size(stop); dn(i)=merge(1.0_dp,0.0_dp,stop(i)==t .and. status(i)==1); end do
  end subroutine event_increment_vector

  subroutine prepare_x(y,x,z)
    integer,intent(in)::y(:)
    real(dp),intent(in)::x(:,:)
    real(dp),intent(out)::z(:,:)
    integer::i
    z(:,1)=real(y,dp)
    do i=1,size(x,2); z(:,i+1)=real(y,dp)*x(:,i); end do
  end subroutine prepare_x

  subroutine fit_ols2(x,dn,y,beta,xminus,ok)
    real(dp),intent(in)::x(:,:),dn(:)
    integer,intent(in)::y(:)
    real(dp),intent(out)::beta(:),xminus(:,:)
    logical,intent(out)::ok
    real(dp),allocatable::xtx(:,:),inv(:,:)
    integer::p
    p=size(x,2); beta=0.0_dp; xminus=0.0_dp; ok=.false.
    if(sum(y)<p) return
    allocate(xtx(p,p),inv(p,p)); xtx=matmul(transpose(x),x)
    call inverse_matrix(xtx,inv,ok); if(.not.ok)return
    xminus=matmul(inv,transpose(x)); beta=matmul(xminus,dn)
  end subroutine fit_ols2

  subroutine aalen_fit(start,stop,status,x,times,result,variance,var_estimator)
    real(dp),intent(in)::start(:),stop(:),x(:,:),times(:)
    integer,intent(in)::status(:)
    type(aalen_result),intent(out)::result
    logical,intent(in),optional::variance
    integer,intent(in),optional::var_estimator
    real(dp),allocatable::popinc(:,:)
    allocate(popinc(size(start),size(times))); popinc=0.0_dp
    call aalen_core(start,stop,status,x,times,popinc,result,variance,var_estimator)
  end subroutine aalen_fit

  subroutine aalen_fit_relative(start,stop,status,x,times,pop_increment,result,variance,var_estimator)
    real(dp),intent(in)::start(:),stop(:),x(:,:),times(:),pop_increment(:,:)
    integer,intent(in)::status(:)
    type(aalen_result),intent(out)::result
    logical,intent(in),optional::variance
    integer,intent(in),optional::var_estimator
    call aalen_core(start,stop,status,x,times,pop_increment,result,variance,var_estimator)
  end subroutine aalen_fit_relative

  subroutine aalen_core(start,stop,status,x,times,popinc,result,variance,var_estimator)
    real(dp),intent(in)::start(:),stop(:),x(:,:),times(:),popinc(:,:)
    integer,intent(in)::status(:)
    type(aalen_result),intent(out)::result
    logical,intent(in),optional::variance
    integer,intent(in),optional::var_estimator
    integer,allocatable::yrisk(:)
    real(dp),allocatable::dn(:),adj(:),z(:,:),beta(:),xm(:,:),vinc(:,:),diagv(:)
    logical::ok,dovar
    integer::j,p,n,vest,i
    n=size(start); p=size(x,2)+1; dovar=.false.; if(present(variance))dovar=variance
    vest=1; if(present(var_estimator))vest=var_estimator
    allocate(result%time(size(times)),result%increment(size(times),p),result%cumulative(size(times),p))
    allocate(result%variance(size(times),p,p)); result%time=times; result%increment=0.0_dp
    result%cumulative=0.0_dp; result%variance=0.0_dp
    allocate(yrisk(n),dn(n),adj(n),z(n,p),beta(p),xm(p,n),vinc(p,p),diagv(n))
    do j=1,size(times)
      call at_risk_vector(start,stop,times(j),yrisk)
      call event_increment_vector(stop,status,times(j),dn)
      adj=dn-popinc(:,j)
      call prepare_x(yrisk,x,z)
      call fit_ols2(z,adj,yrisk,beta,xm,ok)
      if(ok) result%increment(j,:)=beta
      if(j==1) then
        result%cumulative(j,:)=result%increment(j,:)
      else
        result%cumulative(j,:)=result%cumulative(j-1,:)+result%increment(j,:)
      end if
      if(dovar .and. ok) then
        if(vest==2) then
          diagv=matmul(z,beta)
        else
          diagv=adj
        end if
        vinc=0.0_dp
        do i=1,n
          vinc=vinc+diagv(i)*spread(xm(:,i),2,p)*spread(xm(:,i),1,p)
        end do
        if(j==1) then
          result%variance(j,:,:)=vinc
        else
          result%variance(j,:,:)=result%variance(j-1,:,:)+vinc
        end if
      else if(j>1) then
        result%variance(j,:,:)=result%variance(j-1,:,:)
      end if
    end do
  end subroutine aalen_core

  subroutine aalen_fit_const(start,stop,status,x,zconst,times,result)
    real(dp),intent(in)::start(:),stop(:),x(:,:),zconst(:,:),times(:)
    integer,intent(in)::status(:)
    type(aalen_result),intent(out)::result
    integer,allocatable::yrisk(:)
    real(dp),allocatable::dn(:),xx(:,:),zz(:,:),xm(:,:),xtx(:,:),inv(:,:),h(:,:),a(:,:),b(:),aint(:,:),bsum(:)
    real(dp),allocatable::gamma(:),beta(:),rhs(:)
    logical::ok
    integer::n,p,q,j
    real(dp)::dt,prev
    n=size(start); p=size(x,2)+1; q=size(zconst,2)
    allocate(yrisk(n),dn(n),xx(n,p),zz(n,q),xm(p,n),xtx(p,p),inv(p,p),h(n,n),a(q,q),b(q))
    allocate(aint(q,q),bsum(q),gamma(q),beta(p),rhs(n)); aint=0.0_dp; bsum=0.0_dp
    prev=minval(start)
    do j=1,size(times)
      call at_risk_vector(start,stop,times(j),yrisk); call event_increment_vector(stop,status,times(j),dn)
      call prepare_x(yrisk,x,xx); zz=spread(real(yrisk,dp),2,q)*zconst
      xtx=matmul(transpose(xx),xx); call inverse_matrix(xtx,inv,ok)
      if(ok) then
        xm=matmul(inv,transpose(xx)); h=0.0_dp; h=-matmul(xx,xm)
        h=h+identity_matrix(n)
        a=matmul(transpose(zz),matmul(h,zz)); b=matmul(transpose(zz),matmul(h,dn))
        dt=times(j)-prev; aint=aint+a*dt; bsum=bsum+b
      end if
      prev=times(j)
    end do
    call solve_linear(aint,bsum,gamma,ok); if(.not.ok)gamma=0.0_dp
    allocate(result%gamma(q)); result%gamma=gamma
    allocate(result%time(size(times)), result%increment(size(times),p), &
      result%cumulative(size(times),p), result%variance(size(times),p,p))
    result%time=times; result%increment=0.0_dp; result%cumulative=0.0_dp
    result%variance=0.0_dp
    prev=minval(start)
    do j=1,size(times)
      call at_risk_vector(start,stop,times(j),yrisk); call event_increment_vector(stop,status,times(j),dn)
      call prepare_x(yrisk,x,xx); xtx=matmul(transpose(xx),xx); call inverse_matrix(xtx,inv,ok)
      if(ok) then
        xm=matmul(inv,transpose(xx)); dt=times(j)-prev; rhs=dn-matmul(zconst,gamma)*dt
        beta=matmul(xm,rhs); result%increment(j,:)=beta
      end if
      if(j==1) then; result%cumulative(j,:)=result%increment(j,:)
      else; result%cumulative(j,:)=result%cumulative(j-1,:)+result%increment(j,:); end if
      prev=times(j)
    end do
  contains
    function identity_matrix(n) result(a)
      integer,intent(in)::n; real(dp)::a(n,n); integer::i
      a=0.0_dp; do i=1,n; a(i,i)=1.0_dp; end do
    end function identity_matrix
  end subroutine aalen_fit_const
end module relsurv_aalen
