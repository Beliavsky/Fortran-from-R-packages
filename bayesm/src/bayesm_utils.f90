module bayesm_utils
  use bayesm_kinds, only: dp
  use bayesm_linalg, only: correlation_matrix
  implicit none
  private
  public :: cond_mom, nmat, num_eff, log_marg_den_nr, cget_c, median_value
contains
  subroutine cond_mom(x,mu,sigi,i,cmean,cvar)
    real(dp), intent(in) :: x(:),mu(:),sigi(:,:)
    integer, intent(in) :: i
    real(dp), intent(out) :: cmean,cvar
    integer :: j
    cvar=1.0_dp/sigi(i,i)
    cmean=mu(i)
    do j=1,size(x)
      if (j/=i) cmean=cmean-(x(j)-mu(j))*sigi(j,i)*cvar
    end do
  end subroutine cond_mom

  pure function nmat(sigma) result(cor)
    real(dp), intent(in) :: sigma(:,:)
    real(dp) :: cor(size(sigma,1),size(sigma,2))
    cor=correlation_matrix(sigma)
  end function nmat

  subroutine num_eff(x,stderr,f,m_used,m)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: stderr,f
    integer, intent(out) :: m_used
    integer, intent(in), optional :: m
    integer :: n,k,lag,mm
    real(dp) :: meanx,varx,covlag,w
    n=size(x)
    if (present(m)) then
      mm=min(max(0,m),n-1)
    else
      mm=min(n-1,int(min(real(n,dp),100.0_dp/sqrt(5000.0_dp)*sqrt(real(n,dp)))))
    end if
    m_used=mm
    meanx=sum(x)/real(max(1,n),dp)
    if (n>1) then
      varx=sum((x-meanx)**2)/real(n-1,dp)
    else
      varx=0.0_dp
    end if
    f=1.0_dp
    if (varx>0.0_dp) then
      do lag=1,mm
        covlag=0.0_dp
        do k=1,n-lag
          covlag=covlag+(x(k)-meanx)*(x(k+lag)-meanx)
        end do
        covlag=covlag/real(n,dp)
        w=real(mm-lag+1,dp)/real(mm+1,dp)
        f=f+2.0_dp*w*covlag/varx
      end do
    end if
    f=max(0.0_dp,f)
    stderr=sqrt(varx*f/real(max(1,n),dp))
  end subroutine num_eff

  real(dp) function median_value(x) result(med)
    real(dp), intent(in) :: x(:)
    real(dp), allocatable :: y(:)
    real(dp) :: tmp
    integer :: i,j,n
    n=size(x); allocate(y(n)); y=x
    do i=2,n
      tmp=y(i); j=i-1
      do while (j>=1)
        if (y(j)<=tmp) exit
        y(j+1)=y(j); j=j-1
      end do
      y(j+1)=tmp
    end do
    if (mod(n,2)==1) then
      med=y((n+1)/2)
    else
      med=0.5_dp*(y(n/2)+y(n/2+1))
    end if
  end function median_value

  real(dp) function log_marg_den_nr(ll) result(v)
    real(dp), intent(in) :: ll(:)
    real(dp) :: med
    med=median_value(ll)
    v=med-log(sum(exp(-ll+med))/real(size(ll),dp))
  end function log_marg_den_nr

  subroutine cget_c(e,k,c,info)
    real(dp), intent(in) :: e
    integer, intent(in) :: k
    real(dp), intent(out) :: c(k+1)
    integer, intent(out) :: info
    real(dp) :: m1,m2,s0,s1,s2,s3,s4,aq,bq,cq,det,a,b,t
    integer :: i,j
    m1=0.0_dp; m2=0.0_dp
    do i=1,k-1
      t=real(i,dp)+0.5_dp; m1=m1+t; m2=m2+t*t
    end do
    s0=real(k-1,dp); s1=0.0_dp; s2=0.0_dp; s3=0.0_dp; s4=0.0_dp
    do i=1,k-1
      t=real(i,dp); s1=s1+t; s2=s2+t*t; s3=s3+t**3; s4=s4+t**4
    end do
    aq=s0*s2-s1*s1
    bq=2.0_dp*e*s0*s3-2.0_dp*e*s1*s2
    cq=m1*m1-m2*s0+e*e*s0*s4-e*e*s2*s2
    det=bq*bq-4.0_dp*aq*cq
    if (det<0.0_dp) then
      info=1; c=0.0_dp; return
    end if
    b=(-bq+sqrt(det))/(2.0_dp*aq)
    a=(m1-b*s1-e*s2)/s0
    c(1)=-1000.0_dp; c(k+1)=1000.0_dp
    do i=1,k-1
      c(i+1)=a+b*real(i,dp)+e*real(i*i,dp)
    end do
    ! insertion sort as upstream returns sort(c)
    do i=2,k+1
      t=c(i); j=i-1
      do while(j>=1)
        if (c(j)<=t) exit
        c(j+1)=c(j); j=j-1
      end do
      c(j+1)=t
    end do
    info=0
  end subroutine cget_c
end module bayesm_utils
