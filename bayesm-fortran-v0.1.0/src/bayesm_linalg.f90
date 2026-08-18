module bayesm_linalg
  use bayesm_kinds, only: dp
  implicit none
  private
  public :: chol_upper, inverse_upper, solve_upper, solve_lower, solve_spd, inverse_spd
  public :: inverse_general, logdet_spd, outer_product, symmetrize, identity_matrix, trace_matrix
  public :: covariance_matrix, correlation_matrix, matrix_sqrt_spd
contains
  pure function identity_matrix(n) result(a)
    integer, intent(in) :: n
    real(dp) :: a(n,n)
    integer :: i
    a = 0.0_dp
    do i=1,n
      a(i,i)=1.0_dp
    end do
  end function identity_matrix

  pure subroutine chol_upper(a,r,info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: r(size(a,1),size(a,2))
    integer, intent(out) :: info
    integer :: n,i,j,k
    real(dp) :: s
    n=size(a,1)
    r=0.0_dp
    info=0
    if (size(a,2) /= n) then
      info=-1
      return
    end if
    do j=1,n
      s=a(j,j)
      do k=1,j-1
        s=s-r(k,j)*r(k,j)
      end do
      if (s <= 0.0_dp) then
        info=j
        return
      end if
      r(j,j)=sqrt(s)
      do i=j+1,n
        s=a(j,i)
        do k=1,j-1
          s=s-r(k,j)*r(k,i)
        end do
        r(j,i)=s/r(j,j)
      end do
    end do
  end subroutine chol_upper

  pure subroutine solve_upper(r,b,x,info)
    real(dp), intent(in) :: r(:,:), b(:)
    real(dp), intent(out) :: x(size(b))
    integer, intent(out) :: info
    integer :: n,i
    n=size(b)
    x=b
    info=0
    do i=n,1,-1
      if (abs(r(i,i)) <= tiny(1.0_dp)) then
        info=i
        return
      end if
      if (i<n) x(i)=x(i)-dot_product(r(i,i+1:n),x(i+1:n))
      x(i)=x(i)/r(i,i)
    end do
  end subroutine solve_upper

  pure subroutine solve_lower(l,b,x,info)
    real(dp), intent(in) :: l(:,:), b(:)
    real(dp), intent(out) :: x(size(b))
    integer, intent(out) :: info
    integer :: n,i
    n=size(b)
    x=b
    info=0
    do i=1,n
      if (abs(l(i,i)) <= tiny(1.0_dp)) then
        info=i
        return
      end if
      if (i>1) x(i)=x(i)-dot_product(l(i,1:i-1),x(1:i-1))
      x(i)=x(i)/l(i,i)
    end do
  end subroutine solve_lower

  pure subroutine inverse_upper(r,ri,info)
    real(dp), intent(in) :: r(:,:)
    real(dp), intent(out) :: ri(size(r,1),size(r,2))
    integer, intent(out) :: info
    integer :: n,j,istat
    real(dp) :: e(size(r,1)),x(size(r,1))
    n=size(r,1)
    ri=0.0_dp
    info=0
    do j=1,n
      e=0.0_dp
      e(j)=1.0_dp
      call solve_upper(r,e,x,istat)
      if (istat /= 0) then
        info=istat
        return
      end if
      ri(:,j)=x
    end do
  end subroutine inverse_upper

  pure subroutine solve_spd(a,b,x,info)
    real(dp), intent(in) :: a(:,:),b(:)
    real(dp), intent(out) :: x(size(b))
    integer, intent(out) :: info
    real(dp) :: r(size(a,1),size(a,2)),y(size(b))
    call chol_upper(a,r,info)
    if (info/=0) return
    call solve_lower(transpose(r),b,y,info)
    if (info/=0) return
    call solve_upper(r,y,x,info)
  end subroutine solve_spd

  pure subroutine inverse_spd(a,ai,info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: ai(size(a,1),size(a,2))
    integer, intent(out) :: info
    real(dp) :: r(size(a,1),size(a,2)),ri(size(a,1),size(a,2))
    call chol_upper(a,r,info)
    if (info/=0) return
    call inverse_upper(r,ri,info)
    if (info/=0) return
    ai=matmul(ri,transpose(ri))
    call symmetrize(ai)
  end subroutine inverse_spd

  pure subroutine inverse_general(a,ai,info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: ai(size(a,1),size(a,2))
    integer, intent(out) :: info
    integer :: n,i,j,p
    real(dp) :: aug(size(a,1),2*size(a,1)),tmp(2*size(a,1)),pivotv,f
    n=size(a,1)
    info=0
    if (size(a,2)/=n) then
      info=-1
      ai=0.0_dp
      return
    end if
    aug(:,1:n)=a
    aug(:,n+1:2*n)=identity_matrix(n)
    do i=1,n
      p=i
      do j=i+1,n
        if (abs(aug(j,i))>abs(aug(p,i))) p=j
      end do
      if (abs(aug(p,i))<=tiny(1.0_dp)) then
        info=i
        ai=0.0_dp
        return
      end if
      if (p/=i) then
        tmp=aug(i,:)
        aug(i,:)=aug(p,:)
        aug(p,:)=tmp
      end if
      pivotv=aug(i,i)
      aug(i,:)=aug(i,:)/pivotv
      do j=1,n
        if (j==i) cycle
        f=aug(j,i)
        aug(j,:)=aug(j,:)-f*aug(i,:)
      end do
    end do
    ai=aug(:,n+1:2*n)
  end subroutine inverse_general

  real(dp) function logdet_spd(a,info) result(v)
    real(dp), intent(in) :: a(:,:)
    integer, intent(out) :: info
    real(dp) :: r(size(a,1),size(a,2))
    integer :: i
    call chol_upper(a,r,info)
    if (info/=0) then
      v=-huge(1.0_dp)
      return
    end if
    v=0.0_dp
    do i=1,size(a,1)
      v=v+2.0_dp*log(r(i,i))
    end do
  end function logdet_spd

  pure function outer_product(x,y) result(a)
    real(dp), intent(in) :: x(:),y(:)
    real(dp) :: a(size(x),size(y))
    integer :: i
    do i=1,size(x)
      a(i,:)=x(i)*y
    end do
  end function outer_product

  pure subroutine symmetrize(a)
    real(dp), intent(inout) :: a(:,:)
    a=0.5_dp*(a+transpose(a))
  end subroutine symmetrize

  pure real(dp) function trace_matrix(a) result(t)
    real(dp), intent(in) :: a(:,:)
    integer :: i
    t=0.0_dp
    do i=1,min(size(a,1),size(a,2))
      t=t+a(i,i)
    end do
  end function trace_matrix

  pure function covariance_matrix(x) result(cov)
    real(dp), intent(in) :: x(:,:)
    real(dp) :: cov(size(x,2),size(x,2)),mu(size(x,2))
    real(dp) :: xc(size(x,1),size(x,2))
    integer :: n,j
    n=size(x,1)
    mu=sum(x,dim=1)/real(max(1,n),dp)
    xc=x
    do j=1,size(x,2)
      xc(:,j)=xc(:,j)-mu(j)
    end do
    if (n>1) then
      cov=matmul(transpose(xc),xc)/real(n-1,dp)
    else
      cov=0.0_dp
    end if
  end function covariance_matrix

  pure function correlation_matrix(cov) result(cor)
    real(dp), intent(in) :: cov(:,:)
    real(dp) :: cor(size(cov,1),size(cov,2)),s(size(cov,1))
    integer :: i,j
    do i=1,size(cov,1)
      s(i)=sqrt(max(0.0_dp,cov(i,i)))
    end do
    cor=0.0_dp
    do j=1,size(cov,2)
      do i=1,size(cov,1)
        if (s(i)>0.0_dp .and. s(j)>0.0_dp) cor(i,j)=cov(i,j)/(s(i)*s(j))
      end do
    end do
  end function correlation_matrix

  subroutine matrix_sqrt_spd(a,s,info)
    ! Symmetric square root by cyclic Jacobi eigen decomposition.
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: s(size(a,1),size(a,2))
    integer, intent(out) :: info
    integer :: n,p,q,i,j,iter
    real(dp) :: d(size(a,1),size(a,2)),v(size(a,1),size(a,2))
    real(dp) :: app,aqq,apq,tau,t,c,ss,tmp,maxoff
    n=size(a,1)
    d=a
    v=identity_matrix(n)
    info=0
    do iter=1,100*n*n
      maxoff=0.0_dp; p=1; q=min(2,n)
      do j=2,n
        do i=1,j-1
          if (abs(d(i,j))>maxoff) then
            maxoff=abs(d(i,j)); p=i; q=j
          end if
        end do
      end do
      if (maxoff <= 100.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(d)))) exit
      app=d(p,p); aqq=d(q,q); apq=d(p,q)
      tau=(aqq-app)/(2.0_dp*apq)
      if (tau>=0.0_dp) then
        t=1.0_dp/(tau+sqrt(1.0_dp+tau*tau))
      else
        t=-1.0_dp/(-tau+sqrt(1.0_dp+tau*tau))
      end if
      c=1.0_dp/sqrt(1.0_dp+t*t); ss=t*c
      do i=1,n
        if (i/=p .and. i/=q) then
          tmp=d(i,p); d(i,p)=c*tmp-ss*d(i,q); d(p,i)=d(i,p)
          d(i,q)=ss*tmp+c*d(i,q); d(q,i)=d(i,q)
        end if
      end do
      d(p,p)=c*c*app-2.0_dp*c*ss*apq+ss*ss*aqq
      d(q,q)=ss*ss*app+2.0_dp*c*ss*apq+c*c*aqq
      d(p,q)=0.0_dp; d(q,p)=0.0_dp
      do i=1,n
        tmp=v(i,p); v(i,p)=c*tmp-ss*v(i,q); v(i,q)=ss*tmp+c*v(i,q)
      end do
    end do
    if (iter>100*n*n) info=1
    s=0.0_dp
    do i=1,n
      if (d(i,i)<-1.0e-10_dp) then
        info=max(info,2)
      end if
      s=s+sqrt(max(0.0_dp,d(i,i)))*outer_product(v(:,i),v(:,i))
    end do
    call symmetrize(s)
  end subroutine matrix_sqrt_spd
end module bayesm_linalg
