module fastmatrix_base
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
  private
  integer, parameter, public :: dp = real64
  real(dp), parameter, public :: pi = acos(-1.0_dp)
  public :: eye, outer_product, solve_linear, inverse_matrix, chol_lower, normal_rand, gamma_p, normal_cdf
contains
  pure function eye(n) result(a)
    integer, intent(in) :: n
    real(dp) :: a(n,n)
    integer :: i
    a = 0.0_dp
    do i=1,n
    a(i,i)=1.0_dp
    end do
  end function

  pure function outer_product(x,y) result(a)
    real(dp), intent(in) :: x(:), y(:)
    real(dp) :: a(size(x),size(y))
    integer :: i
    do i=1,size(x)
    a(i,:) = x(i)*y
    end do
  end function

  subroutine solve_linear(a,b,x,info)
    real(dp), intent(in) :: a(:,:), b(:)
    real(dp), intent(out) :: x(:)
    integer, intent(out), optional :: info
    real(dp), allocatable :: m(:,:), rhs(:), row(:)
    real(dp) :: fac, piv
    integer :: n,i,k,p
    n=size(b)
    allocate(m(n,n),rhs(n),row(n))
    m=a
    rhs=b
    if (present(info)) info=0
    do k=1,n-1
      p=k
      piv=abs(m(k,k))
      do i=k+1,n
        if (abs(m(i,k))>piv) then
        p=i
        piv=abs(m(i,k))
        end if
      end do
      if (piv <= epsilon(1.0_dp)*max(1.0_dp,maxval(abs(m)))) then
        x=0.0_dp
        if(present(info)) info=k
        return
      end if
      if(p/=k) then
        row=m(k,:)
        m(k,:)=m(p,:)
        m(p,:)=row
        fac=rhs(k)
        rhs(k)=rhs(p)
        rhs(p)=fac
      end if
      do i=k+1,n
        fac=m(i,k)/m(k,k)
        m(i,k)=0.0_dp
        m(i,k+1:n)=m(i,k+1:n)-fac*m(k,k+1:n)
        rhs(i)=rhs(i)-fac*rhs(k)
      end do
    end do
    if(abs(m(n,n))<=epsilon(1.0_dp)) then
    x=0.0_dp
    if(present(info)) info=n
    return
    end if
    x(n)=rhs(n)/m(n,n)
    do i=n-1,1,-1
    x(i)=(rhs(i)-dot_product(m(i,i+1:n),x(i+1:n)))/m(i,i)
    end do
  end subroutine

  subroutine inverse_matrix(a,ainv,info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: ainv(:,:)
    integer, intent(out), optional :: info
    real(dp), allocatable :: e(:),x(:)
    integer :: n,j,ier
    n=size(a,1)
    allocate(e(n),x(n))
    ainv=0.0_dp
    ier=0
    do j=1,n
      e=0.0_dp
      e(j)=1.0_dp
      call solve_linear(a,e,x,ier)
      if(ier/=0) exit
      ainv(:,j)=x
    end do
    if(present(info)) info=ier
  end subroutine

  subroutine chol_lower(a,l,info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: l(:,:)
    integer, intent(out), optional :: info
    integer :: n,i,j,k
    real(dp) :: s
    n=size(a,1)
    l=0.0_dp
    if(present(info)) info=0
    do i=1,n
      do j=1,i
        s=a(i,j)
        do k=1,j-1
        s=s-l(i,k)*l(j,k)
        end do
        if(i==j) then
          if(s<=0.0_dp) then
          if(present(info)) info=i
          return
          end if
          l(i,j)=sqrt(s)
        else
          l(i,j)=s/l(j,j)
        end if
      end do
    end do
  end subroutine

  function normal_rand() result(z)
    real(dp) :: z,u1,u2
    call random_number(u1)
    call random_number(u2)
    u1=max(u1,tiny(1.0_dp))
    z=sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
  end function

  pure function normal_cdf(x) result(p)
    real(dp), intent(in) :: x
    real(dp) :: p
    p=0.5_dp*erfc(-x/sqrt(2.0_dp))
  end function

  pure function gamma_p(a,x) result(p)
    real(dp), intent(in) :: a,x
    real(dp) :: p,ap,del,sumv,b,c,d,h,an
    integer :: n
    if(x<=0.0_dp) then
    p=0.0_dp
    return
    end if
    if(x<a+1.0_dp) then
      ap=a
      del=1.0_dp/a
      sumv=del
      do n=1,500
        ap=ap+1.0_dp
        del=del*x/ap
        sumv=sumv+del
        if(abs(del)<abs(sumv)*1.0e-14_dp) exit
      end do
      p=sumv*exp(-x+a*log(x)-log_gamma(a))
    else
      b=x+1.0_dp-a
      c=1.0_dp/tiny(1.0_dp)
      d=1.0_dp/b
      h=d
      do n=1,500
        an=-real(n,dp)*(real(n,dp)-a)
        b=b+2.0_dp
        d=an*d+b
        if(abs(d)<tiny(1.0_dp)) d=tiny(1.0_dp)
        c=b+an/c
        if(abs(c)<tiny(1.0_dp)) c=tiny(1.0_dp)
        d=1.0_dp/d
        del=d*c
        h=h*del
        if(abs(del-1.0_dp)<1.0e-14_dp) exit
      end do
      p=1.0_dp-exp(-x+a*log(x)-log_gamma(a))*h
    end if
    p=max(0.0_dp,min(1.0_dp,p))
  end function
end module
