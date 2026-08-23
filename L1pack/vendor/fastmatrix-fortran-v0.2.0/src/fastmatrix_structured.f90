module fastmatrix_structured
  use fastmatrix_base, only: dp, eye
  implicit none
  private
  public :: cor_ar1, cor_cs, circulant, hankel, frank, helmert, duplication, commutation, symmetrizer
contains
  pure function cor_ar1(rho,p) result(a)
    real(dp), intent(in)::rho
    integer,intent(in)::p
    real(dp)::a(p,p)
    integer::i,j
    do j=1,p
    do i=1,p
    a(i,j)=rho**abs(i-j)
    end do
    end do
  end function
  pure function cor_cs(rho,p) result(a)
    real(dp),intent(in)::rho
    integer,intent(in)::p
    real(dp)::a(p,p)
    integer::i
    a=rho
    do i=1,p
    a(i,i)=1.0_dp
    end do
  end function
  pure function circulant(x) result(a)
    real(dp),intent(in)::x(:)
    real(dp)::a(size(x),size(x))
    integer::n,i,j,k
    n=size(x)
    do j=1,n
    do i=1,n
    k=mod(i-j,n)+1
    a(i,j)=x(k)
    end do
    end do
  end function
  pure function hankel(x,y) result(a)
    real(dp),intent(in)::x(:),y(:)
    real(dp)::a(size(x),size(y))
    integer::i,j,k,n
    n=size(x)
    do j=1,size(y)
    do i=1,n
    k=i+j-1
    if(k<=n) then
    a(i,j)=x(k)
    else
    a(i,j)=y(k-n+1)
    end if
    end do
    end do
  end function
  pure function frank(n) result(a)
    integer,intent(in)::n
    real(dp)::a(n,n)
    integer::i
    a=0.0_dp
    do i=1,n
    a(i,i)=real(n-i+1,dp)
    if(i<n) then
    a(i,i+1)=real(n-i,dp)
    a(i+1,i)=real(n-i,dp)
    end if
    end do
  end function
  pure function helmert(n) result(h)
    integer,intent(in)::n
    real(dp)::h(n-1,n)
    integer::i
    h=0.0_dp
    do i=1,n-1
      h(i,1:i)=1.0_dp/sqrt(real(i*(i+1),dp))
      h(i,i+1)=-real(i,dp)/sqrt(real(i*(i+1),dp))
    end do
  end function
  pure function commutation(m,n) result(km)
    integer,intent(in)::m,n
    real(dp)::km(m*n,m*n)
    integer::i,j,r,c
    km=0.0_dp
    do j=1,n
    do i=1,m
    c=i+(j-1)*m
    r=j+(i-1)*n
    km(r,c)=1.0_dp
    end do
    end do
  end function
  pure function symmetrizer(n) result(s)
    integer,intent(in)::n
    real(dp)::s(n*n,n*n)
    s=0.5_dp*(eye(n*n)+commutation(n,n))
  end function
  pure function duplication(n) result(d)
    integer,intent(in)::n
    real(dp)::d(n*n,n*(n+1)/2)
    integer::i,j,k,r1,r2
    d=0.0_dp
    k=0
    do j=1,n
    do i=j,n
    k=k+1
    r1=i+(j-1)*n
    r2=j+(i-1)*n
    d(r1,k)=1.0_dp
    d(r2,k)=1.0_dp
    end do
    end do
  end function
end module
