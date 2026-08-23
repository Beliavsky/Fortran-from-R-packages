module fastmatrix_products
  use fastmatrix_base, only: dp, outer_product
  implicit none
  private
  public :: vec, vech, hadamard_prod, kronecker_prod, matrix_inner, bracket_prod, symm_prod, house, house_prod
contains
  pure function vec(a) result(v)
    real(dp),intent(in)::a(:,:)
    real(dp)::v(size(a))
    v=reshape(a,[size(a)])
  end function
  pure function vech(a) result(v)
    real(dp),intent(in)::a(:,:)
    real(dp)::v(size(a,1)*(size(a,1)+1)/2)
    integer::i,j,k,n
    n=size(a,1)
    k=0
    do j=1,n
    do i=j,n
    k=k+1
    v(k)=a(i,j)
    end do
    end do
  end function
  pure function hadamard_prod(a,b) result(c)
    real(dp),intent(in)::a(:,:),b(:,:)
    real(dp)::c(size(a,1),size(a,2))
    c=a*b
  end function
  pure function kronecker_prod(a,b) result(c)
    real(dp),intent(in)::a(:,:),b(:,:)
    real(dp)::c(size(a,1)*size(b,1),size(a,2)*size(b,2))
    integer::i,j,m,n
    m=size(b,1)
    n=size(b,2)
    do j=1,size(a,2)
    do i=1,size(a,1)
    c((i-1)*m+1:i*m,(j-1)*n+1:j*n)=a(i,j)*b
    end do
    end do
  end function
  pure function matrix_inner(a,b) result(x)
    real(dp),intent(in)::a(:,:),b(:,:)
    real(dp)::x
    x=sum(a*b)
  end function
  pure function bracket_prod(a,b,c) result(d)
    real(dp),intent(in)::a(:,:),b(:,:),c(:,:)
    real(dp)::d(size(a,1),size(c,2))
    d=matmul(matmul(a,b),c)
  end function
  pure function symm_prod(a) result(c)
    real(dp),intent(in)::a(:,:)
    real(dp)::c(size(a,1),size(a,1))
    c=matmul(a,transpose(a))
  end function
  pure function house(x) result(h)
    real(dp),intent(in)::x(:)
    real(dp)::h(size(x),size(x)),u(size(x)),nu
    u=x
    nu=sqrt(sum(u*u))
    if(nu<=tiny(1.0_dp)) then
    h=0.0_dp
    return
    end if
    if(u(1)>=0.0_dp) then
    u(1)=u(1)+nu
    else
    u(1)=u(1)-nu
    end if
    nu=sqrt(sum(u*u))
    u=u/nu
    h=-2.0_dp*outer_product(u,u)
    h=h+diag_eye(size(x))
  contains
    pure function diag_eye(n) result(a)
    integer,intent(in)::n
    real(dp)::a(n,n)
    integer::i
    a=0
    do i=1,n
    a(i,i)=1
    end do
    end function
  end function
  pure function house_prod(h,a) result(c)
    real(dp),intent(in)::h(:,:),a(:,:)
    real(dp)::c(size(h,1),size(a,2))
    c=matmul(h,a)
  end function
end module
