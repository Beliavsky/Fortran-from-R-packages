module fastmatrix_linalg
  use fastmatrix_base, only: dp, eye, inverse_matrix, solve_linear, chol_lower
  implicit none
  private
  public :: lu_decomp, lu_solve, lu_inverse, ldl_decomp, jacobi_eigen, power_method, cg_solve, seidel_solve
  public :: sherman_morrison, rank1_update, cholupdate, sweep_operator, matrix_norm, scaled_condition
  public :: matrix_sqrt, matrix_polynomial, whitening, equilibrate
contains
  subroutine lu_decomp(a,l,u,piv,info)
    real(dp),intent(in)::a(:,:)
    real(dp),intent(out)::l(:,:),u(:,:)
    integer,intent(out)::piv(:)
    integer,intent(out),optional::info
    real(dp)::m(size(a,1),size(a,2)),row(size(a,2))
    integer::n,i,j,k,p,itmp
    n=size(a,1)
    m=a
    piv=[(i,i=1,n)]
    if(present(info)) info=0
    do k=1,n-1
      p=k
      do i=k+1,n
      if(abs(m(i,k))>abs(m(p,k))) p=i
      end do
      if(abs(m(p,k))<=epsilon(1.0_dp)) then
      if(present(info)) info=k
      exit
      end if
      if(p/=k) then
      row=m(k,:)
      m(k,:)=m(p,:)
      m(p,:)=row
      itmp=piv(k)
      piv(k)=piv(p)
      piv(p)=itmp
      end if
      do i=k+1,n
      m(i,k)=m(i,k)/m(k,k)
      m(i,k+1:n)=m(i,k+1:n)-m(i,k)*m(k,k+1:n)
      end do
    end do
    l=0.0_dp
    u=0.0_dp
    do i=1,n
    l(i,i)=1.0_dp
    do j=1,i-1
    l(i,j)=m(i,j)
    end do
    do j=i,n
    u(i,j)=m(i,j)
    end do
    end do
  end subroutine

  subroutine lu_solve(l,u,piv,b,x)
    real(dp),intent(in)::l(:,:),u(:,:),b(:)
    integer,intent(in)::piv(:)
    real(dp),intent(out)::x(:)
    real(dp)::y(size(b)),pb(size(b))
    integer::n,i
    n=size(b)
    do i=1,n
    pb(i)=b(piv(i))
    end do
    y(1)=pb(1)
    do i=2,n
    y(i)=pb(i)-dot_product(l(i,1:i-1),y(1:i-1))
    end do
    x(n)=y(n)/u(n,n)
    do i=n-1,1,-1
    x(i)=(y(i)-dot_product(u(i,i+1:n),x(i+1:n)))/u(i,i)
    end do
  end subroutine

  subroutine lu_inverse(a,ainv,info)
    real(dp),intent(in)::a(:,:)
    real(dp),intent(out)::ainv(:,:)
    integer,intent(out),optional::info
    call inverse_matrix(a,ainv,info)
  end subroutine

  subroutine ldl_decomp(a,l,d,info)
    real(dp),intent(in)::a(:,:)
    real(dp),intent(out)::l(:,:),d(:)
    integer,intent(out),optional::info
    integer::n,i,j,k
    real(dp)::s
    n=size(a,1)
    l=0.0_dp
    d=0.0_dp
    if(present(info))info=0
    do i=1,n
      do j=1,i-1
        s=a(i,j)
        do k=1,j-1
        s=s-l(i,k)*d(k)*l(j,k)
        end do
        l(i,j)=s/d(j)
      end do
      s=a(i,i)
      do k=1,i-1
      s=s-l(i,k)*l(i,k)*d(k)
      end do
      d(i)=s
      l(i,i)=1.0_dp
      if(abs(d(i))<=epsilon(1.0_dp)) then
      if(present(info))info=i
      return
      end if
    end do
  end subroutine

  subroutine jacobi_eigen(a,eval,evec,tol,maxiter,info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: eval(:), evec(:,:)
    real(dp), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    integer, intent(out), optional :: info
    real(dp) :: b(size(a,1),size(a,2)), jmat(size(a,1),size(a,2))
    real(dp) :: tmp(size(a,1),size(a,2)), eps, mx, phi, c, sn
    integer :: n,p,q,i,k,it,mi
    n=size(a,1)
    b=a
    evec=eye(n)
    eps=1.0e-12_dp
    if(present(tol)) eps=tol
    mi=100*n*n
    if(present(maxiter)) mi=maxiter
    do it=1,mi
      mx=0.0_dp
      p=1
      q=min(2,n)
      do i=1,n-1
        do k=i+1,n
          if(abs(b(i,k))>mx) then
            mx=abs(b(i,k))
            p=i
            q=k
          end if
        end do
      end do
      if(mx<eps) exit
      phi=0.5_dp*atan2(2.0_dp*b(p,q), b(q,q)-b(p,p))
      c=cos(phi)
      sn=sin(phi)
      jmat=eye(n)
      jmat(p,p)=c
      jmat(q,q)=c
      jmat(p,q)=sn
      jmat(q,p)=-sn
      tmp=matmul(b,jmat)
      b=matmul(transpose(jmat),tmp)
      evec=matmul(evec,jmat)
    end do
    do i=1,n
      eval(i)=b(i,i)
    end do
    if(present(info)) then
      if(it<=mi) then
        info=0
      else
        info=1
      end if
    end if
  end subroutine

  subroutine power_method(a,value,vector,tol,maxiter,iter)
    real(dp),intent(in)::a(:,:)
    real(dp),intent(out)::value,vector(:)
    real(dp),intent(in),optional::tol
    integer,intent(in),optional::maxiter
    integer,intent(out),optional::iter
    real(dp)::v(size(vector)),w(size(vector)),old,eps
    integer::k,mi
    eps=1e-10_dp
    if(present(tol))eps=tol
    mi=1000
    if(present(maxiter))mi=maxiter
    v=1.0_dp/sqrt(real(size(v),dp))
    old=0
    do k=1,mi
    w=matmul(a,v)
    if(sqrt(sum(w*w))<=tiny(1.0_dp))exit
    v=w/sqrt(sum(w*w))
    value=dot_product(v,matmul(a,v))
    if(abs(value-old)<eps*max(1.0_dp,abs(value)))exit
    old=value
    end do
    vector=v
    if(present(iter))iter=k
  end subroutine

  subroutine cg_solve(a,b,x,tol,maxiter,iter)
    real(dp),intent(in)::a(:,:),b(:)
    real(dp),intent(out)::x(:)
    real(dp),intent(in),optional::tol
    integer,intent(in),optional::maxiter
    integer,intent(out),optional::iter
    real(dp)::r(size(b)),p(size(b)),ap(size(b)),rs,rsn,eps
    integer::k,mi
    eps=1e-10_dp
    if(present(tol))eps=tol
    mi=10*size(b)
    if(present(maxiter))mi=maxiter
    x=0
    r=b
    p=r
    rs=dot_product(r,r)
    do k=1,mi
    ap=matmul(a,p)
    if(abs(dot_product(p,ap))<=tiny(1.0_dp))exit
    x=x+(rs/dot_product(p,ap))*p
    r=r-(rs/dot_product(p,ap))*ap
    rsn=dot_product(r,r)
    if(sqrt(rsn)<eps)exit
    p=r+(rsn/rs)*p
    rs=rsn
    end do
    if(present(iter))iter=k
  end subroutine

  subroutine seidel_solve(a,b,x,tol,maxiter,iter)
    real(dp),intent(in)::a(:,:),b(:)
    real(dp),intent(out)::x(:)
    real(dp),intent(in),optional::tol
    integer,intent(in),optional::maxiter
    integer,intent(out),optional::iter
    real(dp)::old(size(b)),eps
    integer::i,k,n,mi
    n=size(b)
    eps=1e-10_dp
    if(present(tol))eps=tol
    mi=1000
    if(present(maxiter))mi=maxiter
    x=0
    do k=1,mi
    old=x
    do i=1,n
    x(i)=(b(i)-dot_product(a(i,1:i-1),x(1:i-1))-dot_product(a(i,i+1:n),old(i+1:n)))/a(i,i)
    end do
    if(maxval(abs(x-old))<eps)exit
    end do
    if(present(iter))iter=k
  end subroutine

  pure function sherman_morrison(a_inv,u,v) result(b)
    real(dp),intent(in)::a_inv(:,:),u(:),v(:)
    real(dp)::b(size(a_inv,1),size(a_inv,2)),au(size(u)),va(size(v)),den
    au=matmul(a_inv,u)
    va=matmul(v,a_inv)
    den=1.0_dp+dot_product(v,au)
    b=a_inv-spread(au,2,size(va))*spread(va,1,size(au))/den
  end function
  pure function rank1_update(a,alpha,x,y) result(b)
    real(dp),intent(in)::a(:,:),alpha,x(:),y(:)
    real(dp)::b(size(a,1),size(a,2))
    integer::i
    b=a
    do i=1,size(x)
    b(i,:)=b(i,:)+alpha*x(i)*y
    end do
  end function

  subroutine cholupdate(l,x,sign_update,info)
    real(dp),intent(inout)::l(:,:)
    real(dp),intent(in)::x(:)
    logical,intent(in),optional::sign_update
    integer,intent(out),optional::info
    real(dp)::a(size(l,1),size(l,1)),xx(size(x),size(x)),newl(size(l,1),size(l,1))
    logical::up
    integer::i
    up=.true.
    if(present(sign_update))up=sign_update
    a=matmul(l,transpose(l))
    do i=1,size(x)
    xx(i,:)=x(i)*x
    end do
    if(up)then
    a=a+xx
    else
    a=a-xx
    end if
    call chol_lower(a,newl,info)
    l=newl
  end subroutine

  subroutine sweep_operator(a,k)
    real(dp),intent(inout)::a(:,:)
    integer,intent(in)::k
    real(dp)::akk
    integer::i,j,n
    n=size(a,1)
    akk=a(k,k)
    do i=1,n
    if(i/=k)then
    do j=1,n
    if(j/=k)a(i,j)=a(i,j)-a(i,k)*a(k,j)/akk
    end do
    end if
    end do
    do i=1,n
    if(i/=k)then
    a(i,k)=a(i,k)/akk
    a(k,i)=a(k,i)/akk
    end if
    end do
    a(k,k)=-1.0_dp/akk
  end subroutine

  pure function matrix_norm(a,kind) result(v)
    real(dp),intent(in)::a(:,:)
    character(len=*),intent(in),optional::kind
    real(dp)::v
    character(len=8)::k
    k='fro'
    if(present(kind))k=kind
    select case(trim(k))
    case('1','one')
    v=maxval(sum(abs(a),dim=1))
    case('inf')
    v=maxval(sum(abs(a),dim=2))
    case default
    v=sqrt(sum(a*a))
    end select
  end function

  function scaled_condition(a) result(c)
    real(dp),intent(in)::a(:,:)
    real(dp)::c,ainv(size(a,1),size(a,2))
    integer::info
    call inverse_matrix(a,ainv,info)
    if(info/=0)then
    c=huge(1.0_dp)
    else
    c=matrix_norm(a,'1')*matrix_norm(ainv,'1')
    end if
  end function

  subroutine matrix_sqrt(a,s,info)
    real(dp),intent(in)::a(:,:)
    real(dp),intent(out)::s(:,:)
    integer,intent(out),optional::info
    real(dp)::eval(size(a,1)),evec(size(a,1),size(a,1))
    integer::i,ier
    call jacobi_eigen(a,eval,evec,info=ier)
    s=0
    do i=1,size(eval)
    if(eval(i)<-1e-10_dp)then
    if(present(info))info=i
    return
    end if
    s=s+sqrt(max(0.0_dp,eval(i)))*spread(evec(:,i),2,size(eval))*spread(evec(:,i),1,size(eval))
    end do
    if(present(info))info=0
  end subroutine

  pure function matrix_polynomial(a,coef) result(p)
    real(dp),intent(in)::a(:,:),coef(:)
    real(dp)::p(size(a,1),size(a,2))
    integer::k,n
    n=size(a,1)
    p=0
    do k=size(coef),1,-1
    p=matmul(p,a)
    p=p+coef(k)*eye(n)
    end do
  end function

  subroutine whitening(sigma,w,info)
    real(dp),intent(in)::sigma(:,:)
    real(dp),intent(out)::w(:,:)
    integer,intent(out),optional::info
    real(dp)::eval(size(sigma,1)),evec(size(sigma,1),size(sigma,1))
    integer::i,ier
    call jacobi_eigen(sigma,eval,evec,info=ier)
    w=0
    do i=1,size(eval)
    if(eval(i)<=0)then
    if(present(info))info=i
    return
    end if
    w=w+(1/sqrt(eval(i)))*spread(evec(:,i),2,size(eval))*spread(evec(:,i),1,size(eval))
    end do
    if(present(info))info=0
  end subroutine

  subroutine equilibrate(a,r,c,b)
    real(dp),intent(in)::a(:,:)
    real(dp),intent(out)::r(:),c(:),b(:,:)
    integer::i,j
    do i=1,size(a,1)
    r(i)=1.0_dp/max(maxval(abs(a(i,:))),tiny(1.0_dp))
    end do
    do j=1,size(a,2)
    c(j)=1.0_dp/max(maxval(abs(a(:,j))*r),tiny(1.0_dp))
    end do
    do i=1,size(a,1)
    do j=1,size(a,2)
    b(i,j)=r(i)*a(i,j)*c(j)
    end do
    end do
  end subroutine
end module
