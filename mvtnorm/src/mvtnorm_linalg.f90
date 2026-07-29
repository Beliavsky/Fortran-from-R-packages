! SPDX-License-Identifier: GPL-2.0-only
module mvtnorm_linalg
  use mvtnorm_kinds, only : dp
  implicit none
  private
  public :: cholesky_lower, solve_lower, solve_upper, solve_spd
  public :: inverse_lower, inverse_spd, covariance_to_correlation
  public :: symmetrize, logdet_cholesky, jacobi_eigen, nearest_psd
  public :: matrix_rank, identity_matrix

contains

  subroutine cholesky_lower(a, l, ok, message, tolerance)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: l(:,:)
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message
    real(dp), intent(in), optional :: tolerance
    integer :: n, i, j
    real(dp) :: s, tol

    n = size(a,1)
    allocate(l(n,n)); l = 0.0_dp
    ok = .false.; message = ''
    if (size(a,2) /= n) then
      message = 'matrix must be square'
      return
    end if
    tol = 100.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(a)))
    if (present(tolerance)) tol = tolerance
    do i = 1, n
      do j = 1, i
        s = a(i,j)
        if (j > 1) s = s-dot_product(l(i,1:j-1),l(j,1:j-1))
        if (i == j) then
          if (s <= tol) then
            message = 'matrix is not positive definite'
            return
          end if
          l(i,j) = sqrt(s)
        else
          l(i,j) = s/l(j,j)
        end if
      end do
    end do
    ok = .true.
  end subroutine cholesky_lower

  subroutine solve_lower(l, b, x, ok, unit_diagonal, transpose)
    real(dp), intent(in) :: l(:,:), b(:,:)
    real(dp), allocatable, intent(out) :: x(:,:)
    logical, intent(out) :: ok
    logical, intent(in), optional :: unit_diagonal, transpose
    logical :: unitd, trans
    integer :: n, m, i, j
    real(dp) :: d

    n = size(l,1); m = size(b,2)
    allocate(x(n,m)); x = b
    ok = .false.
    if (size(l,2) /= n .or. size(b,1) /= n) return
    unitd = .false.; if (present(unit_diagonal)) unitd = unit_diagonal
    trans = .false.; if (present(transpose)) trans = transpose
    if (.not. trans) then
      do i = 1, n
        if (i > 1) x(i,:) = x(i,:)-matmul(l(i,1:i-1),x(1:i-1,:))
        d = merge(1.0_dp,l(i,i),unitd)
        if (abs(d) <= tiny(1.0_dp)) return
        x(i,:) = x(i,:)/d
      end do
    else
      do i = n, 1, -1
        if (i < n) then
          do j = 1, m
            x(i,j) = x(i,j)-dot_product(l(i+1:n,i),x(i+1:n,j))
          end do
        end if
        d = merge(1.0_dp,l(i,i),unitd)
        if (abs(d) <= tiny(1.0_dp)) return
        x(i,:) = x(i,:)/d
      end do
    end if
    ok = .true.
  end subroutine solve_lower

  subroutine solve_upper(u, b, x, ok)
    real(dp), intent(in) :: u(:,:), b(:,:)
    real(dp), allocatable, intent(out) :: x(:,:)
    logical, intent(out) :: ok
    real(dp), allocatable :: lt(:,:)
    lt = transpose(u)
    call solve_lower(lt,b,x,ok,transpose=.true.)
  end subroutine solve_upper

  subroutine solve_spd(a,b,x,ok,message)
    real(dp), intent(in) :: a(:,:), b(:,:)
    real(dp), allocatable, intent(out) :: x(:,:)
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message
    real(dp), allocatable :: l(:,:), y(:,:)
    logical :: ok1
    call cholesky_lower(a,l,ok,message)
    if (.not. ok) then
      allocate(x(size(b,1),size(b,2))); x=0.0_dp
      return
    end if
    call solve_lower(l,b,y,ok1)
    if (.not. ok1) then
      ok=.false.; message='lower triangular solve failed'
      allocate(x(size(b,1),size(b,2))); x=0.0_dp
      return
    end if
    call solve_lower(l,y,x,ok1,transpose=.true.)
    ok=ok1
    if (.not. ok) message='upper triangular solve failed'
  end subroutine solve_spd

  subroutine inverse_lower(l, linv, ok)
    real(dp), intent(in) :: l(:,:)
    real(dp), allocatable, intent(out) :: linv(:,:)
    logical, intent(out) :: ok
    real(dp), allocatable :: eye(:,:)
    integer :: n
    n=size(l,1)
    eye=identity_matrix(n)
    call solve_lower(l,eye,linv,ok)
  end subroutine inverse_lower

  subroutine inverse_spd(a, ainv, ok, message)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: ainv(:,:)
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message
    real(dp), allocatable :: eye(:,:)
    eye=identity_matrix(size(a,1))
    call solve_spd(a,eye,ainv,ok,message)
    if (ok) call symmetrize(ainv)
  end subroutine inverse_spd

  function identity_matrix(n) result(a)
    integer, intent(in) :: n
    real(dp), allocatable :: a(:,:)
    integer :: i
    allocate(a(n,n)); a=0.0_dp
    do i=1,n
      a(i,i)=1.0_dp
    end do
  end function identity_matrix

  subroutine covariance_to_correlation(cov, cor, sd, ok, message)
    real(dp), intent(in) :: cov(:,:)
    real(dp), allocatable, intent(out) :: cor(:,:), sd(:)
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message
    integer :: n,i,j
    n=size(cov,1)
    allocate(cor(n,n),sd(n)); cor=0.0_dp; sd=0.0_dp
    ok=.false.; message=''
    if (size(cov,2)/=n) then
      message='covariance matrix must be square'; return
    end if
    do i=1,n
      if (cov(i,i)<=0.0_dp) then
        message='covariance diagonal must be positive'; return
      end if
      sd(i)=sqrt(cov(i,i))
    end do
    do i=1,n
      do j=1,n
        cor(i,j)=cov(i,j)/(sd(i)*sd(j))
      end do
      cor(i,i)=1.0_dp
    end do
    call symmetrize(cor)
    ok=.true.
  end subroutine covariance_to_correlation

  subroutine symmetrize(a)
    real(dp), intent(inout) :: a(:,:)
    integer :: i,j,n
    real(dp) :: s
    n=min(size(a,1),size(a,2))
    do i=1,n
      do j=1,i-1
        s=0.5_dp*(a(i,j)+a(j,i))
        a(i,j)=s; a(j,i)=s
      end do
    end do
  end subroutine symmetrize

  real(dp) function logdet_cholesky(l) result(v)
    real(dp), intent(in) :: l(:,:)
    integer :: i
    v=0.0_dp
    do i=1,min(size(l,1),size(l,2))
      if (l(i,i)<=0.0_dp) then
        v=-huge(1.0_dp); return
      end if
      v=v+2.0_dp*log(l(i,i))
    end do
  end function logdet_cholesky

  subroutine jacobi_eigen(a, values, vectors, ok)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: values(:), vectors(:,:)
    logical, intent(out) :: ok
    real(dp), allocatable :: b(:,:)
    real(dp) :: app,aqq,apq,tau,t,c,s,bip,biq,vip,viq,maxoff
    integer :: n,p,q,i,iter,maxiter,k
    n=size(a,1)
    allocate(b(n,n),values(n),vectors(n,n))
    b=a; call symmetrize(b); vectors=identity_matrix(n)
    maxiter=max(50,100*n*n); ok=.false.
    do iter=1,maxiter
      maxoff=0.0_dp; p=1; q=min(2,n)
      do i=1,n
        do k=i+1,n
          if (abs(b(i,k))>maxoff) then
            maxoff=abs(b(i,k)); p=i; q=k
          end if
        end do
      end do
      if (maxoff <= 100.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(b)))) then
        ok=.true.; exit
      end if
      app=b(p,p); aqq=b(q,q); apq=b(p,q)
      tau=(aqq-app)/(2.0_dp*apq)
      if (tau>=0.0_dp) then
        t=1.0_dp/(tau+sqrt(1.0_dp+tau*tau))
      else
        t=-1.0_dp/(-tau+sqrt(1.0_dp+tau*tau))
      end if
      c=1.0_dp/sqrt(1.0_dp+t*t); s=t*c
      do i=1,n
        if (i/=p .and. i/=q) then
          bip=b(i,p); biq=b(i,q)
          b(i,p)=c*bip-s*biq; b(p,i)=b(i,p)
          b(i,q)=s*bip+c*biq; b(q,i)=b(i,q)
        end if
        vip=vectors(i,p); viq=vectors(i,q)
        vectors(i,p)=c*vip-s*viq
        vectors(i,q)=s*vip+c*viq
      end do
      b(p,p)=c*c*app-2.0_dp*s*c*apq+s*s*aqq
      b(q,q)=s*s*app+2.0_dp*s*c*apq+c*c*aqq
      b(p,q)=0.0_dp; b(q,p)=0.0_dp
    end do
    do i=1,n
      values(i)=b(i,i)
    end do
    call sort_eigen(values,vectors)
  end subroutine jacobi_eigen

  subroutine sort_eigen(values,vectors)
    real(dp), intent(inout) :: values(:),vectors(:,:)
    integer :: i,j,k,n
    real(dp) :: tv
    real(dp), allocatable :: col(:)
    n=size(values); allocate(col(size(vectors,1)))
    do i=1,n-1
      k=i
      do j=i+1,n
        if (values(j)>values(k)) k=j
      end do
      if (k/=i) then
        tv=values(i); values(i)=values(k); values(k)=tv
        col=vectors(:,i); vectors(:,i)=vectors(:,k); vectors(:,k)=col
      end if
    end do
  end subroutine sort_eigen

  subroutine nearest_psd(a, out, tolerance)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: out(:,:)
    real(dp), intent(in), optional :: tolerance
    real(dp), allocatable :: values(:),vectors(:,:)
    real(dp) :: tol
    logical :: ok
    integer :: i,n
    n=size(a,1); tol=100.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(a)))
    if (present(tolerance)) tol=tolerance
    call jacobi_eigen(a,values,vectors,ok)
    allocate(out(n,n)); out=0.0_dp
    do i=1,n
      values(i)=max(values(i),tol)
      out=out+values(i)*spread(vectors(:,i),2,n)*spread(vectors(:,i),1,n)
    end do
    call symmetrize(out)
  end subroutine nearest_psd

  integer function matrix_rank(a,tolerance) result(r)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(in), optional :: tolerance
    real(dp), allocatable :: ata(:,:),values(:),vectors(:,:)
    real(dp) :: tol
    logical :: ok
    ata=matmul(transpose(a),a)
    call jacobi_eigen(ata,values,vectors,ok)
    tol=100.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(values))
    if (present(tolerance)) tol=tolerance
    r=count(values>tol)
  end function matrix_rank

end module mvtnorm_linalg
