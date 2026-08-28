module gpa_linalg
  use gpa_kinds, only: dp
  implicit none
  private
  public :: eye, inverse_matrix, solve_linear, logabsdet, cholesky_lower
  public :: symmetric_eigen_jacobi, polar_orthogonal, random_orthogonal
  public :: frobenius_norm, max_abs_matrix, normalize_columns

contains

  pure function eye(n) result(a)
    integer, intent(in) :: n
    real(dp) :: a(n,n)
    integer :: i
    a = 0.0_dp
    do i = 1, n
      a(i,i) = 1.0_dp
    end do
  end function eye

  pure function frobenius_norm(a) result(v)
    real(dp), intent(in) :: a(:,:)
    real(dp) :: v
    v = sqrt(sum(a*a))
  end function frobenius_norm

  pure function max_abs_matrix(a) result(v)
    real(dp), intent(in) :: a(:,:)
    real(dp) :: v
    v = maxval(abs(a))
  end function max_abs_matrix

  subroutine inverse_matrix(a, ainv, info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: ainv(size(a,1),size(a,2))
    integer, intent(out) :: info
    real(dp), allocatable :: aug(:,:), tmp(:)
    real(dp) :: piv, factor, scale
    integer :: n, i, k, imax
    n = size(a,1)
    if (size(a,2) /= n) then
      info = -1
      ainv = 0.0_dp
      return
    end if
    allocate(aug(n,2*n), tmp(2*n))
    aug(:,1:n) = a
    aug(:,n+1:2*n) = eye(n)
    info = 0
    do k = 1, n
      imax = k
      piv = abs(aug(k,k))
      do i = k+1, n
        if (abs(aug(i,k)) > piv) then
          piv = abs(aug(i,k))
          imax = i
        end if
      end do
      scale = max(1.0_dp, maxval(abs(aug(:,1:n))))
      if (piv <= sqrt(epsilon(1.0_dp))*scale) then
        info = k
        ainv = 0.0_dp
        return
      end if
      if (imax /= k) then
        tmp = aug(k,:)
        aug(k,:) = aug(imax,:)
        aug(imax,:) = tmp
      end if
      aug(k,:) = aug(k,:) / aug(k,k)
      do i = 1, n
        if (i == k) cycle
        factor = aug(i,k)
        if (abs(factor) > tiny(1.0_dp)) aug(i,:) = aug(i,:) - factor*aug(k,:)
      end do
    end do
    ainv = aug(:,n+1:2*n)
  end subroutine inverse_matrix

  subroutine solve_linear(a, b, x, info)
    real(dp), intent(in) :: a(:,:), b(:,:)
    real(dp), intent(out) :: x(size(a,2),size(b,2))
    integer, intent(out) :: info
    real(dp) :: ainv(size(a,1),size(a,2))
    call inverse_matrix(a, ainv, info)
    if (info == 0) then
      x = matmul(ainv,b)
    else
      x = 0.0_dp
    end if
  end subroutine solve_linear

  function logabsdet(a, sign_det, info) result(v)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out), optional :: sign_det
    integer, intent(out), optional :: info
    real(dp) :: v
    real(dp), allocatable :: m(:,:), tmp(:)
    real(dp) :: piv, fac, sgn, scale
    integer :: n, i, k, imax, ierr
    n = size(a,1)
    ierr = 0
    if (size(a,2) /= n) then
      ierr = -1
      v = -huge(1.0_dp)
      sgn = 0.0_dp
      if (present(sign_det)) sign_det = sgn
      if (present(info)) info = ierr
      return
    end if
    allocate(m(n,n),tmp(n))
    m = a
    v = 0.0_dp
    sgn = 1.0_dp
    do k=1,n
      imax=k
      piv=abs(m(k,k))
      do i=k+1,n
        if (abs(m(i,k))>piv) then
          piv=abs(m(i,k))
          imax=i
        end if
      end do
      scale=max(1.0_dp,maxval(abs(m)))
      if (piv <= sqrt(epsilon(1.0_dp))*scale) then
        ierr=k
        v=-huge(1.0_dp)
        sgn=0.0_dp
        exit
      end if
      if (imax/=k) then
        tmp=m(k,:)
        m(k,:)=m(imax,:)
        m(imax,:)=tmp
        sgn=-sgn
      end if
      if (m(k,k)<0.0_dp) sgn=-sgn
      v=v+log(abs(m(k,k)))
      do i=k+1,n
        fac=m(i,k)/m(k,k)
        m(i,k:n)=m(i,k:n)-fac*m(k,k:n)
      end do
    end do
    if (present(sign_det)) sign_det=sgn
    if (present(info)) info=ierr
  end function logabsdet

  subroutine cholesky_lower(a, l, info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: l(size(a,1),size(a,2))
    integer, intent(out) :: info
    real(dp) :: s
    integer :: n, i, j, k
    n=size(a,1)
    l=0.0_dp
    info=0
    do i=1,n
      do j=1,i
        s=a(i,j)
        do k=1,j-1
          s=s-l(i,k)*l(j,k)
        end do
        if (i==j) then
          if (s <= 0.0_dp) then
            info=i
            return
          end if
          l(i,j)=sqrt(s)
        else
          l(i,j)=s/l(j,j)
        end if
      end do
    end do
  end subroutine cholesky_lower

  subroutine symmetric_eigen_jacobi(a, eval, evec, info, tol, max_sweeps)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: eval(size(a,1)), evec(size(a,1),size(a,1))
    integer, intent(out) :: info
    real(dp), intent(in), optional :: tol
    integer, intent(in), optional :: max_sweeps
    real(dp) :: m(size(a,1),size(a,1)), c, s, tau, t, app, aqq, apq
    real(dp) :: thresh, eps0, mpj, mqj, vjp, vjq
    integer :: n, p, q, j, sweep, ms
    n=size(a,1)
    m=0.5_dp*(a+transpose(a))
    evec=eye(n)
    eps0=1.0e-13_dp
    if(present(tol)) eps0=tol
    ms=100*n*n
    if(present(max_sweeps)) ms=max_sweeps
    info=1
    do sweep=1,ms
      thresh=0.0_dp
      p=1
      q=min(2,n)
      do j=1,n
        if (j<n) then
          if (maxval(abs(m(j,j+1:n)))>thresh) then
            q=j+maxloc(abs(m(j,j+1:n)),dim=1)
            p=j
            thresh=abs(m(p,q))
          end if
        end if
      end do
      if (thresh <= eps0*max(1.0_dp,maxval(abs(m)))) then
        info=0
        exit
      end if
      apq=m(p,q)
      app=m(p,p)
      aqq=m(q,q)
      tau=(aqq-app)/(2.0_dp*apq)
      if (tau>=0.0_dp) then
        t=1.0_dp/(tau+sqrt(1.0_dp+tau*tau))
      else
        t=-1.0_dp/(-tau+sqrt(1.0_dp+tau*tau))
      end if
      c=1.0_dp/sqrt(1.0_dp+t*t)
      s=t*c
      do j=1,n
        if (j/=p .and. j/=q) then
          mpj=m(p,j)
          mqj=m(q,j)
          m(p,j)=c*mpj-s*mqj
          m(j,p)=m(p,j)
          m(q,j)=s*mpj+c*mqj
          m(j,q)=m(q,j)
        end if
      end do
      m(p,p)=c*c*app-2.0_dp*c*s*apq+s*s*aqq
      m(q,q)=s*s*app+2.0_dp*c*s*apq+c*c*aqq
      m(p,q)=0.0_dp
      m(q,p)=0.0_dp
      do j=1,n
        vjp=evec(j,p)
        vjq=evec(j,q)
        evec(j,p)=c*vjp-s*vjq
        evec(j,q)=s*vjp+c*vjq
      end do
    end do
    do j=1,n
      eval(j)=m(j,j)
    end do
  end subroutine symmetric_eigen_jacobi

  subroutine polar_orthogonal(x, q, info)
    real(dp), intent(in) :: x(:,:)
    real(dp), intent(out) :: q(size(x,1),size(x,2))
    integer, intent(out) :: info
    real(dp) :: c(size(x,2),size(x,2)), eval(size(x,2))
    real(dp) :: v(size(x,2),size(x,2)), invsqrt(size(x,2),size(x,2))
    integer :: j
    c=matmul(transpose(x),x)
    call symmetric_eigen_jacobi(c,eval,v,info)
    if (info/=0 .or. minval(eval)<=sqrt(epsilon(1.0_dp))) then
      q=0.0_dp
      info=max(1,info)
      return
    end if
    invsqrt=0.0_dp
    do j=1,size(eval)
      invsqrt=invsqrt+(1.0_dp/sqrt(eval(j)))*outer_col(v(:,j),v(:,j))
    end do
    q=matmul(x,invsqrt)
  contains
    pure function outer_col(a,b) result(m)
      real(dp), intent(in) :: a(:),b(:)
      real(dp) :: m(size(a),size(b))
      integer :: ii,jj
      do jj=1,size(b)
      do ii=1,size(a)
      m(ii,jj)=a(ii)*b(jj)
      end do
      end do
    end function outer_col
  end subroutine polar_orthogonal

  subroutine normalize_columns(x, q, info)
    real(dp), intent(in) :: x(:,:)
    real(dp), intent(out) :: q(size(x,1),size(x,2))
    integer, intent(out) :: info
    real(dp) :: nrm
    integer :: j
    q=x
    info=0
    do j=1,size(x,2)
      nrm=sqrt(sum(x(:,j)*x(:,j)))
      if (nrm<=sqrt(epsilon(1.0_dp))) then
        info=j
        q=0.0_dp
        return
      end if
      q(:,j)=x(:,j)/nrm
    end do
  end subroutine normalize_columns

  subroutine random_orthogonal(k, q, info)
    integer, intent(in) :: k
    real(dp), intent(out) :: q(k,k)
    integer, intent(out) :: info
    real(dp) :: z(k,k), rdiag(k), u1, u2, nrm, proj
    integer :: i,j
    real(dp), parameter :: twopi=6.283185307179586476925286766559_dp
    do j=1,k
      i=1
      do while(i<=k)
        call random_number(u1)
        call random_number(u2)
        u1=max(u1,tiny(1.0_dp))
        z(i,j)=sqrt(-2.0_dp*log(u1))*cos(twopi*u2)
        if(i+1<=k) z(i+1,j)=sqrt(-2.0_dp*log(u1))*sin(twopi*u2)
        i=i+2
      end do
    end do
    q=0.0_dp
    info=0
    do j=1,k
      q(:,j)=z(:,j)
      do i=1,j-1
        proj=dot_product(q(:,i),q(:,j))
        q(:,j)=q(:,j)-proj*q(:,i)
      end do
      nrm=sqrt(sum(q(:,j)*q(:,j)))
      rdiag(j)=nrm
      if(nrm<=sqrt(epsilon(1.0_dp))) then
      info=j
      return
      end if
      q(:,j)=q(:,j)/nrm
      if(rdiag(j)<0.0_dp) q(:,j)=-q(:,j)
    end do
  end subroutine random_orthogonal

end module gpa_linalg
