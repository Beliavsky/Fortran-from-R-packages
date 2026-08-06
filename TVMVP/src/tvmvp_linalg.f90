! SPDX-License-Identifier: MIT
! Modern Fortran translation of computational routines from TVMVP.
module tvmvp_linalg
  use tvmvp_kinds, only : dp
  use tvmvp_status, only : tvmvp_error, clear_error, set_error, tvmvp_invalid_input, &
                           tvmvp_singular, tvmvp_not_positive_definite, tvmvp_convergence_failure
  implicit none
  private
  public :: solve_linear, solve_matrix, invert_matrix, cholesky_factor
  public :: symmetric_eigen, floor_symmetric_eigenvalues, matrix_sqrt_abs
  public :: covariance_matrix, standardize_columns, identity_matrix, outer_product
  public :: sample_mean, sample_variance, correlation
contains
  pure function identity_matrix(n) result(a)
    integer, intent(in) :: n
    real(dp) :: a(n,n)
    integer :: i
    a = 0.0_dp
    do i = 1, n
      a(i,i) = 1.0_dp
    end do
  end function identity_matrix

  pure function outer_product(x, y) result(a)
    real(dp), intent(in) :: x(:), y(:)
    real(dp) :: a(size(x),size(y))
    integer :: i
    do i = 1, size(x)
      a(i,:) = x(i) * y
    end do
  end function outer_product

  pure real(dp) function sample_mean(x)
    real(dp), intent(in) :: x(:)
    if (size(x) > 0) then
      sample_mean = sum(x) / real(size(x),dp)
    else
      sample_mean = 0.0_dp
    end if
  end function sample_mean

  pure real(dp) function sample_variance(x)
    real(dp), intent(in) :: x(:)
    real(dp) :: mu
    if (size(x) > 1) then
      mu = sample_mean(x)
      sample_variance = sum((x-mu)**2) / real(size(x)-1,dp)
    else
      sample_variance = 0.0_dp
    end if
  end function sample_variance

  pure real(dp) function correlation(x, y)
    real(dp), intent(in) :: x(:), y(:)
    real(dp) :: mx, my, sx, sy
    if (size(x) /= size(y) .or. size(x) < 2) then
      correlation = 0.0_dp
      return
    end if
    mx = sample_mean(x)
    my = sample_mean(y)
    sx = sqrt(sum((x-mx)**2))
    sy = sqrt(sum((y-my)**2))
    if (sx <= tiny(1.0_dp) .or. sy <= tiny(1.0_dp)) then
      correlation = 0.0_dp
    else
      correlation = dot_product(x-mx,y-my)/(sx*sy)
    end if
  end function correlation

  subroutine standardize_columns(x, z, err)
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable, intent(out) :: z(:,:)
    type(tvmvp_error), intent(out) :: err
    real(dp) :: mu, sd
    integer :: j
    call clear_error(err)
    allocate(z(size(x,1),size(x,2)))
    do j = 1, size(x,2)
      mu = sample_mean(x(:,j))
      sd = sqrt(sample_variance(x(:,j)))
      if (sd <= sqrt(epsilon(1.0_dp))) then
        call set_error(err, tvmvp_invalid_input, 'cannot standardize a constant column')
        z = 0.0_dp
        return
      end if
      z(:,j) = (x(:,j)-mu)/sd
    end do
  end subroutine standardize_columns

  subroutine covariance_matrix(x, cov, unbiased, err)
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable, intent(out) :: cov(:,:)
    logical, intent(in), optional :: unbiased
    type(tvmvp_error), intent(out) :: err
    real(dp), allocatable :: centered(:,:)
    real(dp) :: denom
    integer :: j, n
    logical :: ub
    call clear_error(err)
    n = size(x,1)
    if (n < 1) then
      allocate(cov(0,0))
      call set_error(err, tvmvp_invalid_input, 'covariance requires observations')
      return
    end if
    ub = .true.
    if (present(unbiased)) ub = unbiased
    if (ub .and. n < 2) then
      allocate(cov(0,0))
      call set_error(err, tvmvp_invalid_input, 'unbiased covariance requires two observations')
      return
    end if
    allocate(centered(n,size(x,2)))
    do j = 1, size(x,2)
      centered(:,j) = x(:,j) - sample_mean(x(:,j))
    end do
    denom = real(merge(n-1,n,ub),dp)
    allocate(cov(size(x,2),size(x,2)))
    cov = matmul(transpose(centered),centered)/denom
    cov = 0.5_dp*(cov+transpose(cov))
  end subroutine covariance_matrix

  subroutine solve_linear(a, b, x, err)
    real(dp), intent(in) :: a(:,:), b(:)
    real(dp), allocatable, intent(out) :: x(:)
    type(tvmvp_error), intent(out) :: err
    real(dp), allocatable :: aa(:,:), bb(:)
    real(dp) :: pivot, factor, scale, tmp
    integer :: n, i, j, k, p
    call clear_error(err)
    n = size(b)
    if (size(a,1) /= n .or. size(a,2) /= n) then
      allocate(x(0))
      call set_error(err,tvmvp_invalid_input,'linear system dimensions do not conform')
      return
    end if
    allocate(aa(n,n),bb(n),x(n))
    aa = a
    bb = b
    scale = max(1.0_dp,maxval(abs(aa)))
    do k = 1, n
      p = k-1+maxloc(abs(aa(k:n,k)),dim=1)
      if (abs(aa(p,k)) <= epsilon(1.0_dp)*scale*real(max(1,n),dp)) then
        x = 0.0_dp
        call set_error(err,tvmvp_singular,'singular linear system')
        return
      end if
      if (p /= k) then
        do j = k, n
          tmp = aa(k,j); aa(k,j)=aa(p,j); aa(p,j)=tmp
        end do
        tmp = bb(k); bb(k)=bb(p); bb(p)=tmp
      end if
      pivot = aa(k,k)
      do i = k+1,n
        factor = aa(i,k)/pivot
        aa(i,k)=0.0_dp
        if (k < n) aa(i,k+1:n)=aa(i,k+1:n)-factor*aa(k,k+1:n)
        bb(i)=bb(i)-factor*bb(k)
      end do
    end do
    do i=n,1,-1
      if (i<n) then
        x(i)=(bb(i)-dot_product(aa(i,i+1:n),x(i+1:n)))/aa(i,i)
      else
        x(i)=bb(i)/aa(i,i)
      end if
    end do
  end subroutine solve_linear

  subroutine solve_matrix(a, b, x, err)
    real(dp), intent(in) :: a(:,:), b(:,:)
    real(dp), allocatable, intent(out) :: x(:,:)
    type(tvmvp_error), intent(out) :: err
    type(tvmvp_error) :: local_err
    real(dp), allocatable :: col(:)
    integer :: j
    call clear_error(err)
    if (size(a,1)/=size(a,2) .or. size(b,1)/=size(a,1)) then
      allocate(x(0,0))
      call set_error(err,tvmvp_invalid_input,'matrix system dimensions do not conform')
      return
    end if
    allocate(x(size(a,2),size(b,2)))
    do j=1,size(b,2)
      call solve_linear(a,b(:,j),col,local_err)
      if (local_err%failed()) then
        x=0.0_dp; err=local_err; return
      end if
      x(:,j)=col
    end do
  end subroutine solve_matrix

  subroutine invert_matrix(a, ainv, err)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: ainv(:,:)
    type(tvmvp_error), intent(out) :: err
    real(dp), allocatable :: eye(:,:)
    if (size(a,1)/=size(a,2)) then
      allocate(ainv(0,0))
      call set_error(err,tvmvp_invalid_input,'matrix inverse requires a square matrix')
      return
    end if
    allocate(eye(size(a,1),size(a,1)))
    eye=identity_matrix(size(a,1))
    call solve_matrix(a,eye,ainv,err)
  end subroutine invert_matrix

  subroutine cholesky_factor(a, l, err)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: l(:,:)
    type(tvmvp_error), intent(out) :: err
    real(dp) :: s, scale
    integer :: n,i,j,k
    call clear_error(err)
    n=size(a,1)
    if (size(a,2)/=n) then
      allocate(l(0,0)); call set_error(err,tvmvp_invalid_input,'Cholesky requires square matrix'); return
    end if
    allocate(l(n,n)); l=0.0_dp
    scale=max(1.0_dp,maxval(abs(a)))
    do i=1,n
      do j=1,i
        s=a(i,j)
        do k=1,j-1
          s=s-l(i,k)*l(j,k)
        end do
        if (i==j) then
          if (s <= epsilon(1.0_dp)*scale*real(max(1,n),dp)) then
            call set_error(err,tvmvp_not_positive_definite,'matrix is not positive definite')
            return
          end if
          l(i,j)=sqrt(s)
        else
          l(i,j)=s/l(j,j)
        end if
      end do
    end do
  end subroutine cholesky_factor

  subroutine symmetric_eigen(a, values, vectors, err, tolerance, max_sweeps)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: values(:), vectors(:,:)
    type(tvmvp_error), intent(out) :: err
    real(dp), intent(in), optional :: tolerance
    integer, intent(in), optional :: max_sweeps
    real(dp), allocatable :: work(:,:)
    real(dp) :: tol, app, aqq, apq, tau, t, c, s, wip, wiq, vip, viq, off
    integer :: n, p, q, i, sweep, nsweep, j, k
    call clear_error(err)
    n=size(a,1)
    if (size(a,2)/=n) then
      allocate(values(0),vectors(0,0)); call set_error(err,tvmvp_invalid_input,'eigen decomposition requires square matrix'); return
    end if
    tol=100.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(a)))
    if (present(tolerance)) tol=tolerance
    nsweep=max(20,10*n*n)
    if (present(max_sweeps)) nsweep=max_sweeps
    allocate(work(n,n),values(n),vectors(n,n))
    work=0.5_dp*(a+transpose(a)); vectors=identity_matrix(n)
    off=0.0_dp
    do sweep=1,nsweep
      off=0.0_dp
      do p=1,n-1
        do q=p+1,n
          off=max(off,abs(work(p,q)))
          apq=work(p,q)
          if (abs(apq)<=tol) cycle
          app=work(p,p); aqq=work(q,q)
          tau=(aqq-app)/(2.0_dp*apq)
          if (tau>=0.0_dp) then
            t=1.0_dp/(tau+sqrt(1.0_dp+tau*tau))
          else
            t=-1.0_dp/(-tau+sqrt(1.0_dp+tau*tau))
          end if
          c=1.0_dp/sqrt(1.0_dp+t*t); s=t*c
          do i=1,n
            if (i/=p .and. i/=q) then
              wip=work(i,p); wiq=work(i,q)
              work(i,p)=c*wip-s*wiq; work(p,i)=work(i,p)
              work(i,q)=s*wip+c*wiq; work(q,i)=work(i,q)
            end if
          end do
          work(p,p)=c*c*app-2.0_dp*s*c*apq+s*s*aqq
          work(q,q)=s*s*app+2.0_dp*s*c*apq+c*c*aqq
          work(p,q)=0.0_dp; work(q,p)=0.0_dp
          do i=1,n
            vip=vectors(i,p); viq=vectors(i,q)
            vectors(i,p)=c*vip-s*viq
            vectors(i,q)=s*vip+c*viq
          end do
        end do
      end do
      if (off<=tol) exit
    end do
    if (sweep>nsweep .and. off>tol) then
      call set_error(err,tvmvp_convergence_failure,'Jacobi eigen decomposition did not converge')
      return
    end if
    do i=1,n
      values(i)=work(i,i)
    end do
    do i=1,n-1
      k=i
      do j=i+1,n
        if (values(j)>values(k)) k=j
      end do
      if (k/=i) then
        app=values(i); values(i)=values(k); values(k)=app
        do j=1,n
          app=vectors(j,i); vectors(j,i)=vectors(j,k); vectors(j,k)=app
        end do
      end if
    end do
  end subroutine symmetric_eigen

  subroutine floor_symmetric_eigenvalues(a, floor_value, out, err)
    real(dp), intent(in) :: a(:,:), floor_value
    real(dp), allocatable, intent(out) :: out(:,:)
    type(tvmvp_error), intent(out) :: err
    real(dp), allocatable :: values(:), vectors(:,:), scaled(:,:)
    integer :: j
    call symmetric_eigen(a,values,vectors,err)
    if (err%failed()) then
      allocate(out(0,0)); return
    end if
    allocate(scaled(size(vectors,1),size(vectors,2)))
    scaled=vectors
    do j=1,size(values)
      scaled(:,j)=scaled(:,j)*max(values(j),floor_value)
    end do
    allocate(out(size(a,1),size(a,2)))
    out=matmul(scaled,transpose(vectors))
    out=0.5_dp*(out+transpose(out))
  end subroutine floor_symmetric_eigenvalues

  subroutine matrix_sqrt_abs(a, root, err)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: root(:,:)
    type(tvmvp_error), intent(out) :: err
    real(dp), allocatable :: values(:), vectors(:,:), scaled(:,:)
    integer :: j
    call symmetric_eigen(a,values,vectors,err)
    if (err%failed()) then
      allocate(root(0,0)); return
    end if
    allocate(scaled(size(vectors,1),size(vectors,2)))
    scaled=vectors
    do j=1,size(values)
      scaled(:,j)=scaled(:,j)*sqrt(abs(values(j)))
    end do
    allocate(root(size(a,1),size(a,2)))
    root=matmul(scaled,transpose(vectors))
    root=0.5_dp*(root+transpose(root))
  end subroutine matrix_sqrt_abs
end module tvmvp_linalg
