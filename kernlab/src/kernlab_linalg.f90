! SPDX-License-Identifier: GPL-2.0-only
module kernlab_linalg
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use kernlab_kinds, only: dp, KL_SUCCESS, KL_INVALID_ARGUMENT, KL_SINGULAR, KL_NOT_CONVERGED
  implicit none
  private
  public :: solve_linear, invert_matrix, cholesky_lower, symmetric_eigen
  public :: center_kernel, standardize_columns, kmeans_dense, quantile_value
  public :: sort_eigenpairs_desc, matrix_inverse_sqrt, vector_ranks, vec_norm

contains

  pure real(dp) function vec_norm(x)
    real(dp), intent(in) :: x(:)
    vec_norm = sqrt(max(0.0_dp, dot_product(x, x)))
  end function vec_norm

  subroutine solve_linear(a, b, x, status)
    real(dp), intent(in) :: a(:,:), b(:,:)
    real(dp), allocatable, intent(out) :: x(:,:)
    integer, intent(out) :: status
    real(dp), allocatable :: aug(:,:)
    real(dp) :: pivot, factor, tmp
    integer :: n, nrhs, i, j, k, p

    n = size(a, 1)
    nrhs = size(b, 2)
    status = KL_INVALID_ARGUMENT
    allocate(x(0,0))
    if (size(a,2) /= n .or. size(b,1) /= n .or. n == 0) return
    if (.not. all(ieee_is_finite(a)) .or. .not. all(ieee_is_finite(b))) return
    allocate(aug(n,n+nrhs))
    aug(:,1:n) = a
    aug(:,n+1:n+nrhs) = b

    do k = 1, n
      p = k
      do i = k + 1, n
        if (abs(aug(i,k)) > abs(aug(p,k))) p = i
      end do
      if (abs(aug(p,k)) <= 100.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(a)))) then
        status = KL_SINGULAR
        return
      end if
      if (p /= k) then
        do j = k, n + nrhs
          tmp = aug(k,j); aug(k,j) = aug(p,j); aug(p,j) = tmp
        end do
      end if
      pivot = aug(k,k)
      aug(k,k:n+nrhs) = aug(k,k:n+nrhs)/pivot
      do i = 1, n
        if (i == k) cycle
        factor = aug(i,k)
        if (abs(factor) > tiny(1.0_dp)) aug(i,k:n+nrhs) = aug(i,k:n+nrhs) - factor*aug(k,k:n+nrhs)
      end do
    end do
    deallocate(x)
    allocate(x(n,nrhs))
    x = aug(:,n+1:n+nrhs)
    status = KL_SUCCESS
  end subroutine solve_linear

  subroutine invert_matrix(a, ainv, status)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: ainv(:,:)
    integer, intent(out) :: status
    real(dp), allocatable :: eye(:,:)
    integer :: n, i
    n = size(a,1)
    allocate(eye(n,n)); eye = 0.0_dp
    do i = 1, n
      eye(i,i) = 1.0_dp
    end do
    call solve_linear(a, eye, ainv, status)
  end subroutine invert_matrix

  subroutine cholesky_lower(a, l, status, jitter)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: l(:,:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: jitter
    real(dp) :: s, jit
    integer :: n, i, j, k
    n = size(a,1)
    status = KL_INVALID_ARGUMENT
    allocate(l(0,0))
    if (size(a,2) /= n .or. n == 0) return
    jit = 0.0_dp
    if (present(jitter)) jit = max(0.0_dp,jitter)
    deallocate(l); allocate(l(n,n)); l = 0.0_dp
    do i = 1, n
      do j = 1, i
        s = a(i,j)
        if (i == j) s = s + jit
        do k = 1, j-1
          s = s - l(i,k)*l(j,k)
        end do
        if (i == j) then
          if (s <= 0.0_dp) then
            status = KL_SINGULAR
            return
          end if
          l(i,j) = sqrt(s)
        else
          l(i,j) = s/l(j,j)
        end if
      end do
    end do
    status = KL_SUCCESS
  end subroutine cholesky_lower

  subroutine symmetric_eigen(a, values, vectors, status, max_sweeps, tol)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: values(:), vectors(:,:)
    integer, intent(out) :: status
    integer, intent(in), optional :: max_sweeps
    real(dp), intent(in), optional :: tol
    real(dp), allocatable :: d(:,:)
    real(dp) :: app, aqq, apq, tau, t, c, s, thresh, off, dpi, dqi, vip, viq
    integer :: n, p, q, i, sweep, ms

    n = size(a,1)
    status = KL_INVALID_ARGUMENT
    allocate(values(0),vectors(0,0))
    if (n == 0 .or. size(a,2) /= n) return
    allocate(d(n,n)); d = 0.5_dp*(a+transpose(a))
    deallocate(values,vectors); allocate(values(n),vectors(n,n)); vectors=0.0_dp
    do i=1,n; vectors(i,i)=1.0_dp; end do
    ms = max(100, 100*n*n); if (present(max_sweeps)) ms=max_sweeps
    thresh = sqrt(epsilon(1.0_dp)); if (present(tol)) thresh=tol
    status = KL_NOT_CONVERGED
    do sweep=1,ms
      off=0.0_dp
      do p=1,n-1
        do q=p+1,n
          off=max(off,abs(d(p,q)))
          if (abs(d(p,q)) <= thresh*max(1.0_dp,abs(d(p,p))+abs(d(q,q)))) cycle
          app=d(p,p); aqq=d(q,q); apq=d(p,q)
          tau=(aqq-app)/(2.0_dp*apq)
          if (tau >= 0.0_dp) then
            t=1.0_dp/(tau+sqrt(1.0_dp+tau*tau))
          else
            t=-1.0_dp/(-tau+sqrt(1.0_dp+tau*tau))
          end if
          c=1.0_dp/sqrt(1.0_dp+t*t); s=t*c
          do i=1,n
            if (i /= p .and. i /= q) then
              dpi=d(i,p); dqi=d(i,q)
              d(i,p)=c*dpi-s*dqi; d(p,i)=d(i,p)
              d(i,q)=s*dpi+c*dqi; d(q,i)=d(i,q)
            end if
            vip=vectors(i,p); viq=vectors(i,q)
            vectors(i,p)=c*vip-s*viq
            vectors(i,q)=s*vip+c*viq
          end do
          d(p,p)=c*c*app-2.0_dp*s*c*apq+s*s*aqq
          d(q,q)=s*s*app+2.0_dp*s*c*apq+c*c*aqq
          d(p,q)=0.0_dp; d(q,p)=0.0_dp
        end do
      end do
      if (off <= thresh*max(1.0_dp,maxval(abs(d)))) then
        status=KL_SUCCESS; exit
      end if
    end do
    do i=1,n; values(i)=d(i,i); end do
    call sort_eigenpairs_desc(values,vectors)
  end subroutine symmetric_eigen

  subroutine sort_eigenpairs_desc(values, vectors)
    real(dp), intent(inout) :: values(:), vectors(:,:)
    real(dp) :: tv
    real(dp), allocatable :: col(:)
    integer :: i,j,k,n
    n=size(values); allocate(col(size(vectors,1)))
    do i=1,n-1
      k=i
      do j=i+1,n
        if(values(j)>values(k)) k=j
      end do
      if(k/=i) then
        tv=values(i); values(i)=values(k); values(k)=tv
        col=vectors(:,i); vectors(:,i)=vectors(:,k); vectors(:,k)=col
      end if
    end do
  end subroutine sort_eigenpairs_desc

  subroutine matrix_inverse_sqrt(a, out, status, floor_value)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: out(:,:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: floor_value
    real(dp), allocatable :: val(:), vec(:,:), scaled(:,:)
    real(dp) :: fl
    integer :: i,n
    call symmetric_eigen(a,val,vec,status)
    if(status/=KL_SUCCESS) then; allocate(out(0,0)); return; end if
    fl=1.0e-12_dp; if(present(floor_value)) fl=floor_value
    n=size(val); allocate(scaled(n,n)); scaled=vec
    do i=1,n
      if(val(i)>fl) then
        scaled(:,i)=scaled(:,i)/sqrt(val(i))
      else
        scaled(:,i)=0.0_dp
      end if
    end do
    allocate(out(n,n)); out=matmul(scaled,transpose(vec))
  end subroutine matrix_inverse_sqrt

  subroutine center_kernel(k, kc, row_mean, grand_mean)
    real(dp), intent(in) :: k(:,:)
    real(dp), allocatable, intent(out) :: kc(:,:), row_mean(:)
    real(dp), intent(out) :: grand_mean
    integer :: n,i,j
    n=size(k,1); allocate(row_mean(n),kc(n,n))
    row_mean=sum(k,dim=2)/real(n,dp)
    grand_mean=sum(row_mean)/real(n,dp)
    do i=1,n; do j=1,n
      kc(i,j)=k(i,j)-row_mean(i)-row_mean(j)+grand_mean
    end do; end do
  end subroutine center_kernel

  subroutine standardize_columns(x, z, means, scales)
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable, intent(out) :: z(:,:), means(:), scales(:)
    integer :: n,p,j
    n=size(x,1); p=size(x,2); allocate(z(n,p),means(p),scales(p))
    means=sum(x,dim=1)/real(n,dp)
    do j=1,p
      if(n>1) then
        scales(j)=sqrt(sum((x(:,j)-means(j))**2)/real(n-1,dp))
      else
        scales(j)=1.0_dp
      end if
      if(scales(j)<=sqrt(epsilon(1.0_dp))) scales(j)=1.0_dp
      z(:,j)=(x(:,j)-means(j))/scales(j)
    end do
  end subroutine standardize_columns

  subroutine kmeans_dense(x, k, labels, centers, withinss, iterations, status, maxiter)
    real(dp), intent(in) :: x(:,:)
    integer, intent(in) :: k
    integer, allocatable, intent(out) :: labels(:)
    real(dp), allocatable, intent(out) :: centers(:,:), withinss(:)
    integer, intent(out) :: iterations, status
    integer, intent(in), optional :: maxiter
    integer :: n,p,i,c,best,mi,limit
    real(dp) :: d,bestd
    integer, allocatable :: old(:), counts(:)
    n=size(x,1); p=size(x,2); status=KL_INVALID_ARGUMENT; iterations=0
    allocate(labels(0),centers(0,0),withinss(0))
    if(k<1 .or. k>n .or. p<1) return
    deallocate(labels,centers,withinss); allocate(labels(n),old(n),centers(k,p),withinss(k),counts(k))
    do c=1,k
      i=1+(c-1)*max(1,(n-1)/max(1,k-1)); i=min(i,n); centers(c,:)=x(i,:)
    end do
    labels=0; limit=200; if(present(maxiter)) limit=maxiter
    do mi=1,limit
      old=labels
      do i=1,n
        best=1; bestd=sum((x(i,:)-centers(1,:))**2)
        do c=2,k
          d=sum((x(i,:)-centers(c,:))**2)
          if(d<bestd) then; bestd=d; best=c; end if
        end do
        labels(i)=best
      end do
      centers=0.0_dp; counts=0
      do i=1,n; centers(labels(i),:)=centers(labels(i),:)+x(i,:); counts(labels(i))=counts(labels(i))+1; end do
      do c=1,k
        if(counts(c)>0) then
          centers(c,:)=centers(c,:)/real(counts(c),dp)
        else
          best=1; bestd=-1.0_dp
          do i=1,n
            d=sum((x(i,:)-centers(labels(i),:))**2)
            if(d>bestd) then; bestd=d; best=i; end if
          end do
          centers(c,:)=x(best,:); labels(best)=c; counts(c)=1
        end if
      end do
      iterations=mi
      if(all(labels==old)) exit
    end do
    withinss=0.0_dp
    do i=1,n; withinss(labels(i))=withinss(labels(i))+sum((x(i,:)-centers(labels(i),:))**2); end do
    status=KL_SUCCESS
  end subroutine kmeans_dense

  real(dp) function quantile_value(x, prob)
    real(dp), intent(in) :: x(:), prob
    real(dp), allocatable :: y(:)
    real(dp) :: h,tmp
    integer :: i,j,n,lo,hi
    n=size(x); if(n==0) then; quantile_value=0.0_dp; return; end if
    allocate(y(n)); y=x
    do i=2,n
      tmp=y(i); j=i-1
      do while(j>=1)
        if(y(j)<=tmp) exit
        y(j+1)=y(j); j=j-1
      end do
      y(j+1)=tmp
    end do
    h=1.0_dp+(real(n-1,dp))*min(1.0_dp,max(0.0_dp,prob))
    lo=floor(h); hi=ceiling(h)
    quantile_value=y(lo)+(h-real(lo,dp))*(y(hi)-y(lo))
  end function quantile_value

  subroutine vector_ranks(x, ranks)
    real(dp), intent(in) :: x(:)
    real(dp), allocatable, intent(out) :: ranks(:)
    integer, allocatable :: idx(:)
    integer :: n,i,j,k,start,finish,it
    real(dp) :: tx
    n=size(x); allocate(ranks(n),idx(n)); idx=[(i,i=1,n)]
    do i=2,n
      it=idx(i); tx=x(it); j=i-1
      do while(j>=1)
        if(x(idx(j))<=tx) exit
        idx(j+1)=idx(j); j=j-1
      end do
      idx(j+1)=it
    end do
    i=1
    do while(i<=n)
      start=i; finish=i
      do while(finish<n)
        if(abs(x(idx(finish+1))-x(idx(start)))>10.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(x(idx(start))))) exit
        finish=finish+1
      end do
      do k=start,finish; ranks(idx(k))=0.5_dp*real(start+finish,dp); end do
      i=finish+1
    end do
  end subroutine vector_ranks

end module kernlab_linalg
