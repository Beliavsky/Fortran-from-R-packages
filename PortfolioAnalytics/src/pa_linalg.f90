! SPDX-License-Identifier: GPL-3.0-only
module pa_linalg
  use pa_kinds, only: dp, pa_pi
  implicit none
  private
  public :: solve_linear, inverse_matrix, cholesky_lower, jacobi_eigen_sym
  public :: project_box_sum, set_random_seed, random_normal, normal_quantile
  public :: normal_pdf, sort_real, clip_value, symmetrize, kronecker_product

contains

  pure real(dp) function clip_value(x, lo, hi) result(y)
    real(dp), intent(in) :: x, lo, hi
    y = min(max(x, lo), hi)
  end function clip_value

  subroutine symmetrize(a)
    real(dp), intent(inout) :: a(:,:)
    a = 0.5_dp * (a + transpose(a))
  end subroutine symmetrize

  subroutine solve_linear(a, b, x, info)
    real(dp), intent(in) :: a(:,:), b(:)
    real(dp), intent(out) :: x(:)
    integer, intent(out) :: info
    real(dp), allocatable :: m(:,:), rhs(:), rowtmp(:)
    real(dp) :: pivot, factor, tmp
    integer :: n, i, j, k, p

    n = size(b)
    info = 0
    x = 0.0_dp
    if (size(a,1) /= n .or. size(a,2) /= n .or. size(x) /= n) then
      info = -1
      return
    end if
    allocate(m(n,n), rhs(n), rowtmp(n))
    m = a
    rhs = b
    do k = 1, n - 1
      p = k
      pivot = abs(m(k,k))
      do i = k + 1, n
        if (abs(m(i,k)) > pivot) then
          pivot = abs(m(i,k))
          p = i
        end if
      end do
      if (pivot <= epsilon(1.0_dp) * max(1.0_dp, maxval(abs(m)))) then
        info = k
        return
      end if
      if (p /= k) then
        rowtmp = m(k,:)
        m(k,:) = m(p,:)
        m(p,:) = rowtmp
        tmp = rhs(k)
        rhs(k) = rhs(p)
        rhs(p) = tmp
      end if
      do i = k + 1, n
        factor = m(i,k) / m(k,k)
        m(i,k) = 0.0_dp
        do j = k + 1, n
          m(i,j) = m(i,j) - factor * m(k,j)
        end do
        rhs(i) = rhs(i) - factor * rhs(k)
      end do
    end do
    if (abs(m(n,n)) <= epsilon(1.0_dp) * max(1.0_dp, maxval(abs(m)))) then
      info = n
      return
    end if
    do i = n, 1, -1
      tmp = rhs(i)
      if (i < n) tmp = tmp - dot_product(m(i,i+1:n), x(i+1:n))
      x(i) = tmp / m(i,i)
    end do
  end subroutine solve_linear

  subroutine inverse_matrix(a, ainv, info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: ainv(:,:)
    integer, intent(out) :: info
    real(dp), allocatable :: e(:), x(:)
    integer :: n, j, istat

    n = size(a,1)
    info = 0
    if (size(a,2) /= n .or. size(ainv,1) /= n .or. size(ainv,2) /= n) then
      info = -1
      return
    end if
    allocate(e(n), x(n))
    do j = 1, n
      e = 0.0_dp
      e(j) = 1.0_dp
      call solve_linear(a, e, x, istat)
      if (istat /= 0) then
        info = istat
        return
      end if
      ainv(:,j) = x
    end do
  end subroutine inverse_matrix

  subroutine cholesky_lower(a, l, info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: l(:,:)
    integer, intent(out) :: info
    integer :: n, i, j, k
    real(dp) :: s

    n = size(a,1)
    info = 0
    l = 0.0_dp
    if (size(a,2) /= n .or. size(l,1) /= n .or. size(l,2) /= n) then
      info = -1
      return
    end if
    do i = 1, n
      do j = 1, i
        s = a(i,j)
        do k = 1, j - 1
          s = s - l(i,k) * l(j,k)
        end do
        if (i == j) then
          if (s <= 0.0_dp) then
            info = i
            return
          end if
          l(i,j) = sqrt(s)
        else
          l(i,j) = s / l(j,j)
        end if
      end do
    end do
  end subroutine cholesky_lower

  subroutine jacobi_eigen_sym(a, eigenvalues, eigenvectors, info, tolerance, max_sweeps)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: eigenvalues(:), eigenvectors(:,:)
    integer, intent(out) :: info
    real(dp), intent(in), optional :: tolerance
    integer, intent(in), optional :: max_sweeps
    real(dp), allocatable :: b(:,:)
    real(dp) :: tol, app, aqq, apq, tau, t, c, s, tmp, maxoff
    integer :: n, sweeps, p, q, i, j, maxit, k

    n = size(a,1)
    info = 0
    if (size(a,2) /= n .or. size(eigenvalues) /= n .or. &
        size(eigenvectors,1) /= n .or. size(eigenvectors,2) /= n) then
      info = -1
      return
    end if
    tol = 1.0e-12_dp
    if (present(tolerance)) tol = tolerance
    maxit = max(50, 20*n*n)
    if (present(max_sweeps)) maxit = max_sweeps
    allocate(b(n,n))
    b = 0.5_dp * (a + transpose(a))
    eigenvectors = 0.0_dp
    do i = 1, n
      eigenvectors(i,i) = 1.0_dp
    end do

    do sweeps = 1, maxit
      maxoff = 0.0_dp
      p = 1
      q = min(2,n)
      do i = 1, n - 1
        do j = i + 1, n
          if (abs(b(i,j)) > maxoff) then
            maxoff = abs(b(i,j))
            p = i
            q = j
          end if
        end do
      end do
      if (maxoff <= tol * max(1.0_dp, maxval(abs(b)))) exit
      app = b(p,p)
      aqq = b(q,q)
      apq = b(p,q)
      tau = (aqq - app) / (2.0_dp * apq)
      if (tau >= 0.0_dp) then
        t = 1.0_dp / (tau + sqrt(1.0_dp + tau*tau))
      else
        t = -1.0_dp / (-tau + sqrt(1.0_dp + tau*tau))
      end if
      c = 1.0_dp / sqrt(1.0_dp + t*t)
      s = t * c
      do k = 1, n
        if (k /= p .and. k /= q) then
          tmp = b(k,p)
          b(k,p) = c*tmp - s*b(k,q)
          b(p,k) = b(k,p)
          b(k,q) = s*tmp + c*b(k,q)
          b(q,k) = b(k,q)
        end if
      end do
      b(p,p) = c*c*app - 2.0_dp*s*c*apq + s*s*aqq
      b(q,q) = s*s*app + 2.0_dp*s*c*apq + c*c*aqq
      b(p,q) = 0.0_dp
      b(q,p) = 0.0_dp
      do k = 1, n
        tmp = eigenvectors(k,p)
        eigenvectors(k,p) = c*tmp - s*eigenvectors(k,q)
        eigenvectors(k,q) = s*tmp + c*eigenvectors(k,q)
      end do
    end do
    if (sweeps > maxit) info = 1
    do i = 1, n
      eigenvalues(i) = b(i,i)
    end do
    call sort_eigenpairs_desc(eigenvalues, eigenvectors)
  end subroutine jacobi_eigen_sym

  subroutine sort_eigenpairs_desc(values, vectors)
    real(dp), intent(inout) :: values(:), vectors(:,:)
    real(dp) :: tv
    real(dp), allocatable :: col(:)
    integer :: i, j, k, n
    n = size(values)
    allocate(col(size(vectors,1)))
    do i = 1, n - 1
      k = i
      do j = i + 1, n
        if (values(j) > values(k)) k = j
      end do
      if (k /= i) then
        tv = values(i)
        values(i) = values(k)
        values(k) = tv
        col = vectors(:,i)
        vectors(:,i) = vectors(:,k)
        vectors(:,k) = col
      end if
    end do
  end subroutine sort_eigenpairs_desc

  subroutine project_box_sum(v, lo, hi, target, w, feasible)
    real(dp), intent(in) :: v(:), lo(:), hi(:), target
    real(dp), intent(out) :: w(:)
    logical, intent(out) :: feasible
    real(dp) :: left, right, mid, s
    integer :: i, iter, n

    n = size(v)
    feasible = .false.
    w = 0.0_dp
    if (size(lo) /= n .or. size(hi) /= n .or. size(w) /= n) return
    if (any(lo > hi)) return
    if (target < sum(lo) - 1.0e-10_dp .or. target > sum(hi) + 1.0e-10_dp) return
    left = minval(v - hi) - 1.0_dp
    right = maxval(v - lo) + 1.0_dp
    do iter = 1, 200
      mid = 0.5_dp * (left + right)
      do i = 1, n
        w(i) = clip_value(v(i) - mid, lo(i), hi(i))
      end do
      s = sum(w)
      if (abs(s - target) <= 1.0e-13_dp * max(1.0_dp, abs(target))) exit
      if (s > target) then
        left = mid
      else
        right = mid
      end if
    end do
    do i = 1, n
      w(i) = clip_value(v(i) - mid, lo(i), hi(i))
    end do
    if (abs(sum(w) - target) > 1.0e-10_dp) then
      if (sum(w) < target) then
        i = maxloc(hi - w, dim=1)
      else
        i = maxloc(w - lo, dim=1)
      end if
      w(i) = w(i) + target - sum(w)
    end if
    feasible = all(w >= lo - 1.0e-9_dp) .and. all(w <= hi + 1.0e-9_dp)
  end subroutine project_box_sum

  subroutine set_random_seed(seed)
    integer, intent(in) :: seed
    integer :: n, i
    integer, allocatable :: put(:)
    call random_seed(size=n)
    allocate(put(n))
    do i = 1, n
      put(i) = modulo(seed + 104729*i + 37*i*i, huge(1) - 1)
      if (put(i) <= 0) put(i) = i + 1
    end do
    call random_seed(put=put)
  end subroutine set_random_seed

  real(dp) function random_normal() result(z)
    real(dp) :: u1, u2
    call random_number(u1)
    call random_number(u2)
    u1 = max(u1, tiny(1.0_dp))
    z = sqrt(-2.0_dp*log(u1)) * cos(2.0_dp*pa_pi*u2)
  end function random_normal

  pure real(dp) function normal_pdf(x) result(y)
    real(dp), intent(in) :: x
    y = exp(-0.5_dp*x*x) / sqrt(2.0_dp*pa_pi)
  end function normal_pdf

  pure real(dp) function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    real(dp), parameter :: a1=-3.969683028665376e+01_dp, a2=2.209460984245205e+02_dp
    real(dp), parameter :: a3=-2.759285104469687e+02_dp, a4=1.383577518672690e+02_dp
    real(dp), parameter :: a5=-3.066479806614716e+01_dp, a6=2.506628277459239e+00_dp
    real(dp), parameter :: b1=-5.447609879822406e+01_dp, b2=1.615858368580409e+02_dp
    real(dp), parameter :: b3=-1.556989798598866e+02_dp, b4=6.680131188771972e+01_dp
    real(dp), parameter :: b5=-1.328068155288572e+01_dp
    real(dp), parameter :: c1=-7.784894002430293e-03_dp, c2=-3.223964580411365e-01_dp
    real(dp), parameter :: c3=-2.400758277161838e+00_dp, c4=-2.549732539343734e+00_dp
    real(dp), parameter :: c5=4.374664141464968e+00_dp, c6=2.938163982698783e+00_dp
    real(dp), parameter :: d1=7.784695709041462e-03_dp, d2=3.224671290700398e-01_dp
    real(dp), parameter :: d3=2.445134137142996e+00_dp, d4=3.754408661907416e+00_dp
    real(dp), parameter :: plow=0.02425_dp, phigh=1.0_dp-plow
    real(dp) :: q, r
    if (p <= 0.0_dp) then
      x = -huge(1.0_dp)
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
    else if (p < plow) then
      q = sqrt(-2.0_dp*log(p))
      x = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
    else if (p <= phigh) then
      q = p - 0.5_dp
      r = q*q
      x = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q / (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.0_dp)
    else
      q = sqrt(-2.0_dp*log(1.0_dp-p))
      x = -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
    end if
  end function normal_quantile

  subroutine sort_real(x)
    real(dp), intent(inout) :: x(:)
    integer :: i, j
    real(dp) :: key
    do i = 2, size(x)
      key = x(i)
      j = i - 1
      do while (j >= 1)
        if (x(j) <= key) exit
        x(j+1) = x(j)
        j = j - 1
      end do
      x(j+1) = key
    end do
  end subroutine sort_real

  subroutine kronecker_product(a, b, c)
    real(dp), intent(in) :: a(:,:), b(:,:)
    real(dp), intent(out) :: c(:,:)
    integer :: i, j, rb, cb
    rb = size(b,1)
    cb = size(b,2)
    if (size(c,1) /= size(a,1)*rb .or. size(c,2) /= size(a,2)*cb) return
    do i = 1, size(a,1)
      do j = 1, size(a,2)
        c((i-1)*rb+1:i*rb, (j-1)*cb+1:j*cb) = a(i,j) * b
      end do
    end do
  end subroutine kronecker_product

end module pa_linalg
