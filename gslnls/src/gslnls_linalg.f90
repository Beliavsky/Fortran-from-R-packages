! SPDX-License-Identifier: LGPL-3.0-only
module gslnls_linalg
  use gslnls_kinds, only : dp
  implicit none
  private
  public :: norm2v, median_value, solve_linear, cholesky_factor, cholesky_solve
  public :: least_squares, symmetric_eigen, pseudo_inverse_sym, matrix_rank
  public :: covariance_from_jacobian, determinant_jtj

contains

  pure real(dp) function norm2v(x) result(v)
    real(dp), intent(in) :: x(:)
    real(dp) :: scale, ssq, ax
    integer :: i
    scale = 0.0_dp
    ssq = 1.0_dp
    do i = 1, size(x)
      if (abs(x(i)) > 0.0_dp) then
        ax = abs(x(i))
        if (scale < ax) then
          ssq = 1.0_dp + ssq * (scale / ax)**2
          scale = ax
        else
          ssq = ssq + (ax / scale)**2
        end if
      end if
    end do
    if (scale <= tiny(1.0_dp)) then
      v = 0.0_dp
    else
      v = scale * sqrt(ssq)
    end if
  end function norm2v

  function median_value(x) result(med)
    real(dp), intent(in) :: x(:)
    real(dp) :: med
    real(dp), allocatable :: a(:)
    real(dp) :: tmp
    integer :: i, j, n
    n = size(x)
    if (n == 0) then
      med = 0.0_dp
      return
    end if
    allocate(a(n)); a = x
    do i = 2, n
      tmp = a(i); j = i - 1
      do while (j >= 1)
        if (a(j) <= tmp) exit
        a(j+1) = a(j); j = j - 1
      end do
      a(j+1) = tmp
    end do
    if (mod(n,2) == 1) then
      med = a((n+1)/2)
    else
      med = 0.5_dp * (a(n/2) + a(n/2+1))
    end if
  end function median_value

  subroutine solve_linear(a, b, x, ok)
    real(dp), intent(in) :: a(:,:), b(:)
    real(dp), intent(out) :: x(:)
    logical, intent(out) :: ok
    real(dp), allocatable :: m(:,:), rhs(:)
    real(dp) :: piv, fac, t
    integer :: n, i, j, k, ip
    n = size(b)
    ok = .false.
    if (size(a,1) /= n .or. size(a,2) /= n .or. size(x) /= n) return
    allocate(m(n,n), rhs(n)); m = a; rhs = b
    do k = 1, n
      ip = k
      piv = abs(m(k,k))
      do i = k + 1, n
        if (abs(m(i,k)) > piv) then
          piv = abs(m(i,k)); ip = i
        end if
      end do
      if (piv <= 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(m)))) return
      if (ip /= k) then
        do j = k, n
          t = m(k,j); m(k,j) = m(ip,j); m(ip,j) = t
        end do
        t = rhs(k); rhs(k) = rhs(ip); rhs(ip) = t
      end if
      do i = k + 1, n
        fac = m(i,k) / m(k,k)
        m(i,k) = 0.0_dp
        m(i,k+1:n) = m(i,k+1:n) - fac * m(k,k+1:n)
        rhs(i) = rhs(i) - fac * rhs(k)
      end do
    end do
    x = 0.0_dp
    do i = n, 1, -1
      if (i < n) then
        x(i) = (rhs(i) - dot_product(m(i,i+1:n), x(i+1:n))) / m(i,i)
      else
        x(i) = rhs(i) / m(i,i)
      end if
    end do
    ok = .true.
  end subroutine solve_linear

  subroutine cholesky_factor(a, l, ok)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: l(:,:)
    logical, intent(out) :: ok
    real(dp) :: s
    integer :: n, i, j, k
    n = size(a,1)
    ok = .false.; l = 0.0_dp
    if (size(a,2) /= n .or. size(l,1) /= n .or. size(l,2) /= n) return
    do i = 1, n
      do j = 1, i
        s = a(i,j)
        do k = 1, j - 1
          s = s - l(i,k) * l(j,k)
        end do
        if (i == j) then
          if (s <= 0.0_dp) return
          l(i,j) = sqrt(s)
        else
          l(i,j) = s / l(j,j)
        end if
      end do
    end do
    ok = .true.
  end subroutine cholesky_factor

  subroutine cholesky_solve(a, b, x, ok)
    real(dp), intent(in) :: a(:,:), b(:)
    real(dp), intent(out) :: x(:)
    logical, intent(out) :: ok
    real(dp), allocatable :: l(:,:), y(:)
    integer :: n, i
    n = size(b)
    allocate(l(n,n), y(n))
    call cholesky_factor(a, l, ok)
    if (.not. ok) return
    do i = 1, n
      if (i == 1) then
        y(i) = b(i) / l(i,i)
      else
        y(i) = (b(i) - dot_product(l(i,1:i-1), y(1:i-1))) / l(i,i)
      end if
    end do
    do i = n, 1, -1
      if (i == n) then
        x(i) = y(i) / l(i,i)
      else
        x(i) = (y(i) - dot_product(l(i+1:n,i), x(i+1:n))) / l(i,i)
      end if
    end do
  end subroutine cholesky_solve

  subroutine least_squares(a, b, x, rank, solver)
    real(dp), intent(in) :: a(:,:), b(:)
    real(dp), intent(out) :: x(:)
    integer, intent(out) :: rank
    integer, intent(in) :: solver
    real(dp), allocatable :: ata(:,:), atb(:), q(:,:), r(:,:), aw(:,:), y(:), norms(:)
    real(dp), allocatable :: eval(:), evec(:,:), z(:)
    real(dp) :: tol, nv, tmp
    integer, allocatable :: perm(:)
    integer :: m, n, i, j, k, jp
    logical :: ok

    m = size(a,1); n = size(a,2)
    x = 0.0_dp; rank = 0
    if (size(b) /= m .or. size(x) /= n) return
    allocate(ata(n,n), atb(n))
    ata = matmul(transpose(a), a)
    atb = matmul(transpose(a), b)
    tol = max(100.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(ata))), tiny(1.0_dp))

    select case (solver)
    case (2)
      call cholesky_solve(ata, atb, x, ok)
      if (ok) then
        rank = n
      else
        block
          real(dp), allocatable :: atai(:,:)
          allocate(atai(n,n))
          call pseudo_inverse_sym(ata, atai, rank)
          x = matmul(atai, atb)
        end block
      end if
    case (3)
      allocate(eval(n), evec(n,n), z(n))
      call symmetric_eigen(ata, eval, evec)
      rank = count(eval > tol)
      z = matmul(transpose(evec), atb)
      do i = 1, n
        if (eval(i) > tol) then
          z(i) = z(i) / eval(i)
        else
          z(i) = 0.0_dp
        end if
      end do
      x = matmul(evec, z)
    case default
      allocate(aw(m,n), q(m,n), r(n,n), y(n), norms(n), perm(n))
      aw = a; q = 0.0_dp; r = 0.0_dp; y = 0.0_dp
      do j = 1, n
        norms(j) = dot_product(aw(:,j), aw(:,j)); perm(j) = j
      end do
      do k = 1, min(m,n)
        jp = k
        do j = k + 1, n
          if (norms(j) > norms(jp)) jp = j
        end do
        if (jp /= k) then
          do i = 1, m
            tmp = aw(i,k); aw(i,k) = aw(i,jp); aw(i,jp) = tmp
          end do
          tmp = norms(k); norms(k) = norms(jp); norms(jp) = tmp
          i = perm(k); perm(k) = perm(jp); perm(jp) = i
          if (k > 1) then
            do i = 1, k - 1
              tmp = r(i,k); r(i,k) = r(i,jp); r(i,jp) = tmp
            end do
          end if
        end if
        nv = norm2v(aw(:,k))
        if (nv <= sqrt(tol)) exit
        rank = rank + 1
        q(:,k) = aw(:,k) / nv
        r(k,k) = nv
        y(k) = dot_product(q(:,k), b)
        do j = k + 1, n
          r(k,j) = dot_product(q(:,k), aw(:,j))
          aw(:,j) = aw(:,j) - r(k,j) * q(:,k)
          norms(j) = dot_product(aw(:,j), aw(:,j))
        end do
      end do
      allocate(z(n)); z = 0.0_dp
      do i = rank, 1, -1
        if (i < rank) then
          z(i) = (y(i) - dot_product(r(i,i+1:rank), z(i+1:rank))) / r(i,i)
        else
          z(i) = y(i) / r(i,i)
        end if
      end do
      do i = 1, n
        x(perm(i)) = z(i)
      end do
    end select
  end subroutine least_squares

  subroutine symmetric_eigen(a, eval, evec)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: eval(:), evec(:,:)
    real(dp), allocatable :: w(:,:)
    real(dp) :: app, aqq, apq, phi, c, s, t1, t2, off, tol, tmp
    integer :: n, i, j, p, q, sweep, k
    n = size(a,1)
    allocate(w(n,n)); w = 0.5_dp * (a + transpose(a))
    evec = 0.0_dp
    do i = 1, n; evec(i,i) = 1.0_dp; end do
    tol = 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(w)))
    do sweep = 1, max(50, 20*n*n)
      off = 0.0_dp; p = 1; q = min(2,n)
      do i = 1, n - 1
        do j = i + 1, n
          if (abs(w(i,j)) > off) then
            off = abs(w(i,j)); p = i; q = j
          end if
        end do
      end do
      if (off <= tol .or. n == 1) exit
      app = w(p,p); aqq = w(q,q); apq = w(p,q)
      phi = 0.5_dp * atan2(2.0_dp*apq, aqq-app)
      c = cos(phi); s = sin(phi)
      do k = 1, n
        if (k /= p .and. k /= q) then
          t1 = w(k,p); t2 = w(k,q)
          w(k,p) = c*t1 - s*t2; w(p,k) = w(k,p)
          w(k,q) = s*t1 + c*t2; w(q,k) = w(k,q)
        end if
      end do
      w(p,p) = c*c*app - 2.0_dp*s*c*apq + s*s*aqq
      w(q,q) = s*s*app + 2.0_dp*s*c*apq + c*c*aqq
      w(p,q) = 0.0_dp; w(q,p) = 0.0_dp
      do k = 1, n
        t1 = evec(k,p); t2 = evec(k,q)
        evec(k,p) = c*t1 - s*t2
        evec(k,q) = s*t1 + c*t2
      end do
    end do
    do i = 1, n; eval(i) = w(i,i); end do
    do i = 1, n - 1
      do j = i + 1, n
        if (eval(j) < eval(i)) then
          tmp = eval(i); eval(i) = eval(j); eval(j) = tmp
          do k = 1, n
            tmp = evec(k,i); evec(k,i) = evec(k,j); evec(k,j) = tmp
          end do
        end if
      end do
    end do
  end subroutine symmetric_eigen

  subroutine pseudo_inverse_sym(a, ainv, rank)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: ainv(:,:)
    integer, intent(out) :: rank
    real(dp), allocatable :: eval(:), evec(:,:), d(:,:)
    real(dp) :: tol
    integer :: n, i
    n = size(a,1)
    allocate(eval(n), evec(n,n), d(n,n)); d = 0.0_dp
    call symmetric_eigen(a, eval, evec)
    tol = max(100.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(eval))), tiny(1.0_dp))
    rank = 0
    do i = 1, n
      if (eval(i) > tol) then
        d(i,i) = 1.0_dp / eval(i); rank = rank + 1
      end if
    end do
    ainv = matmul(evec, matmul(d, transpose(evec)))
  end subroutine pseudo_inverse_sym

  integer function matrix_rank(a) result(rank)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable :: ata(:,:), eval(:), evec(:,:)
    real(dp) :: tol
    integer :: n
    n = size(a,2)
    allocate(ata(n,n), eval(n), evec(n,n))
    ata = matmul(transpose(a), a)
    call symmetric_eigen(ata, eval, evec)
    tol = 100.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(eval)))
    rank = count(eval > tol)
  end function matrix_rank

  subroutine covariance_from_jacobian(j, ssr, dof, cov, rank)
    real(dp), intent(in) :: j(:,:), ssr
    integer, intent(in) :: dof
    real(dp), intent(out) :: cov(:,:)
    integer, intent(out) :: rank
    real(dp), allocatable :: a(:,:), ai(:,:)
    integer :: n
    n = size(j,2)
    allocate(a(n,n), ai(n,n))
    a = matmul(transpose(j),j)
    call pseudo_inverse_sym(a, ai, rank)
    if (dof > 0) then
      cov = ai * ssr / real(dof,dp)
    else
      cov = ai
    end if
  end subroutine covariance_from_jacobian

  real(dp) function determinant_jtj(j) result(det)
    real(dp), intent(in) :: j(:,:)
    real(dp), allocatable :: a(:,:), l(:,:)
    logical :: ok
    integer :: n, i
    n = size(j,2); allocate(a(n,n), l(n,n))
    a = matmul(transpose(j),j)
    call cholesky_factor(a,l,ok)
    if (.not. ok) then
      det = 0.0_dp
    else
      det = 1.0_dp
      do i = 1, n; det = det * l(i,i) * l(i,i); end do
    end if
  end function determinant_jtj

end module gslnls_linalg
