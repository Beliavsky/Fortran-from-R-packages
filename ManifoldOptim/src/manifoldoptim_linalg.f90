! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of the ManifoldOptim/ROPTLIB computational interface.
! See LICENSE and UPSTREAM_PROVENANCE.md for upstream copyright and provenance.
module manifoldoptim_linalg
  use manifoldoptim_kinds, only : dp
  implicit none
  private
  public :: vecnorm, symmetrize, skew_part, mgs_orthonormalize
  public :: symmetric_eigen, symmetric_matrix_function, matrix_inverse
  public :: is_spd, eye_matrix, trace_matrix
  public :: orthogonal_complement, matrix_exponential, solve_linear_system

contains

  pure real(dp) function vecnorm(x) result(v)
    real(dp), intent(in) :: x(:)
    v = sqrt(max(0.0_dp, dot_product(x, x)))
  end function vecnorm

  pure function eye_matrix(n) result(a)
    integer, intent(in) :: n
    real(dp) :: a(n,n)
    integer :: i
    a = 0.0_dp
    do i = 1, n
      a(i,i) = 1.0_dp
    end do
  end function eye_matrix

  pure real(dp) function trace_matrix(a) result(t)
    real(dp), intent(in) :: a(:,:)
    integer :: i
    t = 0.0_dp
    do i = 1, min(size(a,1), size(a,2))
      t = t + a(i,i)
    end do
  end function trace_matrix

  pure subroutine symmetrize(a, s)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: s(:,:)
    s = 0.5_dp * (a + transpose(a))
  end subroutine symmetrize

  pure subroutine skew_part(a, s)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: s(:,:)
    s = 0.5_dp * (a - transpose(a))
  end subroutine skew_part

  subroutine mgs_orthonormalize(a, q, ok)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: q(:,:)
    logical, intent(out) :: ok
    integer :: i, j, n, p
    real(dp) :: nr
    n = size(a,1)
    p = size(a,2)
    q = a
    ok = .true.
    do j = 1, p
      do i = 1, j - 1
        q(:,j) = q(:,j) - dot_product(q(:,i), q(:,j)) * q(:,i)
      end do
      ! one reorthogonalization pass is important for nearly dependent columns
      do i = 1, j - 1
        q(:,j) = q(:,j) - dot_product(q(:,i), q(:,j)) * q(:,i)
      end do
      nr = sqrt(max(0.0_dp, dot_product(q(:,j), q(:,j))))
      if (nr <= 100.0_dp * epsilon(1.0_dp)) then
        ok = .false.
        q(:,j) = 0.0_dp
        if (j <= n) q(j,j) = 1.0_dp
        do i = 1, j - 1
          q(:,j) = q(:,j) - dot_product(q(:,i), q(:,j)) * q(:,i)
        end do
        nr = sqrt(max(0.0_dp, dot_product(q(:,j), q(:,j))))
        if (nr <= epsilon(1.0_dp)) return
      end if
      q(:,j) = q(:,j) / nr
    end do
  end subroutine mgs_orthonormalize

  subroutine symmetric_eigen(a, eigval, eigvec, info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: eigval(:), eigvec(:,:)
    integer, intent(out) :: info
    real(dp), allocatable :: b(:,:)
    real(dp) :: app, aqq, apq, tau, t, c, s, off, tol, bip, biq, vip, viq
    integer :: n, i, p, q, sweep, max_sweep
    n = size(a,1)
    if (size(a,2) /= n .or. size(eigval) /= n .or. any(shape(eigvec) /= [n,n])) then
      info = -1
      return
    end if
    allocate(b(n,n))
    b = 0.5_dp * (a + transpose(a))
    eigvec = eye_matrix(n)
    tol = 50.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(b)))
    max_sweep = max(30, 8*n*n)
    info = 1
    do sweep = 1, max_sweep
      off = 0.0_dp
      p = 1
      q = min(2,n)
      do i = 1, n - 1
        if (maxval(abs(b(i,i+1:n))) > off) then
          q = i + maxloc(abs(b(i,i+1:n)), dim=1)
          p = i
          off = abs(b(p,q))
        end if
      end do
      if (n <= 1 .or. off <= tol) then
        info = 0
        exit
      end if
      app = b(p,p)
      aqq = b(q,q)
      apq = b(p,q)
      if (abs(apq) <= tiny(1.0_dp)) cycle
      tau = (aqq - app) / (2.0_dp * apq)
      if (tau >= 0.0_dp) then
        t = 1.0_dp / (tau + sqrt(1.0_dp + tau*tau))
      else
        t = -1.0_dp / (-tau + sqrt(1.0_dp + tau*tau))
      end if
      c = 1.0_dp / sqrt(1.0_dp + t*t)
      s = t*c
      do i = 1, n
        if (i /= p .and. i /= q) then
          bip = b(i,p)
          biq = b(i,q)
          b(i,p) = c*bip - s*biq
          b(p,i) = b(i,p)
          b(i,q) = s*bip + c*biq
          b(q,i) = b(i,q)
        end if
      end do
      b(p,p) = c*c*app - 2.0_dp*s*c*apq + s*s*aqq
      b(q,q) = s*s*app + 2.0_dp*s*c*apq + c*c*aqq
      b(p,q) = 0.0_dp
      b(q,p) = 0.0_dp
      do i = 1, n
        vip = eigvec(i,p)
        viq = eigvec(i,q)
        eigvec(i,p) = c*vip - s*viq
        eigvec(i,q) = s*vip + c*viq
      end do
    end do
    do i = 1, n
      eigval(i) = b(i,i)
    end do
    call sort_eigenpairs(eigval, eigvec)
  end subroutine symmetric_eigen

  subroutine sort_eigenpairs(d, q)
    real(dp), intent(inout) :: d(:), q(:,:)
    integer :: i, j, k, n
    real(dp) :: td
    real(dp), allocatable :: col(:)
    n = size(d)
    allocate(col(size(q,1)))
    do i = 1, n - 1
      k = i
      do j = i + 1, n
        if (d(j) < d(k)) k = j
      end do
      if (k /= i) then
        td = d(i)
        d(i) = d(k)
        d(k) = td
        col = q(:,i)
        q(:,i) = q(:,k)
        q(:,k) = col
      end if
    end do
  end subroutine sort_eigenpairs

  subroutine symmetric_matrix_function(a, power, mode, result, info, floor_eig)
    real(dp), intent(in) :: a(:,:), power
    character(len=*), intent(in) :: mode
    real(dp), intent(out) :: result(:,:)
    integer, intent(out) :: info
    real(dp), intent(in), optional :: floor_eig
    integer :: n, i
    real(dp) :: fl
    real(dp), allocatable :: d(:), q(:,:), w(:,:)
    n = size(a,1)
    allocate(d(n), q(n,n), w(n,n))
    call symmetric_eigen(a, d, q, info)
    if (info /= 0) return
    fl = 100.0_dp * epsilon(1.0_dp)
    if (present(floor_eig)) fl = floor_eig
    select case (trim(mode))
    case ('power')
      do i = 1, n
        d(i) = max(d(i), fl)**power
      end do
    case ('exp')
      do i = 1, n
        d(i) = exp(d(i))
      end do
    case ('log')
      do i = 1, n
        d(i) = log(max(d(i), fl))
      end do
    case default
      info = -2
      return
    end select
    w = 0.0_dp
    do i = 1, n
      w(:,i) = q(:,i) * d(i)
    end do
    result = matmul(w, transpose(q))
    result = 0.5_dp * (result + transpose(result))
  end subroutine symmetric_matrix_function

  subroutine matrix_inverse(a, ainv, info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: ainv(:,:)
    integer, intent(out) :: info
    integer :: n, i, j, k, piv
    real(dp) :: mx, factor, tmp
    real(dp), allocatable :: aug(:,:)
    n = size(a,1)
    if (size(a,2) /= n .or. any(shape(ainv) /= [n,n])) then
      info = -1
      return
    end if
    allocate(aug(n,2*n))
    aug(:,1:n) = a
    aug(:,n+1:2*n) = eye_matrix(n)
    info = 0
    do k = 1, n
      piv = k
      mx = abs(aug(k,k))
      do i = k + 1, n
        if (abs(aug(i,k)) > mx) then
          piv = i
          mx = abs(aug(i,k))
        end if
      end do
      if (mx <= 100.0_dp*epsilon(1.0_dp)) then
        info = k
        return
      end if
      if (piv /= k) then
        do j = 1, 2*n
          tmp = aug(k,j)
          aug(k,j) = aug(piv,j)
          aug(piv,j) = tmp
        end do
      end if
      tmp = aug(k,k)
      aug(k,:) = aug(k,:) / tmp
      do i = 1, n
        if (i /= k) then
          factor = aug(i,k)
          aug(i,:) = aug(i,:) - factor*aug(k,:)
        end if
      end do
    end do
    ainv = aug(:,n+1:2*n)
  end subroutine matrix_inverse

  subroutine orthogonal_complement(q, qperp, ok)
    real(dp), intent(in) :: q(:,:)
    real(dp), intent(out) :: qperp(:,:)
    logical, intent(out) :: ok
    integer :: n, p, nc, i, j, k
    real(dp) :: nr
    real(dp), allocatable :: v(:)

    n = size(q,1)
    p = size(q,2)
    nc = n - p
    if (nc < 0 .or. any(shape(qperp) /= [n,nc])) then
      ok = .false.
      return
    end if
    if (nc == 0) then
      ok = .true.
      return
    end if

    allocate(v(n))
    qperp = 0.0_dp
    k = 0
    do i = 1, n
      v = 0.0_dp
      v(i) = 1.0_dp
      do j = 1, p
        v = v - dot_product(q(:,j),v) * q(:,j)
      end do
      do j = 1, k
        v = v - dot_product(qperp(:,j),v) * qperp(:,j)
      end do
      ! Reorthogonalize to make the complement deterministic and stable.
      do j = 1, p
        v = v - dot_product(q(:,j),v) * q(:,j)
      end do
      do j = 1, k
        v = v - dot_product(qperp(:,j),v) * qperp(:,j)
      end do
      nr = vecnorm(v)
      if (nr > 100.0_dp * epsilon(1.0_dp)) then
        k = k + 1
        qperp(:,k) = v / nr
        if (k == nc) exit
      end if
    end do
    ok = k == nc
  end subroutine orthogonal_complement

  subroutine matrix_exponential(a, ea, info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: ea(:,:)
    integer, intent(out) :: info
    integer :: n, i, ns, k
    real(dp) :: anorm, termnorm
    real(dp), allocatable :: b(:,:), term(:,:), esum(:,:)

    n = size(a,1)
    if (size(a,2) /= n .or. any(shape(ea) /= [n,n])) then
      info = -1
      return
    end if
    if (n == 0) then
      info = 0
      return
    end if

    anorm = maxval(sum(abs(a),dim=1))
    ns = 0
    if (anorm > 0.5_dp) ns = max(0, ceiling(log(anorm/0.5_dp)/log(2.0_dp)))
    allocate(b(n,n), term(n,n), esum(n,n))
    b = a / (2.0_dp**ns)
    esum = eye_matrix(n)
    term = eye_matrix(n)
    do k = 1, 80
      term = matmul(term,b) / real(k,dp)
      esum = esum + term
      termnorm = maxval(abs(term))
      if (termnorm <= 10.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(esum)))) exit
    end do
    ea = esum
    do i = 1, ns
      ea = matmul(ea,ea)
    end do
    info = 0
  end subroutine matrix_exponential

  subroutine solve_linear_system(a, b, x, info)
    real(dp), intent(in) :: a(:,:), b(:)
    real(dp), intent(out) :: x(:)
    integer, intent(out) :: info
    integer :: n, i, j, k, piv
    real(dp) :: mx, factor, tmp
    real(dp), allocatable :: aug(:,:)

    n = size(b)
    if (any(shape(a) /= [n,n]) .or. size(x) /= n) then
      info = -1
      return
    end if
    allocate(aug(n,n+1))
    aug(:,1:n) = a
    aug(:,n+1) = b
    info = 0
    do k = 1, n
      piv = k
      mx = abs(aug(k,k))
      do i = k + 1, n
        if (abs(aug(i,k)) > mx) then
          piv = i
          mx = abs(aug(i,k))
        end if
      end do
      if (mx <= 100.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(a)))) then
        info = k
        x = 0.0_dp
        return
      end if
      if (piv /= k) then
        do j = k, n + 1
          tmp = aug(k,j)
          aug(k,j) = aug(piv,j)
          aug(piv,j) = tmp
        end do
      end if
      tmp = aug(k,k)
      aug(k,k:n+1) = aug(k,k:n+1) / tmp
      do i = 1, n
        if (i /= k) then
          factor = aug(i,k)
          aug(i,k:n+1) = aug(i,k:n+1) - factor*aug(k,k:n+1)
        end if
      end do
    end do
    x = aug(:,n+1)
  end subroutine solve_linear_system

  logical function is_spd(a, tol) result(ok)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(in), optional :: tol
    real(dp), allocatable :: d(:), q(:,:)
    real(dp) :: t
    integer :: n, info
    n = size(a,1)
    allocate(d(n), q(n,n))
    call symmetric_eigen(0.5_dp*(a+transpose(a)), d, q, info)
    t = 1000.0_dp*epsilon(1.0_dp)
    if (present(tol)) t = tol
    ok = info == 0 .and. minval(d) > t
  end function is_spd

end module manifoldoptim_linalg
