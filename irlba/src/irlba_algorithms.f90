! Computational translation of the R package irlba 2.3.7.
! Upstream package and native core: GPL-3 / GPL-3-or-later.
! See ../LICENSE and ../UPSTREAM.md for provenance and copyright details.
module irlba_algorithms
  use irlba_kinds, only : dp
  use irlba_sparse, only : csc_matrix
  use irlba_operator, only : dense_operator, csc_operator, linear_operator
  use irlba_core, only : irlba_result, irlba_control, irlb_operator, full_svd_operator
  use irlba_linalg, only : svd_complex, thin_qr, fill_normal, vec_norm2, sort_real_ascending
  implicit none
  private

  type, public :: complex_svd_result
    real(dp), allocatable :: d(:)
    complex(dp), allocatable :: u(:, :)
    complex(dp), allocatable :: v(:, :)
    integer :: info = 0
  end type complex_svd_result

  type, public :: eigen_result
    real(dp), allocatable :: values(:)
    real(dp), allocatable :: vectors(:, :)
    integer :: info = 0
  end type eigen_result

  type, public :: pca_result
    real(dp), allocatable :: sdev(:)
    real(dp), allocatable :: rotation(:, :)
    real(dp), allocatable :: scores(:, :)
    real(dp), allocatable :: center(:)
    real(dp), allocatable :: scale(:)
    real(dp), allocatable :: proportion(:)
    real(dp), allocatable :: cumulative(:)
    real(dp) :: total_variance = 0.0_dp
    integer :: info = 0
  end type pca_result

  type, public :: ssvd_result
    real(dp), allocatable :: u(:, :)
    real(dp), allocatable :: v(:, :)
    real(dp), allocatable :: d(:, :)
    real(dp), allocatable :: lambda(:)
    integer, allocatable :: nonzero(:)
    integer :: iter = 0
    integer :: info = 0
  end type ssvd_result

  type, public :: svdr_result
    real(dp), allocatable :: d(:)
    real(dp), allocatable :: u(:, :)
    real(dp), allocatable :: v(:, :)
    real(dp), allocatable :: q(:, :)
    integer :: iter = 0
    integer :: mprod = 0
    integer :: info = 0
  end type svdr_result

  interface irlba
    module procedure irlba_dense
    module procedure irlba_sparse_matrix
  end interface irlba

  interface partial_eigen
    module procedure partial_eigen_dense
    module procedure partial_eigen_sparse
  end interface partial_eigen

  interface svdr
    module procedure svdr_dense
    module procedure svdr_sparse_matrix
  end interface svdr

  public :: irlba, irlba_dense, irlba_sparse_matrix, irlba_complex
  public :: partial_eigen, partial_eigen_dense, partial_eigen_sparse
  public :: prcomp_irlba, ssvd, svdr, svdr_dense, svdr_sparse_matrix

contains

  subroutine setup_dense_operator(a, op)
    real(dp), target, intent(in) :: a(:, :)
    type(dense_operator), intent(out) :: op
    op%nrow = size(a, 1)
    op%ncol = size(a, 2)
    op%a => a
  end subroutine setup_dense_operator

  subroutine setup_sparse_operator(a, op)
    type(csc_matrix), target, intent(in) :: a
    type(csc_operator), intent(out) :: op
    op%nrow = a%nrow
    op%ncol = a%ncol
    op%a => a
  end subroutine setup_sparse_operator

  function irlba_dense(a, nv, nu, control, start_v, scale, center, shift, restart_from, smallest) result(ans)
    real(dp), target, intent(in) :: a(:, :)
    integer, intent(in) :: nv
    integer, intent(in), optional :: nu
    type(irlba_control), intent(in), optional :: control
    real(dp), intent(in), optional :: start_v(:), scale(:), center(:), shift
    type(irlba_result), intent(in), optional :: restart_from
    logical, intent(in), optional :: smallest
    type(irlba_result) :: ans
    type(dense_operator) :: op
    integer :: nleft
    logical :: small

    call setup_dense_operator(a, op)
    nleft = nv
    if (present(nu)) nleft = nu
    small = .false.
    if (present(smallest)) small = smallest
    if (small .or. min(size(a, 1), size(a, 2)) < 6) then
      call full_svd_operator(op, nleft, nv, ans, scale, center, shift, small)
    else
      call irlb_operator(op, nleft, nv, ans, control, start_v, scale, center, shift, restart_from)
    end if
  end function irlba_dense

  function irlba_sparse_matrix(a, nv, nu, control, start_v, scale, center, shift, restart_from, smallest) result(ans)
    type(csc_matrix), target, intent(in) :: a
    integer, intent(in) :: nv
    integer, intent(in), optional :: nu
    type(irlba_control), intent(in), optional :: control
    real(dp), intent(in), optional :: start_v(:), scale(:), center(:), shift
    type(irlba_result), intent(in), optional :: restart_from
    logical, intent(in), optional :: smallest
    type(irlba_result) :: ans
    type(csc_operator) :: op
    integer :: nleft
    logical :: small

    call setup_sparse_operator(a, op)
    nleft = nv
    if (present(nu)) nleft = nu
    small = .false.
    if (present(smallest)) small = smallest
    if (small .or. min(a%nrow, a%ncol) < 6) then
      call full_svd_operator(op, nleft, nv, ans, scale, center, shift, small)
    else
      call irlb_operator(op, nleft, nv, ans, control, start_v, scale, center, shift, restart_from)
    end if
  end function irlba_sparse_matrix

  function irlba_complex(a, nv, nu, smallest) result(ans)
    complex(dp), intent(in) :: a(:, :)
    integer, intent(in) :: nv
    integer, intent(in), optional :: nu
    logical, intent(in), optional :: smallest
    type(complex_svd_result) :: ans
    real(dp), allocatable :: s(:)
    complex(dp), allocatable :: u(:, :), v(:, :)
    integer :: nleft, k, r, i1, ierr
    logical :: small

    nleft = nv
    if (present(nu)) nleft = nu
    k = max(nleft, nv)
    call svd_complex(a, s, u, v, ierr)
    ans%info = ierr
    if (ierr /= 0) return
    r = size(s)
    small = .false.
    if (present(smallest)) small = smallest
    if (small) then
      i1 = r - k + 1
    else
      i1 = 1
    end if
    allocate(ans%d(k), ans%u(size(a, 1), nleft), ans%v(size(a, 2), nv))
    ans%d = s(i1:i1 + k - 1)
    ans%u = u(:, i1:i1 + nleft - 1)
    ans%v = v(:, i1:i1 + nv - 1)
  end function irlba_complex

  function partial_eigen_dense(x, n, control) result(ans)
    real(dp), target, intent(in) :: x(:, :)
    integer, intent(in) :: n
    type(irlba_control), intent(in), optional :: control
    type(eigen_result) :: ans
    type(irlba_result) :: s1, s2
    integer :: j, first_negative
    real(dp) :: sh

    if (size(x, 1) /= size(x, 2)) error stop "partial_eigen: symmetric matrix must be square"
    s1 = irlba_dense(x, n, n, control=control)
    if (s1%info /= 0 .and. s1%info /= -2) then
      ans%info = s1%info
      return
    end if
    first_negative = 0
    do j = 1, n
      if (s1%u(1, j) * s1%v(1, j) < 0.0_dp) then
        first_negative = j
        exit
      end if
    end do
    allocate(ans%values(n), ans%vectors(size(x, 1), n))
    if (first_negative == 0) then
      ans%values = s1%d
      ans%vectors = s1%u
      ans%info = s1%info
    else
      sh = s1%d(first_negative)
      s2 = irlba_dense(x, n, n, control=control, shift=sh)
      ans%values = s2%d - sh
      ans%vectors = s2%u
      ans%info = s2%info
    end if
  end function partial_eigen_dense

  function partial_eigen_sparse(x, n, control) result(ans)
    type(csc_matrix), target, intent(in) :: x
    integer, intent(in) :: n
    type(irlba_control), intent(in), optional :: control
    type(eigen_result) :: ans
    type(irlba_result) :: s1, s2
    integer :: j, first_negative
    real(dp) :: sh

    if (x%nrow /= x%ncol) error stop "partial_eigen: symmetric matrix must be square"
    s1 = irlba_sparse_matrix(x, n, n, control=control)
    if (s1%info /= 0 .and. s1%info /= -2) then
      ans%info = s1%info
      return
    end if
    first_negative = 0
    do j = 1, n
      if (s1%u(1, j) * s1%v(1, j) < 0.0_dp) then
        first_negative = j
        exit
      end if
    end do
    allocate(ans%values(n), ans%vectors(x%nrow, n))
    if (first_negative == 0) then
      ans%values = s1%d
      ans%vectors = s1%u
      ans%info = s1%info
    else
      sh = s1%d(first_negative)
      s2 = irlba_sparse_matrix(x, n, n, control=control, shift=sh)
      ans%values = s2%d - sh
      ans%vectors = s2%u
      ans%info = s2%info
    end if
  end function partial_eigen_sparse

  function prcomp_irlba(x, n, do_center, do_scale, center_values, scale_values, retx, control) result(ans)
    real(dp), target, intent(in) :: x(:, :)
    integer, intent(in) :: n
    logical, intent(in), optional :: do_center, do_scale, retx
    real(dp), intent(in), optional :: center_values(:), scale_values(:)
    type(irlba_control), intent(in), optional :: control
    type(pca_result) :: ans
    type(irlba_result) :: s
    real(dp), allocatable :: ctr(:), scl(:)
    real(dp) :: denom
    integer :: j
    logical :: dc, ds, rx

    dc = .true.; ds = .false.; rx = .true.
    if (present(do_center)) dc = do_center
    if (present(do_scale)) ds = do_scale
    if (present(retx)) rx = retx
    allocate(ctr(size(x, 2)), scl(size(x, 2)))
    ctr = 0.0_dp
    scl = 1.0_dp
    if (present(center_values)) then
      if (size(center_values) /= size(x, 2)) error stop "prcomp_irlba: center length mismatch"
      ctr = center_values
      dc = .true.
    else if (dc) then
      ctr = sum(x, dim=1) / real(size(x, 1), dp)
    end if
    denom = real(max(1, size(x, 1) - 1), dp)
    if (present(scale_values)) then
      if (size(scale_values) /= size(x, 2)) error stop "prcomp_irlba: scale length mismatch"
      scl = scale_values
      ds = .true.
    else if (ds) then
      do j = 1, size(x, 2)
        if (dc) then
          scl(j) = sqrt(sum((x(:, j) - ctr(j)) ** 2) / denom)
        else
          scl(j) = sqrt(sum(x(:, j) ** 2) / denom)
        end if
      end do
    end if

    if (dc .and. ds) then
      s = irlba_dense(x, n, n, control=control, center=ctr, scale=scl)
    else if (dc) then
      s = irlba_dense(x, n, n, control=control, center=ctr)
    else if (ds) then
      s = irlba_dense(x, n, n, control=control, scale=scl)
    else
      s = irlba_dense(x, n, n, control=control)
    end if
    ans%info = s%info
    allocate(ans%sdev(n), ans%rotation(size(x, 2), n), ans%center(size(x, 2)), ans%scale(size(x, 2)))
    ans%sdev = s%d / sqrt(denom)
    ans%rotation = s%v
    ans%center = ctr
    ans%scale = scl
    if (rx) then
      allocate(ans%scores(size(x, 1), n))
      do j = 1, n
        ans%scores(:, j) = s%u(:, j) * s%d(j)
      end do
    end if
    if (ds) then
      ans%total_variance = real(size(x, 2), dp)
      if (.not. dc) then
        ans%total_variance = 0.0_dp
        do j = 1, size(x, 2)
          ans%total_variance = ans%total_variance + sum((x(:, j) / scl(j)) ** 2) / denom
        end do
      end if
    else
      ans%total_variance = 0.0_dp
      do j = 1, size(x, 2)
        if (dc) then
          ans%total_variance = ans%total_variance + sum((x(:, j) - ctr(j)) ** 2) / denom
        else
          ans%total_variance = ans%total_variance + sum(x(:, j) ** 2) / denom
        end if
      end do
    end if
    allocate(ans%proportion(n), ans%cumulative(n))
    ans%proportion = ans%sdev ** 2 / ans%total_variance
    ans%cumulative(1) = ans%proportion(1)
    do j = 2, n
      ans%cumulative(j) = ans%cumulative(j - 1) + ans%proportion(j)
    end do
  end function prcomp_irlba

  function ssvd(x, k, nonzero, maxit, tol, center, scale, alpha, control) result(ans)
    real(dp), target, intent(in) :: x(:, :)
    integer, intent(in) :: k
    integer, intent(in) :: nonzero(:)
    integer, intent(in), optional :: maxit
    real(dp), intent(in), optional :: tol, center(:), scale(:), alpha
    type(irlba_control), intent(in), optional :: control
    type(ssvd_result) :: ans
    type(irlba_result) :: s
    real(dp), allocatable :: vold(:, :), vnew(:, :), xsv(:, :), q(:, :), a(:), z(:), dwork(:, :)
    real(dp) :: al, tolerance, delta, lam, sn
    integer :: itmax, iter, j, p, ierr, first_old, first_new

    if (size(nonzero) /= 1 .and. size(nonzero) /= k) error stop "ssvd: nonzero must have length 1 or k"
    al = 0.0_dp
    if (present(alpha)) al = alpha
    if (al < 0.0_dp .or. al >= 1.0_dp) error stop "ssvd: alpha must satisfy 0 <= alpha < 1"
    itmax = 500
    if (present(maxit)) itmax = maxit
    tolerance = 1.0e-3_dp
    if (present(tol)) tolerance = tol

    if (present(center) .and. present(scale)) then
      s = irlba_dense(x, k, k, control=control, center=center, scale=scale)
    else if (present(center)) then
      s = irlba_dense(x, k, k, control=control, center=center)
    else if (present(scale)) then
      s = irlba_dense(x, k, k, control=control, scale=scale)
    else
      s = irlba_dense(x, k, k, control=control)
    end if
    allocate(vnew(size(x, 2), k), vold(size(x, 1), k), xsv(size(x, 1), k))
    allocate(ans%lambda(k), ans%nonzero(k))
    vnew = s%v
    do j = 1, k
      vnew(:, j) = vnew(:, j) * s%d(j)
      ans%nonzero(j) = nonzero(merge(1, j, size(nonzero) == 1))
      if (ans%nonzero(j) < 1 .or. ans%nonzero(j) >= size(x, 2)) error stop "ssvd: invalid nonzero target"
    end do
    ans%u = s%u
    delta = huge(1.0_dp)
    iter = 0
    allocate(a(size(x, 2)), z(size(x, 2)))
    do while (delta > tolerance .and. iter < itmax)
      vold = ans%u
      do j = 1, k
        vnew(:, j) = matmul(transpose(x), ans%u(:, j))
        if (present(center)) vnew(:, j) = vnew(:, j) - sum(ans%u(:, j)) * center
        if (present(scale)) vnew(:, j) = vnew(:, j) / scale
        a = abs(vnew(:, j))
        z = a
        call sort_real_ascending(z)
        p = size(x, 2) - ans%nonzero(j)
        lam = (1.0_dp - al) * z(p) + al * z(p + 1)
        ans%lambda(j) = lam
        where (a > lam)
          vnew(:, j) = sign(1.0_dp, vnew(:, j)) * (a - lam)
        elsewhere
          vnew(:, j) = 0.0_dp
        end where
        if (present(scale)) vnew(:, j) = vnew(:, j) / scale
      end do
      xsv = matmul(x, vnew)
      if (present(center)) then
        do j = 1, k
          xsv(:, j) = xsv(:, j) - dot_product(center, vnew(:, j))
        end do
      end if
      call thin_qr(xsv, q, ierr)
      if (ierr /= 0) then
        ans%info = ierr
        return
      end if
      ans%u = q(:, 1:k)
      do j = 1, k
        first_old = first_nonzero_index(xsv(:, j))
        first_new = first_nonzero_index(ans%u(:, j))
        if (first_old > 0 .and. first_new > 0) then
          if (xsv(first_old, j) * ans%u(first_new, j) < 0.0_dp) then
            ans%u(:, j) = -ans%u(:, j)
          end if
        end if
      end do
      delta = 0.0_dp
      do j = 1, k
        delta = max(delta, 1.0_dp - abs(dot_product(vold(:, j), ans%u(:, j))))
      end do
      iter = iter + 1
    end do

    do j = 1, k
      sn = vec_norm2(vnew(:, j))
      if (sn > 0.0_dp) vnew(:, j) = vnew(:, j) / sn
    end do
    allocate(ans%v(size(x, 2), k), ans%d(k, k))
    ans%v = vnew
    allocate(dwork(size(x, 2), k))
    dwork = ans%v
    if (present(scale)) then
      do j = 1, k
        dwork(:, j) = dwork(:, j) / scale
      end do
    end if
    xsv = matmul(x, dwork)
    if (present(center)) then
      do j = 1, k
        xsv(:, j) = xsv(:, j) - dot_product(center, dwork(:, j))
      end do
    end if
    ans%d = matmul(transpose(ans%u), xsv)
    ans%iter = iter
    ans%info = merge(0, -2, iter < itmax)
  end function ssvd

  integer function first_nonzero_index(x) result(idx)
    real(dp), intent(in) :: x(:)
    integer :: i
    idx = 0
    do i = 1, size(x)
      if (abs(x(i)) > tiny(1.0_dp)) then
        idx = i
        return
      end if
    end do
  end function first_nonzero_index

  function svdr_dense(x, k, tol, it, extra, center, q_initial, return_q) result(ans)
    real(dp), intent(in) :: x(:, :)
    integer, intent(in) :: k
    real(dp), intent(in), optional :: tol
    integer, intent(in), optional :: it, extra
    real(dp), intent(in), optional :: center(:), q_initial(:, :)
    logical, intent(in), optional :: return_q
    type(svdr_result) :: ans
    real(dp), allocatable :: q(:, :), q2(:, :), b(:, :), s(:), ub(:, :), vb(:, :), dprev(:), ones(:)
    real(dp) :: tolerance, eps2, rel
    integer :: nit, nextra, nsub, j, ierr, r, ii
    logical :: rq

    tolerance = 1.0e-5_dp
    if (present(tol)) tolerance = tol
    nit = 100
    if (present(it)) nit = it
    nextra = min(10, min(size(x, 1), size(x, 2)) - k)
    if (present(extra)) nextra = extra
    nsub = min(size(x, 2), k + max(0, nextra))
    rq = .false.
    if (present(return_q)) rq = return_q
    eps2 = epsilon(1.0_dp) ** (4.0_dp / 5.0_dp)
    allocate(q(size(x, 2), nsub), dprev(k), ones(size(x, 1)))
    ones = 1.0_dp
    dprev = 0.0_dp
    if (present(q_initial)) then
      if (size(q_initial, 1) /= size(x, 2) .or. size(q_initial, 2) < nsub) error stop "svdr: Q dimension mismatch"
      q = q_initial(:, 1:nsub)
    else
      do ii = 1, size(q, 2)
        call fill_normal(q(:, ii))
      end do
    end if

    do j = 1, nit
      q2 = matmul(x, q)
      if (present(center)) q2 = q2 - spread(matmul(center, q), 1, size(x, 1))
      call thin_qr(q2, q, ierr)
      if (ierr /= 0) then
        ans%info = ierr
        return
      end if
      if (present(center)) then
        b = matmul(transpose(q), x) - spread(matmul(transpose(q), ones), 2, size(x, 2)) * &
            spread(center, 1, size(q, 2))
        call thin_qr(transpose(b), q2, ierr)
        q = q2
      else
        b = matmul(transpose(x), q)
        call thin_qr(b, q2, ierr)
        q = q2
      end if
      call svd_real_local_for_svdr(b, s, ierr)
      if (ierr /= 0) then
        ans%info = ierr
        return
      end if
      r = min(k, size(s))
      if (j > 1) then
        rel = 0.0_dp
        do ii = 1, r
          if (s(ii) > eps2 .and. abs(dprev(ii)) > eps2) then
            rel = max(rel, abs((s(ii) - dprev(ii)) / dprev(ii)))
          end if
        end do
        if (rel < tolerance) exit
      end if
      dprev(1:r) = s(1:r)
    end do

    if (rq) then
      allocate(ans%q(size(q, 1), size(q, 2)))
      ans%q = q
    end if
    q2 = matmul(x, q)
    if (present(center)) q2 = q2 - spread(matmul(center, q), 1, size(x, 1))
    call thin_qr(q2, q, ierr)
    if (present(center)) then
      b = matmul(transpose(q), x) - spread(matmul(transpose(q), ones), 2, size(x, 2)) * &
          spread(center, 1, size(q, 2))
    else
      b = matmul(transpose(q), x)
    end if
    call svd_small_matrix(b, s, ub, vb, ierr)
    if (ierr /= 0) then
      ans%info = ierr
      return
    end if
    allocate(ans%d(k), ans%u(size(x, 1), k), ans%v(size(x, 2), k))
    ans%d = s(1:k)
    ans%u = matmul(q, ub(:, 1:k))
    ans%v = vb(:, 1:k)
    ans%iter = j
    ans%mprod = 2 * j + 1
    ans%info = 0
  end function svdr_dense

  function svdr_sparse_matrix(x, k, tol, it, extra, center, q_initial, return_q) result(ans)
    type(csc_matrix), intent(in) :: x
    integer, intent(in) :: k
    real(dp), intent(in), optional :: tol
    integer, intent(in), optional :: it, extra
    real(dp), intent(in), optional :: center(:), q_initial(:, :)
    logical, intent(in), optional :: return_q
    type(svdr_result) :: ans
    real(dp), allocatable :: q(:, :), q2(:, :), b(:, :), xtq(:, :), s(:), ub(:, :), vb(:, :)
    real(dp), allocatable :: dprev(:), ones(:), qtx(:, :)
    real(dp) :: tolerance, eps2, rel
    integer :: nit, nextra, nsub, j, ierr, r, ii
    logical :: rq

    tolerance = 1.0e-5_dp
    if (present(tol)) tolerance = tol
    nit = 100
    if (present(it)) nit = it
    nextra = min(10, min(x%nrow, x%ncol) - k)
    if (present(extra)) nextra = extra
    nsub = min(x%ncol, k + max(0, nextra))
    rq = .false.
    if (present(return_q)) rq = return_q
    eps2 = epsilon(1.0_dp) ** (4.0_dp / 5.0_dp)
    allocate(q(x%ncol, nsub), dprev(k), ones(x%nrow))
    ones = 1.0_dp
    dprev = 0.0_dp
    if (present(q_initial)) then
      if (size(q_initial, 1) /= x%ncol .or. size(q_initial, 2) < nsub) error stop "svdr: Q dimension mismatch"
      q = q_initial(:, 1:nsub)
    else
      do ii = 1, size(q, 2)
        call fill_normal(q(:, ii))
      end do
    end if

    do j = 1, nit
      allocate(q2(x%nrow, nsub))
      call x%matmat(q, q2)
      if (present(center)) q2 = q2 - spread(matmul(center, q), 1, x%nrow)
      call thin_qr(q2, q, ierr)
      deallocate(q2)
      if (ierr /= 0) then
        ans%info = ierr
        return
      end if
      allocate(xtq(x%ncol, size(q, 2)))
      call x%tmatmat(q, xtq)
      if (present(center)) then
        allocate(qtx(size(q, 2), x%ncol))
        qtx = transpose(xtq) - spread(matmul(transpose(q), ones), 2, x%ncol) * &
              spread(center, 1, size(q, 2))
        b = qtx
        call thin_qr(transpose(b), q2, ierr)
        q = q2
        deallocate(qtx)
      else
        b = xtq
        call thin_qr(b, q2, ierr)
        q = q2
      end if
      deallocate(xtq, q2)
      call svd_real_local_for_svdr(b, s, ierr)
      if (ierr /= 0) then
        ans%info = ierr
        return
      end if
      r = min(k, size(s))
      if (j > 1) then
        rel = 0.0_dp
        do ii = 1, r
          if (s(ii) > eps2 .and. abs(dprev(ii)) > eps2) then
            rel = max(rel, abs((s(ii) - dprev(ii)) / dprev(ii)))
          end if
        end do
        if (rel < tolerance) exit
      end if
      dprev(1:r) = s(1:r)
      deallocate(s)
    end do

    if (rq) then
      allocate(ans%q(size(q, 1), size(q, 2)))
      ans%q = q
    end if
    allocate(q2(x%nrow, size(q, 2)))
    call x%matmat(q, q2)
    if (present(center)) q2 = q2 - spread(matmul(center, q), 1, x%nrow)
    call thin_qr(q2, q, ierr)
    deallocate(q2)
    allocate(xtq(x%ncol, size(q, 2)))
    call x%tmatmat(q, xtq)
    b = transpose(xtq)
    if (present(center)) b = b - spread(matmul(transpose(q), ones), 2, x%ncol) * &
                                      spread(center, 1, size(q, 2))
    call svd_small_matrix(b, s, ub, vb, ierr)
    if (ierr /= 0) then
      ans%info = ierr
      return
    end if
    allocate(ans%d(k), ans%u(x%nrow, k), ans%v(x%ncol, k))
    ans%d = s(1:k)
    ans%u = matmul(q, ub(:, 1:k))
    ans%v = vb(:, 1:k)
    ans%iter = j
    ans%mprod = 2 * j + 1
    ans%info = 0
  end function svdr_sparse_matrix

  subroutine svd_real_local_for_svdr(a, s, info)
    use irlba_linalg, only : svd_real
    real(dp), intent(in) :: a(:, :)
    real(dp), allocatable, intent(out) :: s(:)
    integer, intent(out) :: info
    real(dp), allocatable :: u(:, :), v(:, :)
    call svd_real(a, s, u, v, info)
  end subroutine svd_real_local_for_svdr

  subroutine svd_small_matrix(a, s, u, v, info)
    use irlba_linalg, only : svd_real
    real(dp), intent(in) :: a(:, :)
    real(dp), allocatable, intent(out) :: s(:), u(:, :), v(:, :)
    integer, intent(out) :: info
    call svd_real(a, s, u, v, info)
  end subroutine svd_small_matrix

end module irlba_algorithms
