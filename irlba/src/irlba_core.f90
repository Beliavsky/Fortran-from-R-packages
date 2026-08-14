! Computational translation of the R package irlba 2.3.7.
! Upstream package and native core: GPL-3 / GPL-3-or-later.
! See ../LICENSE and ../UPSTREAM.md for provenance and copyright details.
module irlba_core
  use irlba_kinds, only : dp
  use irlba_linalg, only : vec_norm2, orthogonalize, fill_normal, svd_real
  use irlba_operator, only : linear_operator
  implicit none
  private

  type, public :: irlba_result
    real(dp), allocatable :: d(:)
    real(dp), allocatable :: u(:, :)
    real(dp), allocatable :: v(:, :)
    integer :: iter = 0
    integer :: mprod = 0
    integer :: info = 0
    logical :: converged = .false.
  end type irlba_result

  type, public :: irlba_control
    integer :: maxit = 1000
    integer :: work = 0
    real(dp) :: tol = 1.0e-5_dp
    real(dp) :: svtol = -1.0_dp
    real(dp) :: invariant_tol = -1.0_dp
    logical :: reorth = .true.
  end type irlba_control

  public :: irlb_operator, full_svd_operator

contains

  subroutine apply_forward(op, v, y, scale, center, shift)
    class(linear_operator), intent(in) :: op
    real(dp), intent(in) :: v(:)
    real(dp), intent(out) :: y(:)
    real(dp), intent(in), optional :: scale(:), center(:), shift
    real(dp), allocatable :: x(:)
    real(dp) :: beta

    allocate(x(size(v)))
    x = v
    if (present(scale)) x = x / scale
    call op%matvec(x, y)
    if (present(shift)) then
      if (op%nrow /= op%ncol) error stop "irlba: shift requires a square matrix"
      y = y + shift * x
    end if
    if (present(center)) then
      beta = dot_product(x, center)
      y = y - beta
    end if
  end subroutine apply_forward

  subroutine apply_transpose(op, w, y, scale, center, shift)
    class(linear_operator), intent(in) :: op
    real(dp), intent(in) :: w(:)
    real(dp), intent(out) :: y(:)
    real(dp), intent(in), optional :: scale(:), center(:), shift
    real(dp) :: beta

    call op%tmatvec(w, y)
    if (present(shift)) then
      if (op%nrow /= op%ncol) error stop "irlba: shift requires a square matrix"
      y = y + shift * w
    end if
    if (present(scale)) y = y / scale
    if (present(center)) then
      beta = sum(w)
      if (present(scale)) then
        y = y - beta * center / scale
      else
        y = y - beta * center
      end if
    end if
  end subroutine apply_transpose

  subroutine full_svd_operator(op, nu, nv, result, scale, center, shift, smallest)
    class(linear_operator), intent(in) :: op
    integer, intent(in) :: nu, nv
    type(irlba_result), intent(out) :: result
    real(dp), intent(in), optional :: scale(:), center(:), shift
    logical, intent(in), optional :: smallest
    real(dp), allocatable :: a(:, :), x(:), y(:), s(:), u(:, :), v(:, :)
    integer :: j, k, r, i1, i2, ierr
    logical :: small

    allocate(a(op%nrow, op%ncol), x(op%ncol), y(op%nrow))
    do j = 1, op%ncol
      x = 0.0_dp
      x(j) = 1.0_dp
      call apply_forward(op, x, y, scale, center, shift)
      a(:, j) = y
    end do
    call svd_real(a, s, u, v, ierr)
    result%info = ierr
    if (ierr /= 0) return
    k = max(nu, nv)
    r = size(s)
    small = .false.
    if (present(smallest)) small = smallest
    if (small) then
      i1 = r - k + 1
      i2 = r
    else
      i1 = 1
      i2 = k
    end if
    allocate(result%d(k), result%u(op%nrow, nu), result%v(op%ncol, nv))
    result%d = s(i1:i2)
    result%u = u(:, i1:i1 + nu - 1)
    result%v = v(:, i1:i1 + nv - 1)
    result%iter = 0
    result%mprod = 2 * op%ncol
    result%converged = .true.
  end subroutine full_svd_operator

  subroutine irlb_operator(op, nu, nv, result, control, start_v, scale, center, shift, restart_from)
    class(linear_operator), intent(in) :: op
    integer, intent(in) :: nu, nv
    type(irlba_result), intent(out) :: result
    type(irlba_control), intent(in), optional :: control
    real(dp), intent(in), optional :: start_v(:), scale(:), center(:), shift
    type(irlba_result), intent(in), optional :: restart_from

    type(irlba_control) :: ctl
    real(dp), allocatable :: v(:, :), w(:, :), b(:, :), f(:), y(:), tmp(:)
    real(dp), allocatable :: bs(:), ub(:, :), vb(:, :), res(:), svratio(:)
    real(dp), allocatable :: vnew(:, :), wnew(:, :)
    real(dp) :: sstep, rf, rstep, smax, ss, delta
    integer :: m, n, work, kreq, k, j, jj, len_res, ierr, iter
    integer :: restart_k, ncopy_u, ncopy_v
    real(dp) :: eps, eps2, svtol
    logical :: converged

    ctl = irlba_control()
    if (present(control)) ctl = control
    m = op%nrow
    n = op%ncol
    kreq = max(nu, nv)
    if (kreq <= 0) error stop "irlba: requested rank must be positive"
    if (kreq >= min(m, n)) error stop "irlba: requested rank must be less than min(m,n)"
    work = ctl%work
    if (work <= 0) work = kreq + 7
    work = max(work, kreq + 1, 4)
    work = min(work, min(m, n))
    if (work <= kreq) then
      call full_svd_operator(op, nu, nv, result, scale, center, shift, .false.)
      return
    end if
    if (present(scale)) then
      if (size(scale) /= n) error stop "irlba: scale length mismatch"
      if (any(abs(scale) <= tiny(1.0_dp))) error stop "irlba: scale entries must be nonzero"
    end if
    if (present(center)) then
      if (size(center) /= n) error stop "irlba: center length mismatch"
    end if

    eps = epsilon(1.0_dp)
    eps2 = ctl%invariant_tol
    if (eps2 <= 0.0_dp) eps2 = eps ** (4.0_dp / 5.0_dp)
    svtol = ctl%svtol
    if (svtol < 0.0_dp) svtol = ctl%tol

    allocate(v(n, work), w(m, work), b(work, work), f(n), y(m), tmp(max(m, n)))
    allocate(res(work), svratio(work))
    v = 0.0_dp
    w = 0.0_dp
    b = 0.0_dp
    svratio = 0.0_dp
    restart_k = 0

    if (present(restart_from)) then
      if (allocated(restart_from%d) .and. allocated(restart_from%u) .and. allocated(restart_from%v)) then
        restart_k = min(size(restart_from%d), work - 1)
        ncopy_u = min(restart_k, size(restart_from%u, 2))
        ncopy_v = min(restart_k, size(restart_from%v, 2))
        restart_k = min(restart_k, ncopy_u, ncopy_v)
        if (restart_k > 0) then
          v(:, 1:restart_k) = restart_from%v(:, 1:restart_k)
          w(:, 1:restart_k) = restart_from%u(:, 1:restart_k)
          do j = 1, restart_k
            b(j, j) = restart_from%d(j)
          end do
          call fill_normal(f)
          call orthogonalize(v, f, restart_k)
          ss = vec_norm2(f)
          if (ss <= eps2) error stop "irlba: failed to construct restart vector"
          v(:, restart_k + 1) = f / ss
        end if
      end if
    end if

    if (restart_k == 0) then
      if (present(start_v)) then
        if (size(start_v) /= n) error stop "irlba: starting vector length mismatch"
        v(:, 1) = start_v
      else
        call fill_normal(v(:, 1))
      end if
    end if

    iter = 0
    result%mprod = 0
    k = restart_k
    converged = .false.
    smax = 0.0_dp

    do while (iter < ctl%maxit)
      if (restart_k == 0 .and. iter == 0) then
        ss = vec_norm2(v(:, 1))
        if (ss < eps2) then
          result%info = -4
          return
        end if
        v(:, 1) = v(:, 1) / ss
        j = 1
      else
        j = k + 1
      end if

      call apply_forward(op, v(:, j), y, scale, center, shift)
      result%mprod = result%mprod + 1
      w(:, j) = y
      if (iter > 0 .and. j > 1 .and. ctl%reorth) call orthogonalize(w, w(:, j), j - 1)
      sstep = vec_norm2(w(:, j))
      if (sstep < eps2 .and. j == 1) then
        result%info = -4
        return
      end if
      if (sstep < eps2) then
        call fill_normal(w(:, j))
        if (j > 1) call orthogonalize(w, w(:, j), j - 1)
        ss = vec_norm2(w(:, j))
        if (ss <= eps2) then
          result%info = -5
          return
        end if
        w(:, j) = w(:, j) / ss
        sstep = 0.0_dp
      else
        w(:, j) = w(:, j) / sstep
      end if

      do while (j <= work)
        call apply_transpose(op, w(:, j), f, scale, center, shift)
        result%mprod = result%mprod + 1
        f = f - sstep * v(:, j)
        call orthogonalize(v, f, j)

        if (j < work) then
          rf = vec_norm2(f)
          rstep = rf
          if (rf < eps2) then
            call fill_normal(f)
            call orthogonalize(v, f, j)
            ss = vec_norm2(f)
            if (ss <= eps2) then
              result%info = -5
              return
            end if
            v(:, j + 1) = f / ss
            rstep = 0.0_dp
          else
            v(:, j + 1) = f / rf
          end if
          b(j, j) = sstep
          b(j, j + 1) = rstep

          call apply_forward(op, v(:, j + 1), y, scale, center, shift)
          result%mprod = result%mprod + 1
          w(:, j + 1) = y - rstep * w(:, j)
          if (ctl%reorth) call orthogonalize(w, w(:, j + 1), j)
          sstep = vec_norm2(w(:, j + 1))
          if (sstep < eps2) then
            call fill_normal(w(:, j + 1))
            call orthogonalize(w, w(:, j + 1), j)
            ss = vec_norm2(w(:, j + 1))
            if (ss <= eps2) then
              result%info = -5
              return
            end if
            w(:, j + 1) = w(:, j + 1) / ss
            sstep = 0.0_dp
          else
            w(:, j + 1) = w(:, j + 1) / sstep
          end if
        else
          b(j, j) = sstep
        end if
        j = j + 1
      end do

      call svd_real(b, bs, ub, vb, ierr)
      if (ierr /= 0) then
        result%info = ierr
        return
      end if
      rf = vec_norm2(f)
      if (rf > eps2) f = f / rf
      if (rf < eps2) rf = 0.0_dp
      smax = max(smax, maxval(bs))
      smax = max(eps ** (2.0_dp / 3.0_dp), smax)

      do jj = 1, work
        if (bs(jj) > 0.0_dp) then
          delta = abs(svratio(jj) - bs(jj)) / bs(jj)
        else
          delta = 0.0_dp
        end if
        svratio(jj) = delta
        res(jj) = rf * ub(work, jj)
      end do

      len_res = 0
      do jj = 1, work
        if (abs(res(jj)) < ctl%tol * smax .and. svratio(jj) < svtol) len_res = len_res + 1
      end do
      if (len_res >= kreq .or. sstep <= 0.0_dp) then
        converged = .true.
        iter = iter + 1
        exit
      end if

      if (k < kreq + len_res) k = kreq + len_res
      k = min(k, work - 3)
      k = max(k, 1)
      if (k < kreq) k = min(kreq, work - 1)
      svratio = bs

      allocate(vnew(n, k), wnew(m, k))
      vnew = matmul(v(:, 1:work), vb(:, 1:k))
      wnew = matmul(w(:, 1:work), ub(:, 1:k))
      v = 0.0_dp
      w = 0.0_dp
      v(:, 1:k) = vnew
      w(:, 1:k) = wnew
      if (k + 1 <= work) v(:, k + 1) = f
      b = 0.0_dp
      do jj = 1, k
        b(jj, jj) = bs(jj)
        if (k + 1 <= work) b(jj, k + 1) = res(jj)
      end do
      deallocate(vnew, wnew, bs, ub, vb)
      iter = iter + 1
      restart_k = 1
    end do

    result%iter = iter
    result%converged = converged
    if (converged) then
      result%info = 0
    else
      result%info = -2
    end if

    if (.not. allocated(bs)) then
      call svd_real(b, bs, ub, vb, ierr)
      if (ierr /= 0) then
        result%info = ierr
        return
      end if
    end if
    allocate(result%d(kreq), result%u(m, nu), result%v(n, nv))
    result%d = bs(1:kreq)
    result%u = matmul(w(:, 1:work), ub(:, 1:nu))
    result%v = matmul(v(:, 1:work), vb(:, 1:nv))
  end subroutine irlb_operator

end module irlba_core
