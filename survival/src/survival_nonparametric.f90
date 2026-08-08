! SPDX-License-Identifier: LGPL-2.0-or-later
module survival_nonparametric
  use survival_kinds, only : dp
  use survival_types, only : survfit_result
  implicit none
  private
  public :: kaplan_meier, kaplan_meier_counting, aalen_johansen

contains

  subroutine kaplan_meier(time, status, fit, weights)
    real(dp), intent(in) :: time(:)
    integer, intent(in) :: status(:)
    type(survfit_result), intent(out) :: fit
    real(dp), intent(in), optional :: weights(:)

    real(dp), allocatable :: w(:), ut(:)
    real(dp) :: nrisk, deaths, cens, km, ch, vkm, vch
    integer :: n, m, i, k

    n = size(time)
    allocate(w(n), ut(n))
    if (present(weights)) then
      w = weights
    else
      w = 1.0_dp
    end if

    m = 0
    do i = 1, n
      if (.not. contains_time(ut, m, time(i))) then
        m = m + 1
        ut(m) = time(i)
      end if
    end do
    call sort_real(ut(1:m))
    call allocate_fit(fit, m)
    fit%time = ut(1:m)

    km = 1.0_dp
    ch = 0.0_dp
    vkm = 0.0_dp
    vch = 0.0_dp
    do k = 1, m
      nrisk = 0.0_dp
      deaths = 0.0_dp
      cens = 0.0_dp
      do i = 1, n
        if (time(i) >= fit%time(k)) nrisk = nrisk + w(i)
        if (same_time(time(i), fit%time(k))) then
          if (status(i) /= 0) deaths = deaths + w(i)
          if (status(i) == 0) cens = cens + w(i)
        end if
      end do
      call update_fit_row(fit, k, nrisk, deaths, cens, km, ch, vkm, vch)
    end do
  end subroutine kaplan_meier

  subroutine kaplan_meier_counting(start, stop, status, fit, weights)
    real(dp), intent(in) :: start(:), stop(:)
    integer, intent(in) :: status(:)
    type(survfit_result), intent(out) :: fit
    real(dp), intent(in), optional :: weights(:)

    real(dp), allocatable :: w(:), ut(:)
    real(dp) :: nrisk, deaths, cens, km, ch, vkm, vch
    integer :: n, m, i, k

    n = size(stop)
    allocate(w(n), ut(n))
    if (present(weights)) then
      w = weights
    else
      w = 1.0_dp
    end if

    m = 0
    do i = 1, n
      if (.not. contains_time(ut, m, stop(i))) then
        m = m + 1
        ut(m) = stop(i)
      end if
    end do
    call sort_real(ut(1:m))
    call allocate_fit(fit, m)
    fit%time = ut(1:m)

    km = 1.0_dp
    ch = 0.0_dp
    vkm = 0.0_dp
    vch = 0.0_dp
    do k = 1, m
      nrisk = 0.0_dp
      deaths = 0.0_dp
      cens = 0.0_dp
      do i = 1, n
        if (start(i) < fit%time(k) .and. stop(i) >= fit%time(k)) then
          nrisk = nrisk + w(i)
        end if
        if (same_time(stop(i), fit%time(k))) then
          if (status(i) /= 0) deaths = deaths + w(i)
          if (status(i) == 0) cens = cens + w(i)
        end if
      end do
      call update_fit_row(fit, k, nrisk, deaths, cens, km, ch, vkm, vch)
    end do
  end subroutine kaplan_meier_counting

  subroutine aalen_johansen(time, from_state, to_state, nstate, prob, utime, initial)
    real(dp), intent(in) :: time(:)
    integer, intent(in) :: from_state(:), to_state(:), nstate
    real(dp), allocatable, intent(out) :: prob(:,:), utime(:)
    real(dp), intent(in), optional :: initial(:)

    real(dp), allocatable :: u(:), risk(:), d(:,:), p(:), tmp(:)
    integer :: n, m, i, k, a, b

    n = size(time)
    allocate(u(n))
    m = 0
    do i = 1, n
      if (.not. contains_time(u, m, time(i))) then
        m = m + 1
        u(m) = time(i)
      end if
    end do
    call sort_real(u(1:m))
    allocate(utime(m), prob(nstate,m), risk(nstate), d(nstate,nstate))
    allocate(p(nstate), tmp(nstate))
    utime = u(1:m)
    p = 0.0_dp
    if (present(initial)) then
      p = initial
    else
      p(1) = 1.0_dp
    end if

    do k = 1, m
      risk = 0.0_dp
      d = 0.0_dp
      do i = 1, n
        a = from_state(i)
        b = to_state(i)
        if (time(i) >= utime(k) .and. a >= 1 .and. a <= nstate) then
          risk(a) = risk(a) + 1.0_dp
        end if
        if (same_time(time(i), utime(k)) .and. b > 0 .and. a /= b) then
          d(a,b) = d(a,b) + 1.0_dp
        end if
      end do
      tmp = p
      do a = 1, nstate
        if (risk(a) <= 0.0_dp) cycle
        do b = 1, nstate
          if (b /= a) tmp(b) = tmp(b) + p(a) * d(a,b) / risk(a)
        end do
        tmp(a) = tmp(a) - p(a) * sum(d(a,:)) / risk(a)
      end do
      p = tmp
      prob(:,k) = p
    end do
  end subroutine aalen_johansen

  subroutine allocate_fit(fit, m)
    type(survfit_result), intent(out) :: fit
    integer, intent(in) :: m
    allocate(fit%time(m), fit%n_risk(m), fit%n_event(m), fit%n_censor(m))
    allocate(fit%survival(m), fit%cumhaz(m), fit%std_err(m), fit%std_chaz(m))
  end subroutine allocate_fit

  subroutine update_fit_row(fit, k, nrisk, deaths, cens, km, ch, vkm, vch)
    type(survfit_result), intent(inout) :: fit
    integer, intent(in) :: k
    real(dp), intent(in) :: nrisk, deaths, cens
    real(dp), intent(inout) :: km, ch, vkm, vch

    fit%n_risk(k) = nrisk
    fit%n_event(k) = deaths
    fit%n_censor(k) = cens
    if (nrisk > 0.0_dp) then
      ch = ch + deaths / nrisk
      vch = vch + deaths / (nrisk * nrisk)
      if (deaths > 0.0_dp) then
        km = km * max(0.0_dp, 1.0_dp - deaths / nrisk)
        if (nrisk > deaths) vkm = vkm + deaths / (nrisk * (nrisk - deaths))
      end if
    end if
    fit%survival(k) = km
    fit%cumhaz(k) = ch
    fit%std_err(k) = km * sqrt(max(0.0_dp, vkm))
    fit%std_chaz(k) = sqrt(max(0.0_dp, vch))
  end subroutine update_fit_row

  logical function contains_time(x, nused, value) result(found)
    real(dp), intent(in) :: x(:), value
    integer, intent(in) :: nused
    integer :: i
    found = .false.
    do i = 1, nused
      if (same_time(x(i), value)) then
        found = .true.
        return
      end if
    end do
  end function contains_time

  pure logical function same_time(a, b) result(equal)
    real(dp), intent(in) :: a, b
    real(dp) :: scale
    scale = max(1.0_dp, abs(a), abs(b))
    equal = abs(a - b) <= 8.0_dp * epsilon(1.0_dp) * scale
  end function same_time

  subroutine sort_real(x)
    real(dp), intent(inout) :: x(:)
    integer :: i, j
    real(dp) :: value
    do i = 2, size(x)
      value = x(i)
      j = i - 1
      do while (j >= 1)
        if (x(j) <= value) exit
        x(j+1) = x(j)
        j = j - 1
      end do
      x(j+1) = value
    end do
  end subroutine sort_real

end module survival_nonparametric
