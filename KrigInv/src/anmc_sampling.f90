! SPDX-License-Identifier: GPL-3.0-only
module anmc_sampling
  use anmc_kinds, only : dp
  use anmc_types, only : simulation_control
  use anmc_math, only : random_normals, cholesky_lower, pmvnorm, probability_control, probability_result, genz_bretz
  implicit none
  private
  public :: mvrnorm_arma, trmvrnorm_rej_cpp

contains

  function mvrnorm_arma(n, mu, sigma, chol, ok, message) result(x)
    integer, intent(in) :: n
    real(dp), intent(in) :: mu(:), sigma(:,:)
    integer, intent(in), optional :: chol
    logical, intent(out), optional :: ok
    character(len=*), intent(out), optional :: message
    real(dp), allocatable :: x(:,:)
    real(dp), allocatable :: l(:,:), z(:)
    logical :: lok
    character(len=256) :: msg
    integer :: i, p, ichol

    p = size(mu)
    ichol = 0
    if (present(chol)) ichol = chol
    allocate(x(p, max(0,n)))
    if (n <= 0) then
      if (present(ok)) ok = .true.
      if (present(message)) message = ''
      return
    end if
    if (size(sigma,1) /= p .or. size(sigma,2) /= p) then
      x = 0.0_dp
      if (present(ok)) ok = .false.
      if (present(message)) message = 'non-conforming covariance/Cholesky dimensions'
      return
    end if

    if (ichol /= 0) then
      ! Match upstream Armadillo implementation: sigma is the upper
      ! Cholesky factor U and each draw is mu + transpose(U) z.
      l = transpose(sigma)
      lok = .true.
      msg = ''
    else
      call cholesky_lower(sigma, l, lok, msg)
    end if
    if (.not. lok) then
      x = 0.0_dp
      if (present(ok)) ok = .false.
      if (present(message)) message = trim(msg)
      return
    end if

    allocate(z(p))
    do i = 1, n
      call random_normals(z)
      x(:,i) = mu + matmul(l, z)
    end do
    if (present(ok)) ok = .true.
    if (present(message)) message = ''
  end function mvrnorm_arma

  function trmvrnorm_rej_cpp(n, mu, sigma, lower, upper, verb, sim_control, &
                             prob_control, ok, total_draws) result(y)
    integer, intent(in) :: n
    real(dp), intent(in) :: mu(:), sigma(:,:), lower(:), upper(:)
    integer, intent(in), optional :: verb
    type(simulation_control), intent(in), optional :: sim_control
    type(probability_control), intent(in), optional :: prob_control
    logical, intent(out), optional :: ok
    integer, intent(out), optional :: total_draws
    real(dp), allocatable :: y(:,:)

    type(simulation_control) :: sctl
    type(probability_control) :: pctl
    type(probability_result) :: pres
    real(dp), allocatable :: x(:,:)
    real(dp) :: alpha
    integer :: p, remaining, accepted, batch, i, j, v, draws, good
    logical :: lok
    character(len=256) :: msg

    sctl = simulation_control()
    if (present(sim_control)) sctl = sim_control
    pctl = genz_bretz()
    if (present(prob_control)) pctl = prob_control
    v = 0
    if (present(verb)) v = verb
    p = size(mu)
    allocate(y(p, max(0,n)))
    y = 0.0_dp

    if (n <= 0) then
      if (present(ok)) ok = .true.
      if (present(total_draws)) total_draws = 0
      return
    end if
    if (size(sigma,1) /= p .or. size(sigma,2) /= p .or. &
        size(lower) /= p .or. size(upper) /= p) then
      if (present(ok)) ok = .false.
      if (present(total_draws)) total_draws = 0
      return
    end if

    pres = pmvnorm(lower, upper, mu, sigma, pctl)
    alpha = max(0.0_dp, min(1.0_dp, pres%value))
    if (v >= 3) write(*,'(a,es12.4)') 'Acceptance rate: ', alpha
    if (alpha <= 0.0_dp) then
      if (present(ok)) ok = .false.
      if (present(total_draws)) total_draws = 0
      return
    end if

    remaining = n
    accepted = 0
    draws = 0
    do while (remaining > 0)
      if (alpha > tiny(1.0_dp)) then
        if (real(remaining,dp)/alpha >= real(sctl%max_rejection_batch,dp)) then
          batch = sctl%max_rejection_batch
        else
          batch = max(10, ceiling(real(remaining,dp) / alpha))
        end if
      else
        batch = sctl%max_rejection_batch
      end if
      batch = min(batch, sctl%max_rejection_batch)
      if (draws + batch > sctl%max_rejection_draws) then
        batch = sctl%max_rejection_draws - draws
      end if
      if (batch <= 0) exit

      x = mvrnorm_arma(batch, mu, sigma, 0, lok, msg)
      if (.not. lok) exit
      draws = draws + batch
      good = 0
      do i = 1, batch
        if (all(x(:,i) >= lower) .and. all(x(:,i) <= upper)) good = good + 1
      end do
      if (good == 0) cycle
      alpha = real(good,dp) / real(batch,dp)
      if (v >= 4) write(*,'(a,es12.4)') 'Current acceptance rate: ', alpha
      j = 0
      do i = 1, batch
        if (all(x(:,i) >= lower) .and. all(x(:,i) <= upper)) then
          j = j + 1
          if (accepted < n) then
            accepted = accepted + 1
            y(:,accepted) = x(:,i)
            remaining = remaining - 1
          end if
          if (remaining == 0) exit
        end if
      end do
    end do

    lok = (accepted == n)
    if (v >= 3) then
      write(*,'(a,i0)') 'Total samples run: ', draws
      write(*,'(a,i0)') 'Total samples accepted: ', accepted
      if (draws > 0) write(*,'(a,es12.4)') 'Accepted/requested draw ratio: ', real(accepted,dp)/real(draws,dp)
    end if
    if (present(ok)) ok = lok
    if (present(total_draws)) total_draws = draws
  end function trmvrnorm_rej_cpp

end module anmc_sampling
