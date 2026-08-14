module dice_design_lhs_opt
  use iso_fortran_env, only : int64
  use dice_design_kinds, only : dp
  use dice_design_rng, only : rng_state
  use dice_design_criteria, only : discrepancy_value, phi_p
  implicit none
  private

  type, public :: lhs_optimization_result
    real(dp), allocatable :: initial_design(:, :)
    real(dp), allocatable :: design(:, :)
    real(dp), allocatable :: crit_values(:)
    real(dp), allocatable :: temp_values(:)
    real(dp), allocatable :: proba_values(:)
    integer :: steps = 0
  end type lhs_optimization_result

  public :: discrep_sa_lhs, discrep_ese_lhs, maximin_sa_lhs, maximin_ese_lhs

contains

  subroutine swap_random(m, rng, candidate, column, row1, row2)
    real(dp), intent(in) :: m(:, :)
    type(rng_state), intent(inout) :: rng
    real(dp), allocatable, intent(out) :: candidate(:, :)
    integer, intent(in), optional :: column
    integer, intent(out), optional :: row1, row2
    integer :: i1, i2, k
    real(dp) :: tmp

    allocate(candidate(size(m,1),size(m,2)))
    candidate = m
    i1 = rng%integer(1,size(m,1))
    i2 = rng%integer(1,size(m,1))
    if (present(column)) then
      k = column
    else
      k = rng%integer(1,size(m,2))
    end if
    tmp = candidate(i1,k)
    candidate(i1,k) = candidate(i2,k)
    candidate(i2,k) = tmp
    if (present(row1)) row1 = i1
    if (present(row2)) row2 = i2
  end subroutine swap_random

  subroutine replace_random_row(m, rng, candidate)
    real(dp), intent(in) :: m(:, :)
    type(rng_state), intent(inout) :: rng
    real(dp), allocatable, intent(out) :: candidate(:, :)
    integer :: row, j

    allocate(candidate(size(m,1),size(m,2)))
    candidate = m
    row = rng%integer(1,size(m,1))
    do j = 1, size(m,2)
      candidate(row,j) = rng%uniform()
    end do
  end subroutine replace_random_row

  subroutine append_value(x, n, value)
    real(dp), allocatable, intent(inout) :: x(:)
    integer, intent(inout) :: n
    real(dp), intent(in) :: value
    real(dp), allocatable :: tmp(:)
    integer :: newcap

    if (.not. allocated(x)) allocate(x(64))
    if (n >= size(x)) then
      newcap = max(2*size(x),1)
      allocate(tmp(newcap))
      if (n>0) tmp(1:n) = x(1:n)
      call move_alloc(tmp,x)
    end if
    n = n+1
    x(n) = value
  end subroutine append_value

  subroutine finish_history(x, n)
    real(dp), allocatable, intent(inout) :: x(:)
    integer, intent(in) :: n
    real(dp), allocatable :: tmp(:)
    if (.not. allocated(x)) then
      allocate(x(0))
    else if (size(x) /= n) then
      allocate(tmp(n))
      if (n>0) tmp = x(1:n)
      call move_alloc(tmp,x)
    end if
  end subroutine finish_history

  pure function acceptance_probability(oldv, newv, temperature) result(prob)
    real(dp), intent(in) :: oldv, newv, temperature
    real(dp) :: prob
    if (temperature <= 0.0_dp) then
      if (newv <= oldv) then
        prob = 1.0_dp
      else
        prob = 0.0_dp
      end if
    else
      prob = min(1.0_dp,exp((oldv-newv)/temperature))
    end if
  end function acceptance_probability

  subroutine discrep_sa_lhs(initial_design, result, t0, cooling, iterations, criterion, profile, imax, seed)
    real(dp), intent(in) :: initial_design(:, :)
    type(lhs_optimization_result), intent(out) :: result
    real(dp), intent(in), optional :: t0, cooling
    integer, intent(in), optional :: iterations, imax
    character(len=*), intent(in), optional :: criterion, profile
    integer(int64), intent(in), optional :: seed
    type(rng_state) :: rng
    real(dp), allocatable :: m(:, :), cand(:, :), best_design(:, :)
    real(dp) :: temp0, c, t, v, g, prob, ref
    integer :: itmax, imaxv, i, idle, ncrit, ntemp, nprob
    logical :: accepted_any
    character(len=16) :: crit, prof

    temp0 = 10.0_dp
    if (present(t0)) temp0 = t0
    c = 0.95_dp
    if (present(cooling)) c = cooling
    itmax = 2000
    if (present(iterations)) itmax = iterations
    imaxv = 100
    if (present(imax)) imaxv = imax
    crit = 'C2'
    if (present(criterion)) crit = criterion
    prof = 'GEOM'
    if (present(profile)) prof = profile
    call rng%seed(1_int64)
    if (present(seed)) call rng%seed(seed)
    allocate(m(size(initial_design,1),size(initial_design,2)),result%initial_design(size(initial_design,1),size(initial_design,2)))
    m = initial_design
    result%initial_design = initial_design
    v = discrepancy_value(m,trim(crit))
    ncrit = 0
    ntemp = 0
    nprob = 0
    call append_value(result%crit_values,ncrit,v)

    select case (trim(prof))
    case ('GEOM','geom','LINEAR','linear','MC','mc')
      if ((trim(prof)=='MC' .or. trim(prof)=='mc') .and. .not. (trim(crit)=='C2' .or. trim(crit)=='c2')) then
        error stop 'discrep_sa_lhs: upstream MC profile is defined only for C2'
      end if
      t = temp0
      do i = 1, itmax
        if (t <= 0.0_dp) exit
        if (trim(prof)=='MC' .or. trim(prof)=='mc') then
          call replace_random_row(m,rng,cand)
        else
          call swap_random(m,rng,cand)
        end if
        g = discrepancy_value(cand,trim(crit))
        prob = acceptance_probability(v,g,t)
        if (rng%uniform() < prob) then
          m = cand
          v = g
        end if
        call append_value(result%crit_values,ncrit,v)
        call append_value(result%temp_values,ntemp,t)
        call append_value(result%proba_values,nprob,prob)
        if (trim(prof)=='LINEAR' .or. trim(prof)=='linear') then
          t = temp0*(1.0_dp-real(i,dp)/real(itmax,dp))
        else
          t = temp0*c**i
        end if
        deallocate(cand)
      end do
      allocate(result%design(size(m,1),size(m,2)))
      result%design = m
    case ('GEOM_MORRIS','geom_morris')
      t = temp0
      allocate(best_design(size(m,1),size(m,2)))
      best_design = m
      ref = v
      do i = 1, itmax
        accepted_any = .false.
        idle = 0
        do while (idle < imaxv)
          call swap_random(m,rng,cand)
          g = discrepancy_value(cand,trim(crit))
          prob = acceptance_probability(v,g,t)
          if (rng%uniform() < prob) then
            m = cand
            v = g
            accepted_any = .true.
            if (v < ref) then
              ref = v
              best_design = m
              idle = 0
            else
              idle = idle+1
            end if
          else
            idle = idle+1
          end if
          call append_value(result%crit_values,ncrit,ref)
          call append_value(result%temp_values,ntemp,t)
          call append_value(result%proba_values,nprob,prob)
          deallocate(cand)
        end do
        if (accepted_any) then
          t = t*c
        else
          exit
        end if
      end do
      allocate(result%design(size(m,1),size(m,2)))
      result%design = best_design
    case default
      error stop 'discrep_sa_lhs: unknown temperature profile'
    end select
    result%steps = nprob
    call finish_history(result%crit_values,ncrit)
    call finish_history(result%temp_values,ntemp)
    call finish_history(result%proba_values,nprob)
  end subroutine discrep_sa_lhs

  subroutine maximin_sa_lhs(initial_design, result, t0, cooling, iterations, p, profile, imax, seed)
    real(dp), intent(in) :: initial_design(:, :)
    type(lhs_optimization_result), intent(out) :: result
    real(dp), intent(in), optional :: t0, cooling, p
    integer, intent(in), optional :: iterations, imax
    character(len=*), intent(in), optional :: profile
    integer(int64), intent(in), optional :: seed
    type(rng_state) :: rng
    real(dp), allocatable :: m(:, :), cand(:, :), best_design(:, :)
    real(dp) :: temp0, c, t, v, g, prob, ref, pp
    integer :: itmax, imaxv, i, idle, ncrit, ntemp, nprob
    logical :: accepted_any
    character(len=16) :: prof

    temp0 = 10.0_dp
    if (present(t0)) temp0 = t0
    c = 0.95_dp
    if (present(cooling)) c = cooling
    itmax = 2000
    if (present(iterations)) itmax = iterations
    pp = 50.0_dp
    if (present(p)) pp = p
    imaxv = 100
    if (present(imax)) imaxv = imax
    prof = 'GEOM'
    if (present(profile)) prof = profile
    call rng%seed(1_int64)
    if (present(seed)) call rng%seed(seed)
    allocate(m(size(initial_design,1),size(initial_design,2)),result%initial_design(size(initial_design,1),size(initial_design,2)))
    m = initial_design
    result%initial_design = initial_design
    v = phi_p(m,pp)
    ncrit = 0
    ntemp = 0
    nprob = 0
    call append_value(result%crit_values,ncrit,v)

    select case (trim(prof))
    case ('GEOM','geom','LINEAR','linear','MC','mc')
      t = temp0
      do i = 1, itmax
        if (t <= 0.0_dp) exit
        if (trim(prof)=='MC' .or. trim(prof)=='mc') then
          call replace_random_row(m,rng,cand)
        else
          call swap_random(m,rng,cand)
        end if
        g = phi_p(cand,pp)
        prob = acceptance_probability(v,g,t)
        if (rng%uniform() < prob) then
          m = cand
          v = g
        end if
        call append_value(result%crit_values,ncrit,v)
        call append_value(result%temp_values,ntemp,t)
        call append_value(result%proba_values,nprob,prob)
        if (trim(prof)=='LINEAR' .or. trim(prof)=='linear') then
          t = temp0*(1.0_dp-real(i,dp)/real(itmax,dp))
        else
          t = temp0*c**i
        end if
        deallocate(cand)
      end do
      allocate(result%design(size(m,1),size(m,2)))
      result%design = m
    case ('GEOM_MORRIS','geom_morris')
      t = temp0
      allocate(best_design(size(m,1),size(m,2)))
      best_design = m
      ref = v
      do i = 1, itmax
        accepted_any = .false.
        idle = 0
        do while (idle < imaxv)
          call swap_random(m,rng,cand)
          g = phi_p(cand,pp)
          prob = acceptance_probability(v,g,t)
          if (rng%uniform() < prob) then
            m = cand
            v = g
            accepted_any = .true.
            if (v < ref) then
              ref = v
              best_design = m
              idle = 0
            else
              idle = idle+1
            end if
          else
            idle = idle+1
          end if
          call append_value(result%crit_values,ncrit,ref)
          call append_value(result%temp_values,ntemp,t)
          call append_value(result%proba_values,nprob,prob)
          deallocate(cand)
        end do
        if (accepted_any) then
          t = t*c
        else
          exit
        end if
      end do
      allocate(result%design(size(m,1),size(m,2)))
      result%design = best_design
    case default
      error stop 'maximin_sa_lhs: unknown temperature profile'
    end select
    result%steps = nprob
    call finish_history(result%crit_values,ncrit)
    call finish_history(result%temp_values,ntemp)
    call finish_history(result%proba_values,nprob)
  end subroutine maximin_sa_lhs

  subroutine discrep_ese_lhs(initial_design, result, t0, inner_iterations, candidates, outer_iterations, criterion, seed)
    real(dp), intent(in) :: initial_design(:, :)
    type(lhs_optimization_result), intent(out) :: result
    real(dp), intent(in), optional :: t0
    integer, intent(in), optional :: inner_iterations, candidates, outer_iterations
    character(len=*), intent(in), optional :: criterion
    integer(int64), intent(in), optional :: seed
    type(rng_state) :: rng
    real(dp), allocatable :: m(:, :), cand(:, :), bestcand(:, :), best_design(:, :)
    real(dp) :: temperature, best, current, a, value, delta, v1, v2, bold, temp0
    integer :: inner, nj, outer, q, count, j, col, na, ni, ncrit, ntemp, nprob, start_hist
    character(len=16) :: crit

    crit = 'C2'
    if (present(criterion)) crit = criterion
    inner = 100
    if (present(inner_iterations)) inner = inner_iterations
    nj = 50
    if (present(candidates)) nj = candidates
    if (inner < 1 .or. nj < 1) error stop 'discrep_ese_lhs: iteration counts must be positive'
    outer = 2
    if (present(outer_iterations)) outer = outer_iterations
    temp0 = 0.005_dp*discrepancy_value(initial_design,'C2')
    if (present(t0)) temp0 = t0
    call rng%seed(1_int64)
    if (present(seed)) call rng%seed(seed)
    allocate(m(size(initial_design,1),size(initial_design,2)),best_design(size(initial_design,1),size(initial_design,2)), &
      bestcand(size(initial_design,1),size(initial_design,2)), &
      result%initial_design(size(initial_design,1),size(initial_design,2)))
    m = initial_design
    best_design = m
    bestcand = m
    result%initial_design = initial_design
    current = discrepancy_value(m,trim(crit))
    best = current
    temperature = temp0
    ncrit = 0
    ntemp = 0
    nprob = 0
    call append_value(result%crit_values,ncrit,current)

    do q = 1, outer
      bold = best
      na = 0
      ni = 0
      start_hist = ntemp+1
      do count = 1, inner
        col = modulo(count-1,size(m,2))+1
        a = huge(1.0_dp)
        do j = 1, nj
          call swap_random(m,rng,cand,column=col)
          value = discrepancy_value(cand,trim(crit))
          if (value < a) then
            a = value
            bestcand = cand
          end if
          deallocate(cand)
        end do
        delta = a-current
        if (delta <= temperature*rng%uniform()) then
          m = bestcand
          current = a
          na = na+1
          if (a <= best) then
            best = a
            best_design = m
            ni = ni+1
          end if
        end if
        call append_value(result%crit_values,ncrit,best)
        call append_value(result%temp_values,ntemp,temperature)
        call append_value(result%proba_values,nprob,0.0_dp)
      end do
      v1 = real(na,dp)/real(max(inner,1),dp)
      v2 = real(ni,dp)/real(max(inner,1),dp)
      result%proba_values(start_hist:nprob) = v1
      if (best-bold < 0.0_dp) then
        if (v1 >= 0.1_dp .and. v2 <= v1) then
          temperature = 0.8_dp*temperature
        else
          temperature = temperature/0.8_dp
        end if
      else
        if (v1 <= 0.1_dp) then
          temperature = temperature/0.7_dp
        else
          temperature = 0.9_dp*temperature
        end if
      end if
    end do
    allocate(result%design(size(m,1),size(m,2)))
    result%design = best_design
    result%steps = nprob
    call finish_history(result%crit_values,ncrit)
    call finish_history(result%temp_values,ntemp)
    call finish_history(result%proba_values,nprob)
  end subroutine discrep_ese_lhs

  subroutine maximin_ese_lhs(initial_design, result, t0, inner_iterations, candidates, outer_iterations, p, seed)
    real(dp), intent(in) :: initial_design(:, :)
    type(lhs_optimization_result), intent(out) :: result
    real(dp), intent(in), optional :: t0, p
    integer, intent(in), optional :: inner_iterations, candidates, outer_iterations
    integer(int64), intent(in), optional :: seed
    type(rng_state) :: rng
    real(dp), allocatable :: m(:, :), cand(:, :), bestcand(:, :), best_design(:, :)
    real(dp) :: temperature, best, current, a, value, delta, v1, v2, bold, temp0, pp
    integer :: inner, nj, outer, q, count, j, col, na, ni, ncrit, ntemp, nprob, start_hist

    inner = 100
    if (present(inner_iterations)) inner = inner_iterations
    nj = 50
    if (present(candidates)) nj = candidates
    if (inner < 1 .or. nj < 1) error stop 'maximin_ese_lhs: iteration counts must be positive'
    outer = 1
    if (present(outer_iterations)) outer = outer_iterations
    pp = 50.0_dp
    if (present(p)) pp = p
    temp0 = 0.005_dp*phi_p(initial_design,50.0_dp)
    if (present(t0)) temp0 = t0
    call rng%seed(1_int64)
    if (present(seed)) call rng%seed(seed)
    allocate(m(size(initial_design,1),size(initial_design,2)),best_design(size(initial_design,1),size(initial_design,2)), &
      bestcand(size(initial_design,1),size(initial_design,2)), &
      result%initial_design(size(initial_design,1),size(initial_design,2)))
    m = initial_design
    best_design = m
    bestcand = m
    result%initial_design = initial_design
    current = phi_p(m,pp)
    best = current
    temperature = temp0
    ncrit = 0
    ntemp = 0
    nprob = 0
    call append_value(result%crit_values,ncrit,current)

    do q = 1, outer
      bold = best
      na = 0
      ni = 0
      start_hist = ntemp+1
      do count = 1, inner
        col = modulo(count-1,size(m,2))+1
        a = huge(1.0_dp)
        do j = 1, nj
          call swap_random(m,rng,cand,column=col)
          value = phi_p(cand,pp)
          if (value < a) then
            a = value
            bestcand = cand
          end if
          deallocate(cand)
        end do
        delta = a-current
        if (delta <= temperature*rng%uniform()) then
          m = bestcand
          current = a
          na = na+1
          if (a <= best) then
            best = a
            best_design = m
            ni = ni+1
          end if
        end if
        call append_value(result%crit_values,ncrit,best)
        call append_value(result%temp_values,ntemp,temperature)
        call append_value(result%proba_values,nprob,0.0_dp)
      end do
      v1 = real(na,dp)/real(max(inner,1),dp)
      v2 = real(ni,dp)/real(max(inner,1),dp)
      result%proba_values(start_hist:nprob) = v1
      if (best-bold < 0.0_dp) then
        if (v1 >= 0.1_dp .and. v2 <= v1) then
          temperature = 0.8_dp*temperature
        else
          temperature = temperature/0.8_dp
        end if
      else
        if (v1 <= 0.1_dp) then
          temperature = temperature/0.7_dp
        else
          temperature = 0.9_dp*temperature
        end if
      end if
    end do
    allocate(result%design(size(m,1),size(m,2)))
    result%design = best_design
    result%steps = nprob
    call finish_history(result%crit_values,ncrit)
    call finish_history(result%temp_values,ntemp)
    call finish_history(result%proba_values,nprob)
  end subroutine maximin_ese_lhs

end module dice_design_lhs_opt
