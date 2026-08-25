! SPDX-License-Identifier: LGPL-2.0-or-later
module survival_stats
  use survival_kinds, only : dp
  use survival_types, only : concordance_result, survdiff_result
  use survival_linalg, only : solve_sym
  implicit none
  private
  public :: survdiff, concordance_right
contains

  subroutine survdiff(time, status, group, ngroup, result, rho, strata)
    real(dp), intent(in) :: time(:)
    integer, intent(in) :: status(:), group(:), ngroup
    type(survdiff_result), intent(out) :: result
    real(dp), intent(in), optional :: rho
    integer, intent(in), optional :: strata(:)

    real(dp) :: r, wt, nrisk, deaths, km, tmp, chisq
    real(dp), allocatable :: ut(:), risk(:), obs(:), ex(:), var(:,:)
    real(dp), allocatable :: kmleft(:), diff(:), vv(:,:), sol(:)
    integer, allocatable :: st(:)
    integer :: n, i, j, g, h, s, ns, m
    logical :: ok

    n = size(time)
    r = 0.0_dp
    if (present(rho)) r = rho
    allocate(st(n))
    if (present(strata)) then
      st = strata
    else
      st = 1
    end if
    ns = maxval(st)
    allocate(obs(ngroup), ex(ngroup), var(ngroup,ngroup))
    obs = 0.0_dp
    ex = 0.0_dp
    var = 0.0_dp

    do s = 1, ns
      allocate(ut(count(st==s)), kmleft(count(st==s)))
      m = 0
      do i = 1, n
        if (st(i) == s) then
          m = m + 1
          ut(m) = time(i)
        end if
      end do
      call sort_real(ut)
      km = 1.0_dp
      do i = 1, size(ut)
        kmleft(i) = km
        nrisk = 0.0_dp
        deaths = 0.0_dp
        do j = 1, n
          if (st(j) /= s) cycle
          if (time(j) >= ut(i)) nrisk = nrisk + 1.0_dp
          if (same_time(time(j),ut(i)) .and. status(j) /= 0) then
            deaths = deaths + 1.0_dp
          end if
        end do
        if (nrisk > 0.0_dp .and. deaths > 0.0_dp) then
          km = km * (nrisk-deaths)/nrisk
        end if
      end do

      do i = 1, size(ut)
        if (i > 1) then
          if (same_time(ut(i),ut(i-1))) cycle
        end if
        allocate(risk(ngroup))
        risk = 0.0_dp
        nrisk = 0.0_dp
        deaths = 0.0_dp
        do j = 1, n
          if (st(j) /= s) cycle
          if (time(j) >= ut(i)) then
            nrisk = nrisk + 1.0_dp
            risk(group(j)) = risk(group(j)) + 1.0_dp
          end if
          if (same_time(time(j),ut(i)) .and. status(j) /= 0) then
            deaths = deaths + 1.0_dp
          end if
        end do
        if (deaths <= 0.0_dp .or. nrisk <= 0.0_dp) then
          deallocate(risk)
          cycle
        end if
        wt = 1.0_dp
        if (abs(r) > tiny(1.0_dp)) wt = kmleft(i)**r
        do g = 1, ngroup
          do j = 1, n
            if (st(j)==s .and. same_time(time(j),ut(i)) .and. &
                status(j)/=0 .and. group(j)==g) then
              obs(g) = obs(g) + wt
            end if
          end do
          ex(g) = ex(g) + wt*deaths*risk(g)/nrisk
        end do
        if (nrisk > 1.0_dp) then
          do g = 1, ngroup
            tmp = wt*wt*deaths*risk(g)*(nrisk-deaths) / &
                  (nrisk*(nrisk-1.0_dp))
            var(g,g) = var(g,g) + tmp*(1.0_dp-risk(g)/nrisk)
            do h = 1, ngroup
              if (h /= g) var(g,h) = var(g,h) - tmp*risk(h)/nrisk
            end do
          end do
        end if
        deallocate(risk)
      end do
      deallocate(ut,kmleft)
    end do

    allocate(result%observed(ngroup), result%expected(ngroup))
    allocate(result%variance(ngroup,ngroup))
    result%observed = obs
    result%expected = ex
    result%variance = var
    if (ngroup > 1) then
      allocate(diff(ngroup-1), vv(ngroup-1,ngroup-1), sol(ngroup-1))
      diff = obs(1:ngroup-1) - ex(1:ngroup-1)
      vv = var(1:ngroup-1,1:ngroup-1)
      call solve_sym(vv,diff,sol,ok)
      if (ok) then
        chisq = dot_product(diff,sol)
        result%chisq = max(0.0_dp,chisq)
      end if
    end if
  end subroutine survdiff

  subroutine concordance_right(time,status,risk,result,weights)
    real(dp), intent(in) :: time(:), risk(:)
    integer, intent(in) :: status(:)
    type(concordance_result), intent(out) :: result
    real(dp), intent(in), optional :: weights(:)

    real(dp), allocatable :: w(:)
    real(dp) :: pair_weight, denom
    integer :: n, i, j

    n = size(time)
    allocate(w(n))
    if (present(weights)) then
      w = weights
    else
      w = 1.0_dp
    end if
    do i = 1, n-1
      do j = i+1, n
        pair_weight = w(i)*w(j)
        if (same_time(time(i),time(j)) .and. status(i)/=0 .and. status(j)/=0) then
          result%tied_time = result%tied_time + pair_weight
        else if (time(i) < time(j) .and. status(i)/=0) then
          call add_pair(risk(i),risk(j),pair_weight,result)
        else if (time(j) < time(i) .and. status(j)/=0) then
          call add_pair(risk(j),risk(i),pair_weight,result)
        end if
      end do
    end do
    denom = result%concordant + result%discordant + result%tied_risk
    if (denom > 0.0_dp) then
      result%cindex = (result%concordant + 0.5_dp*result%tied_risk)/denom
    end if
  end subroutine concordance_right

  subroutine add_pair(early,late,w,result)
    real(dp), intent(in) :: early, late, w
    type(concordance_result), intent(inout) :: result
    if (early > late) then
      result%concordant = result%concordant + w
    else if (early < late) then
      result%discordant = result%discordant + w
    else
      result%tied_risk = result%tied_risk + w
    end if
  end subroutine add_pair

  pure logical function same_time(a,b) result(equal)
    real(dp), intent(in) :: a,b
    real(dp) :: scale
    scale = max(1.0_dp,abs(a),abs(b))
    equal = abs(a-b) <= 8.0_dp*epsilon(1.0_dp)*scale
  end function same_time

  subroutine sort_real(x)
    real(dp), intent(inout) :: x(:)
    integer :: i,j
    real(dp) :: value
    do i=2,size(x)
      value=x(i)
      j=i-1
      do while(j>=1)
        if(x(j)<=value) exit
        x(j+1)=x(j)
        j=j-1
      end do
      x(j+1)=value
    end do
  end subroutine sort_real
end module survival_stats
