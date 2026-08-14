! Pareto ranking and diversity utilities translated from rmoo.
! rmoo original package: GPL-2.0-or-later.
module rmoo_pareto
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf
  use ga_kinds, only : dp
  implicit none
  private
  public :: dominates, non_dominated_sort, crowding_distance
  public :: generational_distance, gd_plus, igd, igd_plus
  public :: calc_norm_pref_distance, perpendicular_similarity
contains

  pure logical function dominates(a, b) result(ans)
    real(dp), intent(in) :: a(:), b(:)
    if (size(a) /= size(b)) error stop "dominates: dimension mismatch"
    ans = all(a <= b) .and. any(a < b)
  end function dominates

  subroutine non_dominated_sort(fitness, rank)
    real(dp), intent(in) :: fitness(:,:)
    integer, intent(out) :: rank(size(fitness,1))
    integer :: n, i, j, r, nq, nnew
    integer, allocatable :: n_dom(:), q(:), qnew(:)
    logical, allocatable :: dom(:,:)

    n = size(fitness,1)
    rank = 0
    if (n == 0) return
    allocate(n_dom(n), dom(n,n), q(n), qnew(n))
    n_dom = 0
    dom = .false.
    do i = 1, n - 1
      do j = i + 1, n
        if (dominates(fitness(i,:), fitness(j,:))) then
          dom(i,j) = .true.
          n_dom(j) = n_dom(j) + 1
        else if (dominates(fitness(j,:), fitness(i,:))) then
          dom(j,i) = .true.
          n_dom(i) = n_dom(i) + 1
        end if
      end do
    end do

    nq = 0
    do i = 1, n
      if (n_dom(i) == 0) then
        nq = nq + 1
        q(nq) = i
      end if
    end do
    r = 1
    do while (nq > 0)
      do i = 1, nq
        rank(q(i)) = r
      end do
      nnew = 0
      do i = 1, nq
        do j = 1, n
          if (dom(q(i),j)) then
            n_dom(j) = n_dom(j) - 1
            if (n_dom(j) == 0) then
              nnew = nnew + 1
              qnew(nnew) = j
            end if
          end if
        end do
      end do
      if (nnew > 0) q(1:nnew) = qnew(1:nnew)
      nq = nnew
      r = r + 1
    end do
  end subroutine non_dominated_sort

  subroutine crowding_distance(fitness, rank, crowding)
    real(dp), intent(in) :: fitness(:,:)
    integer, intent(in) :: rank(:)
    real(dp), intent(out) :: crowding(size(rank))
    integer :: n, m, max_rank, r, nf, j, k
    integer, allocatable :: idx(:), ord(:)
    real(dp) :: rangev

    n = size(fitness,1); m = size(fitness,2)
    if (size(rank) /= n) error stop "crowding_distance: rank size mismatch"
    crowding = 0.0_dp
    if (n == 0) return
    max_rank = maxval(rank)
    allocate(idx(n), ord(n))
    do r = 1, max_rank
      nf = 0
      do k = 1, n
        if (rank(k) == r) then
          nf = nf + 1
          idx(nf) = k
        end if
      end do
      if (nf == 0) cycle
      if (nf <= 2) then
        crowding(idx(1:nf)) = ieee_value(1.0_dp,ieee_positive_inf)
        cycle
      end if
      do j = 1, m
        call argsort_values(fitness(:,j), idx(1:nf), ord(1:nf))
        crowding(ord(1)) = ieee_value(1.0_dp,ieee_positive_inf)
        crowding(ord(nf)) = ieee_value(1.0_dp,ieee_positive_inf)
        rangev = fitness(ord(nf),j) - fitness(ord(1),j)
        if (abs(rangev) <= tiny(1.0_dp)) cycle
        do k = 2, nf - 1
          if (crowding(ord(k)) < ieee_value(1.0_dp,ieee_positive_inf)) then
            crowding(ord(k)) = crowding(ord(k)) + &
              (fitness(ord(k+1),j) - fitness(ord(k-1),j))/abs(rangev)
          end if
        end do
      end do
    end do
  contains
    subroutine argsort_values(x, subset, out)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: subset(:)
      integer, intent(out) :: out(size(subset))
      integer :: ii, jj, t
      out = subset
      do ii = 2, size(out)
        t = out(ii); jj = ii - 1
        do while (jj >= 1)
          if (x(out(jj)) <= x(t)) exit
          out(jj+1) = out(jj)
          jj = jj - 1
        end do
        out(jj+1) = t
      end do
    end subroutine argsort_values
  end subroutine crowding_distance

  real(dp) function generational_distance(front, true_front, p, inverted, plus) result(gd)
    real(dp), intent(in) :: front(:,:), true_front(:,:)
    real(dp), intent(in), optional :: p
    logical, intent(in), optional :: inverted, plus
    real(dp) :: pv
    logical :: inv, pl
    pv = 1.0_dp; if (present(p)) pv = p
    inv = .false.; if (present(inverted)) inv = inverted
    pl = .false.; if (present(plus)) pl = plus
    if (pv <= 0.0_dp) error stop "generational_distance: p must be positive"
    if (size(front,2) /= size(true_front,2)) error stop "generational_distance: dimension mismatch"
    if (inv) then
      gd = gd_core(true_front, front, pv, pl, .true.)
    else
      gd = gd_core(front, true_front, pv, pl, .false.)
    end if
  contains
    real(dp) function gd_core(a, b, pp, ppplus, was_inverted) result(v)
      real(dp), intent(in) :: a(:,:), b(:,:), pp
      logical, intent(in) :: ppplus, was_inverted
      integer :: ia, ib
      real(dp) :: ss, dd, bb
      if (size(a,1) == 0 .or. size(b,1) == 0) then
        v = huge(1.0_dp); return
      end if
      ss = 0.0_dp
      do ia = 1, size(a,1)
        bb = huge(1.0_dp)
        do ib = 1, size(b,1)
          if (ppplus) then
            if (was_inverted) then
              dd = sqrt(sum(max(b(ib,:) - a(ia,:), 0.0_dp)**2))
            else
              dd = sqrt(sum(max(a(ia,:) - b(ib,:), 0.0_dp)**2))
            end if
          else
            dd = sqrt(sum((b(ib,:) - a(ia,:))**2))
          end if
          bb = min(bb, dd)
        end do
        ss = ss + bb**pp
      end do
      v = ss**(1.0_dp/pp)/real(size(a,1),dp)
    end function gd_core
  end function generational_distance

  real(dp) function gd_plus(front, true_front, p) result(v)
    real(dp), intent(in) :: front(:,:), true_front(:,:)
    real(dp), intent(in), optional :: p
    real(dp) :: pv
    pv=1.0_dp; if(present(p))pv=p
    v=generational_distance(front,true_front,pv,.false.,.true.)
  end function gd_plus

  real(dp) function igd(front, true_front, p) result(v)
    real(dp), intent(in) :: front(:,:), true_front(:,:)
    real(dp), intent(in), optional :: p
    real(dp) :: pv
    pv=1.0_dp; if(present(p))pv=p
    v=generational_distance(front,true_front,pv,.true.,.false.)
  end function igd

  real(dp) function igd_plus(front, true_front, p) result(v)
    real(dp), intent(in) :: front(:,:), true_front(:,:)
    real(dp), intent(in), optional :: p
    real(dp) :: pv
    pv=1.0_dp; if(present(p))pv=p
    v=generational_distance(front,true_front,pv,.true.,.true.)
  end function igd_plus

  subroutine calc_norm_pref_distance(fitness, ref_points, weight, ideal_point, nadir_point, dist)
    real(dp), intent(in) :: fitness(:,:), ref_points(:,:), weight(:), ideal_point(:), nadir_point(:)
    real(dp), intent(out) :: dist(size(fitness,1),size(ref_points,1))
    real(dp) :: denom(size(weight)), d2
    integer :: i,j,k,p
    p=size(fitness,2)
    if(size(ref_points,2)/=p .or. size(weight)/=p .or. size(ideal_point)/=p .or. size(nadir_point)/=p) &
      error stop "calc_norm_pref_distance: dimension mismatch"
    denom=nadir_point-ideal_point
    where(abs(denom)<=tiny(1.0_dp)) denom=1.0e-12_dp
    do i=1,size(fitness,1)
      do j=1,size(ref_points,1)
        d2=0.0_dp
        do k=1,p
          d2=d2+weight(k)*(fitness(i,k)-ref_points(j,k))**2/denom(k)**2
        end do
        dist(i,j)=sqrt(max(0.0_dp,d2*real(p,dp)))
      end do
    end do
  end subroutine calc_norm_pref_distance

  subroutine perpendicular_similarity(x, dirs, cosine)
    real(dp), intent(in) :: x(:,:), dirs(:,:)
    real(dp), intent(out) :: cosine(size(x,1),size(dirs,1))
    real(dp) :: nx, nd
    integer :: i,j
    if(size(x,2)/=size(dirs,2))error stop "perpendicular_similarity: dimension mismatch"
    do i=1,size(x,1)
      nx=sqrt(sum(x(i,:)**2))
      do j=1,size(dirs,1)
        nd=sqrt(sum(dirs(j,:)**2))
        if(nx<=tiny(1.0_dp).or.nd<=tiny(1.0_dp))then
          cosine(i,j)=0.0_dp
        else
          cosine(i,j)=dot_product(x(i,:),dirs(j,:))/(nx*nd)
          cosine(i,j)=min(1.0_dp,max(-1.0_dp,cosine(i,j)))
        end if
      end do
    end do
  end subroutine perpendicular_similarity
end module rmoo_pareto
