! SPDX-License-Identifier: GPL-3.0-or-later
! Modern Fortran translation of computational code from R package adagio 0.9.2.
module adagio_discrete
  use adagio_kinds, only : dp
  use adagio_types, only : assignment_result, subset_result, knapsack_result, &
       mknapsack_result, setcover_result, change_result, binpack_result, hamiltonian_result
  use lpsolve, only : lp_result, solve_lp, lp_assign, LP_MIN, LP_MAX, LP_LE, LP_GE, LP_EQ, &
       LP_OPTIMAL, LP_SUBOPTIMAL
  implicit none
  private
  public :: assignment, subsetsum, sss_test, knapsack, mknapsack
  public :: setcover, change_making, bpp_approx, hamiltonian

contains

  function assignment(cmat, maximize) result(res)
    real(dp), intent(in) :: cmat(:,:)
    logical, intent(in), optional :: maximize
    type(assignment_result) :: res
    type(lp_result) :: lpres
    real(dp), allocatable :: smat(:,:)
    integer :: n, i, direction

    n = size(cmat, 1)
    if (size(cmat, 2) /= n .or. n <= 1) then
       res%err = 1
       allocate(res%perm(0))
       return
    end if

    direction = LP_MIN
    if (present(maximize)) then
       if (maximize) direction = LP_MAX
    end if

    allocate(smat(n,n))
    smat = 0.0_dp
    call lp_assign(cmat, lpres, direction=direction, assignment=smat)
    if (lpres%status /= LP_OPTIMAL .and. lpres%status /= LP_SUBOPTIMAL) then
       res%err = lpres%status
       allocate(res%perm(0))
       return
    end if

    allocate(res%perm(n))
    do i = 1, n
       res%perm(i) = maxloc(smat(i,:), dim=1)
    end do
    res%value = sum(cmat * smat)
    res%err = 0
  end function assignment

  function sss_test(s, target) result(best)
    integer, intent(in) :: s(:), target
    integer :: best, i, v
    logical, allocatable :: reach(:)
    allocate(reach(0:target)); reach=.false.; reach(0)=.true.
    do i=1,size(s)
       do v=target,s(i),-1
          reach(v)=reach(v) .or. reach(v-s(i))
       end do
       if(reach(target)) exit
    end do
    best=0
    do v=target,0,-1
       if(reach(v)) then; best=v; exit; end if
    end do
  end function sss_test

  function subsetsum(s, target, method) result(res)
    integer, intent(in) :: s(:), target
    character(len=*), intent(in), optional :: method
    type(subset_result) :: res
    character(len=16) :: meth
    logical, allocatable :: reach(:), x(:)
    integer, allocatable :: from_item(:)
    integer :: i, v, t, k, n, j, nsel
    meth='greedy'; if(present(method)) meth=adjustl(method)
    n=size(s); allocate(x(n)); x=.false.
    if(index(meth,'dynamic')==1) then
       allocate(reach(0:target),from_item(0:target)); reach=.false.; from_item=0; reach(0)=.true.
       do k=1,n
          do v=target,s(k),-1
             if(.not.reach(v) .and. reach(v-s(k))) then
                reach(v)=.true.; from_item(v)=k
             end if
          end do
          if(reach(target)) exit
       end do
       t=target
       if(.not.reach(t)) then
          do v=target,0,-1; if(reach(v)) then; t=v; exit; end if; end do
       end if
       res%val=t; res%found=(t==target)
       do while(t>0 .and. from_item(t)>0)
          k=from_item(t); x(k)=.true.; t=t-s(k)
       end do
    else
       ! Source-compatible repeated prefix reachability reconstruction.
       t=target
       do while(t>0)
          allocate(reach(0:t)); reach=.false.; reach(0)=.true.; j=0
          do i=1,n
             do v=t,s(i),-1
                reach(v)=reach(v) .or. reach(v-s(i))
             end do
             if(reach(t)) then; j=i; exit; end if
          end do
          deallocate(reach)
          if(j==0) exit
          x(j)=.true.; t=t-s(j)
       end do
       res%found=(t==0); res%val=target-t
    end if
    nsel=count(x); allocate(res%inds(nsel)); k=0
    do i=1,n; if(x(i)) then; k=k+1; res%inds(k)=i; end if; end do
  end function subsetsum

  function knapsack(w, p, cap) result(res)
    integer, intent(in) :: w(:), cap
    real(dp), intent(in) :: p(:)
    type(knapsack_result) :: res
    real(dp), allocatable :: dpv(:,:), g(:)
    logical, allocatable :: take(:)
    integer :: n, k, c, j, nsel
    real(dp) :: f
    n=size(w); allocate(dpv(0:cap,n),g(0:cap)); g=0.0_dp
    do k=1,n
       dpv(:,k)=g
       do c=cap,w(k),-1
          g(c)=max(g(c),g(c-w(k))+p(k))
       end do
    end do
    allocate(take(n)); take=.false.; c=cap; f=g(cap)
    do k=n,1,-1
       if(dpv(c,k)<f) then
          take(k)=.true.; c=c-w(k); f=dpv(c,k)
       end if
    end do
    nsel=count(take); allocate(res%indices(nsel)); j=0
    do k=1,n; if(take(k)) then; j=j+1; res%indices(j)=k; end if; end do
    res%capacity=sum(w,mask=take); res%profit=sum(p,mask=take)
  end function knapsack

  function change_making(items, value) result(res)
    integer, intent(in) :: items(:), value
    type(change_result) :: res
    type(lp_result) :: lpres
    real(dp), allocatable :: obj(:), amat(:,:), rhs(:)
    integer, allocatable :: sense(:), ints(:)
    integer :: n, i

    n = size(items)
    if (n == 0) then
       allocate(res%solution(0))
       res%feasible = (value == 0)
       return
    end if

    allocate(obj(n), amat(1,n), rhs(1), sense(1), ints(n))
    obj = 1.0_dp
    amat(1,:) = real(items, dp)
    rhs(1) = real(value, dp)
    sense(1) = LP_EQ
    ints = [(i, i=1,n)]
    call solve_lp(LP_MIN, obj, amat, sense, rhs, lpres, integer_variables=ints)

    res%feasible = lpres%status == LP_OPTIMAL .or. lpres%status == LP_SUBOPTIMAL
    if (.not. res%feasible) then
       allocate(res%solution(0))
       res%count = 0
       return
    end if

    res%count = nint(lpres%objective)
    if (res%count == 0) then
       allocate(res%solution(0))
    else
       allocate(res%solution(n))
       res%solution = nint(lpres%solution)
    end if
  end function change_making

  function bpp_approx(s, cap, method) result(res)
    real(dp), intent(in) :: s(:), cap
    character(len=*), intent(in), optional :: method
    type(binpack_result) :: res
    real(dp), allocatable :: b(:)
    character(len=16) :: meth
    integer :: n, i, j, k, m
    meth='firstfit'; if(present(method)) meth=adjustl(method)
    n=size(s); allocate(b(n),res%xbins(n)); b=cap; res%xbins=0
    do i=1,n
       j=0
       if(index(meth,'first')==1) then
          do k=1,n; if(b(k)>=s(i)) then; j=k; exit; end if; end do
       else if(index(meth,'best')==1) then
          do k=1,n
             if(b(k)>=s(i)) then
                if (j == 0) then
                   j = k
                else if (b(k) < b(j)) then
                   j = k
                end if
             end if
          end do
       else
          do k=1,n
             if(b(k)>=s(i) .and. b(k)<cap) then
                if (j == 0) then
                   j = k
                else if (b(k) > b(j)) then
                   j = k
                end if
             end if
          end do
          if(j==0) then
             do k=1,n
                if(b(k)>=s(i)) then
                   if (j == 0) then
                   j = k
                else if (b(k) > b(j)) then
                   j = k
                end if
                end if
             end do
          end if
       end if
       if(j==0) cycle
       res%xbins(i)=j; b(j)=b(j)-s(i)
    end do
    m=maxval(res%xbins); res%nbins=m
    allocate(res%sbins(m)); res%sbins=cap-b(1:m)
    if(m>0) res%filled=sum(s)/(real(m,dp)*cap)
  end function bpp_approx

  function setcover(sets, weights) result(res)
    integer, intent(in) :: sets(:,:)
    real(dp), intent(in), optional :: weights(:)
    type(setcover_result) :: res
    type(lp_result) :: lpres
    real(dp), allocatable :: obj(:), amat(:,:), rhs(:)
    integer, allocatable :: sense(:), bins(:)
    integer :: n, m, i, k

    n = size(sets,1)
    m = size(sets,2)
    if (n == 0 .or. m == 0 .or. any((sets /= 0) .and. (sets /= 1))) then
       allocate(res%sets(0))
       return
    end if
    if (any(sum(sets, dim=1) == 0)) then
       allocate(res%sets(0))
       return
    end if

    allocate(obj(n), amat(m,n), rhs(m), sense(m), bins(n))
    obj = 1.0_dp
    if (present(weights)) then
       if (size(weights) /= n) then
          allocate(res%sets(0))
          return
       end if
       obj = weights
    end if
    amat = transpose(real(sets, dp))
    rhs = 1.0_dp
    sense = LP_GE
    bins = [(i, i=1,n)]
    call solve_lp(LP_MIN, obj, amat, sense, rhs, lpres, binary_variables=bins)

    res%feasible = lpres%status == LP_OPTIMAL .or. lpres%status == LP_SUBOPTIMAL
    if (.not. res%feasible) then
       allocate(res%sets(0))
       return
    end if
    res%objective = lpres%objective
    k = count(lpres%solution > 0.5_dp)
    allocate(res%sets(k))
    k = 0
    do i = 1, n
       if (lpres%solution(i) > 0.5_dp) then
          k = k + 1
          res%sets(k) = i
       end if
    end do
  end function setcover

  function mknapsack(w,p,cap) result(res)
    real(dp), intent(in) :: w(:), p(:), cap(:)
    type(mknapsack_result) :: res
    type(lp_result) :: lpres
    real(dp), allocatable :: obj(:), amat(:,:), rhs(:)
    integer, allocatable :: sense(:), bins(:)
    integer :: n, m, nvar, ncon, i, k, j, col

    n = size(w)
    m = size(cap)
    res%bs = 0
    if (size(p) /= n .or. n == 0 .or. m == 0) then
       allocate(res%ksack(0))
       return
    end if

    if (n == 1) then
       allocate(res%ksack(1))
       res%ksack = 0
       do k = 1, m
          if (cap(k) >= w(1)) then
             res%ksack(1) = k
             res%value = p(1)
             return
          end if
       end do
       return
    end if

    if (m == 1) then
       allocate(obj(n), amat(1,n), rhs(1), sense(1), bins(n))
       obj = p
       amat(1,:) = w
       rhs(1) = cap(1)
       sense(1) = LP_LE
       bins = [(i, i=1,n)]
       call solve_lp(LP_MAX, obj, amat, sense, rhs, lpres, binary_variables=bins)
       allocate(res%ksack(n))
       res%ksack = 0
       if (lpres%status == LP_OPTIMAL .or. lpres%status == LP_SUBOPTIMAL) then
          res%ksack = nint(lpres%solution)
          res%value = lpres%objective
       end if
       return
    end if

    nvar = n * m
    ncon = m + n
    allocate(obj(nvar), amat(ncon,nvar), rhs(ncon), sense(ncon), bins(nvar))
    amat = 0.0_dp
    do k = 1, m
       do i = 1, n
          col = (k-1)*n + i
          obj(col) = p(i)
          amat(k,col) = w(i)
          amat(m+i,col) = 1.0_dp
       end do
    end do
    rhs(1:m) = cap
    rhs(m+1:ncon) = 1.0_dp
    sense = LP_LE
    bins = [(j, j=1,nvar)]
    call solve_lp(LP_MAX, obj, amat, sense, rhs, lpres, binary_variables=bins)

    allocate(res%ksack(n))
    res%ksack = 0
    if (lpres%status /= LP_OPTIMAL .and. lpres%status /= LP_SUBOPTIMAL) return
    res%value = lpres%objective
    do k = 1, m
       do i = 1, n
          col = (k-1)*n + i
          if (lpres%solution(col) > 0.5_dp) res%ksack(i) = k
       end do
    end do
  end function mknapsack

  function hamiltonian(edges,start,cycle) result(res)
    integer, intent(in) :: edges(:)
    integer, intent(in), optional :: start
    logical, intent(in), optional :: cycle
    type(hamiltonian_result) :: res
    logical, allocatable :: adj(:,:),used(:)
    integer, allocatable :: path(:)
    integer :: n, i, s
    logical :: cyc, ok
    if(mod(size(edges),2)/=0) then; allocate(res%path(0)); return; end if
    n=maxval(edges); s=1; if(present(start)) s=start; cyc=.true.; if(present(cycle)) cyc=cycle
    allocate(adj(n,n),used(n),path(n)); adj=.false.; used=.false.; path=0
    do i=1,size(edges),2
       adj(edges(i),edges(i+1))=.true.; adj(edges(i+1),edges(i))=.true.
    end do
    path(1)=s; used(s)=.true.; ok=search(2)
    res%found=ok
    if(ok) then; allocate(res%path(n)); res%path=path; else; allocate(res%path(0)); end if
  contains
    recursive logical function search(k) result(found)
      integer, intent(in) :: k
      integer :: v
      found=.false.
      if(k>n) then; found=(.not.cyc .or. adj(path(n),path(1))); return; end if
      do v=1,n
         if(adj(path(k-1),v) .and. .not.used(v)) then
            path(k)=v; used(v)=.true.
            if(search(k+1)) then; found=.true.; return; end if
            used(v)=.false.; path(k)=0
         end if
      end do
    end function search
  end function hamiltonian

end module adagio_discrete
