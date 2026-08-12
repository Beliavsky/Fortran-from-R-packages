module flsss_arbitrary
  use flsss_kinds, only : dp, i8
  use flsss_types, only : subset_solutions, arb_flsss_object, arb_flsss_decomposition, ksum_table
  use flsss_util, only : timer_type, append_solution
  use flsss_bigint, only : big_add, scaled_integer, canonical_bigint
  use flsss_parallel, only : openmp_enabled
  implicit none
  private
  public :: arb_flsss, decompose_arb_flsss, arb_flsss_obj_run, arb_flsss_decomp_run, build_ksum_hash

contains

  function arb_flsss(len, v, target, solution_need, tlimit, given_ksum) result(r)
    integer, intent(in) :: len
    character(len=*), intent(in) :: v(:,:), target(:)
    integer, intent(in), optional :: solution_need
    real(dp), intent(in), optional :: tlimit
    type(ksum_table), intent(in), optional :: given_ksum
    type(subset_solutions) :: r
    logical :: can_hash

    can_hash = .false.
    if (present(given_ksum)) then
      if (given_ksum%k < 0) error stop "arb_flsss: invalid k-sum table"
      can_hash = hash_table_compatible(given_ksum, len, v, target)
    end if
    if (can_hash) then
      call arb_search_hash(len, v, target, given_ksum, r, solution_need, tlimit)
    else
      call arb_search(len, v, target, r, solution_need, tlimit, 0)
    end if
  end function arb_flsss

  function decompose_arb_flsss(len, v, target, approx_ninstance) result(dec)
    integer, intent(in) :: len
    character(len=*), intent(in) :: v(:,:), target(:)
    integer, intent(in), optional :: approx_ninstance
    type(arb_flsss_decomposition) :: dec
    integer :: n, m, i, ng, lo, hi
    n = size(v,1)
    m = max(0, n-len+1)
    ng = m
    if (present(approx_ninstance)) then
      if (approx_ninstance < 1) error stop "decompose_arb_flsss: approx_ninstance must be positive"
      ng = min(m, approx_ninstance)
    end if
    allocate(dec%object(ng), dec%solutions_found%sol(0))
    do i = 1, ng
      lo = 1 + (i-1)*m/max(1,ng)
      hi = i*m/max(1,ng)
      dec%object(i)%len = len
      dec%object(i)%v = v
      dec%object(i)%target = target
      dec%object(i)%prefix_lo = lo
      dec%object(i)%prefix_hi = hi
      if (lo == hi) dec%object(i)%prefix = lo
    end do
  end function decompose_arb_flsss

  function arb_flsss_obj_run(x, solution_need, tlimit) result(r)
    type(arb_flsss_object), intent(in) :: x
    integer, intent(in), optional :: solution_need
    real(dp), intent(in), optional :: tlimit
    type(subset_solutions) :: r, tmp
    type(timer_type) :: timer
    integer :: need, p, plo, phi
    real(dp) :: lim, remain
    need = 1
    if (present(solution_need)) need = max(1, solution_need)
    lim = huge(1.0_dp)
    if (present(tlimit)) lim = tlimit
    call timer%start(lim)
    allocate(r%sol(0))
    plo=x%prefix_lo; phi=x%prefix_hi
    if (plo <= 0 .or. phi < plo) then
      if (x%prefix > 0) then
        plo=x%prefix; phi=x%prefix
      else
        plo=1; phi=max(0,size(x%v,1)-x%len+1)
      end if
    end if
    do p=plo,phi
      remain=max(0.0_dp,lim-timer%elapsed())
      if (remain <= 0.0_dp) then
        r%timed_out=.true.
        exit
      end if
      call arb_search(x%len,x%v,x%target,tmp,need-r%size(),remain,p)
      call merge_subset_results(r,tmp,need)
      if (r%timed_out .or. r%size() >= need) exit
    end do
  end function arb_flsss_obj_run

  function arb_flsss_decomp_run(dec,solution_need,tlimit,parallel,max_threads) result(r)
    type(arb_flsss_decomposition),intent(in)::dec
    integer,intent(in),optional::solution_need,max_threads
    real(dp),intent(in),optional::tlimit
    logical,intent(in),optional::parallel
    type(subset_solutions)::r
    type(subset_solutions),allocatable::part(:)
    integer::i,need,nthreads
    real(dp)::lim,remain
    logical::usepar
    type(timer_type)::timer
    need=1;if(present(solution_need))need=max(1,solution_need)
    lim=huge(1.0_dp);if(present(tlimit))lim=tlimit
    usepar=.false.;if(present(parallel))usepar=parallel
    nthreads=max(1,size(dec%object));if(present(max_threads))nthreads=max(1,max_threads)
    allocate(part(size(dec%object)),r%sol(0))
    if(usepar) then
!$omp parallel do schedule(dynamic) num_threads(nthreads)
      do i=1,size(dec%object)
        part(i)=arb_flsss_obj_run(dec%object(i),need,lim)
      end do
!$omp end parallel do
    else
      call timer%start(lim)
      do i=1,size(dec%object)
        remain=max(0.0_dp,lim-timer%elapsed())
        if(remain<=0.0_dp) then
          part(i)%timed_out=.true.
          cycle
        end if
        part(i)=arb_flsss_obj_run(dec%object(i),need,remain)
      end do
    end if
    do i=1,size(part)
      call merge_subset_results(r,part(i),need)
      if(r%size()>=need)exit
    end do
    r%partitions_run=int(size(part),i8)
    if(usepar .and. openmp_enabled()) then
      r%engine='arb-decomp-omp'
    else
      r%engine='arb-decomp-serial'
    end if
  end function arb_flsss_decomp_run

  function build_ksum_hash(k, v, max_entries) result(tab)
    integer, intent(in) :: k
    character(len=*), intent(in) :: v(:,:)
    integer, intent(in), optional :: max_entries
    type(ksum_table) :: tab
    integer(i8) :: nc
    integer :: n, d, cap, row, maxlen, i, j
    integer, allocatable :: chosen(:)
    character(len=:), allocatable :: sv(:,:), sums(:)

    n = size(v,1)
    d = size(v,2)
    if (k < 0 .or. k > n) error stop "build_ksum_hash: invalid k"
    tab%k = k
    nc = nchoosek_count(n,k)
    cap = 1000000
    if (present(max_entries)) cap = max_entries
    if (nc > cap) error stop "build_ksum_hash: table exceeds max_entries"
    maxlen = 32
    do j = 1, d
      do i = 1, n
        maxlen = max(maxlen, len_trim(v(i,j)) * max(1,k) + 8)
      end do
    end do
    allocate(character(len=maxlen) :: sv(n,d), sums(d), tab%sum(int(nc),d))
    allocate(tab%frac_digits(d))
    call table_scale_problem(v, sv, tab%frac_digits)
    allocate(tab%index(int(nc),k), chosen(k), tab%hash(int(nc)), tab%order(int(nc)))
    row = 0
    call enumerate(1,1)
    tab%order = [(i, i=1,int(nc))]
    call sort_hash_order(tab%hash, tab%order)

  contains
    recursive subroutine enumerate(pos,start)
      integer, intent(in) :: pos,start
      integer :: q, hi, j, u
      character(len=maxlen) :: s
      if (pos > k) then
        row = row + 1
        tab%index(row,:) = chosen
        do j = 1, d
          s = '0'
          do u = 1, k
            s = big_add(trim(s), trim(sv(chosen(u),j)))
          end do
          sums(j) = canonical_bigint(trim(s))
          tab%sum(row,j) = sums(j)
        end do
        tab%hash(row) = tuple_hash(sums)
        return
      end if
      hi = n - (k-pos)
      do q = start, hi
        chosen(pos) = q
        call enumerate(pos+1,q+1)
      end do
    end subroutine enumerate
  end function build_ksum_hash

  subroutine arb_search_hash(len, v, target, tab, result, solution_need, tlimit)
    integer, intent(in) :: len
    character(len=*), intent(in) :: v(:,:), target(:)
    type(ksum_table), intent(in) :: tab
    type(subset_solutions), intent(out) :: result
    integer, intent(in), optional :: solution_need
    real(dp), intent(in), optional :: tlimit
    integer :: n, d, leftk, need, maxlen, i, j
    integer, allocatable :: chosen(:)
    character(len=:), allocatable :: sv(:,:), st(:), sums(:), wanted(:)
    type(timer_type) :: timer
    real(dp) :: lim

    n = size(v,1)
    d = size(v,2)
    leftk = len - tab%k
    need = 1
    if (present(solution_need)) need = max(1, solution_need)
    lim = huge(1.0_dp)
    if (present(tlimit)) lim = tlimit
    maxlen = max(len_trim(tab%sum(1,1)), len_trim(target(1))) + 32
    do j = 1, d
      maxlen = max(maxlen, len_trim(target(j)) + 32)
      do i = 1, n
        maxlen = max(maxlen, len_trim(v(i,j))*max(2,len) + 32)
      end do
    end do
    allocate(character(len=maxlen) :: sv(n,d), st(d), sums(d), wanted(d))
    call scale_problem_with_digits(v, target, tab%frac_digits, sv, st)
    sums = '0'
    allocate(chosen(len), result%sol(0))
    call timer%start(lim)
    if (leftk == 0) then
      call probe_table(0)
    else
      call dfs_left(1,1)
    end if

  contains
    recursive subroutine dfs_left(pos,start)
      integer, intent(in) :: pos,start
      integer :: q, hi, u
      character(len=:), allocatable :: old(:)
      result%nodes = result%nodes + 1_i8
      if (iand(result%nodes, 511_i8) == 0_i8) then
        if (timer%expired()) then
          result%timed_out = .true.
          return
        end if
      end if
      if (size(result%sol) >= need) return
      if (pos > leftk) then
        call probe_table(chosen(leftk))
        return
      end if
      ! Leave at least tab%k entries after the left prefix.  This unique split
      ! makes the last k indices come from the lookup table and avoids duplicates.
      hi = n - (leftk-pos) - tab%k
      allocate(character(len=maxlen) :: old(d))
      old = sums
      do q = start, hi
        chosen(pos) = q
        do u = 1, d
          sums(u) = big_add(trim(old(u)), trim(sv(q,u)))
        end do
        call dfs_left(pos+1,q+1)
        if (result%timed_out .or. size(result%sol) >= need) exit
      end do
      sums = old
    end subroutine dfs_left

    subroutine probe_table(last_left)
      integer, intent(in) :: last_left
      integer(i8) :: h
      integer :: lo, hi, mid, p, row, u
      logical :: match
      do u = 1, d
        wanted(u) = big_add(trim(st(u)), negate_big(trim(sums(u))))
        wanted(u) = canonical_bigint(trim(wanted(u)))
      end do
      h = tuple_hash(wanted)
      result%hash_lookups = result%hash_lookups + 1_i8
      lo = 1
      hi = size(tab%order)
      do while (lo <= hi)
        mid = lo + (hi-lo)/2
        if (tab%hash(tab%order(mid)) < h) then
          lo = mid + 1
        else
          hi = mid - 1
        end if
      end do
      p = lo
      do while (p <= size(tab%order))
        row = tab%order(p)
        if (tab%hash(row) /= h) exit
        result%hash_candidates = result%hash_candidates + 1_i8
        if (tab%k > 0) then
          if (tab%index(row,1) <= last_left) then
            p = p + 1
            cycle
          end if
        end if
        match = .true.
        do u = 1, d
          if (canonical_bigint(trim(tab%sum(row,u))) /= canonical_bigint(trim(wanted(u)))) then
            match = .false.
            exit
          end if
        end do
        if (match) then
          if (leftk > 0) chosen(1:leftk) = chosen(1:leftk)
          if (tab%k > 0) chosen(leftk+1:len) = tab%index(row,:)
          call append_solution(result, chosen)
          if (size(result%sol) >= need) return
        end if
        p = p + 1
      end do
    end subroutine probe_table
  end subroutine arb_search_hash

  subroutine arb_search(len, v, target, result, solution_need, tlimit, prefix)
    integer, intent(in) :: len, prefix
    character(len=*), intent(in) :: v(:,:), target(:)
    type(subset_solutions), intent(out) :: result
    integer, intent(in), optional :: solution_need
    real(dp), intent(in), optional :: tlimit
    integer :: n, d, need, i, j, maxlen
    integer, allocatable :: chosen(:)
    character(len=:), allocatable :: sv(:,:), st(:), sums(:)
    type(timer_type) :: timer
    real(dp) :: lim

    n = size(v,1)
    d = size(v,2)
    if (size(target) /= d) error stop "arb_search: target dimension mismatch"
    need = 1
    if (present(solution_need)) need = max(1,solution_need)
    lim = huge(1.0_dp)
    if (present(tlimit)) lim = tlimit
    maxlen = 64
    do j = 1, d
      do i = 1, n
        maxlen = max(maxlen, len_trim(v(i,j)) * max(2,len) + 32)
      end do
      maxlen = max(maxlen, len_trim(target(j)) * max(2,len) + 32)
    end do
    allocate(character(len=maxlen) :: sv(n,d), st(d), sums(d))
    call scale_problem(v, target, sv, st, .true.)
    sums = '0'
    allocate(chosen(len), result%sol(0))
    call timer%start(lim)
    if (prefix > 0) then
      if (prefix > n-len+1) return
      chosen(1) = prefix
      sums = sv(prefix,:)
      call dfs(2,prefix+1)
    else
      call dfs(1,1)
    end if

  contains
    recursive subroutine dfs(pos,start)
      integer, intent(in) :: pos,start
      integer :: q, hi, u
      character(len=:), allocatable :: old(:)
      result%nodes = result%nodes + 1_i8
      if (iand(result%nodes, 511_i8) == 0_i8) then
        if (timer%expired()) then
          result%timed_out = .true.
          return
        end if
      end if
      if (size(result%sol) >= need) return
      if (pos > len) then
        do u = 1, d
          if (canonical_bigint(trim(sums(u))) /= canonical_bigint(trim(st(u)))) return
        end do
        call append_solution(result, chosen)
        return
      end if
      hi = n - (len-pos)
      allocate(character(len=maxlen) :: old(d))
      old = sums
      do q = start, hi
        chosen(pos) = q
        do u = 1, d
          sums(u) = big_add(trim(old(u)), trim(sv(q,u)))
        end do
        call dfs(pos+1,q+1)
        if (result%timed_out .or. size(result%sol) >= need) exit
      end do
      sums = old
    end subroutine dfs
  end subroutine arb_search

  logical function hash_table_compatible(tab, len, v, target) result(ok)
    type(ksum_table), intent(in) :: tab
    integer, intent(in) :: len
    character(len=*), intent(in) :: v(:,:), target(:)
    integer :: j, p, frac
    ok = .false.
    if (tab%k <= 0 .or. tab%k >= len) return
    if (.not. allocated(tab%hash) .or. .not. allocated(tab%order)) return
    if (.not. allocated(tab%frac_digits)) return
    if (size(tab%sum,2) /= size(v,2)) return
    if (size(target) /= size(v,2)) return
    if (size(tab%frac_digits) /= size(v,2)) return
    do j = 1, size(target)
      p = index(trim(target(j)),'.')
      frac = 0
      if (p > 0) frac = len_trim(target(j)) - p
      if (frac > tab%frac_digits(j)) return
    end do
    ok = .true.
  end function hash_table_compatible

  subroutine table_scale_problem(v, sv, frac_digits)
    character(len=*), intent(in) :: v(:,:)
    character(len=*), intent(out) :: sv(:,:)
    integer, intent(out) :: frac_digits(:)
    integer :: n, d, i, j, p, frac
    n = size(v,1)
    d = size(v,2)
    do j = 1, d
      frac = 0
      do i = 1, n
        p = index(trim(v(i,j)),'.')
        if (p > 0) frac = max(frac,len_trim(v(i,j))-p)
      end do
      frac_digits(j) = frac
      do i = 1, n
        sv(i,j) = scaled_integer(trim(v(i,j)),frac)
      end do
    end do
  end subroutine table_scale_problem

  subroutine scale_problem_with_digits(v, target, frac_digits, sv, st)
    character(len=*), intent(in) :: v(:,:), target(:)
    integer, intent(in) :: frac_digits(:)
    character(len=*), intent(out) :: sv(:,:), st(:)
    integer :: i, j
    do j = 1, size(v,2)
      do i = 1, size(v,1)
        sv(i,j) = scaled_integer(trim(v(i,j)),frac_digits(j))
      end do
      st(j) = scaled_integer(trim(target(j)),frac_digits(j))
    end do
  end subroutine scale_problem_with_digits

  subroutine scale_problem(v, target, sv, st, use_target)
    character(len=*), intent(in) :: v(:,:), target(:)
    character(len=*), intent(out) :: sv(:,:), st(:)
    logical, intent(in) :: use_target
    integer :: n,d,i,j,p,frac
    n=size(v,1)
    d=size(v,2)
    do j=1,d
      frac=0
      do i=1,n
        p=index(trim(v(i,j)),'.')
        if(p>0) frac=max(frac,len_trim(v(i,j))-p)
      end do
      if (use_target .and. j <= size(target)) then
        p=index(trim(target(j)),'.')
        if(p>0) frac=max(frac,len_trim(target(j))-p)
      end if
      do i=1,n
        sv(i,j)=scaled_integer(trim(v(i,j)),frac)
      end do
      if (use_target) then
        st(j)=scaled_integer(trim(target(j)),frac)
      else
        st(j)='0'
      end if
    end do
  end subroutine scale_problem

  function negate_big(s) result(r)
    character(len=*), intent(in) :: s
    character(len=:), allocatable :: r, x
    x = canonical_bigint(s)
    if (x == '0') then
      r = '0'
    else if (x(1:1) == '-') then
      r = x(2:)
    else
      r = '-'//x
    end if
  end function negate_big

  pure integer(i8) function tuple_hash(x) result(h)
    character(len=*), intent(in) :: x(:)
    integer :: i, j
    h = int(z'6A09E667F3BCC909',i8)
    do j = 1, size(x)
      do i = 1, len_trim(x(j))
        h = ieor(ishftc(h, 7), int(iachar(x(j)(i:i)),i8))
        h = ieor(h, ishft(h,-11))
      end do
      h = ieor(ishftc(h,13), int(j,i8))
    end do
  end function tuple_hash

  subroutine sort_hash_order(hash, order)
    integer(i8), intent(in) :: hash(:)
    integer, intent(inout) :: order(:)
    integer, allocatable :: tmp(:)
    allocate(tmp(size(order)))
    call merge_sort(1,size(order))
  contains
    recursive subroutine merge_sort(lo,hi)
      integer, intent(in) :: lo,hi
      integer :: mid,i,j,k
      if (lo >= hi) return
      mid = lo + (hi-lo)/2
      call merge_sort(lo,mid)
      call merge_sort(mid+1,hi)
      i=lo; j=mid+1; k=lo
      do while (i<=mid .and. j<=hi)
        if (hash(order(i)) <= hash(order(j))) then
          tmp(k)=order(i); i=i+1
        else
          tmp(k)=order(j); j=j+1
        end if
        k=k+1
      end do
      do while(i<=mid); tmp(k)=order(i); i=i+1; k=k+1; end do
      do while(j<=hi); tmp(k)=order(j); j=j+1; k=k+1; end do
      order(lo:hi)=tmp(lo:hi)
    end subroutine merge_sort
  end subroutine sort_hash_order

  subroutine merge_subset_results(a,b,need)
    type(subset_solutions), intent(inout) :: a
    type(subset_solutions), intent(in) :: b
    integer, intent(in) :: need
    integer :: i
    do i=1,b%size()
      if(a%size()>=need) exit
      call append_solution(a,b%sol(i)%idx)
    end do
    a%nodes=a%nodes+b%nodes
    a%pruned=a%pruned+b%pruned
    a%hash_lookups=a%hash_lookups+b%hash_lookups
    a%hash_candidates=a%hash_candidates+b%hash_candidates
    a%bound_states=a%bound_states+b%bound_states
    a%timed_out=a%timed_out .or. b%timed_out
  end subroutine merge_subset_results

  integer(i8) function nchoosek_count(n,k) result(c)
    integer,intent(in)::n,k
    integer::i,kk
    if(k<0 .or. k>n) then
      c=0
      return
    end if
    kk=min(k,n-k)
    c=1_i8
    do i=1,kk
      c=c*int(n-kk+i,i8)/int(i,i8)
    end do
  end function nchoosek_count

end module flsss_arbitrary
