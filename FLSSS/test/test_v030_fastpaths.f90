program test_v030_fastpaths
  use flsss_mod
  use flsss_fast_search, only : search_md_i8_packed, search_md_i8_pat
  use flsss_search, only : search_md_i8
  implicit none

  call test_packed_search()
  call test_pat_search()
  call test_public_engines()
  call test_auto_pat_boundary()
  call test_decomposition_runners()
  print *, 'test_v030_fastpaths: PASS'

contains

  subroutine test_packed_search()
    integer(i8) :: v(10,10), target(10), me(10)
    integer :: i, j, k, idx(3)
    type(subset_solutions) :: a, b
    allocate(a%sol(0),b%sol(0))
    do j = 1, 10
      do i = 1, 10
        v(i,j) = int(mod(7*i + 11*j + i*j,31) + j, i8)
      end do
    end do
    idx = [2,6,9]
    target = 0_i8
    do k = 1, 3
      target = target + v(idx(k),:)
    end do
    me = [(int(mod(j,3)+1,i8), j=1,10)]
    call search_md_i8(v,3,target,me,a,10000,10.0_dp)
    call search_md_i8_packed(v,3,target,me,b,10000,10.0_dp)
    call require_same_sets(a,b,'packed versus v0.2 DFS')
    if (b%packed_lanes >= 10 .or. b%packed_lanes < 1) error stop 'packed lane count'
    if (trim(b%engine) /= 'packed-dfs') error stop 'packed engine tag'
  end subroutine test_packed_search

  subroutine test_pat_search()
    integer(i8) :: v(12,6), target(6), me(6)
    integer :: i, j, k, idx(4)
    type(subset_solutions) :: a, b
    allocate(a%sol(0),b%sol(0))
    do j = 1, 6
      do i = 1, 12
        v(i,j) = int(i*(j+2) + (i*i)/(15+j),i8)
      end do
    end do
    idx = [2,5,8,11]
    target = 0_i8
    do k = 1, 4
      target = target + v(idx(k),:)
    end do
    me = 0_i8
    call search_md_i8(v,4,target,me,a,10000,10.0_dp)
    call search_md_i8_pat(v,4,target,me,b,10000,10.0_dp)
    call require_same_sets(a,b,'PAT versus v0.2 DFS')
    if (b%bound_states <= 0_i8) error stop 'PAT state counter'
    if (trim(b%engine) /= 'pat-packed') error stop 'PAT engine tag'
  end subroutine test_pat_search

  subroutine test_public_engines()
    real(dp) :: v(10,8), target(8), me(8)
    integer :: i, j, k, idx(3)
    type(integerized_search_result) :: a, b, c
    do j = 1, 8
      do i = 1, 10
        v(i,j) = real(mod(3*i + 5*j + i*j,19),dp) / 4.0_dp - 2.0_dp
      end do
    end do
    idx = [2,5,9]
    target = 0.0_dp
    do k = 1, 3
      target = target + v(idx(k),:)
    end do
    me = 0.30_dp
    a = mflsss_par_integerized(3,v,target,me,solution_need=10000,engine='dfs')
    b = mflsss_par_integerized(3,v,target,me,solution_need=10000,engine='packed')
    c = mflsss_par_integerized(3,v,target,me,solution_need=10000,engine='auto')
    call require_same_sets(a%solution,b%solution,'public packed engine')
    call require_same_sets(a%solution,c%solution,'public auto engine')
    if (b%solution%packed_lanes >= 8) error stop 'public API did not pack dimensions'
  end subroutine test_public_engines

  subroutine test_auto_pat_boundary()
    real(dp)::v(20,4),target(4),me(4)
    integer::i,j,k,idx(4)
    type(integerized_search_result)::r
    do j=1,4
      do i=1,20
        v(i,j)=real(i*(j+1),dp)
      end do
    end do
    idx=[17,18,19,20]
    target=0.0_dp
    do k=1,4
      target=target+v(idx(k),:)
    end do
    me=0.05_dp
    r=mflsss_par_integerized(4,v,target,me,solution_need=1,engine='auto')
    if(r%solution%size()/=1) error stop 'auto PAT boundary solution'
    if(trim(r%solution%engine)/='mpat-packed') error stop 'auto mPAT boundary engine'
  end subroutine test_auto_pat_boundary

  subroutine test_decomposition_runners()
    real(dp) :: v(9,2), target(2), me(2)
    character(len=8) :: av(7,1), at(1)
    integer :: i
    type(subset_solutions) :: direct, serial, parallel, adirect, aserial, aparallel
    type(mflsss_decomposition) :: dec
    type(arb_flsss_decomposition) :: adec

    do i = 1, 9
      v(i,1) = real(i,dp)
      v(i,2) = real(i*i,dp)
    end do
    target = v(2,:) + v(4,:) + v(7,:)
    me = 0.0_dp
    direct = mflsss_par(3,v,target,me,solution_need=10000)
    dec = decompose_mflsss(3,v,target,me,approx_ninstance=4)
    serial = mflsss_decomp_run(dec,10000,10.0_dp,.false.,2)
    parallel = mflsss_decomp_run(dec,10000,10.0_dp,.true.,2)
    call require_same_sets(direct,serial,'serial decomposition')
    call require_same_sets(direct,parallel,'parallel decomposition')
    if (serial%partitions_run /= 4_i8) error stop 'serial partition count'
    if (parallel%partitions_run /= 4_i8) error stop 'parallel partition count'

    av(:,1) = ['1       ','2       ','3       ','4       ','5       ','6       ','7       ']
    at(1) = '9'
    adirect = arb_flsss(2,av,at,solution_need=10000)
    adec = decompose_arb_flsss(2,av,at,approx_ninstance=3)
    aserial = arb_flsss_decomp_run(adec,10000,10.0_dp,.false.,2)
    aparallel = arb_flsss_decomp_run(adec,10000,10.0_dp,.true.,2)
    call require_same_sets(adirect,aserial,'arbitrary serial decomposition')
    call require_same_sets(adirect,aparallel,'arbitrary parallel decomposition')
  end subroutine test_decomposition_runners

  subroutine require_same_sets(a,b,label)
    type(subset_solutions), intent(in) :: a,b
    character(len=*), intent(in) :: label
    integer :: i
    if (a%size() /= b%size()) then
      print *, trim(label), ': sizes ', a%size(), b%size()
      error stop 'solution-set size mismatch'
    end if
    do i = 1, a%size()
      if (.not. has_solution(b,a%sol(i)%idx)) then
        print *, trim(label), ': missing solution ', a%sol(i)%idx
        error stop 'solution-set mismatch'
      end if
    end do
  end subroutine require_same_sets

  logical function has_solution(r,idx) result(found)
    type(subset_solutions), intent(in) :: r
    integer, intent(in) :: idx(:)
    integer :: i
    found = .false.
    do i = 1, r%size()
      if (size(r%sol(i)%idx) /= size(idx)) cycle
      if (all(r%sol(i)%idx == idx)) then
        found = .true.
        return
      end if
    end do
  end function has_solution

end program test_v030_fastpaths
