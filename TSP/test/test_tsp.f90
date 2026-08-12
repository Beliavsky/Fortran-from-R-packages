program test_tsp
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use tsp
    implicit none

    call set_seed()
    call test_core()
    call test_heuristics()
    call test_atsp_two_opt()
    call test_transformations()
    call test_tsplib_distances()
    call test_tsplib_roundtrip()
    call test_random_two_opt()
    call test_random_insertion()
    print '(A)', 'All TSP tests passed.'

contains


    subroutine set_seed()
        integer, allocatable :: seed(:)
        integer :: n, i
        call random_seed(size=n)
        allocate(seed(n))
        do i = 1, n
            seed(i) = 1234567 + 7919*i
        end do
        call random_seed(put=seed)
    end subroutine set_seed

    subroutine check(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message
        if (.not. condition) then
            write(*,'(A)') 'FAIL: '//trim(message)
            error stop 1
        end if
    end subroutine check

    subroutine check_close(a, b, tol, message)
        real(dp), intent(in) :: a, b, tol
        character(len=*), intent(in) :: message
        call check(abs(a-b) <= tol, message)
    end subroutine check_close

    subroutine test_core()
        real(dp) :: m(4,4), inf
        integer :: order(4)
        real(dp), allocatable :: delta(:), replaced(:,:), coords(:,:), d(:,:)
        integer :: ierr

        inf = positive_infinity()
        m = reshape([ &
            0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
            1.0_dp, 0.0_dp, 1.0_dp, inf, &
            0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
            1.0_dp, inf, 1.0_dp, 0.0_dp ], [4,4], order=[2,1])
        order = [1,2,3,4]
        call check_close(tour_length(m, order), 4.0_dp, 1.0e-12_dp, 'tour length')
        call insertion_cost(m, [1,3], 2, delta)
        call check(size(delta) == 2, 'insertion cost size')
        call replace_infinite(m, replaced, ierr=ierr)
        call check(ierr == 0 .and. all(ieee_is_finite(replaced)), 'replace Inf')

        allocate(coords(3,2))
        coords(1,:) = [0.0_dp,0.0_dp]
        coords(2,:) = [3.0_dp,4.0_dp]
        coords(3,:) = [3.0_dp,0.0_dp]
        call euclidean_distance_matrix(coords, d)
        call check_close(d(1,2), 5.0_dp, 1.0e-12_dp, 'euclidean distance')
        call check_close(etsp_tour_length(coords,[1,2,3]), 12.0_dp, 1.0e-12_dp, 'ETSP tour length')
    end subroutine test_core

    subroutine test_heuristics()
        real(dp) :: m(4,4), inf, before, after
        integer, allocatable :: t(:), t2(:)
        type(tsp_tour) :: sol
        type(tsp_control) :: ctl
        integer :: meth

        inf = positive_infinity()
        m(1,:) = [0.0_dp,1.0_dp,0.0_dp,1.0_dp]
        m(2,:) = [1.0_dp,0.0_dp,1.0_dp,inf]
        m(3,:) = [0.0_dp,1.0_dp,0.0_dp,1.0_dp]
        m(4,:) = [1.0_dp,inf,1.0_dp,0.0_dp]

        do meth = tsp_identity, tsp_repetitive_nn
            if (meth == tsp_random) cycle
            sol = solve_tsp(m, meth)
            call check(size(sol%order) == 4, 'solver tour size')
            call check(is_valid_tour(sol%order,4), 'solver valid tour')
            call check(abs(sol%length-4.0_dp) < 1.0e-12_dp .or. .not. ieee_is_finite(sol%length), &
                'solver expected length')
        end do

        call random_permutation(4,t)
        before = tour_length(m,t)
        call replace_infinite_for_test(m)
        call two_opt(m,t,t2)
        after = tour_length(m,t2)
        call check(after <= before .or. .not. ieee_is_finite(before), 'two opt does not worsen')
        call two_opt_symmetric(m,t,t2)
        after = tour_length(m,t2)
        call check(after <= before .or. .not. ieee_is_finite(before), 'symmetric two opt does not worsen')

        ctl = tsp_control()
        ctl%maxit = 100
        ctl%tmax = 5
        sol = solve_tsp(m, tsp_sa_method, ctl)
        call check(is_valid_tour(sol%order,4), 'SA valid tour')
    end subroutine test_heuristics

    subroutine replace_infinite_for_test(m)
        real(dp), intent(inout) :: m(:,:)
        real(dp), allocatable :: work(:,:)
        integer :: ierr
        call replace_infinite(m,work,ierr=ierr)
        if (ierr == 0) m = work
    end subroutine replace_infinite_for_test

    subroutine test_atsp_two_opt()
        real(dp) :: a(5,5), before, after
        integer :: initial(5)
        integer, allocatable :: improved(:)
        integer :: i, j
        real(dp) :: candidate
        integer, allocatable :: rev(:)

        a = reshape([ &
            0.139303529169410_dp,0.897691324818879_dp,0.509101516567171_dp,0.430898967897519_dp,0.141799068776891_dp, &
            0.033456290373579_dp,0.902805947931483_dp,0.203576791565865_dp,0.435874363640323_dp,0.064170722616836_dp, &
            0.101683554705232_dp,0.631239329231903_dp,0.555331876967102_dp,0.082961557200179_dp,0.272443652851507_dp, &
            0.215095571940765_dp,0.532841097796336_dp,0.795302660670131_dp,0.432568762451410_dp,0.582661165855825_dp, &
            0.250269076088443_dp,0.164849652675912_dp,0.638499777996913_dp,0.857200765516609_dp,0.013439181726426_dp ], &
            [5,5])
        initial = [1,2,3,4,5]
        before = tour_length(a,initial)
        call two_opt(a,initial,improved)
        after = tour_length(a,improved)
        call check(after <= before + 1.0e-10_dp, 'ATSP two opt improves')

        do i = 2, 4
            do j = i+1, 5
                allocate(rev,source=improved)
                rev(i:j) = improved(j:i:-1)
                candidate = tour_length(a,rev)
                call check(candidate >= after - 1.0e-7_dp, 'two opt local optimality')
                deallocate(rev)
            end do
        end do
    end subroutine test_atsp_two_opt

    subroutine test_transformations()
        real(dp) :: a(3,3)
        real(dp), allocatable :: d(:,:), t(:,:)
        integer, allocatable :: filtered(:), path(:)
        integer :: ierr
        type(tsp_path_collection) :: paths

        a = reshape([0.0_dp,4.0_dp,2.0_dp, 1.0_dp,0.0_dp,3.0_dp, 5.0_dp,2.0_dp,0.0_dp],[3,3])
        call insert_dummy(a,2,d,const=0.0_dp)
        call check(all(shape(d) == [5,5]), 'dummy shape')
        call check(.not. ieee_is_finite(d(4,5)), 'dummy separation')

        call reformulate_atsp_as_tsp(a,t,infeasible=100.0_dp,cheap=-10.0_dp,ierr=ierr)
        call check(ierr == 0 .and. all(shape(t) == [6,6]), 'ATSP reformulation size')
        call check(is_symmetric_matrix(t), 'ATSP reformulation symmetric')
        call filter_atsp_tour([1,4,2,5,3,6],a,filtered)
        call check(size(filtered)==3 .and. is_valid_tour(filtered,3), 'filter ATSP dummies')

        call cut_tour_single([1,2,3,4],1,path)
        call check(all(path == [2,3,4]), 'single cut')
        call cut_tour_multiple([1,2,3,4,5,6],[2,5],paths)
        call check(size(paths%path)==2, 'multiple cuts count')
        call check(size(paths%path(1)%city)+size(paths%path(2)%city)==4, 'multiple cuts total')
    end subroutine test_transformations

    subroutine test_tsplib_distances()
        real(dp) :: att(4,2), geo(3,2)
        real(dp), allocatable :: d(:,:)

        att(:,1) = [0.0_dp,3.0_dp,6.0_dp,9.0_dp]
        att(:,2) = [0.0_dp,4.0_dp,8.0_dp,12.0_dp]
        call tsplib_att_distance(att,d)
        call check(all(nint(d) == reshape([0,2,4,5, 2,0,2,4, 4,2,0,2, 5,4,2,0],[4,4],order=[2,1])), 'ATT')

        geo(1,:) = [48.12_dp,16.22_dp]
        geo(2,:) = [46.38_dp,14.18_dp]
        geo(3,:) = [48.18_dp,14.17_dp]
        call tsplib_geo_distance(geo,d)
        call check(all(nint(d) == reshape([0,234,155, 234,0,186, 155,186,0],[3,3],order=[2,1])), 'GEO')
    end subroutine test_tsplib_distances

    subroutine test_tsplib_roundtrip()
        real(dp) :: m(4,4)
        type(tsplib_instance) :: inst
        integer :: ierr
        character(len=*), parameter :: fname='test-roundtrip.tsp'

        m(1,:)=[0.0_dp,1.0_dp,2.0_dp,3.0_dp]
        m(2,:)=[1.0_dp,0.0_dp,4.0_dp,5.0_dp]
        m(3,:)=[2.0_dp,4.0_dp,0.0_dp,6.0_dp]
        m(4,:)=[3.0_dp,5.0_dp,6.0_dp,0.0_dp]
        call write_tsplib_tsp(fname,m,precision=3,ierr=ierr)
        call check(ierr==0,'write TSPLIB')
        call read_tsplib(fname,inst,precision=3,ierr=ierr)
        call check(ierr==0,'read TSPLIB')
        call check(maxval(abs(inst%cost-m))<1.0e-12_dp,'TSPLIB roundtrip')
        open(unit=99,file=fname,status='old',iostat=ierr)
        if (ierr==0) close(99,status='delete')
    end subroutine test_tsplib_roundtrip


    pure real(dp) function subset_cycle_length(cost, order) result(total)
        real(dp), intent(in) :: cost(:,:)
        integer, intent(in) :: order(:)
        integer :: i
        total = 0.0_dp
        if (size(order) <= 1) return
        do i = 1, size(order)-1
            total = total + cost(order(i),order(i+1))
        end do
        total = total + cost(order(size(order)),order(1))
    end function subset_cycle_length

    subroutine test_random_two_opt()
        integer, parameter :: cases = 100
        integer :: c, n, i, j
        real(dp), allocatable :: a(:,:)
        integer, allocatable :: p(:), q(:), z(:)
        real(dp) :: l0, l1, lz

        do c = 1, cases
            n = 3 + mod(c,7)
            allocate(a(n,n))
            call random_number(a)
            do i = 1, n
                a(i,i) = 0.0_dp
            end do
            call random_permutation(n,p)
            l0 = tour_length(a,p)
            call two_opt(a,p,q)
            l1 = tour_length(a,q)
            call check(l1 <= l0 + 1.0e-7_dp, 'random two-opt non-worsening')
            do i = 2, n - 1
                do j = i + 1, n
                    allocate(z, source=q)
                    z(i:j) = q(j:i:-1)
                    lz = tour_length(a,z)
                    call check(lz >= l1 - 1.0001e-7_dp, 'random two-opt local optimality')
                    deallocate(z)
                end do
            end do
            deallocate(a,p,q)
        end do
    end subroutine test_random_two_opt

    subroutine test_random_insertion()
        integer, parameter :: cases = 100
        integer :: c, n, m, k, pos, i
        real(dp), allocatable :: a(:,:), delta(:)
        integer, allocatable :: p(:), newp(:), out(:)
        real(dp) :: base, actual

        do c = 1, cases
            n = 4 + mod(c,6)
            allocate(a(n,n))
            call random_number(a)
            do i = 1, n
                a(i,i) = 0.0_dp
            end do
            call random_permutation(n,p)
            m = 2 + mod(c,n-2)
            k = p(m+1)
            call insertion_cost(a,p(:m),k,delta)
            base = subset_cycle_length(a,p(:m))
            do pos = 1, m
                allocate(newp(m+1))
                newp(:pos) = p(:pos)
                newp(pos+1) = k
                if (pos < m) newp(pos+2:) = p(pos+1:m)
                actual = subset_cycle_length(a,newp) - base
                call check(abs(actual-delta(pos)) <= 1.0e-10_dp, 'random insertion delta')
                deallocate(newp)
            end do
            call insertion_heuristic(a,tsp_nearest_insertion,out,start=1)
            call check(is_valid_tour(out,n),'random nearest insertion valid')
            call insertion_heuristic(a,tsp_farthest_insertion,out,start=1)
            call check(is_valid_tour(out,n),'random farthest insertion valid')
            call insertion_heuristic(a,tsp_cheapest_insertion,out,start=1)
            call check(is_valid_tour(out,n),'random cheapest insertion valid')
            call arbitrary_insertion(a,out)
            call check(is_valid_tour(out,n),'random arbitrary insertion valid')
            deallocate(a,p,delta,out)
        end do
    end subroutine test_random_insertion


end program test_tsp
