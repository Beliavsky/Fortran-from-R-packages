program test_core
    use, intrinsic :: iso_fortran_env, only: int64
    use caramel, only: dp, pareto, dominate, dominated, val2rank, boxes, vol_splx, seed_random, new_xval
    use caramel_delaunay, only: delaunay_nd
    implicit none
    real(dp) :: x(5,2), ranks(4), s(3,2), v(4)
    real(dp), allocatable :: xn(:,:), pc(:,:)
    real(dp) :: param(8,2), crit(8,2), sp(2), bounds(2,2)
    integer :: front(5), onion(5), i
    integer, allocatable :: simp(:,:)
    logical :: d(4), ok
    logical :: sense(2)
    integer(int64) :: box_id(3)

    call seed_random(12345)

    x = reshape([1.0_dp,0.0_dp, 0.0_dp,1.0_dp, 0.5_dp,0.5_dp, 0.0_dp,0.0_dp, 1.0_dp,0.0_dp], [5,2], order=[2,1])
    call pareto(x, front)
    call check(all(front == [1,1,1,0,0]), "pareto")
    call dominate(x, onion)
    call check(all(onion == [1,1,1,3,2]), "dominate")
    call dominated([1.0_dp,1.0_dp], reshape([0.5_dp,0.5_dp, 1.0_dp,1.0_dp, 2.0_dp,0.0_dp, 0.0_dp,2.0_dp], [4,2], order=[2,1]), d)
    call check(all(d .eqv. [.true.,.false.,.false.,.false.]), "dominated")

    v = [3.0_dp,1.0_dp,1.0_dp,2.0_dp]
    call val2rank(v, 1, ranks)
    call check(maxval(abs(ranks-[4.0_dp,1.5_dp,1.5_dp,3.0_dp])) < 1.0e-12_dp, "val2rank average")
    call val2rank(v, 2, ranks)
    call check(maxval(abs(ranks-[3.0_dp,1.0_dp,1.0_dp,2.0_dp])) < 1.0e-12_dp, "val2rank unique")
    call val2rank(v, 3, ranks)
    call check(maxval(abs(ranks-[4.0_dp,2.0_dp,2.0_dp,3.0_dp])) < 1.0e-12_dp, "val2rank max")

    call boxes(reshape([0.0_dp,0.0_dp, 0.1_dp,0.0_dp, 0.0_dp,0.2_dp], [3,2], order=[2,1]), &
               [0.1_dp,0.1_dp], box_id)
    call check(all(box_id == [4_int64,5_int64,7_int64]), "boxes")

    s = reshape([0.0_dp,0.0_dp, 1.0_dp,0.0_dp, 0.0_dp,1.0_dp], [3,2], order=[2,1])
    call check(abs(vol_splx(s)-0.5_dp) < 1.0e-12_dp, "simplex volume")

    call delaunay_nd(reshape([0.0_dp,0.0_dp, 1.0_dp,0.0_dp, 1.0_dp,1.0_dp, 0.0_dp,1.0_dp], [4,2], order=[2,1]), simp, ok)
    call check(ok, "delaunay status")
    call check(size(simp,2) == 3 .and. size(simp,1) >= 2, "delaunay shape")
    do i = 1, size(simp,1)
        call check(all(simp(i,:) >= 1 .and. simp(i,:) <= 4), "delaunay indices")
    end do

    do i = 1, 8
        param(i,1) = -1.0_dp + 2.0_dp * real(i-1,dp) / 7.0_dp
        param(i,2) = sin(real(i,dp))
        crit(i,1) = -(param(i,1)-0.2_dp)**2 - 0.1_dp*param(i,2)**2
        crit(i,2) = -(param(i,1)+0.3_dp)**2 - 0.2_dp*(param(i,2)-0.1_dp)**2
    end do
    sp = [0.2_dp,0.2_dp]
    bounds(1,:) = [-2.0_dp,2.0_dp]
    bounds(2,:) = [-2.0_dp,2.0_dp]
    sense = [.true.,.true.]
    call new_xval(param, crit, sense, sp, bounds, [2,2,2,2], .true., xn, pc)
    call check(size(xn,1) > 0 .and. size(xn,2) == 2, "new_xval shape")
    call check(all(xn(:,1) >= -2.0_dp .and. xn(:,1) <= 2.0_dp), "new_xval bounds 1")
    call check(all(xn(:,2) >= -2.0_dp .and. xn(:,2) <= 2.0_dp), "new_xval bounds 2")

    print '(a)', 'test_core: PASS'
contains
    subroutine check(condition, name)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: name
        if (.not. condition) then
            write(*,'(a)') 'FAIL: '//trim(name)
            error stop 1
        end if
    end subroutine check
end program test_core
