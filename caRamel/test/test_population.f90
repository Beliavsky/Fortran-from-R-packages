program test_population
    use caramel, only: dp, decrease_pop, seed_random
    implicit none
    integer, parameter :: n = 60, d = 3
    real(dp) :: obj(n,d), prec(d)
    logical :: sense(d)
    integer, allocatable :: ia(:), ip(:)
    integer :: trial, i, j

    call seed_random(8675309)
    prec = [0.05_dp,0.05_dp,0.05_dp]
    sense = [.true.,.false.,.true.]
    do trial = 1, 200
        call random_number(obj)
        call decrease_pop(obj, sense, prec, 15, 20, ia, ip)
        call check(size(ia) <= 15, 'archive cap')
        call check(size(ip) <= 20, 'population cap')
        call check(all(ia >= 1 .and. ia <= n), 'archive indices')
        call check(all(ip >= 1 .and. ip <= n), 'population indices')
        do i = 1, size(ia)
            do j = i + 1, size(ia)
                call check(ia(i) /= ia(j), 'unique archive')
            end do
            call check(.not. any(ip == ia(i)), 'archive/pop disjoint')
        end do
        do i = 1, size(ip)
            do j = i + 1, size(ip)
                call check(ip(i) /= ip(j), 'unique population')
            end do
        end do
    end do
    print '(a)', 'test_population: PASS'
contains
    subroutine check(condition, name)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: name
        if (.not. condition) then
            write(*,'(a)') 'FAIL: '//trim(name)
            error stop 1
        end if
    end subroutine check
end program test_population
