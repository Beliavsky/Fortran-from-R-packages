! SPDX-License-Identifier: Artistic-2.0
! Derived from the R matlab package; see COPYRIGHTS and upstream/.
module matlab_timer
    use matlab_kinds, only : dp
    implicit none
    private

    integer, save :: saved_count = 0
    integer, save :: saved_rate = 0

    public :: tic
    public :: toc

contains

    subroutine tic()
        call system_clock(saved_count, saved_rate)
    end subroutine tic

    function toc(echo) result(seconds)
        logical, intent(in), optional :: echo
        real(dp) :: seconds
        integer :: now, rate
        logical :: do_echo

        call system_clock(now, rate)
        if (saved_rate > 0) rate = saved_rate
        seconds = real(now - saved_count, dp) / real(rate, dp)
        do_echo = .true.
        if (present(echo)) do_echo = echo
        if (do_echo) write(*, '(a,f0.6,a)') 'elapsed time is ', seconds, ' seconds'
    end function toc
end module matlab_timer
