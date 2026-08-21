! SPDX-License-Identifier: GPL-2.0-or-later
module adequacy_interfaces
    use adequacy_kinds, only: dp
    implicit none
    private

    public :: objective_fn, density_fn, cdf_fn

    abstract interface
        function objective_fn(par, data) result(value)
            import dp
            real(dp), intent(in) :: par(:)
            real(dp), intent(in) :: data(:)
            real(dp) :: value
        end function objective_fn

        function density_fn(par, x) result(value)
            import dp
            real(dp), intent(in) :: par(:)
            real(dp), intent(in) :: x
            real(dp) :: value
        end function density_fn

        function cdf_fn(par, x) result(value)
            import dp
            real(dp), intent(in) :: par(:)
            real(dp), intent(in) :: x
            real(dp) :: value
        end function cdf_fn
    end interface
end module adequacy_interfaces
