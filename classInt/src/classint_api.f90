! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Modern Fortran translation of R package classInt computational code.
module classint_api
    use classint_kinds, only: dp
    use classint_types, only: classint_options, class_intervals
    use classint_core, only: classint_fit
    use classint_metrics, only: find_cols
    implicit none
    private

    interface classify_intervals
        module procedure classify_intervals_fit
        module procedure classify_intervals_n
        module procedure classify_intervals_auto
    end interface classify_intervals

    public :: classify_intervals

contains

    subroutine classify_intervals_fit(fit, cols)
        type(class_intervals), intent(in) :: fit !! Existing interval fit whose retained observations are assigned to classes.
        integer, allocatable, intent(out) :: cols(:) !! One-based class labels; zero denotes a retained non-finite observation.

        call find_cols(fit, cols)
    end subroutine classify_intervals_fit

    subroutine classify_intervals_n(var, n, style, cols, fit, options, interval_closure)
        real(dp), intent(in) :: var(:) !! Numeric observations to fit and classify; non-finite values receive class label zero.
        integer, intent(in) :: n !! Requested class count for styles that use n; values above the unique count are reduced.
        character(len=*), intent(in) :: style !! Interval style accepted by classint_fit, such as fisher, quantile, or equal.
        integer, allocatable, intent(out) :: cols(:) !! One-based fitted class label for every input observation.
        type(class_intervals), intent(out), optional :: fit !! Optional fitted interval object retained by the caller for reuse.
        type(classint_options), intent(in), optional :: options !! Optional style-specific controls; defaults match classint_fit.
        character(len=*), intent(in), optional :: interval_closure !! Optional "left" or "right" closure; Jenks forces right.
        type(class_intervals) :: local_fit

        if (present(options)) then
            if (present(interval_closure)) then
                call classint_fit(var, n, style, local_fit, options, interval_closure)
            else
                call classint_fit(var, n, style, local_fit, options)
            end if
        else if (present(interval_closure)) then
            call classint_fit(var, n, style, local_fit, interval_closure=interval_closure)
        else
            call classint_fit(var, n, style, local_fit)
        end if
        call find_cols(local_fit, cols)
        if (present(fit)) fit = local_fit
    end subroutine classify_intervals_n

    subroutine classify_intervals_auto(var, style, cols, fit, options, interval_closure)
        real(dp), intent(in) :: var(:) !! Numeric observations classified using a Sturges-derived count where the style needs n.
        character(len=*), intent(in) :: style !! Interval style accepted by the n-free classint_fit convenience overload.
        integer, allocatable, intent(out) :: cols(:) !! One-based fitted class label for every input observation.
        type(class_intervals), intent(out), optional :: fit !! Optional fitted interval object retained by the caller for reuse.
        type(classint_options), intent(in), optional :: options !! Optional style-specific controls; defaults match classint_fit.
        character(len=*), intent(in), optional :: interval_closure !! Optional "left" or "right" closure; Jenks forces right.
        type(class_intervals) :: local_fit

        if (present(options)) then
            if (present(interval_closure)) then
                call classint_fit(var, style, local_fit, options, interval_closure)
            else
                call classint_fit(var, style, local_fit, options)
            end if
        else if (present(interval_closure)) then
            call classint_fit(var, style, local_fit, interval_closure=interval_closure)
        else
            call classint_fit(var, style, local_fit)
        end if
        call find_cols(local_fit, cols)
        if (present(fit)) fit = local_fit
    end subroutine classify_intervals_auto
end module classint_api
