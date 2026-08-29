! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from KFAS 1.6.0 src/kfilter.f90 by Jouni Helske.
! Algorithm retained and mechanically modernized for free-form FPM packaging.
! See NOTICE.md and provenance/upstream-sha256.txt.

! Compatibility entry point for filtering without saved filtered states.
! kfilter2 implements the common algorithm and returns the additional states.
subroutine kfilter(yt, ymiss, timevar, zt, ht, tt, rt, qt, a1, p1, p1inf, p, n, m, r, d, j, &
    at, pt, vt, ft, kt, pinf, finf, kinf, lik, tol, rankp, theta, thetavar, filtersignal)
    use kfas_kinds, only: dp

    implicit none

    integer, intent(in) :: p, n, m, r, filtersignal
    integer, intent(inout) :: d, j, rankp
    integer, intent(in), dimension(5) :: timevar
    integer, intent(in), dimension(n, p) :: ymiss
    real(dp), intent(in), dimension(n, p) :: yt
    real(dp), intent(in), dimension(p, m, (n - 1) * timevar(1) + 1) :: zt
    real(dp), intent(in), dimension(p, p, (n - 1) * timevar(2) + 1) :: ht
    real(dp), intent(in), dimension(m, m, (n - 1) * timevar(3) + 1) :: tt
    real(dp), intent(in), dimension(m, r, (n - 1) * timevar(4) + 1) :: rt
    real(dp), intent(in), dimension(r, r, (n - 1) * timevar(5) + 1) :: qt
    real(dp), intent(in), dimension(m) :: a1
    real(dp), intent(in), dimension(m, m) :: p1, p1inf
    real(dp), intent(in) :: tol
    real(dp), intent(inout), dimension(m, n + 1) :: at
    real(dp), intent(inout), dimension(m, m, n + 1) :: pt, pinf
    real(dp), intent(inout), dimension(p, n) :: vt, ft, finf
    real(dp), intent(inout), dimension(m, p, n) :: kt, kinf
    real(dp), intent(inout) :: lik
    real(dp), intent(inout), dimension(p, p, n) :: thetavar
    real(dp), intent(inout), dimension(n, p) :: theta
    real(dp), dimension(m, n) :: filtered_state
    real(dp), dimension(m, m, n) :: filtered_covariance

    external :: kfilter2

    call kfilter2(yt, ymiss, timevar, zt, ht, tt, rt, qt, a1, p1, p1inf, p, n, m, r, d, j, &
        at, pt, vt, ft, kt, pinf, finf, kinf, lik, tol, rankp, theta, thetavar, filtersignal, &
        filtered_state, filtered_covariance)

end subroutine kfilter
