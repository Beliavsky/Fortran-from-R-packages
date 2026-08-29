! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from KFAS 1.6.0 src/smoothonestep.f90 by Jouni Helske.
! Algorithm retained and mechanically modernized for free-form FPM packaging.
! See NOTICE.md and provenance/upstream-sha256.txt.

!non-diffuse smoothing for single time point
subroutine smooth1step(ymiss, zt, ht, tt, rtv, qt, vt, ft,kt,&
im,p,m,r,j,rt,etahat,epshat,needeps)
    use kfas_kinds, only: dp

    implicit none

    logical, intent(in) :: needeps
    integer, intent(in) :: p, m,j,r
    integer :: i
    integer, intent(in), dimension(p) :: ymiss
    real(dp), intent(in), dimension(m,p) :: zt
    real(dp), intent(in), dimension(p,p) :: ht
    real(dp), intent(in), dimension(m,m) :: tt
    real(dp), intent(in), dimension(m,r) :: rtv
    real(dp), intent(in), dimension(r,r) :: qt
    real(dp), intent(inout), dimension(m) :: rt
    real(dp), intent(inout), dimension(r) :: etahat
    real(dp), intent(inout), dimension(p) :: epshat
    real(dp), intent(in), dimension(p) :: vt,ft
    real(dp), intent(in), dimension(m,p) :: kt
    real(dp), intent(in), dimension(m,m) :: im
    real(dp), dimension(m,m) :: l0
    real(dp), dimension(m) :: mhelp
    real(dp), dimension(r) :: rhelp
    real(dp) :: finv

    real(dp), external :: ddot

    external dgemv, dsymv, dger

    call dgemv('t',m,r,1.0_dp,rtv,m,rt,1,0.0_dp,rhelp,1)
    call dsymv('l',r,1.0_dp,qt,r,rhelp,1,0.0_dp,etahat,1)
    call dgemv('t',m,m,1.0_dp,tt,m,rt,1,0.0_dp,mhelp,1)
    rt = mhelp
    do i = p, j, -1
        if(ymiss(i) == 0) then
            if(ft(i) > 0.0_dp) then
                finv = 1.0_dp / ft(i)
                if(needeps) then
                    epshat(i) = ht(i,i) * (vt(i) - ddot(m,kt(:,i),1,rt,1)) * finv
                end if
                l0 = im
                call dger(m,m,-finv,kt(:,i),1,zt(:,i),1,l0,m)
                call dgemv('t',m,m,1.0_dp,l0,m,rt,1,0.0_dp,mhelp,1)
                rt = mhelp + vt(i) * finv * zt(:,i)
            end if
        end if
    end do

end subroutine smooth1step

subroutine dsmooth1step(ymiss, zt, ht, tt, rtv, qt, vt, ft,kt,&
im,p,m,r,j,rt,rt1,finf,kinf,etahat,epshat,needeps)
    use kfas_kinds, only: dp

    implicit none

    logical, intent(in) :: needeps
    integer, intent(in) :: p, m,j,r
    integer :: i
    integer, intent(in), dimension(p) :: ymiss
    real(dp), intent(in), dimension(m,p) :: zt
    real(dp), intent(in), dimension(p,p) :: ht
    real(dp), intent(in), dimension(m,m) :: tt
    real(dp), intent(in), dimension(m,r) :: rtv
    real(dp), intent(in), dimension(r,r) :: qt
    real(dp), intent(inout), dimension(m) :: rt,rt1
    real(dp), intent(inout), dimension(r) :: etahat
    real(dp), intent(inout), dimension(p) :: epshat
    real(dp), intent(in), dimension(p) :: vt,ft,finf
    real(dp), intent(in), dimension(m,p) :: kt,kinf
    real(dp), intent(in), dimension(m,m) :: im
    real(dp), dimension(m,m) :: l0,linf
    real(dp), dimension(m) :: mhelp
    real(dp), dimension(r) :: rhelp
    real(dp) :: finv

    real(dp), external :: ddot

    external dgemv, dsymv, dger

    if(j == p) then
        call dgemv('t',m,r,1.0_dp,rtv,m,rt,1,0.0_dp,rhelp,1)
        call dsymv('l',r,1.0_dp,qt,r,rhelp,1,0.0_dp,etahat,1)
        call dgemv('t',m,m,1.0_dp,tt,m,rt,1,0.0_dp,mhelp,1)
        rt = mhelp
        call dgemv('t',m,m,1.0_dp,tt,m,rt1,1,0.0_dp,mhelp,1)
        rt1 = mhelp
    end if

    do i = j, 1, -1
        if(ymiss(i) == 0) then
            if(finf(i) > 0.0_dp) then
                finv = 1.0_dp / finf(i)
                if(needeps) then
                    epshat(i) = -ht(i,i) * ddot(m,kinf(:,i),1,rt,1) * finv
                end if
                linf = im
                call dger(m,m,-finv,kinf(:,i),1,zt(:,i),1,linf,m)

                mhelp = kinf(:,i) * ft(i) * finv - kt(:,i)
                l0 = 0.0_dp
                call dger(m,m,finv,mhelp,1,zt(:,i),1,l0,m)

                call dgemv('t',m,m,1.0_dp,linf,m,rt1,1,0.0_dp,mhelp,1)
                rt1 = mhelp
                call dgemv('t',m,m,1.0_dp,l0,m,rt,1,1.0_dp,rt1,1)
                rt1 = rt1 + vt(i) * finv * zt(:,i)

                call dgemv('t',m,m,1.0_dp,linf,m,rt,1,0.0_dp,mhelp,1)
                rt = mhelp
            else
                if(ft(i) > 0.0_dp) then
                    finv = 1.0_dp / ft(i)
                    if(needeps) then
                        epshat(i) = ht(i,i) * (vt(i) - ddot(m,kt(:,i),1,rt,1)) * finv
                    end if

                    l0 = im
                    call dger(m,m,-finv,kt(:,i),1,zt(:,i),1,l0,m)

                    call dgemv('t',m,m,1.0_dp,l0,m,rt,1,0.0_dp,mhelp,1)
                    rt = mhelp + vt(i) * finv * zt(:,i)

                    call dgemv('t',m,m,1.0_dp,l0,m,rt1,1,0.0_dp,mhelp,1)
                    rt1 = mhelp
                end if
            end if
        end if
    end do

end subroutine dsmooth1step
