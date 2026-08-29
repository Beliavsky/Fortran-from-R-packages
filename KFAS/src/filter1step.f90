! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from KFAS 1.6.0 src/filter1step.f90 by Jouni Helske.
! Algorithm retained and mechanically modernized for free-form FPM packaging.
! See NOTICE.md and provenance/upstream-sha256.txt.


!diffuse filtering for single time point
subroutine dfilter1step(ymiss, yt, zt, ht, tt, rqr, at, pt, vt, ft,kt,&
pinf,finf,kinf,rankp,lik,basetol,c,p,m,i)
    use kfas_kinds, only: dp

    implicit none

    integer, intent(in) :: p, m
    integer, intent(inout) :: i,rankp
    integer, intent(in), dimension(p) :: ymiss
    real(dp), intent(in), dimension(p) :: yt
    real(dp), intent(in), dimension(m,p) :: zt
    real(dp), intent(in), dimension(p,p) :: ht
    real(dp), intent(in), dimension(m,m) :: tt
    real(dp), dimension(m,m) :: rqr
    real(dp), intent(in) :: basetol,c
    real(dp), intent(inout) :: lik
    real(dp), intent(inout), dimension(m) :: at
    real(dp), intent(inout), dimension(p) :: vt,ft,finf
    real(dp), intent(inout), dimension(m,p) :: kt,kinf
    real(dp), intent(inout), dimension(m,m) :: pt,pinf
    real(dp), dimension(m) :: filtered_state
    real(dp), dimension(m,m) :: filtered_covariance

    external :: dfilter1step2

    call dfilter1step2(ymiss, yt, zt, ht, tt, rqr, at, pt, vt, ft, kt, pinf, finf, kinf, &
        rankp, lik, basetol, c, p, m, i, filtered_state, filtered_covariance)

end subroutine dfilter1step

!non-diffuse filtering for single time point
subroutine filter1step(ymiss, yt, zt, ht, tt, rqr, at, pt, vt, &
    ft, kt, lik, basetol, c, p, m, j)
    use kfas_kinds, only: dp

    implicit none

    integer, intent(in) :: p, m,j
    integer, intent(in), dimension(p) :: ymiss
    real(dp), intent(in), dimension(p) :: yt
    real(dp), intent(in), dimension(m,p) :: zt
    real(dp), intent(in), dimension(p,p) :: ht
    real(dp), intent(in), dimension(m,m) :: tt
    real(dp), dimension(m,m) :: rqr
    real(dp), intent(in) :: basetol,c
    real(dp), intent(inout) :: lik
    real(dp), intent(inout), dimension(m) :: at
    real(dp), intent(inout), dimension(p) :: vt,ft
    real(dp), intent(inout), dimension(m,p) :: kt
    real(dp), intent(inout), dimension(m,m) :: pt
    real(dp), dimension(m) :: filtered_state
    real(dp), dimension(m,m) :: filtered_covariance

    external :: filter1step2

    call filter1step2(ymiss, yt, zt, ht, tt, rqr, at, pt, vt, ft, kt, lik, basetol, c, &
        p, m, j, filtered_state, filtered_covariance)

end subroutine filter1step

!! as above, but returns also att and ptt


!diffuse filtering for single time point
subroutine dfilter1step2(ymiss, yt, zt, ht, tt, rqr, at, pt, vt, ft,kt,&
pinf,finf,kinf,rankp,lik,basetol,c,p,m,i, att, ptt)
    use kfas_kinds, only: dp

    implicit none

    integer, intent(in) :: p, m
    integer :: k1, k2
    integer, intent(inout) :: i,rankp
    integer, intent(in), dimension(p) :: ymiss
    real(dp), intent(in), dimension(p) :: yt
    real(dp), intent(in), dimension(m,p) :: zt
    real(dp), intent(in), dimension(p,p) :: ht
    real(dp), intent(in), dimension(m,m) :: tt
    real(dp), dimension(m,m) :: rqr
    real(dp), intent(in) :: basetol,c
    real(dp), intent(inout) :: lik
    real(dp), intent(inout), dimension(m) :: at, att
    real(dp), intent(inout), dimension(p) :: vt,ft,finf
    real(dp), intent(inout), dimension(m,p) :: kt,kinf
    real(dp), intent(inout), dimension(m,m) :: pt,pinf, ptt
    real(dp), dimension(m,m) :: mm
    real(dp), dimension(m) :: ahelp
    real(dp) :: finv, tol

    real(dp), external :: ddot

    external dgemm, dsymm, dgemv, dsymv, dsyr, dsyr2
    tol = basetol * minval(abs(zt), mask = abs(zt) > 0.0_dp)**2

    do i = 1, p
        call dsymv('u',m,1.0_dp,pt,m,zt(:,i),1,0.0_dp,kt(:,i),1)
        ft(i) = ddot(m,zt(:,i),1,kt(:,i),1) + ht(i,i)
        if(ymiss(i) == 0) then
            call dsymv('u',m,1.0_dp,pinf,m,zt(:,i),1,0.0_dp,kinf(:,i),1)
            finf(i) = ddot(m,zt(:,i),1,kinf(:,i),1)
            vt(i) = yt(i) - ddot(m,zt(:,i),1,at,1)
            if (finf(i) > tol) then
                finv = 1.0_dp / finf(i)
                at = at + vt(i) * finv * kinf(:,i)
                call dsyr('u',m,ft(i) * finv**2,kinf(:,i),1,pt,m)
                call dsyr2('u',m,-finv,kt(:,i),1,kinf(:,i),1,pt,m)
                call dsyr('u',m,-finv,kinf(:,i),1,pinf,m)
                lik = lik - 0.5_dp * log(finf(i))
                rankp = rankp - 1
            else
                finf(i) = 0.0_dp
                if(ft(i) > tol) then
                    finv = 1.0_dp / ft(i)
                    at = at + vt(i) * finv * kt(:,i)
                    call dsyr('u',m,-finv,kt(:,i),1,pt,m)
                    lik = lik - c - 0.5_dp * (log(ft(i)) + vt(i)**2 * finv)
                end if
            end if
            if (ft(i) <= tol) then
                ft(i) = 0.0_dp
            end if
            if(rankp == 0 .and. i < p) then
                return
            end if
        end if
    end do
    att = at
    do k1 = 1,m
      do k2 = k1,m
        ptt(k1, k2) = pt(k1, k2)
        ptt(k2, k1) = ptt(k1, k2)
      end do
    end do
    call dgemv('n',m,m,1.0_dp,tt,m,at,1,0.0_dp,ahelp,1)
    at = ahelp
    call dsymm('r','u',m,m,1.0_dp,pt,m,tt,m,0.0_dp,mm,m)
    call dgemm('n','t',m,m,m,1.0_dp,mm,m,tt,m,0.0_dp,pt,m)
    pt = pt + rqr

    call dsymm('r','u',m,m,1.0_dp,pinf,m,tt,m,0.0_dp,mm,m)
    call dgemm('n','t',m,m,m,1.0_dp,mm,m,tt,m,0.0_dp,pinf,m)

    do i = 1, m
      if (pt(i,i) <= 0) then
        pt(i,:) = 0.0_dp
        pt(:,i) = 0.0_dp
      end if
      if (ptt(i,i) <= 0) then
        ptt(i,:) = 0.0_dp
        ptt(:,i) = 0.0_dp
      end if
      if (pinf(i,i) <= 0) then
        pinf(i,:) = 0.0_dp
        pinf(:,i) = 0.0_dp
      end if
    end do

end subroutine dfilter1step2

!non-diffuse filtering for single time point
subroutine filter1step2(ymiss, yt, zt, ht, tt, rqr, at, pt, vt, &
    ft, kt, lik, basetol, c, p, m, j, att, ptt)
    use kfas_kinds, only: dp

    implicit none

    integer, intent(in) :: p, m,j
    integer :: i, k1, k2
    integer, intent(in), dimension(p) :: ymiss
    real(dp), intent(in), dimension(p) :: yt
    real(dp), intent(in), dimension(m,p) :: zt
    real(dp), intent(in), dimension(p,p) :: ht
    real(dp), intent(in), dimension(m,m) :: tt
    real(dp), dimension(m,m) :: rqr
    real(dp), intent(in) :: basetol,c
    real(dp), intent(inout) :: lik
    real(dp), intent(inout), dimension(m) :: at, att
    real(dp), intent(inout), dimension(p) :: vt,ft
    real(dp), intent(inout), dimension(m,p) :: kt
    real(dp), intent(inout), dimension(m,m) :: pt, ptt
    real(dp), dimension(m,m) :: mm
    real(dp), dimension(m) :: ahelp
    real(dp) :: finv, tol

    real(dp), external :: ddot

    external dgemm, dsymm, dgemv, dsymv, dsyr
    tol = basetol * minval(abs(zt), mask = abs(zt) > 0.0_dp)**2

    do i = j + 1, p
        call dsymv('u',m,1.0_dp,pt,m,zt(:,i),1,0.0_dp,kt(:,i),1)
        ft(i) = ddot(m,zt(:,i),1,kt(:,i),1) + ht(i,i)
        if(ymiss(i) == 0) then
            vt(i) = yt(i) - ddot(m,zt(:,i),1,at,1)
            if (ft(i) > tol) then
                finv = 1.0_dp / ft(i)
                at = at + vt(i) * finv * kt(:,i)
                call dsyr('u',m,-finv,kt(:,i),1,pt,m)
                lik = lik - c - 0.5_dp * (log(ft(i)) + vt(i)**2 * finv)
            else
                ft(i) = 0.0_dp
            end if
        end if
    end do

    att = at
    do k1 = 1,m
      do k2 = k1,m
        ptt(k1, k2) = pt(k1, k2)
        ptt(k2, k1) = ptt(k1, k2)
      end do
    end do
    call dgemv('n',m,m,1.0_dp,tt,m,at,1,0.0_dp,ahelp,1)
    at = ahelp

    call dsymm('r','u',m,m,1.0_dp,pt,m,tt,m,0.0_dp,mm,m)
    call dgemm('n','t',m,m,m,1.0_dp,mm,m,tt,m,0.0_dp,pt,m)
    pt = pt + rqr

    do i = 1, m
      if (pt(i,i) <= 0) then
        pt(i,:) = 0.0_dp
        pt(:,i) = 0.0_dp
      end if
      if (ptt(i,i) <= 0) then
        ptt(i,:) = 0.0_dp
        ptt(:,i) = 0.0_dp
      end if
    end do

end subroutine filter1step2
