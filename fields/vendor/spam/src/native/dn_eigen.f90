!
subroutine dn_eigen_f(maxnev, ncv, maxitr, &
& n, iwhich, &
& na, a, ja, ia, &
& v, dr, di, iparam)
!
implicit none
!
integer           maxnev, ncv, n, na, &
& iwhich, maxitr
!     %--------------%
!     | Local Arrays |
!     %--------------%
!
integer           iparam(8), ipntr(14), &
& ja(*), ia(na+1)
!
logical           select(ncv)
!
Double precision &
& dr(maxnev+1), di(maxnev+1), resid(n), &
& v(n, ncv), workd(3*n), &
& workev(3*ncv), &
& workl(3*ncv**2+6*ncv), &
& a(*)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
character         bmat*1, which*2
integer           ido, lworkl, info, &
& ierr, ishfts, mode
Double precision &
& tol, sigmar, sigmai
!
!     %------------%
!     | Parameters |
!     %------------%
!
Double precision &
& zero
parameter         (zero = 0.0D+0)
!
bmat   = 'I'
!
lworkl = 3*ncv*ncv+6*ncv
tol    = zero
ido    = 0
info   = 0
!
ishfts = 1
mode   = 1
!
iparam(1) = ishfts
iparam(3) = maxitr
iparam(7) = mode
!
if (iwhich .eq. 1) then
which = 'LM'
else if (iwhich .eq. 2) then
which = 'SM'
else if (iwhich .eq. 3) then
which = 'LR'
else if (iwhich .eq. 4) then
which = 'SR'
else if (iwhich .eq. 5) then
which = 'LI'
else if (iwhich .eq. 6) then
which = 'SI'
else
!
goto 9000
end if
!
!
10 continue
!
!        %---------------------------------------------%
!        | Repeatedly call the routine DNAUPD and take |
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
!
call dnaupd ( ido, bmat, n, which, maxnev, tol, resid, &
& ncv, v, n, iparam, ipntr, workd, workl, lworkl, &
& info )
!
if (ido .eq. -1 .or. ido .eq. 1) then
!
call d_ope (na, workd(ipntr(1)), workd(ipntr(2)), &
& a, ja, ia)
!
go to 10
!
end if
!
if ( info .lt. 0 ) then
!
goto 9000
!
else
!
call dneupd ( .true., 'A', select, dr, di, v, n, &
& sigmar, sigmai, workev, bmat, n, which, maxnev, tol, &
& resid, ncv, v, n, iparam, ipntr, workd, workl, &
& lworkl, ierr )
!
if ( ierr .lt. 0 ) then
!
goto 9000
!
end if
!
endif
!
9000 continue
!
end
!
