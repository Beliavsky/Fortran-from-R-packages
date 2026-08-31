!
subroutine d_ope ( n, x, y, a, ja, ia )
!
! ope computes A * x for a sparse matrix A.
!
implicit none
!
integer n
!
double precision  a(*)
integer i
integer ia(n+1)
integer ja(*)
integer k1
integer k2
double precision  x(*)
double precision  y(*)
!
! spasrse matrix * vector multiplication
!
do i=1,n
k1 = ia(i)
k2 = ia(i+1) -1
y(i) = dot_product ( a(k1:k2), x(ja(k1:k2)) )
end do
!
return
end
!
!-----------------------------------------------------------------------
!
subroutine ds_eigen_f (maxnev, ncv, maxitr, &
& n, iwhich, &
& na, a, ja, ia, &
& v, d, &
& iparam)
!
implicit none
!
!     %--------------------%
!     | Input Declarations |
!     %--------------------%
!
integer          maxnev, ncv, maxitr, n, na, ja(*), ia(na+1), &
& iparam(8), iwhich
!
Double precision &
& a(*), &
& v(n, ncv), d(maxnev), &
& workl(ncv*(ncv+8)), workd(3*n), resid(n)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
Double precision &
& tol, sigma
integer          ipntr(11)
!
character        bmat*1, which*2
integer          ido, info, lworkl, &
& ishfts, mode1, ierr
logical          rvec, select(ncv)
!
!     %------------%
!     | Parameters |
!     %------------%
!
Double precision &
& zero
parameter        (zero = 0.0D+0)
!
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
!     %-------------------------------------------------%
!     | The following include statement and assignments |
!     | initiate trace output from the internal         |
!     | actions of ARPACK.  See debug.doc in the        |
!     | DOCUMENTS directory for usage.  Initially, the  |
!     | most useful information will be a breakdown of  |
!     | time spent in the various stages of computation |
!     | given by setting msaupd = 1.                    |
!     %-------------------------------------------------%
!
!     %-------------------------------------------------%
!     | The following sets dimensions for this problem. |
!     %-------------------------------------------------%
!
bmat  = 'I'
!
if (iwhich .eq. 1) then
which = 'LM'
else if (iwhich .eq. 2) then
which = 'SM'
else if (iwhich .eq. 7) then
which = 'LA'
else if (iwhich .eq. 8) then
which = 'SA'
else if (iwhich .eq. 9) then
which = 'BE'
else
!
goto 9000
end if
!
!     %-----------------------------------------------------%
!     |                                                     |
!     | Specification of stopping rules and initial         |
!     | conditions before calling DSAUPD                    |
!     |                                                     |
!     | TOL  determines the stopping criterion.             |
!     |                                                     |
!     |      Expect                                         |
!     |           abs(lambdaC - lambdaT) < TOL*abs(lambdaC) |
!     |               computed   true                       |
!     |                                                     |
!     |      If TOL .le. 0,  then TOL <- macheps            |
!     |           (machine precision) is used.              |
!     |                                                     |
!     | IDO  is the REVERSE COMMUNICATION parameter         |
!     |      used to specify actions to be taken on return  |
!     |      from DSAUPD. (See usage below.)                |
!     |                                                     |
!     |      It MUST initially be set to 0 before the first |
!     |      call to DSAUPD.                                |
!     |                                                     |
!     | INFO on entry specifies starting vector information |
!     |      and on return indicates error codes            |
!     |                                                     |
!     |      Initially, setting INFO=0 indicates that a     |
!     |      random starting vector is requested to         |
!     |      start the ARNOLDI iteration.  Setting INFO to  |
!     |      a nonzero value on the initial call is used    |
!     |      if you want to specify your own starting       |
!     |      vector (This vector must be placed in RESID.)  |
!     |                                                     |
!     | The work array WORKL is used in DSAUPD as           |
!     | workspace.  Its dimension LWORKL is set as          |
!     | illustrated below.                                  |
!     |                                                     |
!     %-----------------------------------------------------%
!
lworkl = ncv * (ncv + 8)
tol = zero
info = 0
ido = 0
!
!     %---------------------------------------------------%
!     | Specification of Algorithm Mode:                  |
!     |                                                   |
!     | This program uses the exact shift strategy        |
!     | (indicated by setting PARAM(1) = 1).              |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 1 of DSAUPD is used     |
!     | (IPARAM(7) = 1). All these options can be changed |
!     | by the user. For details see the documentation in |
!     | DSAUPD.                                           |
!     %---------------------------------------------------%
!
ishfts = 1
mode1 = 1
!
iparam(1) = ishfts
iparam(3) = maxitr
iparam(7) = mode1
!
!     %------------------------------------------------%
!     | M A I N   L O O P (Reverse communication loop) |
!     %------------------------------------------------%
!
10 continue
!
!        %---------------------------------------------%
!        | Repeatedly call the routine DSAUPD and take |
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
!
!
call dsaupd ( ido, bmat, n, which, maxnev, tol, resid, &
& ncv, v, n, iparam, ipntr, workd, workl, &
& lworkl, info )
!
if (ido .eq. -1 .or. ido .eq. 1) then
!
!           %--------------------------------------%
!           | Perform matrix vector multiplication |
!           |              y <--- OP*x             |
!           | The user should supply his/her own   |
!           | matrix vector multiplication routine |
!           | here that takes workd(ipntr(1)) as   |
!           | the input, and return the result to  |
!           | workd(ipntr(2)).                     |
!           %--------------------------------------%
!
call d_ope (na, workd(ipntr(1)), workd(ipntr(2)), &
& a, ja, ia)
!
!           %-----------------------------------------%
!           | L O O P   B A C K to call DSAUPD again. |
!           %-----------------------------------------%
!
go to 10
!
end if
!
!     %----------------------------------------%
!     | Either we have convergence or there is |
!     | an error.                              |
!     %----------------------------------------%
!
if ( info .lt. 0 ) then
!
!        %--------------------------%
!        | Error message. Check the |
!        | documentation in DSAUPD. |
!        %--------------------------%
!
goto 9000
!
else
!
!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using DSEUPD.                |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |
!        |                                           |
!        | Eigenvectors may be also computed now if  |
!        | desired.  (indicated by rvec = .true.)    |
!        |                                           |
!        | The routine DSEUPD now called to do this  |
!        | post processing (Other modes may require  |
!        | more complicated post processing than     |
!        | mode1.)                                   |
!        |                                           |
!        %-------------------------------------------%
!
rvec = .true.
!
call dseupd ( rvec, 'A', select, d, v, n, sigma, &
& bmat, n, which, maxnev, tol, resid, ncv, v, n, &
& iparam, ipntr, workd, workl, lworkl, ierr )
!
!         %----------------------------------------------%
!         | Eigenvalues are returned in the first column |
!         | of the two dimensional array D and the       |
!         | corresponding eigenvectors are returned in   |
!         | the first NCONV (=IPARAM(5)) columns of the  |
!         | two dimensional array V if requested.        |
!         | Otherwise, an orthogonal basis for the       |
!         | invariant subspace corresponding to the      |
!         | eigenvalues in D is returned in V.           |
!         %----------------------------------------------%
!
if ( ierr .ne. 0) then
!
!            %------------------------------------%
!            | Error condition:                   |
!            | Check the documentation of DSEUPD. |
!            %------------------------------------%
!
goto 9000
!
end if
!
end if
!
9000 continue

end
!
!
