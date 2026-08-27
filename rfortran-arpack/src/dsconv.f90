!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: dsconv
!
!\Description:
!  Convergence testing for the symmetric Arnoldi eigenvalue routine.
!
!\Usage:
!  call dsconv
!     ( N, RITZ, BOUNDS, TOL, NCONV )
!
!\Arguments
!  N       Integer.  (INPUT)
!          Number of Ritz values to check for convergence.
!
!  RITZ    Double precision array of length N.  (INPUT)
!          The Ritz values to be checked for convergence.
!
!  BOUNDS  Double precision array of length N.  (INPUT)
!          Ritz estimates associated with the Ritz values in RITZ.
!
!  TOL     Double precision scalar.  (INPUT)
!          Desired relative accuracy for a Ritz value to be considered
!          "converged".
!
!  NCONV   Integer scalar.  (OUTPUT)
!          Number of "converged" Ritz values.
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Routines called:
!     arscnd  ARPACK utility routine for timing.
!     dlamch  LAPACK routine that determines machine constants.
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\SCCS Information: @(#)
! FILE: sconv.F   SID: 2.4   DATE OF SID: 4/19/96   RELEASE: 2
!
!\Remarks
!     1. Starting with version 2.4, this routine no longer uses the
!        Parlett strategy using the gap conditions.
!
!\EndLib
!
!-----------------------------------------------------------------------
!
subroutine dsconv (n, ritz, bounds, tol, nconv)
!
!     %----------------------------------------------------%
!     | Include files for debugging and timing information |
!     %----------------------------------------------------%
!
include   'debug.h'
include   'stat.h'
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
integer    n, nconv
Double precision &
  &           tol
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
Double precision &
  &           ritz(n), bounds(n)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
integer    i
Double precision &
  &           temp, eps23
!
!     %-------------------%
!     | External routines |
!     %-------------------%
!
Double precision &
  &           dlamch
external   dlamch

!     %---------------------%
!     | Intrinsic Functions |
!     %---------------------%
!
intrinsic    abs
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
call arscnd (t0)
!
eps23 = dlamch('Epsilon-Machine')
eps23 = eps23**(2.0D+0 / 3.0D+0)
!
nconv  = 0
do 10 i = 1, n
!
!        %-----------------------------------------------------%
!        | The i-th Ritz value is considered "converged"       |
!        | when: bounds(i) .le. TOL*max(eps23, abs(ritz(i)))   |
!        %-----------------------------------------------------%
!
   temp = max( eps23, abs(ritz(i)) )
   if ( bounds(i) .le. tol*temp ) then
      nconv = nconv + 1
   end if
!
10 continue
!
call arscnd (t1)
tsconv = tsconv + (t1 - t0)
!
return
!
!     %---------------%
!     | End of dsconv |
!     %---------------%
!
end
