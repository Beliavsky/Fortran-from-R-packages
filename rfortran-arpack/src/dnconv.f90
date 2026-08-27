!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: dnconv
!
!\Description:
!  Convergence testing for the nonsymmetric Arnoldi eigenvalue routine.
!
!\Usage:
!  call dnconv
!     ( N, RITZR, RITZI, BOUNDS, TOL, NCONV )
!
!\Arguments
!  N       Integer.  (INPUT)
!          Number of Ritz values to check for convergence.
!
!  RITZR,  Double precision arrays of length N.  (INPUT)
!  RITZI   Real and imaginary parts of the Ritz values to be checked
!          for convergence.

!  BOUNDS  Double precision array of length N.  (INPUT)
!          Ritz estimates for the Ritz values in RITZR and RITZI.
!
!  TOL     Double precision scalar.  (INPUT)
!          Desired backward error for a Ritz value to be considered
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
!\Local variables:
!     xxxxxx  real
!
!\Routines called:
!     arscnd  ARPACK utility routine for timing.
!     dlamch  LAPACK routine that determines machine constants.
!     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\Revision history:
!     xx/xx/92: Version ' 2.1'
!
!\SCCS Information: @(#)
! FILE: nconv.F   SID: 2.3   DATE OF SID: 4/20/96   RELEASE: 2
!
!\Remarks
!     1. xxxx
!
!\EndLib
!
!-----------------------------------------------------------------------
!
subroutine dnconv (n, ritzr, ritzi, bounds, tol, nconv)
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

Double precision &
  &           ritzr(n), ritzi(n), bounds(n)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
integer    i
Double precision &
  &           temp, eps23
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
Double precision &
  &           dlapy2, dlamch
external   dlapy2, dlamch

!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
!     %-------------------------------------------------------------%
!     | Convergence test: unlike in the symmetric code, I am not    |
!     | using things like refined error bounds and gap condition    |
!     | because I don't know the exact equivalent concept.          |
!     |                                                             |
!     | Instead the i-th Ritz value is considered "converged" when: |
!     |                                                             |
!     |     bounds(i) .le. ( TOL * | ritz | )                       |
!     |                                                             |
!     | for some appropriate choice of norm.                        |
!     %-------------------------------------------------------------%
!
call arscnd (t0)
!
!     %---------------------------------%
!     | Get machine dependent constant. |
!     %---------------------------------%
!
eps23 = dlamch('Epsilon-Machine')
eps23 = eps23**(2.0D+0 / 3.0D+0)
!
nconv  = 0
do 20 i = 1, n
   temp = max( eps23, dlapy2( ritzr(i), ritzi(i) ) )
   if (bounds(i) .le. tol*temp)   nconv = nconv + 1
20 continue
!
call arscnd (t1)
tnconv = tnconv + (t1 - t0)
!
return
!
!     %---------------%
!     | End of dnconv |
!     %---------------%
!
end
