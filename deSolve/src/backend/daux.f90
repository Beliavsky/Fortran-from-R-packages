!  The code in this file is based on ODEPACK from netlib
!    https://www.netlib.org/odepack/
! 
!  Adapted for use in R package deSolve by the deSolve authors.


DOUBLE PRECISION FUNCTION D1MACH (IDUM)
INTEGER IDUM
!-----------------------------------------------------------------------
! THIS ROUTINE COMPUTES THE UNIT ROUNDOFF OF THE MACHINE IN DOUBLE
! PRECISION.  THIS IS DEFINED AS THE SMALLEST POSITIVE MACHINE NUMBER
! U SUCH THAT  1.0D0 + U .NE. 1.0D0 (IN DOUBLE PRECISION).
!-----------------------------------------------------------------------
DOUBLE PRECISION U, COMP
U = 1.0D0
10 U = U*0.5D0
COMP = 1.0D0 + U
IF (COMP .NE. 1.0D0) GO TO 10
D1MACH = U*2.0D0
RETURN
!----------------------- END OF FUNCTION D1MACH ------------------------
END
