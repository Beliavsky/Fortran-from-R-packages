!  The code in this file is was taken from 
!    https://www.netlib.org/odepack/
!  Original author: Hindmarsh, Alan C. (LLNL)
!  Rewritten to be used with R by Karline Soetaert.
!

subroutine rprintd1(msg, d1)
character (len=*) msg
double precision DBL(1), d1
  DBL(1) = d1
  call dblepr(msg, -1, DBL, 1)
end subroutine

subroutine rprintd2(msg, d1, d2)
character (len=*) msg
double precision DBL(2), d1, d2
  DBL(1) = d1
  DBL(2) = d2
  call dblepr(msg, -1, DBL, 2)
end subroutine

subroutine rprinti1(msg, i1)
character (len=*) msg
integer I(1), i1
  I(1) = i1
  call intpr(msg, -1, I, 1)
end subroutine


!DECK XERRWD
SUBROUTINE XERRWD (MSG, NMES, NERR, LEVEL, NI, I1, I2, NR, R1, R2)

!***PURPOSE  Write error message with values.
!***original AUTHOR  Hindmarsh, Alan C., (LLNL)
! Rewritten to be used with R by Karline Soetaert
!
!  All arguments are input arguments.
!
!  MSG    = The message (character array).
!  NMES   = The length of MSG (number of characters).
!  NERR   = The error number (not used).
!  LEVEL  = The error level..
!           0 or 1 means recoverable (control returns to caller).
!           2 means fatal (run is aborted--see note below).
!  NI     = Number of integers (0, 1, or 2) to be printed with message.
!  I1,I2  = Integers to be printed, depending on NI.
!  NR     = Number of reals (0, 1, or 2) to be printed with message.
!  R1,R2  = Reals to be printed, depending on NR.
!
!-----------------------------------------------------------------------
!
!  Declare arguments.
!
DOUBLE PRECISION R1, R2, RVEC(2), Dummy
INTEGER NMES, NERR, LEVEL, NI, I1, I2, NR, Ivec(2)
CHARACTER(LEN=*) MSG

dummy = 0.d0
MSG = MSG(1:NMES) // char(0)
call rprintf(MSG)


IF (NI .EQ. 1) THEN
  MSG = 'In above message, I1 = %d' // char(0)
  call rprintfi1(MSG, I1)
  MSG = ' ' // char(0)
  call rprintf(MSG)
ENDIF

IF (NI .EQ. 2) THEN
 IVEC(1) = I1
 IVEC(2) = I2
  MSG = 'In above message, I1 = %d, I2 = %d' // char(0)
  call rprintfi2(MSG, I1, I2)
  MSG = ' ' // char(0)
  call rprintf(MSG)
ENDIF

IF (NR .EQ. 1) THEN
  MSG = 'In above message, R1 = %g' // char(0)
  call rprintfd1(MSG, R1)
  MSG = ' ' // char(0)
  call rprintf(MSG)
ENDIF

IF (NR .EQ. 2) THEN
 RVEC(1) = R1
 RVEC(2) = R2
  MSG = 'In above message, R1 = %g, R2 = %g' // char(0)
  call rprintfd2(MSG, R1, R2)
  MSG = ' ' // char(0)
  call rprintf(MSG)
ENDIF

!  Abort the run if LEVEL = 2.
 if (LEVEL .EQ. 2) call rexit ("fatal error")
RETURN

END
